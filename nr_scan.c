/*
 * nr_scan.c — NR 频段频谱扫描工具 (USRP B210, UHD C API + FFTW)
 *
 * 与 nr_scan.py 逻辑一致:
 *   分段调谐扫过 n41/n78(或自定义频段), Welch 平均功率谱,
 *   多遍扫描按频点取 max-hold(抓 TDD 突发下行),
 *   以 底噪20%分位+threshold 判占用, 列出已占用频段(含峰值/归属),
 *   滑窗找空闲信道, 候选中心对齐合法 NR-ARFCN, 并给出最近 SSB 栅格点。
 *
 * 输出: 终端报告 + CSV; 若系统装有 gnuplot 则自动生成 PNG(同时保留 .gp 脚本)。
 *
 * 编译: make  (依赖 libuhd-dev, libfftw3-dev)
 */

#define _POSIX_C_SOURCE 199309L
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <complex.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <uhd.h>
#include <fftw3.h>

#define MAX_CHAN_BW 8
#define MAX_BANDS 8

typedef struct { double start, stop; const char *label; } alloc_t;
static const alloc_t ALLOC[] = {
    {2515e6, 2675e6, "移动 n41"},
    {3300e6, 3400e6, "3.3G 室内共享"},
    {3400e6, 3500e6, "电信 n78"},
    {3500e6, 3600e6, "联通 n78"},
    {4800e6, 4900e6, "移动 n79"},
    {4900e6, 4960e6, "广电 n79"},
};

typedef struct { const char *name; double start, stop; } band_t;
static const band_t PRESETS[] = {
    {"n41", 2515e6, 2675e6},
    {"n78", 3300e6, 3800e6},
    {"n79", 4400e6, 5000e6},
};

typedef struct {
    char preset[128];
    double start, stop;
    bool custom_range;
    char device[256];
    char antenna[64];
    double gain, rate, dwell, settle, pause, threshold, free_min;
    int fft, sweeps, top;
    double chan_bw[MAX_CHAN_BW];
    int n_chan_bw;
    char out[256];
} config_t;

typedef struct {
    double ch_db, center;
    long arfcn;
} cand_t;

typedef struct { size_t a, b; } seg_t;  /* 占用段 [a,b], fr/db 下标 */

static void ms_sleep(double sec) {
    if (sec <= 0) return;
    struct timespec ts = { (time_t)sec, (long)((sec - (time_t)sec) * 1e9) };
    nanosleep(&ts, NULL);
}

/* ---- NR 频点换算 (38.104) ---- */
static long freq_to_arfcn(double f) {
    return f < 3000e6 ? lround(f / 5e3) : 600000 + lround((f - 3000e6) / 15e3);
}
static double arfcn_to_freq(long n) {
    return n < 600000 ? n * 5e3 : 3000e6 + (n - 600000) * 15e3;
}
static void nearest_ssb(double f, double *ss_ref, long *gscn) {
    if (f < 3000e6) {
        double best_d = 1e18; long best_n = 1; int best_m = 1;
        for (int mi = 0; mi < 3; mi++) {
            int m = 1 + 2 * mi;
            long n = lround((f - m * 50e3) / 1.2e6);
            if (n < 1) n = 1;
            if (n > 2499) n = 2499;
            double ss = n * 1.2e6 + m * 50e3;
            double d = fabs(ss - f);
            if (d < best_d) { best_d = d; best_n = n; best_m = m; }
        }
        *ss_ref = best_n * 1.2e6 + best_m * 50e3;
        *gscn = 3 * best_n + (best_m - 3) / 2;
    } else {
        long n = lround((f - 3000e6) / 1.44e6);
        if (n < 0) n = 0;
        if (n > 14756) n = 14756;
        *ss_ref = 3000e6 + n * 1.44e6;
        *gscn = 7499 + n;
    }
}
static const char *alloc_label(double f) {
    for (size_t i = 0; i < sizeof(ALLOC) / sizeof(ALLOC[0]); i++)
        if (f >= ALLOC[i].start && f < ALLOC[i].stop) return ALLOC[i].label;
    return "未标注";
}

static void usage(const char *prog) {
    printf(
        "用法: %s [选项]\n"
        "  --preset LIST     预设频段, 逗号分隔: n41,n78,n79 (默认 n41,n78)\n"
        "  --start HZ        自定义起始频率 (与 --stop 一起覆盖 preset)\n"
        "  --stop HZ         自定义结束频率\n"
        "  --device ARGS     UHD 设备参数, 如 'serial=3129082'\n"
        "  --antenna ANT     接收天线口 (默认 RX2)\n"
        "  --gain DB         RX 增益 (默认 40)\n"
        "  --rate SPS        采样率 (默认 23.04e6, USB2 降到 ~8e6)\n"
        "  --fft N           FFT 点数 (默认 2048)\n"
        "  --dwell SEC       每个调谐点驻留时间 (默认 0.08)\n"
        "  --settle SEC      调谐稳定等待 (默认 0.02)\n"
        "  --sweeps N        每频段扫描遍数, max-hold 合并 (默认 3)\n"
        "  --pause SEC       两遍扫描间隔 (默认 1.0)\n"
        "  --threshold DB    占用门限: 底噪+N (默认 10)\n"
        "  --chan-bw LIST    候选带宽 Hz (默认 10e6,20e6,40e6)\n"
        "  --free-min R      候选内空闲 bin 比例下限 (默认 0.95)\n"
        "  --top N           每种带宽输出前 N 个候选 (默认 5)\n"
        "  --out PREFIX      输出文件前缀 (默认 scan_时间戳)\n"
        "  --help            显示本帮助\n", prog);
    exit(0);
}

static void parse_doubles(const char *s, double *out, int *n, int max) {
    char buf[256];
    snprintf(buf, sizeof(buf), "%s", s);
    for (char *t = strtok(buf, ","); t && *n < max; t = strtok(NULL, ","))
        out[(*n)++] = strtod(t, NULL);
}

static void parse_args(int argc, char **argv, config_t *c) {
    memset(c, 0, sizeof(*c));
    snprintf(c->preset, sizeof(c->preset), "n41,n78");
    snprintf(c->antenna, sizeof(c->antenna), "RX2");
    c->gain = 40; c->rate = 23.04e6; c->dwell = 0.08; c->settle = 0.02;
    c->pause = 1.0; c->threshold = 10; c->free_min = 0.95;
    c->fft = 2048; c->sweeps = 3; c->top = 5;

    static const struct option opts[] = {
        {"preset", required_argument, 0, 'p'}, {"start", required_argument, 0, 1},
        {"stop", required_argument, 0, 2},     {"device", required_argument, 0, 'd'},
        {"antenna", required_argument, 0, 'a'},{"gain", required_argument, 0, 'g'},
        {"rate", required_argument, 0, 'r'},   {"fft", required_argument, 0, 'F'},
        {"dwell", required_argument, 0, 'D'},  {"settle", required_argument, 0, 'S'},
        {"sweeps", required_argument, 0, 's'}, {"pause", required_argument, 0, 'P'},
        {"threshold", required_argument, 0, 't'},{"chan-bw", required_argument, 0, 'b'},
        {"free-min", required_argument, 0, 'f'},{"top", required_argument, 0, 'T'},
        {"out", required_argument, 0, 'o'},    {"help", no_argument, 0, 'h'}, {0, 0, 0, 0}
    };
    int idx, ch;
    bool have_start = false, have_stop = false, have_bw = false;
    while ((ch = getopt_long(argc, argv, "h", opts, &idx)) != -1) {
        switch (ch) {
        case 'p': snprintf(c->preset, sizeof(c->preset), "%s", optarg); break;
        case 1: c->start = strtod(optarg, NULL); have_start = true; break;
        case 2: c->stop = strtod(optarg, NULL); have_stop = true; break;
        case 'd': snprintf(c->device, sizeof(c->device), "%s", optarg); break;
        case 'a': snprintf(c->antenna, sizeof(c->antenna), "%s", optarg); break;
        case 'g': c->gain = strtod(optarg, NULL); break;
        case 'r': c->rate = strtod(optarg, NULL); break;
        case 'F': c->fft = atoi(optarg); break;
        case 'D': c->dwell = strtod(optarg, NULL); break;
        case 'S': c->settle = strtod(optarg, NULL); break;
        case 's': c->sweeps = atoi(optarg); break;
        case 'P': c->pause = strtod(optarg, NULL); break;
        case 't': c->threshold = strtod(optarg, NULL); break;
        case 'b': parse_doubles(optarg, c->chan_bw, &c->n_chan_bw, MAX_CHAN_BW); have_bw = true; break;
        case 'f': c->free_min = strtod(optarg, NULL); break;
        case 'T': c->top = atoi(optarg); break;
        case 'o': snprintf(c->out, sizeof(c->out), "%s", optarg); break;
        case 'h': default: usage(argv[0]);
        }
    }
    if (have_start && have_stop) c->custom_range = true;
    if (!have_bw) parse_doubles("10e6,20e6,40e6", c->chan_bw, &c->n_chan_bw, MAX_CHAN_BW);
    if (c->out[0] == '\0') {
        time_t now = time(NULL);
        struct tm tmv;
        localtime_r(&now, &tmv);
        strftime(c->out, sizeof(c->out), "scan_%Y%m%d_%H%M%S", &tmv);
    }
}

static uhd_usrp_handle g_usrp;
static void uhd_check(uhd_error err, const char *what) {
    if (err != UHD_ERROR_NONE) {
        char msg[512] = {0};
        if (g_usrp) uhd_usrp_last_error(g_usrp, msg, sizeof(msg));
        fprintf(stderr, "UHD 错误 @ %s: %s\n", what, msg);
        exit(1);
    }
}

/* 调谐到 freq, 采集 n_samps 个采样写入 out, 返回实际采样数, *ovf 累计溢出次数 */
static size_t capture_hop(uhd_rx_streamer_handle st, uhd_rx_metadata_handle md,
                          double freq, size_t n_samps, double settle,
                          float complex *out, float complex *pkt, size_t max_samps, int *ovf) {
    uhd_tune_request_t tr = {
        .target_freq = freq,
        .rf_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO,
        .dsp_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO,
        .args = "",
    };
    uhd_tune_result_t tres;
    uhd_check(uhd_usrp_set_rx_freq(g_usrp, &tr, 0, &tres), "set_rx_freq");
    ms_sleep(settle);

    uhd_stream_cmd_t cmd = {
        .stream_mode = UHD_STREAM_MODE_NUM_SAMPS_AND_DONE,
        .num_samps = n_samps,
        .stream_now = true,
    };
    uhd_check(uhd_rx_streamer_issue_stream_cmd(st, &cmd), "issue_stream_cmd");

    size_t got = 0;
    while (got < n_samps) {
        void *buffs[] = { pkt };
        size_t rx = 0;
        uhd_rx_streamer_recv(st, buffs, max_samps, &md, 1.0, false, &rx);
        uhd_rx_metadata_error_code_t ec;
        uhd_rx_metadata_error_code(md, &ec);
        if (ec == UHD_RX_METADATA_ERROR_CODE_OVERFLOW) { (*ovf)++; continue; }
        if (ec == UHD_RX_METADATA_ERROR_CODE_TIMEOUT) break;
        if (ec != UHD_RX_METADATA_ERROR_CODE_NONE) continue;
        size_t m = rx < n_samps - got ? rx : n_samps - got;
        memcpy(out + got, pkt, m * sizeof(float complex));
        got += m;
    }
    return got;
}

/* Welch 平均周期图, ps[k] 对应偏移 (k<nfft/2 ? k : k-nfft)*bin_hz, sum(ps)=平均功率 */
static void welch_psd(const float complex *iq, size_t len, int nfft,
                      const float *win, double w2, fftwf_plan plan,
                      float complex *fin, float complex *fout, double *ps) {
    size_t hop = (size_t)nfft / 2;
    size_t nseg = 1 + (len - (size_t)nfft) / hop;
    for (int i = 0; i < nfft; i++) ps[i] = 0;
    for (size_t s = 0; s < nseg; s++) {
        for (int i = 0; i < nfft; i++) fin[i] = iq[s * hop + i] * win[i];
        fftwf_execute(plan);
        for (int i = 0; i < nfft; i++) {
            double re = crealf(fout[i]), im = cimagf(fout[i]);
            ps[i] += re * re + im * im;
        }
    }
    double norm = 1.0 / ((double)nseg * nfft * nfft * w2);
    for (int i = 0; i < nfft; i++) ps[i] *= norm;
}

typedef struct {
    double *sum, *comb_db;  /* sum: 单遍线性累加; comb_db: 跨遍 max-hold (dBFS/Hz) */
    long *cnt;
    size_t nbins;
    double f0, bin_hz;
} grid_t;

static int cmp_double(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0;
}
static int cmp_cand(const void *a, const void *b) {
    double d = ((const cand_t *)a)->ch_db - ((const cand_t *)b)->ch_db;
    return d < 0 ? -1 : d > 0;
}

static double percentile(const double *x, size_t n, double q) {
    double *tmp = malloc(n * sizeof(double));
    memcpy(tmp, x, n * sizeof(double));
    qsort(tmp, n, sizeof(double), cmp_double);
    double r = tmp[(size_t)(q * (n - 1))];
    free(tmp);
    return r;
}

int main(int argc, char **argv) {
    config_t cfg;
    parse_args(argc, argv, &cfg);

    band_t bands[MAX_BANDS];
    int n_bands = 0;
    if (cfg.custom_range) {
        bands[n_bands++] = (band_t){"custom", cfg.start, cfg.stop};
    } else {
        char buf[128];
        snprintf(buf, sizeof(buf), "%s", cfg.preset);
        for (char *t = strtok(buf, ","); t && n_bands < MAX_BANDS; t = strtok(NULL, ",")) {
            bool found = false;
            for (size_t i = 0; i < sizeof(PRESETS) / sizeof(PRESETS[0]); i++)
                if (strcmp(t, PRESETS[i].name) == 0) { bands[n_bands++] = PRESETS[i]; found = true; }
            if (!found) { fprintf(stderr, "未知预设 %s\n", t); return 1; }
        }
    }

    /* ---- 设备初始化 ---- */
    uhd_check(uhd_usrp_make(&g_usrp, cfg.device), "usrp_make");
    char mboard[128] = {0};
    uhd_usrp_get_mboard_name(g_usrp, 0, mboard, sizeof(mboard));
    uhd_meta_range_handle gr;
    uhd_meta_range_make(&gr);
    double glo = 0, ghi = 0;
    if (uhd_usrp_get_rx_gain_range(g_usrp, "", 0, gr) == UHD_ERROR_NONE) {
        uhd_meta_range_start(gr, &glo);
        uhd_meta_range_stop(gr, &ghi);
    }
    uhd_meta_range_free(&gr);
    printf("设备: %s  增益范围: %.0f-%.0f dB\n", mboard, glo, ghi);

    uhd_check(uhd_usrp_set_rx_rate(g_usrp, cfg.rate, 0), "set_rx_rate");
    uhd_check(uhd_usrp_set_rx_gain(g_usrp, cfg.gain, 0, ""), "set_rx_gain");
    uhd_check(uhd_usrp_set_rx_antenna(g_usrp, cfg.antenna, 0), "set_rx_antenna");
    uhd_usrp_set_rx_bandwidth(g_usrp, cfg.rate, 0);

    size_t channel = 0;
    uhd_stream_args_t sargs = {
        .cpu_format = "fc32", .otw_format = "sc16", .args = "",
        .channel_list = &channel, .n_channels = 1,
    };
    uhd_rx_streamer_handle st;
    uhd_rx_streamer_make(&st);
    uhd_check(uhd_usrp_get_rx_stream(g_usrp, &sargs, st), "get_rx_stream");
    uhd_rx_metadata_handle md;
    uhd_rx_metadata_make(&md);
    size_t max_samps = 0;
    uhd_rx_streamer_max_num_samps(st, &max_samps);

    /* ---- FFT / 采集缓冲 ---- */
    int nfft = cfg.fft;
    float *win = malloc(nfft * sizeof(float));
    double w2 = 0;
    for (int i = 0; i < nfft; i++) {
        win[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (nfft - 1)));
        w2 += (double)win[i] * win[i];
    }
    w2 /= nfft;
    float complex *fin = fftwf_malloc(nfft * sizeof(float complex));
    float complex *fout = fftwf_malloc(nfft * sizeof(float complex));
    fftwf_plan plan = fftwf_plan_dft_1d(nfft, (fftwf_complex *)fin,
                                        (fftwf_complex *)fout, FFTW_FORWARD, FFTW_ESTIMATE);
    size_t n_samps = (size_t)(cfg.rate * cfg.dwell);
    if (n_samps < (size_t)nfft * 4) n_samps = (size_t)nfft * 4;
    float complex *iq = malloc(n_samps * sizeof(float complex));
    float complex *pkt = malloc(max_samps * sizeof(float complex));
    double *ps = malloc(nfft * sizeof(double));

    double bin_hz = cfg.rate / nfft;
    double usable = 0.7 * cfg.rate;
    double notch = cfg.rate * 0.05 < 0.6e6 ? cfg.rate * 0.05 : 0.6e6;

    /* ---- 逐频段扫描 ---- */
    for (int bi = 0; bi < n_bands; bi++) {
        double fs = bands[bi].start, fe = bands[bi].stop;
        /* 调谐点 */
        double centers[512];
        int nc = 0;
        if (fe - fs <= usable) {
            centers[nc++] = (fs + fe) / 2;
        } else {
            double step = 0.75 * usable;
            for (double c = fs + usable / 2; c + usable / 2 < fe && nc < 510; c += step)
                centers[nc++] = c;
            centers[nc++] = fe - usable / 2;
        }
        double est = nc * (cfg.dwell + cfg.settle + 0.05) * cfg.sweeps;
        printf("\n频段 %.0f-%.0f MHz: %d 个调谐点 x %d 遍, 预计 %.0f 秒\n",
               fs / 1e6, fe / 1e6, nc, cfg.sweeps, est);

        /* 频率网格 */
        grid_t g;
        g.bin_hz = bin_hz;
        g.f0 = floor(fs / bin_hz) * bin_hz;
        g.nbins = (size_t)ceil((fe - g.f0) / bin_hz) + 1;
        g.sum = malloc(g.nbins * sizeof(double));
        g.cnt = malloc(g.nbins * sizeof(long));
        g.comb_db = malloc(g.nbins * sizeof(double));
        bool *covered = calloc(g.nbins, sizeof(bool));
        for (size_t i = 0; i < g.nbins; i++) g.comb_db[i] = -INFINITY;

        int total_ovf = 0;
        for (int sw = 0; sw < cfg.sweeps; sw++) {
            memset(g.sum, 0, g.nbins * sizeof(double));
            memset(g.cnt, 0, g.nbins * sizeof(long));
            for (int hi = 0; hi < nc; hi++) {
                size_t got = capture_hop(st, md, centers[hi], n_samps, cfg.settle,
                                         iq, pkt, max_samps, &total_ovf);
                if (got < (size_t)nfft) continue;
                size_t sat = 0;
                for (size_t i = 0; i < got; i++)
                    if (cabsf(iq[i]) >= 0.99f) sat++;
                if ((double)sat / got > 0.001)
                    printf("  [警告] %.1f MHz 处 %.1f%% 采样饱和, 建议降低 --gain\n",
                           centers[hi] / 1e6, 100.0 * sat / got);
                welch_psd(iq, got, nfft, win, w2, plan, fin, fout, ps);
                for (int k = 0; k < nfft; k++) {
                    double off = k < nfft / 2 ? k * bin_hz : (k - nfft) * bin_hz;
                    if (fabs(off) > usable / 2 || fabs(off) < notch) continue;
                    long idx2 = llrint((centers[hi] + off - g.f0) / bin_hz);
                    if (idx2 < 0 || (size_t)idx2 >= g.nbins) continue;
                    g.sum[idx2] += ps[k];
                    g.cnt[idx2]++;
                }
            }
            for (size_t i = 0; i < g.nbins; i++) {
                if (!g.cnt[i]) continue;
                covered[i] = true;
                double db = 10 * log10(g.sum[i] / g.cnt[i] / bin_hz + 1e-20);
                if (db > g.comb_db[i]) g.comb_db[i] = db;
            }
            printf("  第 %d/%d 遍完成\n", sw + 1, cfg.sweeps);
            if (sw < cfg.sweeps - 1 && cfg.pause > 0) ms_sleep(cfg.pause);
        }
        if (total_ovf)
            printf("  [警告] 共 %d 次 USB 溢出, 建议降低 --rate\n", total_ovf);

        /* 提取有效区间 */
        size_t n = 0;
        for (size_t i = 0; i < g.nbins; i++) n += covered[i];
        if (n == 0) { printf("%s: 无有效数据\n", bands[bi].name); goto next_band; }
        double *fr = malloc(n * sizeof(double));
        double *db = malloc(n * sizeof(double));
        size_t j = 0;
        for (size_t i = 0; i < g.nbins; i++)
            if (covered[i]) { fr[j] = g.f0 + i * bin_hz; db[j] = g.comb_db[i]; j++; }

        /* ---- 占用分析 ---- */
        double floor_db = percentile(db, n, 0.2);
        double thr = floor_db + cfg.threshold;
        double *sm = malloc(n * sizeof(double));
        bool *occ = malloc(n * sizeof(bool));
        size_t n_occ = 0;
        for (size_t i = 0; i < n; i++) {
            double a = 0; int m = 0;
            for (int k = -2; k <= 2; k++) {
                long ii = (long)i + k;
                if (ii >= 0 && (size_t)ii < n) { a += db[ii]; m++; }
            }
            sm[i] = a / m;
            occ[i] = sm[i] > thr;
            n_occ += occ[i];
        }
        printf("\n===== %s: %.0f-%.0f MHz =====\n", bands[bi].name, fr[0] / 1e6, fr[n - 1] / 1e6);
        printf("底噪(20%%分位) %.1f dBFS/Hz, 占用门限 %.1f dBFS/Hz, 全段占用率 %.1f%%\n",
               floor_db, thr, 100.0 * n_occ / n);
        for (size_t ai = 0; ai < sizeof(ALLOC) / sizeof(ALLOC[0]); ai++) {
            size_t m = 0, o = 0;
            for (size_t i = 0; i < n; i++)
                if (fr[i] >= ALLOC[ai].start && fr[i] < ALLOC[ai].stop) { m++; o += occ[i]; }
            if (m)
                printf("  %-14s %7.0f-%-7.0f MHz 占用率 %5.1f%%\n", ALLOC[ai].label,
                       ALLOC[ai].start / 1e6, ALLOC[ai].stop / 1e6, 100.0 * o / m);
        }

        /* ---- 已占用频段: occ[] 连续段, 间隔 <=2 个空闲 bin 的缝隙合并 ---- */
        seg_t *segs = malloc((n / 2 + 1) * sizeof(seg_t));
        size_t nseg = 0;
        for (size_t i = 0; i < n; ) {
            if (!occ[i]) { i++; continue; }
            size_t a = i, b = i;
            while (b + 1 < n) {
                size_t nxt = b + 1;
                while (nxt < n && !occ[nxt]) nxt++;
                if (nxt < n && nxt - b <= 3) b = nxt;   /* 缝隙 <=2 bin */
                else break;
            }
            segs[nseg++] = (seg_t){ a, b };
            i = b + 1;
        }
        printf("\n-- 已占用频段 (%zu 段, 缝隙<=2 bin 已合并) --\n", nseg);
        if (!nseg) {
            printf("  (未检测到高于门限的占用)\n");
        } else {
            size_t show = nseg < 64 ? nseg : 64;
            for (size_t si = 0; si < show; si++) {
                size_t a = segs[si].a, b = segs[si].b, pk = a;
                for (size_t k = a + 1; k <= b; k++)
                    if (db[k] > db[pk]) pk = k;
                printf("  %2zu. %9.2f-%9.2f MHz (%6.2f MHz)  峰值 %6.1f dBFS/Hz @ %.2f MHz  [%s]\n",
                       si + 1, fr[a] / 1e6, fr[b] / 1e6, (fr[b] - fr[a]) / 1e6,
                       db[pk], fr[pk] / 1e6, alloc_label(fr[pk]));
            }
            if (nseg > show)
                printf("  ... 其余 %zu 段省略, 详见 CSV occupied 列\n", nseg - show);
        }
        free(segs);

        /* ---- 空闲候选 ---- */
        cand_t *picked_all = malloc((size_t)cfg.n_chan_bw * cfg.top * sizeof(cand_t));
        double *picked_bw = malloc((size_t)cfg.n_chan_bw * cfg.top * sizeof(double));
        int *picked_n = calloc((size_t)cfg.n_chan_bw, sizeof(int));
        for (int cbi = 0; cbi < cfg.n_chan_bw; cbi++) {
            double bw = cfg.chan_bw[cbi];
            size_t win = (size_t)llround(bw / bin_hz);
            printf("\n-- %.0f MHz 带宽最空闲的候选 --\n", bw / 1e6);
            if (win >= n) { printf("  (带宽超过扫描范围)\n"); continue; }
            size_t stride = (size_t)llround(1e6 / bin_hz);
            if (stride < 1) stride = 1;
            cand_t *cands = malloc((n / stride + 2) * sizeof(cand_t));
            size_t ncd = 0;
            for (size_t i = 0; i + win <= n; i += stride) {
                size_t o = 0;
                double p = 0;
                for (size_t k = i; k < i + win; k++) {
                    o += occ[k];
                    p += pow(10.0, db[k] / 10.0);
                }
                if (1.0 - (double)o / win < cfg.free_min) continue;
                cands[ncd].ch_db = 10 * log10(p * bin_hz + 1e-20);
                cands[ncd].center = fr[i + win / 2];
                ncd++;
            }
            qsort(cands, ncd, sizeof(cand_t), cmp_cand);
            int cnt_pick = 0;
            for (size_t i = 0; i < ncd && cnt_pick < cfg.top; i++) {
                long narf = freq_to_arfcn(cands[i].center);
                double cf = arfcn_to_freq(narf);
                if (cf - bw / 2 < fr[0] || cf + bw / 2 > fr[n - 1]) continue;
                bool overlap = false;
                for (int k = 0; k < cnt_pick; k++)
                    if (fabs(cf - picked_all[cbi * cfg.top + k].center) < bw) overlap = true;
                if (overlap) continue;
                cand_t pk = { cands[i].ch_db, cf, narf };
                picked_all[cbi * cfg.top + cnt_pick] = pk;
                picked_bw[cbi * cfg.top + cnt_pick] = bw;
                picked_n[cbi]++;
                double ss; long gscn;
                nearest_ssb(cf, &ss, &gscn);
                printf("  %d. %9.2f MHz  NR-ARFCN %-7ld 信道功率 %6.1f dBFS  [%s]  "
                       "最近SSB %.2f MHz (GSCN %ld, ARFCN %ld)\n",
                       cnt_pick + 1, cf / 1e6, narf, pk.ch_db, alloc_label(cf),
                       ss / 1e6, gscn, freq_to_arfcn(ss));
                cnt_pick++;
            }
            if (!cnt_pick) printf("  (无满足条件的空闲信道)\n");
            free(cands);
        }

        /* ---- CSV ---- */
        char path[512];
        snprintf(path, sizeof(path), "%s_%s.csv", cfg.out, bands[bi].name);
        FILE *fp = fopen(path, "w");
        if (fp) {
            fprintf(fp, "freq_hz,power_dbfs_per_hz,occupied\n");
            for (size_t i = 0; i < n; i++)
                fprintf(fp, "%.1f,%.3f,%d\n", fr[i], db[i], occ[i] ? 1 : 0);
            fclose(fp);
            printf("已保存: %s\n", path);
        }

        /* ---- gnuplot 脚本 (装有 gnuplot 时自动生成 PNG) ---- */
        snprintf(path, sizeof(path), "%s_%s.gp", cfg.out, bands[bi].name);
        fp = fopen(path, "w");
        if (fp) {
            double p1 = percentile(db, n, 0.01) - 8;
            double pmax = db[0];
            for (size_t i = 1; i < n; i++) if (db[i] > pmax) pmax = db[i];
            char png[512], csv[512];
            snprintf(png, sizeof(png), "%s_%s.png", cfg.out, bands[bi].name);
            snprintf(csv, sizeof(csv), "%s_%s.csv", cfg.out, bands[bi].name);
            char now[64];
            time_t t = time(NULL);
            struct tm tmv;
            localtime_r(&t, &tmv);
            strftime(now, sizeof(now), "%Y-%m-%d %H:%M", &tmv);
            fprintf(fp,
                    "set terminal pngcairo size 1400,500\nset output '%s'\n"
                    "set title '%s  gain=%.0fdB rate=%.2fMS/s sweeps=%d  %s'\n"
                    "set xlabel 'MHz'\nset ylabel 'dBFS/Hz (max-hold)'\n"
                    "set datafile separator ','\n"
                    "set grid\nset yrange [%.1f:%.1f]\n",
                    png, bands[bi].name, cfg.gain, cfg.rate / 1e6, cfg.sweeps, now,
                    p1, pmax + 12);
            int obj = 0;
            for (size_t i = 0; i < n && obj < 400; i++) {
                if (!occ[i]) continue;
                size_t k = i;
                while (k + 1 < n && occ[k + 1]) k++;
                fprintf(fp, "set object %d rect from %.3f,graph 0 to %.3f,graph 1 "
                            "fc rgb 'red' fs transparent solid 0.12 noborder\n",
                        ++obj, fr[i] / 1e6, fr[k] / 1e6);
                i = k;
            }
            const char *cc[] = {"green", "orange", "purple", "brown"};
            for (int cbi = 0; cbi < cfg.n_chan_bw; cbi++)
                for (int k = 0; k < picked_n[cbi] && obj < 460; k++) {
                    cand_t *pk = &picked_all[cbi * cfg.top + k];
                    double bw = picked_bw[cbi * cfg.top + k];
                    fprintf(fp, "set object %d rect from %.3f,graph 0 to %.3f,graph 1 "
                                "fc rgb '%s' fs transparent solid 0.2 noborder\n",
                            ++obj, (pk->center - bw / 2) / 1e6, (pk->center + bw / 2) / 1e6,
                            cc[cbi % 4]);
                    fprintf(fp, "set label %d '%.1f\\nARFCN %ld' at %.3f,graph 0.97 center font ',7' "
                                "tc rgb '%s' front\n",
                            obj, pk->center / 1e6, pk->arfcn, pk->center / 1e6, cc[cbi % 4]);
                }
            fprintf(fp, "plot '%s' skip 1 using ($1/1e6):2 with lines lw 0.5 lc rgb 'royalblue' notitle, "
                        "%.2f with lines dt 2 lc rgb 'red' title 'threshold'\n", csv, thr);
            fclose(fp);
            printf("已保存: %s\n", path);
            if (system("command -v gnuplot >/dev/null 2>&1") == 0) {
                char cmdline[600];
                snprintf(cmdline, sizeof(cmdline), "gnuplot '%s'", path);
                if (system(cmdline) == 0)
                    printf("已保存: %s\n", png);
            } else {
                printf("(未检测到 gnuplot, 未生成 PNG; 安装后可手动执行: gnuplot %s)\n", path);
            }
        }

        free(fr); free(db); free(sm); free(occ);
        free(picked_all); free(picked_bw); free(picked_n);
next_band:
        free(g.sum); free(g.cnt); free(g.comb_db); free(covered);
    }

    printf("\n提示: 结果为多次扫描 max-hold; TDD 小区可能间歇发射, "
           "重要频点建议用 --sweeps 10 --pause 5 复测确认长期空闲。\n");

    free(win); free(ps); free(iq); free(pkt);
    fftwf_destroy_plan(plan);
    fftwf_free(fin); fftwf_free(fout);
    uhd_rx_metadata_free(&md);
    uhd_rx_streamer_free(&st);
    uhd_usrp_free(&g_usrp);
    return 0;
}
