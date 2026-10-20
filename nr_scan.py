#!/usr/bin/env python3
"""NR 频段频谱扫描工具 (USRP B210 / UHD Python API)

扫过 n41/n78(或自定义频段), 输出:
  - 终端: 底噪、占用率、空闲候选频点(含 NR-ARFCN、最近 SSB 栅格)
  - PNG : 全频段功率谱 + 占用区域/候选频点标注
  - CSV : 频点-功率数据

多次扫描默认按 bin 取 max-hold, 用于捕捉 TDD 小区的突发下行。
"""

import argparse
import datetime
import math
import sys
import time

import numpy as np
import uhd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 预设频段 (Hz)
PRESETS = {
    "n41": (2515e6, 2675e6),
    "n78": (3300e6, 3800e6),
    "n79": (4400e6, 5000e6),
}

# 国内 5G 频谱分配标注
ALLOC = [
    (2515e6, 2675e6, "移动 n41"),
    (3300e6, 3400e6, "3.3G 室内共享"),
    (3400e6, 3500e6, "电信 n78"),
    (3500e6, 3600e6, "联通 n78"),
    (4800e6, 4900e6, "移动 n79"),
    (4900e6, 4960e6, "广电 n79"),
]


def freq_to_arfcn(f_hz):
    """3GPP 38.104 全局频率栅格 -> NR-ARFCN"""
    if f_hz < 3000e6:
        return int(round(f_hz / 5e3))
    return 600000 + int(round((f_hz - 3000e6) / 15e3))


def arfcn_to_freq(n):
    if n < 600000:
        return n * 5e3
    return 3000e6 + (n - 600000) * 15e3


def nearest_ssb(f_hz):
    """最近 SSB 栅格点, 返回 (SS_REF Hz, GSCN)"""
    if f_hz < 3000e6:
        best = None
        for m in (1, 3, 5):
            n = int(round((f_hz - m * 50e3) / 1.2e6))
            n = min(max(n, 1), 2499)
            ss = n * 1.2e6 + m * 50e3
            gscn = 3 * n + (m - 3) // 2
            d = abs(ss - f_hz)
            if best is None or d < best[0]:
                best = (d, ss, gscn)
        return best[1], best[2]
    n = int(round((f_hz - 3000e6) / 1.44e6))
    n = min(max(n, 0), 14756)
    return 3000e6 + n * 1.44e6, 7499 + n


def alloc_label(f_hz):
    for a, b, label in ALLOC:
        if a <= f_hz < b:
            return label
    return "未标注"


def parse_floats(s):
    return [float(x) for x in s.split(",") if x.strip()]


def parse_args():
    p = argparse.ArgumentParser(
        description="USRP B210 NR 频段频谱扫描, 寻找空闲频点",
        epilog="示例: python3 nr_scan.py --preset n41,n78 --sweeps 5 --pause 2\n"
               "      python3 nr_scan.py --start 3400e6 --stop 3600e6 --gain 30",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--preset", default="n41,n78",
                   help="预设频段, 逗号分隔, 可选: " + ",".join(PRESETS))
    p.add_argument("--start", type=float, default=None, help="自定义起始频率 Hz (覆盖 preset)")
    p.add_argument("--stop", type=float, default=None, help="自定义结束频率 Hz")
    p.add_argument("--device", default="", help="UHD 设备参数, 如 'type=b200'")
    p.add_argument("--antenna", default="RX2", help="接收天线口 (默认 RX2)")
    p.add_argument("--gain", type=float, default=40.0, help="RX 增益 dB (默认 40)")
    p.add_argument("--rate", type=float, default=23.04e6,
                   help="采样率 S/s (默认 23.04e6, USB2 下需降到 ~8e6)")
    p.add_argument("--fft", type=int, default=2048, help="FFT 点数 (默认 2048)")
    p.add_argument("--dwell", type=float, default=0.08, help="每个调谐点驻留秒数 (默认 0.08)")
    p.add_argument("--settle", type=float, default=0.02, help="调谐后稳定等待秒数 (默认 0.02)")
    p.add_argument("--sweeps", type=int, default=3, help="每个频段扫描遍数 (默认 3)")
    p.add_argument("--pause", type=float, default=1.0, help="两遍扫描之间的间隔秒数 (默认 1.0)")
    p.add_argument("--threshold", type=float, default=10.0,
                   help="占用判决门限: 底噪+N dB (默认 10)")
    p.add_argument("--chan-bw", default="10e6,20e6,40e6",
                   help="候选信道带宽列表 Hz (默认 10e6,20e6,40e6)")
    p.add_argument("--free-min", type=float, default=0.95,
                   help="候选信道内空闲 bin 比例下限 (默认 0.95)")
    p.add_argument("--top", type=int, default=5, help="每种带宽输出前 N 个候选 (默认 5)")
    p.add_argument("--out", default=None, help="输出文件前缀 (默认 scan_时间戳)")
    return p.parse_args()


def capture_hop(usrp, streamer, freq, n_samps, settle):
    """调谐到 freq, 采集 n_samps 个采样, 返回 (iq, 溢出次数)"""
    usrp.set_rx_freq(uhd.types.TuneRequest(freq), 0)
    if settle > 0:
        time.sleep(settle)
    max_samps = streamer.get_max_num_samps()
    buf = np.empty((1, max_samps), dtype=np.complex64)
    out = np.empty(n_samps, dtype=np.complex64)
    cmd = uhd.types.StreamCMD(uhd.types.StreamMode.num_done)
    cmd.num_samps = int(n_samps)
    cmd.stream_now = True
    md = uhd.types.RXMetadata()
    streamer.issue_stream_cmd(cmd)
    got, ovf = 0, 0
    while got < n_samps:
        n = streamer.recv(buf, md, 1.0)
        err = md.error_code
        if err == uhd.types.RXMetadataErrorCode.overflow:
            ovf += 1
            continue
        if err == uhd.types.RXMetadataErrorCode.timeout:
            break
        if err != uhd.types.RXMetadataErrorCode.none:
            continue
        if n > 0:
            m = min(n, n_samps - got)
            out[got:got + m] = buf[0, :m]
            got += m
    return out[:got], ovf


def welch_psd(iq, nfft):
    """Welch 平均周期图, 返回 (偏移 Hz 数组, 每 bin 功率, 满足 sum(功率)=信号平均功率)"""
    hop = nfft // 2
    win = np.hanning(nfft)
    w2 = float(np.mean(win ** 2))
    nseg = 1 + (len(iq) - nfft) // hop
    acc = np.zeros(nfft)
    for i in range(nseg):
        seg = iq[i * hop:i * hop + nfft] * win
        acc += np.abs(np.fft.fft(seg)) ** 2
    acc /= (nseg * nfft ** 2 * w2)
    return np.fft.fftshift(acc)


def hop_centers(f_start, f_stop, usable):
    if f_stop - f_start <= usable:
        return [(f_start + f_stop) / 2]
    step = 0.75 * usable
    centers = []
    c = f_start + usable / 2
    while c + usable / 2 < f_stop:
        centers.append(c)
        c += step
    centers.append(f_stop - usable / 2)
    return centers


def sweep_band(usrp, streamer, args, f_start, f_stop):
    """对 [f_start, f_stop] 扫 args.sweeps 遍, 按 bin max-hold, 返回 (freqs, psd_dbfs_hz, bin_hz)"""
    rate, nfft = args.rate, args.fft
    bin_hz = rate / nfft
    usable = 0.7 * rate
    notch = min(0.6e6, 0.05 * rate)  # 去掉每个 hop 中心 LO 直流尖峰附近
    n_samps = max(nfft * 4, int(rate * args.dwell))
    centers = hop_centers(f_start, f_stop, usable)

    f0 = math.floor(f_start / bin_hz) * bin_hz
    nbins = int(math.ceil((f_stop - f0) / bin_hz)) + 1
    grid_f = f0 + np.arange(nbins) * bin_hz

    est = len(centers) * (args.dwell + args.settle + 0.05) * args.sweeps
    print(f"  频段 {f_start/1e6:.0f}-{f_stop/1e6:.0f} MHz: "
          f"{len(centers)} 个调谐点 x {args.sweeps} 遍, 预计 {est:.0f} 秒")

    combined = np.full(nbins, -np.inf)
    total_ovf = 0
    for s in range(args.sweeps):
        sum_lin = np.zeros(nbins)
        cnt = np.zeros(nbins)
        for c in centers:
            iq, ovf = capture_hop(usrp, streamer, c, n_samps, args.settle)
            total_ovf += ovf
            if len(iq) < nfft:
                continue
            sat = float(np.mean(np.abs(iq) >= 0.99))
            if sat > 0.001:
                print(f"  [警告] {c/1e6:.1f} MHz 处 {sat*100:.1f}% 采样饱和, 建议降低 --gain")
            ps = welch_psd(iq, nfft)
            offs = (np.arange(nfft) - nfft // 2) * bin_hz
            keep = (np.abs(offs) <= usable / 2) & (np.abs(offs) >= notch)
            idx = np.rint((c + offs[keep] - f0) / bin_hz).astype(int)
            np.add.at(sum_lin, idx, ps[keep])
            np.add.at(cnt, idx, 1)
        valid = cnt > 0
        sweep_db = np.full(nbins, -np.inf)
        sweep_db[valid] = 10 * np.log10(sum_lin[valid] / cnt[valid] / bin_hz + 1e-20)
        combined = np.maximum(combined, sweep_db)
        print(f"  第 {s + 1}/{args.sweeps} 遍完成")
        if s < args.sweeps - 1 and args.pause > 0:
            time.sleep(args.pause)
    if total_ovf:
        print(f"  [警告] 共 {total_ovf} 次 USB 溢出, 频谱可能有缺口, 建议降低 --rate")
    ok = np.isfinite(combined)
    return grid_f[ok], combined[ok], bin_hz


def find_candidates(freqs, psd_db, occupied, bin_hz, bw_hz, args):
    win = int(round(bw_hz / bin_hz))
    if win >= len(freqs):
        return []
    stride = max(1, int(round(1e6 / bin_hz)))
    cands = []
    for i in range(0, len(freqs) - win, stride):
        if 1.0 - occupied[i:i + win].mean() < args.free_min:
            continue
        ch_db = 10 * np.log10(np.sum(10 ** (psd_db[i:i + win] / 10.0)) * bin_hz + 1e-20)
        cands.append((ch_db, freqs[i + win // 2]))
    cands.sort()
    picked = []
    for ch_db, center in cands:
        n = freq_to_arfcn(center)
        cf = arfcn_to_freq(n)  # 对齐到合法 NR-ARFCN 频点
        if cf - bw_hz / 2 < freqs[0] or cf + bw_hz / 2 > freqs[-1]:
            continue
        if any(abs(cf - p[1]) < bw_hz for p in picked):
            continue
        picked.append((ch_db, cf, n))
        if len(picked) >= args.top:
            break
    return picked


def occupied_runs(mask):
    idx = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(np.int8)))
    return idx.reshape(-1, 2)


def analyze(name, freqs, psd_db, bin_hz, args):
    floor = float(np.percentile(psd_db, 20))
    thr = floor + args.threshold
    sm = np.convolve(psd_db, np.ones(5) / 5, mode="same")
    occupied = sm > thr
    occ_pct = 100.0 * occupied.mean()
    print(f"\n===== {name}: {freqs[0]/1e6:.0f}-{freqs[-1]/1e6:.0f} MHz =====")
    print(f"底噪(20%分位) {floor:.1f} dBFS/Hz, 占用门限 {thr:.1f} dBFS/Hz, "
          f"全段占用率 {occ_pct:.1f}%")
    for a, b, label in ALLOC:
        m = (freqs >= a) & (freqs < b)
        if m.any():
            print(f"  {label:<14} {a/1e6:>7.0f}-{b/1e6:<7.0f} MHz 占用率 {100*occupied[m].mean():5.1f}%")

    # 已占用频段: 连续占用段, 间隔 <=2 个空闲 bin 的缝隙合并
    segs = []
    i, nb = 0, len(freqs)
    while i < nb:
        if not occupied[i]:
            i += 1
            continue
        a = b = i
        while b + 1 < nb:
            nxt = b + 1
            while nxt < nb and not occupied[nxt]:
                nxt += 1
            if nxt < nb and nxt - b <= 3:
                b = nxt
            else:
                break
        segs.append((a, b))
        i = b + 1
    print(f"\n-- 已占用频段 ({len(segs)} 段, 缝隙<=2 bin 已合并) --")
    if not segs:
        print("  (未检测到高于门限的占用)")
    for j, (a, b) in enumerate(segs[:64], 1):
        pk = a + int(np.argmax(psd_db[a:b + 1]))
        print(f"  {j:2d}. {freqs[a]/1e6:9.2f}-{freqs[b]/1e6:9.2f} MHz "
              f"({(freqs[b]-freqs[a])/1e6:6.2f} MHz)  峰值 {psd_db[pk]:6.1f} dBFS/Hz "
              f"@ {freqs[pk]/1e6:.2f} MHz  [{alloc_label(freqs[pk])}]")
    if len(segs) > 64:
        print(f"  ... 其余 {len(segs)-64} 段省略, 详见 CSV occupied 列")

    all_picked = []
    for bw in args.chan_bw:
        picked = find_candidates(freqs, psd_db, occupied, bin_hz, bw, args)
        all_picked.append((bw, picked))
        print(f"\n-- {bw/1e6:.0f} MHz 带宽最空闲的前 {len(picked)} 个候选 --")
        if not picked:
            print("  (无满足条件的空闲信道)")
        for i, (ch_db, cf, n) in enumerate(picked, 1):
            ss, gscn = nearest_ssb(cf)
            print(f"  {i}. {cf/1e6:9.2f} MHz  NR-ARFCN {n:<7} "
                  f"信道功率 {ch_db:6.1f} dBFS  [{alloc_label(cf)}]  "
                  f"最近SSB {ss/1e6:.2f} MHz (GSCN {gscn}, ARFCN {freq_to_arfcn(ss)})")
    return occupied, thr, all_picked


def plot_band(name, freqs, psd_db, occupied, thr, all_picked, out_png, args):
    f_mhz = freqs / 1e6
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(f_mhz, psd_db, lw=0.6, color="tab:blue")
    for a, b in occupied_runs(occupied):
        ax.axvspan(f_mhz[a], f_mhz[min(b, len(f_mhz) - 1)], color="red", alpha=0.15, lw=0)
    colors = ["tab:green", "tab:orange", "tab:purple", "tab:brown"]
    for j, (bw, picked) in enumerate(all_picked):
        for ch_db, cf, n in picked:
            ax.axvspan((cf - bw / 2) / 1e6, (cf + bw / 2) / 1e6,
                       color=colors[j % len(colors)], alpha=0.25, lw=0)
            ax.text(cf / 1e6, psd_db.max(), f"{cf/1e6:.1f}\nARFCN {n}",
                    ha="center", va="top", fontsize=7, color=colors[j % len(colors)])
    ax.axhline(thr, color="red", ls="--", lw=0.8, label=f"threshold {thr:.1f} dBFS/Hz")
    for a, b, label in ALLOC:
        if b > freqs[0] and a < freqs[-1]:
            ax.axvline(a / 1e6, color="gray", ls=":", lw=0.5)
    ax.set_xlabel("MHz")
    ax.set_ylabel("dBFS/Hz (max-hold)")
    ax.set_title(f"{name}  gain={args.gain}dB  rate={args.rate/1e6:.2f}MS/s  "
                 f"sweeps={args.sweeps}  {datetime.datetime.now():%Y-%m-%d %H:%M}")
    ax.set_ylim(np.percentile(psd_db, 1) - 8, psd_db.max() + 12)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def main():
    args = parse_args()
    args.chan_bw = parse_floats(args.chan_bw)
    if args.start is not None and args.stop is not None:
        bands = [("custom", (args.start, args.stop))]
    else:
        bands = []
        for name in args.preset.split(","):
            name = name.strip()
            if name not in PRESETS:
                sys.exit(f"未知预设 {name}, 可选: {list(PRESETS)}")
            bands.append((name, PRESETS[name]))
    if args.out is None:
        args.out = "scan_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    usrp = uhd.usrp.MultiUSRP(args.device)
    print(f"设备: {usrp.get_mboard_name(0)}  增益范围: {usrp.get_rx_gain_range(0)}")
    usrp.set_rx_rate(args.rate, 0)
    usrp.set_rx_gain(args.gain, 0)
    usrp.set_rx_antenna(args.antenna, 0)
    try:
        usrp.set_rx_bandwidth(args.rate, 0)
    except Exception:
        pass
    stream_args = uhd.usrp.StreamArgs("fc32", "sc16")
    stream_args.channels = [0]
    streamer = usrp.get_rx_stream(stream_args)

    for name, (f_start, f_stop) in bands:
        try:
            freqs, psd_db, bin_hz = sweep_band(usrp, streamer, args, f_start, f_stop)
        except KeyboardInterrupt:
            print("\n用户中断")
            break
        if len(freqs) == 0:
            print(f"{name}: 无有效数据")
            continue
        occupied, thr, all_picked = analyze(name, freqs, psd_db, bin_hz, args)
        out_png = f"{args.out}_{name}.png"
        out_csv = f"{args.out}_{name}.csv"
        plot_band(name, freqs, psd_db, occupied, thr, all_picked, out_png, args)
        np.savetxt(out_csv, np.column_stack([freqs, psd_db, occupied.astype(int)]),
                   delimiter=",", header="freq_hz,power_dbfs_per_hz,occupied",
                   comments="", fmt=["%.1f", "%.3f", "%d"])
        print(f"已保存: {out_png} , {out_csv}")
    print("\n提示: 结果为多次扫描 max-hold; TDD 小区可能间歇发射, "
          "重要频点建议用 --sweeps 10 --pause 5 复测确认长期空闲。")


if __name__ == "__main__":
    main()
