# USRP NR 频谱扫描工具使用说明

基于 USRP B210 的 NR 频段频谱扫描工具，用于观测周围环境的频谱占用情况，在 n41 / n78（或任意自定义频段）中寻找空闲频点，为自建 ocudu (srsRAN) + Open5GS 实验网选址频点。提供算法完全一致的 Python / C 两种实现。

## 项目文件

| 文件 | 说明 |
|---|---|
| `nr_scan.py` | Python 版扫描工具（UHD Python API + numpy + matplotlib） |
| `nr_scan.c` | C 版扫描工具（UHD C API + FFTW），独立二进制 |
| `Makefile` | C 版构建脚本，`make` 编译、`make clean` 清理 |
| `README.md` | 使用说明 + 开发进度记录 |

## 1. 环境依赖

| 依赖 | 说明 |
|---|---|
| USRP B210 | USB3 连接（USB2 可用但需降低采样率，见 4.2） |
| libuhd + python3-uhd | 本机已装 UHD 4.6 |
| Python 3 | Python 版需 numpy、matplotlib（本机已装） |
| gcc、libuhd-dev、libfftw3-dev、gnuplot | C 版编译与出图（本机已装） |

运行前先确认设备识别正常：

```bash
uhd_find_devices     # 应能看到 product: B210
```

## 2. 快速开始

```bash
cd ~/usrp-scanner

# === Python 版 ===
python3 nr_scan.py                            # 默认扫描: n41 + n78, 每频段 3 遍 max-hold
python3 nr_scan.py --preset n78               # 只扫 n78
python3 nr_scan.py --start 3400e6 --stop 3600e6   # 自定义频率范围

# === C 版（参数与 Python 版完全相同） ===
make                                          # 编译, 生成 ./nr_scan
./nr_scan
./nr_scan --preset n78
```

扫描过程中终端会打印进度；每个频段扫描结束后输出占用统计、空闲候选频点表，并生成 PNG 图和 CSV 数据文件。

## 3. 输出说明

### 3.1 终端输出

```
===== n78: 3300-3800 MHz =====
底噪(20%分位) -134.1 dBFS/Hz, 占用门限 -124.1 dBFS/Hz, 全段占用率 0.0%
  电信 n78            3400-3500    MHz 占用率   0.0%
  ...

-- 已占用频段 (1 段, 缝隙<=2 bin 已合并) --
   1.   3411.96-  3411.96 MHz (  0.00 MHz)  峰值 -109.5 dBFS/Hz @ 3411.96 MHz  [电信 n78]

-- 20 MHz 带宽最空闲的前 5 个候选 --
  1.   3452.61 MHz  NR-ARFCN 630174  信道功率  -61.1 dBFS  [电信 n78]  最近SSB 3452.16 MHz (GSCN 7813, ARFCN 630144)
```

已占用频段每行含义：

- **频率范围 / 宽度**：连续高于门限的频段（缝隙 ≤2 bin 的相邻段已合并为一段）
- **峰值**：段内最大功率及所在频点；**[所属分配]**：峰值频点落在哪个运营商的牌照频段内
- 段数超过 64 时只显示前 64 段，完整占用情况见 CSV 的 `occupied` 列；全部低于门限时显示"未检测到高于门限的占用"

候选频点每行含义：

- **中心频率 / NR-ARFCN**：已对齐到 3GPP 合法频点，可直接用于 gNB 配置
- **信道功率**：该带宽内积分功率（dBFS），越小越安静
- **[所属分配]**：该频点落在哪个运营商的牌照频段内
- **最近 SSB**：最近的 SSB 栅格点及其 GSCN / ARFCN（ocudu 配置 SSB 频点时用）

> 注意：功率单位是 dBFS（相对 ADC 满量程），未做绝对功率校准，只适合相对比较"哪个频点更安静"。

### 3.2 PNG 图（`scan_时间戳_<band>.png`）

- 蓝线：全频段功率谱（多遍扫描 max-hold）
- 红色阴影：判定为"已占用"的区域；红色虚线：占用判决门限
- 绿/橙/紫色阴影：不同带宽的空闲候选信道，标注了中心频率和 NR-ARFCN

### 3.3 CSV（`scan_时间戳_<band>.csv`）

三列：`freq_hz, power_dbfs_per_hz, occupied`（occupied 为 0/1 占用判决结果），可用 pandas / Excel 导入做进一步分析。

## 4. 参数说明

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--preset` | `n41,n78` | 预设频段，逗号分隔，可选 `n41` `n78` `n79` |
| `--start` / `--stop` | 无 | 自定义起止频率 (Hz)，设置后覆盖 `--preset` |
| `--gain` | 40 | RX 增益 dB（B210 范围 0–76）。信号弱调高，出现饱和警告调低 |
| `--rate` | 23.04e6 | 采样率 S/s，决定每个调谐点覆盖的宽度 |
| `--fft` | 2048 | FFT 点数，频率分辨率 = rate/fft ≈ 11.25 kHz |
| `--dwell` | 0.08 | 每个调谐点驻留采集时间（秒），越大平均越充分、底噪越平 |
| `--settle` | 0.02 | 每次调谐后的稳定等待时间（秒） |
| `--sweeps` | 3 | 每个频段扫描遍数（多遍按频点取 max-hold） |
| `--pause` | 1.0 | 两遍扫描之间的间隔（秒） |
| `--threshold` | 10 | 占用判决门限：底噪 + N dB |
| `--chan-bw` | `10e6,20e6,40e6` | 候选信道带宽列表 (Hz) |
| `--free-min` | 0.95 | 候选信道内允许的最大占用比例（0.95 = 95% 的频点必须空闲） |
| `--top` | 5 | 每种带宽输出前 N 个最空闲候选 |
| `--antenna` | RX2 | 接收天线口 |
| `--device` | 空 | UHD 设备参数，多设备时用如 `serial=3129082` |
| `--out` | `scan_时间戳` | 输出文件前缀 |

### 4.1 关键概念：为什么要多遍扫描（max-hold）

n41 / n78 都是 TDD 频段，基站只在下行时隙发射。单遍扫描可能恰好赶上空闲时隙，误判为"无人使用"。`--sweeps N --pause T` 会在 N 遍扫描间拉开时间间隔，并按频点取功率最大值（max-hold），只要基站在任何一遍的任一时刻发过射，就会被记录下来。**部署前对选定频点务必加长复测**：

```bash
python3 nr_scan.py --preset n41 --sweeps 10 --pause 5 --top 3
```

### 4.2 常见警告处理

| 警告 | 原因 | 处理 |
|---|---|---|
| `USB 溢出` | 采样率超过 USB 实际带宽（多见于 USB2 口） | 降低 `--rate`，如 `--rate 15.36e6` 或 USB2 下 `8e6` |
| `采样饱和` | 增益过高，ADC 削波 | 降低 `--gain` 5–10 dB 重扫 |
| 频谱全空、只有等间隔细尖峰 | 尖峰是电脑/USB3 的本地杂散，不是基站 | 属正常现象；可换位置/天线对比 |

### 4.3 提高灵敏度

环境信号弱（室内、天线差）时：

```bash
python3 nr_scan.py --gain 60 --dwell 0.2 --sweeps 5
```

提高增益 + 加长驻留可以把底噪压得更平、更容易看到弱信号。也建议换到窗边或使用增益更好的天线对比扫描。

## 5. 与 ocudu 配置对接

选定候选频点后，gNB 配置（如 `gnb_ru_ran_config.yml` 风格）需要的关键参数都能从候选行直接拿到：

- **band**：n41 或 n78
- **dl_arfcn**：候选行的 NR-ARFCN
- **SSB ARFCN**：候选行"最近 SSB"括号里的 ARFCN（SSB 必须落在 GSCN 栅格上）
- **channel_bandwidth**：候选对应的带宽（10/20/40 MHz）
- **subcarrier_spacing**：n41/n78 用 30 kHz
- **TDD pattern**：国内 n41/n78 宏网常用 DDDSU（2.5ms 双周期）类配比，自建网可按需配置

## 6. C 版本（nr_scan.c）

与 Python 版算法完全一致（同一套扫频 / max-hold / 占用判决 / 候选搜索逻辑），基于 UHD C API + FFTW 实现，编译后是不依赖 Python 的独立二进制：

```bash
cd ~/usrp-scanner
make                                  # 依赖 libuhd-dev、libfftw3-dev（本机已装）
./nr_scan --preset n41,n78            # 参数与 Python 版完全相同
./nr_scan --help                      # 查看全部参数
```

输出格式与 Python 版一致（终端报告 + CSV + PNG）。PNG 通过自动调用 gnuplot 生成（本机已装 gnuplot 6.0），同时保留 `.gp` 脚本可手动重绘；系统没有 gnuplot 时只出 CSV 并在终端提示。两版结果已互相对照验证（底噪、占用率、候选频点一致）。

选型建议：改算法、做数据分析用 Python 版方便；要部署到无 Python 环境、或集成进其他 C 项目用 C 版。

## 7. 法律与安全提醒

- 本工具只做**被动接收**，不发射任何信号，扫描本身安全。
- n41 / n78 / n79 在国内均为运营商 licensed 频段。后续用 ocudu **起站发射**时，务必：使用最小发射功率、串接衰减器、尽量在屏蔽箱/法拉第环境内操作，避免干扰现网（违法且会被监测定位）。
- 扫描结果反映的是**测量时刻、测量位置**的频谱状态；换地点、换时间应重新扫描。


### 2026-10-18

**目标**：基于 B210 + libuhd 开发频谱扫描程序，在 n41/n78 中寻找空闲频点，为 ocudu + Open5GS 自建网选址。

已完成：

- [x] 环境确认：UHD 4.6 + python3-uhd + numpy/matplotlib + libuhd-dev + libfftw3-dev，B210 (serial 3129082) 识别正常
- [x] Python 版 `nr_scan.py`：分段扫频（段间重叠拼接、LO 直流缺口剔除）、Welch PSD、多遍 max-hold（抓 TDD 突发）、底噪+门限占用判决、滑窗空闲信道搜索、候选对齐 NR-ARFCN、最近 SSB 栅格计算、PNG+CSV 输出
- [x] 真机首次实测（n41+n78 各 3 遍）：本环境两频段占用率 ≈0%，只见本地杂散尖峰；结果存于 `scan_20261018_223327_*.png/csv`
- [x] 编写使用说明（本文档）
- [x] C 版 `nr_scan.c` + `Makefile`：与 Python 版算法一致，编译无警告；修复 gnuplot 不识别逗号分隔 CSV 导致空图的问题（加 `set datafile separator ','`）
- [x] 两版对照验证：底噪（-133.9 vs -134.0 dBFS/Hz）、占用率、候选频点一致

待办（TODO）：

- [ ] 选定候选频点后做长时间复测确认（`--sweeps 10 --pause 5`），必要时换位置/天线对比
- [ ] 绝对功率校准：dBFS → dBm（需已知功率的信号源做标定）
- [ ] SSB 突发检测/解调解码：直接识别现网小区（PCI、SSB 频点），比纯能量检测更可靠
- [ ] 连续监测模式：定时循环扫描、结果追加日志，观察频谱占用随时间变化
- [ ] 扫描结果与 ocudu gNB 配置联动：自动生成/校验 gnb 配置中的 band、dl_arfcn、SSB 频点

> 记录格式约定：每次开发在对应日期下追加"已完成 / 待办"事项；待办完成后移入已完成并注明日期。
