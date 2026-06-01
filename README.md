# Ungerboeck M-BCJR FTN 论文复现

本项目复现 IEEE TCOM 2018 论文 “Reduced-Complexity Equalization for Faster-Than-Nyquist Signaling: New Methods Based on Ungerboeck Observation Model” 的主要仿真链路，并扩展实现两个文献对比算法。

当前 README 是项目的唯一进度入口。旧的阶段计划文档已合并到这里，过时的“未开始”表述和重复的实现细节不再保留。

## 当前状态

### 主论文复现

| 内容 | 状态 | 产物 |
| --- | --- | --- |
| Fig. 7: BPSK, tau=0.5, `(7,5)` 卷积码, Ungerboeck M-BCJR | 已完成 | `results/fig7_tau_05/fig7_results.csv`, `fig7_ber.png` |
| Fig. 8: BPSK, tau=0.35, `(7,5)` 卷积码, Ungerboeck M-BCJR | 已完成 | `results/fig8_tau_035/fig8_results.csv`, `fig8_ber.png` |
| Fig. 10: BPSK, R=1/3 Turbo 码, Turbo equalization | 已完成 | `results/fig10/fig10_results.csv`, `fig10_ber.png` |
| Fig. 11: 16-QAM, simplified M-BCJR, `(7,5)` 卷积码 | 已完成，有已知限制 | `results/fig11/fig11_results.csv`, `fig11_ber.png` |

主论文链路已经包含：

- RRC 脉冲生成和 FTN 自相关 tap 计算
- Ungerboeck matched-filter 观测模型和有色噪声生成
- Full BCJR oracle
- Proposed M-BCJR
- Simplified M-BCJR，支持 BPSK 和 16-QAM
- `(7,5)` 卷积码 SISO BCJR
- R=1/3 非对称 Turbo 码
- BER 结果保存为 CSV/NPZ/config JSON，并生成主论文图像

### 文献对比复现

文献对比代码已经实现并生成了 CSV/NPZ 数据，但还没有形成最终分析结论。

| 内容 | 状态 | 产物 |
| --- | --- | --- |
| Paper [14] Prlja & Anderson: Forney/WMF 模型、M-BCJR、Backup M-BCJR | 已实现；uncoded 可参考，Turbo 外信息策略仍待改进 | `ftn/baselines/forney_model.py`, `prlja_mbcjr.py` |
| Paper [26] Rusek & Prlja: channel shortening、shortened BCJR | 已实现，已有对比数据 | `ftn/baselines/channel_shortening.py`, `shortened_bcjr.py` |
| Uncoded comparison, tau=0.35/0.5 | 已生成 CSV/NPZ | `results/comparison/uncoded_tau0350`, `uncoded_tau0500` |
| Turbo comparison, tau=0.35/0.5 | 已重新生成 CSV/NPZ；Paper [14] 已完成有限 LLR 修复但仍明显落后 | `results/comparison/turbo_tau0350`, `turbo_tau0500` |
| Comparison plots | 未完成 | `plot_comparison.py` 当前受 Matplotlib/NumPy ABI 问题影响 |

需要注意：旧计划里把 [14]/[26] 对比算法标成“未开始”，这是过时信息。以本 README 和当前代码/结果目录为准。

## 目录结构

```text
ftn/
  channel.py                 # FTN 信道、有色噪声、复信道
  pulse.py                   # RRC pulse 和 g-taps
  modulation.py              # BPSK、16-QAM
  equalizers/
    full_bcjr.py             # Full BCJR oracle
    ungerboeck_mbcjr.py      # 主论文 proposed M-BCJR
    simplified_mbcjr.py      # 主论文 simplified M-BCJR
  coding/
    conv.py                  # (7,5) 卷积码 SISO BCJR
    turbo.py                 # R=1/3 Turbo 码
  baselines/
    forney_model.py          # Paper [14] Forney/WMF 模型
    prlja_mbcjr.py           # Paper [14] M-BCJR variants
    channel_shortening.py    # Paper [26] 信道缩短
    shortened_bcjr.py        # Paper [26] reduced-memory BCJR
scripts/
  reproduce_fig7.py
  reproduce_fig8.py
  reproduce_fig10.py
  reproduce_fig11.py
  run_parallel.py
  plot_results.py
  run_parallel_comparison.py
  plot_comparison.py
tests/
  test_*.py
results/
  fig7_tau_05/
  fig8_tau_035/
  fig10/
  fig11/
  comparison/
docs/
  Paper.md and paper notes     # 原始论文资料和摘录，非进度入口
```

## 运行与验证

测试需要把项目根目录加入 Python 导入路径：

```bash
PYTHONPATH=. pytest -q
```

当前验证结果：`81 passed, 1 warning`。直接运行 `pytest -q` 可能因为导入路径未设置而报 `No module named 'ftn'`。

主论文快速冒烟：

```bash
PYTHONPATH=. python scripts/run_smoke.py
PYTHONPATH=. python scripts/reproduce_fig7.py --quick
PYTHONPATH=. python scripts/reproduce_fig8.py --quick
PYTHONPATH=. python scripts/reproduce_fig10.py --quick
PYTHONPATH=. python scripts/reproduce_fig11.py --quick
```

主论文完整复现：

```bash
PYTHONPATH=. python scripts/reproduce_fig7.py
PYTHONPATH=. python scripts/reproduce_fig8.py
PYTHONPATH=. python scripts/reproduce_fig10.py
PYTHONPATH=. python scripts/reproduce_fig11.py
```

已有主论文 PNG 保存在对应 `results/fig*` 目录。若需要重新出图，先修复当前环境里的 Matplotlib/NumPy ABI 兼容性，然后运行：

```bash
PYTHONPATH=. python scripts/plot_results.py both
```

文献对比实验：

```bash
PYTHONPATH=. python scripts/run_parallel_comparison.py --mode all --tau 0.5 0.35
PYTHONPATH=. python scripts/plot_comparison.py all
```

`plot_results.py` 和 `plot_comparison.py` 当前可能因环境里的 Matplotlib 与 NumPy 2.x ABI 不兼容而失败。最近日志里的错误是 `_ARRAY_API not found` / `numpy.core.multiarray failed to import`。

## 主要结果摘要

### Fig. 7

`tau=0.5` 下，`M=2,L=5` 在中高 SNR 区间整体优于 `M=2,L=3`，并在 4.5-5.5 dB 附近接近 no-ISI 基线。低 BER 区受 Monte Carlo 错误数较少影响，会有轻微波动。

### Fig. 8

`tau=0.35` 下 ISI 很强，`M=8,L=7` 在瀑布区比 `M=8,L=5` 更有效，并在约 4.5 dB 后接近 no-ISI 基线。低 SNR 区 BER 较高，符合强 ISI 下检测难度增加的预期。

### Fig. 10

Turbo-coded BPSK 链路已修正为“全部 coded bits 通过 FTN 信道，再做 detector/decoder 外信息交换”。`tau=1.0` 曲线与 no-ISI 基线基本吻合，`tau=2/3` 随 SNR 呈明显下降，`tau=0.5` 因 ISI 更强下降较慢。

### Fig. 11

16-QAM 链路和 no-ISI 基线正常，`tau=0.8` 曲线可进入低 BER 区。`tau=2/3` 在当前 `M=4,L=3` 设置下存在 error floor，主要来自 16-QAM 强 ISI 场景下保留状态数过少；这一点应作为已知限制，而不是当作论文完全复现失败。

### Comparison

`results/comparison` 下已有 tau=0.35 和 tau=0.5 的 uncoded/turbo BER 数据。最近一次完整重跑命令为：

```bash
PYTHONPATH=. python scripts/run_parallel_comparison.py --mode all --tau 0.5 0.35
```

这次重跑已写回：

- `results/comparison/uncoded_tau0500/uncoded_ber.csv`
- `results/comparison/turbo_tau0500/turbo_ber.csv`
- `results/comparison/uncoded_tau0350/uncoded_ber.csv`
- `results/comparison/turbo_tau0350/turbo_ber.csv`

Uncoded comparison 的总体趋势正常：所有算法随 SNR 增大 BER 下降。在 `tau=0.5`、6.0 dB 时，`cs_nu5=2.58e-2`、`prlja_M8=2.74e-2`、`full_bcjr=3.90e-2`、`ung_M8=7.77e-2`。在 `tau=0.35`、6.5 dB 时，`cs_nu5=1.05e-1`、`cs_nu3=1.10e-1`、`full_bcjr=1.38e-1`、`prlja_M8=1.66e-1`、`ung_M8=2.08e-1`。这些数值说明 Paper [26] shortened BCJR 在当前 comparison 设置下表现最好，Paper [14] uncoded M=8 也有可解释的下降趋势。

Turbo comparison 中，Ungerboeck M-BCJR 和 Paper [26] shortened BCJR 可进入低 BER 区。在 `tau=0.5`、6.0 dB 时，`turbo_ung_M4=5.0e-6`、`turbo_ung_M8=5.0e-6`、`turbo_cs_nu3=5.0e-6`、`turbo_cs_nu2=2.0e-5`，接近 no-ISI 基线 `1.17e-5`。在 `tau=0.35`、6.5 dB 时，`turbo_cs_nu3=1.67e-6`、`turbo_cs_nu2=5.0e-6`，优于 no-ISI 基线 `3.33e-6` 的差异属于有限 Monte Carlo 误差范围；`turbo_ung_M8=3.83e-2`、`turbo_ung_M4=8.62e-2`，强 ISI 下仍有明显残余误差。

Paper [14] turbo 曲线已从“接近随机判决”改善，但仍不能作为最终性能结论。修复前 `turbo_prlja_M4/M8` 长期在约 0.44-0.50；修复 `PrljaMbcjrResult.llr` 的非有限值后，`tau=0.5`、6.0 dB 降到 `turbo_prlja_M4=2.85e-1`、`turbo_prlja_M8=2.59e-1`，`tau=0.35`、6.5 dB 仍为 `turbo_prlja_M4=4.25e-1`、`turbo_prlja_M8=4.15e-1`。短帧诊断显示保守反馈增益有帮助：`tau=0.5`、5.0 dB、K=1000、3 个 seed 下，默认反馈 final BER 约 `0.320`，反馈增益 0.25 的 best BER 约 `0.078`，但这还不是正式曲线参数。

### Diagnostics

本轮推进修复了 Paper [14] baseline 的一个数值问题：真实 FTN minimum-phase Forney 前端下，`prlja_backup_mbcjr_bpsk` 会输出 `+/-inf`，这些 LLR 在 turbo loop 中被裁剪为饱和值后会产生过强错误先验。现在 public API 会把 `NaN` 转成 0、把 `+/-inf` 限制到 reserve LLR `+/-5`，并新增回归测试覆盖该场景。

诊断结果：

- `scripts/diag_prlja_turbo.py` 修复前：Prlja uncoded BER 约 `0.1455`，但 detector LLR 含 `+/-inf`；turbo 迭代从约 `0.136` 恶化到 `0.449`。
- `scripts/diag_prlja_turbo.py` 修复后：Prlja uncoded BER 约 `0.1565`，LLR 全部有限，turbo 首轮约 `0.114`，但默认全强度反馈后仍振荡并恶化到约 `0.362`。
- `scripts/diag_minphase.py` 显示 `tau=0.35` 下 scipy `minimum_phase(..., method="hilbert")` 会产生 NaN，当前 `min_phase_from_pulse` 已回退到 `homomorphic`，输出有限。
- `scripts/diag_tau035.py` 显示 channel shortening 的 LLR 符号方向正确，反转 LLR 会显著变差。

## 已知问题

- `plot_results.py` 和 `plot_comparison.py` 在当前环境下可能无法导入 Matplotlib，需修复 NumPy/Matplotlib ABI 兼容性。
- Paper [14] 的 turbo comparison 曲线经有限 LLR 修复后仍明显落后，下一步应系统调参外信息反馈增益、smoothing 和 backup 递推策略。
- Fig. 11 的 `tau=2/3` 在 `M=4,L=3` 下存在 error floor，可尝试更大的 `M`、更长 future length 或原始 M-BCJR 进行验证。
- `scripts/plot_results.py` 中 `MARKERS` 字典存在重复 `no_isi` key，当前已有 PNG 不受影响，但后续整理绘图脚本时应修正。

## 下一步

1. 修复绘图环境或更新 Matplotlib，使 comparison PNG 能稳定生成。
2. 为 Paper [14] turbo comparison 增加显式 feedback gain 参数并重跑低增益正式曲线。
3. 对 Fig. 11 `tau=2/3` 做更大 `M`/`L` 的确认实验。
4. 修正 `scripts/plot_results.py` 的重复 marker key 和 Fig.10/Fig.11 label 兼容性。
