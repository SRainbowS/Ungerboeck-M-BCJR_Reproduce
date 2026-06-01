# FTN 库指南

本文档介绍为复现 Ungerboeck M-BCJR 而新增的 `ftn/` 包。它面向未来的实验脚本和调试工作编写，而非最终的论文成果报告。

## 范围

当前库覆盖第一阶段复现工作：

- 实数 BPSK 符号，映射关系为 `0 -> +1`，`1 -> -1`；
- rRC 脉冲生成以及匹配滤波器自相关抽头；
- 有限记忆 Ungerboeck 匹配滤波器信道模型；
- 对数域 BPSK 分支度量；
- 用于短帧验证的全状态 BCJR；
- 所提出的带有局部未来符号搜索的 M-BCJR；
- `(7,5)` 4 态非递归卷积编码器与 SISO 译码器；
- 一个无交织的最小化检测器/译码器 turbo 均衡循环；
- 一个未编码冒烟脚本。

该库目前尚未覆盖 16-QAM、简化 M-BCJR 关键尾部搜索、完整的图 7/8 蒙特卡洛自动化、带交织的 turbo 均衡或面向长帧的优化色噪声生成。

## 依赖流程

各模块按以下顺序使用：

```text
ftn.config
  |
  v
ftn.pulse  ->  ftn.channel  ->  ftn.metrics
                                  |
                                  v
                         ftn.equalizers.full_bcjr
                         ftn.equalizers.ungerboeck_mbcjr
                                  |
                                  v
ftn.modulation  ->  ftn.coding.conv  ->  ftn.turbo_equalization
                                  |
                                  v
                            scripts.run_smoke
```

`anoma/` 保留作为旧版 ANOMA/PCMA 代码路径的参考实现。`ftn/` 包应保持独立，因为论文模型是单用户 FTN ISI，而 `anoma/` 是两用户异步 NOMA/PCMA 模型。

## 模块参考

### `ftn.config`

`FtnConfig` 存储第一阶段实验参数：

- `tau`：FTN 压缩因子。
- `rolloff`：rRC 滚降系数，默认 `0.3`。
- `pulse_span`：以符号周期为单位的半跨度，默认 `15`，与论文一致。
- `sps`：用于数值积分的每符号样点数。
- `isi_len`：接收机 ISI 记忆长度，同时也是网格状态长度。
- `future_len`：M-BCJR 的局部未来符号搜索深度。
- `m_states`：M-BCJR 保留的幸存状态数。
- `snr_db`：单位能量 BPSK 的 `Es/N0`，单位 dB。
- `frame_bits`、`seed`、`turbo_iters`、`llr_clip`：仿真控制参数。

`FtnConfig.n0` 针对单位能量 BPSK 符号返回 `1 / 10^(snr_db/10)`。

### `ftn.modulation`

请使用这些辅助函数，而非手动编写 BPSK 转换：

- `bpsk_modulate(bits)`：映射 `0 -> +1`，`1 -> -1`。
- `bpsk_hard_bits_from_llr(llr)`：从 LLR `log P(+1) / P(-1)` 得到硬判决。

不要直接在 `uint8` 数组上使用 `1 - 2 * bits`。对于比特 `1`，它会下溢并在类型转换前产生 `255`。

### `ftn.pulse`

- `rrc_pulse(t, beta=0.3, T=1.0)` 计算根升余弦脉冲。
- `generate_rrc(beta=0.3, span=15, sps=128, T=1.0)` 在 `[-span*T, span*T]` 上返回 `(t, h)`，并归一化使得 `integral |h(t)|^2 dt = 1`。
- `compute_g(t, h, tau, isi_len, T=1.0)` 计算匹配滤波器自相关抽头：

```text
g_l = integral h(t) h(t - l tau T) dt， -isi_len <= l <= isi_len
```

返回对象是一个以有符号滞后为键的字典，例如 `g[-1]`，`g[0]`，`g[1]`。

### `ftn.channel`

- `build_isi_matrix(g, n)` 构建有限托普利兹矩阵 `G`，满足 `G[row, col] = g[row - col]`。
- `ftn_filter_output(x, g)` 使用有限抽头计算无噪声匹配滤波器输出 `y = Gx`。
- `sample_colored_noise(g, n, n0, rng=None, size=None)` 采样实色噪声，协方差为 `(N0/2) * G`。
- `ftn_awgn_channel(x, g, n0, rng=None)` 返回 `Gx + eta`。

当前噪声采样器使用密集协方差分解。它对于短帧和测试是正确的，但对于 `K=6000` 的完整蒙特卡洛扫描来说太慢了。长序列运行需要基于 FFT/滤波器的色噪声生成方法。

### `ftn.metrics`

- `log_prior_bpsk_symbol(symbol, la)` 返回 BPSK 的对数先验概率，其中 `la = log P(+1) / P(-1)`。
- `log_phi_ungerboeck(y_n, x_n, state_prev, g, n0)` 实现论文中的 Ungerboeck 分支度量：

```text
(2/N0) * x_n * (y_n - 0.5*g_0*x_n - sum_l g_l*x_{n-l})
```

- `log_gamma_bpsk(...)` 将 BPSK 先验加到 `log_phi_ungerboeck` 上。

`state_prev` 的顺序是从最旧到最新。对于 `isi_len=3`，`state_prev[-1]` 是 `x_{n-1}`，`state_prev[-3]` 是 `x_{n-3}`。

### `ftn.equalizers.full_bcjr`

此模块是短 BPSK 序列的正确性基准。

- `full_bcjr_bpsk(y, g, n0, la=None, isi_len=None, initial_state=None)` 在所有 `2^isi_len` 个状态上执行精确的全状态对数域 BCJR。
- `brute_force_bpsk_llr(...)` 枚举每个 BPSK 序列并计算 MAP LLR。它仅适用于极短序列。

两者返回/使用的 LLR 约定均为 `log P(+1|y) / P(-1|y)`。

### `ftn.equalizers.ungerboeck_mbcjr`

`ungerboeck_mbcjr_bpsk(...)` 实现所提出的 M-BCJR 状态剪枝：

1. 将每个幸存状态扩展为 `+1` 和 `-1`；
2. 使用 `logaddexp` 合并到达同一目的地的状态；
3. 通过枚举 `2^future_len` 条尾部路径计算局部未来度量；
4. 按 `log_alpha + log_beta` 对状态排序；
5. 最多保留 `m_states` 个幸存状态。

当 `m_states >= 2^isi_len` 且 `future_len` 覆盖剩余帧时，在短帧测试中该算法与全 BCJR 一致。

### `ftn.coding.conv`

此模块实现论文第一阶段的卷积码：

- `conv_encode_75(bits, initial_state=0)` 使用速率 1/2 的 4 态非递归 `(7,5)` 卷积码进行编码。
- `conv_bcjr_decode(code_llr, initial_state=0, final_state_known=None)` 执行 SISO log-MAP 译码。

重要的约定差异：

- 检测器符号 LLR 为 `log P(+1) / P(-1)`；
- 卷积码比特 LLR 为 `log P(bit=0) / P(bit=1)`。

由于 BPSK 映射 `0 -> +1` 和 `1 -> -1`，这两种约定在检测器/码比特边界上数值上是相容的。

### `ftn.turbo_equalization`

`turbo_equalize_conv_bpsk(...)` 运行一个最小化的检测器/译码器信息交换：

```text
M-BCJR 检测器 LLR
  -> 检测器外信息 = 检测器后验 - 检测器先验
  -> 卷积 SISO 译码器
  -> 译码器外信息 = 译码器码字后验 - 检测器外信息
  -> 下一次检测器先验
```

这目前是一个最小化的无交织 BPSK 链。它足以验证有限 LLR 和基本迭代行为，但图 7/图 8 的复现仍需显式的交织和 BER 曲线脚本。

## 典型用法

### 构建 FTN 信道抽头

```python
from ftn.config import FtnConfig
from ftn.pulse import generate_rrc, compute_g

config = FtnConfig(tau=0.8, isi_len=3, m_states=8, future_len=4)
t, h = generate_rrc(
    beta=config.rolloff,
    span=config.pulse_span,
    sps=config.sps,
)
g = compute_g(t, h, tau=config.tau, isi_len=config.isi_len)
```

### 运行未编码 BPSK M-BCJR 检测

```python
import numpy as np

from ftn.channel import ftn_awgn_channel
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_hard_bits_from_llr, bpsk_modulate

rng = np.random.default_rng(config.seed)
bits = rng.integers(0, 2, size=200, dtype=np.uint8)
x = bpsk_modulate(bits)
y = ftn_awgn_channel(x, g, n0=config.n0, rng=rng)

det = ungerboeck_mbcjr_bpsk(
    y,
    g,
    n0=config.n0,
    isi_len=config.isi_len,
    m_states=config.m_states,
    future_len=config.future_len,
    initial_state=tuple(1.0 for _ in range(config.isi_len)),
)
hard_bits = bpsk_hard_bits_from_llr(det.llr)
```

### 在短序列上对比 M-BCJR 与全 BCJR

```python
from ftn.equalizers.full_bcjr import full_bcjr_bpsk
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk

full = full_bcjr_bpsk(y[:8], g, n0=config.n0, isi_len=config.isi_len)
mbcjr = ungerboeck_mbcjr_bpsk(
    y[:8],
    g,
    n0=config.n0,
    isi_len=config.isi_len,
    m_states=2 ** config.isi_len,
    future_len=8,
)
```

### 运行冒烟脚本

```powershell
python scripts\run_smoke.py --output-dir results\smoke --frame-symbols 80 --frames-per-snr 3 --snr-db 0 2 4
```

该脚本会输出：

- `uncoded_smoke.csv`
- `uncoded_smoke_config.json`

## 测试映射

每个模块都有针对性的测试：

- `tests/test_pulse.py`：rRC 能量、自相关对称性、奈奎斯特抽头。
- `tests/test_channel.py`：托普利兹矩阵构建和色噪声协方差。
- `tests/test_metrics.py`：先验 LLR 和分支度量状态索引。
- `tests/test_modulation.py`：BPSK 映射和硬判决。
- `tests/test_full_bcjr.py`：全 BCJR 对比穷举 MAP。
- `tests/test_mbcjr.py`：保留全部状态时 M-BCJR 对比全 BCJR。
- `tests/test_conv.py`：`(7,5)` 编码器和 SISO 译码器。
- `tests/test_turbo_equalization.py`：有限 LLR 和基本迭代行为。
- `tests/test_smoke.py` 和 `tests/test_smoke_cli.py`：冒烟脚本 API 和命令行接口。

运行所有测试：

```powershell
python -m pytest tests
```

如果 pytest 发出关于 `.pytest_cache` 权限的警告，测试可能仍然有效。该项目之前在受限 ACL 下产生过根级 `pytest-cache-files-*` 临时目录；当它们出现时可以安全删除。

## 当前限制与后续工作

1. 在 `turbo_equalization.py` 中添加显式交织。
2. 使用 `FtnConfig` 添加 `scripts/reproduce_fig7.py` 和 `scripts/reproduce_fig8.py`。
3. 替换长帧的密集色噪声采样。
4. 为 `K=6000` 扫描添加优化的 M-BCJR 数据结构。
5. 实现简化的 M-BCJR 关键尾部搜索。
6. 为 16-QAM 扩展调制和 LLR 汇聚。
