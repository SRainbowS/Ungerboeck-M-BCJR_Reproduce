下面这份计划按 **Codex 可执行工程任务** 来写，目标是复现 Prlja & Anderson 2012 这篇论文的核心低复杂度 ISI 检测器：**minimum/super-minimum phase 前端 + simple-detection M-BCJR + backup/smoothed backup M-BCJR**。论文网页摘要确认其贡献是提出新的 **M-algorithm BCJR, M-BCJR**，用于 FTN 引入的 severe ISI；在迭代译码场景中，为了获得准确 LLR，还提出 **3-recursion M-BCJR**。([Lund University](https://portal.research.lu.se/en/publications/reduced-complexity-receivers-for-strongly-narrowband-intersymbol-/?utm_source=chatgpt.com)) 你上传的论文第 1 页摘要、第 5 页 BCJR 公式、第 7–8 页 M-BCJR/backup M-BCJR 描述、第 9–10 页 BER 实验图，是复现的主要依据。

## 1. 先明确复现范围

这篇论文不是把完整 BCJR “无损降复杂度”，而是做 **reduced-search 近似最优检测**。完整 BCJR 对长度 $m_T$ 的二进制 ISI 信道需要 $2^{m_T}$ 个状态；论文用 M-BCJR 每个时刻只保留 $M$ 个主路径，把复杂度近似降到 $O(M)$ 每 trellis stage。论文自己的说法是：reduced complexity 可以通过 **reduced-trellis** 或 **reduced-search** 两条路实现，本文重点是后者，offset VA/BCJR 和 truncated BCJR 作为 benchmark。

你的复现建议分两层：

**第一层：simple ISI detection 复现。**
复现 Fig. 7 思路：给定 FTN 等效 ISI taps，比较 full BCJR / full VA、offset VA/BCJR、simple M-BCJR 的 BER/EER。这个是最适合先做的核心低复杂度检测器。

**第二层：turbo equalization 复现。**
复现 Fig. 9–10 思路：加入 $(7,5)$ 卷积码、交织器、外层 BCJR 译码器，使用 smoothed backup M-BCJR 做 ISI equalizer。这个工作量明显更大，建议作为第二阶段。

---

## 2. 推荐 Python 环境

```bash
conda create -n ftn-mbcjr python=3.11 -y
conda activate ftn-mbcjr

pip install -U pip setuptools wheel

pip install \
  numpy==1.26.4 \
  scipy==1.11.4 \
  matplotlib==3.8.4 \
  pandas==2.2.2 \
  tqdm==4.66.5 \
  pytest==8.3.3 \
  jupyterlab==4.2.5 \
  ipykernel==6.29.5
```

Codex 官方资料说明 Codex CLI 是能在本地代码库运行的 coding agent；建议在项目根目录放 `AGENTS.md`，把复现目标、运行命令、数学约定和测试标准写清楚，这样 Codex 每轮修改更稳定。([GitHub](https://github.com/openai/codex?utm_source=chatgpt.com))

---

## 3. 建议仓库结构

```text
ftn-mbcjr-repro/
├── AGENTS.md
├── README.md
├── pyproject.toml
├── environment.yml
├── src/
│   └── ftn_mbcjr/
│       ├── __init__.py
│       ├── ftn_models.py          # 论文给出的 FTN 等效 ISI tap
│       ├── modulation.py          # BPSK 映射
│       ├── channel.py             # ISI 卷积 + AWGN
│       ├── trellis.py             # 状态枚举、分支、branch label
│       ├── full_bcjr.py           # 完整 BCJR baseline
│       ├── full_viterbi.py        # 完整 VA baseline
│       ├── m_bcjr.py              # simple detection M-BCJR
│       ├── backup_m_bcjr.py       # backup / smoothed backup M-BCJR
│       ├── offset_receivers.py    # offset VA / offset BCJR benchmark
│       ├── convcode.py            # (7,5) 卷积码，第二阶段
│       ├── turbo_equalizer.py     # 第二阶段
│       ├── metrics.py             # BER/EER/LLR/归一化
│       └── plotting.py
├── experiments/
│   ├── fig7_simple_detection.py
│   ├── fig9_turbo_tau05.py
│   ├── fig10_turbo_tau035.py
│   └── sanity_full_vs_mbcjr.py
├── tests/
│   ├── test_ftn_models.py
│   ├── test_channel.py
│   ├── test_full_bcjr.py
│   ├── test_m_bcjr.py
│   ├── test_backup_m_bcjr.py
│   └── test_end_to_end.py
└── results/
    ├── data/
    └── figures/
```

---

## 4. `AGENTS.md` 给 Codex 的项目说明

```markdown
# AGENTS.md

## Goal

Reproduce the reduced-complexity ISI receivers from:
A. Prlja and J. B. Anderson,
"Reduced-Complexity Receivers for Strongly Narrowband Intersymbol Interference Introduced by Faster-than-Nyquist Signaling",
IEEE Transactions on Communications, 2012.

Primary target:
- simple-detection M-BCJR for severe ISI caused by FTN signaling.
- reproduce the behavior of Fig. 7.

Secondary target:
- backup/smoothed backup M-BCJR for turbo equalization.
- reproduce the behavior of Fig. 9 and Fig. 10.

## Coding rules

- Use Python 3.11.
- Use numpy/scipy only for core algorithms.
- Implement BCJR and M-BCJR in the log domain.
- Binary symbols are +1 and -1.
- AWGN variance convention: real-valued channel noise has variance N0/2 if matching the paper's real BPSK model.
- Every public function must include shape conventions.
- Add pytest tests before large simulations.

## Mathematical conventions

- ISI model: y[n] = sum_{k=0}^{mT} a[n-k] v[k] + w[n].
- Branch label: ell = sum_{k=0}^{mT} a[n-k] v[k].
- Full BCJR has 2^mT states.
- M-BCJR keeps only M dominant paths per trellis stage.
- In M-BCJR, beta paths should prioritize overlap with stored alpha paths.
- Backup M-BCJR uses an additional small recursion with MB, typically MB=2.

## Validation

- For M >= 2^mT, M-BCJR should match full BCJR.
- For memoryless channel, full BCJR should match symbol-by-symbol MAP.
- BER/EER should decrease with SNR.
- Larger M should not perform worse than smaller M except Monte Carlo noise.
```

---

## 5. 第一阶段：直接使用论文给出的离散 ISI 模型

不要一开始就从连续 rRC 脉冲生成 FTN 模型。论文第 4 页已经给出主要测试用的 unit-energy discrete-time channel models，尤其是 $\tau=0.703,0.5,0.35,0.25$ 的 taps。它们已经包含 minimum/super-minimum phase 处理后的模型，论文 Fig. 3 展示了这些模型的 tap 分布。

`ftn_models.py` 先写死这些模型：

```python
import numpy as np

def ftn_model_tau_0703() -> np.ndarray:
    return np.array([
        .553, .793, -.084, -.171, .154, -.064, .006,
        .010, -.012, .015, -.016, .013, -.008
    ], dtype=np.float64)

def ftn_model_tau_05() -> tuple[np.ndarray, int]:
    """
    Returns (v, Kp).
    Kp is the number of low-energy precursor taps ignored by the detector.
    """
    v = np.array([
        -.005, -.003, .007, -.011, -.001, .034, -.019, .003,
        .375, .741, .499, -.070, -.214, .019, .087, -.020,
        -.028, .017
    ], dtype=np.float64)
    Kp = 8
    return v, Kp

def ftn_model_tau_035() -> tuple[np.ndarray, int]:
    v = np.array([
        .025, .012, -.024, .008, .191, .464, .623, .506,
        .176, -.123, -.196, -.075, .060, .080, .013,
        -.035, -.022
    ], dtype=np.float64)
    Kp = 4
    return v, Kp

def ftn_model_tau_025() -> tuple[np.ndarray, int]:
    v = np.array([
        -.010, -.013, -.007, .005, .011, .004, -.008, .001,
        .060, .181, .339, .473, .520, .443, .262, .047,
        -.120, -.182, -.138, -.037, .055, .092, .070, .018,
        -.025, -.037, -.021, .003, .016, .012, .0004, -.008
    ], dtype=np.float64)
    Kp = 8
    return v, Kp
```

注意：论文中 $\tau=0.5,0.35,0.25$ 模型有 italic precursor，所有检测器在形成 branch label 时忽略这些 precursor，并以 $K_p$ delay 工作。这个细节非常重要，否则复现结果会偏。

---

## 6. 第二阶段：基础模块

### 6.1 BPSK 调制

```
modulation.py
def bits_to_bpsk(bits):
    # paper uses 0 -> +1, 1 -> -1 in turbo section
    return 1.0 - 2.0 * bits.astype(np.float64)

def bpsk_to_bits(symbols):
    return (symbols < 0).astype(np.int8)
```

### 6.2 ISI 信道

```
channel.py
def isi_filter(a: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    y_clean[n] = sum_k a[n-k] v[k].
    Use np.convolve(a, v) with consistent indexing.
    """

def add_awgn_real(y: np.ndarray, esn0_db: float, es: float = 1.0, rng=None):
    """
    Paper uses real BPSK. For Es/N0:
    N0 = Es / 10^(EsN0_dB/10)
    real AWGN variance = N0/2.
    """
```

### 6.3 Branch label

论文第 4 页公式写的是：

$$
\ell=\sum_{k=0}^{m_T} a_{n-k}v_k
$$

```
trellis.py
def branch_label(state: int, input_bit: int, v_eff: np.ndarray) -> float:
    """
    state stores previous mT binary symbols.
    input_bit maps to current symbol +1 or -1.
    v_eff excludes ignored precursor taps.
    """
```

---

## 7. 第三阶段：完整 BCJR baseline

完整 BCJR 是验收基准。论文第 5 页列出 BCJR 的 $\alpha$、$\beta$、branch metric 和 LLR 公式。

```
full_bcjr.py
def full_bcjr(y, v_eff, n0, priors=None, terminate=True):
    """
    Full-state BCJR for binary ISI.

    State memory = len(v_eff) - 1.
    Number of states = 2**memory.

    Return:
    - posterior_llr: shape (N,)
    - posterior_prob: shape (N, 2)
    - hard_bits
    """
```

建议全部在 log domain：

```python
from scipy.special import logsumexp
```

branch metric：

$$
\log \Gamma_n(i,j) = \log P(a) - \frac{(y_n-\ell_{i,j})^2}{N_0} + \text{constant}
$$

常数项可省略，因为 LLR 比值会抵消。

必须先通过这些测试：

```text
1. v=[1] 时，full_bcjr 等价于 AWGN symbol MAP。
2. 高 SNR 下 BER 接近 0。
3. 每个时刻 posterior 概率归一化。
4. 小状态模型下 full_bcjr 与 brute-force block MAP 一致。
```

---

## 8. 第四阶段：simple-detection M-BCJR

这是论文最核心的低复杂度检测器。论文第 7 页描述了 simple detection M-BCJR 的 forward recursion、backward recursion 和路径保留策略：每层从 $M$ 个非零 $\alpha$ 或 $\beta$ 扩展到最多 $2M$ 条分支，合并相同状态，只保留最大的 $M$ 条；后向递推时，$\beta$ 路径必须优先保留与已存储 $\alpha$ 路径重叠的状态。

```
m_bcjr.py
@dataclass
class SparseTrellisStage:
    states: np.ndarray       # shape (<=M,)
    log_values: np.ndarray   # alpha or beta log values

def m_bcjr_simple(y, v_eff, n0, M, priors=None, terminate=True):
    """
    Reduced-search BCJR using M-algorithm.

    Steps:
    1. Forward alpha M-search.
    2. Store alpha states and alpha values at every stage.
    3. Backward beta M-search.
    4. During beta search, prioritize states that overlap stored alpha states.
    5. Compute LLR using alpha-beta overlaps.
    6. If L+ or L- is empty, use reserve value Lambda.
    """
```

实现细节：

```text
alpha stage:
    candidates = []
    for each retained state:
        for input symbol in [+1, -1]:
            next_state = transition(state, input_symbol)
            metric = alpha + branch_metric
            append candidate
    merge candidates with same next_state by logsumexp
    keep top M by log_value

beta stage:
    similarly expand backward
    merge same previous_state
    first keep candidates whose state overlaps alpha_states[n]
    then fill remaining slots with largest non-overlap candidates
```

LLR：

$$
LLR(a_n)= \log \frac{ \sum_{j\in L_{+1}}\alpha_n[j]\beta_n[j] }{ \sum_{j\in L_{-1}}\alpha_n[j]\beta_n[j] }
$$

如果 $L_{+1}$ 或 $L_{-1}$ 为空，论文 simple M-BCJR 用 reserve value $\Lambda$ 替代。这个适合 simple detection，但不适合 turbo equalization。

---

## 9. 第五阶段：backup / smoothed backup M-BCJR

论文第 8 页指出 simple M-BCJR 做硬判决还可以，但在 turbo equalization 中 LLR 幅度必须准确；当 $L_{+1}$ 或 $L_{-1}$ 为空时，LLR 幅度没有可靠估计。作者提出第三个低复杂度 recursion，称为 backup M-BCJR；图 8 展示了 $\alpha$、$\beta$、hard path 和 backup search 的关系。

```
backup_m_bcjr.py
def backup_m_bcjr(y, v_eff, n0, M, MB=2, smooth=False, priors=None):
    """
    3-recursion M-BCJR:
    - alpha M-search
    - beta M-search with overlap priority
    - backup recursion of size MB to estimate missing LLR magnitude
    """
```

backup recursion 的工程实现建议：

```text
1. 先运行 simple M-BCJR，得到 hard decision path。
2. 找到 L+ 或 L- 为空的位置。
3. 对这些位置，从 hard path 对应节点出发做一个小规模 forward M-search，大小 MB。
4. 将沿 hard branch 的概率作为已判符号概率；
5. 将 incorrect subset 的概率作为反向符号概率；
6. 用这两个概率生成 backup LLR。
```

平滑版本：

```python
def smooth_backup_llrs(llr, mask):
    """
    Only smooth backup LLRs, not all LLRs.
    Paper uses a simple smoother like {1,3,1}/5 in first iteration.
    """
```

论文 turbo 部分还使用 loop gain $g\le 1$，并将 extrinsic LLR 按 $\sqrt g$ 缩放；$\tau=0.5,0.35$ 附近推荐 $g\approx 0.4$，$\tau=0.25$ 附近推荐 $g\approx 0.25$。这些是第二阶段复现 Fig. 9–10 时要加入的参数。

---

## 10. 第六阶段：benchmark 接收机

为了证明 M-BCJR 的优势，论文比较了三类 benchmark：

### 10.1 Truncated BCJR

只使用前 $m_{\text{tr}}+1$ 个 dominant taps 形成 branch label，不补偿尾巴：

```python
def truncated_bcjr(y, v_eff, n0, m_tr):
    v_tr = v_eff[:m_tr+1]
    return full_bcjr(y, v_tr, n0)
```

这个方法简单，但强 ISI 下性能会明显差，因为丢失了长尾能量。

### 10.2 Offset VA

论文第 5–6 页给出 offset branch label：

$$
\ell = \sum_{k=0}^{m} a_{n-k}v_k + \sum_{k=m+1}^{m_T} a_{n-k}v_k
$$

前一项进入 main state，后一项通过 offset state 修正 branch label。

### 10.3 Single-offset BCJR

这是论文为了公平比较 BCJR benchmark 做的改进：不同 main state 不使用不同 offset，而是使用 single soft offset；软估计来自 $\alpha$ 的概率和，例如论文公式：

$$
\hat a_{n-m}=\hat p_{+1}-\hat p_{-1}
$$

实现上可先简化：只做 truncated BCJR 和 simple M-BCJR；等 Fig. 7 跑通后，再加 offset VA/BCJR。

---

## 11. 实验设计：先复现 Fig. 7

```
experiments/fig7_simple_detection.py
```

目标：复现论文 simple ISI detection 的趋势：$\tau=0.5,0.35,0.25$ 下，M-BCJR 用较小 $M$ 接近 ML/Q-function benchmark，并与 offset VA 做比较。论文第 8 页 Fig. 7 显示 simple M-BCJR 在 $\tau=0.5,0.35,0.25$ 三种 ISI 强度下的 EER 曲线。

最小命令：

```bash
python experiments/fig7_simple_detection.py \
  --tau 0.5 \
  --M-list 2 3 4 5 7 10 \
  --esn0-start 6 \
  --esn0-stop 14 \
  --esn0-step 1 \
  --frame-len 800 \
  --min-error-events 50 \
  --seed 1234
```

对 $\tau=0.35$：

```bash
python experiments/fig7_simple_detection.py \
  --tau 0.35 \
  --M-list 4 5 6 8 10 20 \
  --esn0-start 8 \
  --esn0-stop 18 \
  --esn0-step 1 \
  --frame-len 800 \
  --min-error-events 50 \
  --seed 1234
```

对 $\tau=0.25$：

```bash
python experiments/fig7_simple_detection.py \
  --tau 0.25 \
  --M-list 10 20 40 \
  --esn0-start 12 \
  --esn0-stop 22 \
  --esn0-step 1 \
  --frame-len 800 \
  --min-error-events 50 \
  --seed 1234
```

输出：

```text
results/data/fig7_tau05.csv
results/data/fig7_tau035.csv
results/data/fig7_tau025.csv
results/figures/fig7_tau05.png
results/figures/fig7_tau035.png
results/figures/fig7_tau025.png
```

指标：

```python
def ber(bits_true, bits_hat): ...
def eer(bits_true, bits_hat):
    """
    Paper uses error event rate.
    Simple implementation:
    count contiguous error runs as one event.
    EER = number_error_events / number_bits
    """
```

---

## 12. 第二阶段实验：复现 Fig. 9 和 Fig. 10

这一步加入 turbo equalization。论文第 9 页说明 setup：信息比特先经过 $(7,5)$ rate-1/2 feed-forward convolutional code，再交织，映射到 BPSK；ISI equalizer 和外层卷积码 BCJR 交换 LLR。每个完整 loop 是一次 iteration。论文实验中 block length $N=12000$，一般 20 次迭代，$\tau=0.25$ 用 60 次。

### 12.1 卷积码模块

```
convcode.py
class ConvCode75:
    """
    Rate-1/2 feed-forward convolutional code with generators (7,5) octal.
    Constraint length 3, 4 states.
    """

def conv_encode_75(bits): ...
def conv_bcjr_75(llr_channel, apriori=None): ...
```

### 12.2 Turbo equalizer

```
turbo_equalizer.py
def turbo_equalize(
    y, v_eff, n0, interleaver, M, MB=2, n_iter=20, gain=0.4
):
    """
    Iterative loop:
    1. ISI smoothed backup M-BCJR produces extrinsic LLR.
    2. Deinterleave.
    3. Conv-code BCJR produces decoder extrinsic LLR.
    4. Interleave and feed back to ISI M-BCJR.
    5. Apply sqrt(gain) scaling before each BCJR as in the paper.
    """
```

### 12.3 Fig. 9：$\tau=0.5$

```bash
python experiments/fig9_turbo_tau05.py \
  --M-list 1 2 4 10 \
  --ebn0-start 0 \
  --ebn0-stop 6 \
  --ebn0-step 0.5 \
  --block-len 12000 \
  --n-iter 20 \
  --gain 0.4 \
  --seed 1234
```

### 12.4 Fig. 10：$\tau=0.35$

```bash
python experiments/fig10_turbo_tau035.py \
  --M-list 4 5 6 8 12 16 \
  --ebn0-start 3 \
  --ebn0-stop 8 \
  --ebn0-step 0.5 \
  --block-len 12000 \
  --n-iter 20 \
  --gain 0.4 \
  --seed 1234
```

---

## 13. 必须写的单元测试

### `test_ftn_models.py`

```text
1. 每个 FTN tap 模型能量接近 1。
2. tau=0.5/0.35/0.25 返回合理的 Kp。
3. 去掉 precursor 后 v_eff 的强 tap 应在前部。
```

### `test_channel.py`

```text
1. isi_filter(a, v) 与 np.convolve(a, v) 一致。
2. AWGN 方差在大样本下接近期望值。
```

### `test_full_bcjr.py`

```text
1. v=[1] 时 full BCJR 等价于 AWGN MAP。
2. 高 SNR 下 full BCJR BER 接近 0。
3. posterior probabilities 归一化。
```

### `test_m_bcjr.py`

```text
1. M >= 2^mT 时，M-BCJR 与 full BCJR 一致。
2. M 增大时，BER 不应系统性变差。
3. alpha/beta 每层最多保留 M 个状态。
4. beta overlap priority 确实保留 alpha-overlap states。
```

### `test_backup_m_bcjr.py`

```text
1. 当 L+ 或 L- 为空时，backup recursion 返回有限 LLR。
2. MB=2 时能正常运行。
3. smoothing 只作用于 backup LLR mask。
```

### `test_end_to_end.py`

```text
1. tau=0.5, M=10 的 simple M-BCJR BER 明显优于 M=1。
2. tau=0.35 需要比 tau=0.5 更大的 M。
3. BER/EER 随 Es/N0 增大而下降。
```

---

## 14. 给 Codex 的逐步任务提示词

### Task 1：创建项目骨架

```text
Create a Python package named ftn_mbcjr for reproducing the reduced-complexity M-BCJR receivers from Prlja and Anderson 2012.

Set up:
- pyproject.toml
- src/ftn_mbcjr/
- tests/
- experiments/
- README.md

Use Python 3.11, numpy, scipy, matplotlib, pandas, tqdm, pytest.
Do not implement algorithms yet. Add docstrings and TODOs.
```

### Task 2：实现 FTN taps、BPSK 和 ISI 信道

```text
Implement:
- ftn_models.py with the tau=0.703, 0.5, 0.35, 0.25 tap models from the paper.
- modulation.py with BPSK mapping 0 -> +1, 1 -> -1.
- channel.py with ISI convolution and real AWGN.

Add pytest tests verifying tap energy, convolution correctness, and noise variance.
```

### Task 3：实现完整 BCJR

```text
Implement a full-state log-domain BCJR for binary ISI channels.

Use state memory len(v_eff)-1 and 2^memory states.
Use branch metric:
log_gamma = log_prior - (y[n] - branch_label)^2 / n0

Return posterior LLRs and hard decisions.

Add tests for:
- memoryless AWGN channel
- posterior normalization
- high-SNR BER near zero
```

### Task 4：实现 simple M-BCJR

```text
Implement the simple detection M-BCJR from the paper.

Requirements:
- forward M-search alpha recursion
- merge identical states with logsumexp
- keep top M states per stage
- backward M-search beta recursion
- beta paths must prioritize overlap with stored alpha states
- compute LLR from alpha-beta overlaps
- use finite reserve LLR when L+ or L- is empty

Add tests comparing M-BCJR to full BCJR when M >= full state count.
```

### Task 5：实现 EER/BER 和 Fig. 7 实验

```text
Create experiments/fig7_simple_detection.py.

The script should:
- choose tau in {0.5, 0.35, 0.25}
- remove/ignore precursor taps using Kp
- generate random BPSK frames
- terminate frames with +1 symbols
- pass through ISI + AWGN
- run M-BCJR for multiple M
- compute BER and EER
- save CSV and plot EER vs Es/N0

Add a README command to run tau=0.5 first.
```

### Task 6：实现 backup M-BCJR

```text
Implement backup_m_bcjr.py.

Start from simple M-BCJR.
When L+ or L- is empty, estimate missing LLR magnitude using a third recursion with MB=2 along the decided hard path.

Add optional smoothing with kernel [1,3,1]/5 applied only to backup LLRs.

Add tests that backup M-BCJR returns finite LLRs when simple M-BCJR has empty L sets.
```

### Task 7：实现 turbo equalization

```text
Implement:
- rate-1/2 feed-forward convolutional code with generators (7,5) octal
- full BCJR decoder for this convolutional code
- random interleaver
- iterative turbo equalization loop using smoothed backup M-BCJR as ISI equalizer

Create experiments/fig9_turbo_tau05.py and fig10_turbo_tau035.py.
```

---

## 15. 最小验收标准

先不要追求和论文曲线完全重合；先达到这些结果：

| 验收项              | 期望                                     |
| ------------------- | ---------------------------------------- |
| full BCJR           | 小记忆信道下工作正常                     |
| M-BCJR with large M | 与 full BCJR 近似一致                    |
| $\tau=0.5$          | 小 $M$ 即有较好性能                      |
| $\tau=0.35$         | 需要更大 $M$，但 M-BCJR 明显优于很小 $M$ |
| $\tau=0.25$         | 需要更大 $M$，体现 severe ISI            |
| backup M-BCJR       | LLR 不再出现无穷或空集合失败             |
| turbo equalization  | smoothed backup M-BCJR 能随迭代改善 BER  |

---

## 16. 复现优先级

最推荐的执行顺序是：

1. **只做论文给出的离散 ISI taps**，不要先做 rRC 连续脉冲。
2. **先做 full BCJR**，否则没有基准。
3. **再做 simple M-BCJR**，目标是 Fig. 7。
4. **再做 backup M-BCJR**，解决 LLR 空集合。
5. **最后做 turbo equalization**，目标是 Fig. 9–10。

这样你能在第一周内得到可验证结果，而不是陷入 FTN 前端建模、卷积码 turbo 细节和 M-BCJR 三递推同时调不通的问题。
```