下面是一份适合交给 **Codex** 执行的 Python 复现计划，目标是复现 Rusek & Prlja 的 **低复杂度 ISI 信道最佳检测算法**：先实现完整复杂度 BCJR/Viterbi 基线，再实现论文的 **optimal channel shortening / reduced-memory trellis detector**，最后复现 BER/AIR/GMI 曲线。

这篇论文的核心是：构造 **optimal channel shortening**，也叫 **combined linear Viterbi detection**，用于 ISI/MIMO 信道；优化目标不是传统 MMSE，而是从信息论角度最大化 shortened model 的 achievable information rate，并给出 optimal detector 各组成部分的闭式表达式。论文发表于 *IEEE Transactions on Wireless Communications*, Vol. 11, No. 2, pp. 810–818, 2012。([Lund University](https://www.lu.se/lup/publication/254ef85c-a631-4305-8357-418e90892c97)) 这正好对应你要复现的“低复杂度最佳检测”：原始 VA/BCJR 对 ISI 信道的复杂度随信道记忆指数增长，而 channel shortening 的思想是先滤波缩短信道记忆，再在短记忆 trellis 上运行 VA/BCJR。

---

## 0. 复现目标定义

你的 Python 项目最终应该能完成四件事：

1. **构造 ISI 信道矩阵**
   $$
   \mathbf y=\mathbf H\mathbf x+\mathbf n
   $$

2. **实现完整复杂度最优检测器**
   作为 ground truth：完整 BCJR / Viterbi，trellis memory 等于原始 ISI 信道记忆 $L$。

3. **实现论文的 optimal channel shortening detector** 选择 reduced memory $\nu<L$，构造 shortened metric：
   $$
   q(\mathbf y|\mathbf x) = \exp\left( 2\operatorname{Re}\{\mathbf x^\mathrm H\mathbf H_r\mathbf y\} - \mathbf x^\mathrm H\mathbf G_r\mathbf x \right)
   $$
   并强制 $\mathbf G_r$ 为带状矩阵，使检测复杂度从 $|\mathcal X|^L$ 降为 $|\mathcal X|^\nu$。

4. **复现实验曲线**
   至少复现：
   - full-complexity BCJR；
   - $\nu=0$ 线性检测；
   - $\nu=1,2,\dots$ reduced-memory BCJR；
   - BER/SER vs SNR；
   - AIR/GMI vs SNR；
   - EPR4 或 5-tap ISI channel 场景。

---

## 1. 建议仓库结构

```text
optimal-channel-shortening/
├── AGENTS.md
├── README.md
├── pyproject.toml
├── environment.yml
├── src/
│   └── csdet/
│       ├── __init__.py
│       ├── channels.py
│       ├── modulation.py
│       ├── metrics.py
│       ├── shortening.py
│       ├── bcjr.py
│       ├── viterbi.py
│       ├── information_rate.py
│       ├── simulation.py
│       └── plotting.py
├── experiments/
│   ├── reproduce_epr4.py
│   ├── reproduce_5tap.py
│   ├── sweep_memory.py
│   └── validate_against_full_bcjr.py
├── tests/
│   ├── test_channels.py
│   ├── test_bcjr.py
│   ├── test_shortening.py
│   ├── test_information_rate.py
│   └── test_end_to_end.py
├── notebooks/
│   ├── 01_channel_shortening_derivation.ipynb
│   ├── 02_bcjr_validation.ipynb
│   └── 03_reproduce_figures.ipynb
└── results/
    ├── figures/
    ├── logs/
    └── data/
```

Codex CLI 可以在本地终端读取、修改并运行当前目录中的代码；官方文档也建议用 `AGENTS.md` 给 Codex 提供项目级规则。([OpenAI开发者](https://developers.openai.com/codex/cli))

---

## 2. Python 环境

建议先用干净环境，不引入深度学习框架。

```bash
conda create -n csdet python=3.11 -y
conda activate csdet

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

可选：先用 `cvxpy` 验证 $\mathbf G_r$ 的优化，再写高效实现。

```bash
pip install cvxpy==1.5.3
```

---

## 3. `AGENTS.md`：给 Codex 的项目规则

在仓库根目录写入：

```markdown
# AGENTS.md

## Project goal

Reproduce the reduced-complexity ISI channel detection algorithm from:
"Optimal Channel Shortening for MIMO and ISI Channels", Rusek and Prlja, IEEE TWC 2012.

Focus first on the SISO ISI case, not MIMO.

## Coding rules

- Use Python 3.11.
- Use numpy/scipy only for the core algorithm.
- Use complex128 for all complex-valued linear algebra.
- Implement all probabilistic trellis algorithms in the log domain.
- Every public function must have a docstring with equations and shape conventions.
- Add pytest tests for each module.
- Do not use deep learning libraries.
- Do not optimize prematurely; first write a correct reference implementation.

## Mathematical conventions

- The ISI model is y = H x + n.
- x has unit average symbol energy unless explicitly stated.
- n is circular complex Gaussian with covariance N0 I.
- SNR_dB means Es/N0 in dB.
- Full-complexity BCJR uses original ISI memory L.
- Reduced-complexity BCJR uses shortened memory nu.

## Required validation

- For nu >= original channel memory L, the reduced detector must agree with the full detector.
- Increasing nu should not decrease the optimized information-rate objective except for numerical tolerance.
- BCJR probabilities must normalize at every time index.
- BER must decrease as SNR increases in smoke tests.
```

---

## 4. 第一阶段：实现 ISI 信道和调制

### 4.1 `modulation.py`

先支持 BPSK 和 QPSK。

```python
def bpsk_alphabet() -> np.ndarray:
    return np.array([-1.0, 1.0], dtype=np.complex128)

def qpsk_alphabet() -> np.ndarray:
    return np.array([1+1j, 1-1j, -1+1j, -1-1j], dtype=np.complex128) / np.sqrt(2)
```

需要实现：

```text
bits_to_symbols()
symbols_to_bits_hard()
random_symbols()
symbol_energy()
```

### 4.2 `channels.py`

实现 Toeplitz ISI 矩阵：

```python
def isi_convolution_matrix(h: np.ndarray, n_symbols: int) -> np.ndarray:
    """
    Build H such that y = H x + n.
    h[0] is the current-symbol tap.
    H has shape (n_symbols + L, n_symbols) or optionally truncated to n_symbols.
    """
```

建议支持两种边界模式：

```text
mode="full"      y length = N + L
mode="same"      y length = N
```

先用 `full`，避免边界项混乱。

推荐初始信道：

```python
EPR4 = np.array([1, 1, -1, -1], dtype=np.complex128) / 2
FIVE_TAP = np.ones(5, dtype=np.complex128) / np.sqrt(5)
```

EPR4 和 5-tap uniform-power channel 都适合作为论文级 ISI channel shortening 的复现实验场景。

---

## 5. 第二阶段：完整复杂度 BCJR / Viterbi 基线

这是复现的基准。没有它，无法判断 shortened detector 是否正确。

### 5.1 完整 ISI trellis

若信道记忆为 $L$，调制阶数为 $M$，完整 trellis 状态数：

$$
M^L
$$

实现：

```python
def enumerate_states(alphabet: np.ndarray, memory: int) -> np.ndarray:
    """
    Return all states of length memory.
    shape: (M**memory, memory)
    """
```

transition：

```python
def next_state(prev_state, current_symbol):
    return np.concatenate([[current_symbol], prev_state[:-1]])
```

### 5.2 Full BCJR

`bcjr.py`：

```python
def full_isi_bcjr(
    y: np.ndarray,
    h: np.ndarray,
    alphabet: np.ndarray,
    n0: float,
    priors: np.ndarray | None = None,
) -> dict:
    """
    Full-complexity MAP symbol detector for SISO ISI channel.
    Returns posterior symbol probabilities, LLRs if BPSK, and log-likelihood.
    """
```

branch metric：

$$
\gamma_k(s,a) = -\frac{1}{N_0} \left| y_k-\sum_{\ell=0}^{L} h_\ell x_{k-\ell} \right|^2 + \log P(a)
$$

注意全部使用 log-domain：

```python
from scipy.special import logsumexp
```

必须测试：

```text
L=0 时 BCJR 等价于 AWGN symbol-by-symbol MAP。
高 SNR 下 BER 接近 0。
BPSK posterior 概率归一化。
```

### 5.3 Viterbi

`viterbi.py`：

```python
def full_isi_viterbi(y, h, alphabet, n0):
    """
    ML sequence detection baseline.
    """
```

它用于硬判决 BER 对比；BCJR 用于软输出和信息率估计。

---

## 6. 第三阶段：实现论文的 optimal channel shortening

这部分是复现的核心。

论文的 detector class 是 reduced trellis detector：不是在原 trellis 上剪枝，而是构造一个 reduced trellis，再完整处理这个 reduced trellis；网页镜像中的论文摘要和引言也强调了 reduced-complexity trellis detection 的这一区分。([Studocu](https://www.studocu.com/row/document/sadjad-university-of-technology/data-mining/optimal-channel-shortening-of-mimo-and-isi-channels/81821707))

### 6.1 矩阵模型

给定：

$$
\mathbf y=\mathbf H\mathbf x+\mathbf n,\quad \mathbf n\sim\mathcal{CN}(0,N_0\mathbf I)
$$

先计算 LMMSE / Wiener 相关矩阵：

$$
\mathbf B = \mathbf I - \mathbf H^\mathrm H (\mathbf H\mathbf H^\mathrm H+N_0\mathbf I)^{-1} \mathbf H
$$

实现：

```python
def compute_b_matrix(H: np.ndarray, n0: float) -> np.ndarray:
    """
    B = I - H^H (H H^H + N0 I)^-1 H.
    Use solve(), never explicit inverse.
    """
```

数值实现不要写：

```python
np.linalg.inv(A)
```

而是：

```python
X = scipy.linalg.solve(A, H, assume_a="pos")
B = I - H.conj().T @ X
```

### 6.2 求 optimal banded matrix

Reduced memory $\nu$ 对应：

$$
(\mathbf G_r)_{ij}=0,\quad |i-j|>\nu
$$

更方便写：

$$
\mathbf G=\mathbf I+\mathbf G_r
$$

优化目标可由论文中的 GMI/AIR 公式得到。工程上建议分两步：

#### Step A：CVXPY 参考实现

先写一个慢但清楚的版本：

```python
def solve_banded_g_cvxpy(B: np.ndarray, nu: int) -> np.ndarray:
    """
    Reference solver:
    maximize log_det(G) - trace(G B)
    subject to G Hermitian positive definite
               G[i,j] = 0 for |i-j| > nu
    """
```

这个版本只用于小尺寸单元测试，例如 $N=8,12,20$。

#### Step B：闭式/高效 banded 实现

最终版：

```python
def solve_banded_g_optimal(B: np.ndarray, nu: int) -> np.ndarray:
    """
    Efficient implementation of the paper's closed-form optimal
    banded G = I + G_r.

    Validation target:
    - matches solve_banded_g_cvxpy for small N
    - for nu >= N-1, G ≈ inv(B)
    """
```

如果你先不确定论文 Theorem 2 的索引形式，可以让 Codex 执行一个“公式抽取任务”：

```text
Open the uploaded paper and extract the exact Theorem 1/Theorem 2 formulas for H_opt, G_opt, B, and the GMI objective. Then implement solve_banded_g_optimal(B, nu) exactly as stated. Add a CVXPY reference solver and unit tests comparing both solvers for small matrices.
```

### 6.3 构造 channel-shortening filter

论文中的最优前端滤波器形式可以按：

$$
\mathbf H_r = \mathbf G \mathbf H^\mathrm H (\mathbf H\mathbf H^\mathrm H+N_0\mathbf I)^{-1}
$$

其中：

$$
\mathbf G=\mathbf I+\mathbf G_r
$$

实现：

```python
def compute_hr(H: np.ndarray, n0: float, G: np.ndarray) -> np.ndarray:
    """
    H_r = G H^H (H H^H + N0 I)^-1
    """
```

然后：

$$
\mathbf z=\mathbf H_r\mathbf y
$$

Reduced BCJR 不再直接用 $\mathbf y$，而是使用 $\mathbf z$ 和 banded $\mathbf G$。

---

## 7. 第四阶段：Reduced-memory BCJR

Reduced metric：

$$
q(\mathbf y|\mathbf x) = \exp \left( 2\operatorname{Re}\{\mathbf x^\mathrm H\mathbf z\} - \mathbf x^\mathrm H\mathbf G\mathbf x \right)
$$

其中 $\mathbf G$ 是带宽 $\nu$ 的 Hermitian banded matrix。

将二次型展开成 branch metric：

$$
\mathbf x^\mathrm H\mathbf G\mathbf x = \sum_k G_{k,k}|x_k|^2 + 2\operatorname{Re}\left\{ \sum_{i=1}^{\nu} G_{k,k-i}x_k^\ast x_{k-i} \right\}
$$

所以 transition metric 可写成：

$$
\gamma_k(s,a) = 2\operatorname{Re}\{a^\ast z_k\} - G_{k,k}|a|^2 - 2\operatorname{Re}\left\{ \sum_{i=1}^{\nu} G_{k,k-i}a^\ast x_{k-i} \right\} + \log P(a)
$$

实现：

```python
def shortened_bcjr(
    z: np.ndarray,
    G: np.ndarray,
    alphabet: np.ndarray,
    nu: int,
    priors: np.ndarray | None = None,
) -> dict:
    """
    Reduced-memory BCJR using the optimal channel-shortening metric.
    State size is len(alphabet)**nu.
    """
```

必须做 log-domain forward-backward：

```text
alpha[k+1, next_state] = logsumexp(alpha[k, prev_state] + gamma)
beta[k, prev_state]    = logsumexp(gamma + beta[k+1, next_state])
posterior[k, a]        = logsumexp(alpha + gamma + beta)
```

---

## 8. 第五阶段：信息率 / GMI 复现

实现两个层级。

### 8.1 Gaussian-input AIR/GMI

`information_rate.py`：

```python
def capacity_gaussian(H: np.ndarray, n0: float) -> float:
    """
    C = log det(I + H^H H / N0)
    nats per block.
    """
def gmi_channel_shortening(B: np.ndarray, G: np.ndarray) -> float:
    """
    Implement the exact GMI expression from the paper.
    Include a unit test:
    if G = inv(B), result equals capacity_gaussian(H, n0)
    up to numerical tolerance.
    """
```

这里的关键验证是：

```text
nu = N-1 full memory:
G should equal B^{-1}
GMI should equal full Gaussian information rate.
```

### 8.2 Finite-alphabet Monte Carlo GMI

可选但很有价值：

```python
def monte_carlo_gmi_shortened(
    h, alphabet, n0, nu, n_blocks, block_len, seed
):
    """
    Estimate E[ log q(Y|X) / sum_X' P(X') q(Y|X') ].
    For small block length only.
    """
```

这个用于验证 shortened metric 的信息论排序：

$$
I_{\nu=0}\le I_{\nu=1}\le I_{\nu=2}\le I_\text{full}
$$

允许极小数值误差。

---

## 9. 第六阶段：实验复现路线

### 9.1 Smoke test

先小规模：

```bash
python experiments/validate_against_full_bcjr.py \
  --channel epr4 \
  --mod bpsk \
  --block-len 64 \
  --snr-db 6 \
  --nu 0 1 2 3 \
  --nblocks 20 \
  --seed 1
```

期待：

```text
nu=3 与 full BCJR 基本一致，因为 EPR4 memory = 3。
nu=0 最差，但复杂度最低。
nu=1,2 位于二者之间。
```

### 9.2 BER vs SNR

```bash
python experiments/reproduce_epr4.py \
  --mod bpsk \
  --snr-db-start -2 \
  --snr-db-stop 12 \
  --snr-db-step 1 \
  --block-len 256 \
  --nblocks 500 \
  --nu-list 0 1 2 3 \
  --seed 1234
```

输出：

```text
results/data/epr4_ber.csv
results/figures/epr4_ber_vs_snr.png
```

图中应包含：

```text
Full BCJR
CS nu=0
CS nu=1
CS nu=2
CS nu=3
```

### 9.3 AIR/GMI vs SNR

```bash
python experiments/sweep_memory.py \
  --channel epr4 \
  --snr-db-start -5 \
  --snr-db-stop 20 \
  --snr-db-step 1 \
  --block-len 128 \
  --nu-list 0 1 2 3 \
  --metric gmi
```

输出：

```text
results/data/epr4_gmi.csv
results/figures/epr4_gmi_vs_snr.png
```

---

## 10. 必须写的单元测试

### `test_channels.py`

```text
Toeplitz convolution matrix output equals np.convolve().
H shape is correct for mode="full".
EPR4 energy normalization is correct.
```

### `test_bcjr.py`

```text
For L=0 AWGN, BCJR equals symbol-by-symbol MAP.
For noiseless high-SNR BPSK, BER = 0.
Posterior probabilities sum to 1.
```

### `test_shortening.py`

```text
B is Hermitian positive definite.
G is Hermitian positive definite.
G has zero entries outside bandwidth nu.
For nu=N-1, G ≈ inv(B).
CVXPY G and closed-form/numerical G match for N<=12.
```

### `test_information_rate.py`

```text
Full-memory GMI equals Gaussian capacity.
GMI is nondecreasing with nu.
GMI is nondecreasing with SNR.
```

### `test_end_to_end.py`

```text
For EPR4:
nu=3 shortened detector agrees with full BCJR.
nu=2 outperforms nu=1 on average.
BER decreases with SNR.
```

---

## 11. 给 Codex 的分阶段任务提示词

### Task 1：建立工程骨架

```text
Create a Python package named csdet for reproducing the channel-shortening detector from Rusek and Prlja 2012.

Set up:
- pyproject.toml
- src/csdet/
- tests/
- experiments/
- README.md

Use Python 3.11, numpy, scipy, matplotlib, pandas, tqdm, pytest.

Implement only empty module stubs and docstrings first. Do not implement algorithms yet.
```

### Task 2：实现 ISI 信道和调制

```text
Implement modulation.py and channels.py.

Requirements:
- BPSK and QPSK alphabets with unit average symbol energy.
- random symbol generation with reproducible seed.
- Toeplitz ISI convolution matrix H such that y = H x + n.
- Add tests verifying H @ x equals np.convolve(x, h).
```

### Task 3：实现 full BCJR

```text
Implement a log-domain full-complexity BCJR detector for SISO ISI channels.

Inputs:
- y, h, alphabet, n0
- optional symbol priors

Outputs:
- posterior probabilities with shape (N, M)
- hard symbol estimates
- optional BPSK LLRs

Use scipy.special.logsumexp.
Add tests:
- L=0 equals AWGN MAP
- posterior probabilities normalize
- BER decreases at high SNR
```

### Task 4：实现 optimal channel shortening 矩阵

```text
Implement shortening.py.

Functions:
- compute_b_matrix(H, n0)
- solve_banded_g_cvxpy(B, nu) as a reference solver
- solve_banded_g_optimal(B, nu) following the paper's closed-form theorem
- compute_hr(H, n0, G)

Validation:
- G is Hermitian positive definite
- G[i,j] = 0 for |i-j| > nu
- for nu=N-1, G equals inv(B)
- CVXPY and optimal solver match for small random positive definite B
```

### Task 5：实现 shortened BCJR

```text
Implement shortened_bcjr(z, G, alphabet, nu).

Use the metric:
2 Re{x_k^* z_k}
- G[k,k] |x_k|^2
- 2 Re{sum_i G[k,k-i] x_k^* x_{k-i}}

Use log-domain forward-backward.
Add tests comparing shortened_bcjr with full_bcjr when nu equals original ISI memory.
```

### Task 6：实现实验脚本

```text
Create experiments/reproduce_epr4.py and experiments/sweep_memory.py.

The scripts should:
- generate random BPSK symbols
- pass through EPR4 ISI channel
- add complex AWGN
- run full BCJR and shortened BCJR for nu in [0,1,2,3]
- compute BER/SER
- save CSV
- plot BER vs SNR and GMI vs SNR
```

### Task 7：写 README 复现说明

```text
Write README.md with:
- paper citation
- mathematical model
- environment setup
- commands to reproduce EPR4 BER figure
- commands to reproduce GMI figure
- explanation of nu and complexity M^nu
- known limitations
```

---

## 12. 关键实现细节提醒

**第一，边界处理要统一。**
论文公式通常是 block matrix 形式；你的 Python 仿真如果使用 `full` convolution，则 $y$ 长度是 $N+L$。BCJR 和 $H$ 的构造必须保持一致。

**第二，噪声方差定义要固定。**
建议统一：

$$
N_0 = E_s / 10^{\mathrm{SNR}_{dB}/10}
$$

复高斯噪声：

```python
n = np.sqrt(n0 / 2) * (rng.normal(size=y.shape) + 1j * rng.normal(size=y.shape))
```

**第三，BCJR 必须 log-domain。**
普通概率域很快下溢。

**第四，$\nu=0$ 是重要 sanity check。**
这时 trellis 没有记忆，检测器应退化为逐符号检测，复杂度最低。

**第五，$\nu=L$ 是最重要 sanity check。**
当 reduced memory 等于原始 ISI memory，shortened detector 应接近 full detector。

**第六，先不要追求大规模。**
建议先用：

```text
block_len = 64 或 128
nblocks = 20
```

验证正确性，再扩到：

```text
block_len = 256 或 512
nblocks = 500+
```

---

## 13. 推荐的开发顺序

| 阶段  | 目标                            | 完成标准                             |
| ----- | ------------------------------- | ------------------------------------ |
| Day 1 | 仓库、环境、ISI矩阵、调制       | `pytest tests/test_channels.py` 通过 |
| Day 2 | full BCJR / Viterbi             | AWGN 与短 ISI 测试通过               |
| Day 3 | $B$、$\mathbf G$、$\mathbf H_r$ | CVXPY 与闭式/数值解一致              |
| Day 4 | shortened BCJR                  | $\nu=L$ 与 full BCJR 对齐            |
| Day 5 | BER 曲线                        | EPR4 BER vs SNR 可画图               |
| Day 6 | AIR/GMI 曲线                    | GMI 随 $\nu$ 单调增加                |
| Day 7 | README + 论文图复现             | 一条命令生成 CSV 和 figures          |

---

## 14. 最小可复现命令设计

最终你希望 README 里有：

```bash
conda env create -f environment.yml
conda activate csdet

pytest -q

python experiments/reproduce_epr4.py \
  --mod bpsk \
  --block-len 256 \
  --nblocks 500 \
  --snr-db-start -2 \
  --snr-db-stop 12 \
  --snr-db-step 1 \
  --nu-list 0 1 2 3 \
  --seed 1234

python experiments/sweep_memory.py \
  --channel epr4 \
  --block-len 128 \
  --snr-db-start -5 \
  --snr-db-stop 20 \
  --snr-db-step 1 \
  --nu-list 0 1 2 3
```

---

## 15. 最终验收标准

你可以用下面这组标准判断复现是否成功：

1. `pytest -q` 全部通过。
2. $\nu=0$ 的复杂度为 $M^0=1$，表现接近线性逐符号检测。
3. $\nu=1,2$ 的 BER/GMI 介于 $\nu=0$ 和 full BCJR 之间。
4. $\nu=L$ 与 full BCJR 结果几乎一致。
5. GMI/AIR 随 $\nu$ 增加不下降。
6. BER 随 SNR 增加下降。
7. 结果图中能清楚显示：**用少量 trellis memory 可以逼近 full-complexity optimal detection**。

这就是这篇论文最值得复现的结论。
```