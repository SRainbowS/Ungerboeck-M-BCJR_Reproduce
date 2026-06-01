基于你提供的论文，给出**可用于复现的完整技术梳理、系统模型、算法公式、CPL 规避思想、简化算法、仿真参数、BER 复现路线和代码实现计划**。下面是一份面向代码复现的中文方案。

---

# 1. 论文目标概述

论文：  
**Reduced-Complexity Equalization for Faster-Than-Nyquist Signaling: New Methods Based on Ungerboeck Observation Model**  
Shuangyang Li 等，IEEE Transactions on Communications, 2018.

核心问题：

FTN 信号通过令符号间隔为

$$
\tau T,\quad 0<\tau<1
$$

来突破 Nyquist 速率，从而提高频谱效率。但这样会引入严重 ISI。传统 Forney 模型需要白化滤波器，但当 FTN 符号率超过信号带宽时，频域会出现零区，可能违反 Paley-Wiener 条件，导致白化滤波器不能直接构造。因此论文基于 **Ungerboeck observation model / matched-filter model** 设计低复杂度 M-BCJR 检测器。

论文贡献主要包括：

1. 基于 Ungerboeck 观测模型提出一种新的 reduced-complexity MAP M-BCJR。
2. 在选择 M 个幸存状态时引入若干个“future symbols”，避免传统 Ungerboeck M-BCJR 中的 **correct-path-loss, CPL**。
3. 提出进一步简化的算法：只从一个参考状态搜索“关键尾路径 / correct tail path”，然后用该尾路径计算其他状态的 β 度量。
4. 给出 BPSK、16-QAM、卷积码、Turbo 码下的 BER 仿真结果。
5. 显示相对 Nyquist 系统可获得最高约 186% 频谱效率提升，或最高约 4.5 dB SNR 增益。

---

# 2. FTN 系统模型

## 2.1 编码、交织与调制

信息比特：

$$
\mathbf{u}=[u_1,u_2,\ldots,u_K]^T
$$

经过信道编码得到码字：

$$
\mathbf{c}=[c_1,c_2,\ldots,c_N]^T,\quad c_i\in\{0,1\}
$$

经交织与调制得到符号序列：

$$
\mathbf{x}=[x_1,x_2,\ldots,x_N]^T
$$

BPSK 情况下：

$$
x_n\in\{-1,+1\}
$$

高阶调制如 16-QAM 时：

$$
x_n\in\mathcal{X},\quad |\mathcal{X}|=q
$$

---

## 2.2 FTN 发射信号

FTN 基带发射信号为：

$$
s(t)=\sqrt{\frac{E_s}{T}}\sum_n x_n h(t-n\tau T)
$$

其中：

- $h(t)$：成形脉冲，论文使用 root raised cosine, rRC；
- $E_s$：平均符号能量；
- $T$：Nyquist 符号间隔；
- $\tau<1$：FTN 压缩因子；
- $\tau=1$ 时为 Nyquist 传输；
- $\tau<1$ 时符号更密集，产生 ISI。

论文主要参数：

$$
\beta_{\text{roll-off}}=0.3
$$

脉冲截断范围：

$$
[-15T,15T]
$$

这点很重要。论文特别指出，如果脉冲截断过早，会人为改善最小欧氏距离，导致 BER 结果虚假偏好。

---

## 2.3 AWGN 信道与匹配滤波输出

接收信号：

$$
r(t)=s(t)+w(t)
$$

其中 $w(t)$ 为 AWGN，单边功率谱密度为 $N_0$。

接收端经过匹配滤波器和以 $1/(\tau T)$ 采样后得到：

$$
\mathbf{y}=G\mathbf{x}+\boldsymbol{\eta}
$$

其中 $G$ 是由脉冲自相关采样构成的 Toeplitz 矩阵。

自相关系数：

$$
g_n=\int_{-\infty}^{\infty}h(t)h^*(t-n\tau T)\,dt
$$

接收端考虑有限 ISI 长度：

$$
-L_I\le n\le L_I
$$

噪声向量：

$$
\boldsymbol{\eta}
$$

其协方差为：

$$
E[\boldsymbol{\eta}\boldsymbol{\eta}^H]=\frac{N_0}{2}G
$$

因此匹配滤波输出的噪声是**有色噪声**，这正是 Ungerboeck 模型直接处理的对象。

---

# 3. Ungerboeck Observation Model

Ungerboeck 观测模型不使用白化滤波器，而是基于匹配滤波输出直接构造分支度量。

似然函数可分解为：

$$
P(\mathbf{y}|\mathbf{x})\propto \prod_n \phi(y_n,\mathbf{x})
$$

其中分支度量：

$$
\phi(y_n,S_n,S_{n-1})
=
\exp\left\{
\frac{2}{N_0}
\operatorname{Re}
\left[
x_n^*
\left(
y_n-\frac{1}{2}g_0x_n-\sum_{l=1}^{L_I}g_lx_{n-l}
\right)
\right]
\right\}
$$

如果 $S_{n-1}$ 和 $S_n$ 在 trellis 中合法连接，则使用上式；否则分支度量为 0。

状态定义为：

$$
S_n=(x_{n-L_I+1},x_{n-L_I+2},\ldots,x_n)
$$

对于 BPSK，如果完整 trellis 长度为 $L_I$，总状态数为：

$$
2^{L_I}
$$

高阶调制时为：

$$
q^{L_I}
$$

这会导致最优 BCJR/Viterbi 复杂度迅速爆炸，因此需要 M-algorithm 进行状态裁剪。

---

# 4. 传统 Ungerboeck M-BCJR 的问题：CPL

传统 M-BCJR 在每个时刻只保留度量最大的 $M$ 个状态或路径。

在 Forney 模型下，噪声消失时正确路径通常具有最大度量。但在 Ungerboeck 模型下，由于噪声有色且分支度量包含非因果 ISI 影响，可能出现：

> 即使在无噪声情况下，正确路径的度量也不是最大。

这就是 **correct-path-loss, CPL**。

一旦正确路径在早期被 M-algorithm 裁掉，后续即使信噪比很高也无法恢复，导致误码性能恶化。

论文的核心改进是：

> 选择第 $n$ 时刻的 $M$ 个幸存状态时，不只看过去路径的 α 度量，而是用近似后验概率  
> $$
> P(S_n=s|\mathbf{y})
> $$
> 来选择状态。  
> 为了近似该后验概率，需要额外考虑未来 $L$ 个符号。

---

# 5. Proposed Ungerboeck M-BCJR

## 5.1 MAP 状态选择准则

在第 $n$ 个 trellis section，应保留后验概率最大的 $M$ 个状态：

$$
S_n^{\text{surv}}=\operatorname{TopM}_{s}\ P(S_n=s|\mathbf{y})
$$

因为

$$
P(S_n=s|\mathbf{y})\propto P(S_n=s,\mathbf{y})
$$

论文将其近似为：

$$
P(S_n=s,\mathbf{y})
\approx
\alpha_n(s)\beta_n(s)
$$

其中：

- $\alpha_n(s)$：从起点到当前状态的前向路径概率；
- $\beta_n(s)$：从当前状态继续向后 $L$ 个符号的局部未来概率；
- $L$：future-symbol search depth。

---

## 5.2 分支度量

定义：

$$
\gamma_n(S_{n-1}=s',S_n=s)
=
\phi(y_n,s,s')P(x_n)
$$

其中 $P(x_n)$ 来自先验 LLR。若无 turbo 迭代或首次迭代：

$$
P(x_n=+1)=P(x_n=-1)=\frac{1}{2}
$$

BPSK 下，如果 detector 输入的先验 LLR 为 $L_a(x_n)$，则：

$$
P(x_n=+1)=\frac{e^{L_a(x_n)}}{1+e^{L_a(x_n)}}
$$

$$
P(x_n=-1)=\frac{1}{1+e^{L_a(x_n)}}
$$

实际实现建议在 log 域中使用：

$$
\log P(x_n=+1)= -\log(1+e^{-L_a(x_n)})
$$

$$
\log P(x_n=-1)= -\log(1+e^{L_a(x_n)})
$$

---

## 5.3 前向递推 α

$$
\alpha_n(s)
=
\sum_{s'\in S_{n-1}^{\text{surv}}}
\alpha_{n-1}(s')\gamma_n(s',s)
$$

其中只从上一个时刻的 $M$ 个幸存状态扩展。

对于 BPSK，一个状态可以扩展出 2 个新状态；对于 16-QAM，一个状态扩展出 16 个新状态。

扩展后可能出现多个路径合并到同一状态，此时应将概率相加：

$$
\alpha_n(s)=\sum_{\text{merged paths}}\alpha_n^{(i)}(s)
$$

log 域实现：

$$
\log\alpha_n(s)
=
\operatorname{logsumexp}_{s'}
\left[
\log\alpha_{n-1}(s')+\log\gamma_n(s',s)
\right]
$$

---

## 5.4 局部未来 β

论文不是做传统全局 backward recursion，而是对每个候选状态 $s$，向未来展开 $L$ 个符号：

$$
\beta_n(s)
=
\sum_{\mathbf{x}_{n+1}^{n+L}}
\prod_{k=1}^{L}
\gamma_{n+k}(S_{n+k-1},S_{n+k})
$$

其中枚举所有未来符号组合：

$$
\mathbf{x}_{n+1}^{n+L}\in \mathcal{X}^L
$$

因此复杂度随 $q^L$ 指数增长。

log 域：

$$
\log\beta_n(s)
=
\operatorname{logsumexp}_{\mathbf{x}_{n+1}^{n+L}}
\sum_{k=1}^{L}
\log\gamma_{n+k}(S_{n+k-1},S_{n+k})
$$

---

## 5.5 后验状态概率

$$
P(S_n=s|\mathbf{y})\propto \alpha_n(s)\beta_n(s)
$$

log 域：

$$
\log P(S_n=s|\mathbf{y})
=
\log\alpha_n(s)+\log\beta_n(s)+C
$$

排序选择最大 $M$ 个状态。

---

## 5.6 BPSK LLR 输出

对于 BPSK：

$$
L(x_n)
=
\ln
\frac{P(x_n=+1|\mathbf{y})}
{P(x_n=-1|\mathbf{y})}
$$

论文近似为：

$$
L(x_n)
\approx
\ln
\frac{
\sum_{s\in\mathcal{S}_{+1}}
P(S_n=s,\mathbf{y})
}
{
\sum_{s\in\mathcal{S}_{-1}}
P(S_n=s,\mathbf{y})
}
$$

其中：

- $\mathcal{S}_{+1}$：由当前输入 $x_n=+1$ 导致的状态集合；
- $\mathcal{S}_{-1}$：由当前输入 $x_n=-1$ 导致的状态集合。

log 域实现：

$$
L(x_n)
=
\operatorname{logsumexp}_{s:x_n=+1}
\left[
\log\alpha_n(s)+\log\beta_n(s)
\right]
-
\operatorname{logsumexp}_{s:x_n=-1}
\left[
\log\alpha_n(s)+\log\beta_n(s)
\right]
$$

外信息：

$$
L_e(x_n)=L(x_n)-L_a(x_n)
$$

用于 turbo equalization。

---

# 6. CPL 规避思想

传统 M-BCJR 裁剪依据过去路径度量：

$$
\alpha_n(s)
$$

而论文提出裁剪依据：

$$
\alpha_n(s)\beta_n(s)
$$

即加入未来 $L$ 个符号。

直观解释：

- CPL 发生的原因是当前分支度量不足以体现完整 ISI 影响；
- 某些“错误路径”短期内度量较大，但未来符号展开后会暴露错误；
- 通过 β 搜索未来符号，可以补偿 Ungerboeck 模型的非因果影响；
- 只要 $M$ 和 $L$ 合理，正确路径不会被过早裁掉。

论文中的结论可以理解为：

1. 增大 $L$ 会改善 CPL 规避能力；
2. 增大 $M$ 也能降低正确路径被裁掉的概率；
3. $M$ 和 $L$ 存在替代关系：可以用更大的 $M$ 替代更大的 $L$，以降低 $q^L$ 带来的复杂度。

---

# 7. 简化算法

原始算法的主要瓶颈是：

$$
\beta_n(s)
=
\sum_{\mathcal{X}^L}
\prod_{k=1}^{L}\gamma_{n+k}
$$

复杂度为：

$$
O(q^L)
$$

每个候选状态都要做一次，非常重。

论文提出简化思想：

> 对某一个参考状态 $S_n=s_0$，枚举所有 $q^L$ 条未来路径，找到度量最大的“关键尾路径”。然后对其他候选状态，不再枚举所有未来路径，而是只沿这条关键尾路径计算 β。

---

## 7.1 关键尾路径搜索

对参考状态 $s_0$，计算所有未来路径：

$$
v_i=\mathbf{x}_{n+1}^{n+L},\quad i=1,\ldots,q^L
$$

路径度量：

$$
J(v_i)=
\prod_{k=1}^{L}
\phi(y_{n+k},S_{n+k},S_{n+k-1})
P(x_{n+k})
$$

实际 log 域：

$$
\log J(v_i)=
\sum_{k=1}^{L}
\log\gamma_{n+k}
$$

选最大者：

$$
\hat{v}=\arg\max_i \log J(v_i)
$$

---

## 7.2 简化 β 计算

对于所有候选状态 $s$，使用同一个未来符号序列 $\hat{v}$：

$$
\log\beta_n(s)
=
\sum_{k=1}^{L}
\log\gamma_{n+k}(S_{n+k-1},S_{n+k})
$$

其中 $S_{n+k}$ 由当前状态 $s$ 和未来序列 $\hat{v}$ 推导。

这样复杂度从“每个状态都枚举 $q^L$”降为：

1. 对一个参考状态枚举 $q^L$；
2. 对其他状态只计算 1 条路径。

适合高阶调制，例如 16-QAM。

---

# 8. 频谱效率公式

论文定义归一化频谱效率：

$$
\eta=
\frac{(K/N)\log_2(q)}
{(1+\beta_{\text{roll-off}})\tau}
\cdot
\frac{2}{D}
$$

其中：

- $K/N$：编码率；
- $q$：星座大小；
- $\beta_{\text{roll-off}}$：滚降因子，论文为 0.3；
- $\tau$：FTN 压缩因子；
- $D$：调制维度。

相对于正交系统的频谱效率增益：

$$
\text{gain}
=
\frac{\eta_{\text{FTN}}-\eta_{\text{ORTH}}}
{\eta_{\text{ORTH}}}
\times 100\%
$$

例如：

- $\tau=0.5$：相比 $\tau=1$ 频谱效率提升约 100%；
- $\tau=0.35$：提升约 $1/0.35-1\approx 186\%$。

---

# 9. 论文主要仿真参数整理

## 9.1 公共参数

| 项目     | 参数                                |
| -------- | ----------------------------------- |
| 信道     | AWGN                                |
| 成形脉冲 | root raised cosine, rRC             |
| roll-off | $\beta_{\text{roll-off}}=0.3$       |
| 脉冲截断 | $\pm 15T$                           |
| 接收模型 | Ungerboeck matched-filter model     |
| 检测器   | Proposed M-BCJR / Simplified M-BCJR |
| 实现域   | log domain                          |
| 外迭代   | turbo equalization iterations       |

---

## 9.2 BPSK + 卷积码

论文使用：

| 项目       | 参数                                                     |
| ---------- | -------------------------------------------------------- |
| 调制       | BPSK                                                     |
| 信道码     | $(7,5)$ 4-state rate-1/2 nonrecursive convolutional code |
| 信息长度   | $K=6000$                                                 |
| 压缩因子 1 | $\tau=0.5$                                               |
| 压缩因子 2 | $\tau=0.35$                                              |

### $\tau=0.5$

| 参数                          | 值                                       |
| ----------------------------- | ---------------------------------------- |
| Turbo equalization iterations | 5                                        |
| Proposed 参数                 | $M=2,L=3$ 和 $M=2,L=5$                   |
| 参考结果                      | $M=2,L=5$ 约在 4 dB 附近接近 no-ISI 曲线 |
| 频谱效率增益                  | 约 100%                                  |

### $\tau=0.35$

| 参数                          | 值                                    |
| ----------------------------- | ------------------------------------- |
| Turbo equalization iterations | 15                                    |
| Proposed 参数                 | $M=8,L=5$ 和 $M=8,L=7$                |
| 参考结果                      | $M=8,L=7$ 约在 4.5 dB 附近接近 no-ISI |
| 频谱效率增益                  | 约 186%                               |

---

## 9.3 Turbo coded FTN

论文构造非对称 Turbo code。

| 项目       | 参数                                                         |
| ---------- | ------------------------------------------------------------ |
| 信息长度   | $K=21842$                                                    |
| 码率       | $R=1/3$                                                      |
| 码字长度   | $N=65536$，含终止比特                                        |
| 最大迭代   | $I_{\max}=50$                                                |
| 检测器参数 | $M=8,L=5$                                                    |
| 压缩因子   | $\tau=1,2/3,1/2$ 等                                          |
| 结论       | $\tau=2/3$ 相对正交 $R=1/2$ Turbo 码约有 0.3 dB 增益；$\tau=1/2$ 约有 0.4 dB 增益 |

Turbo 码生成多项式：

$$
g_1(D)=
\left[
1,\frac{1+D+D^2}{1+D^2}
\right]
$$

$$
g_2(D)=
\left[
1,\frac{1+D+D^3}{1+D^2+D^3}
\right]
$$

其中一个分量码对 FTN 更友好，另一个类似 W-CDMA Turbo 码分量。

---

## 9.4 16-QAM + 卷积码 + 简化算法

| 项目             | 参数                                                         |
| ---------------- | ------------------------------------------------------------ |
| 调制             | 16-QAM                                                       |
| 外码             | 与 BPSK 卷积码实验相同，$(7,5)$，$K=6000$                    |
| 检测器           | Simplified algorithm                                         |
| Turbo iterations | 10                                                           |
| $\tau$           | $2/3$、0.8 等                                                |
| 论文结果         | $\tau=2/3$ 的 16-QAM FTN 与同编码 64-QAM Nyquist 频谱效率相当，但可获得约 4.5 dB SNR 增益 |
| 论文比较         | $\tau=0.8,M=4,L=3$ 与文献 [27] $L_E=4$、50 次自迭代性能类似，但总复杂度更低 |

---

# 10. 代码复现总体结构

建议用 Python + NumPy/SciPy 实现，必要时用 Numba/Cython 加速。

项目结构建议：

```text
ftn_ungerboeck_reproduce/
├── main_ber_bpsk_conv.py
├── main_ber_16qam_conv.py
├── main_ber_turbo.py
├── config.py
├── ftn/
│   ├── pulse.py
│   ├── channel.py
│   ├── autocorr.py
│   ├── mapper.py
│   └── metrics.py
├── codes/
│   ├── conv_code.py
│   ├── bcjr_decoder.py
│   ├── turbo_code.py
│   └── interleaver.py
├── equalizers/
│   ├── ungerboeck_mbcjr.py
│   ├── ungerboeck_mbcjr_simplified.py
│   └── full_bcjr.py
├── utils/
│   ├── logsumexp.py
│   ├── ber.py
│   └── seed.py
└── results/
    ├── fig7_tau_05/
    ├── fig8_tau_035/
    ├── fig10_turbo/
    └── fig11_16qam/
```

---

# 11. 核心代码模块设计

## 11.1 rRC 脉冲与自相关

函数：

```python
def rrc_pulse(t, T=1.0, beta=0.3):
    ...
```

需要注意奇点：

- $t=0$
- $t=\pm T/(4\beta)$

然后生成截断脉冲：

```python
def generate_rrc(beta=0.3, span=15, sps=64, T=1.0):
    t = np.arange(-span*T, span*T + T/sps, T/sps)
    h = rrc_pulse(t, T, beta)
    h = h / np.sqrt(np.trapz(np.abs(h)**2, t))
    return t, h
```

计算 $g_l$：

```python
def autocorr_samples(t, h, tau, L_I, T=1.0):
    # g_l = integral h(t) h*(t-l tau T) dt
    # 用插值或离散相关实现
    ...
    return g  # dict or array indexed from -L_I to L_I
```

由于 rRC 与其匹配滤波后的总响应为 raised cosine，自相关采样也可直接使用 raised-cosine 公式近似：

$$
g_l = p(l\tau T)
$$

其中 $p(t)$ 是 raised cosine pulse。

---

## 11.2 FTN 信道生成

两种方法：

### 方法 A：矩阵法

$$
\mathbf{y}=G\mathbf{x}+\boldsymbol{\eta}
$$

生成 Toeplitz 矩阵 $G$，再生成有色高斯噪声：

$$
\boldsymbol{\eta}\sim \mathcal{N}(0,\frac{N_0}{2}G)
$$

代码：

```python
def build_G(g, N, L_I):
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            d = i - j
            if abs(d) <= L_I:
                G[i, j] = g[d]
    return G
```

有色噪声：

```python
eta = np.random.multivariate_normal(
    mean=np.zeros(N),
    cov=(N0/2)*G
)
y = G @ x + eta
```

缺点：$N=6000$ 时矩阵太大。

### 方法 B：卷积法，推荐

利用有限 ISI：

$$
y_n = \sum_{l=-L_I}^{L_I} g_l x_{n-l} + \eta_n
$$

噪声 $\eta_n$ 需要按 $G$ 的协方差生成。可用滤波器法近似，或者先做长序列 FFT 频域生成有色噪声。

简化初版复现可先用：

```python
eta_white = np.sqrt(N0/2) * np.random.randn(N)
```

但这不是严格 Ungerboeck 模型。正式复现应生成协方差为 $(N_0/2)G$ 的有色噪声。

---

## 11.3 Ungerboeck 分支度量

log 域：

$$
\log\phi_n
=
\frac{2}{N_0}
\operatorname{Re}
\left[
x_n^*
\left(
y_n-\frac{1}{2}g_0x_n-\sum_{l=1}^{L_I}g_lx_{n-l}
\right)
\right]
$$

代码：

```python
def log_phi(y_n, x_n, past_symbols, g, L_I, N0):
    isi = 0.0
    for l in range(1, L_I + 1):
        isi += g[l] * past_symbols[-l]
    val = np.real(np.conj(x_n) * (y_n - 0.5 * g[0] * x_n - isi))
    return 2.0 / N0 * val
```

---

# 12. Proposed M-BCJR 伪代码

```python
def ungerboeck_mbcjr(y, g, N0, constellation, La, M, L, L_I):
    """
    y: matched filter outputs
    g: ISI coefficients g[0...L_I]
    La: a priori LLR or symbol prior
    M: number of retained states
    L: future search length
    L_I: ISI memory used by receiver
    """

    # 初始化
    survivors = {zero_state: 0.0}  # log alpha
    L_post = np.zeros(len(y))

    for n in range(N):

        # 1. 从上一时刻幸存状态扩展
        candidates = {}

        for state_prev, log_alpha_prev in survivors.items():
            for x_n in constellation:

                state_new = shift_state(state_prev, x_n)

                log_gamma = compute_log_gamma(
                    y, n, state_prev, state_new,
                    x_n, g, N0, La, constellation
                )

                log_alpha_new = log_alpha_prev + log_gamma

                # 合并相同状态
                candidates[state_new] = logsumexp_pair(
                    candidates.get(state_new, -np.inf),
                    log_alpha_new
                )

        # 2. 对每个候选状态计算 beta
        log_post_state = {}

        for state, log_alpha in candidates.items():
            log_beta = compute_log_beta_full(
                y, n, state, g, N0, La,
                constellation, L, L_I
            )
            log_post_state[state] = log_alpha + log_beta

        # 3. 计算当前符号 LLR
        L_post[n] = compute_symbol_llr_from_states(
            log_post_state, constellation
        )

        # 4. 选 Top-M 状态
        survivors = dict(
            sorted(
                candidates.items(),
                key=lambda item: log_post_state[item[0]],
                reverse=True
            )[:M]
        )

        # 可选：归一化避免数值漂移
        norm = logsumexp(list(survivors.values()))
        for s in survivors:
            survivors[s] -= norm

    return L_post
```

---

# 13. 简化 M-BCJR 伪代码

```python
def compute_log_beta_simplified(y, n, candidate_states, g, N0, La,
                                constellation, L, L_I):
    """
    对一个参考状态枚举所有未来路径，得到关键尾路径；
    对其他状态只使用该关键尾路径。
    """

    ref_state = list(candidate_states)[0]

    # 1. 从参考状态枚举 q^L 条尾路径
    best_path = None
    best_metric = -np.inf

    for future_symbols in product(constellation, repeat=L):
        metric = path_metric(
            y, n, ref_state, future_symbols,
            g, N0, La, L_I
        )
        if metric > best_metric:
            best_metric = metric
            best_path = future_symbols

    # 2. 用 best_path 计算每个候选状态的 beta
    log_beta = {}

    for state in candidate_states:
        log_beta[state] = path_metric(
            y, n, state, best_path,
            g, N0, La, L_I
        )

    return log_beta
```

---

# 14. Turbo Equalization 流程

以 BPSK + 卷积码为例：

```text
信息比特 u
  ↓
卷积编码
  ↓
交织
  ↓
BPSK 映射
  ↓
FTN 调制 + AWGN
  ↓
匹配滤波采样得到 y
  ↓
初始化 La = 0
  ↓
for iter = 1...I:
    Ungerboeck M-BCJR detector:
        输入 y, La
        输出 L_post
        Le_det = L_post - La
    解交织
    BCJR channel decoder:
        输入 Le_det
        输出 decoder extrinsic
    交织
    La = interleaved decoder extrinsic
  ↓
硬判决
  ↓
BER
```

---

# 15. BER 复现步骤

## 15.1 Fig. 7 复现：BPSK，$\tau=0.5$

配置：

```python
K = 6000
code = ConvCode(generator=(7, 5), rate=1/2)
mod = BPSK
tau = 0.5
rolloff = 0.3
pulse_span = 15
M_list = [2]
L_list = [3, 5]
turbo_iters = 5
SNR_dB_range = np.arange(2.0, 6.1, 0.25)
```

流程：

1. 随机生成 $K=6000$ 比特；
2. 使用 $(7,5)$ 卷积码编码；
3. 交织；
4. BPSK 映射；
5. 构造 rRC FTN 信道；
6. 使用 $\tau=0.5$ 生成 $g_l$；
7. 通过 FTN 信道生成 $y$；
8. 执行 5 次 turbo equalization；
9. 统计 BER；
10. 对比 no-ISI BPSK 卷积码曲线。

目标现象：

- $M=2,L=5$ 优于 $M=2,L=3$；
- BER 曲线在约 4 dB 附近接近 no-ISI。

---

## 15.2 Fig. 8 复现：BPSK，$\tau=0.35$

配置：

```python
K = 6000
code = ConvCode(generator=(7, 5), rate=1/2)
mod = BPSK
tau = 0.35
M = 8
L_list = [5, 7]
turbo_iters = 15
SNR_dB_range = np.arange(2.5, 6.5, 0.25)
```

目标现象：

- $\tau=0.35$ ISI 更严重；
- 需要更大的 $M$ 和 $L$；
- $M=8,L=7$ 优于 $M=8,L=5$；
- 约 4.5 dB 左右接近 no-ISI。

---

## 15.3 Fig. 10 复现：Turbo coded FTN

配置：

```python
K = 21842
R = 1/3
N = 65536
tau_list = [1.0, 2/3, 0.5]
M = 8
L = 5
Imax = 50
```

Turbo code：

$$
g_1(D)=
\left[
1,\frac{1+D+D^2}{1+D^2}
\right]
$$

$$
g_2(D)=
\left[
1,\frac{1+D+D^3}{1+D^2+D^3}
\right]
$$

停止准则：

```text
如果两个 RSC decoder 给出的信息比特估计完全一致，则提前停止；
否则达到 Imax=50 停止。
```

目标现象：

- $\tau=2/3$ 时 FTN TC 相比同频谱效率的正交 Turbo code 有约 0.3 dB 增益；
- $\tau=1/2$ 时约 0.4 dB 增益；
- 与二进制输入 Shannon limit 的差距约 0.25–0.3 dB。

---

## 15.4 Fig. 11 复现：16-QAM + 简化算法

配置：

```python
K = 6000
code = ConvCode(generator=(7, 5), rate=1/2)
mod = 16QAM
tau_list = [2/3, 0.8]
M = 4
L = 3
turbo_iters = 10
```

重点：

- 使用 simplified M-BCJR；
- 16-QAM 需要 bit-level LLR 输出；
- 可采用 Gray mapping；
- detector 输出每个调制 bit 的 LLR。

16-QAM 下某一 bit 的 LLR：

$$
L(b_m)
=
\log
\frac{
\sum_{x\in\mathcal{X}_{m,0}}
P(x|\mathbf{y})
}{
\sum_{x\in\mathcal{X}_{m,1}}
P(x|\mathbf{y})
}
$$

或根据映射定义 0/1 的符号集合。

目标现象：

- $\tau=2/3$ 的 16-QAM FTN 与同码率 64-QAM Nyquist 频谱效率相当；
- 但可获得约 4.5 dB 的 SNR 优势。

---

# 16. 实现注意事项

## 16.1 必须使用 log 域

普通概率域会出现严重下溢，必须使用：

```python
scipy.special.logsumexp
```

所有 α、β、γ 均在 log 域计算。

---

## 16.2 LLR 裁剪或等效噪声调节

论文提到 M-algorithm 会恶化 LLR 质量，可通过限制 LLR 幅度或调节噪声方差缓解。

建议实现：

```python
LLR_MAX = 20
L_ext = np.clip(L_ext, -LLR_MAX, LLR_MAX)
```

或者在 detector 内使用等效噪声：

$$
\sigma_{\text{eff}}^2
=
\frac{N_0}{2}
+
E_s\sum_{l=L+1}^{L_I}|g_l|^2
$$

---

## 16.3 ISI 长度 $L_I$

论文没有在所有图中直接给出唯一固定的 $L_I$，实际复现应根据 $g_l$ 衰减选择。

建议：

```python
选择最小 L_I，使得
sum_{|l|>L_I} |g_l|^2 / sum_l |g_l|^2 < 1e-4
```

对于严重 FTN，如 $\tau=0.35$，$L_I$ 需要更大。

---

## 16.4 Trellis 终止

论文要求 trellis 起止状态为全 $+1$ 状态。

BPSK 中可以在数据前后补 $L_I$ 个 $+1$ 符号用于终止 ISI 状态。

---

## 16.5 Monte Carlo 终止条件

建议每个 SNR 点：

```python
min_errors = 200
max_bits = 1e7
```

或：

```python
while bit_errors < 200 and total_bits < 1e7:
    run_one_frame()
```

对于低 BER 区域，至少要有足够错误数，否则曲线抖动严重。

---

# 17. 最小可行复现路线

建议按以下顺序实现，不要一开始就做完整 Turbo code。

## 阶段 1：无编码 FTN + 全状态 BCJR

目标：

- 验证 Ungerboeck 分支度量正确；
- 验证 $g_l$、FTN 信道、噪声生成正确。

配置：

```python
BPSK
N = 200
tau = 0.8
L_I = 3 or 4
full trellis states = 2^L_I
```

---

## 阶段 2：无编码 FTN + Proposed M-BCJR

目标：

- 实现 α、β、Top-M；
- 观察 $M,L$ 对 BER 的影响；
- 验证 $L$ 增大时性能改善。

---

## 阶段 3：卷积码 + Turbo Equalization

目标：

- 复现 Fig. 7 和 Fig. 8；
- 使用 $(7,5)$ 码；
- BPSK；
- 参数：
  - $\tau=0.5,M=2,L=3/5,I=5$
  - $\tau=0.35,M=8,L=5/7,I=15$

---

## 阶段 4：Simplified M-BCJR

目标：

- 实现关键尾路径；
- 先在 BPSK 上验证 simplified 与 original 的差距；
- 再扩展到 16-QAM。

---

## 阶段 5：Turbo Code

目标：

- 实现论文中的非对称 Turbo code；
- 复现 Fig. 10。

---

# 18. 预期复现结果表

| 图      | 系统                  | 参数                       | 目标现象                                      |
| ------- | --------------------- | -------------------------- | --------------------------------------------- |
| Fig. 7  | BPSK + $(7,5)$ 卷积码 | $\tau=0.5,M=2,L=3/5,I=5$   | $L=5$ 接近 no-ISI，约 4 dB                    |
| Fig. 8  | BPSK + $(7,5)$ 卷积码 | $\tau=0.35,M=8,L=5/7,I=15$ | $L=7$ 更好，约 4.5 dB 接近 no-ISI             |
| Fig. 10 | Turbo coded FTN       | $M=8,L=5,\tau=2/3,1/2$     | 比同频谱效率正交 Turbo 码约好 0.3–0.4 dB      |
| Fig. 11 | 16-QAM + 简化算法     | $M=4,L=3,\tau=0.8,2/3$     | $\tau=2/3$ 相比 64-QAM Nyquist 约 4.5 dB 优势 |

---

# 19. 关键实现风险

1. **噪声协方差错误**  
   Ungerboeck 模型下匹配滤波噪声是有色噪声，协方差为 $(N_0/2)G$。若用白噪声，结果可能偏差。

2. **脉冲截断过短**  
   必须使用足够长的 rRC 截断，例如论文的 $\pm15T$。

3. **$E_b/N_0$、$E_s/N_0$ 换算**  
   对编码调制系统需统一定义。常见换算：

   $$
   E_s = R\log_2(q)E_b
   $$
   FTN 下若保持发射 PSD，$\tau$ 对能量归一化也要一致处理。
   
4. **LLR 方向与映射一致性**  
   BPSK 若定义 $0\mapsto +1$，则 decoder 的 LLR 符号必须匹配。

5. **状态定义一致性**  
   论文状态为：

   $$
   S_n=(x_{n-L_I+1},\ldots,x_n)
   $$

   写代码时必须保证 shift 方向、past symbol 索引与 $g_lx_{n-l}$ 一致。

---

# 20. 推荐优先交付代码清单

如果你要真正复现，建议先写以下函数：

```python
# pulse.py
rrc_pulse()
raised_cosine_autocorr()
compute_g()

# mapper.py
bpsk_mod()
qam16_mod()
qam16_llr_sets()

# channel.py
ftn_filter_output()
generate_colored_noise()
ftn_awgn_channel()

# metrics.py
log_prior_symbol()
log_phi_ungerboeck()
log_gamma()

# equalizers/ungerboeck_mbcjr.py
shift_state()
expand_states()
compute_log_alpha_candidates()
compute_log_beta_full()
select_top_m()
compute_llr_bpsk()
ungerboeck_mbcjr()

# equalizers/ungerboeck_mbcjr_simplified.py
find_key_tail_path()
compute_log_beta_simplified()
ungerboeck_mbcjr_simplified()

# codes/conv_code.py
conv_encode_75()
bcjr_conv_decode()

# main
run_ber_curve()
```

---

以上内容可以作为“论文复现说明书”。如果你愿意，我下一步可以继续给你生成一版**Python 代码骨架**，包括：

1. `rrc_pulse()`  
2. `compute_g()`  
3. BPSK FTN 信道  
4. Ungerboeck log-branch metric  
5. Proposed M-BCJR 主循环  
6. BER 仿真脚本模板。