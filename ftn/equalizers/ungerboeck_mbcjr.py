from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from ftn.equalizers.full_bcjr import _prepare_la


LOG_ZERO = -np.inf


@dataclass(frozen=True)
class MbcjrResult:
    llr: np.ndarray
    survivor_counts: list[int]


# ---------------------------------------------------------------------------
# Numba JIT kernel — integer-state representation
#
# BPSK state (s_oldest, ..., s_newest) with each s_i in {-1, +1} is encoded
# as an integer where bit k = 1 iff the (isi_len-1-k)-th element is +1.
# Equivalently, bit 0 = newest symbol, bit (isi_len-1) = oldest.
# Shifting in a new symbol: new = ((old << 1) | bit) & mask.
# ---------------------------------------------------------------------------


@njit(cache=True)
def _nb_logaddexp(a: float, b: float) -> float:
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    if a >= b:
        return a + np.log1p(np.exp(b - a))
    return b + np.log1p(np.exp(a - b))


@njit(cache=True)
def _nb_logsumexp(arr, n: int) -> float:
    if n == 0:
        return -np.inf
    mx = arr[0]
    for i in range(1, n):
        if arr[i] > mx:
            mx = arr[i]
    if not np.isfinite(mx):
        return -np.inf
    s = 0.0
    for i in range(n):
        s += np.exp(arr[i] - mx)
    return mx + np.log(s)


@njit(cache=True)
def _nb_log_gamma(y_n: float, x_n: float, state_int: int,
                  g_taps: np.ndarray, isi_len: int,
                  inv_n0: float, la_n: float) -> float:
    # ISI contribution from past symbols
    isi = 0.0
    taps_len = len(g_taps)
    for lag in range(1, isi_len + 1):
        if lag < taps_len:
            sym = 1.0 if (state_int >> (lag - 1)) & 1 else -1.0
            isi += g_taps[lag] * sym
    g0 = g_taps[0] if taps_len > 0 else 0.0
    centered = y_n - 0.5 * g0 * x_n - isi
    phi = 2.0 * inv_n0 * x_n * centered
    # prior
    if x_n > 0:
        prior = -_nb_logaddexp(0.0, -la_n)
    else:
        prior = -_nb_logaddexp(0.0, la_n)
    return phi + prior


@njit(cache=True)
def _nb_beta_local(y: np.ndarray, g_taps: np.ndarray, isi_len: int,
                   inv_n0: float, la: np.ndarray, state_int: int,
                   idx: int, future_len: int,
                   future_table: np.ndarray) -> float:
    n = y.size
    remaining = n - idx - 1
    horizon = future_len if future_len < remaining else remaining
    if horizon <= 0:
        return 0.0

    n_fut = 1 << horizon
    values = np.empty(n_fut)
    mask = (1 << isi_len) - 1

    for fi in range(n_fut):
        st = state_int
        metric = 0.0
        for off in range(horizon):
            symbol = future_table[fi, off]
            j = idx + off + 1
            metric += _nb_log_gamma(y[j], symbol, st, g_taps, isi_len,
                                    inv_n0, la[j])
            bit = 1 if symbol > 0.0 else 0
            st = ((st << 1) | bit) & mask
        values[fi] = metric

    return _nb_logsumexp(values, n_fut)


@njit(cache=True)
def _nb_mbcjr_kernel(y: np.ndarray, g_taps: np.ndarray, n0: float,
                     la: np.ndarray, isi_len: int, m_states: int,
                     future_len: int, future_table: np.ndarray,
                     init_state_int: int):
    n = y.size
    mask = (1 << isi_len) - 1
    inv_n0 = 1.0 / n0

    max_cand = 2 * m_states + 2

    surv_s = np.empty(m_states, dtype=np.int64)
    surv_m = np.empty(m_states, dtype=np.float64)
    surv_s[0] = init_state_int
    surv_m[0] = 0.0
    n_surv = 1

    llr = np.zeros(n, dtype=np.float64)
    surv_counts = np.zeros(n, dtype=np.int64)

    cand_s = np.empty(max_cand, dtype=np.int64)
    cand_m = np.empty(max_cand, dtype=np.float64)
    lpost = np.empty(max_cand, dtype=np.float64)
    plus_v = np.empty(max_cand, dtype=np.float64)
    minus_v = np.empty(max_cand, dtype=np.float64)

    for idx in range(n):
        # --- expand candidates ---
        n_cand = 0
        for si in range(n_surv):
            sp = surv_s[si]
            am = surv_m[si]
            for sym_idx in range(2):
                sym = -1.0 if sym_idx == 0 else 1.0
                sn = ((sp << 1) | sym_idx) & mask
                g_val = _nb_log_gamma(y[idx], sym, sp, g_taps, isi_len,
                                      inv_n0, la[idx])
                metric = am + g_val
                found = False
                for ci in range(n_cand):
                    if cand_s[ci] == sn:
                        cand_m[ci] = _nb_logaddexp(cand_m[ci], metric)
                        found = True
                        break
                if not found:
                    cand_s[n_cand] = sn
                    cand_m[n_cand] = metric
                    n_cand += 1

        # --- beta + posterior ---
        for ci in range(n_cand):
            lpost[ci] = cand_m[ci] + _nb_beta_local(
                y, g_taps, isi_len, inv_n0, la, cand_s[ci],
                idx, future_len, future_table,
            )

        # --- LLR ---
        n_plus = 0
        n_minus = 0
        for ci in range(n_cand):
            if (cand_s[ci] & 1) == 1:
                plus_v[n_plus] = lpost[ci]
                n_plus += 1
            else:
                minus_v[n_minus] = lpost[ci]
                n_minus += 1
        llr[idx] = _nb_logsumexp(plus_v, n_plus) - _nb_logsumexp(minus_v, n_minus)

        # --- prune top-m ---
        keep = m_states if m_states < n_cand else n_cand
        if keep < n_cand:
            order = np.argsort(lpost[:n_cand])[::-1]
        else:
            order = np.arange(n_cand)
        for ki in range(keep):
            surv_s[ki] = cand_s[order[ki]]
            surv_m[ki] = cand_m[order[ki]]
        n_surv = keep

        # --- normalise ---
        norm = _nb_logsumexp(surv_m, n_surv)
        for si in range(n_surv):
            surv_m[si] -= norm

        surv_counts[idx] = n_surv

    return llr, surv_counts


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _g_dict_to_array(g: dict[int, float], isi_len: int) -> np.ndarray:
    taps = np.zeros(isi_len + 1, dtype=np.float64)
    for lag, val in g.items():
        if 0 <= lag <= isi_len:
            taps[lag] = float(val)
    return taps


def _build_future_table(future_len: int) -> np.ndarray:
    if future_len <= 0:
        return np.empty((1, 0), dtype=np.float64)
    n_seq = 1 << future_len
    table = np.empty((n_seq, future_len), dtype=np.float64)
    for i in range(n_seq):
        for j in range(future_len):
            table[i, j] = 1.0 if (i >> j) & 1 else -1.0
    return table


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def ungerboeck_mbcjr_bpsk(
    y: np.ndarray,
    g: dict[int, float],
    n0: float,
    la: np.ndarray | None = None,
    isi_len: int | None = None,
    m_states: int = 4,
    future_len: int = 3,
    initial_state: tuple[float, ...] | None = None,
) -> MbcjrResult:
    """Proposed M-BCJR state pruning for the BPSK Ungerboeck model."""
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    n = y.size
    if m_states <= 0:
        raise ValueError("m_states must be positive.")
    if future_len < 0:
        raise ValueError("future_len must be non-negative.")
    if isi_len is None:
        isi_len = max((abs(int(k)) for k in g.keys()), default=0)
    if initial_state is None:
        init_state_int = (1 << isi_len) - 1  # all +1
    else:
        if len(initial_state) != isi_len:
            raise ValueError("initial_state length must equal isi_len.")
        init_state_int = 0
        for i, s in enumerate(initial_state):
            if s > 0:
                init_state_int |= (1 << (isi_len - 1 - i))
    la_arr = _prepare_la(la, n)
    g_taps = _g_dict_to_array(g, isi_len)
    future_table = _build_future_table(future_len)

    llr, surv_counts = _nb_mbcjr_kernel(
        y, g_taps, n0, la_arr, isi_len, m_states, future_len,
        future_table, init_state_int,
    )
    return MbcjrResult(llr=llr, survivor_counts=surv_counts.tolist())
