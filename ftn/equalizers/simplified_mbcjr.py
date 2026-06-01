from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


LOG_ZERO = -np.inf


@dataclass(frozen=True)
class SimplifiedMbcjrResult:
    bit_llr: np.ndarray
    survivor_counts: list[int]


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
def _nb_log_gamma_bpsk(y_n: float, x_n: float, state_int: int,
                       g_taps: np.ndarray, isi_len: int,
                       inv_n0: float, la_n: float) -> float:
    isi = 0.0
    taps_len = len(g_taps)
    for lag in range(1, isi_len + 1):
        if lag < taps_len:
            sym = 1.0 if (state_int >> (lag - 1)) & 1 else -1.0
            isi += g_taps[lag] * sym
    g0 = g_taps[0] if taps_len > 0 else 0.0
    centered = y_n - 0.5 * g0 * x_n - isi
    phi = 2.0 * inv_n0 * x_n * centered
    if x_n > 0:
        prior = -_nb_logaddexp(0.0, -la_n)
    else:
        prior = -_nb_logaddexp(0.0, la_n)
    return phi + prior


# ---------------------------------------------------------------------------
# Key tail path search — BPSK
# ---------------------------------------------------------------------------


@njit(cache=True)
def _nb_find_key_tail_bpsk(y: np.ndarray, g_taps: np.ndarray, isi_len: int,
                           inv_n0: float, la: np.ndarray, state_int: int,
                           idx: int, future_len: int,
                           future_table: np.ndarray) -> tuple:
    """Find the best future path from *state_int* at position *idx*.

    Returns (best_metric, best_path_int).
    """
    n = y.size
    remaining = n - idx - 1
    horizon = future_len if future_len < remaining else remaining
    if horizon <= 0:
        return 0.0, 0

    n_fut = 1 << horizon
    mask = (1 << isi_len) - 1
    best_metric = -np.inf
    best_path = 0

    for fi in range(n_fut):
        st = state_int
        metric = 0.0
        for off in range(horizon):
            symbol = future_table[fi, off]
            j = idx + off + 1
            metric += _nb_log_gamma_bpsk(y[j], symbol, st, g_taps, isi_len,
                                         inv_n0, la[j])
            bit = 1 if symbol > 0.0 else 0
            st = ((st << 1) | bit) & mask
        if metric > best_metric:
            best_metric = metric
            best_path = fi

    return best_metric, best_path


@njit(cache=True)
def _nb_beta_single_path_bpsk(y: np.ndarray, g_taps: np.ndarray, isi_len: int,
                              inv_n0: float, la: np.ndarray, state_int: int,
                              idx: int, future_len: int,
                              key_path: np.ndarray) -> float:
    """Compute beta using the single *key_path* instead of enumerating all."""
    n = y.size
    remaining = n - idx - 1
    horizon = future_len if future_len < remaining else remaining
    if horizon <= 0:
        return 0.0

    mask = (1 << isi_len) - 1
    st = state_int
    metric = 0.0
    for off in range(horizon):
        symbol = key_path[off]
        j = idx + off + 1
        metric += _nb_log_gamma_bpsk(y[j], symbol, st, g_taps, isi_len,
                                     inv_n0, la[j])
        bit = 1 if symbol > 0.0 else 0
        st = ((st << 1) | bit) & mask
    return metric


@njit(cache=True)
def _nb_simplified_mbcjr_kernel_bpsk(
    y: np.ndarray, g_taps: np.ndarray, n0: float,
    la: np.ndarray, isi_len: int, m_states: int,
    future_len: int, future_table: np.ndarray,
    init_state_int: int,
) -> tuple:
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
                g_val = _nb_log_gamma_bpsk(y[idx], sym, sp, g_taps, isi_len,
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

        # --- find reference state (best posterior) ---
        # Use the first candidate as reference to find key tail path
        ref_idx = 0
        ref_metric = -np.inf
        for ci in range(n_cand):
            # Quick estimate: alpha metric alone
            if cand_m[ci] > ref_metric:
                ref_metric = cand_m[ci]
                ref_idx = ci

        # Find key tail from reference state
        _, key_fi = _nb_find_key_tail_bpsk(
            y, g_taps, isi_len, inv_n0, la, cand_s[ref_idx],
            idx, future_len, future_table,
        )

        # Extract key path symbols
        remaining = n - idx - 1
        horizon = future_len if future_len < remaining else remaining
        key_path = np.empty(horizon, dtype=np.float64)
        for off in range(horizon):
            key_path[off] = future_table[key_fi, off]

        # --- beta via single key path ---
        for ci in range(n_cand):
            beta_val = _nb_beta_single_path_bpsk(
                y, g_taps, isi_len, inv_n0, la, cand_s[ci],
                idx, future_len, key_path,
            )
            lpost[ci] = cand_m[ci] + beta_val

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


def simplified_mbcjr_bpsk(
    y: np.ndarray,
    g: dict[int, float],
    n0: float,
    la: np.ndarray | None = None,
    isi_len: int | None = None,
    m_states: int = 4,
    future_len: int = 3,
    initial_state: tuple[float, ...] | None = None,
) -> SimplifiedMbcjrResult:
    """Simplified M-BCJR with key tail path search for BPSK."""
    from ftn.equalizers.full_bcjr import _prepare_la

    y = np.asarray(y, dtype=np.float64).reshape(-1)
    n = y.size
    if m_states <= 0:
        raise ValueError("m_states must be positive.")
    if future_len < 0:
        raise ValueError("future_len must be non-negative.")
    if isi_len is None:
        isi_len = max((abs(int(k)) for k in g.keys()), default=0)
    if initial_state is None:
        init_state_int = (1 << isi_len) - 1
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

    llr, surv_counts = _nb_simplified_mbcjr_kernel_bpsk(
        y, g_taps, n0, la_arr, isi_len, m_states, future_len,
        future_table, init_state_int,
    )
    return SimplifiedMbcjrResult(bit_llr=llr, survivor_counts=surv_counts.tolist())


# ---------------------------------------------------------------------------
# 16-QAM simplified M-BCJR
# ---------------------------------------------------------------------------


@njit(cache=True)
def _nb_log_gamma_qam16(
    y_re: float, y_im: float,
    sym_re: float, sym_im: float,
    state_int: int,
    g_taps: np.ndarray, isi_len: int,
    inv_n0: float,
    constellation_re: np.ndarray, constellation_im: np.ndarray,
    bit_map: np.ndarray, la: np.ndarray, la_offset: int,
) -> float:
    """Ungerboeck branch metric for 16-QAM complex symbol."""
    taps_len = len(g_taps)
    # ISI from past symbols stored in state
    isi_re = 0.0
    isi_im = 0.0
    base = 16
    for lag in range(1, isi_len + 1):
        if lag < taps_len:
            sym_idx = (state_int >> (4 * (lag - 1))) & 0xF
            isi_re += g_taps[lag] * constellation_re[sym_idx]
            isi_im += g_taps[lag] * constellation_im[sym_idx]
    g0 = g_taps[0] if taps_len > 0 else 0.0
    # centered = y - 0.5*g0*x - isi
    cent_re = y_re - 0.5 * g0 * sym_re - isi_re
    cent_im = y_im - 0.5 * g0 * sym_im - isi_im
    # phi = (2/N0) * Re(conj(x) * centered)
    phi = 2.0 * inv_n0 * (sym_re * cent_re + sym_im * cent_im)
    # Prior from bit LLRs
    sym_idx_local = -1
    for si in range(16):
        if constellation_re[si] == sym_re and constellation_im[si] == sym_im:
            sym_idx_local = si
            break
    if sym_idx_local < 0:
        return -np.inf
    prior = 0.0
    for m in range(4):
        bit_val = bit_map[sym_idx_local, m]
        llr_m = la[la_offset + m]
        if bit_val == 0:
            prior += -_nb_logaddexp(0.0, -llr_m)
        else:
            prior += -_nb_logaddexp(0.0, llr_m)
    return phi + prior


@njit(cache=True)
def _nb_find_key_tail_qam16(
    y_re: np.ndarray, y_im: np.ndarray,
    g_taps: np.ndarray, isi_len: int,
    inv_n0: float,
    constellation_re: np.ndarray, constellation_im: np.ndarray,
    bit_map: np.ndarray, la: np.ndarray,
    state_int: int, idx: int, future_len: int,
    future_table_idx: np.ndarray,
) -> tuple:
    n = y_re.size
    remaining = n - idx - 1
    horizon = future_len if future_len < remaining else remaining
    if horizon <= 0:
        return 0.0, 0

    n_fut = future_table_idx.shape[0]
    best_metric = -np.inf
    best_fi = 0

    for fi in range(n_fut):
        st = state_int
        metric = 0.0
        for off in range(horizon):
            si = future_table_idx[fi, off]
            j = idx + off + 1
            metric += _nb_log_gamma_qam16(
                y_re[j], y_im[j],
                constellation_re[si], constellation_im[si],
                st, g_taps, isi_len, inv_n0,
                constellation_re, constellation_im, bit_map,
                la, 4 * j,
            )
            st = ((st << 4) | si) & ((1 << (4 * isi_len)) - 1)
        if metric > best_metric:
            best_metric = metric
            best_fi = fi

    return best_metric, best_fi


@njit(cache=True)
def _nb_beta_single_path_qam16(
    y_re: np.ndarray, y_im: np.ndarray,
    g_taps: np.ndarray, isi_len: int,
    inv_n0: float,
    constellation_re: np.ndarray, constellation_im: np.ndarray,
    bit_map: np.ndarray, la: np.ndarray,
    state_int: int, idx: int, future_len: int,
    key_path_idx: np.ndarray,
) -> float:
    n = y_re.size
    remaining = n - idx - 1
    horizon = future_len if future_len < remaining else remaining
    if horizon <= 0:
        return 0.0

    st = state_int
    metric = 0.0
    for off in range(horizon):
        si = key_path_idx[off]
        j = idx + off + 1
        metric += _nb_log_gamma_qam16(
            y_re[j], y_im[j],
            constellation_re[si], constellation_im[si],
            st, g_taps, isi_len, inv_n0,
            constellation_re, constellation_im, bit_map,
            la, 4 * j,
        )
        st = ((st << 4) | si) & ((1 << (4 * isi_len)) - 1)
    return metric


@njit(cache=True)
def _nb_simplified_mbcjr_kernel_qam16(
    y_re: np.ndarray, y_im: np.ndarray,
    g_taps: np.ndarray, n0: float,
    la: np.ndarray, isi_len: int, m_states: int,
    future_len: int,
    future_table_idx: np.ndarray,
    constellation_re: np.ndarray, constellation_im: np.ndarray,
    bit_map: np.ndarray,
    init_state_int: int,
) -> tuple:
    n = y_re.size
    state_mask = (1 << (4 * isi_len)) - 1
    inv_n0 = 1.0 / n0

    max_cand = 16 * m_states + 16

    surv_s = np.empty(m_states, dtype=np.int64)
    surv_m = np.empty(m_states, dtype=np.float64)
    surv_s[0] = init_state_int
    surv_m[0] = 0.0
    n_surv = 1

    n_bits = n * 4
    bit_llr = np.zeros(n_bits, dtype=np.float64)
    surv_counts = np.zeros(n, dtype=np.int64)

    cand_s = np.empty(max_cand, dtype=np.int64)
    cand_m = np.empty(max_cand, dtype=np.float64)
    cand_sym = np.empty(max_cand, dtype=np.int64)
    lpost = np.empty(max_cand, dtype=np.float64)

    # Per-bit accumulator arrays (separate for bit=0 and bit=1)
    bit0_vals = np.empty(max_cand, dtype=np.float64)
    bit1_vals = np.empty(max_cand, dtype=np.float64)

    for idx in range(n):
        la_offset = 4 * idx

        # --- expand candidates ---
        n_cand = 0
        for si in range(n_surv):
            sp = surv_s[si]
            am = surv_m[si]
            for sym_idx in range(16):
                sn = ((sp << 4) | sym_idx) & state_mask
                g_val = _nb_log_gamma_qam16(
                    y_re[idx], y_im[idx],
                    constellation_re[sym_idx], constellation_im[sym_idx],
                    sp, g_taps, isi_len, inv_n0,
                    constellation_re, constellation_im, bit_map,
                    la, la_offset,
                )
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
                    cand_sym[n_cand] = sym_idx
                    n_cand += 1

        # --- find key tail path ---
        ref_idx = 0
        ref_metric = -np.inf
        for ci in range(n_cand):
            if cand_m[ci] > ref_metric:
                ref_metric = cand_m[ci]
                ref_idx = ci

        _, key_fi = _nb_find_key_tail_qam16(
            y_re, y_im, g_taps, isi_len, inv_n0,
            constellation_re, constellation_im, bit_map, la,
            cand_s[ref_idx], idx, future_len, future_table_idx,
        )

        remaining = n - idx - 1
        horizon = future_len if future_len < remaining else remaining
        key_path_idx = np.empty(horizon, dtype=np.int64)
        for off in range(horizon):
            key_path_idx[off] = future_table_idx[key_fi, off]

        # --- beta ---
        for ci in range(n_cand):
            beta_val = _nb_beta_single_path_qam16(
                y_re, y_im, g_taps, isi_len, inv_n0,
                constellation_re, constellation_im, bit_map, la,
                cand_s[ci], idx, future_len, key_path_idx,
            )
            lpost[ci] = cand_m[ci] + beta_val

        # --- bit LLR ---
        for bit_pos in range(4):
            n_b0 = 0
            n_b1 = 0
            for ci in range(n_cand):
                si = cand_sym[ci]
                if bit_map[si, bit_pos] == 0:
                    bit0_vals[n_b0] = lpost[ci]
                    n_b0 += 1
                else:
                    bit1_vals[n_b1] = lpost[ci]
                    n_b1 += 1
            s0 = _nb_logsumexp(bit0_vals, n_b0)
            s1 = _nb_logsumexp(bit1_vals, n_b1)
            bit_llr[4 * idx + bit_pos] = s0 - s1

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

        norm = _nb_logsumexp(surv_m, n_surv)
        for si in range(n_surv):
            surv_m[si] -= norm

        surv_counts[idx] = n_surv

    return bit_llr, surv_counts


def _build_future_table_qam16(future_len: int) -> np.ndarray:
    if future_len <= 0:
        return np.empty((1, 0), dtype=np.int64)
    n_seq = 16 ** future_len
    table = np.empty((n_seq, future_len), dtype=np.int64)
    for i in range(n_seq):
        val = i
        for j in range(future_len):
            table[i, j] = val % 16
            val //= 16
    return table


def simplified_mbcjr_qam16(
    y_complex: np.ndarray,
    g: dict[int, float],
    n0: float,
    la: np.ndarray | None = None,
    isi_len: int = 3,
    m_states: int = 4,
    future_len: int = 3,
) -> SimplifiedMbcjrResult:
    """Simplified M-BCJR with key tail path search for 16-QAM."""
    from ftn.equalizers.full_bcjr import _prepare_la
    from ftn.modulation import qam16_constellation

    y = np.asarray(y_complex, dtype=complex).reshape(-1)
    n = y.size
    n_bits = n * 4
    if la is None:
        la_arr = np.zeros(n_bits, dtype=np.float64)
    else:
        la_arr = np.asarray(la, dtype=np.float64).reshape(-1)
        if la_arr.size != n_bits:
            raise ValueError(f"la size {la_arr.size} != expected {n_bits}")

    symbols, bit_map = qam16_constellation()
    constellation_re = symbols.real.astype(np.float64)
    constellation_im = symbols.imag.astype(np.float64)
    bit_map_arr = bit_map.astype(np.uint8)

    g_taps = _g_dict_to_array(g, isi_len)
    future_table_idx = _build_future_table_qam16(future_len)

    init_state_int = 0

    bit_llr, surv_counts = _nb_simplified_mbcjr_kernel_qam16(
        y.real.copy(), y.imag.copy(),
        g_taps, n0, la_arr, isi_len, m_states, future_len,
        future_table_idx, constellation_re, constellation_im, bit_map_arr,
        init_state_int,
    )
    return SimplifiedMbcjrResult(bit_llr=bit_llr, survivor_counts=surv_counts.tolist())
