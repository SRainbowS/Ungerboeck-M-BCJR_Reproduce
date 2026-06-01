from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


LOG_ZERO = -np.inf


@dataclass(frozen=True)
class ConvDecodeResult:
    info_llr: np.ndarray
    code_llr: np.ndarray


def _parity(value: int) -> int:
    return int(value.bit_count() & 1)


def _next_state_and_outputs(state: int, bit: int) -> tuple[int, tuple[int, int]]:
    """Four-state nonrecursive convolutional code with generators 7 and 5."""
    if state < 0 or state > 3:
        raise ValueError("state must be in [0, 3].")
    bit = int(bit)
    if bit not in (0, 1):
        raise ValueError("bit must be 0 or 1.")

    prev1 = (state >> 1) & 1
    prev2 = state & 1
    reg = (bit << 2) | (prev1 << 1) | prev2
    out_7 = _parity(reg & 0b111)
    out_5 = _parity(reg & 0b101)
    next_state = ((bit << 1) | prev1) & 0b11
    return next_state, (out_7, out_5)


def conv_encode_75(bits: np.ndarray, initial_state: int = 0) -> np.ndarray:
    """Encode bits with the nonrecursive 4-state ``(7,5)`` convolutional code."""
    bit_array = np.asarray(bits, dtype=np.uint8).reshape(-1)
    state = int(initial_state)
    encoded = np.zeros(2 * bit_array.size, dtype=np.uint8)
    for idx, bit in enumerate(bit_array):
        state, outputs = _next_state_and_outputs(state, int(bit))
        encoded[2 * idx : 2 * idx + 2] = outputs
    return encoded


def _log_prob_from_llr(bit: int, llr_zero_over_one: float) -> float:
    if bit == 0:
        return float(-np.logaddexp(0.0, -llr_zero_over_one))
    return float(-np.logaddexp(0.0, llr_zero_over_one))


# ---------------------------------------------------------------------------
# Numba-optimised BCJR kernel
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
def _nb_conv_bcjr_kernel(
    code_llr: np.ndarray,
    init_state: int,
    fin_state: int,
) -> tuple:
    n = code_llr.size // 2

    # ---- precompute trellis tables ----
    ns_tbl = np.empty((4, 2), dtype=np.int64)
    out0_tbl = np.empty((4, 2), dtype=np.int64)
    out1_tbl = np.empty((4, 2), dtype=np.int64)
    for s in range(4):
        for b in range(2):
            p1 = (s >> 1) & 1
            p2 = s & 1
            reg = (b << 2) | (p1 << 1) | p2
            r7 = reg & 7;  p7 = (r7 ^ (r7 >> 1) ^ (r7 >> 2)) & 1
            r5 = reg & 5;  p5 = (r5 ^ (r5 >> 2)) & 1
            ns_tbl[s, b] = ((b << 1) | p1) & 3
            out0_tbl[s, b] = p7
            out1_tbl[s, b] = p5

    # ---- branch metrics (gamma) ----
    gamma = np.empty((n, 4, 2))
    for idx in range(n):
        llr0 = code_llr[2 * idx]
        llr1 = code_llr[2 * idx + 1]
        for s in range(4):
            for b in range(2):
                o0 = out0_tbl[s, b]
                o1 = out1_tbl[s, b]
                lp0 = -_nb_logaddexp(0.0, -llr0) if o0 == 0 else -_nb_logaddexp(0.0, llr0)
                lp1 = -_nb_logaddexp(0.0, -llr1) if o1 == 0 else -_nb_logaddexp(0.0, llr1)
                gamma[idx, s, b] = lp0 + lp1

    # ---- forward recursion (alpha) ----
    alpha = np.full((n + 1, 4), -np.inf)
    alpha[0, init_state] = 0.0
    for idx in range(n):
        for s in range(4):
            if np.isfinite(alpha[idx, s]):
                a_val = alpha[idx, s]
                for b in range(2):
                    ns = ns_tbl[s, b]
                    alpha[idx + 1, ns] = _nb_logaddexp(
                        alpha[idx + 1, ns], a_val + gamma[idx, s, b]
                    )
        mx = alpha[idx + 1, 0]
        for s in range(1, 4):
            if alpha[idx + 1, s] > mx:
                mx = alpha[idx + 1, s]
        for s in range(4):
            alpha[idx + 1, s] -= mx

    # ---- backward recursion (beta) ----
    beta = np.full((n + 1, 4), -np.inf)
    if fin_state < 0:
        for s in range(4):
            beta[n, s] = 0.0
    else:
        beta[n, fin_state] = 0.0
    for idx in range(n - 1, -1, -1):
        for s in range(4):
            v0 = gamma[idx, s, 0] + beta[idx + 1, ns_tbl[s, 0]]
            v1 = gamma[idx, s, 1] + beta[idx + 1, ns_tbl[s, 1]]
            beta[idx, s] = _nb_logaddexp(v0, v1)
        mx = beta[idx, 0]
        for s in range(1, 4):
            if beta[idx, s] > mx:
                mx = beta[idx, s]
        for s in range(4):
            beta[idx, s] -= mx

    # ---- a-posteriori LLRs ----
    info_llr = np.empty(n)
    flat_code = np.empty(2 * n)
    for idx in range(n):
        sum_b0 = -np.inf
        sum_b1 = -np.inf
        p0_0 = -np.inf;  p0_1 = -np.inf
        p1_0 = -np.inf;  p1_1 = -np.inf
        for s in range(4):
            for b in range(2):
                ns = ns_tbl[s, b]
                joint = alpha[idx, s] + gamma[idx, s, b] + beta[idx + 1, ns]
                if b == 0:
                    sum_b0 = _nb_logaddexp(sum_b0, joint)
                else:
                    sum_b1 = _nb_logaddexp(sum_b1, joint)
                o0 = out0_tbl[s, b]
                o1 = out1_tbl[s, b]
                if o0 == 0:
                    p0_0 = _nb_logaddexp(p0_0, joint)
                else:
                    p0_1 = _nb_logaddexp(p0_1, joint)
                if o1 == 0:
                    p1_0 = _nb_logaddexp(p1_0, joint)
                else:
                    p1_1 = _nb_logaddexp(p1_1, joint)
        info_llr[idx] = sum_b0 - sum_b1
        flat_code[2 * idx] = p0_0 - p0_1
        flat_code[2 * idx + 1] = p1_0 - p1_1

    return info_llr, flat_code


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def conv_bcjr_decode(
    code_llr: np.ndarray,
    initial_state: int = 0,
    final_state_known: int | None = None,
) -> ConvDecodeResult:
    """SISO log-MAP decoder for the ``(7,5)`` convolutional code.

    ``code_llr`` uses the convention ``log P(c=0) / P(c=1)``.
    """
    llr = np.asarray(code_llr, dtype=np.float64).reshape(-1)
    if llr.size % 2 != 0:
        raise ValueError("code_llr length must be even.")
    fin = -1 if final_state_known is None else int(final_state_known)
    info_llr, flat_code = _nb_conv_bcjr_kernel(llr, int(initial_state), fin)
    return ConvDecodeResult(info_llr=info_llr, code_llr=flat_code)
