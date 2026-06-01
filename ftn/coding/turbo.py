from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


LOG_ZERO = -np.inf


@dataclass(frozen=True)
class TurboEncodeResult:
    coded_bits: np.ndarray
    n_info: int
    n_tail1: int
    n_tail2: int


@dataclass(frozen=True)
class TurboDecodeResult:
    info_llr: np.ndarray
    n_iterations: int
    converged: bool
    sys_extrinsic: np.ndarray | None = None


# ---------------------------------------------------------------------------
# RSC Encoder helpers (non-Numba, small and called once per frame)
# ---------------------------------------------------------------------------


def _rsc_encode(info_bits: np.ndarray, feedback_poly: list[int],
                fwd_poly: list[int], memory: int,
                initial_state: int = 0) -> tuple[np.ndarray, np.ndarray, int]:
    b = np.asarray(info_bits, dtype=np.uint8).reshape(-1)
    n = b.size
    state = initial_state
    parity = np.zeros(n, dtype=np.uint8)
    for k in range(n):
        fb = int(b[k])
        for j in range(memory):
            if feedback_poly[j]:
                fb ^= (state >> j) & 1
        p = fb
        for j in range(memory):
            if fwd_poly[j]:
                p ^= (state >> j) & 1
        parity[k] = p
        state = ((state << 1) | fb) & ((1 << memory) - 1)
    return b.copy(), parity, state


def rsc_encode_4state(info_bits: np.ndarray, initial_state: int = 0):
    """RSC encoder with g = [1+D+D^2] / [1+D^2]. Returns (sys, parity, final_state)."""
    # feedback poly (denominator): 1 + D^2 = [1,0,1]
    # feedforward poly (numerator): 1 + D + D^2 = [1,1,1]
    fb = [1, 0, 1]
    ff = [1, 1, 1]
    sys, parity, fin = _rsc_encode(info_bits, fb, ff, 2, initial_state)
    return sys, parity, fin


def rsc_encode_8state(info_bits: np.ndarray, initial_state: int = 0):
    """RSC encoder with g = [1+D+D^3] / [1+D^2+D^3]. Returns (sys, parity, final_state)."""
    fb = [1, 0, 1, 1]
    ff = [1, 1, 0, 1]
    sys, parity, fin = _rsc_encode(info_bits, fb, ff, 3, initial_state)
    return sys, parity, fin


def turbo_encode(info_bits: np.ndarray,
                 interleaver: np.ndarray) -> TurboEncodeResult:
    """Rate-1/3 turbo encoder.

    Output layout: [sys(K+mu1), p1(K+mu1), p2(K+mu2)]
    where mu1=2, mu2=3. Total = 3K + 2*mu1 + 2*mu2 = 3K + 10.

    For the tail sections:
    - sys[K:K+mu2] = interleaved tail2 info (zeros in original domain)
    - p1[K:K+mu1] = encoder 1 parity tail
    - p1[K+mu1:] = zeros (no encoder 1 during encoder 2 tail)
    - p2[K:K+mu1] = zeros (no encoder 2 during encoder 1 tail)
    - p2[K+mu1:] = encoder 2 parity tail
    """
    info = np.asarray(info_bits, dtype=np.uint8).reshape(-1)
    K = info.size
    pi = np.asarray(interleaver, dtype=np.int64).reshape(-1)
    mu1 = 2
    mu2 = 3
    n_tail = mu1 + mu2

    # Encoder 1: info bits + mu1 tail bits
    sys1_full, p1_full, _ = rsc_encode_4state(
        np.concatenate([info, np.zeros(mu1, dtype=np.uint8)])
    )

    # Encoder 2: interleaved info bits + mu2 tail bits
    interleaved = info[pi]
    sys2_full, p2_full, _ = rsc_encode_8state(
        np.concatenate([interleaved, np.zeros(mu2, dtype=np.uint8)])
    )

    # Assemble output: sys, p1, p2 — all padded to K + n_tail
    sys_out = np.zeros(K + n_tail, dtype=np.uint8)
    sys_out[:K] = sys1_full[:K]
    sys_out[K:K + mu1] = sys1_full[K:]
    sys_out[K + mu1:] = sys2_full[K:]

    p1_out = np.zeros(K + n_tail, dtype=np.uint8)
    p1_out[:K] = p1_full[:K]
    p1_out[K:K + mu1] = p1_full[K:]

    p2_out = np.zeros(K + n_tail, dtype=np.uint8)
    p2_out[:K] = p2_full[:K]
    p2_out[K + mu1:] = p2_full[K:]

    coded = np.concatenate([sys_out, p1_out, p2_out])
    return TurboEncodeResult(
        coded_bits=coded, n_info=K, n_tail1=mu1, n_tail2=mu2,
    )


# ---------------------------------------------------------------------------
# S-random interleaver
# ---------------------------------------------------------------------------


def s_random_interleaver(k: int, s: int | None = None,
                         seed: int = 0) -> np.ndarray:
    """Generate an S-random interleaver of length ``k``."""
    if s is None:
        s = max(1, int(np.sqrt(k / 4.0)))
    rng = np.random.default_rng(seed)
    pi = -np.ones(k, dtype=np.int64)
    available = list(range(k))
    for i in range(k):
        candidates = list(available)
        rng.shuffle(candidates)
        chosen = None
        for c in candidates:
            ok = True
            lo = max(0, i - s)
            for j in range(lo, i):
                if pi[j] >= 0 and abs(c - pi[j]) < s:
                    ok = False
                    break
            if ok:
                chosen = c
                break
        if chosen is None:
            chosen = candidates[0]
        pi[i] = chosen
        available.remove(chosen)
    return pi


# ---------------------------------------------------------------------------
# Numba helpers
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
def _nb_rsc_bcjr_4state(sys_llr: np.ndarray, parity_llr: np.ndarray,
                        init_state: int, fin_state: int):
    """Log-MAP BCJR for 4-state RSC: g=[1+D+D^2]/[1+D^2].

    Trellis: 4 states, 2 transitions each.
    State encoding: bits (s1, s0) where s0 is LSB (most recent input to shift reg).
    """
    n = sys_llr.size
    n_states = 4

    # Precompute trellis tables
    # next_state[s, b], output_parity[s, b]
    ns_tbl = np.empty((n_states, 2), dtype=np.int64)
    par_tbl = np.empty((n_states, 2), dtype=np.int64)
    for s in range(n_states):
        for b in range(2):
            # feedback poly 1+D^2: taps 0,2 — in state encoding bit0=D^2, bit1=D^1
            fb = b ^ (s & 1)
            # feedforward poly 1+D+D^2: taps 0,1,2
            p = fb ^ (s & 1) ^ ((s >> 1) & 1)
            ns_tbl[s, b] = ((s << 1) | fb) & 3
            par_tbl[s, b] = p

    # Gamma
    gamma = np.empty((n, n_states, 2))
    for k in range(n):
        sl = sys_llr[k]
        pl = parity_llr[k]
        for s in range(n_states):
            for b in range(2):
                lp_sys = -_nb_logaddexp(0.0, -sl) if b == 0 else -_nb_logaddexp(0.0, sl)
                p = par_tbl[s, b]
                lp_par = -_nb_logaddexp(0.0, -pl) if p == 0 else -_nb_logaddexp(0.0, pl)
                gamma[k, s, b] = lp_sys + lp_par

    # Alpha
    alpha = np.full((n + 1, n_states), -np.inf)
    alpha[0, init_state] = 0.0
    for k in range(n):
        for s in range(n_states):
            if np.isfinite(alpha[k, s]):
                av = alpha[k, s]
                for b in range(2):
                    ns = ns_tbl[s, b]
                    alpha[k + 1, ns] = _nb_logaddexp(alpha[k + 1, ns], av + gamma[k, s, b])
        mx = alpha[k + 1, 0]
        for s in range(1, n_states):
            if alpha[k + 1, s] > mx:
                mx = alpha[k + 1, s]
        for s in range(n_states):
            alpha[k + 1, s] -= mx

    # Beta
    beta = np.full((n + 1, n_states), -np.inf)
    if fin_state < 0:
        for s in range(n_states):
            beta[n, s] = 0.0
    else:
        beta[n, fin_state] = 0.0
    for k in range(n - 1, -1, -1):
        for s in range(n_states):
            v0 = gamma[k, s, 0] + beta[k + 1, ns_tbl[s, 0]]
            v1 = gamma[k, s, 1] + beta[k + 1, ns_tbl[s, 1]]
            beta[k, s] = _nb_logaddexp(v0, v1)
        mx = beta[k, 0]
        for s in range(1, n_states):
            if beta[k, s] > mx:
                mx = beta[k, s]
        for s in range(n_states):
            beta[k, s] -= mx

    # LLRs
    info_llr = np.empty(n)
    sys_ext = np.empty(n)
    for k in range(n):
        sum0 = -np.inf
        sum1 = -np.inf
        ext0 = -np.inf
        ext1 = -np.inf
        for s in range(n_states):
            for b in range(2):
                ns = ns_tbl[s, b]
                joint = alpha[k, s] + gamma[k, s, b] + beta[k + 1, ns]
                if b == 0:
                    sum0 = _nb_logaddexp(sum0, joint)
                else:
                    sum1 = _nb_logaddexp(sum1, joint)
                p = par_tbl[s, b]
                if p == 0:
                    ext0 = _nb_logaddexp(ext0, joint)
                else:
                    ext1 = _nb_logaddexp(ext1, joint)
        info_llr[k] = sum0 - sum1
        sys_ext[k] = ext0 - ext1

    return info_llr, sys_ext


@njit(cache=True)
def _nb_rsc_bcjr_8state(sys_llr: np.ndarray, parity_llr: np.ndarray,
                        init_state: int, fin_state: int):
    """Log-MAP BCJR for 8-state RSC: g=[1+D+D^3]/[1+D^2+D^3].

    State encoding: 3 bits (s2, s1, s0).
    """
    n = sys_llr.size
    n_states = 8

    ns_tbl = np.empty((n_states, 2), dtype=np.int64)
    par_tbl = np.empty((n_states, 2), dtype=np.int64)
    for s in range(n_states):
        for b in range(2):
            fb = b ^ (s & 1) ^ ((s >> 2) & 1)
            p = fb ^ (s & 1) ^ ((s >> 1) & 1)
            ns_tbl[s, b] = ((s << 1) | fb) & 7
            par_tbl[s, b] = p

    gamma = np.empty((n, n_states, 2))
    for k in range(n):
        sl = sys_llr[k]
        pl = parity_llr[k]
        for s in range(n_states):
            for b in range(2):
                lp_sys = -_nb_logaddexp(0.0, -sl) if b == 0 else -_nb_logaddexp(0.0, sl)
                p = par_tbl[s, b]
                lp_par = -_nb_logaddexp(0.0, -pl) if p == 0 else -_nb_logaddexp(0.0, pl)
                gamma[k, s, b] = lp_sys + lp_par

    alpha = np.full((n + 1, n_states), -np.inf)
    alpha[0, init_state] = 0.0
    for k in range(n):
        for s in range(n_states):
            if np.isfinite(alpha[k, s]):
                av = alpha[k, s]
                for b in range(2):
                    ns = ns_tbl[s, b]
                    alpha[k + 1, ns] = _nb_logaddexp(alpha[k + 1, ns], av + gamma[k, s, b])
        mx = alpha[k + 1, 0]
        for s in range(1, n_states):
            if alpha[k + 1, s] > mx:
                mx = alpha[k + 1, s]
        for s in range(n_states):
            alpha[k + 1, s] -= mx

    beta = np.full((n + 1, n_states), -np.inf)
    if fin_state < 0:
        for s in range(n_states):
            beta[n, s] = 0.0
    else:
        beta[n, fin_state] = 0.0
    for k in range(n - 1, -1, -1):
        for s in range(n_states):
            v0 = gamma[k, s, 0] + beta[k + 1, ns_tbl[s, 0]]
            v1 = gamma[k, s, 1] + beta[k + 1, ns_tbl[s, 1]]
            beta[k, s] = _nb_logaddexp(v0, v1)
        mx = beta[k, 0]
        for s in range(1, n_states):
            if beta[k, s] > mx:
                mx = beta[k, s]
        for s in range(n_states):
            beta[k, s] -= mx

    info_llr = np.empty(n)
    sys_ext = np.empty(n)
    for k in range(n):
        sum0 = -np.inf
        sum1 = -np.inf
        ext0 = -np.inf
        ext1 = -np.inf
        for s in range(n_states):
            for b in range(2):
                ns = ns_tbl[s, b]
                joint = alpha[k, s] + gamma[k, s, b] + beta[k + 1, ns]
                if b == 0:
                    sum0 = _nb_logaddexp(sum0, joint)
                else:
                    sum1 = _nb_logaddexp(sum1, joint)
                p = par_tbl[s, b]
                if p == 0:
                    ext0 = _nb_logaddexp(ext0, joint)
                else:
                    ext1 = _nb_logaddexp(ext1, joint)
        info_llr[k] = sum0 - sum1
        sys_ext[k] = ext0 - ext1

    return info_llr, sys_ext


# ---------------------------------------------------------------------------
# Turbo decode
# ---------------------------------------------------------------------------


def turbo_decode(
    sys_llr: np.ndarray,
    p1_llr: np.ndarray,
    p2_llr: np.ndarray,
    interleaver: np.ndarray,
    max_iterations: int = 50,
) -> TurboDecodeResult:
    """Iterative turbo decoder with early termination.

    Input LLR layout matches turbo_encode output:
    - sys_llr[K+mu1+mu2], p1_llr[K+mu1+mu2], p2_llr[K+mu1+mu2]
    - sys: [sys(K), sys_tail1(mu1), sys_tail2(mu2)]
    - p1:  [p1(K), p1_tail1(mu1), zeros(mu2)]
    - p2:  [p2(K), zeros(mu1), p2_tail2(mu2)]
    """
    sl = np.asarray(sys_llr, dtype=np.float64).reshape(-1)
    p1l = np.asarray(p1_llr, dtype=np.float64).reshape(-1)
    p2l = np.asarray(p2_llr, dtype=np.float64).reshape(-1)
    pi = np.asarray(interleaver, dtype=np.int64).reshape(-1)
    K = pi.size
    mu1 = 2
    mu2 = 3

    pi_inv = np.empty(K, dtype=np.int64)
    pi_inv[pi] = np.arange(K, dtype=np.int64)

    # RSC1 frame: K info + mu1 tail
    rsc1_sys = np.empty(K + mu1)
    rsc1_sys[:K] = sl[:K]
    rsc1_sys[K:] = sl[K:K + mu1]
    rsc1_par = np.empty(K + mu1)
    rsc1_par[:K] = p1l[:K]
    rsc1_par[K:] = p1l[K:K + mu1]

    # RSC2 frame: K interleaved info + mu2 tail
    rsc2_sys = np.empty(K + mu2)
    rsc2_sys[:K] = sl[pi]
    rsc2_sys[K:] = sl[K + mu1:K + mu1 + mu2]
    rsc2_par = np.empty(K + mu2)
    rsc2_par[:K] = p2l[:K]
    rsc2_par[K:] = p2l[K + mu1:K + mu1 + mu2]

    extrinsic_in = np.zeros(K, dtype=np.float64)
    converged = False
    last_iter = max_iterations
    info1 = np.zeros(K + mu1)
    info2 = np.zeros(K + mu2)

    for it in range(max_iterations):
        # RSC1
        rsc1_input = rsc1_sys.copy()
        rsc1_input[:K] += extrinsic_in
        info1, _ = _nb_rsc_bcjr_4state(rsc1_input, rsc1_par, 0, -1)
        ext1 = np.clip(info1[:K] - rsc1_input[:K], -30.0, 30.0)

        # RSC2
        rsc2_input = rsc2_sys.copy()
        rsc2_input[:K] += ext1[pi]
        info2, _ = _nb_rsc_bcjr_8state(rsc2_input, rsc2_par, 0, -1)
        ext2 = np.clip(info2[:K] - rsc2_input[:K], -30.0, 30.0)

        # De-interleave
        extrinsic_in = ext2[pi_inv]

        # Early termination
        hard1 = (info1[:K] < 0).astype(np.uint8)
        hard2 = (info2[pi_inv] < 0).astype(np.uint8)
        if np.array_equal(hard1, hard2):
            converged = True
            last_iter = it + 1
            break
        last_iter = it + 1

    final_info_llr = info2[pi_inv]
    return TurboDecodeResult(
        info_llr=final_info_llr,
        n_iterations=last_iter,
        converged=converged,
        sys_extrinsic=extrinsic_in.copy(),
    )
