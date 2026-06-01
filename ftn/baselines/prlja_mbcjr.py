"""Paper [14] Prlja & Anderson M-BCJR variants for Forney model FTN signaling.

Implements the Simple M-BCJR and Backup M-BCJR algorithms from:
"Reduced-Complexity Receivers for Strongly Narrowband Intersymbol Interference
Introduced by Faster-than-Nyquist Signaling" by Prlja & Anderson.

The algorithms operate on the Forney (whitened matched filter) model where
observations are y[n] = sum_k v[k]*x[n-k] + w[n] with white Gaussian noise.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PrljaMbcjrResult:
    """Result from a Prlja M-BCJR equalizer run."""
    llr: np.ndarray
    survivor_counts: list[int]


LLR_RESERVE = 5.0


# ---------------------------------------------------------------------------
# Numba JIT helpers
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
def _nb_log_prior(x_n: float, la_n: float) -> float:
    """Log prior for BPSK symbol given LLR la_n."""
    if x_n > 0:
        return -_nb_logaddexp(0.0, -la_n)
    else:
        return -_nb_logaddexp(0.0, la_n)


@njit(cache=True)
def _nb_compute_label(x_n: float, state_int: int, v: np.ndarray, L: int) -> float:
    """Compute expected channel output label = v[0]*x_n + sum_{k=1}^L v[k]*x_{n-k}.

    State encoding: bit i of state_int encodes symbol at position i
    (bit=0 -> -1, bit=1 -> +1). Bit 0 is the newest symbol in the state.
    """
    label = v[0] * x_n
    for k in range(1, L + 1):
        sym = 1.0 if (state_int >> (k - 1)) & 1 else -1.0
        label += v[k] * sym
    return label


@njit(cache=True)
def _nb_branch_metric(y_n: float, label: float, n0: float, log_prior: float) -> float:
    """Forney branch metric: -(y_n - label)^2 / N0 + log_prior."""
    diff = y_n - label
    return -(diff * diff) / n0 + log_prior


# ---------------------------------------------------------------------------
# Simple M-BCJR kernel (Paper [14])
# ---------------------------------------------------------------------------

@njit(cache=True)
def _nb_prlja_simple_kernel(
    y: np.ndarray,
    v: np.ndarray,
    n0: float,
    la: np.ndarray,
    L: int,
    M: int,
) -> tuple:
    """Simple Detection M-BCJR for BPSK Forney model.

    Forward pass: M-algorithm with merge.
    Backward pass: M-algorithm with overlap-first priority.
    LLR: combine alpha + beta for each surviving state.
    """
    N = y.size
    mask = (1 << L) - 1

    max_cand = 2 * M + 2

    # --- Forward pass ---
    fwd_states = np.empty((N + 1, M), dtype=np.int64)
    fwd_metrics = np.empty((N + 1, M), dtype=np.float64)
    fwd_counts = np.empty(N + 1, dtype=np.int64)

    # Initialise: state 0 (all -1 bits), metric 0
    fwd_states[0, 0] = 0
    fwd_metrics[0, 0] = 0.0
    fwd_counts[0] = 1

    cand_s = np.empty(max_cand, dtype=np.int64)
    cand_m = np.empty(max_cand, dtype=np.float64)

    for idx in range(N):
        n_fwd = fwd_counts[idx]
        n_cand = 0

        for si in range(n_fwd):
            sp = fwd_states[idx, si]
            am = fwd_metrics[idx, si]
            for sym_idx in range(2):
                sym = -1.0 if sym_idx == 0 else 1.0
                sn = ((sp << 1) | sym_idx) & mask

                label = _nb_compute_label(sym, sp, v, L)
                lp = _nb_log_prior(sym, la[idx])
                bm = _nb_branch_metric(y[idx], label, n0, lp)
                metric = am + bm

                # Merge: check if destination state already exists
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

        # Keep top-M by metric value
        keep = M if M < n_cand else n_cand
        if keep < n_cand:
            order = np.argsort(cand_m[:n_cand])[::-1]  # descending
        else:
            order = np.arange(n_cand)

        for ki in range(keep):
            fwd_states[idx + 1, ki] = cand_s[order[ki]]
            fwd_metrics[idx + 1, ki] = cand_m[order[ki]]
        fwd_counts[idx + 1] = keep

        # Normalise
        norm = _nb_logsumexp(fwd_metrics[idx + 1, :keep], keep)
        if np.isfinite(norm):
            for si in range(keep):
                fwd_metrics[idx + 1, si] -= norm
        else:
            for si in range(keep):
                fwd_metrics[idx + 1, si] = 0.0

    # --- Backward pass with overlap-first priority ---
    bwd_states = np.empty((N + 1, M), dtype=np.int64)
    bwd_metrics = np.empty((N + 1, M), dtype=np.float64)
    bwd_counts = np.empty(N + 1, dtype=np.int64)

    # Initialise backward: seed with all forward alpha states at time N
    n_alpha_final = fwd_counts[N]
    n_seed = min(n_alpha_final, M)
    for si in range(n_seed):
        bwd_states[N, si] = fwd_states[N, si]
        bwd_metrics[N, si] = 0.0  # uniform terminal metric
    bwd_counts[N] = n_seed

    for idx in range(N - 1, -1, -1):
        n_bwd = bwd_counts[idx + 1]
        n_cand = 0

        # Expand backward: from beta state s' at time idx+1,
        # find source states s at time idx.
        # Forward transition: s' = ((s << 1) | bit) & mask
        # So s bits 0..L-2 = s' bits 1..L-1, and s bit L-1 is free.
        # For each s' and each value of s bit (L-1) (0 or 1):
        #   s = (s' >> 1) | (b << (L-1))
        # The symbol at time idx+1 is bit 0 of s' (the newest bit shifted in).
        for si in range(n_bwd):
            sp = bwd_states[idx + 1, si]
            bm_val = bwd_metrics[idx + 1, si]
            sym_bit = sp & 1  # bit 0 = symbol at time idx+1
            sym = 1.0 if sym_bit == 1 else -1.0

            for old_msb in range(2):
                # Reconstruct source state
                if L > 0:
                    s_src = (sp >> 1) | (old_msb << (L - 1))
                else:
                    s_src = sp

                # Branch metric for the transition s_src -> sp with symbol sym
                label = _nb_compute_label(sym, s_src, v, L)
                lp = _nb_log_prior(sym, la[idx])
                bm_branch = _nb_branch_metric(y[idx], label, n0, lp)
                metric = bm_branch + bm_val

                # Merge: check if source state already exists
                found = False
                for ci in range(n_cand):
                    if cand_s[ci] == s_src:
                        cand_m[ci] = _nb_logaddexp(cand_m[ci], metric)
                        found = True
                        break
                if not found:
                    cand_s[n_cand] = s_src
                    cand_m[n_cand] = metric
                    n_cand += 1

        # Overlap-first priority
        n_alpha = fwd_counts[idx]
        alpha_state_set = set()
        for si in range(n_alpha):
            alpha_state_set.add(int(fwd_states[idx, si]))

        # Separate overlapping and non-overlapping candidates
        overlap_s = np.empty(n_cand, dtype=np.int64)
        overlap_m = np.empty(n_cand, dtype=np.float64)
        other_s = np.empty(n_cand, dtype=np.int64)
        other_m = np.empty(n_cand, dtype=np.float64)
        n_overlap = 0
        n_other = 0

        for ci in range(n_cand):
            if int(cand_s[ci]) in alpha_state_set:
                overlap_s[n_overlap] = cand_s[ci]
                overlap_m[n_overlap] = cand_m[ci]
                n_overlap += 1
            else:
                other_s[n_other] = cand_s[ci]
                other_m[n_other] = cand_m[ci]
                n_other += 1

        # Sort overlapping by metric (descending)
        if n_overlap > 1:
            order_o = np.argsort(overlap_m[:n_overlap])[::-1]
        else:
            order_o = np.arange(n_overlap)

        # Sort non-overlapping by metric (descending)
        if n_other > 1:
            order_n = np.argsort(other_m[:n_other])[::-1]
        else:
            order_n = np.arange(n_other)

        # Fill M slots: overlap first, then non-overlap
        n_keep = 0
        for ki in range(min(n_overlap, M)):
            bwd_states[idx, n_keep] = overlap_s[order_o[ki]]
            bwd_metrics[idx, n_keep] = overlap_m[order_o[ki]]
            n_keep += 1

        remaining = M - n_keep
        for ki in range(min(n_other, remaining)):
            bwd_states[idx, n_keep] = other_s[order_n[ki]]
            bwd_metrics[idx, n_keep] = other_m[order_n[ki]]
            n_keep += 1

        bwd_counts[idx] = n_keep

        # Normalise
        if n_keep > 0:
            temp = np.empty(n_keep, dtype=np.float64)
            for si in range(n_keep):
                temp[si] = bwd_metrics[idx, si]
            norm = _nb_logsumexp(temp, n_keep)
            if np.isfinite(norm):
                for si in range(n_keep):
                    bwd_metrics[idx, si] -= norm
            else:
                for si in range(n_keep):
                    bwd_metrics[idx, si] = 0.0

    # --- LLR computation ---
    LAMBDA = LLR_RESERVE  # reserve value
    llr = np.zeros(N, dtype=np.float64)
    surv_counts = np.zeros(N, dtype=np.int64)

    plus_v = np.empty(M * 2, dtype=np.float64)
    minus_v = np.empty(M * 2, dtype=np.float64)

    for idx in range(N):
        n_alpha = fwd_counts[idx + 1]  # alpha after processing idx
        n_beta = bwd_counts[idx + 1]   # beta at time idx+1

        # The state at time idx+1 has bit 0 = x_idx (newest symbol).
        # Combine alpha[idx+1] and beta[idx+1] for each surviving state.

        n_plus = 0
        n_minus = 0

        for ai in range(n_alpha):
            sa = fwd_states[idx + 1, ai]
            am = fwd_metrics[idx + 1, ai]

            # Find matching beta
            beta_val = -np.inf
            for bi in range(n_beta):
                if bwd_states[idx + 1, bi] == sa:
                    beta_val = bwd_metrics[idx + 1, bi]
                    break

            posterior = am + beta_val

            # bit 0 of state = symbol at time idx
            if (sa & 1) == 1:
                plus_v[n_plus] = posterior
                n_plus += 1
            else:
                minus_v[n_minus] = posterior
                n_minus += 1

        if n_plus > 0 and n_minus > 0:
            p_val = _nb_logsumexp(plus_v, n_plus)
            m_val = _nb_logsumexp(minus_v, n_minus)
            raw = p_val - m_val
            # -inf - (-inf) = NaN when fwd/bwd survivors are disjoint
            llr[idx] = 0.0 if np.isnan(raw) else raw
        elif n_plus > 0:
            llr[idx] = LAMBDA
        elif n_minus > 0:
            llr[idx] = -LAMBDA
        else:
            llr[idx] = 0.0

        surv_counts[idx] = n_alpha

    return llr, surv_counts


# ---------------------------------------------------------------------------
# Backup M-BCJR kernel
# ---------------------------------------------------------------------------

@njit(cache=True)
def _nb_backup_forward(
    y: np.ndarray,
    v: np.ndarray,
    n0: float,
    la: np.ndarray,
    L: int,
    M_B: int,
    start_idx: int,
    start_state: int,
    start_metric: float,
    horizon: int,
) -> float:
    """Run a small forward M-search from a given state and return LLR at start_idx.

    Returns the LLR for the symbol at start_idx based on the local search.
    """
    N = y.size
    mask = (1 << L) - 1
    actual_horizon = min(horizon, N - start_idx)

    if actual_horizon <= 0:
        return 0.0

    max_cand = 2 * M_B + 2

    # Initialise
    surv_s = np.empty(M_B, dtype=np.int64)
    surv_m = np.empty(M_B, dtype=np.float64)
    surv_first_sym = np.empty(M_B, dtype=np.int64)  # 0 for -1, 1 for +1

    # First step: expand from start_state
    n_cand = 0
    cand_s = np.empty(max_cand, dtype=np.int64)
    cand_m = np.empty(max_cand, dtype=np.float64)
    cand_fs = np.empty(max_cand, dtype=np.int64)

    for sym_idx in range(2):
        sym = -1.0 if sym_idx == 0 else 1.0
        sn = ((start_state << 1) | sym_idx) & mask
        label = _nb_compute_label(sym, start_state, v, L)
        lp = _nb_log_prior(sym, la[start_idx])
        bm = _nb_branch_metric(y[start_idx], label, n0, lp)
        metric = start_metric + bm

        cand_s[n_cand] = sn
        cand_m[n_cand] = metric
        cand_fs[n_cand] = np.int64(sym_idx)
        n_cand += 1

    # Keep top M_B
    keep = min(M_B, n_cand)
    if keep < n_cand:
        order = np.argsort(cand_m[:n_cand])[::-1]
    else:
        order = np.arange(n_cand)

    for ki in range(keep):
        surv_s[ki] = cand_s[order[ki]]
        surv_m[ki] = cand_m[order[ki]]
        surv_first_sym[ki] = cand_fs[order[ki]]
    n_surv = keep

    # Subsequent steps
    for step in range(1, actual_horizon):
        t = start_idx + step
        n_cand = 0

        for si in range(n_surv):
            sp = surv_s[si]
            am = surv_m[si]
            fs = surv_first_sym[si]
            for sym_idx in range(2):
                sym = -1.0 if sym_idx == 0 else 1.0
                sn = ((sp << 1) | sym_idx) & mask
                label = _nb_compute_label(sym, sp, v, L)
                lp = _nb_log_prior(sym, la[t])
                bm = _nb_branch_metric(y[t], label, n0, lp)
                metric = am + bm

                found = False
                for ci in range(n_cand):
                    if cand_s[ci] == sn and cand_fs[ci] == fs:
                        cand_m[ci] = _nb_logaddexp(cand_m[ci], metric)
                        found = True
                        break
                if not found:
                    cand_s[n_cand] = sn
                    cand_m[n_cand] = metric
                    cand_fs[n_cand] = fs
                    n_cand += 1

        keep = min(M_B, n_cand)
        if keep < n_cand:
            order = np.argsort(cand_m[:n_cand])[::-1]
        else:
            order = np.arange(n_cand)

        for ki in range(keep):
            surv_s[ki] = cand_s[order[ki]]
            surv_m[ki] = cand_m[order[ki]]
            surv_first_sym[ki] = cand_fs[order[ki]]
        n_surv = keep

    # Compute LLR from first symbol
    plus_v = np.empty(n_surv, dtype=np.float64)
    minus_v = np.empty(n_surv, dtype=np.float64)
    n_plus = 0
    n_minus = 0
    for si in range(n_surv):
        if surv_first_sym[si] == 1:
            plus_v[n_plus] = surv_m[si]
            n_plus += 1
        else:
            minus_v[n_minus] = surv_m[si]
            n_minus += 1

    LAMBDA = LLR_RESERVE
    if n_plus > 0 and n_minus > 0:
        raw = _nb_logsumexp(plus_v, n_plus) - _nb_logsumexp(minus_v, n_minus)
        return 0.0 if np.isnan(raw) else raw
    elif n_plus > 0:
        return LAMBDA
    elif n_minus > 0:
        return -LAMBDA
    return 0.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def prlja_mbcjr_bpsk(
    y: np.ndarray,
    v: np.ndarray,
    n0: float,
    M: int,
    la: np.ndarray | None = None,
) -> PrljaMbcjrResult:
    """Paper [14] Simple Detection M-BCJR for BPSK Forney model.

    Parameters
    ----------
    y : np.ndarray
        Forney model observations (white noise).
    v : np.ndarray
        Minimum-phase ISI sequence from spectral_factorize.
    n0 : float
        One-sided noise power spectral density.
    M : int
        Number of survivors in M-algorithm.
    la : np.ndarray, optional
        A-priori LLRs. Zeros if None.

    Returns
    -------
    PrljaMbcjrResult
        LLR sequence and per-symbol survivor counts.
    """
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    v = np.asarray(v, dtype=np.float64).reshape(-1)
    N = y.size
    L_mem = len(v) - 1  # memory depth from v length

    # Effective state length: need at least 1 bit to encode current symbol
    L = max(L_mem, 1)

    # Pad v with zeros if needed (when L > L_mem)
    if L > L_mem:
        v_padded = np.zeros(L + 1, dtype=np.float64)
        v_padded[:len(v)] = v
        v = v_padded

    if M <= 0:
        raise ValueError("M must be positive.")
    if N == 0:
        return PrljaMbcjrResult(llr=np.array([], dtype=np.float64),
                                survivor_counts=[])

    if la is None:
        la_arr = np.zeros(N, dtype=np.float64)
    else:
        la_arr = np.asarray(la, dtype=np.float64).reshape(-1)
        if la_arr.size != N:
            raise ValueError(f"la must have length {N}, got {la_arr.size}.")

    llr, surv_counts = _nb_prlja_simple_kernel(y, v, n0, la_arr, L, M)
    llr = np.nan_to_num(llr, nan=0.0, posinf=LLR_RESERVE, neginf=-LLR_RESERVE)
    return PrljaMbcjrResult(llr=llr, survivor_counts=surv_counts.tolist())


def prlja_backup_mbcjr_bpsk(
    y: np.ndarray,
    v: np.ndarray,
    n0: float,
    M: int,
    M_B: int = 2,
    smooth: bool = False,
    la: np.ndarray | None = None,
) -> PrljaMbcjrResult:
    """Paper [14] Backup M-BCJR with optional smoothing.

    1. Run Simple M-BCJR, get hard decisions and identify positions where
       L+ or L- is empty (reserve LLR used).
    2. For each such position, run a small forward M-search (M_B survivors)
       starting from the hard-decision state at that position.
    3. Replace reserve values with backup LLRs.
    4. If smooth=True: apply smoothing filter [1,3,1]/5 to backup positions.

    Parameters
    ----------
    y : np.ndarray
        Forney model observations.
    v : np.ndarray
        Minimum-phase ISI sequence.
    n0 : float
        Noise power spectral density.
    M : int
        Main M-BCJR survivor count.
    M_B : int
        Backup search survivor count.
    smooth : bool
        Whether to apply smoothing filter to backup positions.
    la : np.ndarray, optional
        A-priori LLRs.

    Returns
    -------
    PrljaMbcjrResult
        LLR sequence with backup corrections applied.
    """
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    v = np.asarray(v, dtype=np.float64).reshape(-1)
    N = y.size

    L_mem = len(v) - 1
    L = max(L_mem, 1)
    if L > L_mem:
        v_padded = np.zeros(L + 1, dtype=np.float64)
        v_padded[:len(v)] = v
        v = v_padded

    # Prepare la
    if la is None:
        la_arr = np.zeros(N, dtype=np.float64)
    else:
        la_arr = np.asarray(la, dtype=np.float64).reshape(-1)

    # Step 1: Run Simple M-BCJR
    result = prlja_mbcjr_bpsk(y, v, n0, M, la)
    llr = result.llr.copy()

    LAMBDA = LLR_RESERVE
    backup_positions = []

    # Identify reserve positions
    for idx in range(N):
        if abs(abs(llr[idx]) - LAMBDA) < 0.01:
            backup_positions.append(idx)

    if not backup_positions:
        return PrljaMbcjrResult(llr=llr, survivor_counts=result.survivor_counts)

    # Step 2: For each reserve position, run backup search
    hard_decisions = (llr >= 0).astype(np.int64)  # 1 for +1, 0 for -1

    for idx in backup_positions:
        # Reconstruct state from hard decisions: state bits 0..L-1
        # bit 0 = newest = x_{idx-1}, bit 1 = x_{idx-2}, etc.
        state = 0
        for k in range(L):
            pos = idx - 1 - k
            if pos >= 0:
                bit = int(hard_decisions[pos])
            else:
                bit = 1  # assume +1 for positions before start
            state |= (bit << k)

        backup_llr = _nb_backup_forward(
            y, v, n0, la_arr, L, M_B, idx, state, 0.0, 6,
        )
        llr[idx] = backup_llr

    # Step 3: Optional smoothing
    if smooth and len(backup_positions) > 0:
        kernel = np.array([1.0, 3.0, 1.0]) / 5.0
        smoothed = llr.copy()
        for idx in backup_positions:
            if 0 < idx < N - 1:
                v0 = llr[idx - 1]
                v1 = llr[idx]
                v2 = llr[idx + 1]
                if np.isfinite(v0) and np.isfinite(v1) and np.isfinite(v2):
                    smoothed[idx] = kernel[0] * v0 + kernel[1] * v1 + kernel[2] * v2
        for idx in backup_positions:
            if 0 < idx < N - 1:
                if np.isfinite(smoothed[idx]):
                    llr[idx] = smoothed[idx]

    llr = np.nan_to_num(llr, nan=0.0, posinf=LLR_RESERVE, neginf=-LLR_RESERVE)
    return PrljaMbcjrResult(llr=llr, survivor_counts=result.survivor_counts)
