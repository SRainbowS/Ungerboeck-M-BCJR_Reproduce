"""Reduced-memory BCJR equaliser operating on the shortened channel metric.

Paper [26] (Rusek & Prlja) channel shortening produces a banded positive-
definite matrix G (bandwidth *nu*) and a shortening filter h_r.  After
filtering the received signal  z = h_r * y, the reduced-memory BCJR runs on
a trellis with only 2^nu states (instead of 2^L for the full channel).

Branch metric (BPSK, x in {-1, +1}):

    gamma_k(x_k, state) = 2 x_k z_k
                        - G[k,k]
                        - 2 x_k sum_{i=1}^{nu} G[k,k-i] x_{k-i}
                        + log_prior(x_k)

Since G is Toeplitz (LTI channel), G[k,k] is constant and G[k,k-i] depends
only on i.  The G[k,k] term is constant across branches at time k and
cancels in LLR computation, but is included for correctness.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from ftn.baselines.channel_shortening import (
    apply_shortening_filter_fft,
    compute_shortened_params_from_v,
)


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


# ---------------------------------------------------------------------------
# JIT kernel — full forward-backward BCJR on the reduced trellis
# ---------------------------------------------------------------------------

@njit(cache=True)
def _nb_shortened_bcjr_kernel(
    z: np.ndarray,
    g_diag: float,
    g_off: np.ndarray,
    nu: int,
    la: np.ndarray,
) -> np.ndarray:
    """Full forward-backward BCJR on the shortened channel model.

    State encoding: integer *s*, bit i = symbol at position i
    (0 -> -1, 1 -> +1).  Bit 0 is the most recent past symbol x_{k-1}.
    State transition: s' = ((s << 1) | b) & mask  where mask = (1 << nu) - 1.
    """
    n = z.size
    n_states = 1 << nu
    mask = n_states - 1

    # --- Forward pass (alpha) ---
    # alpha[k, s] = log P(z_0..z_{k-1}, state at time k = s)
    alpha = np.full((n + 1, n_states), -np.inf)
    # Start with uniform initial state (all +1)
    init_state = n_states - 1  # all bits set -> all +1
    alpha[0, init_state] = 0.0

    for k in range(n):
        z_k = z[k]
        la_k = la[k]
        for s_prev in range(n_states):
            a_prev = alpha[k, s_prev]
            if a_prev == -np.inf:
                continue
            # Extract past symbols from state for ISI computation
            for bit in range(2):
                x_k = 1.0 if bit == 1 else -1.0
                s_new = ((s_prev << 1) | bit) & mask

                # ISI from past symbols: sum_{i=1}^{nu} g_off[i-1] * x_{k-i}
                isi = 0.0
                for i in range(nu):
                    x_past = 1.0 if (s_prev >> i) & 1 else -1.0
                    isi += g_off[i] * x_past

                # Branch metric
                gamma = (2.0 * x_k * z_k
                         - g_diag
                         - 2.0 * x_k * isi)

                # Prior
                if bit == 1:
                    prior = -_nb_logaddexp(0.0, -la_k)
                else:
                    prior = -_nb_logaddexp(0.0, la_k)
                gamma += prior

                alpha[k + 1, s_new] = _nb_logaddexp(
                    alpha[k + 1, s_new], a_prev + gamma
                )

        # Normalise alpha to prevent overflow
        norm = _nb_logsumexp(alpha[k + 1, :], n_states)
        for s in range(n_states):
            alpha[k + 1, s] -= norm

    # --- Backward pass (beta) ---
    # beta[k, s] = log P(z_k..z_{n-1} | state at time k = s)
    beta = np.full((n + 1, n_states), -np.inf)
    for s in range(n_states):
        beta[n, s] = 0.0

    for k in range(n - 1, -1, -1):
        z_k = z[k]
        la_k = la[k]
        for s_prev in range(n_states):
            values = np.empty(2)
            for bit in range(2):
                x_k = 1.0 if bit == 1 else -1.0
                s_new = ((s_prev << 1) | bit) & mask

                isi = 0.0
                for i in range(nu):
                    x_past = 1.0 if (s_prev >> i) & 1 else -1.0
                    isi += g_off[i] * x_past

                gamma = (2.0 * x_k * z_k
                         - g_diag
                         - 2.0 * x_k * isi)

                if bit == 1:
                    prior = -_nb_logaddexp(0.0, -la_k)
                else:
                    prior = -_nb_logaddexp(0.0, la_k)
                gamma += prior

                values[bit] = gamma + beta[k + 1, s_new]

            beta[k, s_prev] = _nb_logsumexp(values, 2)

        # Normalise beta
        norm = _nb_logsumexp(beta[k, :], n_states)
        for s in range(n_states):
            beta[k, s] -= norm

    # --- LLR computation ---
    llr = np.zeros(n)
    for k in range(n):
        z_k = z[k]
        la_k = la[k]
        plus_vals = np.empty(n_states)
        minus_vals = np.empty(n_states)
        n_plus = 0
        n_minus = 0
        for s_prev in range(n_states):
            a_val = alpha[k, s_prev]
            if a_val == -np.inf:
                continue
            for bit in range(2):
                x_k = 1.0 if bit == 1 else -1.0
                s_new = ((s_prev << 1) | bit) & mask

                isi = 0.0
                for i in range(nu):
                    x_past = 1.0 if (s_prev >> i) & 1 else -1.0
                    isi += g_off[i] * x_past

                gamma = (2.0 * x_k * z_k
                         - g_diag
                         - 2.0 * x_k * isi)

                if bit == 1:
                    prior = -_nb_logaddexp(0.0, -la_k)
                else:
                    prior = -_nb_logaddexp(0.0, la_k)
                gamma += prior

                metric = a_val + gamma + beta[k + 1, s_new]
                if bit == 1:
                    plus_vals[n_plus] = metric
                    n_plus += 1
                else:
                    minus_vals[n_minus] = metric
                    n_minus += 1

        llr[k] = _nb_logsumexp(plus_vals, n_plus) - _nb_logsumexp(
            minus_vals, n_minus
        )

    return llr


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ShortenedBcjrResult:
    """Result from the shortened BCJR equaliser."""
    llr: np.ndarray


def shortened_bcjr_bpsk(
    z: np.ndarray,
    g_diag: float,
    g_off: np.ndarray,
    nu: int,
    la: np.ndarray | None = None,
) -> ShortenedBcjrResult:
    """Paper [26] reduced-memory BCJR on the shortened metric.

    Parameters
    ----------
    z : ndarray, shape (N,)
        Shortened observation (after applying the shortening filter to y).
    g_diag : float
        Diagonal value  G[k,k]  (constant for Toeplitz G).
    g_off : ndarray, shape (nu,)
        Off-diagonal values  G[k, k-i]  for  i = 1 .. nu.
    nu : int
        Memory (bandwidth) of the shortened channel.
    la : ndarray or None
        A-priori LLRs, shape (N,).  Zero if not provided.

    Returns
    -------
    ShortenedBcjrResult with field ``llr``.
    """
    z = np.asarray(z, dtype=np.float64).reshape(-1)
    n = z.size
    g_off = np.asarray(g_off, dtype=np.float64).reshape(-1)
    if g_off.size != nu:
        raise ValueError(f"g_off must have length {nu}, got {g_off.size}.")
    if la is None:
        la_arr = np.zeros(n, dtype=np.float64)
    else:
        la_arr = np.asarray(la, dtype=np.float64).reshape(-1)
        if la_arr.size != n:
            raise ValueError(f"la must have length {n}, got {la_arr.size}.")

    llr = _nb_shortened_bcjr_kernel(z, float(g_diag), g_off, nu, la_arr)
    return ShortenedBcjrResult(llr=llr)


def cs_equalizer_bpsk(
    y: np.ndarray,
    v: np.ndarray,
    n0: float,
    nu: int,
    la: np.ndarray | None = None,
    n_fft: int = 0,
) -> ShortenedBcjrResult:
    """Complete channel-shortening equaliser pipeline.

    Uses the Forney (white noise) model: y = V*x + w, w ~ N(0, n0*I).
    The n0 parameter is the noise variance per dimension (= N0/2 in standard
    notation for single-sided PSD N0).

    1. Precompute shortening parameters from v, n0, nu.
    2. Apply frequency-domain z-filter to y -> z.
    3. Run reduced-memory BCJR.
    """
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    n = len(y)

    # Ensure n_fft is large enough for the signal
    if n_fft < n:
        n_fft = 1
        while n_fft < n + len(v):
            n_fft *= 2

    Z_w, g_diag, g_off, n_fft_used = compute_shortened_params_from_v(
        v, n0, nu, n_fft=n_fft,
    )
    z = apply_shortening_filter_fft(y, Z_w, n_fft_used)

    return shortened_bcjr_bpsk(z, g_diag, g_off, nu, la)
