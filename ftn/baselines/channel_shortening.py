"""Channel shortening optimisation from Paper [26]: Rusek & Prlja.

For ISI channels, the observation model is  y = H x + w  where H is the
convolution matrix of the minimum-phase channel response v and w ~ N(0, N0 I).

The paper derives (Proposition 2) a closed-form banded positive-definite matrix
G_code that maximises the achievable information rate.  For ISI channels, the
B matrix is characterised by the spectral density

    B(omega) = N0 / (|V(omega)|^2 + N0)

and the optimal G_code is obtained via a modified Cholesky decomposition of B.

The z-filter in the frequency domain (Proposition 3) is

    Z(omega) = G_code(omega) * conj(V(omega)) / (|V(omega)|^2 + N0)

and the observation z is computed via circular convolution

    z = IFFT( Z(omega) * FFT(y) ).

Branch metric (BPSK, x in {-1, +1}):

    gamma_k(x_k, state) = 2 x_k z_k
                        - G_code[k,k]
                        - 2 x_k sum_{i=1}^{nu} G_code[k,k-i] x_{k-i}
                        + log_prior(x_k)
"""

from __future__ import annotations

import numpy as np
import scipy.linalg


# ---------------------------------------------------------------------------
# Frequency-domain channel shortening (correct approach from Paper [26])
# ---------------------------------------------------------------------------

def compute_shortened_params_from_v(
    v: np.ndarray,
    n0: float,
    nu: int,
    n_fft: int = 0,
) -> tuple[np.ndarray, float, np.ndarray, int]:
    """Compute channel shortening parameters from min-phase sequence v.

    Uses Proposition 2 (closed-form modified Cholesky) and Proposition 3
    (frequency-domain z-filter) from Paper [26].

    Parameters
    ----------
    v : ndarray
        Minimum-phase channel impulse response (causal, real-valued).
    n0 : float
        Noise variance per dimension (for Forney model: N0/2).
    nu : int
        Desired trellis memory (bandwidth of G_code).
    n_fft : int
        FFT size.  If 0 (default), auto-selected as next power of 2 >= 8*len(v).

    Returns
    -------
    Z_w : ndarray
        Frequency-domain z-filter (rfft coefficients).
    g_diag : float
        Diagonal value G_code[k,k] (constant for Toeplitz interior).
    g_off : ndarray, shape (nu,)
        Off-diagonal values G_code[k,k-i] for i = 1 .. nu.
    n_fft : int
        FFT size used (for apply_shortening_filter_fft).
    """
    v = np.asarray(v, dtype=np.float64).reshape(-1)
    if n_fft <= 0:
        n_fft = 1
        while n_fft < 8 * len(v):
            n_fft *= 2

    # Compute spectra
    v_pad = np.zeros(n_fft)
    v_pad[:len(v)] = v
    V_w = np.fft.rfft(v_pad)
    G_w = np.abs(V_w) ** 2

    # B(omega) = N0 / (G(omega) + N0)
    B_w = n0 / (G_w + n0)

    # b[k] = IDFT(B(omega))[k] for k = 0, ..., nu
    b_full = np.fft.irfft(B_w, n=n_fft)
    b = b_full[:nu + 1]

    # Proposition 2: modified Cholesky decomposition
    # b_vec = [b[1], ..., b[nu]]
    # B_tilde is nu x nu Toeplitz from b[0], ..., b[nu-1]
    b_vec = b[1:nu + 1]

    if nu > 0:
        B_tilde = np.zeros((nu, nu))
        for i in range(nu):
            for j in range(nu):
                B_tilde[i, j] = b[abs(i - j)]

        B_tilde_inv = np.linalg.inv(B_tilde)
        c = b[0] - b_vec @ B_tilde_inv @ b_vec

        u_nn = 1.0 / np.sqrt(c)
        u_off = -u_nn * b_vec @ B_tilde_inv

        # G_code(omega) = |U(omega)|^2  where U row = [u_nn, u_off[0], ...]
        u_full = np.zeros(n_fft)
        u_full[0] = u_nn
        u_full[1:nu + 1] = u_off
        U_w = np.fft.rfft(u_full)
        G_code_w = np.abs(U_w) ** 2
    else:
        c = b[0]
        G_code_w = np.ones(len(V_w)) / c

    # z-filter(omega) = G_code(omega) * conj(V(omega)) / (G(omega) + N0)
    denom = G_w + n0
    denom[denom == 0] = 1e-30
    Z_w = G_code_w * np.conj(V_w) / denom

    # G_code time-domain: extract Toeplitz parameters
    g_code_t = np.fft.irfft(G_code_w, n=n_fft)
    g_diag = g_code_t[0]
    g_off = g_code_t[1:nu + 1]

    return Z_w, g_diag, g_off, n_fft


def apply_shortening_filter_fft(
    y: np.ndarray,
    Z_w: np.ndarray,
    n_fft: int,
) -> np.ndarray:
    """Apply the frequency-domain z-filter to the received signal y.

    Uses circular convolution via FFT, which for zero-padded signals is
    equivalent to linear convolution.
    """
    n = len(y)
    y_pad = np.zeros(n_fft)
    y_pad[:n] = y
    Y_w = np.fft.rfft(y_pad)
    z = np.fft.irfft(Z_w * Y_w, n=n_fft)
    return z[:n]


# ---------------------------------------------------------------------------
# Legacy matrix-based functions (kept for backward compatibility with tests)
# ---------------------------------------------------------------------------

def compute_b_matrix(
    g: dict[int, float],
    n0: float,
    block_len: int = 200,
) -> np.ndarray:
    """Compute B = I - H^T (H H^T + N0 I)^{-1} H for a short block."""
    from ftn.channel import build_isi_matrix
    H = build_isi_matrix(g, block_len)
    A = H @ H.T + n0 * np.eye(block_len)
    L = np.linalg.cholesky(A)
    A_inv_H = scipy.linalg.cho_solve((L, True), H)
    B = np.eye(block_len) - H.T @ A_inv_H
    B = 0.5 * (B + B.T)
    return B


def solve_banded_g(
    B: np.ndarray,
    nu: int,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> np.ndarray:
    """Find banded positive-definite G maximising log det(G) - trace(G B)."""
    N = B.shape[0]
    G = np.eye(N, dtype=np.float64)

    def _project_banded(M: np.ndarray) -> np.ndarray:
        out = M.copy()
        for i in range(N):
            for j in range(N):
                if abs(i - j) > nu:
                    out[i, j] = 0.0
        out = 0.5 * (out + out.T)
        return out

    def _ensure_pd(M: np.ndarray, eps: float = 1e-6) -> np.ndarray:
        eigvals = np.linalg.eigvalsh(M)
        min_eig = eigvals.min()
        if min_eig < eps:
            M = M + (eps - min_eig) * np.eye(N)
        return M

    def _objective(M: np.ndarray) -> float:
        sign, logdet = np.linalg.slogdet(M)
        if sign <= 0:
            return -np.inf
        return logdet - np.trace(M @ B)

    obj_prev = _objective(G)

    for iteration in range(max_iter):
        try:
            G_inv = np.linalg.inv(G)
        except np.linalg.LinAlgError:
            G_inv = np.linalg.pinv(G)

        grad = G_inv - B

        alpha = 1.0
        for _ in range(30):
            G_new = G + alpha * grad
            G_new = _project_banded(G_new)
            G_new = _ensure_pd(G_new)
            obj_new = _objective(G_new)
            if obj_new > obj_prev - 1e-12:
                break
            alpha *= 0.5
        else:
            G_new = G + alpha * grad
            G_new = _project_banded(G_new)
            G_new = _ensure_pd(G_new)
            obj_new = _objective(G_new)

        if abs(obj_new - obj_prev) < tol:
            G = G_new
            break

        G = G_new
        obj_prev = obj_new

    return G


def compute_shortened_params(
    g: dict[int, float],
    n0: float,
    nu: int,
    block_len: int = 200,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Legacy: precompute channel shortening parameters from g-taps."""
    from ftn.channel import build_isi_matrix
    H = build_isi_matrix(g, block_len)
    B = compute_b_matrix(g, n0, block_len)
    G = solve_banded_g(B, nu)

    A = H @ H.T + n0 * np.eye(block_len)
    L = np.linalg.cholesky(A)
    A_inv_H = scipy.linalg.cho_solve((L, True), H)
    H_r = G @ H.T @ A_inv_H

    centre = block_len // 2
    h_r = H_r[centre, :]
    g_diag = G[centre, centre]
    g_off = np.array([G[centre, centre - i] for i in range(1, nu + 1)])

    return h_r, g_diag, g_off


def apply_shortening_filter(y: np.ndarray, h_r: np.ndarray) -> np.ndarray:
    """Legacy: apply FIR shortening filter h_r to y."""
    n = len(y)
    z = np.convolve(y, h_r, mode='full')
    centre = len(h_r) // 2
    return z[centre:centre + n]
