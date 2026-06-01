from __future__ import annotations

import numpy as np
from scipy.signal import lfilter


def build_isi_matrix(g: dict[int, float], n: int) -> np.ndarray:
    """Build the finite Toeplitz matrix for ``y_n = sum_l g_l x_{n-l}``."""
    if n <= 0:
        raise ValueError("n must be positive.")
    G = np.zeros((n, n), dtype=float)
    for row in range(n):
        for col in range(n):
            G[row, col] = float(g.get(row - col, 0.0))
    return G


def ftn_filter_output(x: np.ndarray, g: dict[int, float]) -> np.ndarray:
    """Apply the finite matched-filter ISI response without noise."""
    symbols = np.asarray(x, dtype=float).reshape(-1)
    y = np.zeros_like(symbols, dtype=float)
    for row in range(symbols.size):
        total = 0.0
        for lag, tap in g.items():
            col = row - int(lag)
            if 0 <= col < symbols.size:
                total += float(tap) * symbols[col]
        y[row] = total
    return y


def sample_colored_noise(
    g: dict[int, float],
    n: int,
    n0: float,
    rng: np.random.Generator | None = None,
    size: int | None = None,
) -> np.ndarray:
    """Sample real colored matched-filter noise with covariance ``N0/2 * G``."""
    if n0 < 0.0:
        raise ValueError("n0 must be non-negative.")
    rng = np.random.default_rng() if rng is None else rng
    cov = (n0 / 2.0) * build_isi_matrix(g, n)
    if n0 == 0.0:
        shape = (n,) if size is None else (int(size), n)
        return np.zeros(shape, dtype=float)

    jitter = 1e-12 * np.eye(n)
    try:
        factor = np.linalg.cholesky(cov + jitter)
    except np.linalg.LinAlgError:
        eigvals, eigvecs = np.linalg.eigh(cov)
        eigvals = np.clip(eigvals, 0.0, None)
        factor = eigvecs @ np.diag(np.sqrt(eigvals))

    if size is None:
        white = rng.standard_normal(n)
        return factor @ white
    white = rng.standard_normal((int(size), n))
    return white @ factor.T


def _levinson_durbin(g: dict[int, float]) -> tuple[np.ndarray, float]:
    """Levinson-Durbin spectral factorization of autocorrelation taps ``g``.

    Returns ``(ar_coeffs, sigma2)`` where the AR model
    ``H(z) = sigma / (1 + a_1 z^-1 + ... + a_p z^-p)`` produces noise
    with autocorrelation matching ``g``.
    """
    max_lag = max(abs(int(k)) for k in g.keys())
    r = np.array([g.get(k, 0.0) for k in range(max_lag + 1)])
    a = np.zeros(max_lag + 1)
    a[0] = 1.0
    sigma2 = r[0]
    for p in range(1, max_lag + 1):
        lam = r[p]
        for j in range(1, p):
            lam += a[j] * r[p - j]
        lam /= -sigma2
        a_new = a.copy()
        for j in range(1, p):
            a_new[j] = a[j] + lam * a[p - j]
        a_new[p] = lam
        sigma2 *= 1.0 - lam * lam
        a = a_new
    return a, sigma2


def sample_colored_noise_fast(
    g: dict[int, float],
    n: int,
    n0: float,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Fast colored noise via Levinson-Durbin AR model — O(n*L).

    Uses spectral factorization of the autocorrelation ``g`` to build an AR
    filter, then generates noise with ``scipy.signal.lfilter``.
    """
    if n0 < 0.0:
        raise ValueError("n0 must be non-negative.")
    if n0 == 0.0:
        return np.zeros(n, dtype=float)
    rng = np.random.default_rng() if rng is None else rng
    ar_coeffs, sigma2 = _levinson_durbin(g)
    sigma = np.sqrt(sigma2 * n0 / 2.0)
    w = rng.standard_normal(n)
    return lfilter(np.array([sigma]), ar_coeffs, w)


def ftn_awgn_channel(
    x: np.ndarray,
    g: dict[int, float],
    n0: float,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Generate matched-filter FTN channel output under the Ungerboeck model."""
    symbols = np.asarray(x, dtype=float).reshape(-1)
    n = symbols.size
    if n <= 256:
        noise = sample_colored_noise(g, n, n0, rng)
    else:
        noise = sample_colored_noise_fast(g, n, n0, rng)
    return ftn_filter_output(symbols, g) + noise


# ---------------------------------------------------------------------------
# Complex-valued variants (for QAM)
# ---------------------------------------------------------------------------


def ftn_filter_output_complex(x_complex: np.ndarray, g: dict[int, float]) -> np.ndarray:
    """Apply real-valued ISI taps ``g`` to complex symbols."""
    symbols = np.asarray(x_complex, dtype=complex).reshape(-1)
    y = np.zeros_like(symbols, dtype=complex)
    for row in range(symbols.size):
        total = 0.0 + 0.0j
        for lag, tap in g.items():
            col = row - int(lag)
            if 0 <= col < symbols.size:
                total += float(tap) * symbols[col]
        y[row] = total
    return y


def ftn_awgn_channel_complex(
    x_complex: np.ndarray,
    g: dict[int, float],
    n0: float,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Complex FTN channel: ISI + independent I/Q colored noise."""
    rng = np.random.default_rng() if rng is None else rng
    symbols = np.asarray(x_complex, dtype=complex).reshape(-1)
    n = symbols.size
    signal = ftn_filter_output_complex(symbols, g)
    noise_i = sample_colored_noise_fast(g, n, n0, rng)
    noise_q = sample_colored_noise_fast(g, n, n0, rng)
    return signal + noise_i + 1j * noise_q
