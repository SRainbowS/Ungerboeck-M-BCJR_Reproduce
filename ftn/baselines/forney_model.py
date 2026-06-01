"""Forney (whitened matched filter) channel model via minimum-phase conversion.

Converts the sampled rRC pulse to its minimum-phase version, which serves as
the causal ISI sequence v in the Forney observation model:

    y[n] = sum_{k=0}^{L} v[k] * x[n-k] + w[n],  w ~ N(0, N0/2)

The conversion follows Paper [14] (Prlja & Anderson): sample the rRC pulse at
the FTN signalling rate, then compute the minimum-phase version via the Hilbert
cepstral method.  The result is scaled so that sum(v^2) = g[0], ensuring the
Forney model has the same SNR as the Ungerboeck model.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import lfilter, minimum_phase


def sample_pulse_at_ftn_rate(
    t_pulse: np.ndarray,
    h_pulse: np.ndarray,
    tau: float,
    T: float = 1.0,
) -> np.ndarray:
    """Sample the rRC pulse at the FTN signalling rate 1/(tau*T).

    Returns the sampled sequence c[k] = h(k * tau * T) for k = -K..K,
    where K = span/tau.  The output is a 1-D array of length 2K+1.
    """
    dt = t_pulse[1] - t_pulse[0]
    ftn_step = tau * T
    ftn_step_idx = int(round(ftn_step / dt))
    center_idx = len(h_pulse) // 2
    span = int(round((t_pulse[-1] - t_pulse[0]) / (2 * T)))
    max_k = int(span / tau)

    c = np.zeros(2 * max_k + 1)
    for k in range(-max_k, max_k + 1):
        idx = center_idx + k * ftn_step_idx
        if 0 <= idx < len(h_pulse):
            c[k + max_k] = h_pulse[idx]
    return c


def min_phase_from_pulse(
    t_pulse: np.ndarray,
    h_pulse: np.ndarray,
    tau: float,
    g0: float = 1.0,
    T: float = 1.0,
) -> np.ndarray:
    """Compute the minimum-phase ISI sequence v from the rRC pulse.

    Steps (following Paper [14]):
    1. Sample the pulse at the FTN rate: c[k] = h(k*tau*T)
    2. Compute the minimum-phase version via cepstral method (scipy)
    3. Scale so that sum(v^2) = g0 (the main autocorrelation tap)

    The result v is the causal minimum-phase ISI sequence for the Forney model.
    """
    c = sample_pulse_at_ftn_rate(t_pulse, h_pulse, tau, T)
    v = minimum_phase(c, method="hilbert")
    if np.any(np.isnan(v)):
        v = minimum_phase(c, method="homomorphic")
    # Scale so that sum(v^2) = g0
    energy = np.sum(v ** 2)
    if energy > 0:
        v = v * np.sqrt(g0 / energy)
    return v


def spectral_factorize(g: dict[int, float]) -> np.ndarray:
    """Spectral factorization from autocorrelation g-taps (legacy, less accurate).

    For best results, use ``min_phase_from_pulse`` instead, which operates on
    the pulse directly and avoids numerical issues with high-degree polynomials.
    """
    if not g:
        return np.array([1.0])

    lags = [abs(int(k)) for k in g.keys()]
    L = max(lags)
    r = np.array([float(g.get(lag, 0.0)) for lag in range(L + 1)], dtype=np.float64)

    if L == 0:
        return np.array([np.sqrt(r[0])])

    # Build polynomial P(z) = z^L * G(z)
    poly_coeffs = np.zeros(2 * L + 1, dtype=np.float64)
    poly_coeffs[0] = r[L]
    for l in range(1, L):
        poly_coeffs[l] = r[L - l]
    poly_coeffs[L] = r[0]
    for l in range(1, L + 1):
        poly_coeffs[L + l] = r[l]

    roots = np.roots(poly_coeffs)

    # Classify roots: inside, on-unit-circle, outside
    eps_unit = 1e-6
    inside = []
    on_circle = []
    for root in roots:
        if np.abs(root) < 1.0 - eps_unit:
            inside.append(root)
        elif np.abs(root) > 1.0 + eps_unit:
            pass
        else:
            on_circle.append(root)

    # Select one from each conjugate pair on the unit circle
    selected_circle = []
    used = np.zeros(len(on_circle), dtype=bool)
    for i, root in enumerate(on_circle):
        if used[i]:
            continue
        found_partner = False
        for j in range(i + 1, len(on_circle)):
            if used[j]:
                continue
            if np.abs(root - np.conj(on_circle[j])) < 1e-6:
                used[j] = True
                found_partner = True
                break
        if found_partner:
            selected_circle.append(root if root.imag >= 0 else np.conj(root))
        else:
            selected_circle.append(root)

    min_phase_roots = list(inside) + selected_circle
    if len(min_phase_roots) != L:
        abs_roots = np.abs(roots)
        order = np.argsort(abs_roots)
        min_phase_roots = list(roots[order[:L]])

    v_poly = np.real(np.poly(min_phase_roots))
    g_dc = float(np.sum(r) + np.sum(r[1:]))
    v_dc = float(np.sum(v_poly))
    if abs(v_dc) > 1e-15:
        scale = np.sqrt(abs(g_dc)) / abs(v_dc)
    else:
        scale = np.sqrt(r[0] / np.sum(v_poly ** 2))
    v = v_poly * scale
    if v[0] < 0:
        v = -v
    return v


def forney_channel(
    x: np.ndarray,
    v: np.ndarray,
    n0: float,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Generate Forney model observation: y = conv(x, v) + white_noise.

    The convolution is causal: y[n] = sum_{k=0}^{L} v[k] * x[n-k] + w[n],
    where w ~ N(0, N0/2).
    """
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    v = np.asarray(v, dtype=np.float64).reshape(-1)
    signal = lfilter(v, [1.0], x)
    if n0 > 0.0:
        if rng is None:
            rng = np.random.default_rng()
        noise = rng.standard_normal(x.size) * np.sqrt(n0 / 2.0)
    else:
        noise = 0.0
    return signal + noise


def forney_branch_metric(
    y_n: float, label: float, n0: float, log_prior: float = 0.0
) -> float:
    """Forney branch metric: -(y_n - label)^2 / N0 + log_prior."""
    diff = y_n - label
    return -(diff * diff) / n0 + log_prior
