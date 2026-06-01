from __future__ import annotations

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz


def rrc_pulse(t: np.ndarray, beta: float = 0.3, T: float = 1.0) -> np.ndarray:
    """Evaluate a unit-symbol-period root-raised-cosine pulse."""
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must be in [0, 1].")
    if T <= 0.0:
        raise ValueError("T must be positive.")

    t = np.asarray(t, dtype=float)
    u = t / T
    h = np.empty_like(u, dtype=float)

    if beta == 0.0:
        h[:] = np.sinc(u) / np.sqrt(T)
        return h

    zero = np.isclose(u, 0.0, atol=1e-12)
    singular = np.isclose(np.abs(u), 1.0 / (4.0 * beta), atol=1e-12)
    normal = ~(zero | singular)

    h[zero] = (1.0 + beta * (4.0 / np.pi - 1.0)) / np.sqrt(T)
    h[singular] = (
        beta
        / np.sqrt(2.0 * T)
        * (
            (1.0 + 2.0 / np.pi) * np.sin(np.pi / (4.0 * beta))
            + (1.0 - 2.0 / np.pi) * np.cos(np.pi / (4.0 * beta))
        )
    )

    un = u[normal]
    numerator = (
        np.sin(np.pi * un * (1.0 - beta))
        + 4.0 * beta * un * np.cos(np.pi * un * (1.0 + beta))
    )
    denominator = np.pi * un * (1.0 - (4.0 * beta * un) ** 2)
    h[normal] = numerator / denominator / np.sqrt(T)
    return h


def generate_rrc(
    beta: float = 0.3,
    span: int = 15,
    sps: int = 128,
    T: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a continuous-time sampled rRC pulse over ``[-span*T, span*T]``."""
    if span <= 0:
        raise ValueError("span must be positive.")
    if sps <= 0:
        raise ValueError("sps must be positive.")

    step = T / float(sps)
    t = np.arange(-span * T, span * T + 0.5 * step, step, dtype=float)
    h = rrc_pulse(t, beta=beta, T=T)
    energy = _trapz(np.abs(h) ** 2, t)
    if energy <= 0.0:
        raise ValueError("pulse energy must be positive.")
    h = h / np.sqrt(energy)
    return t, h


def compute_g(
    t: np.ndarray,
    h: np.ndarray,
    tau: float,
    isi_len: int,
    T: float = 1.0,
) -> dict[int, float]:
    """Compute sampled matched-filter autocorrelation taps ``g_l``."""
    if tau <= 0.0:
        raise ValueError("tau must be positive.")
    if isi_len < 0:
        raise ValueError("isi_len must be non-negative.")

    t = np.asarray(t, dtype=float)
    h = np.asarray(h, dtype=float)
    if t.ndim != 1 or h.ndim != 1 or t.shape != h.shape:
        raise ValueError("t and h must be one-dimensional arrays with the same shape.")

    taps: dict[int, float] = {}
    for lag in range(-isi_len, isi_len + 1):
        shifted = np.interp(t - lag * tau * T, t, h, left=0.0, right=0.0)
        taps[lag] = float(_trapz(h * shifted, t))
    return taps
