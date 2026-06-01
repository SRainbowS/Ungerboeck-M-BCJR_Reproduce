from __future__ import annotations

import numpy as np


def log_prior_bpsk_symbol(symbol: float, la: float = 0.0) -> float:
    """Return log P(symbol) for BPSK with LLR ``la = log P(+1)/P(-1)``."""
    if symbol > 0:
        return float(-np.logaddexp(0.0, -float(la)))
    if symbol < 0:
        return float(-np.logaddexp(0.0, float(la)))
    raise ValueError("BPSK symbol must be +1 or -1.")


def log_phi_ungerboeck(
    y_n: float,
    x_n: float,
    state_prev: tuple[float, ...],
    g: dict[int, float],
    n0: float,
) -> float:
    """Ungerboeck matched-filter branch metric for real BPSK."""
    if n0 <= 0.0:
        raise ValueError("n0 must be positive.")
    x_n = float(x_n)
    isi = 0.0
    for lag in range(1, len(state_prev) + 1):
        isi += float(g.get(lag, 0.0)) * float(state_prev[-lag])
    centered = float(y_n) - 0.5 * float(g.get(0, 0.0)) * x_n - isi
    return float((2.0 / n0) * x_n * centered)


def log_gamma_bpsk(
    y_n: float,
    x_n: float,
    state_prev: tuple[float, ...],
    g: dict[int, float],
    n0: float,
    la_n: float = 0.0,
) -> float:
    """BPSK branch metric including the symbol prior."""
    return log_phi_ungerboeck(y_n, x_n, state_prev, g, n0) + log_prior_bpsk_symbol(
        x_n, la_n
    )
