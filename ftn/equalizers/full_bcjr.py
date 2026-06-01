from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.special import logsumexp

from ftn.metrics import log_gamma_bpsk


LOG_ZERO = -np.inf
BPSK_ALPHABET = (-1.0, 1.0)


@dataclass(frozen=True)
class BcjrResult:
    llr: np.ndarray
    log_alpha: list[dict[tuple[float, ...], float]]
    log_beta: list[dict[tuple[float, ...], float]]


def _all_states(isi_len: int) -> list[tuple[float, ...]]:
    return [tuple(float(v) for v in state) for state in product(BPSK_ALPHABET, repeat=isi_len)]


def _shift_state(state: tuple[float, ...], symbol: float) -> tuple[float, ...]:
    if not state:
        return ()
    return tuple((*state[1:], float(symbol)))


def _prepare_la(la: np.ndarray | None, n: int) -> np.ndarray:
    if la is None:
        return np.zeros(n, dtype=float)
    arr = np.asarray(la, dtype=float).reshape(-1)
    if arr.size != n:
        raise ValueError(f"la must have length {n}, got {arr.size}.")
    return arr


def full_bcjr_bpsk(
    y: np.ndarray,
    g: dict[int, float],
    n0: float,
    la: np.ndarray | None = None,
    isi_len: int | None = None,
    initial_state: tuple[float, ...] | None = None,
) -> BcjrResult:
    """Exact full-state BCJR for the finite-memory BPSK Ungerboeck model."""
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.size
    if isi_len is None:
        isi_len = max((abs(int(k)) for k in g.keys()), default=0)
    if initial_state is None:
        initial_state = tuple(1.0 for _ in range(isi_len))
    if len(initial_state) != isi_len:
        raise ValueError("initial_state length must equal isi_len.")
    la_arr = _prepare_la(la, n)
    states = _all_states(isi_len)

    alpha_before: list[dict[tuple[float, ...], float]] = []
    current = {tuple(initial_state): 0.0}
    for idx in range(n):
        alpha_before.append(current)
        next_alpha = {state: LOG_ZERO for state in states}
        for state_prev, log_a in current.items():
            for symbol in BPSK_ALPHABET:
                state_new = _shift_state(state_prev, symbol)
                metric = log_a + log_gamma_bpsk(
                    y[idx], symbol, state_prev, g, n0, la_arr[idx]
                )
                next_alpha[state_new] = np.logaddexp(next_alpha[state_new], metric)
        current = {state: value for state, value in next_alpha.items() if np.isfinite(value)}

    beta_after: list[dict[tuple[float, ...], float]] = [
        {state: 0.0 for state in states} for _ in range(n)
    ]
    if n >= 2:
        for idx in range(n - 2, -1, -1):
            beta_here: dict[tuple[float, ...], float] = {}
            for state in states:
                values = []
                for symbol in BPSK_ALPHABET:
                    state_new = _shift_state(state, symbol)
                    values.append(
                        log_gamma_bpsk(
                            y[idx + 1], symbol, state, g, n0, la_arr[idx + 1]
                        )
                        + beta_after[idx + 1].get(state_new, LOG_ZERO)
                    )
                beta_here[state] = float(logsumexp(values))
            beta_after[idx] = beta_here

    llr = np.zeros(n, dtype=float)
    for idx in range(n):
        plus_values = []
        minus_values = []
        for state_prev, log_a in alpha_before[idx].items():
            for symbol in BPSK_ALPHABET:
                state_new = _shift_state(state_prev, symbol)
                value = (
                    log_a
                    + log_gamma_bpsk(y[idx], symbol, state_prev, g, n0, la_arr[idx])
                    + beta_after[idx].get(state_new, LOG_ZERO)
                )
                if symbol > 0:
                    plus_values.append(value)
                else:
                    minus_values.append(value)
        llr[idx] = float(logsumexp(plus_values) - logsumexp(minus_values))

    return BcjrResult(llr=llr, log_alpha=alpha_before, log_beta=beta_after)


def brute_force_bpsk_llr(
    y: np.ndarray,
    g: dict[int, float],
    n0: float,
    la: np.ndarray | None = None,
    isi_len: int | None = None,
    initial_state: tuple[float, ...] | None = None,
) -> np.ndarray:
    """Brute-force MAP LLR by enumerating all BPSK symbol sequences."""
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.size
    if isi_len is None:
        isi_len = max((abs(int(k)) for k in g.keys()), default=0)
    if initial_state is None:
        initial_state = tuple(1.0 for _ in range(isi_len))
    if len(initial_state) != isi_len:
        raise ValueError("initial_state length must equal isi_len.")
    la_arr = _prepare_la(la, n)

    seqs: list[tuple[float, ...]] = []
    metrics: list[float] = []
    for seq in product(BPSK_ALPHABET, repeat=n):
        state = tuple(initial_state)
        metric = 0.0
        for idx, symbol in enumerate(seq):
            metric += log_gamma_bpsk(y[idx], symbol, state, g, n0, la_arr[idx])
            state = _shift_state(state, symbol)
        seqs.append(tuple(float(v) for v in seq))
        metrics.append(metric)

    metrics_arr = np.asarray(metrics, dtype=float)
    llr = np.zeros(n, dtype=float)
    for idx in range(n):
        plus = [metrics_arr[pos] for pos, seq in enumerate(seqs) if seq[idx] > 0]
        minus = [metrics_arr[pos] for pos, seq in enumerate(seqs) if seq[idx] < 0]
        llr[idx] = float(logsumexp(plus) - logsumexp(minus))
    return llr
