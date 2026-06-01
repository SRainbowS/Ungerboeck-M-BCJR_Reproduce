from __future__ import annotations

import numpy as np


def bpsk_modulate(bits: np.ndarray) -> np.ndarray:
    """Map bits to unit-energy BPSK symbols with ``0 -> +1`` and ``1 -> -1``."""
    bit_array = np.asarray(bits, dtype=np.uint8).reshape(-1)
    if np.any((bit_array != 0) & (bit_array != 1)):
        raise ValueError("BPSK bits must be 0 or 1.")
    return 1.0 - 2.0 * bit_array.astype(float)


def bpsk_hard_bits_from_llr(llr: np.ndarray) -> np.ndarray:
    """Hard decisions for LLR ``log P(+1) / P(-1)``."""
    return (np.asarray(llr, dtype=float).reshape(-1) < 0.0).astype(np.uint8)


# ---------------------------------------------------------------------------
# 16-QAM (Gray-mapped, unit average energy)
# ---------------------------------------------------------------------------

_QAM16_SYMBOLS: np.ndarray | None = None
_QAM16_BIT_MAP: np.ndarray | None = None
_QAM16_BITS_TO_INDEX: np.ndarray | None = None


def _init_qam16() -> tuple[np.ndarray, np.ndarray]:
    global _QAM16_SYMBOLS, _QAM16_BIT_MAP, _QAM16_BITS_TO_INDEX
    if _QAM16_SYMBOLS is not None:
        return _QAM16_SYMBOLS, _QAM16_BIT_MAP
    scale = 1.0 / np.sqrt(10.0)
    # Standard Gray mapping for 16-QAM
    # Symbol order: I={-3,-1,1,3}, Q={-3,-1,1,3}
    # Gray: {-3→00, -1→01, 1→11, 3→10} per axis
    bits_arr = np.array([
        [0, 0, 0, 0],  # (-3-3j)/√10
        [0, 0, 0, 1],  # (-3-1j)/√10
        [0, 0, 1, 1],  # (-3+1j)/√10
        [0, 0, 1, 0],  # (-3+3j)/√10
        [0, 1, 0, 0],  # (-1-3j)/√10
        [0, 1, 0, 1],  # (-1-1j)/√10
        [0, 1, 1, 1],  # (-1+1j)/√10
        [0, 1, 1, 0],  # (-1+3j)/√10
        [1, 1, 0, 0],  # (1-3j)/√10
        [1, 1, 0, 1],  # (1-1j)/√10
        [1, 1, 1, 1],  # (1+1j)/√10
        [1, 1, 1, 0],  # (1+3j)/√10
        [1, 0, 0, 0],  # (3-3j)/√10
        [1, 0, 0, 1],  # (3-1j)/√10
        [1, 0, 1, 1],  # (3+1j)/√10
        [1, 0, 1, 0],  # (3+3j)/√10
    ], dtype=np.uint8)
    symbols = np.array([
        -3 - 3j, -3 - 1j, -3 + 1j, -3 + 3j,
        -1 - 3j, -1 - 1j, -1 + 1j, -1 + 3j,
         1 - 3j,  1 - 1j,  1 + 1j,  1 + 3j,
         3 - 3j,  3 - 1j,  3 + 1j,  3 + 3j,
    ], dtype=complex) * scale
    # Build bits→index lookup: 4-bit natural binary → Gray symbol index
    bits_to_index = np.zeros(16, dtype=np.int64)
    for idx in range(16):
        key = int(bits_arr[idx, 0]) << 3 | int(bits_arr[idx, 1]) << 2 | int(bits_arr[idx, 2]) << 1 | int(bits_arr[idx, 3])
        bits_to_index[key] = idx
    _QAM16_SYMBOLS = symbols
    _QAM16_BIT_MAP = bits_arr
    _QAM16_BITS_TO_INDEX = bits_to_index
    return symbols, bits_arr


def qam16_constellation() -> tuple[np.ndarray, np.ndarray]:
    """Return ``(symbols[16], bit_map[16,4])`` for Gray-mapped 16-QAM."""
    return _init_qam16()


def qam16_modulate(bits: np.ndarray) -> np.ndarray:
    """Map bits to 16-QAM complex symbols (4 bits per symbol, unit average energy)."""
    b = np.asarray(bits, dtype=np.uint8).reshape(-1)
    if b.size % 4 != 0:
        raise ValueError("Number of bits must be a multiple of 4.")
    _init_qam16()
    symbols = _QAM16_SYMBOLS
    lut = _QAM16_BITS_TO_INDEX
    n_sym = b.size // 4
    indices = np.empty(n_sym, dtype=np.int64)
    for k in range(n_sym):
        key = int(b[4 * k]) << 3 | int(b[4 * k + 1]) << 2 | int(b[4 * k + 2]) << 1 | int(b[4 * k + 3])
        indices[k] = lut[key]
    return symbols[indices]


def qam16_demod_bit_llr(log_sym_post: np.ndarray) -> np.ndarray:
    """Convert symbol log-posteriors ``[n, 16]`` to per-bit LLR ``[n*4]``.

    LLR convention: ``log P(bit=0) / P(bit=1)``.
    """
    lp = np.asarray(log_sym_post, dtype=float)
    if lp.ndim != 2 or lp.shape[1] != 16:
        raise ValueError("log_sym_post must have shape (n, 16).")
    _, bit_map = _init_qam16()
    n = lp.shape[0]
    bit_llr = np.zeros(n * 4, dtype=float)
    for m in range(4):
        mask0 = (bit_map[:, m] == 0)
        mask1 = ~mask0
        s0 = _logsumexp_rows(lp[:, mask0])
        s1 = _logsumexp_rows(lp[:, mask1])
        bit_llr[m::4] = s0 - s1
    return bit_llr


def _logsumexp_rows(a: np.ndarray) -> np.ndarray:
    mx = np.max(a, axis=1, keepdims=True)
    return np.squeeze(mx, axis=1) + np.log(np.sum(np.exp(a - mx), axis=1))
