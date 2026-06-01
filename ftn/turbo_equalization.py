from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ftn.coding.conv import conv_bcjr_decode, conv_encode_75
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk
from ftn.modulation import bpsk_hard_bits_from_llr


@dataclass(frozen=True)
class TurboIteration:
    index: int
    detector_ber: float
    info_ber: float


@dataclass(frozen=True)
class TurboEqualizationResult:
    iterations: list[TurboIteration]
    final_info_llr: np.ndarray
    final_detector_llr: np.ndarray


def _hard_bits_from_llr(llr: np.ndarray) -> np.ndarray:
    return bpsk_hard_bits_from_llr(llr)


def _ber_from_llr(llr: np.ndarray, bits: np.ndarray) -> float:
    hard = _hard_bits_from_llr(llr)
    return float(np.mean(hard != np.asarray(bits, dtype=np.uint8).reshape(-1)))


def turbo_equalize_conv_bpsk(
    y: np.ndarray,
    g: dict[int, float],
    n0: float,
    info_bits: np.ndarray,
    turbo_iters: int,
    isi_len: int,
    m_states: int,
    future_len: int,
    llr_clip: float = 20.0,
    initial_state: tuple[float, ...] | None = None,
) -> TurboEqualizationResult:
    """Run detector/conv-code SISO exchange for BPSK without interleaving."""
    if turbo_iters < 0:
        raise ValueError("turbo_iters must be non-negative.")
    info = np.asarray(info_bits, dtype=np.uint8).reshape(-1)
    encoded = conv_encode_75(info)
    if len(y) != encoded.size:
        raise ValueError("y length must match encoded code-bit length.")
    if initial_state is None:
        initial_state = tuple(1.0 for _ in range(isi_len))

    detector_prior = np.zeros(encoded.size, dtype=float)
    iterations: list[TurboIteration] = []
    final_info_llr = np.zeros(info.size, dtype=float)
    final_detector_llr = np.zeros(encoded.size, dtype=float)

    for idx in range(turbo_iters + 1):
        det = ungerboeck_mbcjr_bpsk(
            y,
            g,
            n0=n0,
            la=detector_prior,
            isi_len=isi_len,
            m_states=m_states,
            future_len=future_len,
            initial_state=initial_state,
        )
        final_detector_llr = np.clip(det.llr, -llr_clip, llr_clip)
        det_ext = np.clip(final_detector_llr - detector_prior, -llr_clip, llr_clip)

        dec = conv_bcjr_decode(det_ext)
        final_info_llr = np.clip(dec.info_llr, -llr_clip, llr_clip)
        code_post = np.clip(dec.code_llr, -llr_clip, llr_clip)
        decoder_ext = np.clip(code_post - det_ext, -llr_clip, llr_clip)

        iterations.append(
            TurboIteration(
                index=idx,
                detector_ber=_ber_from_llr(final_detector_llr, encoded),
                info_ber=_ber_from_llr(final_info_llr, info),
            )
        )

        detector_prior = decoder_ext

    return TurboEqualizationResult(
        iterations=iterations,
        final_info_llr=final_info_llr,
        final_detector_llr=final_detector_llr,
    )
