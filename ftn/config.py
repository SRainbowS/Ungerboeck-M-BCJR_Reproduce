from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FtnConfig:
    """Configuration for first-phase FTN Ungerboeck M-BCJR reproduction."""

    tau: float = 0.5
    rolloff: float = 0.3
    pulse_span: int = 15
    sps: int = 128
    isi_len: int = 3
    future_len: int = 5
    m_states: int = 2
    snr_db: float = 4.0
    frame_bits: int = 6000
    seed: int = 20260506
    turbo_iters: int = 5
    llr_clip: float = 20.0

    def __post_init__(self) -> None:
        if not 0.0 < self.tau <= 1.0:
            raise ValueError("tau must be in (0, 1].")
        if not 0.0 <= self.rolloff <= 1.0:
            raise ValueError("rolloff must be in [0, 1].")
        if self.pulse_span <= 0:
            raise ValueError("pulse_span must be positive.")
        if self.sps <= 0:
            raise ValueError("sps must be positive.")
        if self.isi_len <= 0:
            raise ValueError("isi_len must be positive.")
        if self.future_len < 0:
            raise ValueError("future_len must be non-negative.")
        if self.m_states <= 0:
            raise ValueError("m_states must be positive.")
        if self.frame_bits <= 0:
            raise ValueError("frame_bits must be positive.")
        if self.turbo_iters < 0:
            raise ValueError("turbo_iters must be non-negative.")
        if self.llr_clip <= 0.0:
            raise ValueError("llr_clip must be positive.")

    @property
    def n0(self) -> float:
        """Use unit-energy BPSK symbols, so ``N0 = 1 / EsN0``."""
        return 1.0 / (10.0 ** (self.snr_db / 10.0))
