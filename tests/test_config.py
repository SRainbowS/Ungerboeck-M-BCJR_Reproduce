from dataclasses import asdict

import pytest

from ftn.config import FtnConfig


def test_ftn_config_defaults_match_first_phase_reproduction_scope():
    config = FtnConfig()

    assert config.rolloff == 0.3
    assert config.pulse_span == 15
    assert config.tau == 0.5
    assert config.m_states == 2
    assert config.future_len == 5
    assert config.llr_clip == 20.0
    assert asdict(config)["turbo_iters"] == 5


def test_ftn_config_rejects_invalid_values():
    with pytest.raises(ValueError):
        FtnConfig(tau=0.0)
    with pytest.raises(ValueError):
        FtnConfig(isi_len=0)
    with pytest.raises(ValueError):
        FtnConfig(m_states=0)
