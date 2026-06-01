import numpy as np

from ftn.metrics import log_phi_ungerboeck, log_prior_bpsk_symbol


def test_log_prior_bpsk_symbol_uses_llr_for_plus_over_minus():
    la = 1.7

    log_p_plus = log_prior_bpsk_symbol(+1.0, la)
    log_p_minus = log_prior_bpsk_symbol(-1.0, la)

    assert np.isclose(log_p_plus - log_p_minus, la)
    assert np.isclose(np.exp(log_p_plus) + np.exp(log_p_minus), 1.0)


def test_log_phi_uses_previous_state_from_most_recent_symbol_backwards():
    g = {0: 1.0, 1: 0.25, 2: 0.10}
    state_prev = (-1.0, +1.0)
    x_n = -1.0
    y_n = 0.7
    n0 = 0.5

    expected_isi = 0.25 * state_prev[-1] + 0.10 * state_prev[-2]
    expected = (2.0 / n0) * x_n * (y_n - 0.5 * g[0] * x_n - expected_isi)

    assert np.isclose(log_phi_ungerboeck(y_n, x_n, state_prev, g, n0), expected)
