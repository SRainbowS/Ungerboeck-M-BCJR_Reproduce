import numpy as np

from ftn.channel import ftn_filter_output
from ftn.equalizers.full_bcjr import full_bcjr_bpsk
from ftn.equalizers.ungerboeck_mbcjr import ungerboeck_mbcjr_bpsk


def test_mbcjr_matches_full_bcjr_when_all_states_and_full_future_are_kept():
    g = {0: 1.0, 1: 0.28, -1: 0.28, 2: 0.07, -2: 0.07}
    x_true = np.array([1.0, -1.0, 1.0, -1.0])
    y = ftn_filter_output(x_true, g) + np.array([0.03, -0.04, 0.06, -0.02])
    la = np.array([0.0, 0.15, -0.2, 0.1])
    n0 = 0.9

    full = full_bcjr_bpsk(y, g, n0=n0, la=la, isi_len=2, initial_state=(1.0, 1.0))
    mbcjr = ungerboeck_mbcjr_bpsk(
        y,
        g,
        n0=n0,
        la=la,
        isi_len=2,
        m_states=4,
        future_len=len(y),
        initial_state=(1.0, 1.0),
    )

    np.testing.assert_allclose(mbcjr.llr, full.llr, atol=1e-10)


def test_mbcjr_never_keeps_more_than_requested_survivors():
    g = {0: 1.0, 1: 0.28, -1: 0.28, 2: 0.07, -2: 0.07}
    y = np.array([0.6, -0.2, 0.1, -0.7, 0.5])

    result = ungerboeck_mbcjr_bpsk(
        y,
        g,
        n0=1.0,
        isi_len=2,
        m_states=2,
        future_len=2,
        initial_state=(1.0, 1.0),
    )

    assert max(result.survivor_counts) <= 2
