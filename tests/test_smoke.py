from pathlib import Path

import numpy as np

from scripts.run_smoke import run_uncoded_smoke


def test_uncoded_smoke_returns_finite_ber_table_and_writes_csv():
    output_dir = "results/_test_smoke"
    rows = run_uncoded_smoke(
        output_dir=output_dir,
        seed=7,
        frame_symbols=32,
        frames_per_snr=3,
        snr_db_values=(0.0, 3.0),
        tau=0.8,
        isi_len=2,
        m_states=4,
        future_len=3,
    )

    assert len(rows) == 2
    assert Path(output_dir, "uncoded_smoke.csv").exists()
    for row in rows:
        assert set(row) >= {"snr_db", "bit_errors", "total_bits", "ber"}
        assert row["total_bits"] == 96
        assert np.isfinite(row["ber"])
        assert 0.0 <= row["ber"] <= 1.0
