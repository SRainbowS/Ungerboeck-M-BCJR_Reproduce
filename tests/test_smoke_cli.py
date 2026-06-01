import subprocess
import sys
from pathlib import Path


def test_run_smoke_script_works_when_executed_by_path():
    output_dir = Path("results/_test_smoke_cli")
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/run_smoke.py",
            "--output-dir",
            str(output_dir),
            "--frame-symbols",
            "16",
            "--frames-per-snr",
            "1",
            "--snr-db",
            "0",
        ],
        cwd=Path(__file__).resolve().parents[1],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert (output_dir / "uncoded_smoke.csv").exists()
