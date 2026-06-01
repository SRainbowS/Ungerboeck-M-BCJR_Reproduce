"""Run all comparison experiments in parallel.

Launches 4 parallel processes:
  1. Uncoded tau=0.5
  2. Uncoded tau=0.35
  3. Turbo  tau=0.5
  4. Turbo  tau=0.35
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path


def main():
    project = Path(__file__).resolve().parent.parent
    scripts = project / "scripts"

    jobs = [
        ("uncoded_tau05", [sys.executable, str(scripts / "comparison_uncoded.py"), "--tau", "0.5"]),
        ("uncoded_tau035", [sys.executable, str(scripts / "comparison_uncoded.py"), "--tau", "0.35"]),
        ("turbo_tau05", [sys.executable, str(scripts / "comparison_turbo.py"), "--tau", "0.5"]),
        ("turbo_tau035", [sys.executable, str(scripts / "comparison_turbo.py"), "--tau", "0.35"]),
    ]

    procs = {}
    for name, cmd in jobs:
        log = open(project / f"results/comparison/{name}.log", "w")
        print(f"Starting {name}: {' '.join(cmd)}")
        proc = subprocess.Popen(cmd, cwd=str(project), stdout=log, stderr=subprocess.STDOUT)
        procs[name] = (proc, log)

    print(f"\n{len(procs)} experiments running in parallel...\n")

    t0 = time.time()
    results = {}
    for name, (proc, log) in procs.items():
        rc = proc.wait()
        elapsed = time.time() - t0
        log.close()
        results[name] = rc
        status = "OK" if rc == 0 else f"FAILED (rc={rc})"
        print(f"  [{name}] {status}  ({elapsed:.0f}s)")

    print(f"\nTotal wall time: {time.time() - t0:.0f}s")

    # Generate plots
    for tau in [0.5, 0.35]:
        for mode in ["uncoded", "turbo"]:
            cmd = [sys.executable, str(scripts / "plot_comparison.py"), mode, "--tau", str(tau)]
            subprocess.run(cmd, cwd=str(project))

    print("Done. All experiments and plots complete.")


if __name__ == "__main__":
    main()
