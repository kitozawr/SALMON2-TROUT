#!/usr/bin/env python3
"""
run_all.py  -  run every tests/test_*.py and report PASS/FAIL.

Each test is a standalone script that prints a final 'PASS'/'FAIL' line and
exits 0 (pass) / nonzero (fail). This runner aggregates them. Exit code is the
number of failed tests (0 = all pass), so it is CI-friendly.

Usage:  python3 tests/run_all.py
"""
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    tests = sorted(glob.glob(os.path.join(HERE, "test_*.py")))
    if not tests:
        print("no tests found"); return 1
    failed = []
    for t in tests:
        name = os.path.basename(t)
        print(f"\n=== {name} ===")
        r = subprocess.run([sys.executable, t], capture_output=True, text=True)
        sys.stdout.write(r.stdout)
        if r.stderr.strip():
            sys.stdout.write(r.stderr)
        if r.returncode != 0:
            failed.append(name)
    print("\n" + "=" * 50)
    print(f"  {len(tests) - len(failed)}/{len(tests)} passed")
    if failed:
        print("  FAILED: " + ", ".join(failed))
    return len(failed)


if __name__ == "__main__":
    sys.exit(main())
