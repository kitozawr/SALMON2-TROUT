#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Self-contained regression test for dft_to_epm (no SALMON build required).

We synthesise an EPM band structure with *known* Si form factors, hand those
bands to the fitter as if they were a DFT result, and check that the fit
recovers the input form factors (and a zero rigid shift) to tight tolerance.
This exercises the forward model + the least-squares inversion for BOTH the
cubic and primitive forward models.

Run:  python3 tools/dft_to_epm/tests/test_recovery.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dft_to_epm import EPMModel, run_fit, run_fit_zunger  # noqa: E402
from deep_pp_adapter import backend_name  # noqa: E402

RY = 0.5  # Ry -> Ha
SI_VS_RY = {3: -0.2258, 8: 0.05698, 11: 0.070709}  # Kunikiyo
A_LATTICE = 10.26
CUTOFF = 11.1
SHELLS = [3, 8, 11]


def random_kpath(model, nk, seed=0):
    """A handful of k-points inside the BZ (reciprocal-lattice fractional)."""
    rng = np.random.default_rng(seed)
    if model.cell == 'cubic':
        b = (2 * np.pi / A_LATTICE) * np.eye(3)
    else:
        b = model.bmat
    frac = rng.uniform(-0.5, 0.5, size=(nk, 3))
    return frac @ b


def check_cell(cell, tol=1e-6):
    model = EPMModel(cell, A_LATTICE, CUTOFF, SHELLS, [])
    kpts = random_kpath(model, nk=12, seed=42)
    nb_fit = 16
    # synthetic "DFT" bands from the known factors (Hartree)
    e_dft = model.bands(kpts, SI_VS_RY, {}, nb_fit)

    vs, va, delta, res, sol = run_fit(model, kpts, e_dft, SHELLS, [], nb_fit)
    rms = float(np.sqrt(np.mean(res ** 2)))

    errs = {s: abs(vs[s] - SI_VS_RY[s]) for s in SHELLS}
    max_ff_err = max(errs.values())
    print("[%s] recovered V^S(Ry): " % cell
          + ", ".join("%d:%.6f" % (s, vs[s]) for s in SHELLS))
    print("[%s] form-factor max error = %.2e Ry, delta = %.2e Ha, "
          "band RMS = %.2e Ha" % (cell, max_ff_err, delta, rms))

    assert max_ff_err < tol, "%s: form-factor error %.2e exceeds %.2e" % (
        cell, max_ff_err, tol)
    assert abs(delta) < tol, "%s: spurious shift %.2e" % (cell, delta)
    assert rms < tol, "%s: band RMS %.2e too large" % (cell, rms)


def check_zunger(cell, tol=1e-5):
    """Method 'zunger' (vendored DeePseudopot analytic form) must also reproduce
    the known Si shell factors via the band fit."""
    model = EPMModel(cell, A_LATTICE, CUTOFF, SHELLS, [])
    kpts = random_kpath(model, nk=12, seed=7)
    nb_fit = 16
    e_dft = model.bands(kpts, SI_VS_RY, {}, nb_fit)

    vs, va, delta, res, sol, zp = run_fit_zunger(
        model, kpts, e_dft, SHELLS, [], nb_fit, nspecies=1)
    rms = float(np.sqrt(np.mean(res ** 2)))
    max_ff_err = max(abs(vs[s] - SI_VS_RY[s]) for s in SHELLS)
    print("[%s/zunger backend=%s] V^S(Ry): %s" % (
        cell, backend_name(),
        ", ".join("%d:%.6f" % (s, vs[s]) for s in SHELLS)))
    print("[%s/zunger] form-factor max error = %.2e Ry, band RMS = %.2e Ha"
          % (cell, max_ff_err, rms))
    assert max_ff_err < tol, "%s/zunger: ff error %.2e" % (cell, max_ff_err)
    assert rms < tol, "%s/zunger: band RMS %.2e" % (cell, rms)


def main():
    for cell in ('cubic', 'primitive'):
        check_cell(cell)
    check_zunger('cubic')
    print("\nALL TESTS PASSED")
    return 0


if __name__ == '__main__':
    sys.exit(main())
