#!/usr/bin/env python3
"""
test_si_epm_gap.py  -  Part A (Silicon EPM) validation.

Validates the Silicon form factors by computing the indirect gap and the CBM
position from the PRIMITIVE sublattice-block bands of the local EPM (the same
form factors and structure factor the Fortran `theory='epm'` Si path uses).

Targets (see wiki/02_constants.md, wiki/00_implementation_status.md):
  - converged indirect gap ~1.06 eV  (Kunikiyo's own calc 1.068 eV; the 1.12 eV
    experimental value is intentionally NOT reached by a 3-parameter local EPM)
  - CBM near X at ~0.85*(2pi/a) along <100>  (six Delta-valleys)

No SALMON build required; reuses epm_gaas_reference.py.
"""
import importlib.util
import os
import sys

import numpy as np
from numpy.linalg import eigvalsh

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_epm():
    spec = importlib.util.spec_from_file_location(
        "epm", os.path.join(ROOT, "epm_gaas_reference.py"))
    epm = importlib.util.module_from_spec(spec)
    saved = sys.argv
    sys.argv = ["epm"]
    try:
        spec.loader.exec_module(epm)
    finally:
        sys.argv = saved
    return epm


def si_indirect_gap(epm, a=10.26, cutoff=27.0, material="Si"):
    """Primitive-band indirect gap [eV] and CBM position [2pi/a] along <100>."""
    Gcart, _ = epm.build_plane_wave_basis_sc(a, cutoff)
    tpa = 2.0 * np.pi / a
    Gi = np.round(Gcart / tpa).astype(int)
    nv = 4  # scalar Si: 8 valence e- per primitive cell / 2 per band

    def prim_bands(qcart):
        G0 = np.round(qcart / tpa).astype(int)
        ksc = qcart - G0 * tpa
        idx = np.where(epm.sublattice_mask(Gi, G0))[0]
        H = epm.build_hamiltonian_sc(material, ksc, Gcart, a)
        return np.sort(eigvalsh(H[np.ix_(idx, idx)]))

    vbm = prim_bands(np.zeros(3))[nv - 1]
    ts = np.linspace(0.0, 1.0, 101)
    cb = np.array([prim_bands(t * tpa * np.array([1.0, 0.0, 0.0]))[nv] for t in ts])
    icbm = int(cb.argmin())
    gap_ev = (cb[icbm] - vbm) * epm.HARTREE_EV
    return gap_ev, ts[icbm]


def main():
    epm = _load_epm()
    ok = True

    # Converged basis: gap should match Kunikiyo's own calc (1.068 eV) closely.
    gap, tcbm = si_indirect_gap(epm, cutoff=27.0)
    print(f"[Si] converged (cutoff 27 Ry): E_g(indirect) = {gap:.4f} eV, "
          f"CBM @ {tcbm:.3f}*(2pi/a) along <100>")
    if not (abs(gap - 1.068) < 0.030):
        print(f"  FAIL: gap {gap:.4f} eV not within 30 meV of Kunikiyo calc 1.068 eV")
        ok = False
    if not (0.80 <= tcbm <= 0.92):
        print(f"  FAIL: CBM at {tcbm:.3f} not near the 0.85 Delta-valley target")
        ok = False

    # V^A must be exactly zero for Silicon (diamond) on every shell.
    for g2 in (3, 4, 8, 11):
        _vs, va = epm.form_factors("Si", g2)
        if va != 0.0:
            print(f"  FAIL: Si V^A({g2}) = {va} != 0 (diamond must have V^A=0)")
            ok = False

    # GaAs must be untouched (regression guard).
    vs3, va3 = epm.form_factors("GaAs", 3)
    if not (abs(vs3 - (-0.23 * 0.5)) < 1e-12 and abs(va3 - (0.07 * 0.5)) < 1e-12):
        print(f"  FAIL: GaAs form factors changed: V^S(3)={vs3}, V^A(3)={va3}")
        ok = False

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
