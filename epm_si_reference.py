#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone Silicon local-EPM reference (theory='epm' Si path).

Silicon is the diamond structure: two IDENTICAL carbon-group atoms per primitive
cell, so the ANTISYMMETRIC structure factor vanishes and V^A == 0 exactly for
every shell (same simplification as graphene). Only the symmetric V^S enters.
The cubic 8-atom supercell + 4-fold FCC band folding / exact unfolding machinery
is shared with the GaAs reference (epm_gaas_reference.py); this module just
configures it for Si and exposes a Si-only entry point so Si can be generated
and validated independently of GaAs.

Two Si form-factor sets (select with `--variant`):
  * 'Si'    (DEFAULT, Kunikiyo): V^S(3)=-0.2258, V^S(8)=+0.05698,
            V^S(11)=+0.070709 Ry. [T. Kunikiyo et al., J. Appl. Phys. 75, 297
            (1994), Table I]. Converged indirect gap ~1.06 eV (= Kunikiyo's own
            calc 1.068 eV; the 1.12 eV experimental value is intentionally NOT
            reached by a 3-parameter local EPM), CBM near X at ~0.85*(2pi/a).
  * 'Si_cb' (Cohen-Bergstresser alt): V^S(3)=-0.21, V^S(8)=+0.04, V^S(11)=+0.08
            Ry. [M. L. Cohen & T. K. Bergstresser, Phys. Rev. 141, 789 (1966)].
            Provided for cross-validation; gives a slightly different gap.
Both sets have V^A(3)=V^A(4)=V^A(11)=0 exactly (diamond). The ONLY difference
between 'Si' and 'Si_cb' is which V^S triplet is used; everything else (lattice
10.26 Bohr, folding, structure factor) is identical.

Lattice constant a = 10.26 Bohr (5.431 Ang). Scalar (no spin-orbit) by default.

Usage:
    python3 epm_si_reference.py                 # Si (Kunikiyo), writes Si GS files
    python3 epm_si_reference.py --variant Si_cb # Cohen-Bergstresser variant
    python3 epm_si_reference.py --gap           # just report the indirect gap

The heavy lifting (plane-wave basis, secular equation, folding/unfold, the
eigen/tm/k.data writers the SBE reads) is imported from epm_gaas_reference.py --
that validated machinery is reused verbatim, not duplicated.
"""

import argparse
import importlib.util
import os
import sys

import numpy as np
from numpy.linalg import eigvalsh

ROOT = os.path.dirname(os.path.abspath(__file__))
SI_A_BOHR = 10.26          # Si lattice constant [Bohr] = 5.431 Ang


def load_epm_machinery():
    """Import the shared local-EPM machinery (epm_gaas_reference.py) without
    running its __main__."""
    spec = importlib.util.spec_from_file_location(
        "epm_gaas_reference", os.path.join(ROOT, "epm_gaas_reference.py"))
    epm = importlib.util.module_from_spec(spec)
    saved = sys.argv
    sys.argv = ["epm_gaas_reference"]
    try:
        spec.loader.exec_module(epm)
    finally:
        sys.argv = saved
    return epm


def configure_for_si(epm, variant="Si", cutoff_ry=11.1, num_kgrid=(4, 4, 4)):
    """Point the shared machinery at Silicon (scalar, diamond, V^A=0)."""
    assert variant in ("Si", "Si_cb"), "variant must be 'Si' (Kunikiyo) or 'Si_cb'"
    epm.MATERIAL = variant
    epm.A_LATTICE_AU = SI_A_BOHR
    epm.PW_CUTOFF_RY = cutoff_ry
    epm.NUM_KGRID = list(num_kgrid)
    epm.INCLUDE_SPIN_ORBIT = False
    epm.SYSNAME = "Si_cubic" if variant == "Si" else "Si_cb_cubic"
    return epm


def si_indirect_gap(epm, variant="Si", cutoff=27.0):
    """Primitive sublattice-block indirect gap [eV] and CBM position [2pi/a]
    along <100> (the same quantity tests/test_si_epm_gap.py checks)."""
    twopi = 2.0 * np.pi / SI_A_BOHR
    Gcart, _ = epm.build_plane_wave_basis_sc(SI_A_BOHR, cutoff)
    G_idx = np.round(Gcart / twopi).astype(int)
    # primitive Gamma-block (FCC sublattice g0=0): exact-folding subblock
    g0 = np.array([0, 0, 0])
    msk = epm.sublattice_mask(G_idx, g0)
    rows = np.where(msk)[0]

    def bands(kred):
        kcart = np.array(kred) * twopi
        H = epm.build_hamiltonian_sc(variant, kcart, Gcart, SI_A_BOHR)
        return eigvalsh(H[np.ix_(rows, rows)])

    nv = 4                                   # 8 valence e / 2 (spinless) per prim cell
    vbm = sorted(bands([0, 0, 0]))[nv - 1]   # VBM = top valence band at Gamma
    best_cb, best_x = 1e9, 0.0
    for f in np.linspace(0.7, 1.0, 31):      # scan Gamma->X along <100>
        cb = sorted(bands([f, 0, 0]))[nv]    # lowest conduction band
        if cb < best_cb:
            best_cb, best_x = cb, f
    return (best_cb - vbm) * 27.211386245988, best_x


def main():
    ap = argparse.ArgumentParser(description="Standalone Silicon local-EPM reference.")
    ap.add_argument("--variant", choices=["Si", "Si_cb"], default="Si",
                    help="Si (Kunikiyo, default) or Si_cb (Cohen-Bergstresser)")
    ap.add_argument("--cutoff-ry", type=float, default=11.1)
    ap.add_argument("--gap", action="store_true",
                    help="only report the converged indirect gap, do not write files")
    args = ap.parse_args()

    epm = load_epm_machinery()
    configure_for_si(epm, args.variant, cutoff_ry=args.cutoff_ry)

    gap, xcbm = si_indirect_gap(epm, args.variant)
    print(f"[{args.variant}] indirect gap = {gap:.4f} eV, "
          f"CBM @ {xcbm:.3f}*(2pi/a) along <100>  "
          f"(target: Si ~1.06 eV / Kunikiyo 1.068; CBM ~0.85)")
    if args.gap:
        return
    epm.main()       # write the Si GS data files (eigen/tm/k + unfold map)


if __name__ == "__main__":
    main()
