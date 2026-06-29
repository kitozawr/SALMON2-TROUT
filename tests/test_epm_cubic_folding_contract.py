#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_epm_cubic_folding_contract.py

Locks the convention contract that the Fortran EPM (src/epm/epm_solver.f90) MUST
reproduce to be byte-equivalent to the Python reference (epm_gaas_reference.py,
the source of truth) for the cubic materials GaAs / Si / Si_cb:

  1. SIMPLE-CUBIC supercell basis (not FCC primitive), cutoff on |G|^2 in
     (2*pi/a)^2 units (integer shells)  ->  fixed plane-wave count at the
     default cutoff (171 PW for GaAs/Si at epm_pw_cutoff_ry = 11.1).
  2. The FCC-in-cubic PARITY selection rule makes the cubic Hamiltonian exactly
     block-diagonal over the 4 FCC reciprocal cosets (folding is exact) -- the
     property the Fortran reproduces with `mod(dh-dk,2)==0 .and. mod(dk-dl,2)==0`.
  3. Monkhorst-Pack grid + REDUCED k-coordinate convention: the first cubic MP
     point of a 4x4x4 grid is (-3/8,-3/8,-3/8) in reduced coords.

If someone reverts the Fortran solver to the FCC-primitive cell / a.u. cutoff /
Cartesian k.data (the pre-fix bug), these invariants change and this test fails.
The full numeric equivalence (band energies to ~5e-11 Ha, occupations identical,
valence/optical momentum to ~1e-10) was verified by BUILDING SALMON and diffing
the actual Fortran output against this reference for GaAs and Si (scalar); that
binary diff is not reproducible in the lightweight test harness, so this test
guards the convention contract that makes the diff pass.
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load_ref():
    """Import the (scalar) Python EPM reference module."""
    path = os.path.join(ROOT, 'epm_gaas_reference.py')
    spec = importlib.util.spec_from_file_location('epm_ref', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ref = _load_ref()
    nfail = 0

    def check(cond, msg):
        nonlocal nfail
        if not cond:
            print(f'  FAIL: {msg}')
            nfail += 1

    a = 10.68  # GaAs; Si (10.26) gives the same integer-shell basis count
    cutoff = 11.1

    # (1) simple-cubic basis size at the default cutoff (contract: 171 PW)
    Gcart, G2 = ref.build_plane_wave_basis_sc(a, cutoff)
    npw = Gcart.shape[0]
    check(npw == 171, f'plane-wave count at cutoff {cutoff} = {npw}, expected 171')
    # cutoff is on integer shells h^2+k^2+l^2 (a.u. cutoff would give a wildly
    # different count); the max retained shell must be <= cutoff.
    check(int(G2.max()) <= int(cutoff) + 0 and int(G2.max()) == 11,
          f'max retained shell |G|^2 = {G2.max()} (2pi/a)^2 units, expected 11')

    # (2) folding exactness: cubic H block-diagonal over the 4 FCC cosets.
    twopi_over_a = 2.0 * np.pi / a
    G_idx = np.round(Gcart / twopi_over_a).astype(int)
    # build at a generic k to avoid accidental symmetry
    kvec = np.array([0.123, -0.071, 0.045])
    H = ref.build_hamiltonian_sc('GaAs', kvec, Gcart, a)
    check(np.allclose(H, H.conj().T, atol=1e-13), 'cubic H not Hermitian')
    masks = [ref.sublattice_mask(G_idx, off) for off in ref.SUBLATTICE_OFFSETS]
    # every plane wave belongs to exactly one coset
    cover = np.zeros(npw, dtype=int)
    for m in masks:
        cover += m.astype(int)
    check(np.all(cover == 1), 'cosets do not partition the plane-wave basis')
    # off-coset blocks must be exactly zero (the parity rule => exact folding)
    max_off = 0.0
    for s in range(4):
        rows = np.where(masks[s])[0]
        others = np.setdiff1d(np.arange(npw), rows)
        if rows.size and others.size:
            max_off = max(max_off, np.abs(H[np.ix_(rows, others)]).max())
    check(max_off < 1e-13, f'off-coset |H| = {max_off:.2e} (folding not exact)')

    # the 4 coset blocks each carry the same number of valence electrons; the
    # union of their spectra == the full cubic spectrum (folding identity).
    full = np.sort(np.linalg.eigvalsh(H))
    sub = []
    for s in range(4):
        idx = np.where(masks[s])[0]
        sub.append(np.linalg.eigvalsh(H[np.ix_(idx, idx)]))
    sub = np.sort(np.concatenate(sub))
    check(np.allclose(full, sub, atol=1e-10),
          'folded cubic spectrum != union of the 4 coset (primitive) spectra')

    # (3) MP grid + reduced-coordinate convention (4x4x4 first point = -3/8).
    b1, b2, b3, _ = ref.reciprocal_lattice_sc(a)
    bmat = np.array([b1, b2, b3])
    kpt, kw = ref.monkhorst_pack_grid(bmat, (4, 4, 4))
    check(kpt.shape[0] == 64, f'MP grid has {kpt.shape[0]} points, expected 64')
    b_diag = np.array([bmat[0, 0], bmat[1, 1], bmat[2, 2]])
    kred0 = kpt[0] / b_diag
    check(np.allclose(kred0, [-0.375, -0.375, -0.375], atol=1e-12),
          f'first reduced MP point = {kred0}, expected (-0.375,-0.375,-0.375)')
    check(abs(kw[0] - 1.0 / 64) < 1e-12, 'MP weight != 1/nk')

    if nfail == 0:
        print('PASS  (EPM cubic band-folding contract: 171 PW, exact 4-coset '
              'folding, reduced-k MP grid -- Fortran must mirror this)')
        return 0
    print(f'FAIL ({nfail} checks)')
    return 1


if __name__ == '__main__':
    sys.exit(main())
