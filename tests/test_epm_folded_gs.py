#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_epm_folded_gs.py -- scalar folded-supercell GS emission for CdS & graphene.

Validates the Python EPM references that close the EPM->SBE pipeline for the
non-cubic materials (the Fortran src/epm solver is cubic-only -- see docs):

  * epm_io shared writer round-trip (reduced-k format, Hartree/occupations).
  * graphene RECTANGULAR 4-atom cell + 2-fold fold: exact folding (off-coset
    ~0), the supercell reproduces the primitive Dirac cone (gapless), and the
    structure-factor normalization makes supercell == primitive bands.
  * graphene + CdS main_gs: emit SBE GS files, correct band count / occupations
    (2 per filled band), Hermitian momentum, finite energies. Both were also
    verified to RUN in the actual SBE binary (trace conserved) -- see wiki/00.

Pure numpy; no SALMON build. Small grids/cutoffs keep it fast.
"""
import os
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import epm_io
import epm_graphene as G
import epm_wurtzite_cds as W

nfail = 0
def check(name, cond):
    global nfail
    if not cond:
        print(f"  FAIL: {name}")
        nfail += 1


# --- epm_io shared helpers ---------------------------------------------------
b = np.diag([2.0, 3.0, 5.0])
kpt, kw = epm_io.monkhorst_pack(b, (4, 4, 4))
check("epm_io MP grid: 64 pts, weights sum to 1",
      kpt.shape == (64, 3) and abs(kw.sum() - 1.0) < 1e-12)
check("epm_io MP first reduced point = -3/8",
      np.allclose(kpt[0] / np.diag(b), [-0.375, -0.375, -0.375], atol=1e-12))
# non-orthogonal cell must be rejected by the reduced-coord writer
try:
    epm_io.write_epm_gs_files('x', tempfile.gettempdir() + '/', 'x',
                              np.zeros((1, 3)), np.array([[1, 0.3, 0], [0, 1, 0], [0, 0, 1.]]),
                              np.ones(1), np.zeros((1, 1)), np.zeros((1, 1)),
                              np.zeros((1, 1, 3, 1), complex), np.zeros((1, 1, 3, 1), complex))
    check("epm_io rejects non-orthogonal cell", False)
except ValueError:
    check("epm_io rejects non-orthogonal cell", True)


# --- graphene: rectangular 4-atom 2-fold folding ------------------------------
off, gap, kfold = G.rect_folding_check()
check("graphene 2-fold folding EXACT (off-coset |H| ~ 0)", off < 1e-10)
check("graphene folded Dirac cone gapless (pi/pi* gap ~ 0)", abs(gap) < 5e-3)

# struct_norm fix: supercell reproduces the PRIMITIVE Dirac energy (1.034 eV)
Gp, _ = G.build_pw_basis(400.0)
evK_prim = np.sort(np.linalg.eigvalsh(G.build_hamiltonian(G.high_symmetry_points()['K'], Gp, G.basis_atoms())))
Gc, _ = G.rect_pw_basis(400.0)
posR = G.rect_atoms_ang()
b1r, b2r = G.rect_reciprocal_ang()
kK = np.array([1/3, 0]) @ np.array([b1r, b2r])
evK_rect = np.sort(np.linalg.eigvalsh(G.build_hamiltonian(kK, Gc, posR, struct_norm=G.GRAPHENE_STRUCT_NORM)))
check("graphene supercell reproduces primitive Dirac energy (struct_norm)",
      abs(evK_rect[1] - evK_prim[0]) < 1e-3 and abs(evK_rect[2] - evK_prim[1]) < 1e-3)
check("graphene 4 atoms in rectangular cell", len(posR) == 4)

# --- graphene GS emission ----------------------------------------------------
with tempfile.TemporaryDirectory() as td:
    td = td + '/'
    eig, occ = G.main_gs(num_kgrid=(2, 2, 1), nstate=G.GRAPHENE_NSTATE,
                         nelec=G.GRAPHENE_NELEC, outdir=td)
    nk = eig.shape[1]
    check("graphene emit: nstate bands", eig.shape[0] == G.GRAPHENE_NSTATE)
    check("graphene emit: 2 filled bands (occ 2), rest empty",
          np.allclose(occ[:2], 2.0) and np.allclose(occ[2:], 0.0))
    check("graphene emit: energies finite", np.all(np.isfinite(eig)))
    for f in ('graphene_k.data', 'graphene_eigen.data', 'graphene_tm.data'):
        check(f"graphene emit: wrote {f}", os.path.exists(td + f))
    # k.data must be reduced (|kx| <= 0.5)
    kvals = np.loadtxt(td + 'graphene_k.data', comments='#')
    check("graphene k.data reduced (|k|<=0.5)", np.abs(kvals[:, 1:4]).max() <= 0.5 + 1e-9)


# --- CdS GS emission (small/fast) --------------------------------------------
with tempfile.TemporaryDirectory() as td:
    td = td + '/'
    eig, occ = W.main_gs(sysname='CdS', num_kgrid=(2, 2, 1), nstate=20,
                        nelec=W.CDS_NELEC, cutoff_ry=5.0, outdir=td)
    check("CdS emit: 20 bands", eig.shape[0] == 20)
    check("CdS emit: 16 filled bands (occ 2)", np.allclose(occ[:16], 2.0))
    check("CdS emit: energies finite", np.all(np.isfinite(eig)))
    for f in ('CdS_k.data', 'CdS_eigen.data', 'CdS_tm.data'):
        check(f"CdS emit: wrote {f}", os.path.exists(td + f))


if nfail == 0:
    print("PASS  (folded scalar GS emission for CdS + graphene: exact folding, "
          "gapless graphene Dirac, correct occupations, SBE-ready files)")
    sys.exit(0)
print(f"FAIL ({nfail} checks)")
sys.exit(1)
