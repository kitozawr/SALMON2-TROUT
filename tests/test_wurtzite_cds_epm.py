#!/usr/bin/env python3
"""
test_wurtzite_cds_epm.py  -  wurtzite CdS EPM validated against BC1967.

Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967). Checks:

  * orthorhombic (sqrt3 x 1 x 1) al(1:3) box and the hexagonal primitive cell;
  * exactly 8 atoms (orthorhombic) / 4 atoms (primitive), 2 Cd + 2 S each;
  * reciprocal-lattice consistency; H(k) Hermitian; broken inversion (complex H);
  * STRUCTURE FACTORS match BC1967 Table II at the gap shells
    (002: |S^S|,|S^A| = 0.71, 0.71; 101: 0.33, 0.80; 102: 0.35, 0.35);
  * the CITED wurtzite form factors are loaded (V^S(002 shell) = -0.26 Ry);
  * the direct gap at Gamma reproduces BC1967 Table I (2.58 eV), converged,
    via the 1/n_total normalization + the cited wurtzite Table II form factors.

Run:  python3 tests/test_wurtzite_cds_epm.py   (pure numpy, no SALMON build).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import epm_wurtzite_cds as W

nfail = 0
def check(name, cond):
    global nfail
    if not cond:
        print(f"  FAIL: {name}")
        nfail += 1

a_au = W.CDS_A_ANG * W.ANG_TO_BOHR
c_au = W.CDS_C_ANG * W.ANG_TO_BOHR

# --- orthorhombic al(1:3) cell (SBE side) ------------------------------------
cell = W.cds_cell_au()
A, B, C = W.orthorhombic_vectors_au(a_au, c_au)
check("al(1:3) = (a, a*sqrt3, c)",
      np.allclose([np.linalg.norm(A), np.linalg.norm(B), np.linalg.norm(C)], cell, atol=1e-6))
check("cell anisotropic (b = a*sqrt3)", abs(cell[1] - cell[0] * np.sqrt(3)) < 1e-6)
check("orthorhombic axes orthogonal",
      abs(A @ B) < 1e-9 and abs(B @ C) < 1e-9 and abs(A @ C) < 1e-9)
posO, specO = W.wurtzite_atoms_orth(a_au, c_au)
check("orthorhombic cell: 8 atoms (4 Cd + 4 S)",
      len(posO) == 8 and int((specO > 0).sum()) == 4 and int((specO < 0).sum()) == 4)

# --- hexagonal primitive cell (band validation) ------------------------------
a1, a2, a3 = W.hexagonal_vectors_au(a_au, c_au)
Brec, Vol = W.reciprocal(a1, a2, a3)
check("primitive reciprocal: b_i . a_j = 2pi delta",
      np.allclose(Brec @ np.array([a1, a2, a3]).T, 2 * np.pi * np.eye(3), atol=1e-9))
pos, spec = W.hex_primitive_atoms(a_au, c_au)
check("primitive cell: 4 atoms (2 Cd + 2 S)",
      len(pos) == 4 and int((spec > 0).sum()) == 2 and int((spec < 0).sum()) == 2)

# --- structure factors vs BC1967 Table II (the load-bearing #3 check) ---------
# (hkl): expected (|S^S|, |S^A|) from BC1967 Table II (wurtzite, u=3/8).
for hkl, (rs, ra) in {(0, 0, 2): (0.71, 0.71), (1, 0, 1): (0.33, 0.80),
                      (1, 0, 2): (0.35, 0.35), (0, 0, 1): (0.0, 0.0)}.items():
    G = hkl[0] * Brec[0] + hkl[1] * Brec[1] + hkl[2] * Brec[2]
    Ss, Sa = W.structure_factors(G, pos, spec)
    check(f"S^S{hkl} = {rs} (BC1967 Table II)", abs(abs(Ss) - rs) < 0.01)
    check(f"S^A{hkl} = {ra} (BC1967 Table II)", abs(abs(Sa) - ra) < 0.01)

# --- Hamiltonian on a 5x5x5 grid: Hermitian, finite, complex (broken inv.) ---
g = (np.arange(5) - 2) / 5.0
kgrid = np.array([(x, y, z) for z in g for y in g for x in g]) @ Brec
Gcart, _ = W.build_pw_basis(Brec, 9.0)
max_nonherm = max_imag = 0.0
all_finite = True
for kc in kgrid[::13]:
    H = W.build_hamiltonian(kc, Gcart, pos, spec)
    max_nonherm = max(max_nonherm, np.abs(H - H.conj().T).max())
    max_imag = max(max_imag, np.abs(H.imag).max())
    all_finite = all_finite and np.all(np.isfinite(np.linalg.eigvalsh(H)))
check("H Hermitian on the 5x5x5 grid", max_nonherm < 1e-10)
check("all eigenvalues finite", all_finite)
check("broken inversion (Im H != 0)", max_imag > 1e-6)

# --- cited wurtzite form factor loaded (002 shell, |G|^2_BC ~ 3.04) -----------
unit = (2 * np.pi / W._azb_bohr()) ** 2
vs, va = W.form_factor_phys(3.04 * unit)
check("cited wurtzite V^S(002) = -0.26 Ry", abs(vs - (-0.26)) < 1e-9)
check("cited wurtzite V^A(002) = +0.23 Ry", abs(va - (+0.23)) < 1e-9)

# --- THE band-structure validation: gap reproduces BC1967 Table I (2.58 eV) --
ok, gap, npw = W.validate_against_paper()
check(f"direct gap at Gamma = {gap:.3f} eV reproduces BC1967 2.58 eV", ok)

# --- the SBE cell: orthorhombic supercell + EXACT 2-fold folding --------------
# The SBE uses the orthorhombic al(1:3) cell (8 atoms). It must reproduce the
# primitive gap via folding, and the folding must be EXACT (block-diagonal over
# the 2 cosets, the analogue of the cubic 4-fold FCC folding).
off, gap_orth, g_cos0, g_cos1 = W.orth_folding_check(9.0)
check("orthorhombic supercell gap == primitive gap (folding consistent)",
      abs(gap_orth - gap) < 0.05)
check("2-fold folding is EXACT (off-coset |H| ~ 0)", off < 1e-10)
check("coset0 carries the direct gap (Gamma_hex)", abs(g_cos0 - gap) < 0.05)
check("coset1 (zone-edge partner) gap is larger (CBM is at Gamma)", g_cos1 > g_cos0 + 1.0)
check("2-fold coset table has 2 reps", W.WZ_COSET.shape == (2, 3))

if nfail == 0:
    print(f"PASS  (primitive gap@Gamma = {gap:.3f} eV vs BC1967 2.58 eV, {npw} PW; "
          f"orthorhombic SC folds exactly to it: gap={gap_orth:.3f} eV, off-block={off:.1e})")
    sys.exit(0)
else:
    print(f"FAIL ({nfail} checks)")
    sys.exit(1)
