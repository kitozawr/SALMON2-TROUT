#!/usr/bin/env python3
"""
test_graphene_epm.py  -  monolayer graphene local EPM (Ramanujam thesis).

Validates the Config-A (strict-2D, three in-plane form factors) local EPM
against the thesis acceptance tests (G1.5 of the task / Ramanujam 2015
Fig 4.3/4.8/4.9):

  * geometry: a = 2.46 Ang, C-C bond = 1.42 Ang = a/sqrt3, 2 carbon atoms;
  * H(k) real-symmetric (centrosymmetric D6h, V_A = 0);
  * pi/pi* touch with ZERO gap at the Dirac point K (symmetry-protected);
  * linear dispersion, v_F in [0.8, 1.0]e6 m/s (thesis 1L 7.65e5);
  * Gamma valence-band bottom in [-8.3, -7.8] eV below the Dirac point;
  * M-point dip in [-3.0, -2.5] eV;
  * the UNIT-TRAP fix: the form-factor subscript is |G|^2 in (2pi/(sqrt3 a))^2
    units, so V4 sits on the first reciprocal shell and couples the Dirac
    triplet (the cone collapses to the free-electron value otherwise).

Run:  python3 tests/test_graphene_epm.py   (pure numpy, no SALMON build).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import epm_graphene as G

nfail = 0
def check(name, cond):
    global nfail
    if not cond:
        print(f"  FAIL: {name}")
        nfail += 1

# --- geometry ----------------------------------------------------------------
check("lattice constant a = 2.46 Ang", abs(G.A_LATT - 2.46) < 1e-9)
check("C-C bond = a/sqrt3 = 1.42 Ang", abs(G.assert_geometry() - 2.46 / np.sqrt(3)) < 1e-6)
a1, a2 = G.lattice_vectors()
check("two lattice vectors of equal length a", abs(np.linalg.norm(a1) - np.linalg.norm(a2)) < 1e-9)
tau = G.basis_atoms()
check("two carbon atoms in the basis", len(tau) == 2)
# nearest-neighbour distance between the two sublattices = bond
check("sublattice separation = bond", abs(np.linalg.norm(tau[1] - tau[0]) - G.D_BOND) < 1e-6)

# --- Hamiltonian is real-symmetric (V_A = 0, centrosymmetric) ----------------
Gc, _ = G.build_pw_basis(400.0)
hs = G.high_symmetry_points()
H = G.build_hamiltonian(hs['M'], Gc, tau)
check("H Hermitian", np.abs(H - H.conj().T).max() < 1e-10)
check("H real (V_A = 0, centrosymmetric)", np.abs(H.imag).max() < 1e-9)

# --- |K| is the canonical Dirac wavevector 4pi/3a ----------------------------
check("|K| = 4pi/3a", abs(np.linalg.norm(hs['K']) - 4 * np.pi / (3 * G.A_LATT)) < 1e-9)

# --- the unit-trap fix: V4/V12/V16 land on the 1st/3rd/4th reciprocal shells -
# |G|^2 = |b1|^2 * m (m = h^2+hk+k^2); in the (2pi/(sqrt3 a))^2 unit m=1,3,4
# map to the subscripts n=4,12,16. (In the WRONG (2pi/a)^2 unit V4 would land
# off-shell at n=4/3 and never couple the Dirac triplet -> no cone.)
b1, _ = G.reciprocal_vectors()
g1 = b1 @ b1
check("V4 active at the 1st reciprocal shell (m=1)", G.form_factor(g1) == -8.23)
check("V12 active at the 3rd reciprocal shell (m=3)", G.form_factor(3 * g1) == 1.5)
check("V16 active at the 4th reciprocal shell (m=4)", G.form_factor(4 * g1) == 0.05)
check("no form factor on the 2nd shell (m=... none) ", G.form_factor(2 * g1) == 0.0)

# --- the thesis acceptance tests (G1.5) --------------------------------------
ok, d = G.validate_against_thesis(400.0)
check(f"zero gap at Dirac point (gap={d['gap']:.5f} eV)", abs(d['gap']) < 1e-3)
check(f"v_F = {d['vF']:.3e} m/s in [0.8,1.0]e6", 0.8e6 <= d['vF'] <= 1.0e6)
check(f"Gamma VB-bottom = {d['g_bottom']:.2f} eV in [-8.4,-7.7]",
      -8.4 <= d['g_bottom'] <= -7.7)
check(f"M-point dip = {d['m_dip']:.2f} eV in [-3.1,-2.4]",
      -3.1 <= d['m_dip'] <= -2.4)
check("all thesis acceptance tests pass", ok)

if nfail == 0:
    print(f"PASS  (Dirac gap={d['gap']:.5f} eV, v_F={d['vF']:.3e} m/s, "
          f"Gamma={d['g_bottom']:.2f} eV, M={d['m_dip']:.2f} eV, {d['npw']} PW)")
    sys.exit(0)
else:
    print(f"FAIL ({nfail} checks)")
    sys.exit(1)
