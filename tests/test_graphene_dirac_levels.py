#!/usr/bin/env python3
"""
test_graphene_dirac_levels.py - exercise x14 LEVEL CHECK for the graphene
Dirac cone the SBE runs on (Ramanujam pi-EPM; epm_graphene.py is the Python
twin of the in-SALMON Fortran EPM: both select plane waves by |G|^2 <= cutoff,
k-independent), plus the k-resolution advisory for resonant interband physics
on a uniform Monkhorst-Pack mesh.

THE TRAP THIS TEST GUARDS (found 2026-09-04): the 7-plane-wave basis of the
original x11 graphene input (epm_pw_cutoff_ry = 2.94, i.e. |G|^2 <= 2.94 a.u.,
G = 0 + the first shell) is NOT closed under the little group C3v of the K
point (the rotation about K maps a first-shell G onto a SECOND-shell vector
b2-b1 that the basis lacks), so the symmetry protection of the Dirac
degeneracy is broken by the truncation: the "gapless" cone acquires a
SPURIOUS 0.21 eV gap at K. A shell-complete basis up to |G|^2 <= 29.4 a.u.
(43 plane waves, the Python 400 eV cutoff) restores the degeneracy to ~1e-5 eV
and meets the thesis acceptance numbers (v_F, Gamma bottom, M dip). x14 uses
the 43-PW basis; this test fails if anyone reverts to 7 PW believing it gapless.

Checks
  1) 43-PW basis: |gap(K)| < 1 meV; v_F in [0.8, 1.0]e6 m/s (thesis window);
     e-h symmetry near K within 5 %; linear cone E(2q)/E(q) = 2 within 5 %
     up to ~0.5 eV.
  2) 7-PW basis: gap(K) > 0.1 eV (the documented truncation artifact).
  3) ratio of the EPM v_F to the 1e8 cm/s the Rana constants assume (the
     balance density n_i ~ 1/v_F^2 -- a 4 % velocity shift is 8 % in n_i).
  4) resolution advisory helper: points per resonance-shell radius on an
     N x N mesh = (hw/2v_F)/(|b|/N); N = 12 gives < 1 at 0.8 eV (the mesh sees
     only the K point), N = 150 gives >= 3. Table printed for the README.
Run:  python3 tests/test_graphene_dirac_levels.py   (pure numpy)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import epm_graphene as G  # noqa: E402

HBAR2_2M = G.HBAR2_2M          # eV Ang^2
BOHR_ANG = 0.52917721067
HBAR_EVS = G.HBAR_EVS
ANG = G.ANG
VF_RANA_MS = 1.0e6             # the R07 constant used by the Rana/Auger channel

nfail = 0


def check(name, cond):
    global nfail
    if not cond:
        print(f"  FAIL: {name}")
        nfail += 1


def cutoff_ev_from_au(g2_au):
    """Fortran epm_pw_cutoff_ry (|G|^2 <= X a.u.) -> Python kinetic cutoff [eV]."""
    return HBAR2_2M * g2_au / BOHR_ANG**2


def dirac_pair_at(kvec, Gc, tau):
    """Energies of the two lowest bands (the pi/pi* pair the nstate=2 SBE uses)."""
    ev = G.bands_at_k(kvec, Gc, tau)
    return ev[0], ev[1]


hs = G.high_symmetry_points()
K, Gam, M = hs['K'], hs['Gamma'], hs['M']
tau = G.basis_atoms()

# --- the two bases -----------------------------------------------------------
cut7 = cutoff_ev_from_au(2.94)      # the x11 input:  7 PW
cut43 = cutoff_ev_from_au(29.4)     # x14 (shell-complete to n=12): 43 PW
G7, _ = G.build_pw_basis(cut7)
G43, _ = G.build_pw_basis(cut43)
print(f"  basis sizes: |G|^2<=2.94 a.u. -> {len(G7)} PW ;  |G|^2<=29.4 a.u. -> {len(G43)} PW")
check("7-PW basis reproduces the x11 count", len(G7) == 7)
check("29.4 a.u. basis is the 43-PW shell-complete set", len(G43) == 43)

# --- (2) the trap: 7 PW is gapped at K ----------------------------------------
e0, e1 = dirac_pair_at(K, G7, tau)
gap7 = e1 - e0
print(f"  7-PW  : E(K) pair = {e0:8.4f} / {e1:8.4f} eV  ->  gap(K) = {gap7:.4f} eV  (spurious: basis breaks C3v at K)")
check("7-PW basis has the documented spurious gap (> 0.1 eV)", gap7 > 0.1)

# --- (1) 43 PW: gapless, v_F, e-h symmetry, linearity -------------------------
e0, e1 = dirac_pair_at(K, G43, tau)
gap43 = e1 - e0
Ed = 0.5 * (e0 + e1)
print(f"  43-PW : E(K) pair = {e0:8.4f} / {e1:8.4f} eV  ->  gap(K) = {gap43:.2e} eV")
check("43-PW basis: Dirac point gapless (< 1 meV)", abs(gap43) < 1e-3)

vfs, asym = [], []
for tgt in (Gam, M):
    kdir = (tgt - K) / np.linalg.norm(tgt - K)
    for q in (0.02, 0.04):                      # Ang^-1  (E ~ 0.13 - 0.5 eV)
        ev0, ev1 = dirac_pair_at(K + q * kdir, G43, tau)
        Ec, Ev = ev1 - Ed, Ed - ev0
        asym.append(abs(Ec - Ev) / max(Ec, 1e-12))
        vfs.append(0.5 * (Ec + Ev) / q * ANG / HBAR_EVS)
vF = float(np.mean(vfs))
print(f"  43-PW : v_F = {vF:.3e} m/s (isotropic mean, q = 0.02-0.04 A^-1 toward Gamma and M);"
      f" max e-h asymmetry = {max(asym):.3f}")
check("43-PW basis: v_F in the thesis window [0.8, 1.0]e6 m/s", 0.8e6 <= vF <= 1.0e6)
check("43-PW basis: electron-hole symmetric near K (< 5 %)", max(asym) < 0.05)
kdir = (M - K) / np.linalg.norm(M - K)
E1 = np.mean(dirac_pair_at(K + 0.04 * kdir, G43, tau) * np.array([-1, 1])) - 0  # (e1 - e0)/2 style below
Eq = lambda q: 0.5 * (lambda p: p[1] - p[0])(dirac_pair_at(K + q * kdir, G43, tau))   # noqa: E731
lin = Eq(0.08) / Eq(0.04)
print(f"  43-PW : linearity E(2q)/E(q) = {lin:.3f} at E ~ {Eq(0.08):.2f} eV (toward M)")
check("43-PW basis: linear cone up to ~0.5 eV (E(2q)/E(q) = 2 within 5 %)", abs(lin - 2.0) < 0.10)
gb = G.bands_at_k(Gam, G43, tau)[0] - Ed
md = G.bands_at_k(M, G43, tau)[0] - Ed
print(f"  43-PW : Gamma VB bottom = {gb:.2f} eV, M dip = {md:.2f} eV (thesis: -8.3..-7.8, -3.0..-2.5)")
check("43-PW basis: Gamma bottom in the thesis window", -8.4 <= gb <= -7.7)
check("43-PW basis: M dip in the thesis window", -3.1 <= md <= -2.4)

# --- (3) v_F vs the Rana constant -------------------------------------------------
r = vF / VF_RANA_MS
print(f"  v_F(EPM)/v_F(Rana 1e8 cm/s) = {r:.3f}  ->  balance density n_i ~ 1/v_F^2 differs by x{1 / r**2:.3f}")
check("EPM v_F within 10 % of the Rana constant (rates consistent to ~20 %)", abs(r - 1) < 0.10)

# --- (4) resolution advisory --------------------------------------------------------
b_len_au = np.linalg.norm(G.reciprocal_vectors()[0]) * BOHR_ANG       # |b1| in a.u.
vF_au = vF / 2.18769126364e6


def shell_points(hw_ev, N):
    k_res = (hw_ev / 27.211386245988) / (2.0 * vF_au)
    return k_res / (b_len_au / N)


print("  resonance-shell resolution (mesh points per shell radius k_res = hw/2v_F):")
print("     N   :  " + "  ".join(f"{N:5d}" for N in (12, 24, 48, 96, 150, 300)))
for hw in (0.1, 0.4, 0.8, 1.5):
    print(f"   {hw:4.1f} eV:  " + "  ".join(f"{shell_points(hw, N):5.2f}" for N in (12, 24, 48, 96, 150, 300)))
check("12x12 cannot resolve the 0.8 eV shell (< 1 point/radius)", shell_points(0.8, 12) < 1.0)
check("150x150 resolves the 0.8 eV shell (>= 3 points/radius)", shell_points(0.8, 150) >= 3.0)
check("advisory scales linearly with N", abs(shell_points(0.8, 300) / shell_points(0.8, 150) - 2.0) < 1e-12)

if nfail == 0:
    print("PASS  (graphene Dirac levels: 43-PW cone gapless/v_F/symmetric/linear; 7-PW trap documented; resolution advisory)")
    sys.exit(0)
print(f"FAIL ({nfail} checks)")
sys.exit(1)
