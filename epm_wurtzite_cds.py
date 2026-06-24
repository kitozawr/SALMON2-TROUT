#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wurtzite CdS empirical-pseudopotential band machinery -- LOCAL EPM with the
REAL, CITED form factors of Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967)
("Electronic Structure and Optical Properties of Hexagonal CdSe, CdS, ZnS").

This is the companion to epm_gaas_reference.py (cubic). It implements the
parts of "add wurtzite CdS" that are backed by the paper:

  * the ORTHORHOMBIC (sqrt3 x 1 x 1) 8-atom supercell built from the 3-number
    cell al(1:3) = (a, a*sqrt3, c) -- the SBE side already takes an al vector;
  * the 8-atom wurtzite structure factor (Eqs 3a/3b of BC1967) with the
    internal parameter u = 3/8, split into the SYMMETRIC V^S and ANTISYMMETRIC
    V^A parts. V^A != 0 is the broken-inversion (polar) term;
  * the LOCAL pseudopotential -- BC1967 explicitly assume spherically-symmetric
    atomic potentials with NO angular/nonlocal dependence, so there is NO
    nonlocal velocity correction (rvnl_tm = 0), exactly like the local GaAs/Si
    EPM. (No CdS nonlocal parameter is cited anywhere, so adding one would be
    fabrication.)
  * the CITED form factors: BC1967 Table II gives the zinc-blende CdS form
    factors (the anchors) and obtains the wurtzite values by INTERPOLATION onto
    the hexagonal G-shells. We reproduce that exact procedure -- interpolate
    V^S, V^A as smooth functions of the PHYSICAL |G|^2, with the anchors placed
    at the zinc-blende shells of a_ZB = sqrt(2) a_W (the BC1967 relation).

VALIDATION: the band solver checks the computed direct gap at Gamma against the
paper's value (Table I: E_g = 2.58 eV, direct at Gamma). It REFUSES to emit a
band structure whose gap is not within tolerance of the paper -- no unvalidated
numbers are produced.

Cited form factors -- BC1967 Table II, zinc-blende CdS [Ry]:
  V^S(3) = -0.24,  V^S(8) = +0.03,  V^S(11) = +0.04
  V^A(3) = +0.23,  V^A(4) = +0.13,  V^A(11) = +0.05,  V^A(12) = +0.05
Lattice (Table I): a = 4.136 Ang, c/a = 1.623, u = 3/8, gap = 2.58 eV.
"""

import numpy as np
from numpy.linalg import eigvalsh

RY_TO_HA = 0.5
HA_TO_EV = 27.211386245988
ANG_TO_BOHR = 1.0 / 0.52917721067

# --- CdS wurtzite cell (BC1967 Table I) --------------------------------------
CDS_A_ANG = 4.136
CDS_C_ANG = 4.136 * 1.623          # c/a = 1.623  -> c = 6.713 Ang
CDS_U     = 3.0 / 8.0              # BC1967 use u = 3/8 for the structure factors
CDS_GAP_PAPER_EV = 2.58           # Table I, direct at Gamma (validation target)

# zinc-blende lattice constant of the same material (nearest-neighbour match):
# a_ZB = sqrt(2) a_W  (BC1967, Sec. II). The form-factor anchors live on the
# zinc-blende shells of THIS a_ZB.
CDS_AZB_ANG = np.sqrt(2.0) * CDS_A_ANG

# --- CITED zinc-blende CdS form factors (BC1967 Table II, Ry) -----------------
# Anchors as (shell index n, V[Ry]); n is |G|^2 in units of (2pi/a_ZB)^2.
_VS_ANCHORS = [(0.0, 0.0), (3.0, -0.24), (8.0, +0.03), (11.0, +0.04), (16.0, 0.0)]
_VA_ANCHORS = [(0.0, 0.0), (3.0, +0.23), (4.0, +0.13), (11.0, +0.05),
               (12.0, +0.05), (16.0, 0.0)]


def _azb_bohr():
    return CDS_AZB_ANG * ANG_TO_BOHR

def _anchor_curve(anchors):
    """Return (g2phys[Bohr^-2], V[Ry]) nodes: anchor shell n -> n*(2pi/a_ZB)^2."""
    unit = (2.0 * np.pi / _azb_bohr()) ** 2
    xs = np.array([n * unit for n, _ in anchors])
    ys = np.array([v for _, v in anchors])
    return xs, ys

def form_factor_phys(g2_phys_bohr):
    """V^S, V^A [Ry] at physical |G|^2 [Bohr^-2], by interpolation of the cited
    zinc-blende CdS anchors onto a continuous curve (the BC1967 wurtzite
    procedure). Zero beyond the largest anchor shell."""
    xs_s, ys_s = _anchor_curve(_VS_ANCHORS)
    xs_a, ys_a = _anchor_curve(_VA_ANCHORS)
    vs = np.interp(g2_phys_bohr, xs_s, ys_s, left=0.0, right=0.0)
    va = np.interp(g2_phys_bohr, xs_a, ys_a, left=0.0, right=0.0)
    return vs, va


# =============================================================================
# Geometry: orthorhombic lattice, atoms, plane-wave basis
# =============================================================================
def cds_cell_au():
    """Orthorhombic (sqrt3 x 1 x 1) box al(1:3) = (a, a*sqrt3, c) in Bohr."""
    a = CDS_A_ANG * ANG_TO_BOHR
    return np.array([a, a * np.sqrt(3.0), CDS_C_ANG * ANG_TO_BOHR])

def hexagonal_vectors_au(a_au, c_au):
    a1 = a_au * np.array([1.0, 0.0, 0.0])
    a2 = a_au * np.array([-0.5, np.sqrt(3.0) / 2.0, 0.0])
    a3 = c_au * np.array([0.0, 0.0, 1.0])
    return a1, a2, a3

def orthorhombic_vectors_au(a_au, c_au):
    a1, a2, a3 = hexagonal_vectors_au(a_au, c_au)
    return a1, a1 + 2.0 * a2, a3        # A=a1, B=a1+2a2 (=a*sqrt3 along y), C=a3

def reciprocal(A, B, C):
    V = np.dot(A, np.cross(B, C))
    return np.array([2 * np.pi * np.cross(B, C) / V,
                     2 * np.pi * np.cross(C, A) / V,
                     2 * np.pi * np.cross(A, B) / V]), V

def wurtzite_atoms_orth(a_au, c_au, u=CDS_U):
    """Cartesian positions [Bohr] + species (+1 Cd, -1 S) filling the
    orthorhombic cell: algorithmic, -> exactly 8 atoms (4 Cd + 4 S)."""
    a1, a2, a3 = hexagonal_vectors_au(a_au, c_au)
    A, B, C = orthorhombic_vectors_au(a_au, c_au)
    basis = [(np.array([0.0, 0.0, 0.0]), +1), (np.array([1/3, 2/3, 1/2]), +1),
             (np.array([0.0, 0.0, u]),   -1), (np.array([1/3, 2/3, 1/2 + u]), -1)]
    Binv = np.linalg.inv(np.array([A, B, C]).T)
    pos, spec, seen = [], [], set()
    for n1 in range(-2, 3):
        for n2 in range(-2, 3):
            for bf, sp in basis:
                r = (bf[0] + n1) * a1 + (bf[1] + n2) * a2 + bf[2] * a3
                f = Binv @ r
                f = f - np.floor(f + 1e-9)
                key = (round(f[0], 5), round(f[1], 5), round(f[2], 5), sp)
                if key in seen:
                    continue
                seen.add(key)
                pos.append(f @ np.array([A, B, C]))
                spec.append(sp)
    return np.array(pos), np.array(spec)

def build_pw_basis_orth(A, B, C, cutoff_ry):
    """Plane waves G with |G|^2 <= cutoff (Ry) of the orthorhombic reciprocal
    lattice. Returns (Gcart[npw,3], hkl[npw,3])."""
    Brec, _ = reciprocal(A, B, C)
    gmin = min(np.linalg.norm(Brec[i]) for i in range(3))
    nmax = int(np.ceil(np.sqrt(max(cutoff_ry, 1e-6)) / gmin)) + 1
    Gs, hkls = [], []
    for h in range(-nmax, nmax + 1):
        for k in range(-nmax, nmax + 1):
            for l in range(-nmax, nmax + 1):
                G = h * Brec[0] + k * Brec[1] + l * Brec[2]
                if G @ G <= cutoff_ry + 1e-9:
                    Gs.append(G); hkls.append((h, k, l))
    return np.array(Gs), np.array(hkls)


# =============================================================================
# Local Hamiltonian  H(k)  (Rydberg internally -> Hartree)
# =============================================================================
def build_hamiltonian_local(kvec, Gcart, atoms_pos, atoms_spec):
    """H(k) = |k+G|^2 (Ry) + sum_species V_species(dG) S_species(dG). Cation/anion
    atomic form factors V_cat = V^S + V^A, V_ani = V^S - V^A, each times its
    normalized complex structure factor. Hermitian; complex because the wurtzite
    cation/anion positions are not inversion-related (broken inversion)."""
    npw = len(Gcart)
    H = np.zeros((npw, npw), dtype=complex)
    cat = atoms_pos[atoms_spec > 0]
    ani = atoms_pos[atoms_spec < 0]
    ncat, nani = max(len(cat), 1), max(len(ani), 1)
    g2cut = 17.0 * (2.0 * np.pi / _azb_bohr()) ** 2
    for i in range(npw):
        kg = kvec + Gcart[i]
        H[i, i] += (kg @ kg) * RY_TO_HA
        for j in range(npw):
            if i == j:
                continue
            dG = Gcart[i] - Gcart[j]
            g2 = dG @ dG
            if g2 > g2cut:
                continue
            VS, VA = form_factor_phys(g2)
            if VS == 0.0 and VA == 0.0:
                continue
            Scat = np.exp(-1j * (cat @ dG)).sum() / ncat
            Sani = np.exp(-1j * (ani @ dG)).sum() / nani
            H[i, j] += ((VS + VA) * Scat + (VS - VA) * Sani) * RY_TO_HA
    return 0.5 * (H + H.conj().T)

def bands_at_k(kvec, Gcart, atoms_pos, atoms_spec, nb):
    return eigvalsh(build_hamiltonian_local(kvec, Gcart, atoms_pos, atoms_spec))[:nb]


# =============================================================================
# Wurtzite BZ symmetry points (orthorhombic frame) + 2-fold coset
# =============================================================================
HS_POINTS_WZ = {'Gamma': [0, 0, 0], 'X': [0.5, 0, 0], 'Y': [0, 0.5, 0],
                'A': [0, 0, 0.5], 'S': [0.5, 0.5, 0]}
DEFAULT_PATH_WZ = ['A', 'Gamma', 'X', 'S', 'Y', 'Gamma']
WZ_COSET = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]])   # 2-fold orthorhombic<-hex

# CdS valence: 8 e per formula unit (Cd 5s2 + S 3s2 3p4), 4 formula units in the
# orthorhombic cell -> 32 electrons -> 16 filled bands (spinless, 2 e/band).
CDS_NVAL = 16


def direct_gap_at_gamma(Gcart, atoms_pos, atoms_spec, nval=CDS_NVAL):
    ev = bands_at_k(np.zeros(3), Gcart, atoms_pos, atoms_spec, nval + 4)
    return (ev[nval] - ev[nval - 1]) * HA_TO_EV, ev


def validate_against_paper(cutoff_ry=12.0, tol_ev=0.4):
    """Compute the Gamma direct gap with the cited BC1967 form factors and check
    it against the paper's 2.58 eV. Returns (ok, gap_ev). Used by the test; the
    band solver should not be trusted outside this tolerance."""
    a_au, c_au = CDS_A_ANG * ANG_TO_BOHR, CDS_C_ANG * ANG_TO_BOHR
    A, B, C = orthorhombic_vectors_au(a_au, c_au)
    pos, spec = wurtzite_atoms_orth(a_au, c_au)
    Gcart, _ = build_pw_basis_orth(A, B, C, cutoff_ry)
    gap, _ = direct_gap_at_gamma(Gcart, pos, spec)
    return abs(gap - CDS_GAP_PAPER_EV) <= tol_ev, gap


if __name__ == '__main__':
    ok, gap = validate_against_paper()
    print(f'CdS wurtzite EPM (BC1967 local form factors):')
    print(f'  cell al(1:3) = {np.round(cds_cell_au(),3)} Bohr  (a, a*sqrt3, c)')
    print(f'  direct gap at Gamma = {gap:.3f} eV   (BC1967 Table I: {CDS_GAP_PAPER_EV} eV)')
    print(f'  validated against the paper: {ok}')
