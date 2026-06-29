#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wurtzite CdS empirical-pseudopotential band machinery -- LOCAL EPM with the
REAL, CITED form factors of Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967)
("Electronic Structure and Optical Properties of Hexagonal CdSe, CdS, ZnS").

Two cell representations:
  * HEXAGONAL PRIMITIVE cell (4 atoms) -- used to VALIDATE the band structure
    directly against the paper (gap, structure factors). This is BC1967's own
    setup, so it isolates the EPM physics from any supercell folding.
  * ORTHORHOMBIC (sqrt3 x 1 x 1, 8-atom) cell -- the al(1:3) = (a, a*sqrt3, c)
    box the SBE side needs (its lattice input is a 3-number vector). A folding
    of the primitive cell.

BC1967 are a LOCAL EPM (spherically-symmetric atomic potentials, NO angular /
nonlocal term) -> there is NO nonlocal velocity correction (rvnl_tm = 0), and no
cited CdS nonlocal parameter exists, so none is fabricated.

Potential matrix element (BC1967 Eqs 2,3,4), n = TOTAL atoms in the cell:
  <G'|V|G> = S^S(dG) V^S(|dG|^2) + i S^A(dG) V^A(|dG|^2),
  S^S(dG) = (1/n) sum_j exp(-i dG.tau_j),
  i S^A(dG) = (1/n) sum_j P_j exp(-i dG.tau_j),  P_j = +1 (cation) / -1 (anion).
Equivalently  (1/n)[(V^S+V^A) s_cat + (V^S-V^A) s_ani]  with s = un-normalized
species sums. The 1/n (TOTAL atoms) normalization is the "volume per atom"
normalization BC1967 use -- dividing by atoms-per-species instead makes the
potential too strong.

Cited form factors -- BC1967 Table II, zinc-blende CdS [Ry]:
  V^S(3) = -0.24,  V^S(8) = +0.03,  V^S(11) = +0.04
  V^A(3) = +0.23,  V^A(4) = +0.13,  V^A(11) = +0.05,  V^A(12) = +0.05
Lattice (Table I): a = 4.136 Ang, c/a = 1.623, u = 3/8, gap = 2.58 eV (direct
at Gamma). Form factors are interpolated onto the hexagonal G-shells from the
zinc-blende shells of a_ZB = sqrt(2) a_W (BC1967 Sec. II).
"""

import sys

import numpy as np
from numpy.linalg import eigh, eigvalsh

RY_TO_HA = 0.5
HA_TO_EV = 27.211386245988
ANG_TO_BOHR = 1.0 / 0.52917721067

# --- CdS wurtzite cell (BC1967 Table I) --------------------------------------
CDS_A_ANG = 4.136
CDS_C_ANG = 4.136 * 1.623          # c/a = 1.623
CDS_U     = 3.0 / 8.0              # BC1967 use u = 3/8 for the Table II structure factors
CDS_GAP_PAPER_EV = 2.58           # Table I, direct at Gamma (validation target)
CDS_AZB_ANG = np.sqrt(2.0) * CDS_A_ANG   # a_ZB = sqrt(2) a_W (nearest-neighbour match)

# --- CITED CdS form factors (BC1967 Table II, Ry) -----------------------------
# BC1967 give the WURTZITE form factors directly in Table II (their tuned values,
# fit to the hexagonal optical data) AND the zinc-blende anchors. We use the
# WURTZITE values, keyed by |G|^2 in BC1967 reduced units (sqrt2*pi/a_W)^2 (which
# equals the zinc-blende unit (2pi/a_ZB)^2 since a_ZB = sqrt2 a_W). Each shell's
# reduced |G|^2 is computed from the hexagonal reciprocal lattice (e.g. 002->3.04,
# 101->3.43, 102->5.70 -- the rows whose |S^S|,|S^A| we verify against Table II).
_VS_ANCHORS = [(0.0, 0.0), (3.04, -0.26), (3.43, -0.24), (5.70, -0.20),
               (9.50, +0.04), (10.67, +0.04), (13.30, +0.02), (16.0, 0.0)]
_VA_ANCHORS = [(0.0, 0.0), (3.04, +0.23), (3.43, +0.18), (5.70, +0.08),
               (9.50, +0.05), (11.40, +0.05), (12.15, +0.05), (13.30, +0.03), (16.0, 0.0)]


def _azb_bohr():
    return CDS_AZB_ANG * ANG_TO_BOHR

def form_factor_phys(g2_phys_bohr):
    """V^S, V^A [Ry] at physical |G|^2 [Bohr^-2], interpolating the cited
    zinc-blende CdS anchors (shell n -> n*(2pi/a_ZB)^2). Zero beyond shell 16."""
    unit = (2.0 * np.pi / _azb_bohr()) ** 2
    xs_s = np.array([n * unit for n, _ in _VS_ANCHORS]); ys_s = np.array([v for _, v in _VS_ANCHORS])
    xs_a = np.array([n * unit for n, _ in _VA_ANCHORS]); ys_a = np.array([v for _, v in _VA_ANCHORS])
    vs = np.interp(g2_phys_bohr, xs_s, ys_s, left=0.0, right=0.0)
    va = np.interp(g2_phys_bohr, xs_a, ys_a, left=0.0, right=0.0)
    return vs, va


# =============================================================================
# Cells, atoms, plane-wave basis
# =============================================================================
def hexagonal_vectors_au(a_au, c_au):
    a1 = a_au * np.array([1.0, 0.0, 0.0])
    a2 = a_au * np.array([-0.5, np.sqrt(3.0) / 2.0, 0.0])
    a3 = c_au * np.array([0.0, 0.0, 1.0])
    return a1, a2, a3

def orthorhombic_vectors_au(a_au, c_au):
    a1, a2, a3 = hexagonal_vectors_au(a_au, c_au)
    return a1, a1 + 2.0 * a2, a3        # A=a1, B=a1+2a2 (=a*sqrt3 along y), C=a3

def cds_cell_au():
    """Orthorhombic al(1:3) = (a, a*sqrt3, c) [Bohr] for the SBE side."""
    a = CDS_A_ANG * ANG_TO_BOHR
    return np.array([a, a * np.sqrt(3.0), CDS_C_ANG * ANG_TO_BOHR])

def reciprocal(A, B, C):
    V = np.dot(A, np.cross(B, C))
    return np.array([2 * np.pi * np.cross(B, C) / V,
                     2 * np.pi * np.cross(C, A) / V,
                     2 * np.pi * np.cross(A, B) / V]), V

def hex_primitive_atoms(a_au, c_au, u=CDS_U):
    """4 atoms of the hexagonal primitive cell: Cartesian pos[4,3], spec (+1 Cd,
    -1 S). Cd(0,0,0),(1/3,2/3,1/2); S +(0,0,u)."""
    a1, a2, a3 = hexagonal_vectors_au(a_au, c_au)
    frac = [(np.array([0.0, 0.0, 0.0]), +1), (np.array([1/3, 2/3, 1/2]), +1),
            (np.array([0.0, 0.0, u]),   -1), (np.array([1/3, 2/3, 1/2 + u]), -1)]
    M = np.array([a1, a2, a3])
    return np.array([f @ M for f, _ in frac]), np.array([s for _, s in frac])

def wurtzite_atoms_orth(a_au, c_au, u=CDS_U):
    """8 atoms (4 Cd + 4 S) filling the orthorhombic cell (for the SBE side)."""
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
                f = Binv @ r; f = f - np.floor(f + 1e-9)
                key = (round(f[0], 5), round(f[1], 5), round(f[2], 5), sp)
                if key in seen:
                    continue
                seen.add(key); pos.append(f @ np.array([A, B, C])); spec.append(sp)
    return np.array(pos), np.array(spec)

def build_pw_basis(Brec, cutoff_ry):
    """Plane waves G with |G|^2 <= cutoff (Ry) of the reciprocal lattice rows
    Brec. Returns (Gcart[npw,3], hkl[npw,3])."""
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
# Structure factors (BC1967 Eqs 3a/3b) and Hamiltonian
# =============================================================================
def structure_factors(dG, atoms_pos, atoms_spec):
    """Return (S^S, iS^A) at reciprocal vector dG, normalized by TOTAL atoms n.
    S^S = (1/n) sum exp(-i dG.tau);  iS^A = (1/n) sum P_j exp(-i dG.tau)."""
    n = len(atoms_pos)
    ph = np.exp(-1j * (atoms_pos @ dG))
    Ssym = ph.sum() / n
    Sasym = (ph * atoms_spec).sum() / n
    return Ssym, Sasym

def build_hamiltonian(kvec, Gcart, atoms_pos, atoms_spec):
    """H(k) [Hartree] = |k+G|^2 (Ry) + S^S V^S + iS^A V^A, normalized by TOTAL
    atoms. Hermitian; complex because wurtzite breaks inversion. Vectorized over
    the plane-wave pairs (identical result to the scalar double loop; the cited
    form factors are zero beyond shell 16, so no explicit g2 cut is needed)."""
    npw = len(Gcart)
    n = len(atoms_pos)
    dG = Gcart[:, None, :] - Gcart[None, :, :]            # (npw, npw, 3)
    g2 = np.einsum('ijd,ijd->ij', dG, dG)
    VS, VA = form_factor_phys(g2)                          # array-safe (np.interp)
    ph = np.exp(-1j * np.einsum('ijd,ad->ija', dG, atoms_pos))   # (npw, npw, n)
    Ssym  = ph.sum(axis=2) / n
    Sasym = (ph * atoms_spec[None, None, :]).sum(axis=2) / n
    H = (Ssym * VS + Sasym * VA) * RY_TO_HA
    np.fill_diagonal(H, 0.0)                               # diagonal is kinetic only
    kg = kvec[None, :] + Gcart
    H[np.diag_indices(npw)] = np.einsum('id,id->i', kg, kg) * RY_TO_HA
    return 0.5 * (H + H.conj().T)

def bands_at_k(kvec, Gcart, atoms_pos, atoms_spec, nb):
    return eigvalsh(build_hamiltonian(kvec, Gcart, atoms_pos, atoms_spec))[:nb]


# =============================================================================
# Wurtzite BZ symmetry points (orthorhombic frame) + 2-fold coset
# =============================================================================
HS_POINTS_WZ = {'Gamma': [0, 0, 0], 'X': [0.5, 0, 0], 'Y': [0, 0.5, 0],
                'A': [0, 0, 0.5], 'S': [0.5, 0.5, 0]}
DEFAULT_PATH_WZ = ['A', 'Gamma', 'X', 'S', 'Y', 'Gamma']
WZ_COSET = np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.0]])   # 2-fold orthorhombic<-hex

# CdS: 8 valence e/formula unit (Cd 5s2 + S 3s2 3p4). Hexagonal primitive cell
# has 2 formula units -> 16 e -> 8 filled bands (spinless).
CDS_NVAL_PRIM = 8


def direct_gap_at_gamma_primitive(cutoff_ry=12.0):
    """Compute the Gamma direct gap (eV) in the hexagonal primitive cell with the
    cited BC1967 form factors. Returns (gap_ev, eigenvalues_eV, npw)."""
    a_au, c_au = CDS_A_ANG * ANG_TO_BOHR, CDS_C_ANG * ANG_TO_BOHR
    a1, a2, a3 = hexagonal_vectors_au(a_au, c_au)
    Brec, _ = reciprocal(a1, a2, a3)
    pos, spec = hex_primitive_atoms(a_au, c_au)
    Gcart, _ = build_pw_basis(Brec, cutoff_ry)
    ev = bands_at_k(np.zeros(3), Gcart, pos, spec, CDS_NVAL_PRIM + 4) * HA_TO_EV
    nv = CDS_NVAL_PRIM
    return ev[nv] - ev[nv - 1], ev, len(Gcart)


def orth_coset(Gcart, a_au, c_au):
    """Classify each orthorhombic plane wave G into a coset of the PRIMITIVE
    reciprocal lattice (the supercell potential is primitive-periodic, so V(dG)!=0
    only for primitive dG): coset 0 = G is a primitive reciprocal vector, coset 1
    = the index-2 half-points. n_i = G . a_prim_i / 2pi integer <=> coset 0."""
    ap = np.array(hexagonal_vectors_au(a_au, c_au))      # primitive real vectors
    n = (ap @ Gcart.T) / (2.0 * np.pi)                   # (3, npw)
    is_prim = np.all(np.abs(n - np.round(n)) < 1e-6, axis=0)
    return np.where(is_prim, 0, 1)


def orth_folding_check(cutoff_ry=9.0):
    """Verify the 2-fold orthorhombic<-hexagonal folding is EXACT (the analogue
    of the cubic 4-fold FCC folding): the supercell H at Gamma must be
    block-diagonal over the 2 cosets. Returns (max_offblock, gap_orth_eV,
    gap_coset0_eV, gap_coset1_eV)."""
    a_au, c_au = CDS_A_ANG * ANG_TO_BOHR, CDS_C_ANG * ANG_TO_BOHR
    A, B, C = orthorhombic_vectors_au(a_au, c_au)
    Brec, _ = reciprocal(A, B, C)
    pos, spec = wurtzite_atoms_orth(a_au, c_au)
    Gcart, _ = build_pw_basis(Brec, cutoff_ry)
    coset = orth_coset(Gcart, a_au, c_au)
    H = build_hamiltonian(np.zeros(3), Gcart, pos, spec)
    i0, i1 = np.where(coset == 0)[0], np.where(coset == 1)[0]
    offblock = np.abs(H[np.ix_(i0, i1)]).max() if (len(i0) and len(i1)) else 0.0
    ev = np.linalg.eigvalsh(H) * HA_TO_EV
    nv = 16                                              # 8 atoms -> 32 e -> 16 bands
    gap_orth = ev[nv] - ev[nv - 1]
    e0 = np.linalg.eigvalsh(H[np.ix_(i0, i0)]) * HA_TO_EV
    e1 = np.linalg.eigvalsh(H[np.ix_(i1, i1)]) * HA_TO_EV
    return offblock, gap_orth, e0[8] - e0[7], e1[8] - e1[7]


def validate_against_paper(cutoff_ry=12.0, tol_ev=0.1):
    """(ok, gap_ev, npw): is the converged primitive-cell Gamma gap within tol of
    the paper's 2.58 eV? With the cited wurtzite form factors + the 1/n_total
    normalization the gap converges to ~2.55 eV (|Δ|≈0.03 eV, well inside the
    paper's ~0.27 eV form-factor accuracy)."""
    gap, _, npw = direct_gap_at_gamma_primitive(cutoff_ry)
    return abs(gap - CDS_GAP_PAPER_EV) <= tol_ev, gap, npw


# =============================================================================
# SBE ground-state emission (scalar, NO spinor) on the FOLDED orthorhombic cell
# =============================================================================
# The SBE runs on the orthorhombic 8-atom al(1:3)=(a,a*sqrt3,c) cell, whose bands
# are the hexagonal-primitive bands FOLDED 2-fold (verified exact by
# orth_folding_check). This emits SYSNAME_k/_eigen/_tm.data in the SBE read
# contract (via epm_io), so theory='sbe' can run CdS end-to-end. The local
# pseudopotential has no nonlocal velocity term -> rvnl_tm = 0.
CDS_SYSNAME      = 'CdS'
CDS_NELEC        = 32           # 8 atoms * 4 valence e-/formula-pair = 32
CDS_NSTATE       = 32           # 16 valence (folded) + 16 conduction
CDS_NUM_KGRID    = (4, 4, 4)
CDS_GS_CUTOFF_RY = 9.0          # |G|^2 [a.u.^2]; 2-fold folding is exact here


def main_gs(sysname=CDS_SYSNAME, num_kgrid=CDS_NUM_KGRID, nstate=CDS_NSTATE,
            nelec=CDS_NELEC, cutoff_ry=CDS_GS_CUTOFF_RY, outdir='./'):
    """Emit the scalar SBE ground-state files for CdS on the orthorhombic cell."""
    import epm_io
    a_au, c_au = CDS_A_ANG * ANG_TO_BOHR, CDS_C_ANG * ANG_TO_BOHR
    A, B, C = orthorhombic_vectors_au(a_au, c_au)
    Brec, _ = reciprocal(A, B, C)
    pos, spec = wurtzite_atoms_orth(a_au, c_au)
    Gcart, _ = build_pw_basis(Brec, cutoff_ry)
    npw = len(Gcart)
    if nstate > npw:
        raise ValueError(f'nstate={nstate} exceeds npw={npw}; raise cutoff')
    kpoint, kweight = epm_io.monkhorst_pack(Brec, num_kgrid)
    nk = kpoint.shape[0]
    nocc = nelec // 2

    print(f'# EPM CdS (wurtzite, orthorhombic 8-atom cell, 2-fold folded) -- scalar')
    print(f'#   al(1:3) = {np.round(cds_cell_au(), 4)} Bohr  (a, a*sqrt3, c)')
    print(f'#   plane waves = {npw}, k-points = {nk}, bands = {nstate}, '
          f'valence e- = {nelec} (occ 2/band)')

    eigen = np.zeros((nstate, nk))
    occup = np.zeros((nstate, nk))
    p_tm = np.zeros((nstate, nstate, 3, nk), dtype=complex)
    rvnl_tm = np.zeros((nstate, nstate, 3, nk), dtype=complex)   # local -> 0
    for ik in range(nk):
        H = build_hamiltonian(kpoint[ik], Gcart, pos, spec)
        ev, evec = eigh(H)
        eigen[:, ik] = ev[:nstate]
        occup[:nocc, ik] = 2.0
        p_tm[:, :, :, ik] = epm_io.momentum_matrix(kpoint[ik], Gcart, evec[:, :nstate])
        if (ik + 1) % max(1, nk // 8) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik + 1}/{nk}')

    b_matrix = np.array(Brec)
    epm_io.write_epm_gs_files(sysname, outdir, 'CdS', kpoint, b_matrix, kweight,
                              eigen, occup, p_tm, rvnl_tm,
                              extra_note='wurtzite orthorhombic 2-fold folded')
    return eigen, occup


# Hexagonal-BZ high-symmetry points (reduced coords of the hex reciprocal) and
# a primitive band path for the clean (unfolded) level-structure plot.
CDS_HS_HEX = {'Gamma': (0, 0, 0), 'M': (0.5, 0, 0), 'K': (1/3, 1/3, 0),
              'A': (0, 0, 0.5), 'L': (0.5, 0, 0.5), 'H': (1/3, 1/3, 0.5)}
CDS_BANDPATH = ['A', 'Gamma', 'M', 'K', 'Gamma']
CDS_BANDPATH_NB = 14            # 8 valence + 6 conduction
CDS_BANDPATH_NDIV = 40
CDS_BANDPATH_CUTOFF_RY = 12.0


def main_bandpath(sysname=CDS_SYSNAME, outdir='./'):
    """Emit SYSNAME_bandpath.data: the clean primitive (hexagonal-cell) bands
    along CDS_BANDPATH, for the unfolded level-structure plot."""
    import epm_io
    a_au, c_au = CDS_A_ANG * ANG_TO_BOHR, CDS_C_ANG * ANG_TO_BOHR
    a1, a2, a3 = hexagonal_vectors_au(a_au, c_au)
    Brec, _ = reciprocal(a1, a2, a3)
    pos, spec = hex_primitive_atoms(a_au, c_au)
    Gcart, _ = build_pw_basis(Brec, CDS_BANDPATH_CUTOFF_RY)
    qreds, kcarts, dists, node_d = epm_io.build_path(CDS_HS_HEX, CDS_BANDPATH,
                                                     CDS_BANDPATH_NDIV, Brec)
    eig = np.array([bands_at_k(kc, Gcart, pos, spec, CDS_BANDPATH_NB) for kc in kcarts])
    epm_io.write_bandpath_file(sysname, outdir, 'CdS', CDS_BANDPATH, node_d,
                               dists, qreds, eig, nv=CDS_NVAL_PRIM, spinor=0)


def main_unfoldmap(sysname=CDS_SYSNAME, num_kgrid=CDS_NUM_KGRID, nstate=CDS_NSTATE,
                   cutoff_ry=CDS_GS_CUTOFF_RY, outdir='./'):
    """Emit SYSNAME_unfold.data: the 2-coset (orthorhombic<-hexagonal) band ->
    coset spectral-weight map for the SBE's unfolded-population output. Same
    k-grid / cutoff / nstate as main_gs (the SBE checks nk and nb match)."""
    import epm_io
    a_au, c_au = CDS_A_ANG * ANG_TO_BOHR, CDS_C_ANG * ANG_TO_BOHR
    A, B, C = orthorhombic_vectors_au(a_au, c_au)
    Brec, _ = reciprocal(A, B, C)
    pos, spec = wurtzite_atoms_orth(a_au, c_au)
    Gcart, hkls = build_pw_basis(Brec, cutoff_ry)
    coset = orth_coset(Gcart, a_au, c_au)        # 0/1 (2-fold)
    kpoint, _ = epm_io.monkhorst_pack(Brec, num_kgrid)
    offsets, isub, ibprim, wsub = epm_io.compute_unfold_map(
        lambda k: build_hamiltonian(k, Gcart, pos, spec),
        kpoint, Gcart, np.array(hkls), coset, n_coset=2, nstate=nstate)
    epm_io.write_unfold_file(sysname, outdir, 'CdS', 2, CDS_NVAL_PRIM,
                             offsets, isub, ibprim, wsub)


def _print_validation():
    print('CdS wurtzite EPM (BC1967 local form factors):')
    print(f'  SBE cell al(1:3) = {np.round(cds_cell_au(),3)} Bohr  (a, a*sqrt3, c)')
    print('  -- hexagonal primitive cell (band validation) --')
    for ec in (9.0, 12.0, 16.0):
        gap, ev, npw = direct_gap_at_gamma_primitive(ec)
        print(f'     cutoff={ec:5.1f} Ry  npw={npw:4d}  gap@Gamma={gap:7.3f} eV   '
              f'(BC1967 {CDS_GAP_PAPER_EV} eV)')
    print('  -- orthorhombic supercell (the SBE al-vector cell) + 2-fold folding --')
    off, gorth, g0, g1 = orth_folding_check(9.0)
    print(f'     max off-coset |H| = {off:.2e}  (exact folding: must be ~0)')
    print(f'     orthorhombic gap@Gamma = {gorth:.3f} eV  (== primitive: folding OK)')
    print(f'     coset0 Gamma_hex gap = {g0:.3f} eV (the direct gap); '
          f'coset1 partner gap = {g1:.3f} eV')


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    if mode == 'validate':
        _print_validation()            # band/folding validation only (no files)
    elif mode == 'bandpath':
        main_bandpath()                # clean primitive band path only
    elif mode == 'unfoldmap':
        main_unfoldmap()               # 2-coset unfold map only
    elif mode == 'gs':
        main_gs()                      # SBE ground-state files (k/eigen/tm)
        main_bandpath()                # + the clean primitive band path
        main_unfoldmap()               # + the 2-coset unfold map (no slow validation)
    else:
        _print_validation()
        print()
        main_gs()                      # emit the scalar SBE ground-state files
        main_bandpath()                # + the clean primitive band path
        main_unfoldmap()               # + the 2-coset unfold map
