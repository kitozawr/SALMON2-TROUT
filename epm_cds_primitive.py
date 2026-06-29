"""
CdS wurtzite PRIMITIVE-CELL EPM ground state -- 4-atom HEXAGONAL primitive,
NON-orthogonal, NO folding. The CdS companion to epm_wurtzite_cds.py (which runs
the SBE on a 2-fold-folded orthorhombic sqrt3 x 1 x 1 supercell); this module
runs the SBE on the genuine hexagonal primitive cell instead -- the same
"unfolded test" we did for GaAs/Si, but on a truly hexagonal (a1,a2 at 120 deg)
non-orthogonal lattice, not the FCC rhombohedron.

Geometry (Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967); reused verbatim
from epm_wurtzite_cds):
  a1 = a (1,0,0),  a2 = a (-1/2, sqrt3/2, 0),  a3 = c (0,0,1)
  4 atoms: Cd(0,0,0),(1/3,2/3,1/2); S +(0,0,u), u=3/8.   a=4.136 A, c/a=1.623.
The hexagonal reciprocal lattice is itself hexagonal (b1,b2 at 120 deg). The MP
mesh on it is a non-orthogonal (triclinic) k-grid; the clean SBE reads Cartesian
p / band energies and is grid-agnostic, the Coulomb kernel is metric-aware, and
the reduced k it stores are only labels. There is NO 2-coset folding, NO unfold
map, NO sublattice projection -- the primitive cell is already the irreducible
problem. Local pseudopotential (BC1967) -> no nonlocal velocity term (rvnl=0).

CdS has 8 valence e- per formula unit; the primitive cell holds 2 formula units
-> 16 valence e- -> 8 filled bands (scalar/spinless). nstate = 16 (8 v + 8 c).

Output: CdS_prim_k.data (reduced hex-reciprocal coords + the b1/b2/b3 reciprocal
vectors in the header), _eigen.data, _tm.data, _bandpath.data. No _unfold.data.

Usage:
    python3 epm_cds_primitive.py            # validate gap + write GS + bandpath
    python3 epm_cds_primitive.py gap        # just the Gamma direct gap (vs 2.58 eV)
    python3 epm_cds_primitive.py gs         # GS files + bandpath
    python3 epm_cds_primitive.py bandpath   # bandpath only
"""
import sys

import numpy as np

import epm_wurtzite_cds as wz       # hexagonal-primitive geometry, H, form factors
import epm_io
import epm_gaas_primitive as fcc    # reuse the proven non-orthogonal GS writer

SYSNAME = 'CdS_prim'
MATERIAL = 'CdS'
NUM_KGRID = (7, 7, 5)              # ODD in-plane -> Gamma is an explicit grid point
NELEC = 16                         # 2 formula units x 8 valence e- -> 8 filled bands
NSTATE = 16                        # 8 valence + 8 conduction
CUTOFF_RY = 12.0                   # primitive Gamma gap converged to ~2.55 eV here


def _cell_au():
    a_au = wz.CDS_A_ANG * wz.ANG_TO_BOHR
    c_au = wz.CDS_C_ANG * wz.ANG_TO_BOHR
    return a_au, c_au


def report_gap():
    """Print the converged Gamma direct gap of the hexagonal primitive cell."""
    for ec in (9.0, 12.0, 16.0):
        gap, ev, npw = wz.direct_gap_at_gamma_primitive(ec)
        print(f'  cutoff={ec:5.1f} Ry  npw={npw:4d}  gap@Gamma = {gap:7.3f} eV   '
              f'(BC1967 {wz.CDS_GAP_PAPER_EV} eV)')


def main_gs(outdir='./'):
    a_au, c_au = _cell_au()
    a1, a2, a3 = wz.hexagonal_vectors_au(a_au, c_au)
    Brec, V = wz.reciprocal(a1, a2, a3)                  # rows b1,b2,b3 [a.u.]
    pos, spec = wz.hex_primitive_atoms(a_au, c_au)
    Gcart, _ = wz.build_pw_basis(Brec, CUTOFF_RY)
    kpoint, kweight = epm_io.monkhorst_pack(Brec, NUM_KGRID)
    kfrac = kpoint @ np.linalg.inv(Brec)
    nk, npw = len(kpoint), len(Gcart)
    if NSTATE > npw:
        raise ValueError(f'NSTATE={NSTATE} exceeds npw={npw}; raise CUTOFF_RY')
    nocc = NELEC // 2
    print(f'# EPM CdS PRIMITIVE wurtzite (4-atom hexagonal, non-orthogonal, NO folding)')
    print(f'#   a = {a_au:.4f} Bohr, c = {c_au:.4f} Bohr (c/a = {c_au/a_au:.3f}), '
          f'V_cell = {V:.3f}')
    print(f'#   plane waves = {npw}, k-points = {nk}, bands = {NSTATE}, '
          f'valence e- = {NELEC} (occ 2/band -> {nocc} filled)')

    eigen = np.zeros((NSTATE, nk))
    occup = np.zeros((NSTATE, nk))
    p_tm = np.zeros((NSTATE, NSTATE, 3, nk), dtype=complex)
    rvnl_tm = np.zeros((NSTATE, NSTATE, 3, nk), dtype=complex)   # local -> 0
    for ik in range(nk):
        H = wz.build_hamiltonian(kpoint[ik], Gcart, pos, spec)
        ev, evec = np.linalg.eigh(H)
        eigen[:, ik] = ev[:NSTATE]
        occup[:nocc, ik] = 2.0
        p_tm[:, :, :, ik] = epm_io.momentum_matrix(kpoint[ik], Gcart, evec[:, :NSTATE])
        if (ik + 1) % max(1, nk // 8) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik+1}/{nk}')

    # Reuse the proven non-orthogonal GS writer (sets sysname/material via globals).
    fcc.SYSNAME, fcc.MATERIAL = SYSNAME, MATERIAL
    fcc._write_gs_files_nonorth(
        outdir, kfrac, kpoint, kweight, eigen, occup, p_tm, rvnl_tm,
        note='PRIMITIVE wurtzite (4-atom hexagonal, non-orthogonal, unfolded -- no cosets)',
        b_matrix=np.array(Brec))
    np.savez(f'{outdir}{SYSNAME}_kcart.npz', kcart=kpoint, b_matrix=np.array(Brec),
             a=a_au, c=c_au)
    print(f'# EPM (CdS primitive): wrote {SYSNAME}_k/_eigen/_tm.data (nk={nk}, nb={NSTATE})')


def main_bandpath(outdir='./'):
    """Clean hexagonal-primitive band path (A-Gamma-M-K-Gamma), reusing the
    wurtzite reference's primitive path generator with the primitive sysname."""
    wz.main_bandpath(sysname=SYSNAME, outdir=outdir)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    if mode == 'gap':
        report_gap()
    elif mode == 'bandpath':
        main_bandpath()
    elif mode == 'gs':
        main_gs()
        main_bandpath()
    else:
        report_gap()
        print()
        main_gs()
        main_bandpath()
