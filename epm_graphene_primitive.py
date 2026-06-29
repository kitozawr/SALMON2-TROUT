"""
Graphene PRIMITIVE-CELL EPM ground state -- 2-atom HEXAGONAL primitive,
NON-orthogonal (a1,a2 at 120 deg), 2D sheet embedded in a 3D vacuum box, NO
folding. The graphene companion to epm_graphene.py (which runs the SBE on the
2-fold-folded orthorhombic 4-atom rectangular supercell); this module runs on
the genuine 2-atom hexagonal primitive cell -- the same unfolded test as
GaAs/Si/CdS, on graphene's gapless Dirac system.

Geometry (reused verbatim from epm_graphene): a1=a(1,0), a2=a(1/2,sqrt3/2),
2 carbon atoms at +/- delta/2 (bond center origin), minimal Ramanujam pi-EPM
(3 form factors). The Dirac cone is the lowest band PAIR (1 pi e-/atom), so the
2-atom primitive cell has nelec=2 (occ 2 -> 1 filled valence pi band) and the
SBE uses nstate=2 (the pi/pi* Dirac pair, the standard 2-band graphene model).
struct_norm=1 (the validated 2-atom primitive convention).

The strict-2D primitive is hexagonal NON-orthogonal in-plane; we embed it in a
3D box with GRAPHENE_VACUUM vacuum along z (no kz dispersion -> num_kgrid z = 1),
write the (non-orthogonal) reciprocal vectors into the k.data header, and reuse
the proven non-orthogonal GS writer. There is NO 2-fold rect folding, NO unfold
map. Local pseudopotential -> rvnl_tm = 0.

Usage:
    python3 epm_graphene_primitive.py            # validate + write GS + bandpath
    python3 epm_graphene_primitive.py gap        # Dirac gap at K (should be ~0)
    python3 epm_graphene_primitive.py gs         # GS files + bandpath
"""
import sys

import numpy as np

import epm_graphene as gr          # hexagonal-primitive geometry, H, pi form factors
import epm_io
import epm_gaas_primitive as fcc   # reuse the proven non-orthogonal GS writer

SYSNAME = 'graphene_prim'
MATERIAL = 'graphene'
VACUUM_ANG = 20.0                  # c-axis vacuum embedding the 2D sheet
NUM_KGRID = (12, 12, 1)            # multiple of 3 -> the Dirac point K=(2/3,1/3) is on the mesh
NELEC = 2                          # 2 atoms x 1 pi e- -> 1 filled valence pi band
NSTATE = 2                         # pi (valence) + pi* (conduction): the Dirac pair
CUTOFF_EV = 40.0


def _cell_3d_au():
    """3D primitive cell [a.u.]: hexagonal a1,a2 in-plane (Bohr), vacuum on z."""
    a1, a2 = gr.lattice_vectors()                       # 2D, Angstrom
    A2B = gr.ANG_TO_BOHR
    a1_3d = np.array([a1[0], a1[1], 0.0]) * A2B
    a2_3d = np.array([a2[0], a2[1], 0.0]) * A2B
    a3_3d = np.array([0.0, 0.0, VACUUM_ANG]) * A2B
    return a1_3d, a2_3d, a3_3d


def report_gap():
    """Dirac gap at K (should be ~0) and the band touching."""
    Gcart2d, _ = gr.build_pw_basis(CUTOFF_EV)
    tau = gr.basis_atoms()
    hs = gr.high_symmetry_points()
    for name in ('Gamma', 'K'):
        k2d = np.array(hs[name]) if name in hs else None
        ev = gr.bands_at_k(hs[name], Gcart2d, tau, nb=4)
        print(f'  {name:6s}: E(pi)={ev[0]:7.3f}  E(pi*)={ev[1]:7.3f}  gap = {ev[1]-ev[0]:7.4f} eV')


def main_gs(outdir='./'):
    a1, a2, a3 = _cell_3d_au()
    b_matrix, V = fcc.reciprocal_vectors(a1, a2, a3)        # general 3D reciprocal (rows)
    Gcart2d, _ = gr.build_pw_basis(CUTOFF_EV)
    tau = gr.basis_atoms()
    npw = len(Gcart2d)
    if NSTATE > npw:
        raise ValueError(f'NSTATE={NSTATE} exceeds npw={npw}; raise CUTOFF_EV')
    # 3D plane waves [a.u.]: in-plane G (Ang^-1 -> Bohr^-1), G_z = 0
    Gcart_au = np.zeros((npw, 3))
    Gcart_au[:, 0:2] = Gcart2d / gr.ANG_TO_BOHR
    kpoint, kweight = epm_io.monkhorst_pack(b_matrix, NUM_KGRID)
    kfrac = kpoint @ np.linalg.inv(b_matrix)
    nk = len(kpoint)
    nocc = NELEC // 2
    print(f'# EPM graphene PRIMITIVE (2-atom hexagonal, non-orthogonal, 2D-in-vacuum, NO folding)')
    print(f'#   plane waves = {npw}, k-points = {nk}, bands = {NSTATE}, '
          f'valence e- = {NELEC} (occ 2/band -> {nocc} filled)')

    eigen = np.zeros((NSTATE, nk))
    occup = np.zeros((NSTATE, nk))
    p_tm = np.zeros((NSTATE, NSTATE, 3, nk), dtype=complex)
    rvnl_tm = np.zeros((NSTATE, NSTATE, 3, nk), dtype=complex)   # local -> 0
    for ik in range(nk):
        k2d = kpoint[ik, 0:2] * gr.ANG_TO_BOHR                  # back to Ang^-1 (in-plane)
        H = gr.build_hamiltonian(k2d, Gcart2d, tau, struct_norm=1.0)   # [eV], primitive
        ev, evec = np.linalg.eigh(H)
        eigen[:, ik] = ev[:NSTATE] / gr.HA_TO_EV                # eV -> Ha
        occup[:nocc, ik] = 2.0
        p_tm[:, :, :, ik] = epm_io.momentum_matrix(kpoint[ik], Gcart_au, evec[:, :NSTATE])
        if (ik + 1) % max(1, nk // 8) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik+1}/{nk}')

    fcc.SYSNAME, fcc.MATERIAL = SYSNAME, MATERIAL
    fcc._write_gs_files_nonorth(
        outdir, kfrac, kpoint, kweight, eigen, occup, p_tm, rvnl_tm,
        note='PRIMITIVE graphene (2-atom hexagonal, non-orthogonal, 2D-in-vacuum, no folding)',
        b_matrix=np.array(b_matrix))
    np.savez(f'{outdir}{SYSNAME}_kcart.npz', kcart=kpoint, b_matrix=np.array(b_matrix))
    print(f'# EPM (graphene primitive): wrote {SYSNAME}_k/_eigen/_tm.data (nk={nk}, nb={NSTATE})')


def main_bandpath(outdir='./'):
    """Clean primitive Dirac band path (Gamma-M-K-Gamma), reusing the reference's
    primitive path generator with the primitive sysname."""
    gr.main_bandpath(sysname=SYSNAME, outdir=outdir)


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
