#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone Python/NumPy/SciPy reference implementation of the local EPM 
for the CONVENTIONAL CUBIC CELL of GaAs (8 atoms, 32 valence electrons).

This script uses a mathematical "band-folding" trick: it builds the Hamiltonian 
in a simple cubic plane-wave basis, but applies a parity selection rule to the 
Cohen-Bergstresser form factors. This naturally and exactly folds the 4 primitive 
FCC Brillouin zones into 1 simple cubic zone, yielding 16 valence bands and 16+ 
conduction bands, perfectly matching a DFT/SBE setup that uses the cubic cell.
"""

import numpy as np
from scipy.linalg import eigh

# =============================================================================
# Hardcoded run parameters (Cubic Cell Setup)
# =============================================================================
MATERIAL            = 'GaAs'
SYSNAME             = 'GaAs_cubic'
OUTPUT_DIR          = './'

A_LATTICE_AU        = 10.68         # Cubic lattice constant a [Bohr]
PW_CUTOFF_RY        = 11.1          # |G|^2 cutoff in (2*pi/a)^2 units

NUM_KGRID           = (4, 4, 4)     # Monkhorst-Pack grid for the CUBIC BZ

NSTATE              = 32            # 16 valence + 16 conduction (folded bands)
NELEC               = 32            # 8 atoms * 4 valence e- = 32 electrons

# =============================================================================
# Cohen-Bergstresser (1966) form factors (Ry -> Ha)
# =============================================================================
RY_TO_HA = 0.5
_CB_FORM_FACTORS_RY = {
    3:  (-0.23,  0.07),
    4:  ( 0.00,  0.05),
    8:  ( 0.01,  0.00),
    11: ( 0.06,  0.01),
}

def form_factors(material, G2):
    if material != 'GaAs':
        raise ValueError("Only 'GaAs' supported")
    vs_ry, va_ry = _CB_FORM_FACTORS_RY.get(G2, (0.0, 0.0))
    return vs_ry * RY_TO_HA, va_ry * RY_TO_HA

# =============================================================================
# Simple Cubic Lattice & Basis Generation
# =============================================================================
def lattice_vectors_sc(a):
    """Conventional simple cubic vectors."""
    return np.array([a, 0, 0]), np.array([0, a, 0]), np.array([0, 0, a])

def reciprocal_lattice_sc(a):
    """Simple cubic reciprocal lattice."""
    twopi = 2.0 * np.pi
    b1 = np.array([twopi/a, 0, 0])
    b2 = np.array([0, twopi/a, 0])
    b3 = np.array([0, 0, twopi/a])
    return b1, b2, b3, a**3

def build_plane_wave_basis_sc(a_lattice, cutoff_ry):
    """Generate simple cubic G-vectors. Cutoff is on |G|^2 in (2pi/a)^2 units."""
    twopi_over_a = 2.0 * np.pi / a_lattice
    nmax = int(np.ceil(np.sqrt(cutoff_ry))) + 1
    
    Gcart_list, G2_list = [], []
    for h in range(-nmax, nmax + 1):
        for k in range(-nmax, nmax + 1):
            for l in range(-nmax, nmax + 1):
                g2_units = h**2 + k**2 + l**2
                if g2_units <= cutoff_ry + 1.0e-8:
                    Gcart_list.append(twopi_over_a * np.array([h, k, l]))
                    G2_list.append(g2_units)
                    
    return np.array(Gcart_list), np.array(G2_list, dtype=int)

def monkhorst_pack_grid(b_matrix, num_kgrid):
    n1, n2, n3 = num_kgrid
    nk = n1 * n2 * n3
    kpoint = np.zeros((nk, 3))
    kweight = np.full(nk, 1.0 / nk)
    ik = 0
    for i1 in range(1, n1 + 1):
        for i2 in range(1, n2 + 1):
            for i3 in range(1, n3 + 1):
                f1 = (2 * i1 - n1 - 1) / (2.0 * n1)
                f2 = (2 * i2 - n2 - 1) / (2.0 * n2)
                f3 = (2 * i3 - n3 - 1) / (2.0 * n3)
                kpoint[ik] = f1 * b_matrix[0] + f2 * b_matrix[1] + f3 * b_matrix[2]
                ik += 1
    return kpoint, kweight

# =============================================================================
# Hamiltonian with Parity Selection Rule (The Band-Folding Trick)
# =============================================================================
def build_hamiltonian_sc(material, kvec, Gcart, a_lattice):
    npw = Gcart.shape[0]
    H = np.zeros((npw, npw), dtype=complex)
    
    twopi_over_a = 2.0 * np.pi / a_lattice
    # Extract integer (h,k,l) indices for all G vectors
    G_indices = np.round(Gcart / twopi_over_a).astype(int)

    kpg = kvec[None, :] + Gcart
    diag = 0.5 * np.einsum('ij,ij->i', kpg, kpg)
    np.fill_diagonal(H, diag)

    for i in range(npw):
        for j in range(i + 1, npw):
            dG_idx = G_indices[i] - G_indices[j]
            h, k, l = dG_idx
            
            # PARITY SELECTION RULE: 
            # For an FCC lattice embedded in a simple cubic supercell, the structure 
            # factor is EXACTLY ZERO unless h, k, l are all even or all odd.
            # This naturally forces the band folding without manual 8-atom summation!
            if (h % 2 == k % 2) and (k % 2 == l % 2):
                dG2 = h**2 + k**2 + l**2
                VS, VA = form_factors(material, dG2)
                if VS == 0.0 and VA == 0.0:
                    continue
                
                # Phase corresponds to tau = (a/8)(1,1,1)
                phase = np.pi / 4.0 * (h + k + l)
                val = complex(VS * np.cos(phase), VA * np.sin(phase))
                H[i, j] = val
                H[j, i] = np.conj(val)
                
    return H

def momentum_matrix(kvec, Gcart, evec):
    npw, nb = evec.shape
    kpg = kvec[None, :] + Gcart
    p_mn = np.zeros((nb, nb, 3), dtype=complex)
    for idir in range(3):
        Dc = kpg[:, idir][:, None] * evec
        p_mn[:, :, idir] = evec.conj().T @ Dc
    return p_mn

# =============================================================================
# Main Execution
# =============================================================================
def main():
    a1, a2, a3 = lattice_vectors_sc(A_LATTICE_AU)
    b1, b2, b3, volume = reciprocal_lattice_sc(A_LATTICE_AU)
    b_matrix = np.array([b1, b2, b3])

    Gcart, G2 = build_plane_wave_basis_sc(A_LATTICE_AU, PW_CUTOFF_RY)
    npw = Gcart.shape[0]

    kpoint, kweight = monkhorst_pack_grid(b_matrix, NUM_KGRID)
    nk = kpoint.shape[0]

    nb = NSTATE
    nocc = NELEC // 2  # 32 / 2 = 16 occupied bands

    print(f'# EPM CUBIC CELL (Band-Folded) -- Python Reference')
    print(f'#   plane waves        = {npw} (Simple Cubic basis)')
    print(f'#   k-points           = {nk} (Cubic BZ)')
    print(f'#   bands requested    = {nb} (16 val + 16 cond) / valence e- = {NELEC}')

    eigen = np.zeros((nb, nk))
    occup = np.zeros((nb, nk))
    p_tm = np.zeros((nb, nb, 3, nk), dtype=complex)

    for ik in range(nk):
        H = build_hamiltonian_sc(MATERIAL, kpoint[ik], Gcart, A_LATTICE_AU)
        evals, evecs = eigh(H)              

        eigen[:, ik] = evals[:nb]
        occup[:nocc, ik] = 2.0
        occup[nocc:nb, ik] = 0.0

        p_mn = momentum_matrix(kpoint[ik], Gcart, evecs[:, :nb])
        p_tm[:, :, :, ik] = p_mn

        if (ik + 1) % max(1, nk // 10) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik + 1}/{nk}')

    write_epm_files(SYSNAME, OUTPUT_DIR, MATERIAL, kpoint, kweight, eigen, occup, p_tm, b_matrix)

# =============================================================================
# Output Writers (Byte-for-byte compatible with SALMON gs_info_ssbe)
# =============================================================================
def write_epm_files(sysname, outdir, material, kpoint, kweight, eigen, occup, p_tm, b_matrix):
    nk = kpoint.shape[0]
    nb = eigen.shape[0]

    # Convert Cartesian k-vectors [a.u.] to reduced (dimensionless) coordinates,
    # matching DFT output: kx_red = kx / b11, etc. (orthogonal-cell convention).
    b_diag = np.array([b_matrix[0, 0], b_matrix[1, 1], b_matrix[2, 2]])
    kpoint_red = kpoint / b_diag

    with open(f'{outdir}{sysname}_k.data', 'w') as f:
        f.write('# k-point data\n')
        f.write('# generated by EPM (Cubic Cell / Band-Folded) -- Python reference\n')
        f.write(f'# material = {material}, nk = {nk}\n')
        f.write('# units: kx,ky,kz [reduced, dimensionless], weight (sums to 1)\n')
        f.write('# ik, kx, ky, kz, weight\n')
        for ik in range(nk):
            f.write('{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                ik + 1, kpoint_red[ik, 0], kpoint_red[ik, 1], kpoint_red[ik, 2], kweight[ik]))

    with open(f'{outdir}{sysname}_eigen.data', 'w') as f:
        f.write('# eigenvalue data\n')
        f.write('# generated by EPM (Cubic Cell / Band-Folded) -- Python reference\n')
        f.write(f'# nk = {nk:6d}, nb = {nb:6d}\n')
        for ik in range(nk):
            f.write(f'# ik = {ik + 1:6d}\n')
            for ib in range(nb):
                f.write('{:6d}{:18.10E}{:18.10E}\n'.format(ib + 1, eigen[ib, ik], occup[ib, ik]))

    with open(f'{outdir}{sysname}_tm.data', 'w') as f:
        f.write('# transition matrix data\n')
        f.write('# generated by EPM (Cubic Cell / Band-Folded) -- Python reference\n')
        f.write('# block 1: p_tm = <u_m|p|u_n>  (ik, ib, jb, Re px, Im px, Re py, Im py, Re pz, Im pz)\n')
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    px, py, pz = p_tm[ib, jb, 0, ik], p_tm[ib, jb, 1, ik], p_tm[ib, jb, 2, ik]
                    f.write('{:6d}{:6d}{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                        ik + 1, ib + 1, jb + 1,
                        px.real, px.imag, py.real, py.imag, pz.real, pz.imag))
        f.write('# block 2: rvnl_tm = -i[r,Vnl]  (all zero: local pseudopotential)\n')
        zeros9 = '{:6d}{:6d}{:6d}' + '{:18.10E}' * 6 + '\n'
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    f.write(zeros9.format(ik + 1, ib + 1, jb + 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    print(f'# EPM (Cubic): wrote ground-state data files for sysname = {sysname}')

if __name__ == '__main__':
    main()
