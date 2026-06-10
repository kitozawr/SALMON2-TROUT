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

Optionally (INCLUDE_SPIN_ORBIT = True) the scalar problem is promoted to the
SPINOR problem of dimension 2*N_PW: the basis |G> is doubled to |G,s>,
s in {up,down} (spin-block ordering: first N_PW components spin-up), and

    H0_spinor(k) = H_loc(k) (x) 1_2  +  H_SO(k),

where H_SO is the projected Weisz / Bloom-Bergstresser spin-orbit operator
(Chelikowsky-Cohen form):

    <k+G,s|H_SO|k+G',s'> = -i [lam_S(G,G') S_S(G-G') + i lam_A(G,G') S_A(G-G')]
                              * [(k+G) x (k+G')] . sigma_ss',

with S_S = cos(dG.tau), S_A = sin(dG.tau), lam_S/A = (lam_cat +/- lam_an)/2,
lam_cat = mu B(K)B(K'), lam_an = alpha mu B(K)B(K') and a single fitted
constant mu, calibrated so that the Gamma8-Gamma7 splitting equals
Delta0 = 0.341 eV (GaAs).  The l=1 radial kernel is modelled as
B(K) = K / (1 + (K/zeta)^2)^3 (normalized so B(K)/K -> 1 as K -> 0).

Because the spin-orbit term is NONLOCAL, the velocity acquires the mandatory
correction v_SO = -i[r, H_SO] = grad_k H_SO(k); it is written to block 2
("rvnl_tm") of the SYSNAME_tm.data file, so the SBE run must use
yn_vnl_correction = 'y' (and yn_sbe_spinor = 'y') for a spinor dataset.
"""

import numpy as np
from scipy.linalg import eigh, eigvalsh

# =============================================================================
# Hardcoded run parameters (Cubic Cell Setup)
# =============================================================================
MATERIAL            = 'GaAs'
OUTPUT_DIR          = './'

A_LATTICE_AU        = 10.68         # Cubic lattice constant a [Bohr]
PW_CUTOFF_RY        = 11.1          # |G|^2 cutoff in (2*pi/a)^2 units

NUM_KGRID           = (4, 4, 4)     # Monkhorst-Pack grid for the CUBIC BZ

NSTATE              = 32            # 16 valence + 16 conduction (folded bands)
NELEC               = 32            # 8 atoms * 4 valence e- = 32 electrons

# -----------------------------------------------------------------------------
# Hardcoded spin-orbit (spinor) switch.
#   False : scalar N_PW problem, 32 bands, occupation 2 per band (legacy).
#   True  : spinor 2*N_PW problem, 64 spin-orbit-split bands, occupation 1
#           per band; mu is auto-calibrated to Delta0 at Gamma.
# -----------------------------------------------------------------------------
INCLUDE_SPIN_ORBIT  = True

SO_ALPHA            = 1.5           # anion/cation SO strength ratio lam_an/lam_cat
                                    # (As 4p splitting > Ga 4p splitting)
SO_ZETA_AU          = 2.0           # inverse-length scale of B(K) [Bohr^-1]
SO_MU_GUESS         = 1.0e-4        # initial guess for mu [Ha * Bohr^4]
SO_DELTA0_TARGET_EV = 0.341         # target Gamma8-Gamma7 splitting (GaAs)

SYSNAME = 'GaAs_cubic_so' if INCLUDE_SPIN_ORBIT else 'GaAs_cubic'

HARTREE_EV = 27.211386245988

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
# Spin-orbit (spinor) machinery
# =============================================================================
def so_radial_kernel(K, zeta=SO_ZETA_AU):
    """B(K) = K / (1 + (K/zeta)^2)^3 -- l=1 kernel with B(K)/K -> 1 as K -> 0.
    (Exact result for a hydrogenic-like radial p function R_1(r) ~ r e^{-zeta r}.)"""
    x = (K / zeta)**2
    return K / (1.0 + x)**3

def so_radial_kernel_deriv(K, zeta=SO_ZETA_AU):
    """dB/dK for the kernel above: (1 - 5x) / (1 + x)^4, x = (K/zeta)^2."""
    x = (K / zeta)**2
    return (1.0 - 5.0 * x) / (1.0 + x)**4

def spin_blocks(Ax, Ay, Az):
    """Assemble a 2npw x 2npw matrix from the sigma-contraction kappa.sigma:
    [[Az, Ax - i Ay], [Ax + i Ay, -Az]] (spin-block ordering: up first)."""
    return np.block([[Az, Ax - 1j * Ay], [Ax + 1j * Ay, -Az]])

def build_so_matrices(kvec, Gcart, a_lattice, with_velocity):
    """Per-unit-mu spin-orbit Hamiltonian H_SO/mu (2npw x 2npw) and, when
    with_velocity, the analytic velocity correction v_SO/mu = grad_k (H_SO/mu)
    as a list of three 2npw x 2npw matrices."""
    npw = Gcart.shape[0]
    twopi_over_a = 2.0 * np.pi / a_lattice
    G_idx = np.round(Gcart / twopi_over_a).astype(int)

    # Structure factors on dG = G - G' with the same FCC-in-cubic parity rule
    # as the local potential (8-atom structure factor vanishes otherwise).
    dG_idx = G_idx[:, None, :] - G_idx[None, :, :]
    parity = ((dG_idx[..., 0] - dG_idx[..., 1]) % 2 == 0) & \
             ((dG_idx[..., 1] - dG_idx[..., 2]) % 2 == 0)
    phase = np.pi / 4.0 * dG_idx.sum(axis=-1)
    SS = np.cos(phase)
    SA = np.sin(phase)

    # lam_S/A = (lam_cat +/- lam_an)/2 with lam_cat = mu B B', lam_an = alpha mu B B'
    c_s = 0.5 * (1.0 + SO_ALPHA)
    c_a = 0.5 * (1.0 - SO_ALPHA)
    # Pair prefactor (without mu, without B B'): -i [c_s S_S + i c_a S_A]
    F0 = -1j * (c_s * SS + 1j * c_a * SA) * parity

    K = kvec[None, :] + Gcart                      # (npw, 3)
    Kmag = np.linalg.norm(K, axis=1)
    B = so_radial_kernel(Kmag)
    BB = F0 * (B[:, None] * B[None, :])            # full pair factor / mu

    kappa = np.cross(K[:, None, :], K[None, :, :]) # (npw, npw, 3): (k+G) x (k+G')

    H_so = spin_blocks(BB * kappa[..., 0], BB * kappa[..., 1], BB * kappa[..., 2])

    if not with_velocity:
        return H_so, None

    # v_SO/mu = grad_k of the above. Two contributions per direction d:
    #   d/dk_d [B(|K|)B(|K'|)] = B'(|K|) Khat_d B(|K'|) + B(|K|) B'(|K'|) Khat'_d
    #   d/dk_d [(K x K')]      = e_d x (K' - K) = e_d x (G' - G)
    Bp = so_radial_kernel_deriv(Kmag)
    with np.errstate(invalid='ignore', divide='ignore'):
        Khat = np.where(Kmag[:, None] > 1e-12, K / Kmag[:, None], 0.0)
    dG_cart = Gcart[None, :, :] - Gcart[:, None, :]   # G' - G, (npw, npw, 3)

    v_so = []
    for d in range(3):
        dBB = F0 * (Bp[:, None] * Khat[:, None, d] * B[None, :]
                    + B[:, None] * Bp[None, :] * Khat[None, :, d])
        e_d = np.zeros(3)
        e_d[d] = 1.0
        dkappa = np.cross(e_d, dG_cart)               # (npw, npw, 3)
        comps = [dBB * kappa[..., c] + BB * dkappa[..., c] for c in range(3)]
        v_so.append(spin_blocks(*comps))
    return H_so, v_so

def build_hamiltonian_spinor(material, kvec, Gcart, a_lattice, mu, with_velocity=False):
    """H0_spinor = H_loc (x) 1_2 + mu * (H_SO/mu); optionally also v_SO."""
    H_loc = build_hamiltonian_sc(material, kvec, Gcart, a_lattice)
    H_so, v_so = build_so_matrices(kvec, Gcart, a_lattice, with_velocity)
    H = np.kron(np.eye(2), H_loc) + mu * H_so
    if with_velocity:
        return H, [mu * v for v in v_so]
    return H

def momentum_matrix_spinor(kvec, Gcart, evec):
    """Local momentum p (x) 1_2 in the spinor band basis (diagonal k+G in PW)."""
    nb = evec.shape[1]
    kpg = kvec[None, :] + Gcart
    p_mn = np.zeros((nb, nb, 3), dtype=complex)
    for idir in range(3):
        diag = np.concatenate([kpg[:, idir], kpg[:, idir]])
        p_mn[:, :, idir] = evec.conj().T @ (diag[:, None] * evec)
    return p_mn

def calibrate_so_mu(Gcart, a_lattice):
    """Fit mu so that the Gamma8-Gamma7 splitting at Gamma equals
    SO_DELTA0_TARGET_EV. Delta0 is (nearly) linear in mu, so a simple
    proportional update converges in a few iterations."""
    kgamma = np.zeros(3)
    target_ha = SO_DELTA0_TARGET_EV / HARTREE_EV
    mu = SO_MU_GUESS
    delta0 = 0.0
    for _ in range(30):
        H = build_hamiltonian_spinor(MATERIAL, kgamma, Gcart, a_lattice, mu)
        evals = eigvalsh(H)
        # Top six valence states at Gamma: Gamma7 doublet + Gamma8 quadruplet
        g7 = evals[NELEC-6:NELEC-4].mean()
        g8 = evals[NELEC-4:NELEC].mean()
        delta0 = g8 - g7
        if abs(delta0 - target_ha) < 1e-9:
            break
        mu *= target_ha / delta0
    eg = evals[NELEC] - evals[NELEC-1]
    print(f'#   SO calibration: mu = {mu:.8E} [a.u.]')
    print(f'#   Delta0 (Gamma8-Gamma7) = {delta0*HARTREE_EV:.4f} eV (target {SO_DELTA0_TARGET_EV} eV)')
    print(f'#   direct gap Eg(Gamma)   = {eg*HARTREE_EV:.4f} eV')
    return mu

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

    if INCLUDE_SPIN_ORBIT:
        nb = 2 * NSTATE     # spinor bands (spin-orbit split)
        nocc = NELEC        # one electron per spinor band
        occ_val = 1.0
    else:
        nb = NSTATE
        nocc = NELEC // 2   # 32 / 2 = 16 occupied bands
        occ_val = 2.0

    print(f'# EPM CUBIC CELL (Band-Folded) -- Python Reference')
    print(f'#   spin-orbit (spinor)= {INCLUDE_SPIN_ORBIT}')
    print(f'#   plane waves        = {npw} (Simple Cubic basis)'
          + (f' -> spinor dim {2*npw}' if INCLUDE_SPIN_ORBIT else ''))
    print(f'#   k-points           = {nk} (Cubic BZ)')
    print(f'#   bands requested    = {nb} / valence e- = {NELEC} (occ {occ_val} per band)')

    mu = calibrate_so_mu(Gcart, A_LATTICE_AU) if INCLUDE_SPIN_ORBIT else 0.0

    eigen = np.zeros((nb, nk))
    occup = np.zeros((nb, nk))
    p_tm = np.zeros((nb, nb, 3, nk), dtype=complex)
    rvnl_tm = np.zeros((nb, nb, 3, nk), dtype=complex)

    for ik in range(nk):
        if INCLUDE_SPIN_ORBIT:
            H, v_so = build_hamiltonian_spinor(MATERIAL, kpoint[ik], Gcart,
                                               A_LATTICE_AU, mu, with_velocity=True)
            assert np.allclose(H, H.conj().T), 'H_spinor not Hermitian'
            evals, evecs = eigh(H)

            evecs_nb = evecs[:, :nb]
            p_tm[:, :, :, ik] = momentum_matrix_spinor(kpoint[ik], Gcart, evecs_nb)
            # Nonlocal SO velocity correction v_SO = -i[r, H_SO] = grad_k H_SO,
            # in the band basis -> block 2 (rvnl_tm) of the tm file.
            for idir in range(3):
                rvnl_tm[:, :, idir, ik] = evecs_nb.conj().T @ (v_so[idir] @ evecs_nb)
        else:
            H = build_hamiltonian_sc(MATERIAL, kpoint[ik], Gcart, A_LATTICE_AU)
            evals, evecs = eigh(H)
            p_tm[:, :, :, ik] = momentum_matrix(kpoint[ik], Gcart, evecs[:, :nb])

        eigen[:, ik] = evals[:nb]
        occup[:nocc, ik] = occ_val
        occup[nocc:nb, ik] = 0.0

        if (ik + 1) % max(1, nk // 10) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik + 1}/{nk}')

    write_epm_files(SYSNAME, OUTPUT_DIR, MATERIAL, kpoint, kweight, eigen, occup,
                    p_tm, rvnl_tm, b_matrix)

# =============================================================================
# Output Writers (Byte-for-byte compatible with SALMON gs_info_ssbe)
# =============================================================================
def write_epm_files(sysname, outdir, material, kpoint, kweight, eigen, occup,
                    p_tm, rvnl_tm, b_matrix):
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
        if INCLUDE_SPIN_ORBIT:
            f.write('# block 2: rvnl_tm = -i[r,Vso] = grad_k H_SO (nonlocal spin-orbit velocity)\n')
        else:
            f.write('# block 2: rvnl_tm = -i[r,Vnl]  (all zero: local pseudopotential)\n')
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    px, py, pz = rvnl_tm[ib, jb, 0, ik], rvnl_tm[ib, jb, 1, ik], rvnl_tm[ib, jb, 2, ik]
                    f.write('{:6d}{:6d}{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                        ik + 1, ib + 1, jb + 1,
                        px.real, px.imag, py.real, py.imag, pz.real, pz.imag))

    print(f'# EPM (Cubic): wrote ground-state data files for sysname = {sysname}')

if __name__ == '__main__':
    main()
