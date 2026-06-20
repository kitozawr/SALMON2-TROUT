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

import sys

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

# -----------------------------------------------------------------------------
# Unfolded (primitive-cell) band path mode:  python3 epm_gaas_reference.py bandpath
#
# The cubic 8-atom cell is a supercell of 4 primitive FCC cells, so its bands
# are the primitive bands FOLDED 4-fold: every cubic k carries the states of
# 4 primitive BZ points, and along any path the conduction manifold shows 4
# overlaid copies of CB1/CB2/CB3 -- dense crossings that are an artifact of
# the supercell representation, not of the physics (cf. band unfolding,
# Quan-Rybin-Scheffler-Carbogno, PRB 113, 085112 (2026); Popescu-Zunger).
#
# Because the folding is EXACT here (the parity selection rule makes H block-
# diagonal over the 4 FCC reciprocal sublattices to machine precision), the
# primitive bands can be recovered exactly: for a primitive k = k_sc + G0 we
# extract the sublattice block {G : G - G0 in the FCC reciprocal set} of the
# cubic Hamiltonian at k_sc and diagonalize it. This mode writes
# SYSNAME_bandpath.data with the clean primitive bands along the path below
# (plotted by plot_sbe_results.py, including a Dresselhaus spin-splitting
# panel in the spinor case). It does NOT touch the MP-grid SBE dataset.
# -----------------------------------------------------------------------------
BANDPATH_LABELS = ['L', 'Gamma', 'X', 'W', 'K', 'Gamma']
BANDPATH_NDIV   = 40      # k-points per path segment
BANDPATH_NB     = 24      # states written per primitive k (12 in scalar mode)

# High-symmetry points in REDUCED coordinates of the FCC primitive reciprocal
# basis b1=(2pi/a)(-1,1,1), b2=(2pi/a)(1,-1,1), b3=(2pi/a)(1,1,-1).
HS_POINTS_FCC_PRIM = {
    'Gamma': (0.000, 0.000, 0.000),
    'X':     (0.000, 0.500, 0.500),
    'L':     (0.500, 0.500, 0.500),
    'W':     (0.250, 0.500, 0.750),
    'K':     (0.375, 0.375, 0.750),
    'U':     (0.250, 0.625, 0.625),
}

HARTREE_EV = 27.211386245988

# =============================================================================
# Cohen-Bergstresser (1966) form factors (Ry -> Ha)
# =============================================================================
RY_TO_HA = 0.5
# GaAs (zincblende): symmetric V^S and antisymmetric V^A form factors, keyed by
# |G|^2 in (2pi/a)^2 units. [M. L. Cohen & T. K. Bergstresser, Phys. Rev. 141,
# 789 (1966)]
_CB_FORM_FACTORS_RY = {
    3:  (-0.23,  0.07),
    4:  ( 0.00,  0.05),
    8:  ( 0.01,  0.00),
    11: ( 0.06,  0.01),
}
# Silicon (diamond): two IDENTICAL atoms per primitive cell -> the antisymmetric
# structure factor vanishes, so V^A == 0 exactly for all shells. Default set is
# Kunikiyo; the Cohen-Bergstresser set is provided for validation.
#   Kunikiyo: V^S(3)=-0.2258, V^S(8)=+0.05698, V^S(11)=+0.070709 Ry
#     [T. Kunikiyo et al., J. Appl. Phys. 75, 297 (1994), Table I]
#   Cohen-Bergstresser (alt): V^S(3)=-0.21, V^S(8)=+0.04, V^S(11)=+0.08 Ry
#     [Cohen & Bergstresser, Phys. Rev. 141, 789 (1966)]
# V^A(3)=V^A(4)=V^A(11)=0 exactly (diamond). [Cohen-Bergstresser 1966]
_SI_FORM_FACTORS_KUNIKIYO_RY = {
    3:  (-0.2258,   0.0),
    8:  (+0.05698,  0.0),
    11: (+0.070709, 0.0),
}
_SI_FORM_FACTORS_CB_RY = {
    3:  (-0.21, 0.0),
    8:  (+0.04, 0.0),
    11: (+0.08, 0.0),
}

def form_factors(material, G2):
    """Local-EPM symmetric/antisymmetric form factors (Hartree) for shell |G|^2.
    GaAs (zincblende): Cohen-Bergstresser V^S, V^A. Si (diamond): V^A == 0
    identically (two identical atoms), Kunikiyo V^S by default."""
    if material == 'GaAs':
        table = _CB_FORM_FACTORS_RY
    elif material in ('Si', 'Si_kunikiyo'):
        table = _SI_FORM_FACTORS_KUNIKIYO_RY
    elif material == 'Si_cb':
        table = _SI_FORM_FACTORS_CB_RY
    else:
        raise ValueError("material must be 'GaAs', 'Si' (Kunikiyo) or 'Si_cb'")
    vs_ry, va_ry = table.get(G2, (0.0, 0.0))
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
# Unfolded primitive-cell band path (exact sublattice-block extraction)
# =============================================================================
def fcc_reciprocal_rows(a_lattice):
    """FCC primitive reciprocal vectors as rows [Cartesian, a.u.]."""
    c = 2.0 * np.pi / a_lattice
    return c * np.array([[-1.0, 1.0, 1.0],
                         [ 1.0, -1.0, 1.0],
                         [ 1.0, 1.0, -1.0]])

def sublattice_mask(G_indices, g0_idx):
    """Plane waves G with G - G0 in the FCC reciprocal set (all components of
    (h,k,l) of equal parity) -- the exact primitive block at k_sc + G0."""
    d = G_indices - g0_idx[None, :]
    return ((d[:, 0] - d[:, 1]) % 2 == 0) & ((d[:, 1] - d[:, 2]) % 2 == 0)

def generate_bandpath(mu, Gcart, a_lattice):
    """Compute the UNFOLDED primitive-cell bands along BANDPATH_LABELS by
    diagonalizing the exact FCC-sublattice blocks of the cubic Hamiltonian.
    Also asserts the block-diagonality (folding correctness) at every point."""
    twopi_over_a = 2.0 * np.pi / a_lattice
    G_indices = np.round(Gcart / twopi_over_a).astype(int)
    B = fcc_reciprocal_rows(a_lattice)
    npw = Gcart.shape[0]

    nb_out = BANDPATH_NB if INCLUDE_SPIN_ORBIT else BANDPATH_NB // 2
    ne_prim = NELEC // 4                      # electrons per primitive 2-atom cell
    nv = ne_prim if INCLUDE_SPIN_ORBIT else ne_prim // 2

    # Build the sampled path (primitive reduced -> Cartesian)
    qpts, dists, node_d = [], [], [0.0]
    cum = 0.0
    for iseg in range(len(BANDPATH_LABELS) - 1):
        qa = np.array(HS_POINTS_FCC_PRIM[BANDPATH_LABELS[iseg]])
        qb = np.array(HS_POINTS_FCC_PRIM[BANDPATH_LABELS[iseg + 1]])
        ka, kb = qa @ B, qb @ B
        seg_len = np.linalg.norm(kb - ka) / twopi_over_a
        last = (iseg == len(BANDPATH_LABELS) - 2)
        for s in range(BANDPATH_NDIV + (1 if last else 0)):
            t = s / BANDPATH_NDIV
            qpts.append(qa + t * (qb - qa))
            dists.append(cum + t * seg_len)
        cum += seg_len
        node_d.append(cum)

    max_offblock = 0.0
    rows_e = []
    for q in qpts:
        kcart = q @ B
        # Fold into the cubic BZ: k_sc = kcart - G0, G0 a simple-cubic G vector
        g0_idx = np.round(kcart / twopi_over_a).astype(int)
        ksc = kcart - twopi_over_a * g0_idx
        msk = sublattice_mask(G_indices, g0_idx)
        if INCLUDE_SPIN_ORBIT:
            H = build_hamiltonian_spinor(MATERIAL, ksc, Gcart, a_lattice, mu)
            rows = np.concatenate([np.where(msk)[0], np.where(msk)[0] + npw])
        else:
            H = build_hamiltonian_sc(MATERIAL, ksc, Gcart, a_lattice)
            rows = np.where(msk)[0]
        others = np.setdiff1d(np.arange(H.shape[0]), rows)
        max_offblock = max(max_offblock, np.abs(H[np.ix_(rows, others)]).max())
        evals = eigvalsh(H[np.ix_(rows, rows)])
        rows_e.append(evals[:nb_out])

    print(f'#   folding self-check: max off-sublattice |H| = {max_offblock:.3e} '
          f'(must be ~0: folding is exact)')
    assert max_offblock < 1e-12, 'parity selection rule violated -- folding broken'

    fname = f'{OUTPUT_DIR}{SYSNAME}_bandpath.data'
    with open(fname, 'w') as f:
        f.write('# unfolded primitive-cell band path -- EPM cubic reference\n')
        f.write(f'# spinor = {1 if INCLUDE_SPIN_ORBIT else 0}\n')
        f.write(f'# nv = {nv}\n')
        f.write(f'# nb = {nb_out}\n')
        f.write('# nodes: ' + '  '.join(
            f'{l} {d:.7f}' for l, d in zip(BANDPATH_LABELS, node_d)) + '\n')
        f.write('# ik, dist [2pi/a], q1, q2, q3 (fcc-primitive reduced), E_1..E_nb [Ha]\n')
        for ik, (q, d, e) in enumerate(zip(qpts, dists, rows_e)):
            f.write('{:6d}{:14.7f}{:10.5f}{:10.5f}{:10.5f}'.format(
                ik + 1, d, q[0], q[1], q[2]))
            f.write(''.join(f'{x:18.10E}' for x in e) + '\n')
    print(f'# EPM bandpath: wrote {fname} '
          f'({len(qpts)} k-points, {nb_out} bands, nv = {nv})')

def main_bandpath():
    """Fast standalone mode: only the unfolded band path, no MP-grid dataset."""
    Gcart, G2 = build_plane_wave_basis_sc(A_LATTICE_AU, PW_CUTOFF_RY)
    print(f'# EPM unfolded band path (primitive cell via exact sublattice blocks)')
    print(f'#   spin-orbit = {INCLUDE_SPIN_ORBIT}, path = {"-".join(BANDPATH_LABELS)}')
    mu = calibrate_so_mu(Gcart, A_LATTICE_AU) if INCLUDE_SPIN_ORBIT else 0.0
    generate_bandpath(mu, Gcart, A_LATTICE_AU)

# =============================================================================
# Unfold map: cubic supercell band -> (FCC sublattice, primitive band index)
# =============================================================================
SUBLATTICE_OFFSETS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

def main_unfoldmap():
    """
    Standalone mode:  python3 epm_gaas_reference.py unfoldmap

    Re-diagonalizes the MP-grid Hamiltonians (deterministic -- identical
    eigenvalues to the SBE dataset) and writes SYSNAME_unfold.data mapping
    every cubic band (ik, ib) to its FCC sublattice isub (1..4, i.e. the
    primitive BZ point k_prim = k_sc + G0(isub)) and its primitive band index
    ibprim (energy rank within that sublattice). Because the folding is exact,
    each eigenvector lives wholly in one sublattice (weight ~1) except at
    accidental cross-sublattice degeneracies, where the dominant weight is
    taken. Much cheaper than the full dataset generation: no momentum/tm
    matrices and no large file output.

    The SBE reads this map (if present next to the GS files) to output the
    population of PHYSICAL primitive bands (spins summed) at the unfolded
    primitive k-points instead of energy-ordered supercell branches.
    """
    a1, a2, a3 = lattice_vectors_sc(A_LATTICE_AU)
    b1, b2, b3, _ = reciprocal_lattice_sc(A_LATTICE_AU)
    b_matrix = np.array([b1, b2, b3])
    Gcart, G2 = build_plane_wave_basis_sc(A_LATTICE_AU, PW_CUTOFF_RY)
    twopi_over_a = 2.0 * np.pi / A_LATTICE_AU
    G_indices = np.round(Gcart / twopi_over_a).astype(int)
    npw = Gcart.shape[0]

    kpoint, kweight = monkhorst_pack_grid(b_matrix, NUM_KGRID)
    nk = kpoint.shape[0]
    nb = 2 * NSTATE if INCLUDE_SPIN_ORBIT else NSTATE
    ne_prim = NELEC // 4
    nv_prim = ne_prim if INCLUDE_SPIN_ORBIT else ne_prim // 2

    print(f'# EPM unfold map: {nk} k-points, {nb} bands -> 4 sublattices, nv_prim = {nv_prim}')
    mu = calibrate_so_mu(Gcart, A_LATTICE_AU) if INCLUDE_SPIN_ORBIT else 0.0

    masks = [sublattice_mask(G_indices, off) for off in SUBLATTICE_OFFSETS]
    # Row/column index sets of the 4 FCC-sublattice blocks (spinor: plane wave
    # G and its spin partner G+npw). H is exactly block-diagonal over these, so
    # diagonalizing a block yields the primitive bands at k_sc + G0(s) and a
    # symmetric per-sublattice energy ranking (used instead of an argmax counter).
    block_idx = []
    for msk in masks:
        idx = np.where(msk)[0]
        if INCLUDE_SPIN_ORBIT:
            idx = np.concatenate([idx, idx + npw])
        block_idx.append(np.sort(idx))

    n_ambig = 0
    fname = f'{OUTPUT_DIR}{SYSNAME}_unfold.data'
    with open(fname, 'w') as f:
        f.write('# unfold map (cubic band -> FCC sublattice spectral weights)\n')
        f.write('# nk, nb, nv_prim\n')
        f.write(f'{nk:8d}{nb:8d}{nv_prim:8d}\n')
        f.write('# isub, offset G0 (sc reduced)\n')
        for isub, off in enumerate(SUBLATTICE_OFFSETS):
            f.write(f'{isub + 1:4d}{off[0]:4d}{off[1]:4d}{off[2]:4d}\n')
        # isub/ibprim fix the physical-band slot (valence rank); w1..w4 are the
        # spectral weights |<psi|P_s|psi>|^2 of the cubic band on each of the 4
        # FCC sublattices (sum = 1). The SBE distributes the band population over
        # the sublattices by these weights -- exact at symmetry degeneracies,
        # where a hard argmax would dump it all on one equivalent primitive point.
        f.write('# ik, ib, isub, ibprim, w1, w2, w3, w4\n')
        for ik in range(nk):
            if INCLUDE_SPIN_ORBIT:
                H = build_hamiltonian_spinor(MATERIAL, kpoint[ik], Gcart,
                                             A_LATTICE_AU, mu)
            else:
                H = build_hamiltonian_sc(MATERIAL, kpoint[ik], Gcart, A_LATTICE_AU)
            evals, evecs = eigh(H)
            w2 = np.abs(evecs[:, :nb])**2
            wsub = np.zeros((4, nb))
            for s, msk in enumerate(masks):
                if INCLUDE_SPIN_ORBIT:
                    wsub[s] = w2[np.concatenate([np.where(msk)[0],
                                                 np.where(msk)[0] + npw])].sum(axis=0)
                else:
                    wsub[s] = w2[msk].sum(axis=0)
            # Normalise columns to exactly 1 (guard against basis incompleteness).
            wsub /= np.maximum(wsub.sum(axis=0), 1e-300)
            isub_b = np.argmax(wsub, axis=0)
            wmax_b = wsub[isub_b, np.arange(nb)]
            n_ambig += int((wmax_b < 0.99).sum())
            # Per-sublattice block spectra for a SYMMETRIC energy ranking: the
            # primitive-band rank of cubic band ib is its position in the energy
            # spectrum of its dominant sublattice block. Unlike an argmax counter
            # this is identical for symmetry-equivalent (degenerate) states, so
            # the physical level (VB-1/VB/CB1/CB2) is assigned consistently.
            block_evals = [np.sort(eigvalsh(H[np.ix_(bi, bi)])) for bi in block_idx]
            for ib in range(nb):
                s = isub_b[ib]
                ibprim = 1 + int(np.searchsorted(block_evals[s], evals[ib] - 1e-9))
                f.write('{:6d}{:6d}{:4d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n'.format(
                    ik + 1, ib + 1, s + 1, ibprim,
                    wsub[0, ib], wsub[1, ib], wsub[2, ib], wsub[3, ib]))
            if (ik + 1) % max(1, nk // 10) == 0 or ik == nk - 1:
                print(f'#   ... mapped k-point {ik + 1}/{nk}')
    print(f'# EPM unfold map: wrote {fname}'
          + (f' ({n_ambig} bands with cross-sublattice degeneracy, '
             f'distributed by spectral weight)'
             if n_ambig else ' (all weights > 0.99)'))

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
    if len(sys.argv) > 1 and sys.argv[1] == 'bandpath':
        main_bandpath()   # fast: unfolded primitive bands only, no MP dataset
    elif len(sys.argv) > 1 and sys.argv[1] == 'unfoldmap':
        main_unfoldmap()  # band -> sublattice map for the existing MP dataset
    else:
        main()
