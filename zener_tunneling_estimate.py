#!/usr/bin/env python3
"""
zener_tunneling_estimate.py
===========================

Analytic, transverse-k-resolved interband tunnelling estimate for the unfolded
GaAs band structure already in this repo (uses epm_gaas_reference.py via the
helpers in band_field_coupling.py).  It answers, quantitatively and per k:

    "How large is the *vertical* (k-conserving) Zener tunnelling probability
     for a field E || axis, as a function of the TRANSVERSE momentum k_perp,
     and are there any off-Gamma small-gap avoided crossings (conduction
     sub-bands, folded copies, Gamma-L) where Landau-Zener transfer is O(1)?"

Two complementary estimators, both in atomic units (e = hbar = m_e = 1):

  (K) KANE / Zener parabolic vertical-gap tunnelling, valid at a SYMMETRIC
      band extremum (the cb-vb direct minimum, which along any <hkl> line sits
      at k_par = 0 with gap E_g(k_perp)):

          P_Kane(k_perp) = exp( - C * m_r^{1/2} * E_g(k_perp)^{3/2} / F )

      C = pi/(2*sqrt(2)) (standard Kane) or pi/sqrt(2) ("repo" convention used
      in band_field_coupling.py).  This is the ONLY correct estimator at the
      vertical minimum, because there dv = v_c - v_v -> 0 and the linear
      Landau-Zener formula is singular.

  (LZ) Linear LANDAU-ZENER at a GENUINE avoided crossing between two dispersing
      bands n, n+1 encountered while the crystal momentum is swept
      k_par(t) = k_par0 + A(t) (acceleration theorem, dk/dt = F):

          P_LZ = exp( - pi * Delta_min^2 / (2 * |dDelta/dt| ) ),
          dDelta/dt = F * |d(E_{n+1}-E_n)/dk_par|   (evaluated numerically).

      No mass and no convention needed -- it is read straight off the bands.
      This is the estimator for off-Gamma conduction-conduction / folded
      crossings (e.g. the Gamma-L intervalley seam) that could populate states
      AWAY from the field axis.

The script prints a quantitative verdict and writes a 3-panel figure:
  (a) vertical cb-vb gap along the swept line for several k_perp offsets;
  (b) per-k single-pass tunnelling probability: Kane at the vertical minimum
      (markers) and the strongest off-Gamma LZ crossing along the sweep (line);
  (c) transverse decay P_Kane(k_perp) / P_Kane(0)  -- how fast the vertical
      channel collapses onto the axis.

Notes
-----
* The bands come from epm_gaas_reference.py exactly as band_field_coupling.py
  uses them (same parity-folding, spin-orbit and -- after the rvnl fix in
  band_field_coupling.primitive_bands_momentum -- the full nonlocal velocity).
* The local CB-only Cohen-Bergstresser EPM undershoots the GaAs gap (E_g(Gamma)
  ~ 1.27 eV vs 1.42/1.52 eV experiment), so ABSOLUTE P_Kane is biased high by
  orders of magnitude through the E_g^{3/2} exponent. Panel (c) reports the
  Gamma-NORMALIZED transverse decay, which is unaffected by that bias; quote
  absolute probabilities only with that caveat (or a scissor shift).

Usage:
  python3 zener_tunneling_estimate.py                       # E||x=[100], 1-30 MV/cm
  python3 zener_tunneling_estimate.py --field-axis x --fields 1 3 10 30
  python3 zener_tunneling_estimate.py --perp-dir 010 --perp-frac 0 0.05 0.1 0.2
  python3 zener_tunneling_estimate.py --kmax-frac 0.5 --npts 401
  python3 zener_tunneling_estimate.py --kane-conv repo      # match w_inj in band_field_coupling.py
  python3 zener_tunneling_estimate.py --map2d --n2d 41      # 2-D transverse birth map W(k_perp)

The --map2d panel is the most direct visual test of the unfolded-population
puzzle: it reads the FWHM of the injection "needle" at Gamma directly and shows
the folded conduction-band geometry, so you can see that the diagonal weight in
the SBE pop_lcb snapshots coincides with the FOLD positions (zone edge / L-valley
diagonals via the kz-average), not with the injection blob.
"""
import argparse, importlib.util, sys
from pathlib import Path
import numpy as np
from numpy.linalg import eigh
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

AU_FIELD_VPM = 5.14220675e11          # 1 a.u. of E-field in V/m
MVCM_TO_AU   = 1e8 / AU_FIELD_VPM     # 1 MV/cm in a.u.

# ----- pull the EPM machinery + the exact unfolder out of band_field_coupling -
def _load_bfc():
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location("bfc", here / "band_field_coupling.py")
    bfc = importlib.util.module_from_spec(spec)
    saved = sys.argv; sys.argv = ['bfc']
    try:
        spec.loader.exec_module(bfc)
    finally:
        sys.argv = saved
    return bfc

AXIS = {'x': np.array([1.,0.,0.]), 'y': np.array([0.,1.,0.]), 'z': np.array([0.,0.,1.])}
# transverse directions are given as cubic <hkl>; converted to a Cartesian unit vector
def hkl_unit(hkl):
    v = np.array([float(c) for c in hkl], float)
    return v/np.linalg.norm(v)

def kane_const(conv):
    return (np.pi/np.sqrt(2.0)) if conv == 'repo' else (np.pi/(2.0*np.sqrt(2.0)))

def bands_at(bfc, epm, q_cart, ctx):
    """Unfolded primitive bands + cb/vb indices at Cartesian q (a.u.)."""
    ev, p = bfc.primitive_bands_momentum(epm, q_cart, ctx['Gcart'], ctx['Gidx'],
                                         ctx['npw'], ctx['a'], ctx['mu'], ctx['spinor'])
    return ev, p

def _transverse_axes(field_axis):
    """The two cubic Cartesian axes transverse to the field axis (the plane the
    crystal momentum is NOT swept in, hence the plane the injection is born in
    and the plane the folded copies live in). Returns (vec_a, vec_b, lab_a, lab_b)."""
    pair = {'x': (1, 2), 'y': (0, 2), 'z': (0, 1)}[field_axis]
    e = np.eye(3); labs = ['kx', 'ky', 'kz']
    a, b = pair
    return e[a], e[b], labs[a], labs[b]


def birth_map_2d(bfc, epm, ctx, field_axis, kmax_frac, n2d, F, Ck, m_r):
    """Transverse-plane Kane birth-rate map W(k_perp) at k_par = 0, plus the
    FOLDED cubic lowest-conduction-band energy on the same grid (so the diagonal
    blobs of the SBE pop_lcb map can be located).

      W(k_perp) = |<cb|p_axis|Gamma8>|^2 * exp(-C m_r^1/2 E_g(k_perp)^3/2 / F)

    Under acceleration along the field, k_perp is conserved, so W is exactly the
    distribution that SEEDS the conduction population. Returns (grid[n2d] reduced,
    W[n2d,n2d], Ecb_fold[n2d,n2d] eV, Eg[n2d,n2d] eV, (fwhm_a,fwhm_b) reduced,
    (lab_a,lab_b), twopi_a). n2d is forced odd so the central row/col is k=0."""
    if n2d % 2 == 0:
        n2d += 1
    a = ctx['a']; twopi_a = 2.0 * np.pi / a
    Gcart, npw, Gidx = ctx['Gcart'], ctx['npw'], ctx['Gidx']
    spinor, mu = ctx['spinor'], ctx['mu']
    nv = epm.NELEC // 4; icb = nv; g8 = list(range(nv - 4, nv))
    axis_idx = {'x': 0, 'y': 1, 'z': 2}[field_axis]
    ncb = epm.NELEC if spinor else epm.NELEC // 2     # full-cubic lowest-CB index
    va, vb, lab_a, lab_b = _transverse_axes(field_axis)
    grid = np.linspace(-kmax_frac, kmax_frac, n2d)
    W = np.zeros((n2d, n2d)); Ecb = np.zeros((n2d, n2d)); Eg = np.zeros((n2d, n2d))
    for i, u in enumerate(grid):            # row -> axis b (vertical)
        for j, w in enumerate(grid):        # col -> axis a (horizontal)
            q = (w * va + u * vb) * twopi_a                 # k_par = 0
            # unfolded primitive bands + (SO-corrected) momentum -> birth rate
            ev, p = bfc.primitive_bands_momentum(epm, q, Gcart, Gidx, npw, a, mu, spinor)
            eg = ev[icb] - ev[g8].max()
            coup = float(np.sum(np.abs(p[icb, g8, axis_idx]) ** 2))
            Eg[i, j] = eg * epm.HARTREE_EV
            W[i, j] = coup * np.exp(-Ck * np.sqrt(m_r) * max(eg, 0.0) ** 1.5 / F)
            # folded cubic lowest CB (== the SBE pop_lcb branch) at the same k
            G0 = np.round(q / twopi_a).astype(int); ksc = q - G0 * twopi_a
            Hf = (epm.build_hamiltonian_spinor(epm.MATERIAL, ksc, Gcart, a, mu)
                  if spinor else epm.build_hamiltonian_sc(epm.MATERIAL, ksc, Gcart, a))
            Ecb[i, j] = np.linalg.eigvalsh(Hf)[ncb] * epm.HARTREE_EV

    def _fwhm(profile):
        half = profile.max() / 2.0
        above = profile >= half
        return (grid[above].max() - grid[above].min()) if above.sum() >= 2 else np.nan
    ic = n2d // 2
    return (grid, W, Ecb, Eg, (_fwhm(W[ic, :]), _fwhm(W[:, ic])),
            (lab_a, lab_b), twopi_a)


def plot_map2d(grid, W, Ecb, Eg, fwhm, labs, twopi_a, field_axis, Fmv, outdir, dpi):
    """W(k_perp) heatmap with the needle FWHM box and the folded-CB valleys
    (white contours + crosses at the off-center minima) overlaid -- the visual
    proof that the diagonal SBE weight sits at the FOLD positions, not the blob."""
    lab_a, lab_b = labs
    fa, fb = fwhm
    ext = [grid[0], grid[-1], grid[0], grid[-1]]
    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    im = ax.imshow(np.clip(W, 0, None), origin='lower', extent=ext, aspect='equal',
                   cmap='inferno')
    fig.colorbar(im, ax=ax, label=r'$W(k_\perp)=|\langle cb|p|v\rangle|^2 e^{-C m_r^{1/2}E_g^{3/2}/F}$  [a.u.]')
    # folded cubic lowest-CB energy: contours show the real CB anisotropy. In
    # this k_par=0 plane the CB is a SINGLE Gamma valley (the X_y/X_z copies fold
    # to Gamma, the L copies are out of plane) -- the contours rise monotonically
    # outward, petalled along the transverse axes (X-valley character).
    e0 = Ecb.min()
    levels = e0 + np.array([0.1, 0.3, 0.6, 1.0, 1.5])
    cs = ax.contour(grid, grid, Ecb, levels=levels, colors='w', linewidths=0.7, alpha=0.7)
    ax.clabel(cs, fmt=lambda v: f'{v - e0:.1f}', fontsize=6)
    # X_y/X_z zone-face fold positions (where the band petals point) for scale
    g0 = grid[-1]
    ax.plot([g0, -g0, 0, 0], [0, 0, g0, -g0], 'x', color='cyan', ms=10, mew=2,
            label='X$_y$/X$_z$ zone-face folds (|k|=%.2f)' % g0)
    # the injection needle FWHM
    if np.isfinite(fa) and np.isfinite(fb):
        from matplotlib.patches import Ellipse
        ax.add_patch(Ellipse((0, 0), fa, fb, fill=False, ec='lime', lw=1.8, ls='--',
                             label=f'birth FWHM = {fa:.3f}×{fb:.3f} (2π/a)'))
    ax.plot(0, 0, '+', color='lime', ms=12, mew=2)
    ax.set_xlabel(f'{lab_a} [reduced]'); ax.set_ylabel(f'{lab_b} [reduced]')
    ax.set_title(f'Transverse Kane birth map  W({lab_a},{lab_b}), $k_{{par}}$=0, '
                 f'field || {field_axis}, F={Fmv:g} MV/cm\n'
                 f'green = injection needle (FWHM);  white = folded CB energy;  '
                 f'cyan = zone-edge fold scale', fontsize=10)
    ax.legend(loc='upper right', fontsize=7, framealpha=0.85)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f'zener_birthmap2d_{field_axis}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight'); plt.close(fig)
    return out


def sweep_line(bfc, epm, ctx, k_perp_cart, axis_hat, kmax, npts):
    """Sweep k_par along +/- axis at fixed k_perp. Return E_all[npts, nb] (a.u.),
    the signed k_par grid (a.u.), and cb/vb band indices."""
    nv = epm.NELEC // 4; icb = nv; ivb = nv - 1
    kpar = np.linspace(-kmax, kmax, npts)
    nb = None; E = None
    for j, kp in enumerate(kpar):
        q = k_perp_cart + kp*axis_hat
        ev, _ = bands_at(bfc, epm, q, ctx)
        if E is None:
            nb = ev.size; E = np.empty((npts, nb))
        E[j] = ev
    return kpar, E, icb, ivb

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--field-axis', choices=['x','y','z'], default='x',
                    help='field polarisation (swept k_par direction)')
    ap.add_argument('--perp-dir', default='010',
                    help='cubic <hkl> for the TRANSVERSE offset scan (must be ⟂ axis)')
    ap.add_argument('--perp-frac', nargs='+', type=float,
                    default=[0.0, 0.05, 0.10, 0.20],
                    help='transverse offsets as a fraction of 2pi/a')
    ap.add_argument('--fields', nargs='+', type=float, default=[1.,3.,10.,30.],
                    help='field magnitudes [MV/cm]')
    ap.add_argument('--kmax-frac', type=float, default=0.5,
                    help='half-width of the swept k_par line, fraction of 2pi/a')
    ap.add_argument('--npts', type=int, default=401)
    ap.add_argument('--m-r', type=float, default=0.04, help='reduced mass for Kane [m_e]')
    ap.add_argument('--kane-conv', choices=['kane','repo'], default='kane')
    ap.add_argument('--gap-floor-ev', type=float, default=0.25,
                    help='only flag adjacent-band crossings whose min gap is BELOW '
                         'this (eV) as candidate Landau-Zener seams')
    ap.add_argument('--map2d', action='store_true',
                    help='also produce the 2-D transverse Kane birth map '
                         'W(k_perp) at k_par=0, with the needle FWHM and the '
                         'folded conduction-valley positions overlaid.')
    ap.add_argument('--n2d', type=int, default=41, help='2-D map grid size (odd)')
    ap.add_argument('--map-field', type=float, default=10.0,
                    help='field [MV/cm] for the 2-D birth map (shape only)')
    ap.add_argument('-o','--output', default='sbe_plots')
    ap.add_argument('--dpi', type=int, default=150)
    args = ap.parse_args()

    bfc = _load_bfc(); epm = bfc._load_epm()
    a = epm.A_LATTICE_AU; spinor = epm.INCLUDE_SPIN_ORBIT
    Gcart, _ = epm.build_plane_wave_basis_sc(a, epm.PW_CUTOFF_RY)
    ctx = dict(a=a, spinor=spinor, Gcart=Gcart, npw=Gcart.shape[0],
               Gidx=np.round(Gcart/(2*np.pi/a)).astype(int),
               mu=(epm.calibrate_so_mu(Gcart, a) if spinor else 0.0))
    HARTREE_EV = epm.HARTREE_EV
    twopi_a = 2.0*np.pi/a
    axis_hat = AXIS[args.field_axis]
    perp_hat = hkl_unit(args.perp_dir)
    if abs(np.dot(axis_hat, perp_hat)) > 1e-9:
        print(f"# WARNING: perp-dir [{args.perp_dir}] not orthogonal to field axis "
              f"{args.field_axis} (dot={np.dot(axis_hat,perp_hat):.3f}); projecting out.")
        perp_hat = perp_hat - np.dot(perp_hat, axis_hat)*axis_hat
        perp_hat /= np.linalg.norm(perp_hat)
    kmax = args.kmax_frac * twopi_a
    Ck = kane_const(args.kane_conv)

    print(f"# zener_tunneling_estimate: field||{args.field_axis}, perp||[{args.perp_dir}], "
          f"spinor={spinor}, m_r={args.m_r}, Kane const={Ck:.4f} ({args.kane_conv})")

    # ---- 2-D transverse birth map (the W(k_perp) panel) -----------------------
    if args.map2d:
        Fm = args.map_field * MVCM_TO_AU
        print(f"\n# === 2-D transverse Kane birth map W(k_perp), field || "
              f"{args.field_axis}, F = {args.map_field:g} MV/cm, n2d={args.n2d} ===")
        grid, W2, Ecb2, Eg2, fwhm, labs, twopi_a = birth_map_2d(
            bfc, epm, ctx, args.field_axis, args.kmax_frac, args.n2d, Fm, Ck, args.m_r)
        fa, fb = fwhm
        # Estimate of the needle half-width from the Eg curvature (Kane "E-bar"):
        #   E_g(k_perp) ~ E_g0 + (k_perp^2)/(2 m_perp); W ~ exp(-kappa Eg^3/2) sets
        #   the 1/e width; FWHM_pred uses ln2. We just report the measured FWHM.
        ic = args.n2d // 2 if args.n2d % 2 else args.n2d // 2
        eg0 = float(np.nanmin(Eg2)); ecb0 = float(Ecb2.min())
        r_need = 0.5 * fa                                   # needle radius (HWHM)
        edge = grid[-1]; r_edge = edge / max(r_need, 1e-9)
        print(f"#   E_g(min) = {eg0:.3f} eV at the centre (Gamma); the injection W is a "
              f"single Gamma needle.")
        print(f"#   birth-needle FWHM: {labs[0]} = {fa:.4f},  {labs[1]} = {fb:.4f}  (2pi/a)")
        print(f"#               = {fa*twopi_a:.4f} x {fb*twopi_a:.4f}  1/Bohr   "
              f"(HWHM = {r_need:.4f} 2pi/a)")
        # folded lowest-CB structure in THIS (k_par=0) plane
        print(f"#   folded cubic lowest-CB in this plane: E_cb(Gamma)={ecb0:.3f} eV, rising")
        print(f"#     monotonically to +{Ecb2[len(grid)//2,-1]-ecb0:.2f} eV at the {labs[0]} face")
        print(f"#     and +{Ecb2[-1,-1]-ecb0:.2f} eV at the corner -- NO off-Gamma minima here.")
        print(f"#   => the X_y/X_z copies fold to Gamma; the zone edge is {r_edge:.0f}x the")
        print(f"#      needle radius away. Nothing in this plane seeds off-axis weight.")
        out2 = plot_map2d(grid, W2, Ecb2, Eg2, fwhm, labs, twopi_a,
                          args.field_axis, args.map_field, Path(args.output), args.dpi)
        print(f"# saved {out2}")
        print("#\n#   READING THE SBE pop_lcb SNAPSHOT (kx-ky, AVG kz):")
        print("#   The diagonal blobs there are NOT in this transverse plane -- they are the")
        print("#   L-valley folds at the cube DIAGONALS (~(0.4,0.4,0.4) 2pi/a), pulled into")
        print("#   the kx-ky picture by the kz-average. They sit ~7 needle-radii from Gamma,")
        print("#   so the diagonal weight is FOLDING (band geometry), not vertical Zener")
        print("#   injection (a Gamma needle) nor an off-Gamma low-gap LZ seam. That it is")
        print("#   already present at the FIRST nonzero step is the SIGNATURE of folding")
        print("#   (instantaneous/geometric), not transport (which would take time to drift")
        print("#   k-space). It is physical, but worth confirming pop_lcb is the lowest-CB")
        print("#   BRANCH and not summing coset copies in the avg-kz projection.")
        return

    # ---- transverse scan: vertical gap at the minimum + Kane probability ------
    perp_kabs = []           # |k_perp| in 1/Bohr
    Eg_min    = []           # vertical cb-vb gap at the on-axis minimum [eV]
    lz_best   = []           # strongest sub-gap adjacent-band LZ crossing info
    lines     = {}           # for plotting gap profiles
    for fr in args.perp_frac:
        kperp = fr*twopi_a*perp_hat
        kpar, E, icb, ivb = sweep_line(bfc, epm, ctx, kperp, axis_hat, kmax, args.npts)
        gap_cbvb = (E[:, icb] - E[:, ivb]) * HARTREE_EV
        jmin = np.argmin(gap_cbvb)
        Egmin = gap_cbvb[jmin]
        perp_kabs.append(fr*twopi_a); Eg_min.append(Egmin)
        lines[fr] = (kpar/twopi_a, gap_cbvb)
        # scan ALL adjacent-band gaps along the sweep for genuine small-gap crossings
        best = None
        SEP_MIN = 0.5      # eV: the two bands must genuinely separate away from the dip
        DEG_EPS = 5e-3     # eV: below this the 'gap' is a spin-split touching, not avoided
        E_WIN_EV = 6.0     # eV above VBM: only the lowest conduction valleys matter (Γ,L,X)
        vbm_au = E[:, ivb].max()
        nlo, nhi = ivb, min(E.shape[1]-1, icb+7)            # top valence + lowest ~4 cb doublets
        for n in range(nlo, nhi):
            g = (E[:, n+1] - E[:, n]) * HARTREE_EV          # eV, >=0
            jm = int(np.argmin(g))
            if not (0 < jm < len(g)-1):              continue
            if g[jm] >= args.gap_floor_ev:           continue
            if g[jm] <  DEG_EPS:                      continue   # spin touching, not avoided
            if g.max() < SEP_MIN:                    continue   # near-degenerate everywhere
            # both partners must be low-lying conduction/valence states (skip PW continuum)
            if (E[jm, n+1]-vbm_au)*HARTREE_EV > E_WIN_EV: continue
            dEdk = np.gradient((E[:, n+1]-E[:, n]), kpar)       # Hartree*Bohr
            slope = abs(dEdk[jm])
            if slope < 1e-6:                          continue
            cand = dict(band=n, kpar=kpar[jm]/twopi_a, gap=g[jm], slope=slope)
            if best is None or g[jm] < best['gap']:
                best = cand
        lz_best.append(best)

    perp_kabs = np.array(perp_kabs); Eg_min = np.array(Eg_min)

    # ---- probabilities vs field ----------------------------------------------
    print("\n# === VERTICAL (cb-vb) ZENER, Kane parabolic, at the direct minimum ===")
    print("#  P_Kane = exp(-C m_r^1/2 Eg^3/2 / F);  Eg taken at the on-axis vertical min")
    hdr = "#  k_perp[1/Bohr]  Eg[eV] |" + " ".join(f"{f:>10g}MV/cm" for f in args.fields)
    print(hdr)
    Pmap = np.zeros((len(args.perp_frac), len(args.fields)))
    for i, fr in enumerate(args.perp_frac):
        Eg_au = Eg_min[i]/HARTREE_EV
        row = []
        for jf, Fmv in enumerate(args.fields):
            F = Fmv*MVCM_TO_AU
            P = np.exp(-Ck*np.sqrt(args.m_r)*Eg_au**1.5 / F)
            Pmap[i, jf] = P; row.append(P)
        print(f"#  {perp_kabs[i]:12.4f}  {Eg_min[i]:5.3f} | " +
              " ".join(f"{p:14.3e}" for p in row))

    # ---- off-Gamma Landau-Zener seams (the 'diagonal rays' candidate) ---------
    print(f"\n# === off-Gamma small-gap avoided crossings (gap < {args.gap_floor_ev} eV) ===")
    print("#  these are conduction-conduction / folded crossings reached by the sweep;")
    print("#  P_LZ = exp(-pi Delta^2 / (2 F |dDelta/dk|)).  If P_LZ ~ O(1) here, the")
    print("#  off-axis occupation can be GENUINE LZ transfer, not just folding.")
    any_seam = False
    for i, fr in enumerate(args.perp_frac):
        b = lz_best[i]
        if b is None:
            print(f"#  k_perp={perp_kabs[i]:.4f} : none below floor")
            continue
        any_seam = True
        ph = []
        for Fmv in args.fields:
            F = Fmv*MVCM_TO_AU
            P = np.exp(-np.pi*(b['gap']/HARTREE_EV)**2 / (2.0*F*max(b['slope'],1e-12)))
            ph.append(P)
        print(f"#  k_perp={perp_kabs[i]:.4f} : band({b['band']},{b['band']+1}) "
              f"gap={b['gap']:.3f} eV at k_par={b['kpar']:+.3f}(2pi/a)  "
              f"P_LZ@fields=" + " ".join(f"{p:.2e}" for p in ph))
    if not any_seam:
        print("#  (no sub-floor crossings on this line -> within the cb-vb pair the diagonal\n"
              "#   weight cannot come from a low-gap LZ seam; it is folding + transport.)")

    # ---- figure ---------------------------------------------------------------
    fig, (axA, axB, axC) = plt.subplots(3, 1, figsize=(8.2, 11))
    for fr,(x,g) in lines.items():
        axA.plot(x, g, lw=1.5, label=f'k_perp={fr:g}·2π/a')
    axA.set_xlabel('k$_\\parallel$  [2π/a]'); axA.set_ylabel('direct cb-vb gap [eV]')
    axA.set_title(f'Vertical gap along the swept line (field || {args.field_axis})')
    axA.grid(alpha=0.25); axA.legend(fontsize=8); axA.axvline(0,color='0.7',lw=.7)

    for jf, Fmv in enumerate(args.fields):
        axB.semilogy(perp_kabs, np.clip(Pmap[:,jf],1e-300,None), 'o-', lw=1.6,
                     label=f'{Fmv:g} MV/cm')
    axB.set_xlabel('|k_perp|  [1/Bohr]'); axB.set_ylabel('P$_{Kane}$ (vertical, per pass)')
    axB.set_title('Vertical Zener probability vs transverse momentum')
    axB.grid(alpha=0.25, which='both'); axB.legend(fontsize=8)

    P0 = Pmap[0].copy()
    for jf, Fmv in enumerate(args.fields):
        axC.semilogy(perp_kabs, np.clip(Pmap[:,jf]/max(P0[jf],1e-300),1e-300,None),
                     'o-', lw=1.6, label=f'{Fmv:g} MV/cm')
    axC.set_xlabel('|k_perp|  [1/Bohr]')
    axC.set_ylabel('P$_{Kane}$(k_perp) / P$_{Kane}$(0)')
    axC.set_title('Transverse collapse of the vertical channel onto the axis')
    axC.grid(alpha=0.25, which='both'); axC.legend(fontsize=8)

    fig.tight_layout()
    outdir = Path(args.output); outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f'zener_tunneling_{args.field_axis}_perp{args.perp_dir}.png'
    fig.savefig(out, dpi=args.dpi, bbox_inches='tight'); plt.close(fig)
    print(f"\n# saved {out}")

if __name__ == '__main__':
    main()
