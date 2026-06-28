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
