#!/usr/bin/env python3
"""The x14 headline figure: transmission of a DOPED and an INTRINSIC graphene sheet
against the peak driving field, on the same mesh and the same pulse, plus the sheet
conductivity of the doped run against the value a measurement implies.

    python3 field_scan_plot.py --doped 'runs_doped/*/graphene_sit_sbe_rt.data' \
                               --intrinsic 'runs_intrinsic/*/graphene_sit_sbe_rt.data' \
                               --t-meas 0.60 0.70 --n-sub 1.65 --out doped_vs_intrinsic.png

Left panel: T(E0) for both occupations. Only the initial occupation differs, so the
two curves isolate what the doping carriers do. The intrinsic sheet darkens
monotonically (Landau-Zener pair creation, wiki/12 sec. 7); the doped sheet is flat
in the linear regime, darkens slightly, then BRIGHTENS once the vector-potential
excursion A_0 exceeds the Fermi radius k_F -- the drift velocity saturates at v_F
and the differential conductivity falls as k_F/A_0 (wiki/12 sec. 4a.3). The
saturation field A_0 = k_F is drawn from the run's own E_F.

Right panel: Re sigma of the doped sheet over the incident band, with the sheet
conductance(s) implied by the measured transmissions passed in --t-meas (the
substrate's Fresnel loss is divided out with --n-sub). This is the panel that
compares calculation and experiment: both are sheet conductances in units of
sigma_univ = e^2/4hbar, no fitting in between.

Both run sets must be the same mesh, pulse and dt -- otherwise the comparison is
between meshes, not between occupations.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import analyze, AU_EV, VF_EPM_AU  # noqa: E402
from drude_check import sheet_from_transmission, Z0, run_variable  # noqa: E402

SIGMA_UNIV = 0.25
A0_PER_KVCM = 6.213e-4      # peak |A| [a.u.] per kV/cm of the scaled DAST transient


def collect(patterns):
    """(E0 [kV/cm], T, R, A, Re sigma/sigma_univ, E_F [eV]) per run, sorted by field."""
    files = sorted(sum((glob.glob(p) for p in patterns), []))
    rows = []
    for f in files:
        r, _ = analyze(f)
        ef = float(run_variable(f, 'sbe_ef_ev', 0.0) or 0.0)
        rows.append((r['E0_kvcm'], r['T'], r['R'], r['A'], r['resig'], ef))
    rows.sort()
    return np.array([r[:5] for r in rows]), (rows[0][5] if rows else 0.0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--doped', nargs='+', required=True, help='*_sbe_rt.data of the doped scan (globs ok)')
    ap.add_argument('--intrinsic', nargs='*', default=[], help='*_sbe_rt.data of the undoped control scan')
    ap.add_argument('--t-meas', nargs='*', type=float, default=[], help='measured transmissions, substrate included')
    ap.add_argument('--series', nargs='*', default=[],
                    help='extra doped series to overlay, each "label:glob" (e.g. "with e-ph ring:diss/runs/*/..._rt.data"); '
                         'plotted on the T and sigma panels beside --doped')
    ap.add_argument('--n-sub', type=float, default=1.65, help='substrate index for --t-meas (PET ~ 1.65)')
    ap.add_argument('--out', default='doped_vs_intrinsic.png')
    ap.add_argument('--title', default=None)
    args = ap.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available -- no figure'); return 0

    dop, ef = collect(args.doped)
    if dop.size == 0:
        print('# no doped runs matched'); return 1
    ins, _ = collect(args.intrinsic) if args.intrinsic else (np.empty((0, 5)), 0.0)
    kF = (abs(ef) / AU_EV) / VF_EPM_AU if ef else 0.0
    e_sat = kF / A0_PER_KVCM if kF else None

    npan = 3 if ins.size else 2
    extra = []
    for spec in args.series:
        lab, _, pat = spec.partition(':')
        arr, _ = collect([pat])
        if arr.size:
            extra.append((lab, arr))
    fig, ax = plt.subplots(1, npan, figsize=(5.5 * npan, 4.0))
    ax[0].semilogx(dop[:, 0], dop[:, 1], 'o-', color='#c0392b',
                   label=f'doped, $E_F$ = {ef:g} eV (metal)')
    for (lab, arr), col, mk in zip(extra, ('#e67e22', '#16a085', '#8e44ad'), ('D-', 'v-', '^-')):
        ax[0].semilogx(arr[:, 0], arr[:, 1], mk, ms=5, color=col, label=lab)
    if ins.size:
        ax[0].semilogx(ins[:, 0], ins[:, 1], 's-', color='#2c3e50', label='intrinsic (semimetal)')
    if e_sat:
        ax[0].axvline(e_sat, ls='--', c='#7f8c8d', lw=1)
        ax[0].axvspan(e_sat, max(dop[:, 0].max(), 1.0), color='#f1c40f', alpha=0.13)
        ax[0].annotate('$A_0 = k_F$\n(current saturation)', xy=(e_sat, dop[:, 1].min()),
                       xytext=(e_sat * 1.8, dop[:, 1].min() - 0.06), fontsize=8,
                       arrowprops=dict(arrowstyle='->', lw=0.8, color='#7f8c8d'))
    zs_meas = []
    for tm in args.t_meas:
        zs, t_bare = sheet_from_transmission(tm, args.n_sub)
        zs_meas.append(zs / Z0 / SIGMA_UNIV)
        ax[0].axhline(tm / t_bare, ls=':', c='#27ae60', lw=1.2)
    if args.t_meas:
        ax[0].text(dop[0, 0] * 1.2, max(t / sheet_from_transmission(t, args.n_sub)[1] for t in args.t_meas) + 0.012,
                   'measured sheet transmission (substrate divided out)', fontsize=7.5, color='#27ae60')
    ax[0].set_xlabel('peak field $E_0$ [kV/cm]'); ax[0].set_ylabel('transmission $T$')
    ax[0].set_title(args.title or 'Same mesh, same pulse: only the initial occupation differs', fontsize=9)
    ax[0].legend(fontsize=8, loc='lower left'); ax[0].grid(alpha=0.25)

    # middle panel: the extinction the DOPING carriers add, 1 - T_doped/T_intrinsic.
    # This is the observable that isolates them: the intrinsic Landau-Zener darkening
    # divides out, and what is left is the Drude response of the doped carriers, which
    # collapses once the drift saturates.
    if ins.size:
        tin = {round(r[0], 3): r[1] for r in ins}
        common = np.array([[r[0], 1.0 - r[1] / tin[round(r[0], 3)]] for r in dop
                           if round(r[0], 3) in tin])
        if common.size:
            ax[1].semilogx(common[:, 0], 100.0 * common[:, 1], 'o-', color='#8e44ad')
            ax[1].axhline(0.0, c='k', lw=0.6)
            if e_sat:
                ax[1].axvline(e_sat, ls='--', c='#7f8c8d', lw=1)
                ax[1].axvspan(e_sat, common[:, 0].max(), color='#f1c40f', alpha=0.13)
            ax[1].set_xlabel('peak field $E_0$ [kV/cm]')
            ax[1].set_ylabel(r'extinction added by the doping, $1-T_{\rm doped}/T_{\rm intrinsic}$  [%]')
            ax[1].set_title('The doping carriers alone: their extinction\npeaks at the saturation field and collapses',
                            fontsize=9)
            ax[1].grid(alpha=0.25)
    ia = npan - 1
    ax[ia].semilogx(dop[:, 0], dop[:, 4], 'o-', color='#c0392b', label=r'Re $\sigma$, doped')
    for (lab, arr), col, mk in zip(extra, ('#e67e22', '#16a085', '#8e44ad'), ('D-', 'v-', '^-')):
        ax[ia].semilogx(arr[:, 0], arr[:, 4], mk, ms=5, color=col, label=f'Re $\\sigma$, {lab}')
    for z, tm in zip(zs_meas, args.t_meas):
        ax[ia].axhline(z, ls=':', c='#27ae60')
        ax[ia].text(dop[0, 0] * 1.2, z * 1.02, f'measured {z:.1f} $\\sigma_{{univ}}$ (T = {100 * tm:.0f} %)',
                    fontsize=7.5, color='#27ae60')
    if e_sat:
        ax[ia].axvline(e_sat, ls='--', c='#7f8c8d', lw=1)
    if len(dop) > 1:
        imax = int(np.argmax(dop[:, 4]))
        drop = 100.0 * (dop[-1, 4] / dop[imax, 4] - 1.0)
        ax[ia].set_title(rf'Sheet conductivity: {drop:+.0f} % from its peak to {dop[-1, 0]:.0f} kV/cm', fontsize=9)
    ax[ia].set_xlabel('peak field $E_0$ [kV/cm]'); ax[ia].set_ylabel(r'Re $\sigma$ / $\sigma_{univ}$')
    ax[ia].legend(fontsize=8); ax[ia].grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')

    print(f'{"E0[kV/cm]":>9} {"T doped":>9} {"T intr":>9} {"R doped":>9} {"A doped":>9} {"Re s/s0":>9}')
    tin = {round(r[0], 3): r[1] for r in ins} if ins.size else {}
    for r in dop:
        ti = tin.get(round(r[0], 3))
        print(f'{r[0]:9.2f} {r[1]:9.5f} ' + (f'{ti:9.5f} ' if ti is not None else f'{"--":>9} ')
              + f'{r[2]:9.2e} {r[3]:9.5f} {r[4]:9.3f}')
    if e_sat:
        print(f'# k_F = {kF:.5f} a.u. (E_F = {ef:g} eV) -> current saturation A_0 = k_F at E_0 = {e_sat:.1f} kV/cm')
    return 0


if __name__ == '__main__':
    sys.exit(main())
