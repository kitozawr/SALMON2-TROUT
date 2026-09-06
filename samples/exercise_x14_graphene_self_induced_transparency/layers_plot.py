#!/usr/bin/env python3
"""One layer against two, at two dopings: is a bilayer's transmission T_1 squared?

    python3 layers_plot.py --set 0.6 ef0.6_L1/runs/*_rt.data ef0.6_L2/runs/*_rt.data \
                           --set 0.4 ef0.4_L1/runs/*_rt.data ef0.4_L2/runs/*_rt.data \
                           [--t-meas 0.60 0.70] [--n-sub 1.65] [--out layers.png]

Why the question is not trivial. For N electronically decoupled sheets inside
d << lambda the driver adds their currents in the same local field
(`sbe_sheet_nlayers`), so the stack is one sheet of conductance N*sigma, NOT N
independent Fresnel interfaces. With z = Z0 sigma complex the sheet BC gives

    T_1 = |2/(2+z)|^2,  T_2 = |2/(2+2z)|^2,
    T_1^2 / T_2 = 16 |1+z|^2 / |2+z|^4 = 1 + [(Im z)^2 - (Re z)^2]/2 + O(z^3),

so the sign of the T*T error is set by the CHARACTER of the sheet: a dissipative
one (z real) makes T_1^2 an underestimate, a reactive one (z imaginary) an
overestimate. A doped graphene sheet at THz is an inductor, so squaring the
monolayer transmission comes out HIGH and inverting a measured bilayer as
sqrt(T_2) understates the monolayer conductance.

Left panel: T(E_0) for one and two layers at each doping, with T_1^2 dotted.
Right panel: the ratio T_1^2/T_2 against field -- the size and sign of the error.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import analyze, TruncatedRun, AU_EV, VF_EPM_AU, C_AU, SIGMA_UNIV  # noqa: E402
from drude_check import sheet_from_transmission  # noqa: E402
from transmission import run_variable  # noqa: E402

A0_PER_KVCM = 6.213e-4
Z0 = 4.0 * np.pi / C_AU


def collect(patterns):
    """(E0 [kV/cm], T, R, A, sigma_c) per run, sorted by field."""
    rows = []
    for f in sorted(sum((glob.glob(p) for p in patterns), [])):
        try:
            r, _ = analyze(f)
        except TruncatedRun as exc:
            print(f'# SKIPPED: {exc}')
            continue
        ring = (run_variable(f, 'yn_sbe_eph', 'n') or 'n').strip().strip("'\"").lower() == 'y'
        rows.append((r['E0_kvcm'], r['T'], r['R'], r['A'], r['sig_c'], r['nk'], ring))
    rows.sort(key=lambda r: r[0])
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--set', nargs='+', action='append', required=True, metavar='EF L1_GLOB L2_GLOB',
                    help='Fermi level in eV, then the one-layer and two-layer globs')
    ap.add_argument('--t-meas', nargs='*', type=float, default=[])
    ap.add_argument('--n-sub', type=float, default=1.65)
    ap.add_argument('--out', default='layers.png')
    args = ap.parse_args(argv)

    sets = []
    for spec in args.set:
        ef = float(spec[0])
        # each of the two may hold several whitespace-separated globs, so a series can
        # be assembled from a subset of a run directory (dropping a contaminated field)
        a, b = collect(spec[1].split()), collect(spec[2].split())
        if a and b:
            sets.append((ef, np.array([r[:4] for r in a]), np.array([r[:4] for r in b]),
                         [r[4] for r in a], a[0][5], a[0][6]))
    if not sets:
        print('# nothing matched'); return 1

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available'); return 0

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    cols = ('#c0392b', '#2980b9', '#16a085', '#8e44ad')
    for isr, ((ef, a, b, sg, nk, ring), c) in enumerate(zip(sets, cols)):
        kF = (abs(ef) / AU_EV) / VF_EPM_AU
        esat = kF / A0_PER_KVCM
        ax[0].semilogx(a[:, 0], a[:, 1], 'o-', color=c, label=f'1 layer, $E_F$ = {ef:g} eV')
        ax[0].semilogx(b[:, 0], b[:, 1], 's-', color=c, alpha=0.65, label=f'2 layers, $E_F$ = {ef:g} eV')
        ax[0].semilogx(a[:, 0], a[:, 1]**2, ':', color=c, lw=1.6,
                       label=f'$T_1^2$ (naive), $E_F$ = {ef:g} eV')
        ax[0].axvline(esat, ls='--', lw=0.9, color=c, alpha=0.5)
        ax[0].annotate(f'$A_0=k_F$ ({ef:g} eV)', xy=(esat, 0.985 - 0.055 * isr),
                       xycoords=('data', 'axes fraction'),
                       fontsize=7.5, color=c, ha='center', va='top')
        ratio = a[:, 1]**2 / b[:, 1]
        ax[1].semilogx(a[:, 0], 100.0 * (ratio - 1.0), 'o-', color=c,
                       label=f'$E_F$ = {ef:g} eV')
        zb = Z0 * SIGMA_UNIV * sg[0]
        cphi = zb.real / abs(zb)
        # decide the branch from the data, not from an assumption: Re/|z| near 0 is an
        # inductor (T_1^2 too high), near 1 a Drude absorber (T_1^2 too low)
        kind = ('reactive (inductor)  -> T_1^2 too HIGH' if cphi < 0.4 else
                'resistive (Drude)    -> T_1^2 too LOW' if cphi > 0.6 else
                'mixed                -> either sign possible')
        print(f'# E_F = {ef:g} eV: band-averaged z = {zb.real:+.3f}{zb.imag:+.3f}j '
              f'(Re/|z| = {cphi:.3f}) -> {kind}')
        print(f'{"E0":>7} {"T1":>8} {"T2":>8} {"T1^2":>8} {"T1^2/T2":>9}')
        for i in range(len(a)):
            print(f'{a[i, 0]:7.1f} {a[i, 1]:8.5f} {b[i, 1]:8.5f} {a[i, 1]**2:8.5f} {ratio[i]:9.4f}')
    for tm in args.t_meas:
        ax[0].axhline(tm / sheet_from_transmission(tm, args.n_sub)[1], ls=':', c='#27ae60', lw=1.2)
    if args.t_meas:
        ax[0].text(0.995, 0.30, 'green: measured monolayer, substrate divided out',
                   transform=ax[0].transAxes, fontsize=7.5, color='#27ae60', ha='right', va='bottom',
                   bbox=dict(fc='white', ec='none', alpha=0.75, pad=1.5))
    ax[0].set_xlabel('peak field $E_0$ [kV/cm]'); ax[0].set_ylabel('transmission $T$')
    ax[0].set_title('One layer, two layers, and the square of one', fontsize=9)
    ax[0].legend(fontsize=7, ncol=2, loc='lower left', framealpha=0.9); ax[0].grid(alpha=0.25)
    ax[1].axhline(0.0, c='k', lw=0.8)
    ax[1].set_xlabel('peak field $E_0$ [kV/cm]')
    ax[1].set_ylabel(r'error of $T\cdot T$:  $100\,(T_1^2/T_2 - 1)$  [%]')
    cph = [(Z0 * SIGMA_UNIV * s[3][0]).real / abs(Z0 * SIGMA_UNIV * s[3][0]) for s in sets]
    ttl = ('A reactive sheet makes $T_1^2$ an OVERestimate' if max(cph) < 0.4 else
           'A resistive (Drude) sheet makes $T_1^2$ an UNDERestimate' if min(cph) > 0.6 else
           'Sign of the $T\\cdot T$ error follows the phase of $\\sigma$')
    ax[1].set_title(ttl, fontsize=9)
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
    nks = ' / '.join(str(int(round(np.sqrt(s[4])))) + '^2' for s in sets if s[4])
    mode = ('with the e-ph + Rana ring' if all(s[5] for s in sets) else
            'coherent (no dissipation)' if not any(s[5] for s in sets) else
            'MIXED coherent/dissipative -- check the run set')
    fig.text(0.008, 0.008,
             'Free-standing calculation (no substrate); 2 layers = sbe_sheet_nlayers = 2, one shared '
             f'local field, NOT two Fresnel interfaces. Mesh {nks}, DAST transient + 100 fs '
             f'ring-down, {mode}.',
             fontsize=7.0, color='#34495e')
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
