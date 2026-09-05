#!/usr/bin/env python3
"""Semiclassical picture of a doped Dirac cone in a velocity-gauge field.

    python3 cone_kinematics.py [--ef-ev 0.6] [--drive DAST_E100kVcm.txt] [--out cone_kinematics.png]

What the solver actually does, drawn. In the velocity gauge the canonical label k does
NOT move; the field enters only through the Hamiltonian

    H_k(A) = v_F sigma . (k + A(t)),      E_pm = +- v_F |k + A|,
    v_pm   = +- v_F (k + A) / |k + A|,    |v| = v_F ALWAYS.

Three things follow, and each is a panel here.

1. WHERE THE FIELD PUSHES.  A(t) = -int E dt' displaces the KINETIC momentum k + A.
   The occupied disc stays put in label space; what moves through it is the point of
   instantaneous degeneracy, k = -A, where the two branches touch.

2. WHERE THE ELECTRON TURNS.  The velocity of a conduction state points radially AWAY
   from that moving degeneracy point. At A = 0 the arrows point outward from the
   centre and cancel. As |A| grows they all swing towards +A, and once |A| >> k_F they
   are essentially parallel: the drift velocity cannot exceed v_F, so the current
   saturates at n e v_F and the sheet brightens (wiki/12 sec. 4a.5.3).

3. WHEN IT CAN JUMP TO THE UPPER CONE, AND WHERE PAULI FORBIDS IT.  An interband
   (Landau-Zener) transition needs the gap 2 v_F |k+A| to close, i.e. the state must
   pass the degeneracy at label k = -A. Whether that is allowed is pure Pauli:

       |A| < k_F : the point -A lies INSIDE the occupied disc. The conduction state
                   there is already filled -- blocked, no transitions, pure drift.
       |A| > k_F : the point -A has left the disc. The conduction state is empty and
                   pairs are created, in a strip of half-width k_perp = sqrt(E/pi v_F)
                   around the path (wiki/12 Eq. 8).

   So the SAME number, A_0 = k_F, both ends the linear regime and opens pair creation:
   E_sat = 27 kV/cm at E_F = 0.2 eV, 81 kV/cm at 0.6 eV.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import AU_EV, VF_EPM_AU  # noqa: E402

A0_PER_KVCM = 6.213e-4


def velocity_field(kx, ky, A, vf=1.0):
    """Conduction-band velocity of the label (kx, ky) at vector potential A (along x)."""
    px, py = kx + A, ky
    n = np.hypot(px, py)
    n = np.where(n > 0, n, 1e-30)
    return vf * px / n, vf * py / n


def disc_panel(ax, u, kF, lz_halfwidth=None, nring=3):
    """The occupied disc, the moving degeneracy point, and the velocity arrows."""
    A = u * kF
    th = np.linspace(0, 2 * np.pi, 400)
    # occupied disc in LABEL space -- it does not move
    ax.fill(kF * np.cos(th), kF * np.sin(th), color='#3498db', alpha=0.18, zorder=0)
    ax.plot(kF * np.cos(th), kF * np.sin(th), color='#2471a3', lw=1.4, zorder=1)
    # velocity arrows on a few rings inside it
    for frac in np.linspace(1.0, 0.35, nring):
        n = 12 if frac > 0.6 else 8
        a = np.linspace(0, 2 * np.pi, n, endpoint=False) + (0.13 if frac < 0.6 else 0.0)
        kx, ky = frac * kF * np.cos(a), frac * kF * np.sin(a)
        vx, vy = velocity_field(kx, ky, A)
        s = 0.32 * kF
        ax.quiver(kx, ky, vx, vy, angles='xy', scale_units='xy', scale=1.0 / s,
                  width=0.008, color='#c0392b', alpha=0.9, zorder=3)
    # the instantaneous degeneracy: label k = -A
    inside = abs(A) < kF
    ax.plot([-A], [0.0], marker='x', ms=11, mew=2.6,
            color='#7f8c8d' if inside else '#27ae60', zorder=5)
    ax.plot([-kF * 3, -A], [0, 0], ls=':', lw=1.1, color='#7f8c8d', zorder=2)
    if not inside and lz_halfwidth:
        ax.fill_between([-A - 0.1 * kF, -A + 0.1 * kF], -lz_halfwidth, lz_halfwidth,
                        color='#27ae60', alpha=0.25, zorder=2)
    ax.annotate('Dirac point,\n$k=-A$', xy=(-A, 0), xytext=(-A, -1.55 * kF),
                ha='center', va='top', fontsize=8,
                color='#7f8c8d' if inside else '#27ae60',
                arrowprops=dict(arrowstyle='->', lw=1.0,
                                color='#7f8c8d' if inside else '#27ae60'))
    ax.annotate('', xy=(0.85 * kF, 1.62 * kF), xytext=(-0.15 * kF, 1.62 * kF),
                arrowprops=dict(arrowstyle='-|>', lw=2.0, color='#8e44ad'))
    ax.text(0.35 * kF, 1.70 * kF, '$A(t)$', color='#8e44ad', fontsize=9, ha='center')
    ax.set_xlim(-2.1 * kF, 2.1 * kF); ax.set_ylim(-2.1 * kF, 2.0 * kF)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    return inside


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ef-ev', type=float, default=0.6)
    ap.add_argument('--drive', default=None, help='a DAST_*.txt to draw A(t)/k_F from')
    ap.add_argument('--out', default='cone_kinematics.png')
    args = ap.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available'); return 0

    kF = (abs(args.ef_ev) / AU_EV) / VF_EPM_AU
    esat = kF / A0_PER_KVCM
    fig = plt.figure(figsize=(13.4, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.85), hspace=0.55, wspace=0.24)

    # --- (a), (b): the disc, below and above saturation -----------------------------
    for j, (u, tag) in enumerate(((0.5, 'a'), (2.0, 'b'))):
        ax = fig.add_subplot(gs[0, j])
        # Landau-Zener strip half-width at the field that gives this u
        E_kv = u * esat
        Eau = E_kv * 1e5 / 5.14220675e11
        kperp = np.sqrt(Eau / (np.pi * VF_EPM_AU))
        inside = disc_panel(ax, u, kF, lz_halfwidth=kperp)
        ax.set_title(f'({tag})  $u=A_0/k_F={u:g}$   ($E_0={E_kv:.0f}$ kV/cm)', fontsize=10)
        if inside:
            ax.set_xlabel('degeneracy sits INSIDE the occupied disc\n'
                          'upper cone already filled there: PAULI BLOCKS IT\n'
                          'pure intraband drift, no pairs',
                          fontsize=8.5, color='#7f8c8d', labelpad=6)
        else:
            ax.set_xlabel('degeneracy has LEFT the disc\n'
                          'upper cone empty there: pairs are created\n'
                          'in a strip $|k_\\perp| \\sim \\sqrt{E/\\pi v_F}$',
                          fontsize=8.5, color='#27ae60', labelpad=6)

    # --- (c): the drift function G(u)/u ---------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    uu = np.logspace(-1.4, 1.1, 240)

    def g_of_u(u, n=40001):
        a, b = abs(1.0 - u), 1.0 + u
        rho = a + (np.arange(n) + 0.5) * (b - a) / n
        f = ((rho + u)**2 - 1.0) * (1.0 - (rho - u)**2)
        return np.sum(np.sqrt(np.maximum(f, 0.0))) * (b - a) / n / (np.pi * u)

    g = np.array([g_of_u(x) for x in uu])
    ax.loglog(uu, g, color='#c0392b', lw=2.0, label='$J/n e v_F = G(u)$')
    ax.loglog(uu, g / uu, color='#2980b9', lw=2.0, label='$\\sigma_{\\rm eff}/\\sigma_{\\rm lin}=G(u)/u$')
    ax.loglog(uu, np.minimum(uu, 1.0), ls=':', color='#7f8c8d', lw=1.2,
              label='linear drift / saturation')
    ax.axvline(1.0, ls='--', color='#27ae60', lw=1.2)
    ax.text(1.15, 0.55, 'pairs\nallowed', fontsize=8, color='#27ae60', va='center')
    ax.set_ylim(3e-2, 1.6)
    ax.set_xlabel('$u=A_0/k_F$'); ax.set_ylabel('current / conductivity, normalised')
    ax.set_title('(c)  the drift saturates: $|v|=v_F$ is a ceiling', fontsize=10)
    ax.legend(fontsize=7.5, loc='lower left'); ax.grid(alpha=0.25, which='both')

    # --- (d): energy cut along the field axis ---------------------------------------
    ax = fig.add_subplot(gs[1, :2])
    kx = np.linspace(-3 * kF, 3 * kF, 800)
    for u, c, ls in ((0.0, '#95a5a6', '--'), (0.5, '#2980b9', '-'), (2.0, '#27ae60', '-')):
        A = u * kF
        e = VF_EPM_AU * np.abs(kx + A) * AU_EV
        ax.plot(kx / kF, e, color=c, ls=ls, lw=1.8, label=f'$u={u:g}$')
        ax.plot(kx / kF, -e, color=c, ls=ls, lw=1.8)
    ef = args.ef_ev
    ax.axhline(ef, color='#c0392b', ls=':', lw=1.4)
    ax.text(2.95, ef, f'  $E_F$ = {ef:g} eV', color='#c0392b', fontsize=8.5, va='bottom', ha='right')
    ax.fill_between([-1, 1], 0, ef, color='#3498db', alpha=0.18)
    ax.text(0.0, ef * 0.45, 'occupied by the doping\n(upper cone, $|k| \\leq k_F$)',
            ha='center', fontsize=8.5, color='#2471a3')
    ax.plot([-0.5], [0], marker='x', ms=10, mew=2.4, color='#7f8c8d')
    ax.plot([-2.0], [0], marker='x', ms=10, mew=2.4, color='#27ae60')
    ax.annotate('Pauli-blocked', xy=(-0.5, 0), xytext=(0.55, -ef * 0.95), ha='left', fontsize=8.5,
                color='#7f8c8d', arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=1.0))
    ax.annotate('jump to the upper cone', xy=(-2.0, 0), xytext=(-2.9, -ef * 1.0), ha='left',
                fontsize=8.5, color='#27ae60',
                arrowprops=dict(arrowstyle='->', color='#27ae60', lw=1.0))
    ax.set_xlabel('label $k_x/k_F$ along the field'); ax.set_ylabel('energy [eV]')
    ax.set_title('(d)  $E_\\pm=\\pm v_F|k+A|$: the cone slides, the occupation does not',
                 fontsize=10)
    ax.set_ylim(-1.35 * ef, 1.55 * ef); ax.set_xlim(-3.2, 3.2)
    ax.legend(fontsize=8, loc='upper left', ncol=3)
    ax.grid(alpha=0.2)

    # --- (e): the drive on the same axis --------------------------------------------
    ax = fig.add_subplot(gs[1, 2])
    if args.drive and os.path.exists(args.drive):
        d = np.loadtxt(args.drive, comments='#')
        t, Ax = d[:, 0], d[:, 1]
        # the file is A [fs*V/Ang]; convert to a.u. via its own peak against A0_PER_KVCM
        peak_kv = np.max(np.abs(np.gradient(Ax, t))) * 1e5 / 1e5   # V/Ang/fs -> arbitrary
        Aau = Ax / np.max(np.abs(Ax)) if np.max(np.abs(Ax)) > 0 else Ax
        for E_kv, c in ((0.5 * esat, '#2980b9'), (2.0 * esat, '#27ae60')):
            ax.plot(t, Aau * E_kv * A0_PER_KVCM / kF, color=c, lw=1.5,
                    label=f'$E_0={E_kv:.0f}$ kV/cm')
        ax.set_xlabel('time [fs]')
    else:
        tt = np.linspace(-1.6, 1.6, 500)
        env = np.exp(-tt**2) * np.sin(np.pi * tt) / 0.62
        for u, c in ((0.5, '#2980b9'), (2.0, '#27ae60')):
            ax.plot(tt, u * env, color=c, lw=1.8, label=f'$A_0/k_F={u:g}$')
        ax.set_xlabel('time (single-cycle transient)')
    ax.axhspan(-1, 1, color='#3498db', alpha=0.15)
    ax.axhline(1, color='#2471a3', lw=1.0); ax.axhline(-1, color='#2471a3', lw=1.0)
    ax.text(0.98, 0.06, 'inside this band the degeneracy stays\nin the disc: Pauli-blocked',
            transform=ax.transAxes, fontsize=8, color='#2471a3', va='bottom', ha='right',
            bbox=dict(fc='white', ec='none', alpha=0.8, pad=1.5))
    ax.set_ylabel('$A(t)/k_F$')
    ax.set_title('(e)  when the excursion leaves the disc', fontsize=10)
    ax.legend(fontsize=8, loc='upper right'); ax.grid(alpha=0.25)

    fig.suptitle(f'A doped Dirac cone in a velocity-gauge field  '
                 f'($E_F$ = {ef:g} eV, $k_F$ = {kF:.4f} a.u., $A_0=k_F$ at {esat:.0f} kV/cm)',
                 fontsize=11.5)
    fig.subplots_adjust(top=0.90, bottom=0.075, left=0.055, right=0.985)
    fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')
    print(f'# k_F = {kF:.5f} a.u.,  E_sat = {esat:.1f} kV/cm,  '
          f'LZ strip half-width at E_sat = {np.sqrt(esat * 1e5 / 5.14220675e11 / (np.pi * VF_EPM_AU)):.5f} a.u.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
