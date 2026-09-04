#!/usr/bin/env python3
"""Drift saturation of a doped Dirac sheet: the universal curve behind the
non-monotonic transmission, and how well a given k-mesh represents it.

    python3 drift_saturation.py [GSDIR/graphene_sit --ef-ev 0.6] [more GS/EF pairs ...]
                                [--out drift_saturation.png]

Physics.  In the velocity gauge a uniform field displaces the whole occupied set by
the vector potential, k -> k + A, without moving the canonical labels.  For a Dirac
cone the band velocity has FIXED MODULUS v_F and only its direction turns, so the
sheet current of a displaced Fermi disc is

    J(A) = n e v_F * G(A/k_F),      G(u) = (1/pi u_disc) \\int_{|x|<=1} (x_par+u)/|x+u| d^2x,

which is linear at small u and saturates at n e v_F as u -> infinity: the carriers
cannot drift faster than v_F.  The differential (measured) conductivity is therefore

    sigma_eff(E_0)/sigma_lin = G'(A_0/k_F)      ->   k_F/A_0   for A_0 >> k_F,

i.e. the sheet becomes transparent as 1/E_0 once the excursion A_0 = |A|_max exceeds
the Fermi radius k_F = E_F/hbar v_F.  This script evaluates that ratio three ways:

  * the continuum integral above at the requested temperature (the physical answer);
  * the same sum on an actual SALMON k-mesh, with the eigenvalues, momentum matrix
    elements and Fermi-Dirac occupations the solver itself uses -- i.e. the adiabatic
    ground-state current of the truncated H_k(A) minus its undoped reference, exactly
    the quantity calc_current_bloch outputs;
  * optionally the current of real runs (--runs), which adds what the adiabatic
    picture leaves out: interband Landau-Zener pair creation at high field.

The mesh curve is the diagnostic.  A Fermi disc holding only one shell of mesh points
reproduces the linear limit to ~10 % but develops a spurious 10-20 % BUMP near
u ~ 0.2-0.5, which shows up in a field scan as a fake darkening before the physical
brightening; the bump shrinks as the disc is filled (wiki/12 sec. 4a.5).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sumrule_check import read_gs  # noqa: E402
from transmission import A_BOHR, AU_EV, VF_EPM_AU  # noqa: E402

KB_AU = 3.166811563e-6
AREA_BOHR2 = A_BOHR**2 * np.sqrt(3.0) / 2.0


def g_continuum(u, T_K, ef_au, vf=VF_EPM_AU, nr=700, nt=560):
    """J(A)/(A E_F/pi) of the continuum displaced Dirac sea at temperature T."""
    kT = KB_AU * max(T_K, 1.0)
    kF = ef_au / vf
    if u <= 0:
        return 1.0
    A = u * kF
    kmax = kF + max(8.0 * kT / vf, 3.0 * A) + 0.02
    r = (np.arange(nr) + 0.5) * kmax / nr
    th = (np.arange(nt) + 0.5) * 2.0 * np.pi / nt
    R, TH = np.meshgrid(r, th, indexing='ij')
    X = R * np.cos(TH); Y = R * np.sin(TH)
    fc = 1.0 / (1.0 + np.exp(np.clip((vf * R - ef_au) / kT, -60, 60)))
    den = np.hypot(X + A, Y)
    vx = np.where(den > 0, (X + A) / np.maximum(den, 1e-14), 0.0) * vf
    J = (4.0 / (4.0 * np.pi**2)) * np.sum(fc * vx * R) * (kmax / nr) * (2.0 * np.pi / nt)
    return J / (A * ef_au / np.pi)


def g_mesh(gs_prefix, ef_ev, T_K, us):
    """The same ratio evaluated on a real k-mesh: the adiabatic ground-state current of
    the truncated H_k(A) = eps + A p, doped occupation minus the undoped reference."""
    k, w, eig, p = read_gs(gs_prefix)
    nk = len(w)
    eD = 0.5 * (eig[:, 0].max() + eig[:, 1].min())
    ef = ef_ev / AU_EV
    kF = ef / VF_EPM_AU
    kT = KB_AU * max(T_K, 1.0)
    mu = eD + ef
    f = np.stack([2.0 / (1.0 + np.exp(np.clip((eig[:, b] - mu) / kT, -60, 60))) for b in range(2)], 1)
    f0 = np.stack([np.where(eig[:, b] < eD, 2.0, 0.0) for b in range(2)], 1)
    df = f - f0
    px = p[:, :, :, 0]
    eye = np.eye(2)
    npart = int(np.sum((f > 0.05) & (f < 1.95)))
    out = []
    for u in us:
        A = u * kF
        H = np.zeros((nk, 2, 2), dtype=complex)
        H[:, 0, 0] = eig[:, 0]; H[:, 1, 1] = eig[:, 1]
        H = H + A * px
        ev, U = np.linalg.eigh(H)
        occ = np.einsum('kin,kij,kjn->kn', U.conj(), px + A * eye, U).real
        J = (np.sum(w[:, None] * df * occ) / w.sum()) / AREA_BOHR2
        out.append(J / (A * ef / np.pi))
    return kF, np.array(out), npart, int(round(np.sqrt(nk)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gs', nargs='*', default=[], help='GS prefixes; pair each with --ef-ev')
    ap.add_argument('--ef-ev', nargs='*', type=float, default=[], help='E_F [eV] per GS prefix')
    ap.add_argument('--temp-k', type=float, default=300.0)
    ap.add_argument('--out', default='drift_saturation.png')
    args = ap.parse_args(argv)
    if len(args.ef_ev) != len(args.gs):
        args.ef_ev = list(args.ef_ev) + [0.2] * (len(args.gs) - len(args.ef_ev))

    us = np.concatenate([np.linspace(0.02, 1.0, 30), np.linspace(1.1, 8.0, 40)])
    ref_ef = args.ef_ev[0] / AU_EV if args.ef_ev else 0.2 / AU_EV
    cont = np.array([g_continuum(u, args.temp_k, ref_ef) for u in us])
    print(f'# continuum displaced Dirac sea, T = {args.temp_k:g} K: '
          f'sigma_eff/sigma_lin = 1 at u -> 0, {cont[np.argmin(abs(us - 1.0))]:.3f} at u = 1, '
          f'{cont[-1]:.3f} at u = {us[-1]:.0f}  (asymptote k_F/A_0)')
    curves = []
    for gsp, ef in zip(args.gs, args.ef_ev):
        kF, c, npart, nside = g_mesh(gsp, ef, args.temp_k, us)
        curves.append((f'{nside}$^2$, $E_F$ = {ef:g} eV ({npart} partial pts)', c))
        i04 = np.argmin(abs(us - 0.4))
        print(f'# {nside}^2, E_F = {ef:g} eV: k_F = {kF:.5f} a.u., {npart} partially occupied points; '
              f'ratio at u = 0.05 / 0.4 / 1 / 5: {c[0]:.3f} / {c[i04]:.3f} '
              f'/ {c[np.argmin(abs(us - 1.0))]:.3f} / {c[np.argmin(abs(us - 5.0))]:.3f}')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available -- numbers only'); return 0
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(us, cont, 'k-', lw=2.2, label='continuum displaced Dirac sea (exact)')
    ax.plot(us, 1.0 / np.maximum(us, 1e-9), 'k:', lw=1.0, label=r'$k_F/A_0$ asymptote')
    for (lab, c), col in zip(curves, ('#c0392b', '#2980b9', '#27ae60', '#8e44ad')):
        ax.plot(us, c, 'o-', ms=3, lw=1.2, color=col, label=lab)
    ax.axvline(1.0, ls='--', c='#7f8c8d', lw=1)
    ax.text(1.06, 0.06, '$A_0 = k_F$', color='#7f8c8d', fontsize=9)
    ax.set_xscale('log'); ax.set_xlabel(r'$A_0/k_F$  (vector-potential excursion / Fermi radius)')
    ax.set_ylabel(r'$\sigma_{\rm eff}/\sigma_{\rm linear}$')
    ax.set_title('Drift saturation on the Dirac cone, and how a k-mesh represents it', fontsize=10)
    ax.set_ylim(0, 1.35); ax.grid(alpha=0.25); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
