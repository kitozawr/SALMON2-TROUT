#!/usr/bin/env python3
"""Initial level occupation of a graphene sheet: what `sbe_ef_ev` / `sbe_temp_init_k`
actually put on the k-mesh, side by side with the undoped filling.

    python3 plot_occupation.py GSDIR/graphene_sit [--ef-ev 0.2] [--temp-init-k 300]
                               [--out occupation.png] [--kmax 0.28]

Reads the ground-state files the SBE reads (`*_k.data`, `*_eigen.data`,
`*_tm.data`) and applies exactly the rule of `gs_info_ssbe`:
f_n(k) = occ_max * f_FD(eps_n(k); mu, T_init), mu = E_D + sbe_ef_ev, so the picture
is the solver's initial density matrix, not a sketch. For each occupation it prints
and draws:

  * the two Dirac bands against |k - K| (obtained from the local gap, 2 v_F |k-K|),
    coloured by occupation -- the Fermi level and k_F are marked;
  * the sheet density n_2D = (sum_nk f - 2)/(N_k A_2D), i.e. the ADDED charge, which
    is the number to quote for a doped run;
  * the count of partially occupied k-points, the mesh criterion of wiki/12 sec. 4a.2
    (below ~20 the Fermi surface is not resolved and the Drude weight is meaningless);
  * the radial profile of the conduction occupation, which shows directly how many
    mesh shells sit inside the Fermi circle.

Use it before a doped production run: if the middle panel shows a Fermi disc with
only one shell in it, raise num_kgrid or E_F (wiki/12 sec. 4a.0, Eq. 4a.3).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sumrule_check import read_gs  # noqa: E402
from transmission import AU_EV, VF_EPM_AU, A_BOHR  # noqa: E402

KB_AU = 3.166811563e-6
BOHR_CM = 0.52917721067e-8
AREA_BOHR2 = A_BOHR**2 * np.sqrt(3.0) / 2.0


def occupation(eig, mu, T, occ_max=2.0):
    if T <= 0.0:
        return np.where(eig < mu - 1e-9, occ_max, np.where(eig > mu + 1e-9, 0.0, 0.5 * occ_max))
    return occ_max / (1.0 + np.exp(np.clip((eig - mu) / (KB_AU * T), -60.0, 60.0)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gs_prefix', help='GSDIR/SYSNAME of the ground state')
    ap.add_argument('--ef-ev', type=float, default=0.2, help='doped Fermi level from the Dirac point [eV]')
    ap.add_argument('--temp-init-k', type=float, default=300.0)
    ap.add_argument('--kmax', type=float, default=0.28, help='|k-K| range of the band panels [a.u.]')
    ap.add_argument('--out', default='occupation.png')
    args = ap.parse_args(argv)

    k, w, eig, p = read_gs(args.gs_prefix)
    nk, nb = eig.shape
    if nb < 2:
        print('# need at least 2 bands'); return 1
    eD = 0.5 * (eig[:, 0].max() + eig[:, 1].min())          # Dirac point (e-h symmetric)
    dk = np.abs(eig[:, 1] - eig[:, 0]) / (2.0 * VF_EPM_AU)  # |k - K| from the local gap
    cases = [('intrinsic (undoped)', 0.0, 0.0), (f'doped, $E_F$ = {args.ef_ev:g} eV', args.ef_ev, args.temp_init_k)]

    stats = []
    for lab, ef, T in cases:
        mu = eD + ef / AU_EV
        f = np.column_stack([occupation(eig[:, b], mu, T) for b in range(2)])
        n2d = (np.sum(f) / nk - 2.0) / (AREA_BOHR2 * BOHR_CM**2)
        npart = int(np.sum((f > 0.05) & (f < 1.95)))
        kF = (abs(ef) / AU_EV) / VF_EPM_AU
        stats.append((lab, ef, T, f, n2d, npart, kF))
        print(f'# {lab}: added charge {np.sum(f) / nk - 2.0:+.5e} e/cell -> n_2D = {n2d:+.3e} cm^-2;'
              f' k_F = {kF:.5f} a.u.; partially occupied k-points = {npart}'
              + ('  <-- Fermi surface UNDER-RESOLVED' if 0 < npart < 20 else ''))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available -- numbers only'); return 0

    E = [(eig[:, b] - eD) * AU_EV for b in range(2)]
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.1))
    for j, (lab, ef, T, f, n2d, npart, kF) in enumerate(stats):
        sc = ax[j].scatter(np.r_[dk, dk], np.r_[E[0], E[1]], c=np.r_[f[:, 0], f[:, 1]],
                           s=7, cmap='viridis', vmin=0, vmax=2)
        ax[j].axhline(ef, ls='--', c='#c0392b', lw=1.2)
        ax[j].text(args.kmax * 0.78, ef + 0.06, '$E_F$', color='#c0392b', fontsize=9)
        if kF > 0:
            ax[j].axvline(kF, ls=':', c='#c0392b', lw=1)
            ax[j].text(kF * 1.06, -1.45, '$k_F$', color='#c0392b', fontsize=8)
        ax[j].set_title(f'{lab}\n$n_{{2D}}$ = {n2d:.2e} cm$^{{-2}}$, partially occupied k-points: {npart}',
                        fontsize=9)
        ax[j].set_xlabel('$|k-K|$ [a.u.]'); ax[j].set_xlim(0, args.kmax); ax[j].set_ylim(-1.8, 1.8)
        ax[j].grid(alpha=0.2)
    ax[0].set_ylabel('$E-E_D$ [eV]')
    cb = fig.colorbar(sc, ax=ax[1]); cb.set_label('occupation per cell (0..2)', fontsize=8)

    o = np.argsort(dk)
    for (lab, ef, T, f, n2d, npart, kF), c in zip(stats, ('#2c3e50', '#c0392b')):
        ax[2].plot(dk[o], f[o, 1], '.', ms=3, color=c, label=f'{lab.split(",")[0]}: conduction')
        if kF > 0:
            ax[2].axvline(kF, ls=':', c=c, lw=1)
            ax[2].text(kF * 1.06, 1.2, '$k_F$', color=c, fontsize=8)
    n_side = int(round(np.sqrt(nk)))
    ax[2].set_xlim(0, max(4.0 * max(s[6] for s in stats), 0.02))
    ax[2].set_xlabel('$|k-K|$ [a.u.]'); ax[2].set_ylabel('conduction occupation')
    ax[2].set_title(f'Fermi surface on the {n_side}$^2$ mesh', fontsize=9)
    ax[2].legend(fontsize=8); ax[2].grid(alpha=0.2)
    fig.tight_layout(); fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
