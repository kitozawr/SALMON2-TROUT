#!/usr/bin/env python3
"""Level-population dynamics of a graphene sheet run: what the k-resolved
diabatic (fixed-basis) populations do while the THz field acts.

    python3 plot_levels.py runs/E100kVcm_mem/graphene_sit [--times 100,150,200,300] [--band 2]
                          [--out plot_levels.png] [--kwin 0.35]

Reads, per prefix P:
    P_sbe_nex_k_lev_real.data   snapshots (one block per "# t = ... fs" line) of the
                                fixed-basis population of every band at every k
                                (per spin, 0..1; VBM = band 1); written every
                                out_projection_k_step steps (make_inputs: --snap-fs)
    P_sbe_nex.data              continuous Houston conduction population (per cell)
    P_sbe_nex_nonad.data        col 3: dressed-reference (ring-visible) density
    P_k.data                    k-points (for the K-point markers)
Prints, for every snapshot, the BZ-averaged population of each band (electrons per
cell per spin) and the electron/hole numbers; with matplotlib it draws (a) the band
populations vs time on top of the continuous Houston trace, (b) k-space maps of the
chosen band's population around the Dirac points at the requested times (the LZ
pair creation shows up as population on the k-points nearest K along the field).
"""
import argparse
import os
import re
import sys

import numpy as np

AU_T_FS = 0.02418884326505
A_BOHR = 4.648726
B_LEN = 4.0 * np.pi / (np.sqrt(3.0) * A_BOHR)         # |b1| = |b2|


def read_snapshots(path):
    """-> times [fs], k (nk,3), pops (nt, nk, nb)"""
    times, blocks, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                m = re.search(r'#\s*t\s*=\s*([-+0-9.eE]+)\s*fs', line)
                if m:
                    if cur:
                        blocks.append(np.array(cur, float)); cur = []
                    times.append(float(m.group(1)))
                continue
            if line.strip():
                cur.append(line.split())
    if cur:
        blocks.append(np.array(cur, float))
    n = min(len(times), len(blocks))
    data = np.array(blocks[:n])                       # (nt, nk, 4 + nb)
    return np.array(times[:n]), data[0, :, 1:4], data[:, :, 4:]


def k_points_special():
    """K and K' of the hexagonal cell used by x14 (a1 = a x, a2 = a(1/2, sqrt3/2)) in Cartesian a.u."""
    b1 = 2 * np.pi / A_BOHR * np.array([1.0, -1.0 / np.sqrt(3.0)])
    b2 = 2 * np.pi / A_BOHR * np.array([0.0, 2.0 / np.sqrt(3.0)])
    K = (b1 + 2 * b2) / 3.0
    Kp = (2 * b1 + b2) / 3.0
    return K, Kp


def load_xy(path):
    d = np.loadtxt(path, comments='#')
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') and re.search(r'1:time\[', line, re.I):
                unit = re.search(r'1:time\[([^\]]*)\]', line, re.I).group(1)
                t = d[:, 0] if unit.startswith('fs') else d[:, 0] * AU_T_FS
                return t, d[:, 1:]
    return d[:, 0] * AU_T_FS, d[:, 1:]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('prefix')
    ap.add_argument('--times', default='', help='snapshot times [fs] to map (nearest available); default: 4 spread over the run')
    ap.add_argument('--band', type=int, default=0, help='band index for the k-maps (default: first conduction band = VBM+1)')
    ap.add_argument('--nval', type=int, default=1, help='number of filled bands in the file (VBM = band nval)')
    ap.add_argument('--kwin', type=float, default=0.35, help='half-width of the k-map window around K [a.u.]')
    ap.add_argument('--out', default='plot_levels.png')
    args = ap.parse_args(argv)
    P = args.prefix
    snap = P + '_sbe_nex_k_lev_real.data'
    if not os.path.exists(snap):
        print(f'no {snap}: set out_projection_k_step (make_inputs.py --snap-fs) to write snapshots'); return 1
    times, k, pops = read_snapshots(snap)
    nt, nk, nb = pops.shape
    band = args.band if args.band > 0 else args.nval + 1
    # file normalisation: a filled band reads 2 (spin-summed) or 1 (per spin) at t = 0
    full = float(np.round(pops[0, :, :args.nval].max())) or 1.0
    spin = 2.0 / full
    print(f'# {snap}: {nt} snapshots, nk = {nk}, nb = {nb}, VBM = band {args.nval}, map band = {band}, filled band = {full:g}')
    print('#   t[fs]   ' + ' '.join(f'<pop_b{b + 1}>' for b in range(nb)) + '   n_e/cell(spin-summed)  n_h/cell')
    ne_t, nh_t = [], []
    for it in range(nt):
        avg = pops[it].mean(axis=0)
        ne = spin * avg[args.nval:].sum(); nh = spin * (full * args.nval - avg[:args.nval].sum())
        ne_t.append(ne); nh_t.append(nh)
        print(f'  {times[it]:8.2f}   ' + ' '.join(f'{a:9.6f}' for a in avg) + f'   {ne:.4e}  {nh:.4e}')
    if abs(ne_t[-1] - nh_t[-1]) > 1e-3 * max(ne_t[-1], 1e-30):
        print('#   NOTE: n_e /= n_h at the end -> population left in bands outside the written set (or dissipators active)')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available -- table only'); return 0
    want = [float(x) for x in args.times.split(',') if x.strip()] or list(np.linspace(times[0], times[-1], 5)[1:])
    idx = sorted(set(int(np.argmin(np.abs(times - tw))) for tw in want))
    fig = plt.figure(figsize=(4.2 * (1 + len(idx)), 4.0))
    ax0 = fig.add_subplot(1, 1 + len(idx), 1)
    for b in range(nb):
        y = spin * pops[:, :, b].mean(axis=1) if b >= args.nval else spin * (full - pops[:, :, b].mean(axis=1))
        ax0.plot(times, y, 'o-', ms=3, label=f'band {b + 1} ({"e" if b >= args.nval else "h"}) snapshots')
    p = P + '_sbe_nex.data'
    if os.path.exists(p):
        t, y = load_xy(p); ax0.plot(t, y[:, 0], 'k-', lw=0.8, label='Houston n_c (continuous)')
    p = P + '_sbe_nex_nonad.data'
    if os.path.exists(p):
        t, y = load_xy(p); ax0.plot(t, y[:, 1], 'k--', lw=0.8, label='dressed-reference n_c')
    ax0.set_yscale('log'); ax0.set_xlabel('t [fs]'); ax0.set_ylabel('carriers per cell'); ax0.legend(fontsize=6)
    ax0.set_title(os.path.basename(os.path.dirname(P)) or P, fontsize=9)
    K, Kp = k_points_special()
    vmax = max(pops[idx][:, :, band - 1].max(), 1e-12)
    for j, it in enumerate(idx):
        ax = fig.add_subplot(1, 1 + len(idx), 2 + j)
        # fold k into the window around K (use the shortest image over the reciprocal lattice)
        b1 = 2 * np.pi / A_BOHR * np.array([1.0, -1.0 / np.sqrt(3.0)]); b2 = 2 * np.pi / A_BOHR * np.array([0.0, 2.0 / np.sqrt(3.0)])
        kk = k[:, :2] - K
        best = kk.copy(); bestd = (kk**2).sum(axis=1)
        for i1 in (-1, 0, 1):
            for i2 in (-1, 0, 1):
                sh = kk + i1 * b1 + i2 * b2; dd = (sh**2).sum(axis=1); m = dd < bestd
                best[m] = sh[m]; bestd[m] = dd[m]
        m = np.abs(best[:, 0]) < args.kwin
        m &= np.abs(best[:, 1]) < args.kwin
        sc = ax.scatter(best[m, 0], best[m, 1], c=pops[it, m, band - 1], s=18, cmap='magma', vmin=0, vmax=vmax, marker='h')
        ax.plot(0, 0, 'c+', ms=10); ax.set_aspect('equal')
        ax.set_title(f'band {band} pop., t = {times[it]:.0f} fs', fontsize=9); ax.set_xlabel('k_x - K_x [a.u.]')
        if j == 0:
            ax.set_ylabel('k_y - K_y [a.u.]')
    fig.colorbar(sc, ax=fig.axes[1:], shrink=0.8, label=f'population (filled band = {full:g})')
    fig.savefig(args.out, dpi=140, bbox_inches='tight')
    print(f'# wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
