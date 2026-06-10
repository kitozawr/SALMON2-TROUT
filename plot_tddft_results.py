#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified plotter for SALMON TDDFT real-time output files.

Drop into a calculation directory and run:

    python3 plot_tddft_results.py                       # auto-detect all
    python3 plot_tddft_results.py --bands 16 17 18 19   # explicit band list
    python3 plot_tddft_results.py --downsample 10       # thin large rt.data

What is plotted
---------------
  *_rt.data          : E-field / vector potential + current density / dipole
  *_rt_energy.data   : total energy (Eall and Eall-Eall0) vs time
  *_nex.data         : excited electrons and holes vs time
  *_ovlp.data        : per-band, per-k occupation projections:
                         band population k-averaged vs time (selected bands)
                         time × k-marginal maps per selected band
                         [optional --snapshots] per-time k-space 2D projections

Note: SALMON writes _ovlp.data (not _ovpl.data).  The file is always
auto-detected by the correct suffix.

k-coordinates
-------------
  If *_k.data is present (reduced coords, written by both DFT and EPM),
  k-projections and time-k maps are produced.  Without it only the
  band-population-vs-time line plot is written.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import argparse

plt.switch_backend('Agg')

# ---------------------------------------------------------------------------
# Hardcoded defaults (override at runtime with CLI flags)
# ---------------------------------------------------------------------------
CMAP_POP       = 'turbo'
CMAP_LOG_SCALE = False   # --log-cmap
SNAP_ENABLED   = False   # --snapshots

# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _make_norm(vmin, vmax, log_scale):
    if log_scale and vmax > 0:
        floor = max(vmax * 1e-6, 1e-30)
        return mcolors.LogNorm(vmin=max(vmin, floor), vmax=vmax)
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def _detect_time_unit(filepath):
    with open(filepath) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            m = re.search(r'[Tt]ime\[(\w+)\]', line)
            if m:
                return m.group(1)
    return 'a.u.'


def _parse_col_names(filepath):
    """Return {col_name: 0-based_index} from '# N:name[unit]' header lines."""
    cols = {}
    with open(filepath) as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            for m in re.finditer(r'(\d+):([\w\-]+)\[', line):
                cols[m.group(2)] = int(m.group(1)) - 1
    return cols


def _load_columns(filepath, downsample=1):
    """Load numeric data (skip # lines) into an (N, ncol) ndarray."""
    rows = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                pass
    if not rows:
        return np.empty((0, 0))
    data = np.array(rows)
    return data[::downsample] if downsample > 1 else data


def _detect_sysname(calc_dir):
    for suffix in ['_rt.data', '_nex.data', '_rt_energy.data', '_ovlp.data']:
        for p in sorted(Path(calc_dir).glob(f'*{suffix}')):
            return p.name[: -len(suffix)]
    return None


def _find_file(calc_dir, sysname, suffix):
    p = calc_dir / f'{sysname}{suffix}'
    return p if p.exists() else None


def _load_k_coords(k_filepath):
    """Return (nk, 3) array of reduced k-coords from *_k.data."""
    pts = []
    with open(k_filepath) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    pts.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except ValueError:
                    pass
    return np.array(pts) if pts else None


# ---------------------------------------------------------------------------
# *_rt.data — external/total fields + current density / dipole moment
# ---------------------------------------------------------------------------

def plot_rt(filepath, output_dir, dpi=150, downsample=1):
    filepath = Path(filepath)
    print(f'# plot_rt: {filepath.name}')
    cols = _parse_col_names(filepath)
    data = _load_columns(filepath, downsample)
    if data.size == 0:
        return
    t      = data[:, 0]
    t_unit = _detect_time_unit(filepath)

    def _get(name):
        i = cols.get(name)
        return data[:, i] if (i is not None and i < data.shape[1]) else None

    def _first(*names):
        for n in names:
            v = _get(n)
            if v is not None:
                return v, n
        return None, None

    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)

    # Panel 1: external field
    ax = axes[0]
    for comp, c in [('x', 'C0'), ('y', 'C1'), ('z', 'C2')]:
        v, lbl = _first(f'E_ext_{comp}', f'Ac_ext_{comp}')
        if v is not None:
            ax.plot(t, v, color=c, lw=0.8, label=comp)
    ylab = 'E_ext' if cols.get('E_ext_x') is not None else 'Ac_ext'
    ax.set_ylabel(ylab); ax.legend(fontsize=7, ncol=3); ax.grid(True, lw=0.3)

    # Panel 2: matter current (3D) or dipole moment (0D)
    ax = axes[1]
    if cols.get('Jm_x') is not None:
        for comp, c in [('x', 'C0'), ('y', 'C1'), ('z', 'C2')]:
            v = _get(f'Jm_{comp}')
            if v is not None:
                ax.plot(t, v, color=c, lw=0.8, label=comp)
        ax.set_ylabel('J_mat')
    elif cols.get('dm_x') is not None:
        for comp, c in [('x', 'C0'), ('y', 'C1'), ('z', 'C2')]:
            v = _get(f'dm_{comp}')
            if v is not None:
                ax.plot(t, v, color=c, lw=0.8, label=comp)
        ax.set_ylabel('dipole')
    elif cols.get('ddm_e_x') is not None:
        for comp, c in [('x', 'C0'), ('y', 'C1'), ('z', 'C2')]:
            v = _get(f'ddm_e_{comp}')
            if v is not None:
                ax.plot(t, v, color=c, lw=0.8, label=comp)
        ax.set_ylabel('d(dipole)/dt')
    ax.legend(fontsize=7, ncol=3); ax.grid(True, lw=0.3)

    # Panel 3: total field
    ax = axes[2]
    for comp, c in [('x', 'C0'), ('y', 'C1'), ('z', 'C2')]:
        v, _ = _first(f'E_tot_{comp}', f'Ac_tot_{comp}')
        if v is not None:
            ax.plot(t, v, color=c, lw=0.8, label=comp)
    ylab = 'E_tot' if cols.get('E_tot_x') is not None else 'Ac_tot'
    ax.set_ylabel(ylab); ax.set_xlabel(f'time [{t_unit}]')
    ax.legend(fontsize=7, ncol=3); ax.grid(True, lw=0.3)

    fig.suptitle(filepath.stem, fontsize=10)
    fig.tight_layout()
    out = Path(output_dir) / f'{filepath.stem}_fields.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'#   -> {out.name}')


# ---------------------------------------------------------------------------
# *_rt_energy.data
# ---------------------------------------------------------------------------

def plot_rt_energy(filepath, output_dir, dpi=150, downsample=1):
    filepath = Path(filepath)
    print(f'# plot_rt_energy: {filepath.name}')
    data = _load_columns(filepath, downsample)
    if data.size == 0:
        return
    t      = data[:, 0]
    t_unit = _detect_time_unit(filepath)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    # Column layout: 1:Time 2:Eall 3:Eall-Eall0  (use col 2, index 2)
    if data.shape[1] >= 3:
        ax.plot(t, data[:, 2], 'C0', lw=0.9, label='Eall - Eall0')
    if data.shape[1] >= 2:
        ax.plot(t, data[:, 1], 'C1--', lw=0.6, alpha=0.6, label='Eall')
    ax.set_xlabel(f'time [{t_unit}]'); ax.set_ylabel('Energy')
    ax.legend(fontsize=8); ax.grid(True, lw=0.3)
    fig.suptitle(filepath.stem, fontsize=10)
    fig.tight_layout()
    out = Path(output_dir) / f'{filepath.stem}_energy.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'#   -> {out.name}')


# ---------------------------------------------------------------------------
# *_nex.data
# ---------------------------------------------------------------------------

def plot_nex(filepath, output_dir, dpi=150, downsample=1):
    filepath = Path(filepath)
    print(f'# plot_nex: {filepath.name}')
    data = _load_columns(filepath, downsample)
    if data.size == 0:
        return
    t      = data[:, 0]
    t_unit = _detect_time_unit(filepath)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.plot(t, data[:, 1], 'C0', lw=0.9, label='N_exc')
    if data.shape[1] >= 3:
        ax.plot(t, data[:, 2], 'C1--', lw=0.9, label='N_hole')
    ax.set_xlabel(f'time [{t_unit}]'); ax.set_ylabel('count')
    ax.legend(fontsize=8); ax.grid(True, lw=0.3)
    fig.suptitle(filepath.stem, fontsize=10)
    fig.tight_layout()
    out = Path(output_dir) / f'{filepath.stem}_nex.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'#   -> {out.name}')


# ---------------------------------------------------------------------------
# *_ovlp.data streaming parser
# ---------------------------------------------------------------------------

def _iter_ovlp_blocks(filepath):
    """
    Yields (it: int, occ: ndarray shape (nk, nb)) for each time block.

    File layout after header lines:
        <it>                          <- integer on its own line
        ik  occ_1  occ_2 ... occ_NB  <- nk lines (ik is 1-based, discarded)
        <it>
        ...
    Values are occupation numbers divided by k-weight (SALMON convention:
    coef/wtk), so for equal-weight grids each value ≈ nk × fractional_occ.
    For spin-unpolarised calculations the range is 0..2*nk/nk = 0..2.
    """
    with open(filepath) as fh:
        # skip header comment lines
        first_data = ''
        for line in fh:
            if not line.startswith('#'):
                first_data = line
                break
        if not first_data.strip():
            return
        current_it = int(first_data.strip())
        rows = []
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) == 1:
                if rows:
                    yield current_it, np.array(rows, dtype=float)
                    rows = []
                current_it = int(parts[0])
            else:
                # first token = ik (discard), rest = occupations
                rows.append([float(x) for x in parts[1:]])
        if rows:
            yield current_it, np.array(rows, dtype=float)


def _guess_nocc(occ_block):
    """Estimate number of occupied bands from the first (ground-state) block."""
    mean_occ = occ_block.mean(axis=0)
    thresh   = 0.5 * mean_occ.max()
    return int(np.sum(mean_occ >= thresh))


def _k_marginal(occ_1d, k_col, k_unique):
    """Mean occupation for each unique k-coordinate value."""
    out = np.empty(len(k_unique))
    for i, kv in enumerate(k_unique):
        mask = np.abs(k_col - kv) < 1e-9
        out[i] = float(occ_1d[mask].mean()) if mask.any() else 0.0
    return out


# ---------------------------------------------------------------------------
# *_ovlp.data  — top-level plot driver
# ---------------------------------------------------------------------------

def plot_ovlp(filepath, nex_filepath, output_dir, dpi=150,
              bands_to_plot=None, log_scale=False, snapshots=False,
              k_filepath=None):
    """
    Three products from *_ovlp.data:

    1. *_ovlp_bands.png       — k-averaged occupation vs time (selected bands)
    2. *_ovlp_kt_band<N>.png  — time × k-marginal heat maps per selected band
    3. snap_band<N>_*.png     — [if snapshots] 2D k-projections per time step
    """
    filepath = Path(filepath)
    print(f'# plot_ovlp: {filepath.name}  '
          f'(log_cmap={log_scale}, snapshots={snapshots})')

    # Physical time axis from nex.data (written in the same subroutine call)
    phys_times = None
    t_unit = 'step'
    if nex_filepath and Path(nex_filepath).exists():
        nex_data = _load_columns(nex_filepath)
        if nex_data.size:
            phys_times = nex_data[:, 0]
            t_unit = _detect_time_unit(nex_filepath)

    # k-coordinates (optional — enables spatial projections)
    kpoints = _load_k_coords(k_filepath) if (k_filepath and Path(k_filepath).exists()) else None
    k_uniq  = [np.unique(kpoints[:, d]) for d in range(3)] if kpoints is not None else None

    # --- single streaming pass ---
    nb = nk = nocc = None
    band_pop_time = []   # (nt,) of (nb,) arrays  — k-averaged occupation
    it_indices    = []
    # time-k marginal accumulators  {band: list of 1-D marginal arrays}
    marg = {d: {} for d in range(3)}  # marg[dim][band] = list of (nk_d,) arrays

    snap_idx = 0
    for it, occ in _iter_ovlp_blocks(filepath):
        if nb is None:
            nk, nb = occ.shape
            nocc   = _guess_nocc(occ)
            if bands_to_plot is None:
                # default: last 2 valence + first 4 conduction bands
                bands_to_plot = list(range(max(1, nocc - 1), min(nb + 1, nocc + 5)))
            for d in range(3):
                for b in bands_to_plot:
                    marg[d][b] = []

        band_pop_time.append(occ.mean(axis=0))
        it_indices.append(it)

        if kpoints is not None:
            for b in bands_to_plot:
                ib0 = b - 1
                if ib0 < nb:
                    pop = occ[:, ib0]
                    for d in range(3):
                        marg[d][b].append(_k_marginal(pop, kpoints[:, d], k_uniq[d]))

        if snapshots and kpoints is not None:
            _save_ovlp_snap(occ, kpoints, k_uniq, it, snap_idx, phys_times,
                            t_unit, bands_to_plot, nb, output_dir, dpi, log_scale)
        snap_idx += 1

    if nb is None:
        print('#   ovlp.data empty, skipping.')
        return

    band_pop = np.array(band_pop_time)   # (nt, nb)
    nt       = band_pop.shape[0]

    if phys_times is not None:
        n_use = min(nt, len(phys_times))
        t_arr = phys_times[:n_use]
        band_pop = band_pop[:n_use]
    else:
        t_arr = np.array(it_indices, dtype=float)
        n_use = nt

    # 1 — band populations vs time
    _plot_band_vs_time(band_pop, t_arr, t_unit, bands_to_plot, nocc, nb,
                       filepath.stem, output_dir, dpi)

    # 2 — time-k maps
    if kpoints is not None:
        for b in bands_to_plot:
            if not marg[0].get(b):
                continue
            mx = np.array(marg[0][b][:n_use]).T   # (nkx, nt)
            my = np.array(marg[1][b][:n_use]).T
            mz = np.array(marg[2][b][:n_use]).T
            _plot_kt_map(mx, my, mz, k_uniq, t_arr, t_unit,
                         b, filepath.stem, output_dir, dpi, log_scale)


# ---------------------------------------------------------------------------
# ovlp sub-plots
# ---------------------------------------------------------------------------

def _plot_band_vs_time(band_pop, t_arr, t_unit, bands, nocc, nb,
                       stem, output_dir, dpi):
    nt = min(len(t_arr), band_pop.shape[0])
    fig, ax = plt.subplots(figsize=(9, 4))
    for b in bands:
        ib0 = b - 1
        if not (0 <= ib0 < nb):
            continue
        tag = ''
        if b == nocc:
            tag = ' (VBM)'
        elif b == nocc + 1:
            tag = ' (LCB)'
        ax.plot(t_arr[:nt], band_pop[:nt, ib0], lw=0.9, label=f'band {b}{tag}')
    ax.set_xlabel(f'time [{t_unit}]')
    ax.set_ylabel('mean occupation (k-avg)')
    ax.legend(fontsize=7, ncol=min(4, len(bands)))
    ax.grid(True, lw=0.3)
    fig.suptitle(f'{stem}  —  band populations', fontsize=10)
    fig.tight_layout()
    out = Path(output_dir) / f'{stem}_ovlp_bands.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'#   -> {out.name}')


def _plot_kt_map(mx, my, mz, k_uniq, t_arr, t_unit,
                 band, stem, output_dir, dpi, log_scale):
    nt = min(len(t_arr), mx.shape[1], my.shape[1], mz.shape[1])
    vmax = float(np.nanmax([mx[:, :nt], my[:, :nt], mz[:, :nt]]))
    norm = _make_norm(0.0, vmax, log_scale)

    dirs = [
        ('kx', k_uniq[0], mx[:, :nt]),
        ('ky', k_uniq[1], my[:, :nt]),
        ('kz', k_uniq[2], mz[:, :nt]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, (lbl, k_ax, mat) in zip(axes, dirs):
        ext = [t_arr[0], t_arr[nt - 1], k_ax[0], k_ax[-1]]
        im  = ax.imshow(mat, origin='lower', aspect='auto', extent=ext,
                        cmap=CMAP_POP, norm=norm, interpolation='bilinear')
        ax.set_xlabel(f'time [{t_unit}]', fontsize=8)
        ax.set_ylabel(f'{lbl} [reduced]', fontsize=8)
        plt.colorbar(im, ax=ax, pad=0.02, label='occ')
    fig.suptitle(f'{stem}  band {band}  time-k map', fontsize=9)
    fig.tight_layout()
    out = Path(output_dir) / f'{stem}_ovlp_kt_band{band:03d}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'#   -> {out.name}')


def _save_ovlp_snap(occ, kpoints, k_uniq, it, snap_idx, phys_times,
                    t_unit, bands, nb, output_dir, dpi, log_scale):
    """2-D k-projection snapshots for each selected band at one time step."""
    if phys_times is not None and snap_idx < len(phys_times):
        t_label = f't{phys_times[snap_idx]:.4g}{t_unit}'
    else:
        t_label = f'it{it:06d}'

    plane_specs = [
        ('kx', 'ky', 0, 1),
        ('kx', 'kz', 0, 2),
        ('ky', 'kz', 1, 2),
    ]
    for b in bands:
        ib0 = b - 1
        if ib0 >= nb:
            continue
        pop  = occ[:, ib0]
        vmax = float(np.nanmax(pop))
        norm = _make_norm(0.0, vmax, log_scale)

        fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
        for ax, (xl, yl, d1, d2) in zip(axes, plane_specs):
            ku1, ku2 = k_uniq[d1], k_uniq[d2]
            grid = np.zeros((len(ku1), len(ku2)))
            cnt  = np.zeros_like(grid)
            for ik in range(len(pop)):
                i1 = int(np.argmin(np.abs(ku1 - kpoints[ik, d1])))
                i2 = int(np.argmin(np.abs(ku2 - kpoints[ik, d2])))
                grid[i1, i2] += pop[ik]
                cnt[i1, i2]  += 1
            with np.errstate(invalid='ignore'):
                grid = np.where(cnt > 0, grid / cnt, 0.0)
            im = ax.imshow(grid.T, origin='lower', aspect='auto',
                           extent=[ku1[0], ku1[-1], ku2[0], ku2[-1]],
                           cmap=CMAP_POP, norm=norm, interpolation='bicubic')
            ax.set_xlabel(f'{xl} [red]', fontsize=8)
            ax.set_ylabel(f'{yl} [red]', fontsize=8)
            plt.colorbar(im, ax=ax, pad=0.02)
        fig.suptitle(f'band {b}  {t_label}', fontsize=9)
        fig.tight_layout()
        # sanitise t_label for use in filename
        safe_lbl = re.sub(r'[^\w\.\-]', '_', t_label)
        out = Path(output_dir) / f'snap_band{b:03d}_{safe_lbl}.png'
        fig.savefig(out, dpi=dpi, bbox_inches='tight')
        plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Plot SALMON TDDFT real-time output files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('directory', nargs='?', default='.',
                        help='Calculation directory')
    parser.add_argument('--sysname',    help='Override auto-detected SYSNAME')
    parser.add_argument('--output-dir', default='tddft_plots',
                        help='Output directory for PNG files')
    parser.add_argument('--dpi',        type=int, default=150)
    parser.add_argument('--downsample', type=int, default=1,
                        help='Keep every N-th line in large rt / energy files')
    parser.add_argument('--bands',      type=int, nargs='+',
                        help='Band indices (1-based) for ovlp plots; '
                             'default: last 2 valence + first 4 conduction')
    parser.add_argument('--log-cmap',  action='store_true', default=CMAP_LOG_SCALE,
                        help=f'Logarithmic colormap for 2D heat maps '
                             f'(hardcoded default: CMAP_LOG_SCALE={CMAP_LOG_SCALE})')
    parser.add_argument('--snapshots', action='store_true', default=SNAP_ENABLED,
                        help=f'Write per-time-step k-space snapshots '
                             f'(hardcoded default: SNAP_ENABLED={SNAP_ENABLED})')
    args = parser.parse_args()

    calc_dir = Path(args.directory).resolve()
    sysname  = args.sysname or _detect_sysname(calc_dir)
    if sysname is None:
        print('ERROR: cannot auto-detect SYSNAME.  Use --sysname.')
        return

    print(f'# SYSNAME     = {sysname}')
    print(f'# output dir  = {args.output_dir}')
    out_dir = calc_dir / args.output_dir
    out_dir.mkdir(exist_ok=True)

    def f(suf):
        return _find_file(calc_dir, sysname, suf)

    rt  = f('_rt.data')
    rte = f('_rt_energy.data')
    nex = f('_nex.data')
    ovlp = f('_ovlp.data')
    kfile = f('_k.data')

    if rt:
        plot_rt(rt, out_dir, args.dpi, args.downsample)
    if rte:
        plot_rt_energy(rte, out_dir, args.dpi, args.downsample)
    if nex:
        plot_nex(nex, out_dir, args.dpi, args.downsample)
    if ovlp:
        plot_ovlp(ovlp, nex, out_dir,
                  dpi=args.dpi,
                  bands_to_plot=args.bands,
                  log_scale=args.log_cmap,
                  snapshots=args.snapshots,
                  k_filepath=kfile)

    print('# Done.')


if __name__ == '__main__':
    main()
