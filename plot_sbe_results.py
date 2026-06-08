#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotter for SALMON SBE real-time output files.

Drop this script into a calculation directory (next to SYSNAME_sbe_rt.data,
SYSNAME_sbe_rt_energy.data, SYSNAME_sbe_nex.data, SYSNAME_sbe_nex_k.data, ...)
and run it.  It scans the directory for these files, plots:

  - SYSNAME_sbe_rt_energy.data  : total energy vs time
  - SYSNAME_sbe_nex.data        : number of excited electrons/holes vs time
  - SYSNAME_sbe_nex_k.data      : per-k Houston-basis LCB population
        * snapshot PNG per saved time: three 2-D projections
          kx-ky, kx-kz, ky-kz (averaged over passive k-direction)
          with bicubic interpolation – one PNG per moment, time in filename
        * time-momentum map PNGs: one for each of kx/ky/kz showing
          population vs (k_axis, time) by marginalising the other two
          k-directions.  All time steps combined in a single image.

Memory strategy: SYSNAME_sbe_nex_k.data can be many GB.  The file is read
line-by-line; only one block's data is live at a time.  Snapshot plots are
saved immediately so the array can be reused.  Time-momentum maps accumulate
only three small 1-D marginals per time step (one number per unique k-value),
which is negligible regardless of nk.

No interactive windows are opened (Agg backend); everything is saved as PNG
into an output directory.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import argparse
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import zoom as ndimage_zoom

plt.switch_backend('Agg')

# ---------------------------------------------------------------------------
# Colourmap: turbo (perceptually-uniform rainbow, good for population maps)
# ---------------------------------------------------------------------------
CMAP_POP   = 'turbo'   # k-space snapshot and time-k maps
CMAP_CURVE = None      # line plots (default matplotlib colour cycle)


# ---------------------------------------------------------------------------
# Helpers shared by energy / nex plots
# ---------------------------------------------------------------------------

def parse_header(header_line):
    """Extract column names from a numbered header, ignoring units in []."""
    return re.findall(r'\d+:([^\[\s]+)(?:\[[^\]]*\])?', header_line)


def find_header(filepath):
    """Return the first comment line that carries column-index 1 (e.g. '# 1:Time...')."""
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if re.match(r'#\s*1\s*:\s*\S', line):
                return line.strip(), i
    raise ValueError(f"Numbered header line not found in {filepath}")


def load_columns(filepath):
    """Load a whitespace-separated SBE .data file into (column_names, 2-D array)."""
    header_line, _ = find_header(filepath)
    column_names = parse_header(header_line)
    rows = []
    with open(filepath, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                continue
    data = np.array(rows) if rows else np.empty((0, len(column_names)))
    return column_names, data


def plot_xy(time, values, time_name, col_name, output_path, dpi=150):
    if len(time) == 0:
        print(f"  (skip) no data for {col_name}")
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time, values, linewidth=1.0)
    ax.set_xlabel(time_name)
    ax.set_ylabel(col_name)
    ax.set_title(f'{col_name} vs {time_name}')
    ax.grid(True, alpha=0.3, linestyle='--')
    fig.tight_layout()
    safe_col  = re.sub(r'[^\w\-]', '_', col_name)
    safe_time = re.sub(r'[^\w\-]', '_', time_name)
    out_file = output_path / f'{safe_col}_vs_{safe_time}.png'
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out_file.name}")


def plot_energy_and_nex(filepath, output_dir, dpi=150):
    print(f"Processing {filepath.name} ...")
    cols, data = load_columns(filepath)
    if data.size == 0:
        print("  (skip) no data")
        return
    time_name, time = cols[0], data[:, 0]
    for j in range(1, len(cols)):
        plot_xy(time, data[:, j], time_name, cols[j], output_dir, dpi=dpi)


# ---------------------------------------------------------------------------
# Streaming nex_k parser
# ---------------------------------------------------------------------------

def _iter_nex_k_blocks(filepath):
    """
    Generator: yield (t_value, t_unit, kpoints[nk,3], pop[nk]) one block at a
    time while reading *filepath* line-by-line.  Only the current block is
    kept in memory at once.
    """
    time_re = re.compile(r'#\s*t\s*=\s*([-+\d.eEdD]+)\s*(\S*)')
    t_value = None
    t_unit  = ''
    kx, ky, kz, pop = [], [], [], []

    def _flush():
        if t_value is not None and kx:
            yield (t_value, t_unit,
                   np.column_stack([kx, ky, kz]),
                   np.asarray(pop, dtype=float))

    with open(filepath, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            m = time_re.match(s)
            if m:
                # flush previous block
                if t_value is not None and kx:
                    yield (t_value, t_unit,
                           np.column_stack([kx, ky, kz]),
                           np.asarray(pop, dtype=float))
                t_value = float(m.group(1))
                t_unit  = m.group(2)
                kx, ky, kz, pop = [], [], [], []
                continue
            if s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 5:
                continue
            try:
                kx.append(float(parts[1]))
                ky.append(float(parts[2]))
                kz.append(float(parts[3]))
                pop.append(float(parts[4]))
            except ValueError:
                continue

    # last block
    if t_value is not None and kx:
        yield (t_value, t_unit,
               np.column_stack([kx, ky, kz]),
               np.asarray(pop, dtype=float))


# ---------------------------------------------------------------------------
# K-grid helpers
# ---------------------------------------------------------------------------

def _build_grid_info(kpoints):
    """
    From a flat list of k-points (nk,3) discover the regular 3-D grid
    structure and return index arrays ix, iy, iz so that pop[ik] belongs
    to cell (ix[ik], iy[ik], iz[ik]).
    """
    rtol = 1e-9
    kx_u = np.unique(np.round(kpoints[:, 0], 9))
    ky_u = np.unique(np.round(kpoints[:, 1], 9))
    kz_u = np.unique(np.round(kpoints[:, 2], 9))
    ix = np.searchsorted(kx_u, np.round(kpoints[:, 0], 9))
    iy = np.searchsorted(ky_u, np.round(kpoints[:, 1], 9))
    iz = np.searchsorted(kz_u, np.round(kpoints[:, 2], 9))
    return kx_u, ky_u, kz_u, ix, iy, iz


def _fill_3d(pop, ix, iy, iz, nx, ny, nz):
    """Scatter flat pop array into a 3-D grid; cells with no data stay NaN."""
    g = np.full((nx, ny, nz), np.nan)
    g[ix, iy, iz] = pop
    return g


def _project(g3d, axis):
    """Average over *axis* ignoring NaN; returns 2-D array."""
    return np.nanmean(g3d, axis=axis)


def _interp2d(grid2d, k_a, k_b, factor=8):
    """
    Bicubic interpolation of a 2-D population grid onto a finer mesh.
    Returns (ka_fine, kb_fine, grid_fine).
    """
    na, nb = grid2d.shape
    if na < 2 or nb < 2:
        return k_a, k_b, grid2d

    # Replace NaN with nearest-neighbour fill before interpolation
    filled = np.where(np.isnan(grid2d), 0.0, grid2d)

    interp = RegularGridInterpolator(
        (k_a, k_b), filled, method='linear', bounds_error=False, fill_value=None)

    ka_fine = np.linspace(k_a[0], k_a[-1], na * factor)
    kb_fine = np.linspace(k_b[0], k_b[-1], nb * factor)
    KA, KB  = np.meshgrid(ka_fine, kb_fine, indexing='ij')
    grid_fine = interp((KA, KB))
    return ka_fine, kb_fine, grid_fine


def _heatmap_ax(ax, k_a, k_b, grid2d, label_a, label_b, title,
                vmin=None, vmax=None, factor=8):
    """
    Draw a projected + bicubic-interpolated heatmap on *ax*.
    k_a → horizontal axis, k_b → vertical axis.
    """
    if grid2d.size == 0 or np.all(np.isnan(grid2d)):
        ax.set_title(title + " (no data)")
        return None

    ka_f, kb_f, gf = _interp2d(grid2d, k_a, k_b, factor=factor)

    # imshow: rows = kb (vertical), cols = ka (horizontal), origin='lower'
    # extent: [ka_min, ka_max, kb_min, kb_max]
    im = ax.imshow(
        gf.T,                       # .T so rows=kb, cols=ka
        origin='lower',
        aspect='auto',
        extent=[ka_f[0], ka_f[-1], kb_f[0], kb_f[-1]],
        cmap=CMAP_POP,
        vmin=vmin, vmax=vmax,
        interpolation='nearest',
    )
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel(f'{label_a} [a.u.]')
    ax.set_ylabel(f'{label_b} [a.u.]')
    ax.set_title(title)
    return im


# ---------------------------------------------------------------------------
# Snapshot: 3 projected planes for one time step
# ---------------------------------------------------------------------------

def _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi):
    """
    Save one PNG with three heatmap panels (projected kx-ky, kx-kz, ky-kz).
    pop3d has shape (nx, ny, nz) with axes in (kx, ky, kz) order.
    """
    vmin = np.nanmin(pop3d)
    vmax = max(np.nanmax(pop3d), vmin + 1e-30)   # avoid zero-range cmap

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    # kx-ky: average over kz (axis 2)
    _heatmap_ax(axes[0], kx_u, ky_u, _project(pop3d, 2),
                'kx', 'ky', 'pop_lcb: kx-ky (avg kz)',
                vmin=vmin, vmax=vmax)

    # kx-kz: average over ky (axis 1)
    _heatmap_ax(axes[1], kx_u, kz_u, _project(pop3d, 1),
                'kx', 'kz', 'pop_lcb: kx-kz (avg ky)',
                vmin=vmin, vmax=vmax)

    # ky-kz: average over kx (axis 0)
    _heatmap_ax(axes[2], ky_u, kz_u, _project(pop3d, 0),
                'ky', 'kz', 'pop_lcb: ky-kz (avg kx)',
                vmin=vmin, vmax=vmax)

    fig.suptitle(f'Houston-basis LCB population,  t = {t_val:.6f} {t_unit}')
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    safe_t = f'{t_val:.6f}'.replace('-', 'm').replace('+', 'p')
    out_file = output_dir / f'nex_k_snap_t{safe_t}{t_unit}.png'
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out_file.name}")


# ---------------------------------------------------------------------------
# Time–k maps (one per axis direction)
# ---------------------------------------------------------------------------

def _save_kt_map(times, t_unit, k_vals, label_k, marginals, output_dir, dpi):
    """
    2-D colour map: horizontal axis = time, vertical axis = k_axis.
    *marginals* is a list of 1-D arrays (one per time step), each of length
    len(k_vals) = population averaged over the other two k-directions.
    """
    if not marginals:
        return
    mat = np.array(marginals).T           # shape (nk_1d, nt)
    nt  = len(times)
    nk  = len(k_vals)

    fig, ax = plt.subplots(figsize=(max(8, nt * 0.08 + 2), 5))

    # Use pcolormesh so each cell is exactly one (time, k) bin
    t_edges = _bin_edges(np.asarray(times))
    k_edges = _bin_edges(k_vals)
    im = ax.pcolormesh(t_edges, k_edges, mat, cmap=CMAP_POP, shading='flat')
    plt.colorbar(im, ax=ax, label='population_lcb (avg)')
    ax.set_xlabel(f'time [{t_unit}]')
    ax.set_ylabel(f'{label_k} [a.u.]')
    ax.set_title(f'LCB population vs time and {label_k}')
    fig.tight_layout()

    out_file = output_dir / f'nex_k_ktmap_{label_k}.png'
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out_file.name}")


def _bin_edges(centers):
    """Convert bin centres to edges (one extra element)."""
    c = np.asarray(centers, dtype=float)
    if len(c) == 1:
        delta = 1.0
    else:
        delta = np.diff(c)
        delta = np.append(delta, delta[-1])
    edges = np.empty(len(c) + 1)
    edges[0]  = c[0]  - 0.5 * (delta[0] if len(c) > 1 else delta)
    edges[1:] = c + 0.5 * delta
    return edges


# ---------------------------------------------------------------------------
# Main nex_k driver (streaming)
# ---------------------------------------------------------------------------

def plot_nex_k(filepath, output_dir, dpi=150):
    print(f"Processing {filepath.name} ...")

    # Grid structure (filled from first block)
    kx_u = ky_u = kz_u = None
    ix = iy = iz = None
    nx = ny = nz = 0
    pop3d = None   # reused each block

    # Marginal accumulators (small, kept for all time steps)
    times = []
    t_unit_last = ''
    marg_kx, marg_ky, marg_kz = [], [], []

    n_blocks = 0
    for t_val, t_unit, kpoints, pop in _iter_nex_k_blocks(filepath):
        t_unit_last = t_unit
        n_blocks   += 1

        if kx_u is None:
            kx_u, ky_u, kz_u, ix, iy, iz = _build_grid_info(kpoints)
            nx, ny, nz = len(kx_u), len(ky_u), len(kz_u)
            pop3d = np.empty((nx, ny, nz))

        pop3d.fill(np.nan)
        pop3d[ix, iy, iz] = pop

        # Snapshot plot – save immediately, reuse pop3d next iteration
        _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi)

        # Accumulate marginals for kt-maps
        times.append(t_val)
        marg_kx.append(np.nanmean(pop3d, axis=(1, 2)))
        marg_ky.append(np.nanmean(pop3d, axis=(0, 2)))
        marg_kz.append(np.nanmean(pop3d, axis=(0, 1)))

    if n_blocks == 0:
        print("  (skip) no data blocks found")
        return

    # Time–k maps
    print(f"  writing time-k maps ({n_blocks} time steps) ...")
    _save_kt_map(times, t_unit_last, kx_u, 'kx', marg_kx, output_dir, dpi)
    _save_kt_map(times, t_unit_last, ky_u, 'ky', marg_ky, output_dir, dpi)
    _save_kt_map(times, t_unit_last, kz_u, 'kz', marg_kz, output_dir, dpi)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Plot SALMON SBE real-time output files.')
    parser.add_argument('-i', '--input-dir', default='.',
                        help='Directory with SYSNAME_sbe_*.data files')
    parser.add_argument('-o', '--output', default='sbe_plots',
                        help='Output directory for PNGs')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Image resolution')
    args = parser.parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    energy_files = sorted(input_dir.glob('*_sbe_rt_energy.data'))
    nex_files    = sorted(input_dir.glob('*_sbe_nex.data'))
    nex_k_files  = sorted(input_dir.glob('*_sbe_nex_k.data'))

    for f in energy_files:
        plot_energy_and_nex(f, output_dir, dpi=args.dpi)
    for f in nex_files:
        plot_energy_and_nex(f, output_dir, dpi=args.dpi)
    for f in nex_k_files:
        plot_nex_k(f, output_dir, dpi=args.dpi)

    if not (energy_files or nex_files or nex_k_files):
        print(f"No *_sbe_rt_energy.data / *_sbe_nex.data / *_sbe_nex_k.data "
              f"files found in {input_dir.resolve()}")
        return

    print(f"\nDone.  Output: {output_dir.resolve()}")


if __name__ == '__main__':
    main()
