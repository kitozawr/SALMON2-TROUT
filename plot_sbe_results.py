#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified plotter for SALMON-SBE real-time output and EPM/DFT ground-state files.

Drop into a calculation directory and run:

    python3 plot_sbe_results.py                     # auto-detect everything
    python3 plot_sbe_results.py --downsample 200    # thin out large RT curves
    python3 plot_sbe_results.py --band-path L Gamma X W K

What is plotted
---------------
  *_sbe_rt.data          : fields + current vs time  (downsampled if requested)
  *_sbe_rt_energy.data   : total energy vs time
  *_sbe_nex.data         : excited electron count vs time
  *_sbe_nex_k.data       : per-k Houston-basis LCB population:
                             snapshot PNGs (3 projected planes) + time-k maps
  *_k.data + *_eigen.data: band structure along the requested path
                             k in reduced coords, energy shifted to VBM = 0 eV

Spinor (spin-orbit split) input files
-------------------------------------
  Datasets generated with spin-orbit coupling (epm_gaas_reference.py with
  INCLUDE_SPIN_ORBIT = True, consumed by SBE with yn_sbe_spinor = 'y') carry
  2*Nb spinor bands with occupation 1 per band instead of 2.  The plotter
  detects this automatically from the occupation column of *_eigen.data
  (max occupation <= 1) and then treats each pair of adjacent (Kramers
  partner) bands as ONE level: occupations of the two spin sub-bands are
  summed (1+1 = 2 per valence level) and the band plot draws one curve per
  level (energy = mean of the spin pair) on top of the faint spin-resolved
  sub-bands, so tiny Dresselhaus splittings don't render as doubled lines
  while the real spin-orbit splittings (e.g. Gamma8/Gamma7) stay visible.
  Control with --spin-sum {auto,on,off}.

Memory strategy
---------------
  Large *_sbe_rt*.data files are read with downsampling (--downsample N keeps
  every N-th data line).
  *_sbe_nex_k.data is read line-by-line; only one time block lives in RAM.
  Band eigenvalues are read fully (typically < 10 MB).

k-vector units
--------------
  After the EPM Python fix (PR #34), *_k.data contains dimensionless reduced
  coordinates kx/(2π/a), BZ ∈ [-0.5, 0.5].  nex_k writes gs%kpoint verbatim,
  so nex_k k-values are also in reduced units.  All k-axis labels say "reduced".
"""

import re
import itertools
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from scipy.interpolate import RegularGridInterpolator
import matplotlib.colors as mcolors

plt.switch_backend('Agg')

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HA_TO_EV = 27.211386245988     # Hartree → eV
CMAP_POP = 'turbo'             # population heat maps

# ---- Hardcoded default for 2-D colormap scaling ----
# False = linear (Normalize),  True = logarithmic (LogNorm).
# Override at runtime with the --log-cmap CLI flag.
CMAP_LOG_SCALE = False

# ---- Hardcoded default for per-time-step k-space snapshots ----
# False = only time-k maps are written (automatic, lightweight).
# True  = additionally write one PNG per time block (can be many files).
# Override at runtime with the --snapshots CLI flag.
SNAP_ENABLED = False

# FCC high-symmetry points in reduced coordinates of the conventional cubic BZ.
# BZ spans [-0.5, 0.5].  Labels correspond to folded FCC points.
HS_POINTS = {
    'Gamma': [ 0.000,  0.000,  0.000],
    'L':     [ 0.500,  0.500,  0.500],
    'X':     [ 0.000,  0.500,  0.500],
    'W':     [ 0.500,  0.250, -0.250],   # [0.5,0.25,0.75] wrapped to BZ
    'K':     [ 0.375,  0.375, -0.250],   # [0.375,0.375,0.75] wrapped
    'U':     [-0.375,  0.250, -0.375],   # [0.625,0.25,0.625] wrapped
}
DEFAULT_BAND_PATH = ['L', 'Gamma', 'X', 'W', 'K']


# ===========================================================================
# Shared column-file helpers
# ===========================================================================

def parse_header(header_line):
    """Extract column names from a numbered header, ignoring units in []."""
    return re.findall(r'\d+:([^\[\s]+)(?:\[[^\]]*\])?', header_line)


def find_header(filepath):
    """Return first comment line that carries column index 1 (e.g. '# 1:Time...')."""
    with open(filepath, 'r') as f:
        for line in f:
            if re.match(r'#\s*1\s*:\s*\S', line):
                return line.strip()
    raise ValueError(f"Numbered header line not found in {filepath}")


def load_columns_streaming(filepath, downsample=1):
    """
    Read a whitespace-separated SBE .data file.
    With downsample > 1 only every N-th data line is kept (cheap for GB files).
    Returns (column_names, 2-D numpy array).
    """
    header_line = find_header(filepath)
    column_names = parse_header(header_line)
    rows = []
    data_count = 0
    with open(filepath, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            data_count += 1
            if data_count % downsample != 0:
                continue
            try:
                rows.append([float(x) for x in s.split()])
            except ValueError:
                continue
    data = np.array(rows) if rows else np.empty((0, len(column_names)))
    return column_names, data


# ===========================================================================
# RT line plots  (*_sbe_rt.data, *_sbe_rt_energy.data, *_sbe_nex.data)
# ===========================================================================

def _plot_xy(time, values, time_name, col_name, output_path, dpi=150):
    if len(time) == 0:
        print(f"  (skip) no data for {col_name}")
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(time, values, linewidth=0.8)
    ax.set_xlabel(time_name)
    ax.set_ylabel(col_name)
    ax.set_title(f'{col_name}  vs  {time_name}')
    ax.grid(True, alpha=0.3, linestyle='--')
    fig.tight_layout()
    safe_col  = re.sub(r'[^\w\-]', '_', col_name)
    safe_time = re.sub(r'[^\w\-]', '_', time_name)
    out = output_path / f'{safe_col}_vs_{safe_time}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def plot_rt_file(filepath, output_dir, downsample=1, dpi=150):
    """Plot all non-time columns in a columnar SBE .data file."""
    print(f"Processing {filepath.name}  (downsample={downsample}) ...")
    cols, data = load_columns_streaming(filepath, downsample=downsample)
    if data.size == 0:
        print("  (skip) no data")
        return
    n_kept = data.shape[0]
    if downsample > 1:
        print(f"  {n_kept:,} lines kept after downsampling")
    time_name, time = cols[0], data[:, 0]
    for j in range(1, min(len(cols), data.shape[1])):
        _plot_xy(time, data[:, j], time_name, cols[j], output_dir, dpi=dpi)


# ===========================================================================
# Streaming nex_k  (*_sbe_nex_k.data)
# ===========================================================================

def _iter_nex_k_blocks(filepath):
    """Yield (t_val, t_unit, kpoints[nk,3], pop[nk]) one block at a time."""
    time_re = re.compile(r'#\s*t\s*=\s*([-+\d.eEdD]+)\s*(\S*)')
    t_value, t_unit = None, ''
    kx, ky, kz, pop = [], [], [], []

    with open(filepath, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            m = time_re.match(s)
            if m:
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
                kx.append(float(parts[1])); ky.append(float(parts[2]))
                kz.append(float(parts[3])); pop.append(float(parts[4]))
            except ValueError:
                continue

    if t_value is not None and kx:
        yield (t_value, t_unit,
               np.column_stack([kx, ky, kz]),
               np.asarray(pop, dtype=float))


def _build_grid_info(kpoints):
    kx_u = np.unique(np.round(kpoints[:, 0], 9))
    ky_u = np.unique(np.round(kpoints[:, 1], 9))
    kz_u = np.unique(np.round(kpoints[:, 2], 9))
    ix = np.searchsorted(kx_u, np.round(kpoints[:, 0], 9))
    iy = np.searchsorted(ky_u, np.round(kpoints[:, 1], 9))
    iz = np.searchsorted(kz_u, np.round(kpoints[:, 2], 9))
    return kx_u, ky_u, kz_u, ix, iy, iz


def _project(g3d, axis):
    return np.nanmean(g3d, axis=axis)


def _interp2d(grid2d, k_a, k_b, factor=8):
    na, nb_ = grid2d.shape
    if na < 2 or nb_ < 2:
        return k_a, k_b, grid2d
    filled = np.where(np.isnan(grid2d), 0.0, grid2d)
    interp = RegularGridInterpolator(
        (k_a, k_b), filled, method='linear', bounds_error=False, fill_value=None)
    ka_f = np.linspace(k_a[0], k_a[-1], na  * factor)
    kb_f = np.linspace(k_b[0], k_b[-1], nb_ * factor)
    KA, KB = np.meshgrid(ka_f, kb_f, indexing='ij')
    return ka_f, kb_f, interp((KA, KB))


def _make_norm(vmin, vmax, log_scale):
    """Return a matplotlib Normalize or LogNorm for colormap scaling."""
    if log_scale and vmax > 0:
        # Floor at 1e-6 × peak so zeros don't break the log scale
        floor = max(vmax * 1e-6, 1e-30)
        return mcolors.LogNorm(vmin=max(vmin, floor), vmax=vmax)
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def _heatmap_ax(ax, k_a, k_b, grid2d, label_a, label_b, title,
                vmin=None, vmax=None, factor=8, log_scale=False):
    if grid2d.size == 0 or np.all(np.isnan(grid2d)):
        ax.set_title(title + " (no data)")
        return None
    ka_f, kb_f, gf = _interp2d(grid2d, k_a, k_b, factor=factor)
    norm = _make_norm(vmin if vmin is not None else np.nanmin(gf),
                      vmax if vmax is not None else np.nanmax(gf),
                      log_scale)
    im = ax.imshow(
        gf.T, origin='lower', aspect='auto',
        extent=[ka_f[0], ka_f[-1], kb_f[0], kb_f[-1]],
        cmap=CMAP_POP, norm=norm, interpolation='nearest')
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel(f'{label_a} [reduced]')
    ax.set_ylabel(f'{label_b} [reduced]')
    ax.set_title(title)
    return im


def _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi,
                   log_scale=False):
    vmin = np.nanmin(pop3d)
    vmax = max(np.nanmax(pop3d), vmin + 1e-30)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))
    _heatmap_ax(axes[0], kx_u, ky_u, _project(pop3d, 2),
                'kx', 'ky', 'pop_lcb: kx-ky (avg kz)',
                vmin=vmin, vmax=vmax, log_scale=log_scale)
    _heatmap_ax(axes[1], kx_u, kz_u, _project(pop3d, 1),
                'kx', 'kz', 'pop_lcb: kx-kz (avg ky)',
                vmin=vmin, vmax=vmax, log_scale=log_scale)
    _heatmap_ax(axes[2], ky_u, kz_u, _project(pop3d, 0),
                'ky', 'kz', 'pop_lcb: ky-kz (avg kx)',
                vmin=vmin, vmax=vmax, log_scale=log_scale)

    fig.suptitle(f'Houston-basis LCB population,  t = {t_val:.6f} {t_unit}')
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    safe_t  = f'{t_val:.6f}'.replace('-', 'm').replace('+', 'p')
    out = output_dir / f'nex_k_snap_t{safe_t}{t_unit}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def _bin_edges(centers):
    c = np.asarray(centers, dtype=float)
    delta = np.append(np.diff(c), np.diff(c)[-1]) if len(c) > 1 else np.array([1.0])
    edges = np.empty(len(c) + 1)
    edges[0]  = c[0]  - 0.5 * delta[0]
    edges[1:] = c + 0.5 * delta
    return edges


def _save_kt_map(times, t_unit, k_vals, label_k, marginals, output_dir, dpi,
                 log_scale=False):
    if not marginals:
        return
    mat = np.array(marginals).T           # (nk_1d, nt)
    vmin_m, vmax_m = float(np.nanmin(mat)), float(np.nanmax(mat))
    norm = _make_norm(vmin_m, max(vmax_m, vmin_m + 1e-30), log_scale)
    fig, ax = plt.subplots(figsize=(max(8, len(times) * 0.08 + 2), 5))
    im = ax.pcolormesh(_bin_edges(np.asarray(times)), _bin_edges(k_vals),
                       mat, cmap=CMAP_POP, norm=norm, shading='flat')
    plt.colorbar(im, ax=ax, label='population_lcb (avg)')
    ax.set_xlabel(f'time [{t_unit}]')
    ax.set_ylabel(f'{label_k} [reduced]')
    ax.set_title(f'LCB population vs time and {label_k}')
    fig.tight_layout()
    out = output_dir / f'nex_k_ktmap_{label_k}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def plot_nex_k(filepath, output_dir, dpi=150, log_scale=False, snapshots=False):
    print(f"Processing {filepath.name}  "
          f"(cmap={'log' if log_scale else 'linear'}, "
          f"snapshots={'on' if snapshots else 'off'}) ...")
    kx_u = ky_u = kz_u = ix = iy = iz = None
    pop3d = None
    times, marg_kx, marg_ky, marg_kz = [], [], [], []
    t_unit_last = ''
    n_blocks = 0

    for t_val, t_unit, kpoints, pop in _iter_nex_k_blocks(filepath):
        t_unit_last = t_unit
        n_blocks   += 1
        if kx_u is None:
            kx_u, ky_u, kz_u, ix, iy, iz = _build_grid_info(kpoints)
            pop3d = np.empty((len(kx_u), len(ky_u), len(kz_u)))
        pop3d.fill(np.nan)
        pop3d[ix, iy, iz] = pop
        if snapshots:
            _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi,
                           log_scale=log_scale)
        times.append(t_val)
        marg_kx.append(np.nanmean(pop3d, axis=(1, 2)))
        marg_ky.append(np.nanmean(pop3d, axis=(0, 2)))
        marg_kz.append(np.nanmean(pop3d, axis=(0, 1)))

    if n_blocks == 0:
        print("  (skip) no data blocks found")
        return

    print(f"  writing time-k maps ({n_blocks} time steps) ...")
    _save_kt_map(times, t_unit_last, kx_u, 'kx', marg_kx, output_dir, dpi,
                 log_scale=log_scale)
    _save_kt_map(times, t_unit_last, ky_u, 'ky', marg_ky, output_dir, dpi,
                 log_scale=log_scale)
    _save_kt_map(times, t_unit_last, kz_u, 'kz', marg_kz, output_dir, dpi,
                 log_scale=log_scale)


# ===========================================================================
# Band structure  (*_k.data + *_eigen.data)
# ===========================================================================

def _load_kpoints(kfile):
    """
    Parse SYSNAME_k.data.
    k-vectors are in reduced (dimensionless) coordinates: kx/(2π/a) ∈ [-0.5, 0.5].
    Returns ndarray (nk, 3).
    """
    pts = []
    with open(kfile, 'r') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) >= 4:
                try:
                    pts.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except ValueError:
                    continue
    if not pts:
        raise ValueError(f"No k-points found in {kfile}")
    return np.array(pts)


def _load_eigenvalues(eigenfile, nk):
    """
    Parse SYSNAME_eigen.data (EPM/SALMON format).
    Block headers are comment lines '# ik = N'.
    Data lines: 'ib  energy_Ha  occup'.
    Returns (eigen[nb, nk] in Ha, occup[nb, nk], vbm_ha).
    """
    ik_re = re.compile(r'#\s*ik\s*=\s*(\d+)')
    eigen_map = {}     # 1-based ik → list of (energy_Ha, occup)
    current_k = None
    vbm = -np.inf

    with open(eigenfile, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                m = ik_re.match(s)
                if m:
                    current_k = int(m.group(1))
                    eigen_map.setdefault(current_k, [])
                continue
            if current_k is None:
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                e, occ = float(parts[1]), float(parts[2])
                eigen_map[current_k].append((e, occ))
                if occ > 0.1 and e > vbm:
                    vbm = e
            except ValueError:
                continue

    if not eigen_map:
        raise ValueError(f"No eigenvalue blocks found in {eigenfile}")

    nb = max(len(v) for v in eigen_map.values())
    eigen = np.full((nb, nk), np.nan)
    occup = np.full((nb, nk), np.nan)
    for ik, entries in eigen_map.items():
        if 1 <= ik <= nk:
            for ib, (e, occ) in enumerate(entries[:nb]):
                eigen[ib, ik - 1] = e
                occup[ib, ik - 1] = occ

    return eigen, occup, vbm


def _detect_spinor(occup):
    """
    True if the occupation column looks like a spinor (spin-orbit split)
    dataset: occupied bands carry 1 electron each instead of 2
    (epm_gaas_reference.py with INCLUDE_SPIN_ORBIT, SBE with yn_sbe_spinor='y').
    """
    occ_max = np.nanmax(occup)
    return bool(0.0 < occ_max <= 1.0 + 1e-6)


def _sym_equivalents(point):
    """
    Generate all distinct cubic symmetry equivalents of *point*
    (all permutations × all sign flips), each wrapped to (-0.5, 0.5].
    """
    p = np.asarray(point, dtype=float)
    seen, result = set(), []
    for perm in itertools.permutations(p):
        for signs in itertools.product([1, -1], repeat=3):
            c = np.array(perm) * np.array(signs)
            c = c - np.floor(c + 0.5)          # wrap to (-0.5, 0.5]
            key = tuple(np.round(c, 8))
            if key not in seen:
                seen.add(key)
                result.append(c)
    return result


def _nearest_k(ideal, k_db, snap_tol):
    """
    Return index (1-based) of the grid k-point nearest to *ideal*,
    accounting for all cubic symmetry equivalents of *ideal*.
    Returns (k_idx, dist) or (None, dist) if dist > snap_tol.
    """
    equivs = _sym_equivalents(ideal)
    best_id, best_dist = None, np.inf
    for kp in k_db:
        d = min(np.linalg.norm(kp['c'] - eq) for eq in equivs)
        if d < best_dist:
            best_dist, best_id = d, kp['id']
    return (best_id if best_dist <= snap_tol else None), best_dist


def plot_band_structure(kfile, eigenfile, output_dir,
                        path_labels=None, hs_points=None,
                        energy_range_ev=(-6, 12), dpi=150, spin_sum='auto'):
    """
    Plot band structure from *_k.data and *_eigen.data.
    k-points must be in reduced (dimensionless) coordinates.
    Eigenvalues are read in Hartree, shifted to VBM = 0, plotted in eV.

    spin_sum: 'auto' — detect a spinor (spin-orbit split) dataset from the
              occupation column (max occ <= 1) and merge adjacent (Kramers
              partner) spin sub-bands into levels: occupations summed,
              level energy = mean of the pair; the spin-resolved sub-bands
              are kept as faint lines underneath.
              'on'  — force the merge, 'off' — plot all bands as-is.
    """
    if path_labels is None:
        path_labels = DEFAULT_BAND_PATH
    if hs_points is None:
        hs_points = HS_POINTS

    # Validate path labels
    unknown = [l for l in path_labels if l not in hs_points]
    if unknown:
        raise ValueError(f"Unknown high-symmetry labels: {unknown}. "
                         f"Available: {list(hs_points)}")

    print(f"Processing band structure: {kfile.name} + {eigenfile.name}")

    # --- k-points -------------------------------------------------------
    kpts = _load_kpoints(kfile)
    nk   = len(kpts)
    k_db = [{'id': i + 1, 'c': kpts[i]} for i in range(nk)]
    print(f"  {nk} k-points")

    # Estimate grid spacing → snap tolerance (use full grid spacing)
    kx_uniq = np.unique(np.round(kpts[:, 0], 8))
    grid_sp  = float(np.min(np.diff(kx_uniq))) if len(kx_uniq) > 1 else 0.5
    snap_tol = grid_sp            # snap within one grid step

    # --- eigenvalues ----------------------------------------------------
    eigen_ha, occup, vbm_ha = _load_eigenvalues(eigenfile, nk)
    nb = eigen_ha.shape[0]
    print(f"  {nb} bands, VBM = {vbm_ha:.6f} Ha = {vbm_ha * HA_TO_EV:.4f} eV")

    # Spinor (spin-orbit split) input: merge spin pairs into levels?
    if spin_sum == 'on':
        spinor_merge = True
    elif spin_sum == 'off':
        spinor_merge = False
    else:  # 'auto'
        spinor_merge = _detect_spinor(occup)
    if spinor_merge and nb % 2 != 0:
        print(f"  WARNING: odd band count ({nb}) — cannot pair spin sub-bands, "
              "plotting bands as-is")
        spinor_merge = False
    if spinor_merge:
        n_lvl = nb // 2
        lvl_occ = occup[0::2, :] + occup[1::2, :]     # spins summed per level
        n_occ_lvl = int(np.nanmax(np.sum(lvl_occ > 0.1, axis=0)))
        print(f"  spinor input detected (occupation <= 1 per band): "
              f"summing spin pairs -> {n_lvl} levels "
              f"({n_occ_lvl} occupied, {2:.0f} e- per occupied level)")

    eigen_ev = (eigen_ha - vbm_ha) * HA_TO_EV      # shift VBM → 0, convert

    # --- build path -----------------------------------------------------
    full_path  = []    # list of {dist, energies[nb]}
    node_dists = [0.0]
    cum_dist   = 0.0
    steps      = 40

    for seg in range(len(path_labels) - 1):
        la, lb  = path_labels[seg], path_labels[seg + 1]
        pa = np.asarray(hs_points[la], dtype=float)
        pb = np.asarray(hs_points[lb], dtype=float)
        seg_len = np.linalg.norm(pb - pa)
        print(f"  {la} → {lb}  (|Δk| = {seg_len:.4f} r.l.u.)")

        for s in range(steps + 1):
            ideal  = pa + (s / steps) * (pb - pa)
            kid, _ = _nearest_k(ideal, k_db, snap_tol)
            if kid is None:
                continue
            if full_path and full_path[-1]['kid'] == kid:
                continue
            full_path.append({
                'dist':     cum_dist + (s / steps) * seg_len,
                'energies': eigen_ev[:, kid - 1],
                'kid':      kid,
            })

        cum_dist += seg_len
        node_dists.append(cum_dist)

    if not full_path:
        print(f"  WARNING: no k-points mapped (snap_tol = {snap_tol:.4f}). "
              "Try a finer k-grid.")
        return

    print(f"  {len(full_path)} unique k-points on path")

    # --- plot -----------------------------------------------------------
    dists = np.array([p['dist'] for p in full_path])
    bands = np.array([p['energies'] for p in full_path])    # (npath, nb)

    fig, ax = plt.subplots(figsize=(7, 6))

    if spinor_merge:
        # Faint spin-resolved sub-bands underneath ...
        for b in range(bands.shape[1]):
            ax.plot(dists, bands[:, b], '-', color='#7799cc', lw=0.5, alpha=0.45)
        # ... one solid curve per level (Kramers pair, spins summed)
        levels = 0.5 * (bands[:, 0::2] + bands[:, 1::2])
        for b in range(levels.shape[1]):
            ax.plot(dists, levels[:, b], 'k-', lw=0.9, alpha=0.75)
        ax.plot([], [], '-', color='#7799cc', lw=0.8,
                label='spin sub-bands')
        ax.plot([], [], 'k-', lw=0.9,
                label=f'levels (spin pairs, {levels.shape[1]})')
    else:
        for b in range(bands.shape[1]):
            ax.plot(dists, bands[:, b], 'k-', lw=0.8, alpha=0.6)

    for pos in node_dists:
        ax.axvline(pos, color='#888888', linestyle='--', lw=0.7)
    ax.axhline(0.0, color='tab:red', linestyle='-', lw=0.8, alpha=0.7,
               label='VBM = 0')

    tick_labels = [r'$\Gamma$' if l == 'Gamma' else f'${l}$'
                   for l in path_labels]
    ax.set_xticks(node_dists)
    ax.set_xticklabels(tick_labels, fontsize=12)
    ax.set_ylabel('Energy (eV)', fontsize=11)
    ax.set_xlim(0.0, node_dists[-1])
    ax.set_ylim(*energy_range_ev)
    title = f'Band structure — {kfile.stem.replace("_k", "")}'
    if spinor_merge:
        title += '  (spinor: spins summed per level)'
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()

    path_str = '-'.join(path_labels)
    out = output_dir / f'band_structure_{path_str}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


# ===========================================================================
# Entry point
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Plot SALMON SBE real-time data and EPM/DFT band structure.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input-dir', default='.',
                        help='Directory containing data files')
    parser.add_argument('-o', '--output', default='sbe_plots',
                        help='Output directory for PNG files')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Image resolution (DPI)')
    parser.add_argument('--downsample', type=int, default=1,
                        help='Keep every N-th line in RT files (1 = all lines)')
    parser.add_argument('--band-path', nargs='+', default=DEFAULT_BAND_PATH,
                        metavar='PT',
                        help='High-symmetry path for band structure plot')
    parser.add_argument('--energy-range', nargs=2, type=float,
                        default=[-6.0, 12.0], metavar=('EMIN', 'EMAX'),
                        help='Energy window for band structure plot (eV)')
    parser.add_argument('--spin-sum', choices=['auto', 'on', 'off'],
                        default='auto',
                        help='Spinor (spin-orbit split) eigen files: merge '
                             'adjacent spin sub-bands into levels (occupations '
                             'summed, energy = pair mean). "auto" detects '
                             'spinor input from occupations <= 1 per band.')

    # Mode shortcuts (each implies the complementary --no-* flag)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--only-sbe', action='store_true',
                      help='Only SBE RT/energy/nex/nex_k plots — skip band structure '
                           '(legacy mode, equivalent to --no-bands)')
    mode.add_argument('--only-bands', action='store_true',
                      help='Only band structure — skip all RT/nex plots '
                           '(equivalent to --no-rt)')
    # Fine-grained toggles kept for backwards compatibility
    parser.add_argument('--no-bands', action='store_true',
                        help='Skip band structure even if *_k.data exists')
    parser.add_argument('--no-rt', action='store_true',
                        help='Skip all RT / nex plots')
    parser.add_argument('--log-cmap', action='store_true',
                        default=CMAP_LOG_SCALE,
                        help='Use logarithmic colormap scaling for 2-D population '
                             'plots (snapshots + time-k maps). '
                             f'Hardcoded default: CMAP_LOG_SCALE = {CMAP_LOG_SCALE}')
    parser.add_argument('--snapshots', action='store_true',
                        default=SNAP_ENABLED,
                        help='Write one k-space snapshot PNG per time block in '
                             'nex_k (3 projected planes). Time-k maps are always '
                             'written regardless of this flag. '
                             f'Hardcoded default: SNAP_ENABLED = {SNAP_ENABLED}')
    args = parser.parse_args()

    # Resolve mode shortcuts
    if args.only_sbe:
        args.no_bands = True
    if args.only_bands:
        args.no_rt = True

    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    found_any = False

    # --- RT line files --------------------------------------------------
    if not args.no_rt:
        for pattern in ('*_sbe_rt.data', '*_sbe_rt_energy.data', '*_sbe_nex.data'):
            for f in sorted(input_dir.glob(pattern)):
                found_any = True
                plot_rt_file(f, output_dir, downsample=args.downsample, dpi=args.dpi)

        for f in sorted(input_dir.glob('*_sbe_nex_k.data')):
            found_any = True
            plot_nex_k(f, output_dir, dpi=args.dpi,
                       log_scale=args.log_cmap, snapshots=args.snapshots)

    # --- Band structure -------------------------------------------------
    if not args.no_bands:
        for kf in sorted(input_dir.glob('*_k.data')):
            stem = kf.name[:-len('_k.data')]
            ef   = kf.parent / f'{stem}_eigen.data'
            if not ef.exists():
                print(f"  (skip bands) {ef.name} not found alongside {kf.name}")
                continue
            found_any = True
            try:
                plot_band_structure(
                    kf, ef, output_dir,
                    path_labels=args.band_path,
                    energy_range_ev=tuple(args.energy_range),
                    dpi=args.dpi, spin_sum=args.spin_sum)
            except Exception as exc:
                print(f"  ERROR in band structure for {kf.name}: {exc}")

    if not found_any:
        print(f"No data files found in {input_dir.resolve()}")
        print("Expected: *_sbe_rt.data, *_sbe_rt_energy.data, *_sbe_nex.data, "
              "*_sbe_nex_k.data, *_k.data + *_eigen.data")
        return

    print(f"\nDone.  Output: {output_dir.resolve()}")


if __name__ == '__main__':
    main()
