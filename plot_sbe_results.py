#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified plotter for SALMON-SBE real-time output and EPM/DFT ground-state files.

Drop into a calculation directory and run:

    python3 plot_sbe_results.py                     # auto-detect everything
    python3 plot_sbe_results.py --downsample 200    # thin out large RT curves
    python3 plot_sbe_results.py --band-path L Gamma X W K
    python3 plot_sbe_results.py --lattice wurtzite          # CdS (orthorhombic)

What is plotted
---------------
  *_sbe_rt.data          : fields + current vs time  (downsampled if requested)
                           + optical conductivity sigma(w) = J(w)/E(w) along the
                             driven axis: a global spectrum (Re and Im on two
                             Y-axes, 0-4 THz) and a strongly-overlapped STFT
                             Re-sigma(w,t) 2-D map (--no-conductivity to skip;
                             --fmax-thz / --stft-window-fs / --stft-hop to tune).
                             True THz resolution needs a ps-scale run.
  *_sbe_rt_energy.data   : total energy vs time
  *_sbe_nex.data         : excited electron count vs time
  *_sbe_nex_k.data       : per-k Houston-basis LCB population:
                             snapshot PNGs (3 projected planes) + time-k maps
  *_k.data + *_eigen.data: band structure along the requested path
                             k in reduced coords, energy shifted to VBM = 0 eV
  band.dat               : band structure from a theory='dft_band' run,
                             plotted vs path distance (energy shifted to the
                             band index given by --band-vbm, default nb//2)

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
import warnings
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from scipy.interpolate import RegularGridInterpolator
import matplotlib.colors as mcolors

# The wrapped FCC primitive-zone grid is irregular, so some 1-D marginal slices
# are all-NaN and nanmean legitimately returns NaN (rendered blank). Silence the
# accompanying numpy warning -- it is expected, not an error.
warnings.filterwarnings('ignore', message='Mean of empty slice')

plt.switch_backend('Agg')

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HA_TO_EV = 27.211386245988     # Hartree → eV
AU_TIME_FS = 0.02418884326505  # 1 atomic unit of time in fs
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

# FCC high-symmetry points in REDUCED coordinates of the FCC PRIMITIVE
# reciprocal basis b1=(2pi/a)(-1,1,1), b2=(2pi/a)(1,-1,1), b3=(2pi/a)(1,1,-1).
# For the folded cubic-cell data they are converted to Cartesian (in 2pi/a
# units = simple-cubic reduced coordinates) via _fcc_prim_to_sc_reduced() and
# wrapped into the cubic BZ before snapping to the MP grid. NOTE: under the
# 4-fold folding several FCC points land on the same cubic star (e.g.
# X = (2pi/a)(1,0,0) wraps onto Gamma), so the folded plot overlays the
# states of up to 4 primitive BZ points at every tick -- use the unfolded
# *_bandpath.data plot for a clean primitive-cell picture.
HS_POINTS_FCC = {
    'Gamma': [ 0.000,  0.000,  0.000],
    'X':     [ 0.000,  0.500,  0.500],
    'L':     [ 0.500,  0.500,  0.500],
    'W':     [ 0.250,  0.500,  0.750],
    'K':     [ 0.375,  0.375,  0.750],
    'U':     [ 0.250,  0.625,  0.625],
}
HS_POINTS = HS_POINTS_FCC          # back-compat alias (cubic GaAs/Si default)

# Wurtzite CdS is carried as an 8-atom ORTHORHOMBIC (sqrt3 x 1 x 1) supercell,
# al(1:3) = (a, a*sqrt3, c). Its *_k.data is in orthorhombic-cell reduced
# coordinates, so the high-symmetry points are the standard orthorhombic-BZ
# points (unambiguous in that frame); the hexagonal correspondence is noted.
#   Z == A  (zone top along the c-axis [0001]: the inversion-breaking direction
#            where the wurtzite even harmonics live); X,Y,S relate to M/K.
HS_POINTS_WZ = {
    'Gamma': [0.0, 0.0, 0.0],
    'X':     [0.5, 0.0, 0.0],   # along a* (the sqrt3-short direction)
    'Y':     [0.0, 0.5, 0.0],   # along b*
    'A':     [0.0, 0.0, 0.5],   # == hexagonal A (Gamma-A along c)
    'S':     [0.5, 0.5, 0.0],   # in-plane corner
    'U':     [0.5, 0.0, 0.5],
    'T':     [0.0, 0.5, 0.5],
    'R':     [0.5, 0.5, 0.5],
}
# Rows of the FCC primitive reciprocal basis in 2pi/a (= sc-reduced) units.
_B_FCC_RED = np.array([[-1.0, 1.0, 1.0],
                       [ 1.0, -1.0, 1.0],
                       [ 1.0, 1.0, -1.0]])
# Wurtzite-CdS orthorhombic reciprocal-length weights (relative |b_i| = 1/L_i,
# L = a : a*sqrt3 : c with c/a = 1.623): used only to space band-path segments
# by physical |dk| on the anisotropic cell.
_WZ_RECIP_W = np.array([1.0, 1.0 / np.sqrt(3.0), 1.0 / 1.623])

def _fcc_prim_to_sc_reduced(q):
    """FCC-primitive reduced coordinates -> simple-cubic reduced (k/(2pi/a))."""
    return np.asarray(q, dtype=float) @ _B_FCC_RED

def _wz_orth_to_cart(q):
    """Orthorhombic reduced coords -> Cartesian (in 2pi/a units) for the
    wurtzite cell, weighting each axis by its physical reciprocal length so
    band-path segment distances are physically proportioned."""
    return np.asarray(q, dtype=float) * _WZ_RECIP_W

DEFAULT_BAND_PATH_FCC = ['L', 'Gamma', 'X', 'W', 'K']
DEFAULT_BAND_PATH_WZ  = ['A', 'Gamma', 'X', 'S', 'Y', 'Gamma']
DEFAULT_BAND_PATH = DEFAULT_BAND_PATH_FCC

# Per-lattice plotting context: symmetry points, default path, reduced->Cartesian
# transform, and the point-group sign/permutation set used to snap a path point
# to the MP grid. FCC allows cubic permutations+signs; orthorhombic (wurtzite
# supercell) allows sign flips only (the three axes are inequivalent).
def _lattice_context(lattice):
    if lattice == 'wurtzite':
        return dict(hs=HS_POINTS_WZ, path=DEFAULT_BAND_PATH_WZ,
                    to_cart=_wz_orth_to_cart, permute=False)
    return dict(hs=HS_POINTS_FCC, path=DEFAULT_BAND_PATH_FCC,
                to_cart=_fcc_prim_to_sc_reduced, permute=True)


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
# Optical conductivity  sigma(omega) = J(omega)/E(omega)   (*_sbe_rt.data)
# ===========================================================================

def _parse_header_units(header_line):
    """Return [(name, unit), ...] from a numbered '# 1:Time[fs] 2:...' header."""
    return re.findall(r'\d+:([^\[\s]+)\[([^\]]*)\]', header_line)


def _time_unit_fs(filepath):
    """fs per time-step-unit of column 1 of an rt file: AU_TIME_FS for an a.u.
    file, 1.0 for an fs file. Detected from the column-1 unit string; defaults
    to a.u. (the SBE solver's native unit)."""
    try:
        nu = _parse_header_units(find_header(filepath))
        if nu:
            u = nu[0][1].lower()
            if 'fs' in u:
                return 1.0
            if 'au' in u or 'a.u' in u:
                return AU_TIME_FS
    except Exception:
        pass
    return AU_TIME_FS


def _rt_drive_axis(cols, data):
    """Pick the Cartesian axis the field drives (largest peak |E_tot|) and return
    (axis_label, t, E_tot[axis], Jm[axis]). Falls back to E_ext if E_tot absent."""
    name_to_col = {n: i for i, n in enumerate(cols)}
    etag = 'E_tot' if 'E_tot_x' in name_to_col else 'E_ext'
    t = data[:, 0]
    best, best_amp, E, J = 'z', -1.0, None, None
    for ax in ('x', 'y', 'z'):
        ecol = name_to_col.get(f'{etag}_{ax}')
        jcol = name_to_col.get(f'Jm_{ax}')
        if ecol is None or jcol is None:
            continue
        amp = float(np.nanmax(np.abs(data[:, ecol])))
        if amp > best_amp:
            best, best_amp, E, J = ax, amp, data[:, ecol], data[:, jcol]
    return best, t, E, J


def _sigma_ratio(Jw, Ew, eta=1e-3):
    """Regularised spectral ratio sigma = J conj(E) / (|E|^2 + (eta*max|E|)^2).
    Avoids blow-ups where the drive spectrum E(omega) has nulls; reduces to
    J/E where |E| is appreciable."""
    floor = (eta * np.nanmax(np.abs(Ew)))**2
    return Jw * np.conj(Ew) / (np.abs(Ew)**2 + floor)


def plot_conductivity(filepath, output_dir, fmax_thz=4.0, dpi=150):
    """Global optical conductivity sigma(omega) = J(omega)/E(omega) along the
    driven axis, from *_sbe_rt.data. Re and Im are drawn on two Y-axes; the
    frequency axis is in THz, restricted to [0, fmax_thz]. A Hann window is
    applied to both J(t) and E(t) before the FFT (it largely cancels in the
    ratio while suppressing spectral leakage)."""
    print(f"Processing {filepath.name}  (optical conductivity sigma(w)=J/E) ...")
    cols, data = load_columns_streaming(filepath, downsample=1)
    if data.shape[0] < 8:
        print("  (skip) too few time samples")
        return
    tu_fs = _time_unit_fs(filepath)
    axis, t, E, J = _rt_drive_axis(cols, data)
    if E is None:
        print("  (skip) no E/Jm columns found")
        return
    n = len(t)
    dt = float(np.mean(np.diff(t)))
    if dt <= 0:
        print("  (skip) non-monotonic time column")
        return
    w = np.hanning(n)
    f = np.fft.rfftfreq(n, d=dt)                 # 1/time-unit
    f_thz = f * 1000.0 / tu_fs                   # -> THz
    Ew = np.fft.rfft(E * w)
    Jw = np.fft.rfft(J * w)
    sigma = _sigma_ratio(Jw, Ew)

    mask = f_thz <= fmax_thz
    if np.count_nonzero(mask) < 3:
        # trace too short to resolve the requested band -> show the lowest bins
        keep = min(max(8, 3), len(f_thz))
        mask = np.zeros_like(f_thz, dtype=bool); mask[:keep] = True
        print(f"  NOTE: trace too short for 0-{fmax_thz:g} THz "
              f"(df = {f_thz[1]-f_thz[0]:.3g} THz); showing the lowest {keep} bins. "
              f"Use a longer run (ps-scale) for true THz resolution.")
    fx = f_thz[mask]
    re = np.real(sigma[mask])
    im = np.imag(sigma[mask])

    fig, ax_re = plt.subplots(figsize=(10, 5))
    c_re, c_im = 'tab:blue', 'tab:red'
    l1, = ax_re.plot(fx, re, '-', color=c_re, lw=1.3, label=r'Re $\sigma$')
    ax_re.set_xlabel('Frequency [THz]')
    ax_re.set_ylabel(r'Re $\sigma(\omega)$  [a.u.]', color=c_re)
    ax_re.tick_params(axis='y', labelcolor=c_re)
    ax_re.axhline(0.0, color='0.7', lw=0.6)
    ax_re.set_xlim(0.0, fmax_thz if np.count_nonzero(mask) >= 3 else fx[-1])
    ax_re.grid(True, alpha=0.3, linestyle='--')

    ax_im = ax_re.twinx()
    l2, = ax_im.plot(fx, im, '-', color=c_im, lw=1.3, label=r'Im $\sigma$')
    ax_im.set_ylabel(r'Im $\sigma(\omega)$  [a.u.]', color=c_im)
    ax_im.tick_params(axis='y', labelcolor=c_im)

    ax_re.set_title(f'Optical conductivity  $\\sigma(\\omega)=J_{axis}(\\omega)/E_{axis}(\\omega)$'
                    f'  —  {filepath.stem}')
    ax_re.legend(handles=[l1, l2], loc='upper right', fontsize=9)
    fig.tight_layout()
    out = output_dir / f'{filepath.stem}_conductivity.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}  (driven axis: {axis})")


def plot_conductivity_stft(filepath, output_dir, fmax_thz=4.0,
                           window_fs=None, hop=1, dpi=150, max_cols=1000):
    """Short-time-Fourier 2-D map of Re sigma(omega, t) along the driven axis.
    Strong overlap: hop defaults to 1 sample (N-1 of N overlap). The window
    length defaults to a quarter of the trace (override with --stft-window-fs).
    The number of rendered time columns is capped (the effective hop is raised
    if a 1-sample hop would exceed it) so the PNG stays reasonable while the
    overlap remains as strong as practical. Frequency axis restricted to
    [0, fmax_thz]; only the real part is mapped (the dissipative conductivity)."""
    print(f"Processing {filepath.name}  (STFT Re sigma(w,t) 2-D map) ...")
    cols, data = load_columns_streaming(filepath, downsample=1)
    if data.shape[0] < 16:
        print("  (skip) too few time samples for an STFT")
        return
    tu_fs = _time_unit_fs(filepath)
    axis, t, E, J = _rt_drive_axis(cols, data)
    if E is None:
        print("  (skip) no E/Jm columns found")
        return
    n = len(t)
    dt = float(np.mean(np.diff(t)))
    if dt <= 0:
        print("  (skip) non-monotonic time column")
        return

    if window_fs is not None and window_fs > 0:
        N = int(round(window_fs / (dt * tu_fs)))
    else:
        N = n // 4
    N = max(8, min(N, n))
    if N >= n:
        print("  (skip) window not shorter than the trace")
        return

    hop = max(1, int(hop))
    n_full = (n - N) // hop + 1
    # cap the rendered columns: raise the effective hop if hop=1 over-resolves.
    eff_hop = hop
    if n_full > max_cols:
        eff_hop = hop * int(np.ceil(n_full / max_cols))
        print(f"  NOTE: {n_full} windows at hop={hop} exceeds {max_cols}; "
              f"rendering with effective hop={eff_hop} (overlap still "
              f"{100*(1-eff_hop/N):.1f}%).")
    starts = list(range(0, n - N + 1, eff_hop))

    f = np.fft.rfftfreq(N, d=dt)
    f_thz = f * 1000.0 / tu_fs
    fmask = f_thz <= fmax_thz
    if np.count_nonzero(fmask) < 2:
        keep = min(8, len(f_thz)); fmask = np.zeros_like(f_thz, dtype=bool)
        fmask[:keep] = True
        print(f"  NOTE: window too short for 0-{fmax_thz:g} THz "
              f"(df = {f_thz[1]-f_thz[0]:.3g} THz); showing the lowest {keep} bins. "
              f"Use a longer --stft-window-fs (ps-scale).")
    f_keep = f_thz[fmask]

    w = np.hanning(N)
    centers_fs = np.empty(len(starts))
    nfk = np.count_nonzero(fmask)
    Smap = np.empty((nfk, len(starts)))
    Emag = np.empty((nfk, len(starts)))    # local drive power per (f, t) pixel
    for i, s in enumerate(starts):
        seg = slice(s, s + N)
        Ew = np.fft.rfft(E[seg] * w)
        Jw = np.fft.rfft(J[seg] * w)
        Smap[:, i] = np.real(_sigma_ratio(Jw, Ew))[fmask]
        Emag[:, i] = np.abs(Ew)[fmask]
        centers_fs[i] = (t[s] + 0.5 * N * dt) * tu_fs

    # sigma = J/E is only meaningful where the local drive has power. Blank the
    # pixels whose |E(omega,t)| is a tiny fraction of the global peak -- this is
    # where the regularised ratio would otherwise produce meaningless spikes.
    emax = float(np.nanmax(Emag)) if Emag.size else 0.0
    if emax > 0:
        Smap[Emag < 1e-2 * emax] = np.nan
    # robust symmetric scale (99th percentile, not the max) so a few residual
    # edge spikes don't wash out the real conductivity.
    finite = Smap[np.isfinite(Smap)]
    vmax = float(np.nanpercentile(np.abs(finite), 99)) if finite.size else 1.0
    vmax = max(vmax, 1e-30)
    fig, ax = plt.subplots(figsize=(max(8, len(starts) * 0.02 + 3), 5))
    im = ax.pcolormesh(_bin_edges(centers_fs), _bin_edges(f_keep), Smap,
                       cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='flat')
    plt.colorbar(im, ax=ax, label=r'Re $\sigma(\omega, t)$  [a.u.]')
    ax.set_xlabel('time [fs]')
    ax.set_ylabel('Frequency [THz]')
    ax.set_ylim(0.0, fmax_thz if np.count_nonzero(fmask) >= 2 else f_keep[-1])
    win_fs = N * dt * tu_fs
    ax.set_title(f'STFT Re $\\sigma(\\omega,t)$ (axis {axis}, window {win_fs:.1f} fs, '
                 f'{len(starts)} cols)  —  {filepath.stem}')
    fig.tight_layout()
    out = output_dir / f'{filepath.stem}_conductivity_stft.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


# ===========================================================================
# Streaming nex_k  (*_sbe_nex_k.data)
# ===========================================================================

# Physical-level columns of the unfolded population file, in file order.
UNFOLD_LEVELS = ('vbm1', 'vb', 'cb1', 'cb2')
_UNFOLD_PRIMARY = 'cb1'   # level used for the dynamics k-maps / baseline


def _iter_nex_k_blocks(filepath, unfold=False):
    """Yield (t_val, t_unit, kpoints[nk,3], pop[nk], levels) one block at a time.

    folded  (*_sbe_nex_k.data):        columns ik, kx, ky, kz, pop
        -> pop is the LCB population, levels is None.
    unfold  (*_sbe_nex_k_unfold.data): columns ik, isub, kx, ky, kz,
        pop_vbm1, pop_vb, pop_cb1, pop_cb2
        -> pop is the CB1 population (for the existing k-maps), levels is a
           dict {name: ndarray[nk]} with all four physical-band populations."""
    time_re = re.compile(r'#\s*t\s*=\s*([-+\d.eEdD]+)\s*(\S*)')
    icol = 2 if unfold else 1
    ncol = 9 if unfold else 5
    t_value, t_unit = None, ''
    kx, ky, kz = [], [], []
    sub = []      # unfold: the sublattice index (column 2), enables the folded view
    lev = {name: [] for name in UNFOLD_LEVELS} if unfold else None
    pop = []      # folded: the single population column

    def _emit():
        kp = np.column_stack([kx, ky, kz])
        sub_arr = np.asarray(sub, dtype=int) if unfold else None
        if unfold:
            levels = {name: np.asarray(lev[name], dtype=float) for name in UNFOLD_LEVELS}
            return (t_value, t_unit, kp, levels[_UNFOLD_PRIMARY], levels, sub_arr)
        return (t_value, t_unit, kp, np.asarray(pop, dtype=float), None, None)

    with open(filepath, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            m = time_re.match(s)
            if m:
                if t_value is not None and kx:
                    yield _emit()
                t_value = float(m.group(1))
                t_unit  = m.group(2)
                kx, ky, kz = [], [], []
                pop = []
                sub = []
                if unfold:
                    lev = {name: [] for name in UNFOLD_LEVELS}
                continue
            if s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < ncol:
                continue
            try:
                kx.append(float(parts[icol])); ky.append(float(parts[icol + 1]))
                kz.append(float(parts[icol + 2]))
                if unfold:
                    sub.append(int(float(parts[1])))
                    for j, name in enumerate(UNFOLD_LEVELS):
                        lev[name].append(float(parts[icol + 3 + j]))
                else:
                    pop.append(float(parts[icol + 3]))
            except (ValueError, IndexError):
                continue

    if t_value is not None and kx:
        yield _emit()


def _wrap_to_fcc_bz(kpoints):
    """Wrap k-points (sc-reduced units, i.e. Cartesian/(2pi/a)) into the
    first FCC primitive BZ: subtract the FCC reciprocal vector (all-equal-
    parity integer triplet) closest to each point."""
    cands = np.array([g for g in itertools.product((-2, -1, 0, 1, 2), repeat=3)
                      if (g[0] - g[1]) % 2 == 0 and (g[1] - g[2]) % 2 == 0],
                     dtype=float)
    d = kpoints[:, None, :] - cands[None, :, :]
    best = np.argmin((d**2).sum(axis=2), axis=1)
    return kpoints - cands[best]


def _read_unfold_offsets(filepath):
    """Parse the '# isub, offset G0 (sc reduced)' table from an unfold file.
    Returns {isub: ndarray[3]}; falls back to the canonical FCC sublattice
    offsets if the table is absent."""
    offs = {}
    with open(filepath, 'r') as f:
        in_block = False
        for line in f:
            s = line.strip()
            if s.startswith('# isub'):
                in_block = True
                continue
            if in_block:
                if s.startswith('#') or not s:
                    break
                p = s.split()
                if len(p) == 4:
                    offs[int(p[0])] = np.array([float(p[1]), float(p[2]), float(p[3])])
                else:
                    break
    if len(offs) != 4:
        offs = {1: np.zeros(3), 2: np.array([1., 0., 0.]),
                3: np.array([0., 1., 0.]), 4: np.array([0., 0., 1.])}
    return offs


def _fold_to_cubic(kpoints_prim, sub, offsets):
    """Recover the cubic supercell k (k_sc = k_prim - G0(isub)) and wrap it
    into [-0.5, 0.5). Summing populations over the four sublattices at a fixed
    k_sc collapses the FCC valleys (Gamma + the three X points) back onto the
    regular cubic grid -- the clean single-zone view."""
    g0 = np.array([offsets[int(s)] for s in sub])
    ksc = kpoints_prim - g0
    return ksc - np.round(ksc)


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
                   log_scale=False, tag='nex_k'):
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
    out = output_dir / f'{tag}_snap_t{safe_t}{t_unit}.png'
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
                 log_scale=False, tag='nex_k'):
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
    out = output_dir / f'{tag}_ktmap_{label_k}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def plot_nex_k(filepath, output_dir, dpi=150, log_scale=False, snapshots=False,
               unfold=False, subtract_baseline=False):
    print(f"Processing {filepath.name}  "
          f"(cmap={'log' if log_scale else 'linear'}, "
          f"snapshots={'on' if snapshots else 'off'}"
          f"{', unfolded primitive BZ' if unfold else ''}"
          f"{', baseline-subtracted' if subtract_baseline else ''}) ...")
    tag = 'nex_k_unfold' if unfold else 'nex_k'
    if subtract_baseline:
        tag += '_db'
    # The unfolded primitive-zone map legitimately shows the CB1 population of
    # every FCC valley (Gamma for sublattice 1, the X points for 2/3/4), so it
    # carries satellite peaks at the zone boundary. The FOLDED view sums the
    # four sublattices back onto the regular cubic grid (k_sc = k_prim - G0),
    # collapsing the valleys into a single clean zone -- the per-cubic-k total
    # of the lowest conduction band.
    offsets = _read_unfold_offsets(filepath) if unfold else None
    ftag = 'nex_k_fold' + ('_db' if subtract_baseline else '')

    kx_u = ky_u = kz_u = ix = iy = iz = None
    pop3d = None
    fkx_u = fky_u = fkz_u = fix = fiy = fiz = None
    fpop3d = None
    pop_baseline = None       # per-k population of the first non-zero frame
    times, marg_kx, marg_ky, marg_kz = [], [], [], []
    fmarg_kx, fmarg_ky, fmarg_kz = [], [], []
    t_unit_last = ''
    n_blocks = 0

    for t_val, t_unit, kpoints, pop, _levels, sub in _iter_nex_k_blocks(filepath, unfold=unfold):
        t_unit_last = t_unit
        n_blocks   += 1
        kpoints_prim = kpoints
        if unfold:
            kpoints = _wrap_to_fcc_bz(kpoints)
        # Optional baseline removal: subtract the first frame that carries any
        # population. This is the gauge/Houston-projection O(A^2) offset the
        # crystal-gauge projection leaves at finite vector potential; removing
        # it isolates the genuine field-driven excitation. Population can't be
        # negative, so the difference is clipped at zero.
        if subtract_baseline:
            if pop_baseline is None and np.nanmax(pop) > 1e-30:
                pop_baseline = pop.copy()
            if pop_baseline is not None:
                pop = np.clip(pop - pop_baseline, 0.0, None)
        if kx_u is None:
            kx_u, ky_u, kz_u, ix, iy, iz = _build_grid_info(kpoints)
            pop3d = np.empty((len(kx_u), len(ky_u), len(kz_u)))
        pop3d.fill(np.nan)
        pop3d[ix, iy, iz] = pop
        if snapshots:
            _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi,
                           log_scale=log_scale, tag=tag)
        times.append(t_val)
        marg_kx.append(np.nanmean(pop3d, axis=(1, 2)))
        marg_ky.append(np.nanmean(pop3d, axis=(0, 2)))
        marg_kz.append(np.nanmean(pop3d, axis=(0, 1)))

        # Folded cubic-zone view: sum the four sublattices at each k_sc.
        if unfold and sub is not None:
            ksc = _fold_to_cubic(kpoints_prim, sub, offsets)
            fpop = pop  # baseline (if any) already applied to `pop` above
            if fkx_u is None:
                fkx_u, fky_u, fkz_u, fix, fiy, fiz = _build_grid_info(ksc)
                fpop3d = np.zeros((len(fkx_u), len(fky_u), len(fkz_u)))
            fpop3d.fill(0.0)
            np.add.at(fpop3d, (fix, fiy, fiz), fpop)
            if snapshots:
                _save_snapshot(fpop3d, fkx_u, fky_u, fkz_u, t_val, t_unit, output_dir,
                               dpi, log_scale=log_scale, tag=ftag)
            fmarg_kx.append(np.nanmean(fpop3d, axis=(1, 2)))
            fmarg_ky.append(np.nanmean(fpop3d, axis=(0, 2)))
            fmarg_kz.append(np.nanmean(fpop3d, axis=(0, 1)))

    if n_blocks == 0:
        print("  (skip) no data blocks found")
        return

    print(f"  writing time-k maps ({n_blocks} time steps) ...")
    _save_kt_map(times, t_unit_last, kx_u, 'kx', marg_kx, output_dir, dpi,
                 log_scale=log_scale, tag=tag)
    _save_kt_map(times, t_unit_last, ky_u, 'ky', marg_ky, output_dir, dpi,
                 log_scale=log_scale, tag=tag)
    _save_kt_map(times, t_unit_last, kz_u, 'kz', marg_kz, output_dir, dpi,
                 log_scale=log_scale, tag=tag)

    if unfold and fmarg_kx:
        print(f"  writing folded cubic-zone time-k maps ...")
        _save_kt_map(times, t_unit_last, fkx_u, 'kx', fmarg_kx, output_dir, dpi,
                     log_scale=log_scale, tag=ftag)
        _save_kt_map(times, t_unit_last, fky_u, 'ky', fmarg_ky, output_dir, dpi,
                     log_scale=log_scale, tag=ftag)
        _save_kt_map(times, t_unit_last, fkz_u, 'kz', fmarg_kz, output_dir, dpi,
                     log_scale=log_scale, tag=ftag)


def _unfold_peak_levels(filepath):
    """Stream the unfold file once and return the frame of peak excitation
    (largest total CB population), as (t_val, t_unit, kpoints[nk,3] wrapped to
    the FCC BZ in sc-reduced units, levels={name: pop[nk]})."""
    best = None
    best_sum = -1.0
    for t_val, t_unit, kpoints, _pop, levels, _sub in _iter_nex_k_blocks(filepath, unfold=True):
        if levels is None:
            return None
        s = float(np.nansum(levels['cb1']) + np.nansum(levels['cb2']))
        if s > best_sum:
            best_sum = s
            best = (t_val, t_unit, _wrap_to_fcc_bz(kpoints),
                    {k: v.copy() for k, v in levels.items()})
    return best


def _bandpath_level_energies(eig_ha, nv, spinor):
    """Return {name: (E_central_ev[N], kinetic_ev[N])} for VB-1, VB, CB1, CB2,
    referenced to the valence-band maximum. The central energy is the mean of
    the spin doublet (spinor). kinetic_ev is the carrier kinetic energy used as
    the line broadening: E - E_CBM for the conduction bands (electron kinetic
    energy, the variable of the Stobbe impact-ionization rate) and E_VBM - E for
    the valence bands (hole kinetic energy)."""
    n_spin = 2 if spinor else 1
    nv_phys = nv // n_spin
    slot_phys = {'vbm1': nv_phys - 1, 'vb': nv_phys,
                 'cb1': nv_phys + 1, 'cb2': nv_phys + 2}

    def _level(lvl):                       # 1-based physical level -> mean energy in Ha
        if spinor:
            return 0.5 * (eig_ha[:, 2 * lvl - 2] + eig_ha[:, 2 * lvl - 1])
        return eig_ha[:, lvl - 1]

    vbm_ha = np.nanmax(_level(nv_phys))
    cbm_ha = np.nanmin(_level(nv_phys + 1))
    out = {}
    for name, lvl in slot_phys.items():
        if lvl < 1 or (spinor and 2 * lvl > eig_ha.shape[1]) or \
           (not spinor and lvl > eig_ha.shape[1]):
            continue
        e_ev = (_level(lvl) - vbm_ha) * HA_TO_EV
        if name in ('cb1', 'cb2'):
            kinetic = np.clip((_level(lvl) - cbm_ha) * HA_TO_EV, 0.0, None)
        else:
            kinetic = np.clip((vbm_ha - _level(lvl)) * HA_TO_EV, 0.0, None)
        out[name] = (e_ev, kinetic)
    return out


def _grid_spacing(grid_kpts, sample=400):
    """Median nearest-neighbour distance of the (wrapped) unfolded k-grid, used
    to size the path<->grid matching tolerance for any grid density."""
    n = grid_kpts.shape[0]
    idx = np.arange(n) if n <= sample else \
        np.random.default_rng(0).choice(n, sample, replace=False)
    sub = grid_kpts[idx]
    d2 = ((sub[:, None, :] - grid_kpts[None, :, :]) ** 2).sum(axis=2)
    d2[d2 < 1e-12] = np.inf            # drop self-distance
    return float(np.median(np.sqrt(d2.min(axis=1))))


def _map_path_population(qred, grid_kpts, grid_pop, max_dist=None):
    """For each band-path point (FCC-reduced qred) return the population of the
    nearest unfolded grid point (sc-reduced, wrapped). Points farther than
    max_dist (default ~1 grid spacing) from any grid point get NaN. On a coarse
    grid the tolerance grows automatically so the path is still coloured."""
    if max_dist is None:
        max_dist = 1.1 * _grid_spacing(grid_kpts)
    q_sc = _wrap_to_fcc_bz(_fcc_prim_to_sc_reduced(qred))
    d2 = ((q_sc[:, None, :] - grid_kpts[None, :, :]) ** 2).sum(axis=2)
    j = np.argmin(d2, axis=1)
    pop = grid_pop[j].astype(float)
    pop[np.sqrt(d2[np.arange(len(j)), j]) > max_dist] = np.nan
    return pop


def plot_unfold_spectral(filepath, bpfile, output_dir, dpi=150, max_frames=150):
    """A(k,E)-style spectral plots of the unfolded population, ONE PER TIME FRAME.

    Skeleton: the clean primitive-cell dispersion from *_bandpath.data (VB-1,
    VB, CB1, CB2; spins summed). Decoration: at EVERY output time, each band is
    coloured by its physical population (mapped from the nearest unfolded
    MP-grid point) and broadened by the per-k carrier kinetic energy. For each
    frame two views are written into a `spectral_frames/` subfolder -- along the
    high-symmetry path and projected onto kx -- so the band dynamics can be
    watched frame by frame (assemble into a movie with e.g. ffmpeg/imagemagick).
    The colour scale is fixed across all frames (global CB peak) so frames are
    directly comparable. Energies come entirely from the band path."""
    from matplotlib.collections import LineCollection
    print(f"Processing {filepath.name}  (spectral A(k,E) per frame, skeleton {bpfile.name}) ...")

    dist, eig_ha, nv, spinor, nodes, qred = _load_bandpath(bpfile)
    if not nv:
        print("  (skip) band path has no nv header — cannot identify levels")
        return
    levels = _bandpath_level_energies(eig_ha, nv, spinor)
    kx_path = _fcc_prim_to_sc_reduced(qred)[:, 0]
    ke_max = max((np.nanmax(levels[n][1]) for n in levels), default=1.0)
    ke_max = max(ke_max, 1e-9)

    # Pass 1: global colour scale (peak CB population over all frames) + count.
    cb_peak, n_frames = 0.0, 0
    for _t, _tu, _kp, _pop, grid_lev, _sub in _iter_nex_k_blocks(filepath, unfold=True):
        if grid_lev is None:
            print("  (skip) unfold file has no per-level population columns")
            return
        n_frames += 1
        cb_peak = max(cb_peak, max((np.nanmax(grid_lev[n])
                                    for n in ('cb1', 'cb2') if n in grid_lev), default=0.0))
    if n_frames == 0:
        print("  (skip) no data blocks found")
        return
    norm = mcolors.Normalize(vmin=0.0, vmax=max(cb_peak, 1e-12))
    stride = max(1, n_frames // max_frames)   # cap the number of rendered frames

    frame_dir = output_dir / 'spectral_frames'
    frame_dir.mkdir(parents=True, exist_ok=True)

    def _draw(ax, xvals, grid_kpts, grid_lev, scatter):
        last_im = None
        for name in UNFOLD_LEVELS:
            if name not in levels or name not in grid_lev:
                continue
            e_ev, kinetic = levels[name]
            pop = _map_path_population(qred, grid_kpts, grid_lev[name])
            lw = 1.0 + 6.0 * np.nan_to_num(kinetic) / ke_max
            if scatter:
                good = np.isfinite(pop)
                last_im = ax.scatter(xvals[good], e_ev[good], c=pop[good],
                                     s=8 + 60 * lw[good], cmap=CMAP_POP, norm=norm,
                                     edgecolors='none', alpha=0.85)
            else:
                pts = np.column_stack([xvals, e_ev])
                segs = np.stack([pts[:-1], pts[1:]], axis=1)
                seg_pop = 0.5 * (np.nan_to_num(pop[:-1]) + np.nan_to_num(pop[1:]))
                seg_lw  = 0.5 * (lw[:-1] + lw[1:])
                ax.plot(xvals, e_ev, color='0.75', lw=0.6, zorder=1)
                lc = LineCollection(segs, cmap=CMAP_POP, norm=norm, zorder=2)
                lc.set_array(seg_pop); lc.set_linewidths(seg_lw)
                ax.add_collection(lc)
                last_im = lc
        return last_im

    # Pass 2: render every (strided) frame.
    n_written = 0
    for iframe, (t_val, t_unit, kpoints, _pop, grid_lev, _sub) in enumerate(
            _iter_nex_k_blocks(filepath, unfold=True)):
        if iframe % stride != 0:
            continue
        grid_kpts = _wrap_to_fcc_bz(kpoints)
        safe_t = f'{t_val:.4f}'.replace('-', 'm').replace('+', 'p')
        tag = f'f{iframe:04d}_t{safe_t}{t_unit}'

        # view 1: along the high-symmetry path
        fig, ax = plt.subplots(figsize=(8, 6))
        im = _draw(ax, dist, grid_kpts, grid_lev, scatter=False)
        for _lbl, d in nodes:
            ax.axvline(d, color='#888888', linestyle='--', lw=0.7)
        ax.set_xticks([d for _, d in nodes])
        ax.set_xticklabels([r'$\Gamma$' if l == 'Gamma' else f'${l}$' for l, _ in nodes])
        ax.set_xlim(dist[0], dist[-1]); ax.axhline(0.0, color='tab:red', lw=0.8, alpha=0.7)
        ax.set_ylabel('Energy [eV]  (VBM = 0)')
        ax.set_title(f'Unfolded A(k,E): VB-1/VB/CB1/CB2,  t = {t_val:.3f} {t_unit}\n'
                     f'colour = population, width = carrier kinetic energy')
        if im is not None:
            plt.colorbar(im, ax=ax, label='physical-band population')
        fig.tight_layout()
        fig.savefig(frame_dir / f'nex_k_unfold_spectral_path_{tag}.png',
                    dpi=dpi, bbox_inches='tight'); plt.close(fig)

        # view 2: projected onto kx
        fig, ax = plt.subplots(figsize=(8, 6))
        im = _draw(ax, kx_path, grid_kpts, grid_lev, scatter=True)
        ax.axhline(0.0, color='tab:red', lw=0.8, alpha=0.7)
        ax.set_xlabel('kx [reduced, FCC BZ]'); ax.set_ylabel('Energy [eV]  (VBM = 0)')
        ax.set_title(f'Unfolded population vs (kx, E),  t = {t_val:.3f} {t_unit}\n'
                     f'(band-path points; colour = population, size = kinetic energy)')
        if im is not None:
            plt.colorbar(im, ax=ax, label='physical-band population')
        fig.tight_layout()
        fig.savefig(frame_dir / f'nex_k_unfold_spectral_kx_{tag}.png',
                    dpi=dpi, bbox_inches='tight'); plt.close(fig)
        n_written += 1

    print(f"  saved {n_written} frame(s) x 2 views into {frame_dir.name}/ "
          f"(of {n_frames} time steps, stride {stride})")


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


def _sym_equivalents(point, permute=True):
    """
    Generate all distinct symmetry equivalents of *point*, each wrapped to
    (-0.5, 0.5]. permute=True (cubic FCC): all axis permutations × sign flips.
    permute=False (orthorhombic wurtzite supercell): sign flips only — the three
    axes are inequivalent, so permuting them is NOT a symmetry.
    """
    p = np.asarray(point, dtype=float)
    perms = itertools.permutations(p) if permute else (tuple(p),)
    seen, result = set(), []
    for perm in perms:
        for signs in itertools.product([1, -1], repeat=3):
            c = np.array(perm) * np.array(signs)
            c = c - np.floor(c + 0.5)          # wrap to (-0.5, 0.5]
            key = tuple(np.round(c, 8))
            if key not in seen:
                seen.add(key)
                result.append(c)
    return result


def _nearest_k(ideal, k_db, snap_tol, permute=True):
    """
    Return index (1-based) of the grid k-point nearest to *ideal*,
    accounting for the lattice symmetry equivalents of *ideal*.
    Returns (k_idx, dist) or (None, dist) if dist > snap_tol.
    """
    equivs = _sym_equivalents(ideal, permute=permute)
    best_id, best_dist = None, np.inf
    for kp in k_db:
        d = min(np.linalg.norm(kp['c'] - eq) for eq in equivs)
        if d < best_dist:
            best_dist, best_id = d, kp['id']
    return (best_id if best_dist <= snap_tol else None), best_dist


def plot_band_structure(kfile, eigenfile, output_dir,
                        path_labels=None, hs_points=None,
                        energy_range_ev=(-6, 12), dpi=150, spin_sum='auto',
                        lattice='fcc'):
    """
    Plot band structure from *_k.data and *_eigen.data.
    k-points must be in reduced (dimensionless) coordinates.
    Eigenvalues are read in Hartree, shifted to VBM = 0, plotted in eV.

    lattice: 'fcc' (cubic GaAs/Si, default) or 'wurtzite' (orthorhombic CdS
             supercell) — selects the high-symmetry points, the default path,
             the reduced->Cartesian transform, and the snap symmetry group.

    spin_sum: 'auto' — detect a spinor (spin-orbit split) dataset from the
              occupation column (max occ <= 1) and merge adjacent (Kramers
              partner) spin sub-bands into levels: occupations summed,
              level energy = mean of the pair; the spin-resolved sub-bands
              are kept as faint lines underneath.
              'on'  — force the merge, 'off' — plot all bands as-is.
    """
    ctx = _lattice_context(lattice)
    to_cart = ctx['to_cart']
    permute = ctx['permute']
    if path_labels is None:
        path_labels = ctx['path']
    if hs_points is None:
        hs_points = ctx['hs']

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
        # Convert the (lattice-specific) reduced labels to the Cartesian/reduced
        # basis of *_k.data; the ideal path points are then wrapped into the BZ
        # by _sym_equivalents() (lattice symmetry group) before snapping.
        pa = to_cart(hs_points[la])
        pb = to_cart(hs_points[lb])
        seg_len = np.linalg.norm(pb - pa)
        print(f"  {la} → {lb}  (|Δk| = {seg_len:.4f} r.l.u.)")

        for s in range(steps + 1):
            ideal  = pa + (s / steps) * (pb - pa)
            kid, _ = _nearest_k(ideal, k_db, snap_tol, permute=permute)
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
# Band structure from a dft_band run  (band.dat)
# ===========================================================================

def _load_band_dat(bandfile):
    """
    Parse the 'band.dat' written by theory='dft_band'.

    Layout:
      Number_of_Bands:            <nb>
      Number_of_kpt_in_each_block:<nk_block>
      Number_of_blocks:           <nblocks>
      <ik  kred(1:3)  kcart(1:3)>      x (nk_block*nblocks)   (7 columns)
      <ik  ib  e(spin1) [e(spin2)]>    eigenvalues, Hartree   (3-4 columns)

    Eigenvalues are written one block at a time with ik restarting at 1 each
    block, so the global path index is reconstructed from ik wrap-arounds.
    Returns (kcart[N,3], eigen_ha[N, nb], nspin).
    """
    header = {}
    coords = []          # (kx, ky, kz) Cartesian, global order
    eig_rows = []        # (global_ik, ib, [energies])
    nb = nk_block = None
    block_offset = 0
    prev_ik = 0

    with open(bandfile, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if ':' in s:                       # header line
                key, _, val = s.partition(':')
                try:
                    header[key.strip()] = int(val.strip())
                except ValueError:
                    pass
                if 'Number_of_Bands' in key:
                    nb = header.get('Number_of_Bands')
                if 'Number_of_kpt_in_each_block' in key:
                    nk_block = header.get('Number_of_kpt_in_each_block')
                continue
            parts = s.split()
            if len(parts) == 7:                # coordinate line
                try:
                    coords.append([float(parts[4]), float(parts[5]), float(parts[6])])
                except ValueError:
                    continue
            elif len(parts) in (3, 4):         # eigenvalue line
                try:
                    ik = int(parts[0]); ib = int(parts[1])
                    energies = [float(x) for x in parts[2:]]
                except ValueError:
                    continue
                if nk_block and ik < prev_ik:  # ik wrapped -> next block
                    block_offset += nk_block
                prev_ik = ik
                eig_rows.append((block_offset + ik, ib, energies))

    if nb is None or not eig_rows:
        raise ValueError(f"No band data parsed from {bandfile}")

    nspin = max(len(e) for _, _, e in eig_rows)
    nk_tot = max(gik for gik, _, _ in eig_rows)
    eigen = np.full((nk_tot, nb, nspin), np.nan)
    for gik, ib, energies in eig_rows:
        if 1 <= gik <= nk_tot and 1 <= ib <= nb:
            for isp, e in enumerate(energies):
                eigen[gik - 1, ib - 1, isp] = e

    kcart = np.array(coords) if coords else np.full((nk_tot, 3), np.nan)
    # Guard against a coordinate/eigen length mismatch (e.g. trailing fill).
    if kcart.shape[0] != nk_tot:
        n = min(kcart.shape[0], nk_tot)
        kcart = kcart[:n]
        eigen = eigen[:n]
    return kcart, eigen, nspin


def _path_distance_and_nodes(kcart):
    """
    Cumulative |Δk| distance along the Cartesian path and the indices where the
    path direction changes (segment end points), used to draw vertical guides.
    """
    n = kcart.shape[0]
    dist = np.zeros(n)
    seg = np.diff(kcart, axis=0)
    seglen = np.linalg.norm(seg, axis=1)
    dist[1:] = np.cumsum(seglen)

    nodes = [0]
    for i in range(1, n - 1):
        a, b = seg[i - 1], seg[i]
        na, nb_ = np.linalg.norm(a), np.linalg.norm(b)
        if na < 1e-10 or nb_ < 1e-10:
            nodes.append(i)
            continue
        cosang = np.dot(a, b) / (na * nb_)
        if cosang < 1.0 - 1e-4:                # direction changed -> node
            nodes.append(i)
    nodes.append(n - 1)
    # De-duplicate consecutive nodes
    nodes = sorted(set(nodes))
    return dist, nodes


def plot_band_dat(bandfile, output_dir, energy_range_ev=(-6, 12), dpi=150,
                  vbm_index=None):
    """
    Plot the band structure produced by theory='dft_band' (band.dat).

    band.dat stores no occupations, so the valence-band-maximum reference is
    taken as band index `vbm_index` (1-based); defaults to nb//2 (half filling,
    spin-degenerate), override with --band-vbm. Energies are in Hartree and
    converted to eV; the path nodes (high-symmetry points) are detected from
    direction changes and annotated with their reduced... here Cartesian-based
    distance only (labels are not stored in band.dat).
    """
    print(f"Processing dft_band output: {bandfile.name}")
    kcart, eigen_ha, nspin = _load_band_dat(bandfile)
    nk, nb = eigen_ha.shape[0], eigen_ha.shape[1]
    print(f"  {nk} path k-points, {nb} bands, {nspin} spin channel(s)")

    if vbm_index is None:
        vbm_index = nb // 2                    # half filling (spin-degenerate)
    vbm_index = max(1, min(vbm_index, nb))
    vbm_ha = np.nanmax(eigen_ha[:, vbm_index - 1, :])
    print(f"  VBM taken at band index {vbm_index}: "
          f"{vbm_ha:.6f} Ha = {vbm_ha * HA_TO_EV:.4f} eV "
          f"(override with --band-vbm)")

    eigen_ev = (eigen_ha - vbm_ha) * HA_TO_EV
    dist, nodes = _path_distance_and_nodes(kcart)

    fig, ax = plt.subplots(figsize=(7, 6))
    spin_colors = ['k', 'tab:blue']
    for isp in range(nspin):
        col = spin_colors[isp % len(spin_colors)]
        for b in range(nb):
            ax.plot(dist, eigen_ev[:, b, isp], '-', color=col, lw=0.8, alpha=0.6)
    if nspin == 2:
        ax.plot([], [], 'k-', label='spin up')
        ax.plot([], [], '-', color='tab:blue', label='spin down')

    for idx in nodes:
        ax.axvline(dist[idx], color='#888888', linestyle='--', lw=0.7)
    ax.axhline(0.0, color='tab:red', linestyle='-', lw=0.8, alpha=0.7,
               label=f'VBM (band {vbm_index}) = 0')

    # Annotate nodes with their reduced k (rounded) since labels aren't stored.
    ax.set_xticks([dist[i] for i in nodes])
    ax.set_xticklabels([str(i + 1) for i in nodes], fontsize=9)
    ax.set_xlabel('k-point index at path nodes', fontsize=10)
    ax.set_ylabel('Energy (eV)', fontsize=11)
    ax.set_xlim(0.0, dist[-1])
    ax.set_ylim(*energy_range_ev)
    ax.set_title(f'Band structure (dft_band) — {bandfile.stem}', fontsize=11)
    if nspin == 2 or True:
        ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()

    out = output_dir / f'band_dat_{bandfile.stem}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


# ===========================================================================
# Unfolded primitive-cell band path  (*_bandpath.data from epm py 'bandpath')
# ===========================================================================

def _load_bandpath(bpfile):
    """
    Parse SYSNAME_bandpath.data written by `epm_gaas_reference.py bandpath`:
      # spinor = 0/1
      # nv = <valence states per primitive k>
      # nb = <states per line>
      # nodes: LBL dist  LBL dist ...
      data: ik dist q1 q2 q3 E_1..E_nb [Ha]
    Returns (dist[N], eigen_ha[N, nb], nv, spinor, nodes=[(label, dist), ...],
             q[N, 3] in FCC-primitive reduced coords).
    """
    spinor, nv, nb, nodes = False, None, None, []
    dist, eig, qred = [], [], []
    with open(bpfile, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                m = re.match(r'#\s*spinor\s*=\s*(\d+)', s)
                if m:
                    spinor = bool(int(m.group(1)))
                m = re.match(r'#\s*nv\s*=\s*(\d+)', s)
                if m:
                    nv = int(m.group(1))
                m = re.match(r'#\s*nb\s*=\s*(\d+)', s)
                if m:
                    nb = int(m.group(1))
                m = re.match(r'#\s*nodes:\s*(.*)', s)
                if m:
                    toks = m.group(1).split()
                    nodes = [(toks[i], float(toks[i + 1]))
                             for i in range(0, len(toks) - 1, 2)]
                continue
            parts = s.split()
            if nb is None or len(parts) < 5 + nb:
                continue
            try:
                dist.append(float(parts[1]))
                qred.append([float(parts[2]), float(parts[3]), float(parts[4])])
                eig.append([float(x) for x in parts[5:5 + nb]])
            except ValueError:
                continue
    if not eig:
        raise ValueError(f"No band data parsed from {bpfile}")
    return (np.array(dist), np.array(eig), nv, spinor, nodes, np.array(qred))


def plot_bandpath(bpfile, output_dir, energy_range_ev=(-6, 12), dpi=150):
    """
    Plot the UNFOLDED primitive-cell band structure (and, for spinor data,
    the Dresselhaus spin splitting of the levels around the gap, in meV).
    Unlike the folded MP-grid plot, conduction bands here are NOT overlaid
    4-fold -- CB1/CB2/CB3 are individually resolved.
    """
    print(f"Processing unfolded band path: {bpfile.name}")
    dist, eig_ha, nv, spinor, nodes, _q = _load_bandpath(bpfile)
    nk, nb = eig_ha.shape
    print(f"  {nk} path k-points, {nb} bands, nv = {nv}, spinor = {spinor}")

    vbm_ha = np.nanmax(eig_ha[:, nv - 1]) if nv else np.nanmax(eig_ha)
    eig_ev = (eig_ha - vbm_ha) * HA_TO_EV

    def _decorate(ax):
        for lbl, d in nodes:
            ax.axvline(d, color='#888888', linestyle='--', lw=0.7)
        ax.set_xticks([d for _, d in nodes])
        ax.set_xticklabels([r'$\Gamma$' if l == 'Gamma' else f'${l}$'
                            for l, _ in nodes], fontsize=12)
        ax.set_xlim(dist[0], dist[-1])

    # --- band structure ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 6))
    for b in range(nb):
        ax.plot(dist, eig_ev[:, b], 'k-', lw=0.9, alpha=0.7)
    ax.axhline(0.0, color='tab:red', linestyle='-', lw=0.8, alpha=0.7,
               label='VBM = 0')
    _decorate(ax)
    ax.set_ylabel('Energy (eV)', fontsize=11)
    ax.set_ylim(*energy_range_ev)
    ax.set_title(f'Unfolded primitive-cell bands — {bpfile.stem}', fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout()
    out = output_dir / f'bandpath_{bpfile.stem}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")

    # --- Dresselhaus spin splitting of levels around the gap (spinor) ------
    if not spinor or nv is None or nv % 2 != 0:
        return
    n_lvl_v = min(3, nv // 2)
    n_lvl_c = min(3, (nb - nv) // 2)
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for j in range(n_lvl_c):       # CB1, CB2, ...
        lvl = nv // 2 + 1 + j      # 1-based level index
        d_mev = (eig_ha[:, 2 * lvl - 1] - eig_ha[:, 2 * lvl - 2]) * HA_TO_EV * 1e3
        axes[0].plot(dist, d_mev, lw=1.1, label=f'CB{j + 1}')
    for j in range(n_lvl_v):       # VB1 = topmost valence level, ...
        lvl = nv // 2 - j
        d_mev = (eig_ha[:, 2 * lvl - 1] - eig_ha[:, 2 * lvl - 2]) * HA_TO_EV * 1e3
        axes[1].plot(dist, d_mev, lw=1.1, label=f'VB{j + 1}')
    axes[0].set_title(f'Spin splitting of bands — {bpfile.stem}', fontsize=11)
    for ax, tag in zip(axes, ('CB', 'VB')):
        _decorate(ax)
        ax.set_ylabel(rf'$\Delta_j(k)$ {tag} (meV)', fontsize=10)
        ax.legend(fontsize=9, loc='upper right')
        ax.set_ylim(bottom=0)
    fig.tight_layout()
    out = output_dir / f'bandpath_spin_splitting_{bpfile.stem}.png'
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
    parser.add_argument('--band-path', nargs='+', default=None,
                        metavar='PT',
                        help='High-symmetry path for band structure plot '
                             '(default: the --lattice default path).')
    parser.add_argument('--energy-range', nargs=2, type=float,
                        default=[-6.0, 12.0], metavar=('EMIN', 'EMAX'),
                        help='Energy window for band structure plot (eV)')
    parser.add_argument('--spin-sum', choices=['auto', 'on', 'off'],
                        default='auto',
                        help='Spinor (spin-orbit split) eigen files: merge '
                             'adjacent spin sub-bands into levels (occupations '
                             'summed, energy = pair mean). "auto" detects '
                             'spinor input from occupations <= 1 per band.')
    parser.add_argument('--band-vbm', type=int, default=None, metavar='IDX',
                        help='1-based band index taken as the valence-band '
                             'maximum (energy zero) when plotting band.dat from '
                             'a dft_band run. Default: nb//2 (half filling).')
    parser.add_argument('--lattice', choices=['fcc', 'wurtzite'], default='fcc',
                        help='Lattice/symmetry-point set for the band-structure '
                             'plot: "fcc" (cubic GaAs/Si, default) or "wurtzite" '
                             '(orthorhombic CdS supercell: Gamma/X/Y/A/S high-'
                             'symmetry points, default path A-Gamma-X-S-Y-Gamma).')

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
    parser.add_argument('--subtract-baseline', action='store_true',
                        help='For nex_k population plots, subtract the first '
                             'non-zero frame (the finite-A gauge/projection '
                             'offset) from every frame to isolate the genuine '
                             'field-driven excitation. Writes *_db_* files.')
    parser.add_argument('--spectral', action='store_true',
                        help='Also write an A(kx,E)-style spectral map of the '
                             'unfolded CB1 population (needs the 7-column '
                             '*_sbe_nex_k_unfold.data with the e_cb1 column).')
    parser.add_argument('--no-conductivity', action='store_true',
                        help='Skip the optical-conductivity sigma(w)=J/E plots '
                             '(global Re/Im spectrum + STFT Re-sigma map) from '
                             '*_sbe_rt.data.')
    parser.add_argument('--fmax-thz', type=float, default=4.0, metavar='F',
                        help='Upper frequency [THz] for the conductivity plots.')
    parser.add_argument('--stft-window-fs', type=float, default=None, metavar='T',
                        help='STFT window length [fs] for the Re-sigma(w,t) map. '
                             'Default: a quarter of the trace. Use a ps-scale '
                             'window to resolve the 0-4 THz band.')
    parser.add_argument('--stft-hop', type=int, default=1, metavar='H',
                        help='STFT hop in samples (1 = maximal N-1 overlap). The '
                             'effective hop is raised automatically if a 1-sample '
                             'hop would render too many columns.')
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

        # Optical conductivity sigma(w)=J(w)/E(w) from the current/field columns
        if not args.no_conductivity:
            for f in sorted(input_dir.glob('*_sbe_rt.data')):
                found_any = True
                try:
                    plot_conductivity(f, output_dir, fmax_thz=args.fmax_thz,
                                      dpi=args.dpi)
                    plot_conductivity_stft(f, output_dir, fmax_thz=args.fmax_thz,
                                           window_fs=args.stft_window_fs,
                                           hop=args.stft_hop, dpi=args.dpi)
                except Exception as exc:
                    print(f"  ERROR in conductivity for {f.name}: {exc}")

        for f in sorted(input_dir.glob('*_sbe_nex_k.data')):
            found_any = True
            plot_nex_k(f, output_dir, dpi=args.dpi,
                       log_scale=args.log_cmap, snapshots=args.snapshots)
            if args.subtract_baseline:
                plot_nex_k(f, output_dir, dpi=args.dpi,
                           log_scale=args.log_cmap, snapshots=args.snapshots,
                           subtract_baseline=True)

        # Physical (unfolded) CB1 populations on the primitive BZ
        for f in sorted(input_dir.glob('*_sbe_nex_k_unfold.data')):
            found_any = True
            plot_nex_k(f, output_dir, dpi=args.dpi,
                       log_scale=args.log_cmap, snapshots=args.snapshots,
                       unfold=True)
            if args.subtract_baseline:
                plot_nex_k(f, output_dir, dpi=args.dpi,
                           log_scale=args.log_cmap, snapshots=args.snapshots,
                           unfold=True, subtract_baseline=True)
            if args.spectral:
                stem = f.name[:-len('_sbe_nex_k_unfold.data')]
                bpfile = f.parent / f'{stem}_bandpath.data'
                if bpfile.exists():
                    plot_unfold_spectral(f, bpfile, output_dir, dpi=args.dpi)
                else:
                    print(f"  (skip spectral) {bpfile.name} not found "
                          f"(generate it with: epm_gaas_reference.py bandpath)")

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
                    dpi=args.dpi, spin_sum=args.spin_sum,
                    lattice=args.lattice)
            except Exception as exc:
                print(f"  ERROR in band structure for {kf.name}: {exc}")

        # dft_band output (band.dat)
        for bf in sorted(input_dir.glob('band.dat')):
            found_any = True
            try:
                plot_band_dat(
                    bf, output_dir,
                    energy_range_ev=tuple(args.energy_range),
                    dpi=args.dpi, vbm_index=args.band_vbm)
            except Exception as exc:
                print(f"  ERROR in dft_band plot for {bf.name}: {exc}")

        # Unfolded primitive-cell band path (epm_gaas_reference.py bandpath)
        for bf in sorted(input_dir.glob('*_bandpath.data')):
            found_any = True
            try:
                plot_bandpath(
                    bf, output_dir,
                    energy_range_ev=tuple(args.energy_range),
                    dpi=args.dpi)
            except Exception as exc:
                print(f"  ERROR in bandpath plot for {bf.name}: {exc}")

    if not found_any:
        print(f"No data files found in {input_dir.resolve()}")
        print("Expected: *_sbe_rt.data, *_sbe_rt_energy.data, *_sbe_nex.data, "
              "*_sbe_nex_k.data, *_k.data + *_eigen.data, band.dat")
        return

    print(f"\nDone.  Output: {output_dir.resolve()}")


if __name__ == '__main__':
    main()
