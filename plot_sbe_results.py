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
  *_sbe_nex.data         : excited electron count vs time (FULL, incl. dressing)
  *_sbe_nex_nonad.data   : NON-ADIABATIC (real) excited density -- dressed-frame
                           conduction population; drops the reversible A^2(t)
                           virtual dressing. A full-vs-non-adiabatic overlay
                           (twin Y-axes + shared log) is emitted automatically.
  *_sbe_nex_k_real.data  : per-k REAL-carrier LCB population (fixed-basis
                             diabatic n_ex, no A^2 breathing) -- the default
                             carrier map: snapshots + time-k maps. With
                             *_unfold_real.data, also the primitive-BZ map.
  *_sbe_nex_k.data       : per-k instantaneous Houston-basis LCB population
                             (carries the reversible virtual breathing; plotted
                             only with --instantaneous or if no _real file)
  *_sbe_intra_current.data: intra-band (Houston) current vs the total current
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
import os
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
    'Gamma': [0.0,       0.0,       0.0],
    'A':     [0.0,       0.0,       0.5],   # hexagonal A: Gamma-A along c* (the c-axis)
    'M':     [0.5,       0.0,       0.0],   # hexagonal M: b1/2 (in-plane zone-edge midpoint)
    'K':     [1.0 / 3.0, 1.0 / 3.0, 0.0],   # hexagonal K: Dirac-type corner (b1+b2)/3
    'H':     [1.0 / 3.0, 1.0 / 3.0, 0.5],   # K + c* (zone-top corner)
    'L':     [0.5,       0.0,       0.5],   # M + c* (zone-top edge)
    # orthorhombic-supercell aliases kept for back-compat with older CdS datasets:
    'X':     [0.5,       0.0,       0.0],
    'Y':     [0.0,       0.5,       0.0],
    'S':     [0.5,       0.5,       0.0],
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
# Hexagonal wurtzite path: matches the EPM-emitted *_bandpath.data nodes
# (A-Gamma-M-K-Gamma) so the with- and without-population CdS plots agree, and
# the Gamma-A c-axis segment leads.
DEFAULT_BAND_PATH_WZ  = ['A', 'Gamma', 'M', 'K', 'Gamma']
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


def _detect_lattice(kfile):
    """Auto-pick the band-path symmetry set from the '# material = X' header of a
    _k.data file: wurtzite for CdS (hexagonal -> the Gamma-A c-axis segment is in
    the default path), fcc otherwise (cubic GaAs/Si). 'fcc' if unreadable."""
    try:
        with open(kfile) as f:
            for line in f:
                if not line.startswith('#'):
                    break
                m = re.search(r'#\s*material\s*=\s*([A-Za-z0-9_]+)', line)
                if m:
                    mat = m.group(1).lower()
                    return 'wurtzite' if ('cds' in mat or 'wurtz' in mat) else 'fcc'
    except OSError:
        pass
    return 'fcc'


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
# C1: per-channel dissipation ledger (*_sbe_channels.data)
# ---------------------------------------------------------------------------
def plot_channels(filepath, output_dir, dpi=150):
    """Cumulative per-channel ledger: which dissipation channel did what.
    Two panels: conduction-population change dN (pairs created > 0 /
    recombined < 0) and energy change dE [Ha], per ring channel
    (e-ph, II, ring Auger, 2D Rana), vs time."""
    d = np.loadtxt(filepath, comments='#')
    if d.ndim == 1:
        d = d.reshape(1, -1)
    if d.shape[0] < 2:
        return
    t = d[:, 0] * AU_TIME_FS
    names  = ['e-ph', 'impact ionization', 'ring Auger', '2D Rana']
    colors = ['#7a5195', '#ef5675', '#003f5c', '#ffa600']
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 7.5), sharex=True)
    for i, (nm, c) in enumerate(zip(names, colors)):
        dn, de = d[:, 1 + 2*i], d[:, 2 + 2*i]
        if np.max(np.abs(dn)) == 0 and np.max(np.abs(de)) == 0:
            continue                       # channel off -> skip the flat line
        a1.plot(t, dn, lw=2, color=c, label=nm)
        a2.plot(t, de, lw=2, color=c, label=nm)
    a1.axhline(0, color='k', lw=.6)
    a2.axhline(0, color='k', lw=.6)
    a1.set_ylabel('cumulative dN$_{CB}$ [e$^-$/cell]')
    a2.set_ylabel('cumulative dE [Ha/cell]')
    a2.set_xlabel('t [fs]')
    a1.set_title('per-channel dissipation ledger (C1): pair creation (+) / recombination (-)')
    a1.grid(alpha=.3); a2.grid(alpha=.3); a1.legend(fontsize=9)
    out = output_dir / 'sbe_channels_ledger.png'
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


# ---------------------------------------------------------------------------
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


def plot_nex_comparison(nex_file, output_dir, dpi=150):
    """Overlay the FULL excited density (_sbe_nex.data, which includes the
    reversible A^2(t) virtual 'dressing') against the two NON-ADIABATIC real
    densities from _sbe_nex_nonad.data: col 2 = dressed-conduction projection,
    col 3 = the Option-A dressed-reference (delta0-subtracted, clamped) density
    the ring dissipators actually see. Two panels: linear with twin Y-axes (the
    curves live on very different scales during the pulse), and a shared log
    axis. The shaded gap is the reversible dressing that is NOT real carriers."""
    nonad_file = nex_file.parent / nex_file.name.replace('_sbe_nex.data',
                                                         '_sbe_nex_nonad.data')
    if not nonad_file.exists():
        return

    try:
        f = np.loadtxt(nex_file, comments='#')
        d = np.loadtxt(nonad_file, comments='#')
        if f.ndim == 1:
            f = f.reshape(1, -1)
        if d.ndim == 1:
            d = d.reshape(1, -1)
        t, full = f[:, 0], f[:, 1]
        tn, proj = d[:, 0], d[:, 1]
        dref = d[:, 2] if d.shape[1] > 2 else d[:, 1]
    except Exception as exc:
        print(f"  (skip) nex comparison: {exc}")
        return
    n = min(len(t), len(tn))
    if n < 2:
        return
    t, full, proj, dref = t[:n], full[:n], proj[:n], dref[:n]

    C_FULL, C_PROJ, C_DREF = '#1f77b4', '#d62728', '#2ca02c'
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 8))

    # -- top: linear, twin Y-axes (full on the left, non-adiabatic on the right) --
    l1, = a1.plot(t, full, color=C_FULL, lw=1.2,
                  label=r'full n$_{ex}$ (incl. dressing)')
    fill = a1.fill_between(t, dref, full, color=C_FULL, alpha=0.10,
                           label='reversible dressing (not real)')
    a1.set_ylabel(r'full n$_{ex}$ [cm$^{-3}$]', color=C_FULL)
    a1.tick_params(axis='y', labelcolor=C_FULL)
    a1b = a1.twinx()
    l2, = a1b.plot(t, proj, color=C_PROJ, lw=1.3,
                   label='non-adiabatic: dressed-conduction')
    l3, = a1b.plot(t, dref, color=C_DREF, lw=1.6,
                   label='non-adiabatic: Option-A (ring sees this)')
    a1b.set_ylabel(r'non-adiabatic n$_{ex}$ [cm$^{-3}$]')
    a1.set_xlabel('time [fs]')
    a1.set_title('Excited density: full vs non-adiabatic (twin Y-axes)')
    a1.grid(True, alpha=0.3, ls='--')
    a1.legend(handles=[l1, fill, l2, l3], loc='upper left', fontsize=8)

    # -- bottom: shared log axis --
    a2.semilogy(t, np.clip(full, 1e-30, None), color=C_FULL, lw=1.2,
                label='full (incl. dressing)')
    a2.semilogy(t, np.clip(proj, 1e-30, None), color=C_PROJ, lw=1.3,
                label='non-adiabatic: dressed-conduction')
    a2.semilogy(t, np.clip(dref, 1e-30, None), color=C_DREF, lw=1.6,
                label='non-adiabatic: Option-A (ring sees this)')
    a2.set_xlabel('time [fs]')
    a2.set_ylabel(r'n$_{ex}$ [cm$^{-3}$]  (log)')
    a2.set_title('Same, shared log scale')
    a2.grid(True, which='both', alpha=0.3, ls='--')
    a2.legend(loc='lower right', fontsize=8)
    ip = int(np.argmax(full))
    ratio = full[ip] / max(dref[ip], 1e-30)
    a2.annotate(f'peak full / Option-A = {ratio:.0f}x',
                xy=(t[ip], max(full[ip], 1e-30)), fontsize=9, ha='center',
                va='bottom', color='0.2')

    fig.tight_layout()
    out = output_dir / (nex_file.name.replace('.data', '') + '_vs_nonad.png')
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


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


def _fold_to_cubic(kpoints_prim):
    """Recover the supercell k (k_sc = k_prim - G0(isub)) and wrap it into
    [-0.5, 0.5). The coset offset G0 is a SUPERCELL reciprocal-lattice vector,
    hence an integer triplet in sc-reduced units, so subtracting it then wrapping
    is identical to wrapping k_prim directly: k_sc = k_prim - round(k_prim). This
    needs neither the offset table (the *_sbe_nex_k_unfold.data output does not
    carry one) nor the coset count, so it folds 2-coset (CdS/graphene), 4-coset
    (cubic FCC) and any N-coset map alike. Summing the populations of the cosets
    sharing a k_sc collapses the folded valleys back onto the single supercell
    zone -- the clean per-k total of the lowest conduction band."""
    return kpoints_prim - np.round(kpoints_prim)


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


def _read_bmatrix(kdata_path):
    """Read the reciprocal vectors (# b1/# b2/# b3 [a.u.]) from a _k.data header,
    written by the non-orthogonal EPM. Returns rows b1,b2,b3 (3x3) or None when
    absent (orthogonal/legacy cubic dataset -> the reduced coords are already
    axis-aligned and the standard heatmap is used)."""
    try:
        rows = {}
        with open(kdata_path) as f:
            for line in f:
                s = line.strip()
                if not s.startswith('#'):
                    break
                m = re.match(r'#\s*b([123])\s*=\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)', s)
                if m:
                    rows[int(m.group(1))] = [float(m.group(2)), float(m.group(3)), float(m.group(4))]
        if len(rows) == 3:
            return np.array([rows[1], rows[2], rows[3]])
    except (OSError, ValueError):
        pass
    return None


def _bmatrix_for(nex_file, suffix):
    """Given a *_sbe_nex_k*.data file, return the reciprocal b_matrix from the
    sibling ground-state {stem}_k.data, or None when it is absent / orthogonal.
    Only the non-orthogonal (primitive) EPM writes the b1/b2/b3 header, so a
    legacy cubic dataset transparently keeps the standard reduced-axis heatmap."""
    stem = nex_file.name[:-len(suffix)] if nex_file.name.endswith(suffix) else nex_file.stem
    kdata = nex_file.parent / f'{stem}_k.data'
    return _read_bmatrix(kdata) if kdata.exists() else None


def _cartesian_bz_grid(kfrac, pop, b_matrix, nbin):
    """Un-shear a triclinic k-grid into a regular CARTESIAN heatmap volume.
    kfrac (nk,3) reduced -> Cartesian k = kfrac @ b_matrix, wrapped into the
    Wigner-Seitz BZ (nearest reciprocal-lattice vector), then averaged into an
    nbin^3 regular Cartesian grid (NaN where no k falls -- the BZ corners). For
    a dense input grid this gives a smooth heatmap of the true Brillouin zone.
    Returns kx_u, ky_u, kz_u (bin centres, a.u.) and pop3d."""
    kc = (kfrac - np.round(kfrac)) @ b_matrix          # into the reciprocal cell
    # Wigner-Seitz wrap: subtract the nearest reciprocal lattice vector
    ijk = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    G = ijk @ b_matrix                                  # (27,3) candidate G's
    d2 = ((kc[:, None, :] - G[None, :, :]) ** 2).sum(axis=2)
    kc = kc - G[np.argmin(d2, axis=1)]
    kmax = np.abs(kc).max() * 1.0001
    edges = np.linspace(-kmax, kmax, nbin + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    ix = np.clip(np.digitize(kc[:, 0], edges) - 1, 0, nbin - 1)
    iy = np.clip(np.digitize(kc[:, 1], edges) - 1, 0, nbin - 1)
    iz = np.clip(np.digitize(kc[:, 2], edges) - 1, 0, nbin - 1)
    ssum = np.zeros((nbin, nbin, nbin)); cnt = np.zeros((nbin, nbin, nbin))
    np.add.at(ssum, (ix, iy, iz), pop)
    np.add.at(cnt, (ix, iy, iz), 1.0)
    pop3d = np.where(cnt > 0, ssum / np.maximum(cnt, 1.0), np.nan)
    return cen, cen, cen, pop3d


def _make_norm(vmin, vmax, log_scale):
    """Return a matplotlib Normalize or LogNorm for colormap scaling."""
    if log_scale and vmax > 0:
        # Floor at 1e-6 × peak so zeros don't break the log scale
        floor = max(vmax * 1e-6, 1e-30)
        return mcolors.LogNorm(vmin=max(vmin, floor), vmax=vmax)
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def _heatmap_ax(ax, k_a, k_b, grid2d, label_a, label_b, title,
                vmin=None, vmax=None, factor=8, log_scale=False, unit='reduced',
                clip_poly=None):
    if grid2d.size == 0 or np.all(np.isnan(grid2d)):
        ax.set_title(title + " (no data)")
        return None
    ka_f, kb_f, gf = _interp2d(grid2d, k_a, k_b, factor=factor)
    if clip_poly is not None:
        # mask interpolated pixels OUTSIDE the BZ silhouette so the map shows the
        # true zone shape (and the fade-out) instead of a filled square.
        from matplotlib.path import Path as _Path
        GA, GB = np.meshgrid(ka_f, kb_f, indexing='ij')
        inside = _Path(clip_poly).contains_points(
            np.column_stack([GA.ravel(), GB.ravel()])).reshape(GA.shape)
        gf = np.where(inside, gf, np.nan)
    cmap = plt.get_cmap(CMAP_POP).copy(); cmap.set_bad(alpha=0.0)
    norm = _make_norm(vmin if vmin is not None else np.nanmin(gf),
                      vmax if vmax is not None else np.nanmax(gf),
                      log_scale)
    im = ax.imshow(
        gf.T, origin='lower', aspect='equal' if unit != 'reduced' else 'auto',
        extent=[ka_f[0], ka_f[-1], kb_f[0], kb_f[-1]],
        cmap=cmap, norm=norm, interpolation='nearest')
    plt.colorbar(im, ax=ax, shrink=0.8)
    ax.set_xlabel(f'{label_a} [{unit}]')
    ax.set_ylabel(f'{label_b} [{unit}]')
    ax.set_title(title)
    return im


def _cubic_valleys(b_matrix, delta_frac=0.85):
    """High-symmetry markers for an FCC/diamond cell in Cartesian a.u.: Gamma,
    the six X points (2*pi/a along the cubic axes; 2*pi/a = |b1|/sqrt3 for FCC),
    and the six Delta-valley minima (delta_frac*X, Si CBM at ~0.85*X). Returns
    a list of (kx,ky,kz, marker, color, ms, label) overlay points."""
    kX = float(np.linalg.norm(b_matrix[0]) / np.sqrt(3.0))
    ov = [(0.0, 0.0, 0.0, '+', 'lime', 12, 'Γ')]
    for ax in range(3):
        for s in (+1.0, -1.0):
            X = [0.0, 0.0, 0.0]; X[ax] = s * kX
            D = [0.0, 0.0, 0.0]; D[ax] = s * delta_frac * kX
            ov.append((X[0], X[1], X[2], 'x', 'cyan', 9, 'X'))
            ov.append((D[0], D[1], D[2], 'o', 'red', 6, 'Δ'))
    return ov, kX


def _bz_outline_2d(b_matrix, ia, ib, nsamp=6000):
    """Closed polygon = the silhouette (projection onto Cartesian axes ia,ib) of
    the Brillouin zone, i.e. the Wigner-Seitz cell of the reciprocal lattice
    `b_matrix`. General: FCC -> truncated octahedron (square-with-cut-corners in
    kx-ky), hexagonal -> hexagon, etc. So the map shows the TRUE zone boundary,
    not just the square plot frame. Returns None if scipy is unavailable."""
    rng = (-2, -1, 0, 1, 2)
    G = np.array([i*b_matrix[0] + j*b_matrix[1] + k*b_matrix[2]
                  for i in rng for j in rng for k in rng if (i, j, k) != (0, 0, 0)])
    half = 0.5 * (G * G).sum(1)                         # Bragg planes k.G = |G|^2/2
    s = np.arange(nsamp) + 0.5                          # Fibonacci-sphere directions
    phi = np.arccos(1.0 - 2.0*s/nsamp); th = np.pi*(1.0 + 5.0**0.5)*s
    u = np.column_stack([np.sin(phi)*np.cos(th), np.sin(phi)*np.sin(th), np.cos(phi)])
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(u @ G.T > 1e-9, half[None, :] / (u @ G.T), np.inf).min(axis=1)
    pts = (u * r[:, None])[:, [ia, ib]]                 # BZ-surface points, projected
    try:
        from scipy.spatial import ConvexHull
        poly = pts[ConvexHull(pts).vertices]
    except Exception:
        return None
    return np.vstack([poly, poly[:1]])


def _bz_wireframe_3d(b_matrix):
    """Edge list [(p0,p1), ...] of the first Brillouin zone (the Wigner-Seitz
    cell of the reciprocal lattice `b_matrix`), via the Voronoi cell of the
    origin: FCC -> truncated octahedron, hexagonal -> hexagonal prism. Each
    ridge polygon shared with the origin contributes its edges. Returns None
    if scipy is unavailable."""
    try:
        from scipy.spatial import Voronoi
    except Exception:
        return None
    rng = (-1, 0, 1)
    pts = np.array([i*b_matrix[0] + j*b_matrix[1] + k*b_matrix[2]
                    for i in rng for j in rng for k in rng])
    i0 = int(np.argmin((pts * pts).sum(1)))            # the origin
    vor = Voronoi(pts)
    edges, seen = [], set()
    for (p1, p2), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if i0 not in (p1, p2) or -1 in rv:
            continue
        poly = vor.vertices[rv]
        for a in range(len(rv)):
            i, j = rv[a], rv[(a + 1) % len(rv)]
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            edges.append((poly[a], poly[(a + 1) % len(rv)]))
    return edges


def _save_bz3d(kfrac, pop, b_matrix, t_val, t_unit, output_dir, dpi,
               tag='nex_k', basis_label='real-carrier', vmax=None,
               valleys=None, pop_floor=0.02):
    """Paper-style 3D Brillouin-zone population plot: the BZ wireframe (Wigner-
    Seitz cell) + the Monkhorst-Pack k-points as a scatter coloured AND sized by
    the population. Weakly-populated points fade out (alpha, size -> small) and
    below pop_floor*max they are dropped entirely, so only the populated valleys
    show -- the 3D analogue of the *_cart_snap_* maps."""
    edges = _bz_wireframe_3d(b_matrix)
    kc = (kfrac - np.round(kfrac)) @ b_matrix           # un-shear
    ijk = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    G = ijk @ b_matrix
    d2 = ((kc[:, None, :] - G[None, :, :]) ** 2).sum(axis=2)
    kc = kc - G[np.argmin(d2, axis=1)]                  # Wigner-Seitz wrap

    vmax = float(vmax if vmax is not None else max(np.nanmax(pop), 1e-30))
    m = np.clip(np.nan_to_num(pop) / vmax, 0.0, 1.0)
    sel = m > pop_floor                                  # transparent below the floor

    fig = plt.figure(figsize=(8.5, 8))
    ax = fig.add_subplot(projection='3d')
    if edges is not None:
        for p0, p1 in edges:
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]],
                    color='cornflowerblue', lw=0.9, alpha=0.9)
    if np.any(sel):
        sc = ax.scatter(kc[sel, 0], kc[sel, 1], kc[sel, 2],
                        s=4.0 + 120.0 * m[sel] ** 1.5, c=pop[sel],
                        cmap=CMAP_POP, vmin=0.0, vmax=vmax,
                        alpha=None, depthshade=True)
        sc.set_alpha(np.clip(0.15 + 0.85 * m[sel], 0.0, 1.0))
        fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.08, label='population')
    if valleys is not None:
        for p in valleys:
            ax.plot([p[0]], [p[1]], [p[2]], marker=p[3], color=p[4],
                    ms=p[5], mew=2, ls='none')
    lim = max(np.linalg.norm(b_matrix[i]) for i in range(3)) * 0.62
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('kx [a.u.]'); ax.set_ylabel('ky [a.u.]'); ax.set_zlabel('kz [a.u.]')
    ax.set_title(f'{basis_label} conduction population in the 3D BZ,  t = {t_val:.3f} {t_unit}\n'
                 f'(point size & opacity ∝ population; < {pop_floor:.0%} of max hidden)')
    safe_t = f'{t_val:.6f}'.replace('-', 'm').replace('+', 'p')
    out = output_dir / f'{tag}_bz3d_t{safe_t}{t_unit}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def _cube_faces(x0, x1, y0, y1, z0, z1):
    """The 6 quad faces of an axis-aligned box, as lists of 4 (x,y,z) corners."""
    p = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    return [[p[0], p[1], p[2], p[3]], [p[4], p[5], p[6], p[7]],   # z0, z1
            [p[0], p[1], p[5], p[4]], [p[3], p[2], p[6], p[7]],   # y0, y1
            [p[1], p[2], p[6], p[5]], [p[0], p[3], p[7], p[4]]]   # x1, x0


def _save_bz3d_voxel(cpop3d, cx, cy, cz, b_matrix, t_val, t_unit, output_dir,
                     dpi, tag='nex_k', basis_label='real-carrier', vmax=None,
                     valleys=None, pop_floor=0.05, smooth_sigma=0.0,
                     gap=0.05, alpha_gamma=2.5):
    """Variant (b) of the 3D BZ population plot: a SMOOTHED SEMI-TRANSPARENT
    VOXEL CLOUD rendered from the un-sheared Cartesian grid `cpop3d` (the same
    array behind the *_cart_snap_* maps), instead of the MP-point scatter of
    _save_bz3d. Gaussian-smoothed (scipy), voxel opacity ∝ population (weak
    voxels below pop_floor*max are fully transparent), BZ wireframe overlaid.
    The volumetric analogue of the paper-style scatter -- shows the SHAPE of
    the populated valleys rather than the sampled points."""
    pop_raw = np.nan_to_num(np.asarray(cpop3d, dtype=float))
    # Colour scale from the TRUE (unsmoothed) peak population, so the colorbar
    # is honest even if we smooth: on a coarse grid a sharp Gamma spike would
    # otherwise read many times too low (Gaussian smearing rescales the peak).
    vmax = float(vmax if vmax is not None else max(pop_raw.max(), 1e-30))
    pop = pop_raw
    if smooth_sigma > 0:
        try:
            from scipy.ndimage import gaussian_filter
            pop = gaussian_filter(pop_raw, sigma=smooth_sigma)
        except Exception:
            pop = pop_raw                               # unsmoothed fallback
    m = np.clip(pop / vmax, 0.0, 1.0)
    filled = m > pop_floor
    if not np.any(filled):
        return                                           # nothing to draw yet

    # per-voxel RGBA: colormap colour, NON-LINEAR opacity in population.
    # (plt.get_cmap, not matplotlib.cm.get_cmap -- the latter was removed in mpl 3.9)
    # alpha = amin + (amax-amin) * t**gamma, t = (m-floor)/(1-floor) in [0,1].
    # gamma > 1 keeps faint cells nearly transparent so that many overlapping
    # weak voxels don't accumulate into an opaque fog, while the LOCALIZED
    # maxima (the interesting part) stay bright and pop out. gamma=1 -> linear.
    rgba = plt.get_cmap(CMAP_POP)(m)
    _amin, _amax = 0.035, 0.90
    _t = np.clip((m - pop_floor) / (1.0 - pop_floor), 0.0, 1.0)
    rgba[..., 3] = np.where(filled, _amin + (_amax - _amin) * _t ** max(alpha_gamma, 0.01), 0.0)

    # voxel corner edges per axis (cx/cy/cz are bin centres, uniform spacing)
    def _edges(c):
        d = c[1] - c[0] if len(c) > 1 else 1.0
        return np.concatenate([c - 0.5 * d, [c[-1] + 0.5 * d]]), abs(d)
    ex, dx = _edges(cx); ey, dy = _edges(cy); ez, dz = _edges(cz)

    fig = plt.figure(figsize=(8.5, 8))
    ax = fig.add_subplot(projection='3d')
    edges = _bz_wireframe_3d(b_matrix)
    if edges is not None:
        for p0, p1 in edges:
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]],
                    color='cornflowerblue', lw=0.9, alpha=0.9)

    # Build every filled cube as 6 explicit faces in ONE Poly3DCollection so
    # matplotlib depth-sorts ALL faces together (correct back-to-front alpha
    # compositing). ax.voxels() instead makes a separate collection per cube and
    # cannot composite transparency across them -- you'd look straight through
    # the outer shell without the inner cubes ever showing. Each cube is shrunk
    # by `gap` so neighbours don't touch and the interior stays visible.
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    gap = float(np.clip(gap, 0.0, 0.9))
    faces, fcolors = [], []
    for i, j, k in zip(*np.where(filled)):
        x0, x1 = ex[i] + 0.5*gap*dx, ex[i+1] - 0.5*gap*dx
        y0, y1 = ey[j] + 0.5*gap*dy, ey[j+1] - 0.5*gap*dy
        z0, z1 = ez[k] + 0.5*gap*dz, ez[k+1] - 0.5*gap*dz
        faces += _cube_faces(x0, x1, y0, y1, z0, z1)
        fcolors += [rgba[i, j, k]] * 6
    pc = Poly3DCollection(faces, facecolors=fcolors, edgecolors='none')
    ax.add_collection3d(pc)
    if valleys is not None:
        for p in valleys:
            ax.plot([p[0]], [p[1]], [p[2]], marker=p[3], color=p[4],
                    ms=p[5], mew=2, ls='none')
    sm = plt.cm.ScalarMappable(cmap=CMAP_POP,
                               norm=plt.Normalize(vmin=0.0, vmax=vmax))
    fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.08, label='population')
    lim = max(np.linalg.norm(b_matrix[i]) for i in range(3)) * 0.62
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('kx [a.u.]'); ax.set_ylabel('ky [a.u.]'); ax.set_zlabel('kz [a.u.]')
    _smtxt = f'Gaussian σ={smooth_sigma:g} bins' if smooth_sigma > 0 else 'raw k-bins (unsmoothed)'
    ax.set_title(f'{basis_label} conduction population voxel cloud,  t = {t_val:.3f} {t_unit}\n'
                 f'(opacity ∝ population^{alpha_gamma:g}, {_smtxt}; '
                 f'< {pop_floor:.0%} of max transparent)')
    safe_t = f'{t_val:.6f}'.replace('-', 'm').replace('+', 'p')
    out = output_dir / f'{tag}_bz3dvox_t{safe_t}{t_unit}.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def _overlay_valleys(ax, ia, ib, valleys, b_matrix=None):
    """Project the 3D valley markers onto panel axes (ia,ib) and plot them, and
    (if b_matrix given) draw the true BZ boundary silhouette + size the axes to it."""
    seen = set()
    for p in valleys:
        a, b, mk, col, ms = p[ia], p[ib], p[3], p[4], p[5]
        key = (round(a, 6), round(b, 6), mk)
        if key in seen:
            continue
        seen.add(key)
        ax.plot(a, b, marker=mk, color=col, ms=ms, mew=2, ls='none', zorder=5)
    if b_matrix is not None:
        poly = _bz_outline_2d(b_matrix, ia, ib)
        if poly is not None:
            ax.plot(poly[:, 0], poly[:, 1], '-', color='magenta', lw=1.4, zorder=4)
            m = 1.08 * float(np.abs(poly).max())
            ax.set_xlim(-m, m); ax.set_ylim(-m, m)      # extend past the coarse-grid data


def _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi,
                   log_scale=False, tag='nex_k', basis_label='Houston-basis',
                   unit='reduced', valleys=None, b_matrix=None):
    vmin = np.nanmin(pop3d)
    vmax = max(np.nanmax(pop3d), vmin + 1e-30)

    # clip the Cartesian heatmaps to the BZ silhouette (so the map shows the true
    # zone shape + fade-out, not a filled square from the interpolation)
    cp = [None, None, None]
    if b_matrix is not None:
        cp = [_bz_outline_2d(b_matrix, 0, 1), _bz_outline_2d(b_matrix, 0, 2),
              _bz_outline_2d(b_matrix, 1, 2)]

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))
    _heatmap_ax(axes[0], kx_u, ky_u, _project(pop3d, 2),
                'kx', 'ky', 'pop_lcb: kx-ky (avg kz)',
                vmin=vmin, vmax=vmax, log_scale=log_scale, unit=unit, clip_poly=cp[0])
    _heatmap_ax(axes[1], kx_u, kz_u, _project(pop3d, 1),
                'kx', 'kz', 'pop_lcb: kx-kz (avg ky)',
                vmin=vmin, vmax=vmax, log_scale=log_scale, unit=unit, clip_poly=cp[1])
    _heatmap_ax(axes[2], ky_u, kz_u, _project(pop3d, 0),
                'ky', 'kz', 'pop_lcb: ky-kz (avg kx)',
                vmin=vmin, vmax=vmax, log_scale=log_scale, unit=unit, clip_poly=cp[2])
    if valleys is not None:
        _overlay_valleys(axes[0], 0, 1, valleys, b_matrix)   # kx-ky
        _overlay_valleys(axes[1], 0, 2, valleys, b_matrix)   # kx-kz
        _overlay_valleys(axes[2], 1, 2, valleys, b_matrix)   # ky-kz

    zone = 'Cartesian BZ (a.u.)' if unit != 'reduced' else 'reduced k'
    fig.suptitle(f'{basis_label} conduction population [{zone}],  t = {t_val:.6f} {t_unit}')
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


def _read_xy_data(filepath):
    """Read a whitespace SALMON .data file -> ndarray (skip '#' comment lines)."""
    rows = []
    for line in open(filepath):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            continue
    return np.array(rows) if rows else np.empty((0, 0))


def plot_intra_current(filepath, rt_filepath, output_dir, dpi=150):
    """Intra-band (Houston-basis) current J_intra(t), overlaid with the total
    gauge-invariant current J_tot(t) from *_sbe_rt.data when available. In the
    velocity gauge only J_tot is physical; J_intra is the physical drift part in
    the Houston basis and J_inter = J_tot - J_intra the interband polarization."""
    print(f"Processing {filepath.name} (intra-band Houston current) ...")
    a = _read_xy_data(filepath)
    if a.size == 0 or a.shape[1] < 4:
        print("  (skip) no data")
        return
    t = a[:, 0]
    Ji = a[:, 1:4]
    Jt = None
    if rt_filepath is not None:
        rt = _read_xy_data(rt_filepath)
        # *_sbe_rt.data total current Jm = columns 14:16 (1-indexed) -> 13:16
        if rt.size and rt.shape[1] >= 16:
            Jt = np.array([np.interp(t, rt[:, 0], rt[:, 13 + d]) for d in range(3)]).T

    labels = ['x', 'y', 'z']
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), sharex=True)
    for d in range(3):
        ax = axes[d]
        if Jt is not None:
            ax.plot(t, Jt[:, d], color='gray', lw=1.4, label='J_total (gauge-inv.)')
            ax.plot(t, Jt[:, d] - Ji[:, d], color='steelblue', lw=1.0,
                    alpha=0.8, label='J_inter (= tot - intra)')
        ax.plot(t, Ji[:, d], color='crimson', lw=1.6, label='J_intra (Houston)')
        ax.axhline(0, color='k', lw=0.5, ls=':')
        ax.set_xlabel('time [a.u.]')
        ax.set_ylabel(f'J_{labels[d]} [a.u.]')
        ax.set_title(f'current {labels[d]}')
        ax.legend(fontsize=8)
    fig.suptitle('Intra-band (Houston-basis) vs total current  '
                 '— intra = Boltzmann drift, vanishes when the field is off')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = output_dir / 'sbe_intra_current.png'
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved {out.name}")


def plot_nex_k(filepath, output_dir, dpi=150, log_scale=False, snapshots=False,
               unfold=False, subtract_baseline=False, real=False, b_matrix=None,
               mark_valleys=False, bz3d=False, bz3d_voxel=False, voxel_smooth=0.0,
               cb_sum=False, voxel_gap=0.05, voxel_gamma=2.5):
    print(f"Processing {filepath.name}  "
          f"(cmap={'log' if log_scale else 'linear'}, "
          f"snapshots={'on' if snapshots else 'off'}"
          f"{', REAL carriers' if real else ''}"
          f"{', unfolded primitive BZ' if unfold else ''}"
          f"{', baseline-subtracted' if subtract_baseline else ''}) ...")
    # real = fixed-basis diabatic (real carriers, no A^2 breathing);
    # else = instantaneous Houston-basis population.
    basis_label = 'real-carrier (diabatic)' if real else 'Houston-basis'
    suff = '_real' if real else ''
    tag = ('nex_k_unfold' if unfold else 'nex_k') + suff
    if subtract_baseline:
        tag += '_db'
    # The unfolded primitive-zone map legitimately shows the CB1 population of
    # every folded valley (each coset's offset G0), so it carries satellite peaks
    # at the zone boundary. The FOLDED view sums the cosets sharing a supercell k
    # (k_sc = k_prim - G0, an integer wrap) back onto the single supercell zone --
    # the per-k total of the lowest conduction band (any coset count).
    ftag = 'nex_k_fold' + suff + ('_db' if subtract_baseline else '')

    # Optional (--bz3d-cb-sum): use the SUM of BOTH recorded conduction bands
    # (CB1+CB2) as the per-k population instead of the lowest-CB-only column.
    # The two CB populations are in the sibling four-level real file
    # SYSNAME_sbe_nex_k_lev_real.data (cols pop_cb1, pop_cb2). Keyed by time and
    # substituted below; the k-order matches (both written in ik order).
    cbsum_by_t = None
    if cb_sum:
        levfile = filepath.parent / filepath.name.replace('_sbe_nex_k', '_sbe_nex_k_lev')
        if levfile.exists():
            cbsum_by_t = {round(lt, 6): pop4[:, 2] + pop4[:, 3]
                          for lt, _lu, _lk, pop4 in _iter_nex_k_lev_blocks(levfile)}
            tag += '_cbsum'
            basis_label += ', CB1+CB2'
            print(f"  --bz3d-cb-sum: population = CB1+CB2 from {levfile.name}")
        else:
            print(f"  (warn) --bz3d-cb-sum: {levfile.name} not found "
                  f"-- falling back to lowest-CB only")

    kx_u = ky_u = kz_u = ix = iy = iz = None
    pop3d = None
    fkx_u = fky_u = fkz_u = fix = fiy = fiz = None
    fpop3d = None
    pop_baseline = None       # per-k population of the first non-zero frame
    times, marg_kx, marg_ky, marg_kz = [], [], [], []
    fmarg_kx, fmarg_ky, fmarg_kz = [], [], []
    cx_u = cy_u = cz_u = None                       # un-sheared Cartesian axes (a.u.)
    cmarg_kx, cmarg_ky, cmarg_kz = [], [], []       # Cartesian-BZ k-t marginals
    t_unit_last = ''
    n_blocks = 0

    for t_val, t_unit, kpoints, pop, _levels, sub in _iter_nex_k_blocks(filepath, unfold=unfold):
        t_unit_last = t_unit
        if cbsum_by_t is not None:                       # --bz3d-cb-sum override
            cs = cbsum_by_t.get(round(t_val, 6))
            if cs is not None and len(cs) == len(pop):
                pop = cs
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
            # For a non-orthogonal (triclinic) grid the reduced axes are sheared,
            # so the per-axis heatmap is geometrically misleading. With the
            # reciprocal vectors (b_matrix) we un-shear into a regular CARTESIAN
            # Wigner-Seitz volume -- a true picture of the Brillouin zone that
            # gets smoother the denser the k-grid is.
            cnbin = max(len(kx_u), len(ky_u), len(kz_u))
        pop3d.fill(np.nan)
        pop3d[ix, iy, iz] = pop
        if snapshots:
            _save_snapshot(pop3d, kx_u, ky_u, kz_u, t_val, t_unit, output_dir, dpi,
                           log_scale=log_scale, tag=tag, basis_label=basis_label)
        # Un-shear into the true Cartesian BZ (when the reciprocal vectors are
        # known): used BOTH for the Cartesian snapshot AND the Cartesian k-t maps,
        # so both the k-k and the k-t views exist in reduced AND Cartesian coords.
        if b_matrix is not None:
            cx, cy, cz, cpop3d = _cartesian_bz_grid(kpoints, pop, b_matrix, cnbin)
            if cx_u is None:
                cx_u, cy_u, cz_u = cx, cy, cz
            cmarg_kx.append(np.nanmean(cpop3d, axis=(1, 2)))
            cmarg_ky.append(np.nanmean(cpop3d, axis=(0, 2)))
            cmarg_kz.append(np.nanmean(cpop3d, axis=(0, 1)))
            if snapshots:
                vlys = _cubic_valleys(b_matrix)[0] if mark_valleys else None
                _save_snapshot(cpop3d, cx, cy, cz, t_val, t_unit, output_dir, dpi,
                               log_scale=log_scale, tag=tag + '_cart',
                               basis_label=basis_label, unit='a.u.', valleys=vlys,
                               b_matrix=b_matrix if mark_valleys else None)
            if bz3d:
                _save_bz3d(kpoints, pop, b_matrix, t_val, t_unit, output_dir, dpi,
                           tag=tag, basis_label=basis_label,
                           valleys=_cubic_valleys(b_matrix)[0] if mark_valleys else None)
            if bz3d_voxel:
                _save_bz3d_voxel(cpop3d, cx, cy, cz, b_matrix, t_val, t_unit,
                                 output_dir, dpi, tag=tag, basis_label=basis_label,
                                 valleys=_cubic_valleys(b_matrix)[0] if mark_valleys else None,
                                 smooth_sigma=voxel_smooth, gap=voxel_gap,
                                 alpha_gamma=voxel_gamma)
        times.append(t_val)
        marg_kx.append(np.nanmean(pop3d, axis=(1, 2)))
        marg_ky.append(np.nanmean(pop3d, axis=(0, 2)))
        marg_kz.append(np.nanmean(pop3d, axis=(0, 1)))

        # Folded supercell-zone view: sum the cosets sharing each k_sc.
        if unfold and sub is not None:
            ksc = _fold_to_cubic(kpoints_prim)
            fpop = pop  # baseline (if any) already applied to `pop` above
            if fkx_u is None:
                fkx_u, fky_u, fkz_u, fix, fiy, fiz = _build_grid_info(ksc)
                fpop3d = np.zeros((len(fkx_u), len(fky_u), len(fkz_u)))
            fpop3d.fill(0.0)
            np.add.at(fpop3d, (fix, fiy, fiz), fpop)
            if snapshots:
                _save_snapshot(fpop3d, fkx_u, fky_u, fkz_u, t_val, t_unit, output_dir,
                               dpi, log_scale=log_scale, tag=ftag, basis_label=basis_label)
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

    # Cartesian (un-sheared, a.u.) k-t maps -- the same evolution against the
    # PHYSICAL k-axis, so it can be read together with the *_cart_snap_* frames.
    if cmarg_kx:
        print(f"  writing Cartesian-BZ time-k maps ...")
        _save_kt_map(times, t_unit_last, cx_u, 'kx', cmarg_kx, output_dir, dpi,
                     log_scale=log_scale, tag=tag + '_cart')
        _save_kt_map(times, t_unit_last, cy_u, 'ky', cmarg_ky, output_dir, dpi,
                     log_scale=log_scale, tag=tag + '_cart')
        _save_kt_map(times, t_unit_last, cz_u, 'kz', cmarg_kz, output_dir, dpi,
                     log_scale=log_scale, tag=tag + '_cart')

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


def _map_primitive_population(qred, gridk, gpop, spacing, tol=1.3):
    """Map the per-k LCB population (gpop on the MP grid gridk) onto each
    band-path point qred by nearest k in the native fractional convention
    (wrapped mod 1 -- works for any, incl. the non-orthogonal FCC, lattice).
    Points farther than tol*spacing from any grid k get NaN."""
    P = np.full(len(qred), np.nan)
    for i, q in enumerate(qred):
        d = gridk - q[None, :]
        d -= np.round(d)                       # wrap mod 1 (umklapp)
        d2 = (d ** 2).sum(axis=1)
        j = int(np.argmin(d2))
        if np.sqrt(d2[j]) <= tol * spacing:
            P[i] = gpop[j]
    return P


def _iter_nex_k_lev_blocks(filepath):
    """Stream SYSNAME_sbe_nex_k_lev_real.data: per time block yield
    (t_val, t_unit, kpoints[nk,3], pop4[nk,4]) where pop4 columns are the diabatic
    populations of VB-1, VB, CB1, CB2."""
    t_val, t_unit, kk, pp = None, '', [], []
    for line in open(filepath):
        s = line.strip()
        if s.startswith('# t ='):
            if t_val is not None and kk:
                yield t_val, t_unit, np.array(kk), np.array(pp)
            parts = s.split('=')[1].split()
            t_val = float(parts[0]); t_unit = parts[1] if len(parts) > 1 else 'a.u.'
            kk, pp = [], []
            continue
        if s.startswith('#') or not s:
            continue
        v = s.split()
        if len(v) >= 8:
            kk.append([float(v[1]), float(v[2]), float(v[3])])
            pp.append([float(v[4]), float(v[5]), float(v[6]), float(v[7])])
    if t_val is not None and kk:
        yield t_val, t_unit, np.array(kk), np.array(pp)


def plot_primitive_spectral(filepath, bpfile, output_dir, dpi=150, max_frames=150,
                            occupation_mode=True, autoscale=False):
    """A(k,E) spectral movie for a PRIMITIVE (unfolded) cell -- ONE FRAME PER STEP.

    Skeleton: the clean primitive bands from *_bandpath.data (thin grey lines).
    Decoration: if the FOUR-level file SYSNAME_sbe_nex_k_lev_real.data is present
    every gap-edge band (VB-1, VB, CB1, CB2) is coloured. Two colouring modes:

    * OCCUPATION (default, occupation_mode=True): each band coloured by its
      FRACTIONAL OCCUPATION f = pop/occ_full in [0,1]. The valence bands start
      FULL (f=1, bright) at t=0 and DEPLETE as holes form; the conduction bands
      start EMPTY (f=0, dark) and FILL. This shows the population *and its
      evolution* directly -- the valence band is visible from t=0 (fixing the
      earlier "valence shows nothing until carriers appear" behaviour, which
      coloured the valence by the hole density occ-pop = 0 at equilibrium).
    * EXCITATION (occupation_mode=False): holes (occ-pop) in the valence bands,
      electrons (pop) in the conduction bands -- both 0 at t=0, growing with
      excitation (the carrier/excitation view; pass --spectral-excitation).

    Otherwise only CB1 is coloured from the LCB file (legacy). Two views/frame go
    into `spectral_frames/` (path + kx projection); assemble with
    `ffmpeg -i nex_k_prim_spectral_path_f%04d*.png movie.mp4`."""
    from matplotlib.collections import LineCollection
    dist, eig_ha, nv, spinor, nodes, qred = _load_bandpath(bpfile)
    if not nv:
        print("  (skip) band path has no nv header — cannot identify levels")
        return
    levels = _bandpath_level_energies(eig_ha, nv, spinor)
    kx_path = qred[:, 0]
    occ_full = 1.0 if spinor else 2.0
    node_lbl = [n[0] for n in nodes]; node_dst = [n[1] for n in nodes]

    # Prefer the 4-level file (colour all bands); else fall back to the LCB file.
    lev_file = filepath.parent / filepath.name.replace('_nex_k_real.data',
                                                       '_nex_k_lev_real.data')
    use_lev = lev_file.exists()
    src = lev_file if use_lev else filepath
    print(f"Processing {src.name}  (primitive spectral A(k,E) per frame, "
          f"{'4-level' if use_lev else 'CB1-only'}, skeleton {bpfile.name}) ...")

    # which bandpath levels to colour, and the carrier sign per level
    if use_lev:
        col_levels = [('vbm1', 0, 'hole'), ('vb', 1, 'hole'),
                      ('cb1', 2, 'elec'), ('cb2', 3, 'elec')]
    else:
        col_levels = [('cb1', None, 'elec')]

    # colour value per band: fractional occupation (default, valence visible from
    # t=0) or fractional carrier/excitation (holes in VB, electrons in CB).
    def colour_value(p, sign):
        if occupation_mode:
            return p / occ_full                       # f in [0,1], VB starts at 1
        return ((occ_full - p) if sign == 'hole' else p) / occ_full   # excitation

    def frames():
        if use_lev:
            for t, tu, kpts, pop4 in _iter_nex_k_lev_blocks(lev_file):
                yield t, tu, kpts, pop4
        else:
            for t, tu, kpts, pop, _l, _s in _iter_nex_k_blocks(filepath, unfold=False):
                yield t, tu, kpts, pop[:, None]

    # Pass 1: global colour scale over the carrier populations + grid spacing.
    # peak_e / peak_h separately: on a fine grid + long runs the thermalized
    # ELECTRONS spread over the whole zone (per-k f ~ 1e-3..1e-4) while the
    # HOLES stay concentrated (per-k ~ 0.1..1) -- on one absolute scale the CB
    # then renders black ("electrons don't change"). --spectral-autoscale
    # normalizes each carrier branch by its own run-peak (annotated below).
    peak, peak_e, peak_h, n_frames, gridk = 0.0, 0.0, 0.0, 0, None
    for _t, _tu, kpts, pcols in frames():
        if gridk is None:
            gridk = kpts
        n_frames += 1
        for name, col, sign in col_levels:
            p = pcols[:, col] if col is not None else pcols[:, 0]
            carrier = colour_value(p, sign)
            pk = float(np.nanmax(carrier)) if carrier.size else 0.0
            peak = max(peak, pk)
            if sign == 'elec':
                # electron-branch peak is the CARRIER content even in
                # occupation mode (f itself, CB starts empty)
                peak_e = max(peak_e, float(np.nanmax(p / occ_full)) if p.size else 0.0)
            else:
                peak_h = max(peak_h, float(np.nanmax((occ_full - p) / occ_full)) if p.size else 0.0)
    if n_frames == 0:
        print("  (skip) no data blocks found"); return
    spacing = _grid_spacing(gridk)
    # occupation: fixed [0,1] scale (full valence = 1); excitation: data peak.
    norm = mcolors.Normalize(vmin=0.0, vmax=(1.0 if occupation_mode else max(peak, 1e-12)))
    if autoscale:
        norm = mcolors.Normalize(vmin=0.0, vmax=1.0)
    stride = max(1, n_frames // max_frames)
    frame_dir = output_dir / 'spectral_frames'
    frame_dir.mkdir(parents=True, exist_ok=True)
    if not use_lev:
        clabel = 'CB1 real population'
    elif occupation_mode:
        clabel = 'occupation  f = pop/occ  (VB full = 1, CB empty = 0)'
    else:
        clabel = 'carrier / excitation (e in CB, h in VB)'
    if autoscale:
        if occupation_mode:
            # VB stays absolute (f in [0,1]); only the CB carrier colour gains
            clabel = (f'occupation; CB autoscaled: colour = f/{max(peak_e,1e-30):.2e} '
                      f'(true CB peak f = {peak_e:.2e})')
        else:
            clabel = (f'per-branch autoscale: h/{max(peak_h,1e-30):.2e}, '
                      f'e/{max(peak_e,1e-30):.2e} (true peaks)')

    # frame-loop colour: per-branch autoscaled when requested (VB stays
    # absolute in occupation mode -- only the CB carrier colour gains).
    def colour_scaled(p, sign):
        v = colour_value(p, sign)
        if not autoscale:
            return v
        if occupation_mode:
            if sign == 'elec':
                return (p / occ_full) / max(peak_e, 1e-30)
            return v
        return v / max(peak_h if sign == 'hole' else peak_e, 1e-30)

    n_written = 0
    for iframe, (t_val, t_unit, kpts, pcols) in enumerate(frames()):
        if iframe % stride != 0:
            continue
        safe_t = f'{t_val:.4f}'.replace('-', 'm').replace('+', 'p')
        tag = f'f{iframe:04d}_t{safe_t}{t_unit}'
        # colour value mapped onto the path, per coloured level
        mapped = {}
        for name, col, sign in col_levels:
            p = pcols[:, col] if col is not None else pcols[:, 0]
            mapped[name] = _map_primitive_population(qred, kpts, colour_scaled(p, sign), spacing)

        # view 1: along the high-symmetry path (thin skeleton + coloured bands)
        fig, ax = plt.subplots(figsize=(8, 6))
        for nm in UNFOLD_LEVELS:
            if nm in levels:
                ax.plot(dist, levels[nm][0], color='0.8', lw=0.5, zorder=1)
        lc_last = None
        for name, col, sign in col_levels:
            if name not in levels:
                continue
            e_b, kin = levels[name]
            kemax = max(np.nanmax(kin), 1e-9)
            lw = 1.0 + 6.0 * np.nan_to_num(kin) / kemax
            P = mapped[name]
            pts = np.column_stack([dist, e_b])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            lc = LineCollection(segs, cmap=CMAP_POP, norm=norm, zorder=3)
            lc.set_array(0.5 * (np.nan_to_num(P[:-1]) + np.nan_to_num(P[1:])))
            lc.set_linewidths(0.5 * (lw[:-1] + lw[1:]))
            ax.add_collection(lc); lc_last = lc
        for d in node_dst:
            ax.axvline(d, color='#888888', linestyle='--', lw=0.7)
        ax.set_xticks(node_dst)
        ax.set_xticklabels([r'$\Gamma$' if l == 'Gamma' else f'${l}$' for l in node_lbl])
        ax.set_xlim(dist[0], dist[-1]); ax.set_ylim(-6, 8)
        ax.axhline(0.0, color='tab:red', lw=0.8, alpha=0.7)
        ax.set_ylabel('Energy [eV]  (VBM = 0)')
        _what = 'occupation' if (use_lev and occupation_mode) else 'carrier population'
        ax.set_title(f'Primitive A(k,E): {"4 gap-edge bands" if use_lev else "CB1"} '
                     f'coloured by {_what},  t = {t_val:.3f} {t_unit}\n'
                     f'colour = {_what}, width = carrier kinetic energy')
        if lc_last is not None:
            plt.colorbar(lc_last, ax=ax, label=clabel)
        fig.tight_layout()
        fig.savefig(frame_dir / f'nex_k_prim_spectral_path_{tag}.png',
                    dpi=dpi, bbox_inches='tight'); plt.close(fig)

        # view 2: projected onto kx (scatter, all coloured levels)
        fig, ax = plt.subplots(figsize=(8, 6))
        for nm in UNFOLD_LEVELS:
            if nm in levels:
                ax.plot(kx_path, levels[nm][0], color='0.85', lw=0.4,
                        zorder=1, marker='.', ms=1.5, ls='none')
        im = None
        for name, col, sign in col_levels:
            if name not in levels:
                continue
            e_b, kin = levels[name]
            kemax = max(np.nanmax(kin), 1e-9)
            lw = 1.0 + 6.0 * np.nan_to_num(kin) / kemax
            P = mapped[name]; good = np.isfinite(P)
            im = ax.scatter(kx_path[good], e_b[good], c=P[good],
                            s=8 + 60 * lw[good], cmap=CMAP_POP, norm=norm,
                            edgecolors='none', alpha=0.85, zorder=3)
        ax.axhline(0.0, color='tab:red', lw=0.8, alpha=0.7)
        ax.set_xlabel('kx [reduced, primitive BZ]')
        ax.set_ylabel('Energy [eV]  (VBM = 0)')
        ax.set_title(f'Primitive carrier population vs (kx, E),  t = {t_val:.3f} {t_unit}')
        if im is not None:
            plt.colorbar(im, ax=ax, label=clabel)
        fig.tight_layout()
        fig.savefig(frame_dir / f'nex_k_prim_spectral_kx_{tag}.png',
                    dpi=dpi, bbox_inches='tight'); plt.close(fig)
        n_written += 1

    print(f"  saved {n_written} frame(s) x 2 views into {frame_dir.name}/ "
          f"(of {n_frames} time steps, stride {stride})")


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
    Parse SYSNAME_eigen.data.
    EPM format: block headers are comment lines '# ik = N', energies in Ha.
    DFT format (theory='dft' + yn_out_tm): block headers are data-ish lines
    'k=     N,  spin=     M', energies in the unit named by the '# 1:io,
    2:esp[eV], 3:occ' header (eV for unit_system='A_eV_fs') — converted to Ha.
    Data lines (both): 'ib  energy  occup'.
    Returns (eigen[nb, nk] in Ha, occup[nb, nk], vbm_ha).
    """
    ik_re = re.compile(r'#\s*ik\s*=\s*(\d+)')
    dft_k_re = re.compile(r'k\s*=\s*(\d+)\s*,\s*spin\s*=\s*(\d+)')
    eigen_map = {}     # 1-based ik → list of (energy_Ha, occup)
    current_k = None
    vbm = -np.inf
    e_conv = 1.0       # → Ha

    with open(eigenfile, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                if 'esp[eV]' in s:
                    e_conv = 1.0 / HA_TO_EV
                m = ik_re.match(s)
                if m:
                    current_k = int(m.group(1))
                    eigen_map.setdefault(current_k, [])
                continue
            m = dft_k_re.match(s)
            if m:
                if int(m.group(2)) != 1:      # spin-2 block: keep spin-1 only
                    current_k = None
                    continue
                current_k = int(m.group(1))
                eigen_map.setdefault(current_k, [])
                continue
            if current_k is None:
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                e, occ = float(parts[1]) * e_conv, float(parts[2])
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

# ===========================================================================
# Auto-assemble the frame series into animations (mp4 via ffmpeg, else gif)
# ===========================================================================

def _anim_group_key(png_path):
    """(group_stem, sort_number) for a frame PNG, or None if it is not part of a
    time/frame series. Strips a trailing _t<time><unit> token and any _f<index>
    token; frames that share the remaining stem animate together. Single images
    (k-t maps, band structures, *_vs_Time) carry neither token -> None."""
    d, base = os.path.split(png_path)
    stem = base[:-4] if base.lower().endswith('.png') else base
    mt = re.search(r'_t([mp]?)(\d+(?:\.\d+)?)', stem)      # _t<time> (m=neg, p=pos)
    mf = re.search(r'_f(\d+)', stem)                       # _f<frame index>
    if mt is None and mf is None:
        return None
    if mt is not None:
        num = float(mt.group(2)) * (-1.0 if mt.group(1) == 'm' else 1.0)
    else:
        num = float(mf.group(1))
    key = re.sub(r'_t[mp]?\d+(?:\.\d+)?[a-zA-Z.]*$', '', stem)   # trailing time+unit
    key = re.sub(r'_f\d+', '', key)                             # frame index (any pos)
    return os.path.join(d, key.rstrip('_')), num


def _ffmpeg_mp4(frames, out, fps):
    """Encode an mp4 from the ordered PNG list via ffmpeg. Returns True on success."""
    import subprocess, tempfile, shutil as _sh
    with tempfile.TemporaryDirectory() as td:
        for i, f in enumerate(frames):                    # sequential names for -i %05d
            dst = os.path.join(td, f'{i:05d}.png')
            try:
                os.symlink(os.path.abspath(f), dst)
            except OSError:
                _sh.copyfile(f, dst)
        cmd = ['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps),
               '-i', os.path.join(td, '%05d.png'),
               # pad to even dimensions (yuv420p requirement) on a white canvas
               '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white',
               '-pix_fmt', 'yuv420p', out]
        try:
            subprocess.run(cmd, check=True)
            return os.path.exists(out)
        except Exception:
            return False


def _gif_pillow(frames, out, fps):
    """Encode an animated gif from the ordered PNG list via Pillow (a matplotlib
    dependency, so always available). Returns True on success."""
    try:
        from PIL import Image
        imgs = [Image.open(f).convert('RGB') for f in frames]
        imgs[0].save(out, save_all=True, append_images=imgs[1:],
                     duration=int(1000 / max(fps, 1)), loop=0, optimize=True)
        return True
    except Exception:
        return False


def _assemble_animations(output_dir, fps=6, fmt='auto', min_frames=2):
    """Group every frame-series PNG under output_dir and write one animation per
    series next to its frames (<stem>_anim.mp4 / .gif). mp4 when ffmpeg is on the
    PATH (fmt='auto'), else gif; 'mp4'/'gif'/'both' force the choice. Single
    images are skipped (they carry no _t/_f token)."""
    import glob as _glob, shutil as _sh
    groups = {}
    for p in _glob.glob(os.path.join(str(output_dir), '**', '*.png'), recursive=True):
        if p.endswith('_anim.png'):
            continue
        gk = _anim_group_key(p)
        if gk is not None:
            groups.setdefault(gk[0], []).append((gk[1], p))
    have_ffmpeg = _sh.which('ffmpeg') is not None
    made = []
    for stem, items in sorted(groups.items()):
        if len(items) < min_frames:
            continue
        frames = [p for _, p in sorted(items, key=lambda x: x[0])]
        want_mp4 = fmt in ('mp4', 'both') or (fmt == 'auto' and have_ffmpeg)
        want_gif = fmt in ('gif', 'both') or (fmt == 'auto' and not have_ffmpeg)
        ok = False
        if want_mp4:
            if _ffmpeg_mp4(frames, stem + '_anim.mp4', fps):
                made.append(stem + '_anim.mp4'); ok = True
            elif fmt == 'auto':
                want_gif = True                            # ffmpeg missing/failed -> gif
        if want_gif and _gif_pillow(frames, stem + '_anim.gif', fps):
            made.append(stem + '_anim.gif'); ok = True
        if not ok:
            print(f"  (warn) could not animate {os.path.basename(stem)} "
                  f"({len(frames)} frames) -- need ffmpeg or Pillow")
    if made:
        print(f"  animations ({fps} fps): "
              + ", ".join(os.path.basename(m) for m in made))
    return made


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
    parser.add_argument('--no-animate', action='store_true',
                        help='Do NOT auto-assemble the per-frame series (band '
                             'maps, k-maps, BZ snapshots, bz3d/voxel, spectral '
                             'frames) into an animation. By default every series '
                             'of >=2 frames is stitched into <stem>_anim.mp4 '
                             '(ffmpeg) or .gif (Pillow fallback).')
    parser.add_argument('--fps', type=int, default=6,
                        help='Frame rate for the auto-assembled animations '
                             '(default 6).')
    parser.add_argument('--anim-format', choices=['auto', 'mp4', 'gif', 'both'],
                        default='gif',
                        help="Animation container (default 'gif'): 'gif' (Pillow, "
                             "no external tool), 'mp4' (ffmpeg), 'both', or 'auto' "
                             "= mp4 if ffmpeg is on PATH else gif.")
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
    parser.add_argument('--lattice', choices=['auto', 'fcc', 'wurtzite'],
                        default='auto',
                        help='Lattice/symmetry-point set for the band-structure '
                             'plot. "auto" (default) reads "# material =" from the '
                             '_k.data header: CdS -> wurtzite (Gamma/X/Y/A/S '
                             'points, default path A-Gamma-X-S-Y-Gamma, so the '
                             'Gamma-A c-axis segment is included), everything else '
                             '-> fcc (cubic GaAs/Si: L-Gamma-X-W-K). Force with '
                             '"fcc"/"wurtzite".')

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
    parser.add_argument('--no-data-copy', action='store_true',
                        help='Do not snapshot the source *.data into output/data/ '
                             '(by default a copy is made so a later solver re-run, '
                             'which overwrites the run-dir *.data, cannot clobber '
                             'the data this plot set was made from)')
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
    parser.add_argument('--bz3d', action='store_true',
                        help='Paper-style 3D Brillouin-zone population plot per '
                             'time step: the BZ wireframe (Wigner-Seitz cell) + '
                             'the MP k-points coloured & sized by population '
                             '(weak points fade/hidden). Needs the reciprocal '
                             'vectors in the k.data header (primitive datasets).')
    parser.add_argument('--bz3d-voxel', action='store_true',
                        help='Variant (b) of --bz3d: a semi-transparent VOXEL '
                             'CLOUD rendered from the un-sheared Cartesian '
                             'population grid (the *_cart_snap_* array) instead '
                             'of the MP-point scatter -- shows the SHAPE of the '
                             'populated valleys. By default UNSMOOTHED (one cube '
                             'per populated k-bin, faithful to --bz3d); see '
                             '--voxel-smooth. Needs the reciprocal vectors in '
                             'the k.data header (primitive datasets).')
    parser.add_argument('--voxel-smooth', type=float, default=0.0, metavar='SIGMA',
                        help='Gaussian smoothing width (in grid bins) for '
                             '--bz3d-voxel. Default 0 = unsmoothed crisp cubes '
                             '(recommended for coarse grids; the colorbar always '
                             'reflects the TRUE peak population). Set e.g. 1.0 for '
                             'a soft cloud on dense grids.')
    parser.add_argument('--voxel-gamma', type=float, default=2.5, metavar='G',
                        help='Opacity non-linearity exponent for --bz3d-voxel: '
                             'alpha ∝ (population)^G. Default 2.5 keeps faint '
                             'overlapping cells transparent so the localized '
                             'maxima stand out. G=1 = linear opacity.')
    parser.add_argument('--voxel-gap', type=float, default=0.05, metavar='FRAC',
                        help='Fractional shrink of each --bz3d-voxel cube (0..0.9) '
                             'so neighbours separate and the interior stays '
                             'visible. Default 0.05 (small, good for dense grids); '
                             'raise for coarse grids.')
    parser.add_argument('--bz3d-cb-sum', action='store_true',
                        help='For --bz3d/--bz3d-voxel and the k-maps, use the SUM '
                             'of BOTH recorded conduction bands (CB1+CB2) as the '
                             'per-k population instead of the lowest CB only. '
                             'Reads the four-level file SYSNAME_sbe_nex_k_lev_real'
                             '.data (REAL carriers, primitive cell); outputs get a '
                             '_cbsum tag. Falls back to lowest-CB if absent.')
    parser.add_argument('--spectral-autoscale', action='store_true',
                        help='Normalize each carrier branch of the --spectral '
                             'colouring by its OWN run-peak (annotated in the '
                             'colorbar). Cures the "electrons look empty" view '
                             'on fine grids/long runs, where thermalized '
                             'electrons spread to per-k f ~ 1e-3..1e-4 while '
                             'holes stay concentrated (~0.1..1): on one '
                             'absolute 0..1 scale the CB renders black even '
                             'though the totals match. In occupation mode the '
                             'VB stays absolute; only the CB colour gains.')
    parser.add_argument('--spectral-excitation', action='store_true',
                        help='Colour the primitive spectral A(k,E) frames by '
                             'EXCITATION (holes in VB, electrons in CB; both 0 at '
                             't=0) instead of the default OCCUPATION (f=pop/occ; '
                             'valence full=1 from t=0, watch it deplete).')
    parser.add_argument('--valleys', action='store_true',
                        help='Overlay the FCC/diamond high-symmetry markers on the '
                             'Cartesian-BZ snapshot maps: Gamma (+), the six X '
                             'points (x), and the six Delta-valley minima at 0.85*X '
                             '(o) -- e.g. to check the Si hot spots sit in the '
                             'correct Delta valleys.')
    parser.add_argument('--instantaneous', action='store_true',
                        help='Also plot the instantaneous Houston-basis nex_k '
                             'maps (*_sbe_nex_k.data / *_unfold.data). These carry '
                             'the reversible A^2(t) virtual-polarization breathing. '
                             'By default only the REAL-carrier maps '
                             '(*_sbe_nex_k_real.data) are plotted when present.')
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

    # Tidy layout: the main line plots (fields A/E, currents, nex, conductivity)
    # stay in the output root; the k-space / Brillouin-zone evolution maps go into
    # a kspace/ subfolder so the root is not swamped by per-time-step frames.
    kdir = output_dir / 'kspace'
    kdir.mkdir(parents=True, exist_ok=True)

    # Snapshot the source *.data into output/data/ so a later re-run of the solver
    # (which overwrites the run-dir *.data in place, and truncates them on a
    # non-restart) does not clobber the results this plot set was made from.
    if not args.no_data_copy:
        ddir = output_dir / 'data'
        ddir.mkdir(parents=True, exist_ok=True)
        import shutil as _shutil
        for df in sorted(input_dir.glob('*.data')):
            try:
                dst = ddir / df.name
                if df.resolve() != dst.resolve():
                    _shutil.copy2(df, dst)
            except Exception as exc:
                print(f"  (skip data copy) {df.name}: {exc}")

    found_any = False

    # --- RT line files --------------------------------------------------
    if not args.no_rt:
        for f in sorted(input_dir.glob('*_sbe_channels.data')):
            found_any = True
            plot_channels(f, output_dir, dpi=args.dpi)
        for pattern in ('*_sbe_rt.data', '*_sbe_rt_energy.data', '*_sbe_nex.data',
                        '*_sbe_nex_nonad.data'):
            for f in sorted(input_dir.glob(pattern)):
                found_any = True
                plot_rt_file(f, output_dir, downsample=args.downsample, dpi=args.dpi)
        # full vs non-adiabatic overlay (twin Y-axes + shared log scale)
        for f in sorted(input_dir.glob('*_sbe_nex.data')):
            plot_nex_comparison(f, output_dir, dpi=args.dpi)

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

        # Intra-band (Houston-basis) current vs the total current
        for f in sorted(input_dir.glob('*_sbe_intra_current.data')):
            found_any = True
            rt = f.parent / (f.name[:-len('_sbe_intra_current.data')] + '_sbe_rt.data')
            plot_intra_current(f, rt if rt.exists() else None, output_dir, dpi=args.dpi)

        # REAL-carrier maps (default): fixed-basis diabatic population, the
        # k-resolved n_ex -- no reversible A^2(t) virtual breathing.
        real_k  = sorted(input_dir.glob('*_sbe_nex_k_real.data'))
        real_uk = sorted(input_dir.glob('*_sbe_nex_k_unfold_real.data'))
        for f in real_k:
            found_any = True
            plot_nex_k(f, kdir, dpi=args.dpi,
                       log_scale=args.log_cmap, snapshots=args.snapshots, real=True,
                       b_matrix=_bmatrix_for(f, '_sbe_nex_k_real.data'),
                       mark_valleys=args.valleys, bz3d=args.bz3d,
                       bz3d_voxel=args.bz3d_voxel, voxel_smooth=args.voxel_smooth,
                       cb_sum=args.bz3d_cb_sum, voxel_gap=args.voxel_gap,
                       voxel_gamma=args.voxel_gamma)
        for f in real_uk:
            found_any = True
            plot_nex_k(f, kdir, dpi=args.dpi,
                       log_scale=args.log_cmap, snapshots=args.snapshots,
                       unfold=True, real=True)

        # Instantaneous Houston-basis maps (carry the virtual breathing): only
        # when --instantaneous is set, or when no REAL file is present (old runs).
        if args.instantaneous or not real_k:
            for f in sorted(input_dir.glob('*_sbe_nex_k.data')):
                found_any = True
                plot_nex_k(f, kdir, dpi=args.dpi,
                           log_scale=args.log_cmap, snapshots=args.snapshots,
                           b_matrix=_bmatrix_for(f, '_sbe_nex_k.data'),
                           mark_valleys=args.valleys, bz3d=args.bz3d,
                       bz3d_voxel=args.bz3d_voxel, voxel_smooth=args.voxel_smooth)
                if args.subtract_baseline:
                    plot_nex_k(f, kdir, dpi=args.dpi,
                               log_scale=args.log_cmap, snapshots=args.snapshots,
                               subtract_baseline=True)

        # Physical (unfolded) CB1 populations on the primitive BZ
        unfold_for_spectral = real_uk or sorted(input_dir.glob('*_sbe_nex_k_unfold.data'))
        if args.instantaneous or not real_uk:
            for f in sorted(input_dir.glob('*_sbe_nex_k_unfold.data')):
                found_any = True
                plot_nex_k(f, kdir, dpi=args.dpi,
                           log_scale=args.log_cmap, snapshots=args.snapshots,
                           unfold=True)
                if args.subtract_baseline:
                    plot_nex_k(f, kdir, dpi=args.dpi,
                               log_scale=args.log_cmap, snapshots=args.snapshots,
                               unfold=True, subtract_baseline=True)

        # Optional spectral A(kx,E) map (one, from whichever unfold file exists)
        if args.spectral and unfold_for_spectral:
            for f in unfold_for_spectral[:1]:
                suffix = ('_sbe_nex_k_unfold_real.data' if f.name.endswith('_real.data')
                          else '_sbe_nex_k_unfold.data')
                stem = f.name[:-len(suffix)]
                bpfile = f.parent / f'{stem}_bandpath.data'
                if bpfile.exists():
                    plot_unfold_spectral(f, bpfile, kdir, dpi=args.dpi)
                else:
                    print(f"  (skip spectral) {bpfile.name} not found "
                          f"(generate it with: epm_gaas_reference.py bandpath)")
        elif args.spectral:
            # PRIMITIVE (unfolded) cell: no unfold file. Use the folded-format
            # LCB population (real preferred) + the primitive bandpath.
            prim = (real_k or sorted(input_dir.glob('*_sbe_nex_k.data')))
            for f in prim[:1]:
                suffix = ('_sbe_nex_k_real.data' if f.name.endswith('_real.data')
                          else '_sbe_nex_k.data')
                stem = f.name[:-len(suffix)]
                bpfile = f.parent / f'{stem}_bandpath.data'
                if bpfile.exists():
                    plot_primitive_spectral(f, bpfile, kdir, dpi=args.dpi,
                                            occupation_mode=not args.spectral_excitation,
                                            autoscale=args.spectral_autoscale)
                else:
                    print(f"  (skip spectral) {bpfile.name} not found")

    # --- Band structure -------------------------------------------------
    if not args.no_bands:
        for kf in sorted(input_dir.glob('*_k.data')):
            stem = kf.name[:-len('_k.data')]
            ef   = kf.parent / f'{stem}_eigen.data'
            if not ef.exists():
                print(f"  (skip bands) {ef.name} not found alongside {kf.name}")
                continue
            found_any = True
            lat = _detect_lattice(kf) if args.lattice == 'auto' else args.lattice
            if args.lattice == 'auto':
                print(f"  band path: lattice='{lat}' auto-detected from "
                      f"{kf.name} header"
                      + ("  (Gamma-A c-axis included)" if lat == 'wurtzite' else ""))
            try:
                plot_band_structure(
                    kf, ef, kdir,
                    path_labels=args.band_path,
                    energy_range_ev=tuple(args.energy_range),
                    dpi=args.dpi, spin_sum=args.spin_sum,
                    lattice=lat)
            except Exception as exc:
                print(f"  ERROR in band structure for {kf.name}: {exc}")

        # dft_band output (band.dat)
        for bf in sorted(input_dir.glob('band.dat')):
            found_any = True
            try:
                plot_band_dat(
                    bf, kdir,
                    energy_range_ev=tuple(args.energy_range),
                    dpi=args.dpi, vbm_index=args.band_vbm)
            except Exception as exc:
                print(f"  ERROR in dft_band plot for {bf.name}: {exc}")

        # Unfolded primitive-cell band path (epm_gaas_reference.py bandpath)
        for bf in sorted(input_dir.glob('*_bandpath.data')):
            found_any = True
            try:
                plot_bandpath(
                    bf, kdir,
                    energy_range_ev=tuple(args.energy_range),
                    dpi=args.dpi)
            except Exception as exc:
                print(f"  ERROR in bandpath plot for {bf.name}: {exc}")

    if not found_any:
        print(f"No data files found in {input_dir.resolve()}")
        print("Expected: *_sbe_rt.data, *_sbe_rt_energy.data, *_sbe_nex.data, "
              "*_sbe_nex_k.data, *_k.data + *_eigen.data, band.dat")
        return

    if not args.no_animate:
        _assemble_animations(output_dir, fps=args.fps, fmt=args.anim_format)

    print(f"\nDone.  Output: {output_dir.resolve()}")


if __name__ == '__main__':
    main()
