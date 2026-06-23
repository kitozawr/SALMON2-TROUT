#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dft_to_epm.py  --  extract EPM-compatible local form factors from a SALMON DFT
                   band structure.

WHY
---
The local Empirical Pseudopotential Method solver in this fork (``theory='epm'``,
``src/epm``) needs a small table of symmetric / antisymmetric local form factors
V^S(|G|^2), V^A(|G|^2) keyed by the reciprocal-shell index |G|^2 (in (2*pi/a)^2
units).  Cohen-Bergstresser / Kunikiyo tabulated those numbers only for a handful
of diamond / zincblende semiconductors.  For every *other* crystal there is no
table -- but SALMON can still produce a perfectly good Kohn-Sham band structure
with ``theory='dft'``.

This tool closes that gap: it reads a SALMON DFT band structure and *fits* the
EPM local form factors so that the EPM band structure reproduces the DFT one in
the least-squares sense.  The fitted numbers are written in three forms:

  1. ``<prefix>_epm_formfactors.data``  -- a table that the patched
     ``epm_cohen_bergstresser`` reader loads directly when you run
     ``theory='epm'`` with ``epm_material='file'`` and
     ``epm_formfactor_file='<prefix>_epm_formfactors.data'``.  No recompilation.
  2. ``<prefix>_fit_report.txt``        -- the fit quality (per-band / overall RMS)
     and the fitted V^S, V^A in both Rydberg and Hartree.
  3. a ready-to-paste Fortran ``case`` block for ``cb_get_form_factors`` (printed
     to stdout), in case you would rather hard-wire a new material.

The EPM forward model used in the fit is a faithful NumPy re-implementation of
the *same* plane-wave Hamiltonian SALMON builds:

  * ``--cell primitive``  mirrors ``src/epm/epm_solver.f90`` (2-atom FCC
    primitive cell, cutoff on |G|^2 in Hartree atomic units, tau = (a/8)(1,1,1)).
    Use this to fit data produced by ``theory='epm'`` itself, or by a DFT run on
    the 2-atom primitive diamond/zincblende cell.
  * ``--cell cubic``      mirrors ``epm_gaas_reference.py`` (conventional 8-atom
    simple-cubic cell with the FCC parity selection rule, cutoff on |G|^2 in
    (2*pi/a)^2 units).  Use this to fit a DFT run on the conventional cubic cell
    (e.g. ``samples/exercise_04_bulkSi_gs``).

Only the absolute energy zero differs between an LDA/GGA DFT run and the
zero-of-(local)-potential EPM convention, so a single rigid energy shift delta
is fitted alongside the form factors (residual = E_epm + delta - E_dft).

This is a *semi-empirical* extraction (band fitting), the standard route to EPM
form factors; it is in the same spirit as the machine-learned pseudopotentials
of the vendored DeePseudopot project (``external/DeePseudopot``).

USAGE (see --help)
------------------
  # fit the conventional cubic Si DFT band structure (band.dat from dft_band):
  python3 dft_to_epm.py --dft band.dat --format band_dat --cell cubic \
      --a-lattice-au 10.26 --material-name Si --shells-s 3,8,11 \
      --nval 16 --nbands-fit 18 --weight-valence 3.0 --out-prefix Si_dft

  # fit data written by theory='epm' (closed-loop self-test, primitive cell):
  python3 dft_to_epm.py --dft Si_epm_eigen.data --kfile Si_epm_k.data \
      --format salmon_eigen --cell primitive --a-lattice-au 10.26 \
      --material-name Si --shells-s 3,8,11 --out-prefix Si_recovered
"""

import argparse
import sys

import numpy as np
from scipy.linalg import eigvalsh
from scipy.optimize import least_squares

try:
    from deep_pp_adapter import zunger_form_factor, backend_name
except ImportError:  # allow `python3 tools/dft_to_epm/dft_to_epm.py` from repo root
    import os as _os
    sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    from deep_pp_adapter import zunger_form_factor, backend_name

RY_TO_HA = 0.5
HARTREE_EV = 27.211386245988
TWO_PI = 2.0 * np.pi


# =============================================================================
# EPM forward model (NumPy mirror of src/epm)
# =============================================================================
class EPMModel:
    """Local-EPM plane-wave Hamiltonian, parameterised by the form factors.

    The form factors enter linearly, so for each k-point we precompute, per
    active shell, the (Hermitian) cosine matrix that multiplies V^S and the
    (anti-Hermitian-times-i) sine matrix that multiplies V^A.  Re-evaluating the
    band structure for a trial parameter vector is then just a few axpy's plus a
    Hermitian eigensolve -- fast enough for least_squares with finite-difference
    Jacobians.
    """

    def __init__(self, cell, a_lattice, cutoff, shells_s, shells_a):
        self.cell = cell
        self.a = float(a_lattice)
        self.cutoff = float(cutoff)
        self.shells_s = list(shells_s)
        self.shells_a = list(shells_a)

        if cell == 'cubic':
            self._build_basis_cubic()
        elif cell == 'primitive':
            self._build_basis_primitive()
        else:
            raise ValueError("cell must be 'cubic' or 'primitive'")

        self._kcache = {}

    # ---- basis: conventional simple-cubic cell (FCC parity selection) -------
    def _build_basis_cubic(self):
        twopi_over_a = 2.0 * np.pi / self.a
        nmax = int(np.ceil(np.sqrt(self.cutoff))) + 1
        idx, gcart = [], []
        for h in range(-nmax, nmax + 1):
            for k in range(-nmax, nmax + 1):
                for l in range(-nmax, nmax + 1):
                    if h * h + k * k + l * l <= self.cutoff + 1e-8:
                        idx.append((h, k, l))
                        gcart.append(twopi_over_a * np.array([h, k, l]))
        self.Gidx = np.array(idx, dtype=int)
        self.Gcart = np.array(gcart, dtype=float)
        self.npw = len(idx)

    # ---- basis: 2-atom FCC primitive cell (mirror epm_solver.f90) -----------
    def _build_basis_primitive(self):
        h = 0.5 * self.a
        a1 = np.array([0.0, h, h]); a2 = np.array([h, 0.0, h]); a3 = np.array([h, h, 0.0])
        vol = np.dot(a1, np.cross(a2, a3))
        twopi = 2.0 * np.pi
        b1 = twopi * np.cross(a2, a3) / vol
        b2 = twopi * np.cross(a3, a1) / vol
        b3 = twopi * np.cross(a1, a2) / vol
        self.bmat = np.array([b1, b2, b3])
        bnorm = np.linalg.norm(self.bmat, axis=1)
        nmax = np.ceil(np.sqrt(self.cutoff) / bnorm).astype(int) + 1
        conv = (self.a / twopi) ** 2  # |G|^2 [a.u.] -> (2*pi/a)^2 units
        idx, gcart = [], []
        for m1 in range(-nmax[0], nmax[0] + 1):
            for m2 in range(-nmax[1], nmax[1] + 1):
                for m3 in range(-nmax[2], nmax[2] + 1):
                    G = m1 * b1 + m2 * b2 + m3 * b3
                    if np.dot(G, G) <= self.cutoff + 1e-8:
                        idx.append((m1, m2, m3))
                        gcart.append(G)
        self.Gidx = np.array(idx, dtype=int)
        self.Gcart = np.array(gcart, dtype=float)
        self.npw = len(idx)
        self.tau = self.a / 8.0 * np.ones(3)
        self._conv = conv

    # ---- per-k precompute ----------------------------------------------------
    def _shell_of(self, dG_cart, dG_idx):
        if self.cell == 'cubic':
            h, k, l = dG_idx
            # FCC parity selection rule: nonzero only if h,k,l all same parity.
            if not ((h % 2 == k % 2) and (k % 2 == l % 2)):
                return None, None
            g2 = int(h * h + k * k + l * l)
            phase = np.pi / 4.0 * (h + k + l)
            return g2, phase
        else:
            g2 = int(round(np.dot(dG_cart, dG_cart) * self._conv))
            phase = float(np.dot(dG_cart, self.tau))
            return g2, phase

    def _prep_k(self, kvec):
        key = (round(kvec[0], 12), round(kvec[1], 12), round(kvec[2], 12))
        if key in self._kcache:
            return self._kcache[key]
        npw = self.npw
        kpg = kvec[None, :] + self.Gcart
        diag = 0.5 * np.einsum('ij,ij->i', kpg, kpg)
        # per-shell cos/sin matrices (full npw x npw, Hermitian assembled later)
        cos_mats = {s: np.zeros((npw, npw)) for s in self.shells_s}
        sin_mats = {s: np.zeros((npw, npw)) for s in self.shells_a}
        for i in range(npw):
            for j in range(npw):
                if i == j:
                    continue
                dG_cart = self.Gcart[i] - self.Gcart[j]
                dG_idx = self.Gidx[i] - self.Gidx[j]
                g2, phase = self._shell_of(dG_cart, dG_idx)
                if g2 is None:
                    continue
                if g2 in cos_mats:
                    cos_mats[g2][i, j] = np.cos(phase)
                if g2 in sin_mats:
                    sin_mats[g2][i, j] = np.sin(phase)
        out = (diag, cos_mats, sin_mats)
        self._kcache[key] = out
        return out

    def bands(self, kpoints, vs_ry, va_ry, nb):
        """Lowest ``nb`` eigenvalues [Ha] of H(k) for each k in ``kpoints``.

        ``vs_ry``/``va_ry`` are dicts {shell: value_in_Rydberg}.
        """
        nk = len(kpoints)
        out = np.empty((nk, nb))
        for ik, kvec in enumerate(kpoints):
            diag, cos_mats, sin_mats = self._prep_k(np.asarray(kvec, dtype=float))
            H = np.diag(diag).astype(complex)
            for s, C in cos_mats.items():
                H += (vs_ry.get(s, 0.0) * RY_TO_HA) * C
            for s, S in sin_mats.items():
                H += 1j * (va_ry.get(s, 0.0) * RY_TO_HA) * S
            ev = eigvalsh(H)
            out[ik] = ev[:nb]
        return out


# =============================================================================
# Readers for SALMON DFT output
# =============================================================================
def read_band_dat(path):
    """Parse a SALMON ``band.dat`` (theory='dft_band').

    Returns (kpoints_cart [nk,3] in 1/Bohr, energies [nk, nb] in the run's energy
    unit).  The file layout (see src/gs/band_dft.f90) is:
        line 1: ' Number_of_Bands:           NB'
        line 2: ' Number_of_kpt_in_each_block:   NK'
        line 3: ' Number_of_blocks:    NBLK'
      then, per block: NK k-point lines  '<ik> kx_red ky_red kz_red kx ky kz'
                       followed by NK*NB energy lines '<ik> <iob> e1 [e2]'.
    """
    with open(path) as fh:
        lines = fh.readlines()

    def _int_after_colon(s):
        return int(s.split(':')[1].split()[0])

    NB = _int_after_colon(lines[0])
    NK = _int_after_colon(lines[1])
    NBLK = _int_after_colon(lines[2])

    kpoints, energies = [], []
    pos = 3
    for _ in range(NBLK):
        block_k = []
        for _ in range(NK):
            tok = lines[pos].split(); pos += 1
            # ik, kred(3), kcart(3)
            block_k.append([float(tok[4]), float(tok[5]), float(tok[6])])
        block_e = [[None] * NB for _ in range(NK)]
        for _ in range(NK * NB):
            tok = lines[pos].split(); pos += 1
            ik = int(tok[0]) - 1
            iob = int(tok[1]) - 1
            block_e[ik][iob] = float(tok[2])  # spin-1 (column 3)
        kpoints.extend(block_k)
        energies.extend(block_e)
    return np.array(kpoints), np.array(energies, dtype=float)


def read_salmon_eigen(eigen_path, k_path):
    """Parse SYSNAME_eigen.data + SYSNAME_k.data (the EPM/SBE GS format written
    by src/epm and read by gs_info_ssbe).  Energies are in Hartree, k in 1/Bohr.

    SYSNAME_k.data:     5 header lines, then 'ik kx ky kz weight'.
    SYSNAME_eigen.data: 3 header lines, then per ik a '# ik = ..' line followed
                        by nb lines 'ib energy occ'.
    """
    kpoints = []
    with open(k_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            tok = s.split()
            kpoints.append([float(tok[1]), float(tok[2]), float(tok[3])])
    kpoints = np.array(kpoints)

    blocks = []
    cur = None
    with open(eigen_path) as fh:
        header = 0
        for line in fh:
            s = line.strip()
            if s.startswith('#'):
                if 'ik' in s and '=' in s:
                    if cur is not None:
                        blocks.append(cur)
                    cur = []
                else:
                    header += 1
                continue
            if not s:
                continue
            if cur is None:
                cur = []
            tok = s.split()
            cur.append(float(tok[1]))  # energy [Ha]
        if cur is not None:
            blocks.append(cur)
    nb = min(len(b) for b in blocks)
    energies = np.array([b[:nb] for b in blocks], dtype=float)
    return kpoints, energies


# =============================================================================
# Fit
# =============================================================================
def _optimal_delta(e_epm, e_dft, weights2d):
    """Rigid energy shift delta* minimising the weighted squared band residual
    (E_epm + delta - E_dft); the closed form is the weighted mean of
    (E_dft - E_epm).  Eliminating delta analytically removes its (otherwise
    badly-conditioned) degeneracy from the nonlinear fit."""
    diff = e_dft - e_epm
    return np.sum(weights2d * diff) / np.sum(weights2d)


def run_fit(model, kpoints, e_dft_ha, shells_s, shells_a, nb_fit,
            band_weights=None):
    ns, na = len(shells_s), len(shells_a)
    e_dft = e_dft_ha[:, :nb_fit]
    if band_weights is None:
        band_weights = np.ones(nb_fit)
    w2d = np.broadcast_to(band_weights[None, :nb_fit], e_dft.shape)

    def unpack(p):
        vs = {s: p[i] for i, s in enumerate(shells_s)}
        va = {s: p[ns + i] for i, s in enumerate(shells_a)}
        return vs, va

    def residuals(p):
        vs, va = unpack(p)
        e_epm = model.bands(kpoints, vs, va, nb_fit)
        delta = _optimal_delta(e_epm, e_dft, w2d)
        return (np.sqrt(w2d) * ((e_epm + delta) - e_dft)).ravel()

    p0 = np.zeros(ns + na)
    # mild physical priors: V^S(3) attractive (negative), others small
    for i, s in enumerate(shells_s):
        p0[i] = -0.2 if s == 3 else 0.05

    sol = least_squares(residuals, p0, method='trf', xtol=1e-13, ftol=1e-13,
                        gtol=1e-13, max_nfev=4000)
    vs, va = unpack(sol.x)
    e_epm = model.bands(kpoints, vs, va, nb_fit)
    delta = _optimal_delta(e_epm, e_dft, w2d)
    res = (e_epm + delta) - e_dft
    return vs, va, delta, res, sol


def _zunger_factors(model, params_by_species, shells_s, shells_a):
    """Sample the Zunger local form V(|G|) at the EPM integer shells and turn it
    into symmetric / antisymmetric form factors.

    Monoatomic basis (diamond): one species -> V^S(s) = V(|G_s|), V^A = 0.
    Two-species basis (zincblende): V^S = (V_cat + V_an)/2, V^A = (V_cat - V_an)/2.
    ``params_by_species`` is a list of 4-vectors [a0,a1,a2,a3] (1 or 2 entries).
    """
    twopi_over_a = TWO_PI / model.a
    shells = sorted(set(shells_s) | set(shells_a))
    vs, va = {}, {}
    for s in shells:
        q = np.sqrt(s) * twopi_over_a  # |G_s| in Bohr^-1
        vals = [zunger_form_factor(q, p) for p in params_by_species]
        if len(vals) == 1:
            vs[s] = vals[0]
            va[s] = 0.0
        else:
            vs[s] = 0.5 * (vals[0] + vals[1])
            va[s] = 0.5 * (vals[0] - vals[1])
    return vs, va


def run_fit_zunger(model, kpoints, e_dft_ha, shells_s, shells_a, nb_fit,
                   band_weights=None, nspecies=1):
    """Fit the DeePseudopot analytic Zunger local form (a0..a3 per species) to
    the DFT bands, then sample it at the EPM shells. Smoother / more constrained
    than the per-shell least squares, so it extrapolates sensibly to shells that
    were not independently free."""
    e_dft = e_dft_ha[:, :nb_fit]
    if band_weights is None:
        band_weights = np.ones(nb_fit)
    w2d = np.broadcast_to(band_weights[None, :nb_fit], e_dft.shape)

    # p = [a0,a1,a2,a3] (* nspecies). a2 > 1 guarantees no real-q pole.
    p0_single = [0.2, 2.2, 2.0, 0.25]
    p0 = np.array(p0_single * nspecies, dtype=float)
    lo = np.array([-5.0, -10.0, 1.001, 1e-3] * nspecies)
    hi = np.array([5.0, 10.0, 50.0, 5.0] * nspecies)

    def split(p):
        return [p[4 * i:4 * i + 4] for i in range(nspecies)]

    def residuals(p):
        vs, va = _zunger_factors(model, split(p), shells_s, shells_a)
        e_epm = model.bands(kpoints, vs, va, nb_fit)
        delta = _optimal_delta(e_epm, e_dft, w2d)
        return (np.sqrt(w2d) * ((e_epm + delta) - e_dft)).ravel()

    sol = least_squares(residuals, p0, bounds=(lo, hi), method='trf',
                        xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=8000)
    vs, va = _zunger_factors(model, split(sol.x), shells_s, shells_a)
    e_epm = model.bands(kpoints, vs, va, nb_fit)
    delta = _optimal_delta(e_epm, e_dft, w2d)
    res = (e_epm + delta) - e_dft
    return vs, va, delta, res, sol, split(sol.x)


def residuals_unweighted(model, kpoints, e_dft_ha, vs, va, delta, nb_fit):
    e_epm = model.bands(kpoints, vs, va, nb_fit)
    return (e_epm + delta) - e_dft_ha[:, :nb_fit]


# =============================================================================
# Output
# =============================================================================
def write_formfactor_file(path, material, a_lat, cell, shells_s, shells_a, vs, va):
    all_shells = sorted(set(shells_s) | set(shells_a))
    with open(path, 'w') as fh:
        fh.write('# EPM local pseudopotential form factors '
                 '(compatible with src/epm)\n')
        fh.write('# generated by tools/dft_to_epm/dft_to_epm.py\n')
        fh.write('# material = %s   a_lattice_au = %.6f   cell = %s\n'
                 % (material, a_lat, cell))
        fh.write('# tau convention: zincblende (a/8)(1,1,1)\n')
        fh.write('# units: VS, VA in Rydberg '
                 '(cb_get_form_factors multiplies by 0.5 -> Hartree)\n')
        fh.write('# columns:  G2   VS_ry   VA_ry\n')
        for s in all_shells:
            fh.write('%-6d %18.10e %18.10e\n'
                     % (s, vs.get(s, 0.0), va.get(s, 0.0)))


def fortran_case_block(material, shells_s, shells_a, vs, va):
    all_shells = sorted(set(shells_s) | set(shells_a))
    lines = []
    lines.append("        case ('%s')" % material)
    lines.append("            select case (G2)")
    for s in all_shells:
        lines.append("            case (%d)" % s)
        lines.append("                VS_ry = %14.8fd0;  VA_ry = %14.8fd0"
                     % (vs.get(s, 0.0), va.get(s, 0.0)))
    lines.append("            end select")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fit EPM-compatible local form factors to a SALMON DFT band "
                    "structure.")
    ap.add_argument('--dft', required=True,
                    help="DFT band file: band.dat (theory='dft_band') or "
                         "SYSNAME_eigen.data (theory='epm'/SBE format).")
    ap.add_argument('--kfile', default=None,
                    help="SYSNAME_k.data (required for --format salmon_eigen).")
    ap.add_argument('--format', choices=['band_dat', 'salmon_eigen'],
                    default='band_dat')
    ap.add_argument('--cell', choices=['cubic', 'primitive'], default='cubic',
                    help="EPM forward-model cell. Must match the DFT cell: "
                         "cubic = conventional 8-atom cell, primitive = 2-atom "
                         "FCC cell / theory='epm' output.")
    ap.add_argument('--a-lattice-au', type=float, required=True,
                    help="(cubic) conventional lattice constant a [Bohr]; "
                         "(primitive) same a [Bohr] that defines the FCC vectors.")
    ap.add_argument('--cutoff-ry', type=float, default=11.1,
                    help="Plane-wave cutoff. cubic: |G|^2 in (2pi/a)^2 units; "
                         "primitive: |G|^2 in a.u. (matches epm_pw_cutoff_ry).")
    ap.add_argument('--material-name', default='DFT',
                    help="Label written into outputs / the 'case' block.")
    ap.add_argument('--shells-s', default='3,8,11',
                    help="Comma list of |G|^2 shells with a symmetric V^S.")
    ap.add_argument('--shells-a', default='',
                    help="Comma list of |G|^2 shells with an antisymmetric V^A "
                         "(empty for diamond / monoatomic basis).")
    ap.add_argument('--dft-energy-unit', choices=['ha', 'ev'], default='ha',
                    help="Energy unit of the --dft file. SALMON writes BOTH "
                         "band.dat and *_eigen.data in Hartree (a.u.) regardless "
                         "of &units, so 'ha' is correct for native SALMON output. "
                         "Use 'ev' only if you pre-converted the file.")
    ap.add_argument('--nval', type=int, default=None,
                    help="Number of valence bands (for VBM reporting / default "
                         "fit window).")
    ap.add_argument('--nbands-fit', type=int, default=None,
                    help="Number of lowest bands included in the fit "
                         "(default: nval+4 if --nval given, else all available).")
    ap.add_argument('--weight-valence', type=float, default=1.0,
                    help="Extra least-squares weight on the valence bands.")
    ap.add_argument('--method', choices=['lsq', 'zunger'], default='lsq',
                    help="Fitting model. 'lsq': free per-shell form factors "
                         "(fast, exact for an EPM-generated reference). 'zunger': "
                         "fit the vendored DeePseudopot analytic Zunger local "
                         "form a0..a3 (per species) and sample it at the shells -- "
                         "smoother / more physically constrained.")
    ap.add_argument('--out-prefix', default='dft_epm')
    args = ap.parse_args(argv)

    shells_s = [int(x) for x in args.shells_s.split(',') if x.strip() != '']
    shells_a = [int(x) for x in args.shells_a.split(',') if x.strip() != '']

    # --- read DFT bands ----------------------------------------------------
    if args.format == 'band_dat':
        kpoints, e_dft = read_band_dat(args.dft)
    else:
        if args.kfile is None:
            ap.error("--format salmon_eigen requires --kfile")
        kpoints, e_dft = read_salmon_eigen(args.dft, args.kfile)

    if args.dft_energy_unit == 'ev':
        e_dft_ha = e_dft / HARTREE_EV
    else:
        e_dft_ha = e_dft.copy()

    nb_avail = e_dft_ha.shape[1]
    if args.nbands_fit is not None:
        nb_fit = min(args.nbands_fit, nb_avail)
    elif args.nval is not None:
        nb_fit = min(args.nval + 4, nb_avail)
    else:
        nb_fit = nb_avail

    # --- build EPM model & check it has enough bands -----------------------
    model = EPMModel(args.cell, args.a_lattice_au, args.cutoff_ry,
                     shells_s, shells_a)
    if model.npw < nb_fit:
        ap.error("EPM basis (%d plane waves) smaller than nb_fit=%d; raise "
                 "--cutoff-ry." % (model.npw, nb_fit))

    band_weights = np.ones(nb_fit)
    if args.nval is not None and args.weight_valence != 1.0:
        band_weights[:min(args.nval, nb_fit)] = args.weight_valence

    print("# dft_to_epm: cell=%s  a=%.4f Bohr  cutoff=%.3f  npw=%d  "
          "nk=%d  nb_fit=%d" % (args.cell, args.a_lattice_au, args.cutoff_ry,
                                model.npw, len(kpoints), nb_fit))
    print("# method=%s   symmetric shells V^S: %s   antisymmetric shells V^A: %s"
          % (args.method, shells_s, shells_a or '(none)'))

    zunger_params = None
    if args.method == 'zunger':
        nspecies = 2 if shells_a else 1
        print("# zunger backend: %s   species=%d" % (backend_name(), nspecies))
        vs, va, delta, res, sol, zunger_params = run_fit_zunger(
            model, kpoints, e_dft_ha, shells_s, shells_a, nb_fit,
            band_weights, nspecies)
    else:
        vs, va, delta, res, sol = run_fit(model, kpoints, e_dft_ha, shells_s,
                                          shells_a, nb_fit, band_weights)

    rms = np.sqrt(np.mean(res ** 2)) * HARTREE_EV
    maxerr = np.max(np.abs(res)) * HARTREE_EV
    per_band_rms = np.sqrt(np.mean(res ** 2, axis=0)) * HARTREE_EV

    # --- report ------------------------------------------------------------
    report = []
    report.append("EPM form-factor fit from SALMON DFT")
    report.append("=" * 52)
    report.append("material            : %s" % args.material_name)
    report.append("cell / forward model : %s" % args.cell)
    report.append("a_lattice [Bohr]    : %.6f" % args.a_lattice_au)
    report.append("cutoff               : %.4f" % args.cutoff_ry)
    report.append("plane waves (npw)    : %d" % model.npw)
    report.append("k-points fitted      : %d" % len(kpoints))
    report.append("bands fitted         : %d" % nb_fit)
    report.append("method               : %s" % args.method)
    if zunger_params is not None:
        report.append("zunger backend       : %s" % backend_name())
        for i, p in enumerate(zunger_params):
            report.append("  species %d a0..a3   : %12.6f %12.6f %12.6f %12.6f"
                          % (i, p[0], p[1], p[2], p[3]))
    report.append("global shift delta   : %.6f Ha (%.4f eV)"
                  % (delta, delta * HARTREE_EV))
    report.append("converged            : %s (nfev=%d, cost=%.3e)"
                  % (sol.success, sol.nfev, sol.cost))
    report.append("")
    report.append("overall band RMS     : %.4f eV" % rms)
    report.append("overall band MAXERR  : %.4f eV" % maxerr)
    report.append("")
    report.append("fitted local form factors:")
    report.append("  shell |G|^2    V^S [Ry]      V^A [Ry]      V^S [Ha]      V^A [Ha]")
    for s in sorted(set(shells_s) | set(shells_a)):
        report.append("  %5d     %12.6f  %12.6f  %12.6f  %12.6f"
                      % (s, vs.get(s, 0.0), va.get(s, 0.0),
                         vs.get(s, 0.0) * RY_TO_HA, va.get(s, 0.0) * RY_TO_HA))
    report.append("")
    report.append("per-band RMS [eV] (band index from 1):")
    for b in range(nb_fit):
        report.append("  band %3d : %.4f" % (b + 1, per_band_rms[b]))
    report_txt = "\n".join(report)

    rep_path = args.out_prefix + "_fit_report.txt"
    with open(rep_path, 'w') as fh:
        fh.write(report_txt + "\n")

    ff_path = args.out_prefix + "_epm_formfactors.data"
    write_formfactor_file(ff_path, args.material_name, args.a_lattice_au,
                          args.cell, shells_s, shells_a, vs, va)

    print(report_txt)
    print("\n# wrote %s" % ff_path)
    print("# wrote %s" % rep_path)

    # Primary consumer: the Python EPM (epm_gaas_reference.py). It can load the
    # table directly via load_form_factor_file('%s').
    print("\n# --- Python EPM (epm_gaas_reference.py, primary) ---")
    print("#   set MATERIAL='file' and FORM_FACTOR_FILE='%s', or paste:" % ff_path)
    print("_FITTED_FORM_FACTORS_RY = {")
    for s in sorted(set(shells_s) | set(shells_a)):
        print("    %2d: (%.6f, %.6f)," % (s, vs.get(s, 0.0), va.get(s, 0.0)))
    print("}")

    print("\n# --- legacy Fortran theory='epm' (deprecated) ---")
    print("#   &epm: epm_material='file', epm_formfactor_file='%s'" % ff_path)
    print(fortran_case_block(args.material_name, shells_s, shells_a, vs, va))
    return 0


if __name__ == '__main__':
    sys.exit(main())
