#!/usr/bin/env python3
"""Velocity-gauge f-sum-rule completeness of a truncated SBE basis (wiki/12 sec. 6).

In the velocity gauge the current J = Tr[(p + A) rho]/V carries the DIAMAGNETIC
term A*N_e/V for every electron. For a FILLED band it must be cancelled by the
paramagnetic (interband) response, which happens exactly only with a COMPLETE
basis: per band n and Cartesian direction a,

    S_n^a(k) = sum_{m != n} 2 |p_nm^a(k)|^2 / (eps_m - eps_n) = 1 - d^2 eps_n / dk_a^2,

and the full-band average of the curvature vanishes, so <S_n^a>_k = 1. With nb
bands only, S^{nb} < 1 and the fraction

    eta_a = 1 - <S_n^{a,nb}>_{k, n occupied}

of the diamagnetic current is NOT cancelled: a spurious reactive sheet current
J_s = -eta (N_e/A_2D) A(t) proportional to A = E/omega -- negligible in the
near-IR, catastrophic at THz (the sheet turns into a plasma mirror). This
script computes eta_a from a GS dataset (SYS_tm.data, SYS_eigen.data,
SYS_k.data) and, if run directories are given, the DYNAMIC eta by fitting
J_s(t) against -(N_e/A_2D) A_tot(t) in *_sbe_rt.data, with T/R.

    python3 sumrule_check.py GSDIR/SYSNAME [--runs runs/*/SYSNAME_sbe_rt.data] [--occ 2] [--nval 1]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from transmission import read_rt, fluence_tra, LZ_BOHR, AREA_BOHR2  # noqa: E402


def read_gs(prefix):
    k = np.loadtxt(prefix + '_k.data', comments='#')
    nk = k.shape[0]
    w = k[:, 4]
    ev, cur = {}, None
    for line in open(prefix + '_eigen.data'):
        m = re.match(r'#\s*ik\s*=\s*(\d+)', line)
        if m:
            cur = int(m.group(1)); ev[cur] = []; continue
        if line.startswith('#') or not line.strip():
            continue
        ev[cur].append(float(line.split()[1]))
    nb = len(ev[1])
    eig = np.array([ev[i + 1] for i in range(nk)])                   # (nk, nb) Ha
    tm = np.loadtxt(prefix + '_tm.data', comments='#')[:nk * nb * nb]   # block 1 only (block 2 = rvnl)
    p = np.zeros((nk, nb, nb, 3), complex)
    ik = tm[:, 0].astype(int) - 1; ib = tm[:, 1].astype(int) - 1; jb = tm[:, 2].astype(int) - 1
    for a in range(3):
        p[ik, ib, jb, a] = tm[:, 3 + 2 * a] + 1j * tm[:, 4 + 2 * a]
    return k, w, eig, p


def static_eta(eig, p, w, nval, de_min=1e-3):
    """eta_a and the captured oscillator strength S_a averaged over the occupied bands."""
    nk, nb = eig.shape
    S = np.zeros(3); norm = 0.0; skipped = 0
    for n in range(nval):
        for m in range(nb):
            if m == n:
                continue
            de = eig[:, m] - eig[:, n]
            ok = np.abs(de) > de_min
            skipped += int((~ok).sum())
            for a in range(3):
                S[a] += np.sum(w[ok] * 2.0 * np.abs(p[ok, n, m, a])**2 / de[ok])
        norm += w.sum()
    S /= norm
    return 1.0 - S, S, skipped


def gs_current_adiabatic(eig, p, w, nval, axis, A_grid, occ=2.0):
    """Adiabatic ground-state current of the truncated H_k(A) = eps + A p_axis per cell area,
    J_gs(A) = sum_k w sum_{n<=nval} occ <phi_n(A)| p + A |phi_n(A)> / sum_k w  (the pure-gauge
    remainder the solver subtracts with yn_sbe_vg_sumrule = 'y'; = eta N A at linear order)."""
    nk, nb = eig.shape
    pa = p[:, :, :, axis]                                     # (nk, nb, nb)
    out = np.zeros(len(A_grid))
    eye = np.eye(nb)
    for i, A in enumerate(A_grid):
        H = np.einsum('kn,nm->knm', eig, eye) + A * pa
        ev, U = np.linalg.eigh(H)                             # ascending, columns
        V = U[:, :, :nval]                                    # (nk, nb, nval)
        cur = np.einsum('kin,kij,kjn->kn', V.conj(), pa + A * eye, V).real.sum(axis=1)
        out[i] = occ * np.sum(w * cur) / w.sum()
    return out


def run_variable(rt_path, key, default=None):
    vlog = os.path.join(os.path.dirname(rt_path) or '.', 'variables.log')
    if not os.path.exists(vlog):
        return default
    pat = re.compile(r'^#\s*' + re.escape(key) + r'\s*=\s*(\S+)')
    with open(vlog) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                return m.group(1).strip("'\"")
    return default


def dynamic_eta(rt_path, ne, area, lz, gs=None, nval=1):
    d = read_rt(rt_path); t = d['Time']
    ax = max('xyz', key=lambda a: np.max(np.abs(d[f'E_ext_{a}'])))
    iax = 'xyz'.index(ax)
    Js = -d[f'Jm_{ax}'] * lz
    Atot = d.get(f'Ac_tot_{ax}', d[f'Ac_ext_{ax}'])
    E_inc = d[f'E_ext_{ax}']; E_t = d.get(f'E_tot_{ax}', E_inc)
    Jdia = -(ne / area) * Atot
    m = np.abs(Atot) > 0.3 * np.abs(Atot).max()
    eta = np.sum(Js[m] * Jdia[m]) / max(np.sum(Jdia[m]**2), 1e-300)
    corrA = np.corrcoef(Js, Atot)[0, 1] if np.std(Js) > 0 else 0.0
    T, R, A = fluence_tra(t, E_inc, E_t, E_t - E_inc, 1.0)
    res = dict(eta=eta, corrA=corrA, T=T, R=R, A=A, Jpeak=np.abs(Js).max(),
               screen=np.abs(Atot).max() / max(np.abs(d[f'Ac_ext_{ax}']).max(), 1e-300),
               corrected=(run_variable(rt_path, 'yn_sbe_vg_sumrule', 'n') == 'y'), eta_res=np.nan)
    if gs is not None:
        # residual A-projection once the adiabatic ground-state current of the truncated
        # basis is removed (already inside Jm if the run had yn_sbe_vg_sumrule = 'y')
        k, w, eig, p = gs
        amax = max(np.abs(Atot).max(), 1e-12)
        Ag = np.linspace(-1.05 * amax, 1.05 * amax, 121)
        Jgs = -gs_current_adiabatic(eig, p, w, nval, iax, Ag) / area         # per-cell current -> sheet current (sign as Js)
        Jres = Js if res['corrected'] else Js - np.interp(Atot, Ag, Jgs)
        res['eta_res'] = np.sum(Jres[m] * Jdia[m]) / max(np.sum(Jdia[m]**2), 1e-300)
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gs_prefix', help='GSDIR/SYSNAME (reads _k/_eigen/_tm.data)')
    ap.add_argument('--runs', nargs='*', default=[], help='*_sbe_rt.data files for the dynamic fit')
    ap.add_argument('--occ', type=float, default=2.0, help='occupation per filled band')
    ap.add_argument('--nval', type=int, default=1, help='number of filled bands')
    ap.add_argument('--area-bohr2', type=float, default=AREA_BOHR2)
    ap.add_argument('--lz-bohr', type=float, default=LZ_BOHR)
    args = ap.parse_args(argv)
    k, w, eig, p = read_gs(args.gs_prefix)
    nk, nb = eig.shape
    eta, S, skipped = static_eta(eig, p, w, args.nval)
    print(f'# {args.gs_prefix}: nk = {nk}, nb = {nb}, filled bands = {args.nval}'
          + (f' ({skipped} (near-)degenerate pairs skipped, |de| < 1e-3 Ha: the exact Dirac point)' if skipped else ''))
    print(f'# captured oscillator strength <S>  x {S[0]:.4f}  y {S[1]:.4f}  z {S[2]:.4f}')
    print(f'# STATIC uncancelled diamagnetic fraction eta:  x {eta[0]:.4f}  y {eta[1]:.4f}  z {eta[2]:.4f}')
    ne = args.occ * args.nval
    files = sorted(sum((glob.glob(f) for f in args.runs), []))
    if files:
        print('# eta_dyn = A-projection of the sheet current in units of the diamagnetic N_e A/area (the reactive'
              ' artifact; ~0 after the pure-gauge restoration); eta_res = the same after subtracting the adiabatic'
              ' ground-state current of this truncated basis (run flag yn_sbe_vg_sumrule auto-detected: corr = y/n).')
        print('# NOTE for a DOPED run (sbe_ef_ev): a partially filled band has a REAL inductive current, so eta_res'
              ' is then the physical Drude weight in units of N_e/area, not an artifact -- eta_res * N_e/A_2D = D/pi'
              ' (use drude_check.py, which separates D and tau).')
        print(f'{"eta_dyn":>8} {"eta_res":>8} {"corr":>4} {"corr(J,A)":>9} {"|J|peak":>9} {"A_tot/A_ext":>11} {"T":>9} {"R":>9} {"A":>9}  file')
        for f in files:
            r = dynamic_eta(f, ne, args.area_bohr2, args.lz_bohr, gs=(k, w, eig, p), nval=args.nval)
            print(f'{r["eta"]:8.4f} {r["eta_res"]:8.4f} {"y" if r["corrected"] else "n":>4} {r["corrA"]:9.4f} {r["Jpeak"]:9.2e}'
                  f' {r["screen"]:11.3f} {r["T"]:9.5f} {r["R"]:9.2e} {r["A"]:9.5f}  {f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
