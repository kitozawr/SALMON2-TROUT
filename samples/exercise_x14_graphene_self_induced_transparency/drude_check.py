#!/usr/bin/env python3
"""Drude sector of a DOPED graphene sheet: conductivity, scattering time, mean free
path, mobility -- and what they do as the driving field grows.

    python3 drude_check.py runs/E100kVcm_diss/graphene_sit_sbe_rt.data [more ...]
            [--ef-ev 0.2] [--temp-k 300] [--n-sub 1.65] [--t-meas 0.60]

Why this script exists.  A doped sheet absorbs at THz through the intraband
(Drude) conductivity

    sigma(omega) = (D/pi) / (1/tau - i omega),     D = (2 e^2 kT/hbar^2) ln[2 cosh(mu/2kT)]

(the finite-temperature Drude weight of the Dirac cone; D -> e^2 mu/hbar^2 when
mu >> kT).  Field-induced transparency of such a sheet can only come from D or
from tau:

  * D at FIXED carrier density falls at most ~5-12 % on heating (mu drops as T_e
    rises), passes a minimum near T_e ~ 1000-1500 K and RISES again once thermal
    pairs are created -- printed below as a table, it cannot explain a large
    bleaching on its own;
  * tau falls hard as soon as carriers are accelerated past the optical-phonon
    threshold (E2g 196 meV, A1' 160 meV): each emission is a momentum-randomising
    event.  v_F A_0 = 0.74 eV at 100 kV/cm, so at THz the whole distribution is
    pushed over the threshold twice per cycle.

The script therefore extracts, per run, Re sigma over the incident band (from the
sheet current, exactly as transmission.py does), divides by the D of the initial
occupation and reports tau = pi Re(sigma)/D, the mean free path l = v_F tau and
the mobility mu_e = e v_F^2 tau / (E_F ...) -- i.e. the numbers a transport
measurement quotes.  With --t-meas it also converts a measured (substrate-included)
transmission into the same sheet conductivity, so calculation and experiment are
compared in one place.
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transmission import (read_rt, band_tra, sheet_fields, fluence_tra, run_variable,  # noqa: E402
                          LZ_BOHR, AREA_BOHR2, AU_E_VM, AU_EV, C_AU, SIGMA_UNIV, VF_EPM_AU)

KB_AU = 3.166811563e-6          # Hartree/K
BOHR_CM = 0.52917721067e-8
AU_T_FS = 0.02418884326505
Z0 = 4.0 * np.pi / C_AU         # vacuum impedance in a.u.


def drude_weight_au(mu_au, kT_au):
    """D = 2 kT ln[2 cosh(mu/2kT)] in atomic units (per the Dirac-cone Drude weight;
    equals |mu| in the degenerate limit, 2 kT ln2 at the Dirac point)."""
    x = abs(mu_au) / (2.0 * kT_au)
    return 2.0 * kT_au * (x + np.log1p(np.exp(-2.0 * x)))      # ln(2 cosh x) = x + ln(1+e^-2x)


_K59 = np.arange(1, 60)


def _fd1(eta):
    """Complete Fermi-Dirac integral F_1(eta) = int_0^inf t dt/(1+e^(t-eta)) = -Li_2(-e^eta),
    from the alternating series for |e^-|eta|| <= 1 plus the inversion
    F_1(eta) = eta^2/2 + pi^2/6 - F_1(-eta). Exact to 1e-13, no quadrature."""
    a = abs(float(eta))
    z = np.exp(-a)
    li = float(np.sum(((-1.0)**(_K59 + 1)) * z**_K59 / _K59**2))     # = F_1(-a)
    return li if eta < 0.0 else 0.5 * a * a + np.pi**2 / 6.0 - li


def density_cm2(mu_au, kT_au, vf=VF_EPM_AU):
    """n - p of the Dirac cone (g = 4) in cm^-2:
       n - p = (2/pi)(kT/hbar v_F)^2 [F_1(mu/kT) - F_1(-mu/kT)]  ->  (mu/hbar v_F)^2/pi at T -> 0."""
    x = mu_au / kT_au
    n2d = (2.0 / np.pi) * (kT_au / vf)**2 * (_fd1(x) - _fd1(-x))
    return n2d / BOHR_CM**2


def mu_of_density(n_cm2, T, vf=VF_EPM_AU):
    """Chemical potential holding n - p fixed at temperature T."""
    kT = KB_AU * T
    lo, hi = 0.0, 2.0
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if density_cm2(mid, kT, vf) < n_cm2:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sheet_from_transmission(T_meas, n_sub, faces=2):
    """Sheet conductance Z0*sigma from a measured transmission that still contains the
    substrate's own Fresnel losses (incoherent, no etalon): T_meas = T_bare * ratio,
    T_bare = [4 n/(1+n)^2] for `faces` air/substrate faces, ratio = [(1+n)/(1+n+Z0 sigma)]^2."""
    t_face = 4.0 * n_sub / (1.0 + n_sub)**2
    t_bare = t_face**faces
    ratio = T_meas / t_bare
    return (1.0 + n_sub) * (1.0 / np.sqrt(max(ratio, 1e-12)) - 1.0), t_bare


def analyze(path, ef_ev, temp_k, lz=LZ_BOHR, area=AREA_BOHR2):
    d = read_rt(path)
    t = d['Time']
    ax = max('xyz', key=lambda a: np.max(np.abs(d[f'E_ext_{a}'])))
    E_inc = d[f'E_ext_{ax}']
    n_layers = int(run_variable(path, 'sbe_sheet_nlayers', 1) or 1)
    J_s = -d[f'Jm_{ax}'] * lz * n_layers
    E_t = d.get(f'E_tot_{ax}', None)
    if E_t is None or np.max(np.abs(E_t - E_inc)) < 1e-12 * max(np.max(np.abs(E_inc)), 1e-300):
        E_t, E_r = sheet_fields(E_inc, J_s)
    else:
        E_r = E_t - E_inc
    T, R, A = fluence_tra(t, E_inc, E_t, E_r)
    Tb, Rb, Ab, resig, w0 = band_tra(t, E_inc, E_t, E_r, J_s)
    # Direct Drude extraction from the trajectory: the intraband sheet current of a
    # metal obeys   dJ_s/dt = (D/pi) E_local - J_s/tau,
    # so a least-squares fit of dJ/dt against (E_tot, J_s) returns the Drude weight D
    # AND the momentum-relaxation time tau of THIS run -- no assumption about which
    # of the two the field changes. Restricted to the driven window (|E| > 5 % of peak)
    # so the field-free tail does not dominate the normal equations.
    dJ = np.gradient(J_s, t)
    m = np.abs(E_t) > 0.05 * np.max(np.abs(E_t))
    D_fit, tau_fit = np.nan, np.nan
    if m.sum() > 20:
        M = np.column_stack([E_t[m], J_s[m]])
        coef, *_ = np.linalg.lstsq(M, dJ[m], rcond=None)
        D_fit = np.pi * coef[0]
        tau_fit = -1.0 / coef[1] if coef[1] < 0 else np.inf
    ef_run = float(run_variable(path, 'sbe_ef_ev', ef_ev) or ef_ev)
    ti_run = float(run_variable(path, 'sbe_temp_init_k', temp_k) or temp_k)
    kT = KB_AU * max(ti_run, 1.0)
    mu = abs(ef_run) / AU_EV
    D = drude_weight_au(mu, kT)                       # a.u.
    n0 = density_cm2(mu, kT)
    sig = resig * SIGMA_UNIV                          # a.u.
    tau = np.pi * sig / D if D > 0 else np.nan        # a.u. of time
    sig_fit = D_fit * tau_fit / np.pi if np.isfinite(D_fit * tau_fit) else np.nan
    return dict(file=path, E0_kvcm=float(np.max(np.abs(E_inc))) * AU_E_VM / 1e5, hw_ev=w0 * AU_EV,
                T=T, R=R, A=A, resig=resig, sigma_au=sig, ef_ev=ef_run, temp_init=ti_run,
                D_au=D, D_ev=D * AU_EV, n0_cm2=n0, tau_fs=tau * AU_T_FS,
                mfp_nm=VF_EPM_AU * tau * BOHR_CM * 1e7, Z0sig=Z0 * sig,
                D_fit_ev=D_fit * AU_EV, tau_fit_fs=tau_fit * AU_T_FS,
                mfp_fit_nm=VF_EPM_AU * tau_fit * BOHR_CM * 1e7,
                sig_fit_s0=sig_fit / SIGMA_UNIV)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='*', default=[], help='*_sbe_rt.data of DOPED runs (globs ok)')
    ap.add_argument('--ef-ev', type=float, default=0.2, help='fallback E_F if the run has no variables.log')
    ap.add_argument('--temp-k', type=float, default=300.0, help='fallback initial temperature')
    ap.add_argument('--n-sub', type=float, default=1.65, help='substrate index for the --t-meas conversion (PET ~ 1.65)')
    ap.add_argument('--t-meas', type=float, nargs='*', default=[],
                    help='measured transmission(s) INCLUDING the substrate, converted to the sheet conductance')
    args = ap.parse_args(argv)

    print('# Drude weight at fixed carrier density vs electron temperature (why heating alone cannot bleach much):')
    print(f'#  {"n [cm^-2]":>11} ' + ' '.join(f'{T:>7.0f}K' for T in (300, 1000, 2000, 3000, 5000)))
    for n_t in (1e12, 3e12, 1e13, 3e13):
        row = []
        d0 = None
        for T in (300, 1000, 2000, 3000, 5000):
            kT = KB_AU * T
            dd = drude_weight_au(mu_of_density(n_t, T), kT)
            d0 = dd if d0 is None else d0
            row.append(f'{dd / d0:8.3f}')
        print(f'#  {n_t:11.1e} ' + ' '.join(row) + '   (D/D(300 K))')

    if args.t_meas:
        print(f'# measured transmission -> sheet conductance (substrate n = {args.n_sub}, 2 faces, incoherent):')
        for tm in args.t_meas:
            zs, tb = sheet_from_transmission(tm, args.n_sub)
            print(f'#   T_meas = {tm:.3f}  ->  bare substrate {tb:.3f}, Z0*sigma = {zs:.4f}'
                  f' = {zs / Z0 / SIGMA_UNIV:6.1f} sigma_univ')

    files = sorted(sum((glob.glob(f) for f in args.files), []))
    if not files:
        return 0
    print('# D_fit, tau_fit: least-squares fit of dJ_s/dt = (D/pi) E_tot - J_s/tau over the driven window'
          ' (the run\'s OWN Drude weight and momentum-relaxation time); D_eq = the equilibrium Drude weight of E_F, T_init.')
    print(f'{"E0[kV/cm]":>9} {"hw[eV]":>7} {"E_F[eV]":>8} {"n0[cm^-2]":>10} {"D_eq[eV]":>9} {"D_fit[eV]":>10}'
          f' {"tau[fs]":>8} {"l[nm]":>7} {"s_fit/s0":>9} {"Re s/s0":>8} {"T":>8} {"R":>8} {"A":>8}  file')
    for f in files:
        r = analyze(f, args.ef_ev, args.temp_k)
        print(f'{r["E0_kvcm"]:9.2f} {r["hw_ev"]:7.4f} {r["ef_ev"]:8.3f} {r["n0_cm2"]:10.2e} {r["D_ev"]:9.4f}'
              f' {r["D_fit_ev"]:10.4f} {r["tau_fit_fs"]:8.2f} {r["mfp_fit_nm"]:7.1f} {r["sig_fit_s0"]:9.2f}'
              f' {r["resig"]:8.3f} {r["T"]:8.5f} {r["R"]:8.5f} {r["A"]:8.5f}  {f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
