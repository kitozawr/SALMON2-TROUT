#!/usr/bin/env python3
"""Population saturation on the Dirac cone -- what density does the ring see,
and is it bounded by the Rana generation/recombination balance?

    python3 saturation_check.py runs/E100kVcm_mem/graphene_sit [more prefixes ...]
            [--lz-bohr 37.794523] [--area-bohr2 18.7154] [--temp-k 300] [--vf-ms 0.96e6] [--plot]

Reads, per prefix P:
    P_sbe_nex.data         Houston conduction population (carries the reversible A^2(t) dressing)
    P_sbe_nex_nonad.data   col 2: dressed-conduction projection; col 3: the Option-A dressed-
                           reference density = what the ring dissipators actually see
    P_sbe_channels.data    cumulative ring ledger; dN_rana > 0 = net carrier multiplication,
                           < 0 = net Auger recombination (per cell)
All 3D cm^-3 values are cell-volume normalized (vacuum included): sheet density
n_2d = n_3d * L_z. The balance density of the two-branch Dirac plasma at the
temperature the Rana rates are evaluated at (the e-ph bath T in the solver) is

    n_i(T) = (pi/6) (k_B T / hbar v_F)^2      (= 8.1e10 cm^-2 at 300 K)

(tests/test_rana_saturation.f90): above it the pair population net-recombines,
below it net-generates -> the population SATURATES at n_i(T). Printed per run:
peak/final Houston and ring-visible densities, the virtual fraction the ring
would have seen without the filters, the Rana ledger phases, and n_final/n_i.
"""
import argparse
import os
import re
import sys

import numpy as np

AU_T_FS = 0.02418884326505
BOHR_CM = 0.52917721067e-8
KB_EV = 8.617333262e-5
HBAR_EVS = 6.582119569e-16
LZ_BOHR = 37.794523
AREA_BOHR2 = 4.648726**2 * np.sqrt(3.0) / 2.0


def n_intrinsic_cm2(T, vf_ms):
    return (np.pi / 6.0) * (KB_EV * T / (HBAR_EVS * vf_ms * 1e2))**2     # (eV / (eV s cm/s))^2 = cm^-2


def _time_to_fs(path, t):
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') and re.search(r'1:time\[', line, re.I):
                unit = re.search(r'1:time\[([^\]]*)\]', line, re.I).group(1)
                return t if unit.startswith('fs') else t * AU_T_FS
    return t * AU_T_FS


def load_xy(path):
    d = np.loadtxt(path, comments='#')
    return _time_to_fs(path, d[:, 0]), d[:, 1:]


def analyze(prefix, lz, area, T, vf_ms, plot=False):
    lz_cm = lz * BOHR_CM
    area_cm2 = area * BOHR_CM**2
    out = {'prefix': prefix}
    t, y = load_xy(prefix + '_sbe_nex.data')
    nH = y[:, 0] * lz_cm                                  # cm^-2
    out.update(t=t, nH=nH, nH_peak=nH.max(), nH_final=nH[-1])
    p = prefix + '_sbe_nex_nonad.data'
    if os.path.exists(p):
        t2, y2 = load_xy(p)
        nproj, ndref = y2[:, 0] * lz_cm, y2[:, 1] * lz_cm
        out.update(t2=t2, nproj=nproj, ndref=ndref, ndref_peak=ndref.max(), ndref_final=ndref[-1],
                   virtual_frac_peak=1.0 - ndref.max() / max(nH.max(), 1e-300))
    p = prefix + '_sbe_channels.data'
    if os.path.exists(p):
        c = np.loadtxt(p, comments='#')
        tc, dn_rana = c[:, 0] * AU_T_FS, c[:, 7] / area_cm2      # cumulative pairs per cm^2
        out.update(tc=tc, dn_rana=dn_rana, rana_final=dn_rana[-1], rana_min=dn_rana.min(), rana_max=dn_rana.max())
    ni = n_intrinsic_cm2(T, vf_ms)
    out['n_i'] = ni
    print(f'--- {prefix}')
    print(f'  Houston n_2d:      peak {out["nH_peak"]:.3e}  final {out["nH_final"]:.3e} cm^-2')
    if 'ndref' in out:
        print(f'  ring-visible n_2d: peak {out["ndref_peak"]:.3e}  final {out["ndref_final"]:.3e} cm^-2 '
              f'(dressed reference)  -> virtual fraction at peak = {out["virtual_frac_peak"]:.3f}')
    if 'dn_rana' in out:
        phase = 'net carrier multiplication' if out['rana_final'] > 0 else 'net Auger recombination'
        print(f'  Rana ledger:       final {out["rana_final"]:+.3e} pairs/cm^2 ({phase}); '
              f'min {out["rana_min"]:+.2e} max {out["rana_max"]:+.2e}')
    print(f'  balance density n_i({T:.0f} K) = {ni:.3e} cm^-2  ->  n_final/n_i = {out["nH_final"] / ni:.3f}'
          f'   [n_i(1000 K) = {n_intrinsic_cm2(1000, vf_ms):.2e}, n_i(3000 K) = {n_intrinsic_cm2(3000, vf_ms):.2e}]')
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('prefixes', nargs='+', help='run prefixes, e.g. runs/E100kVcm_mem/graphene_sit')
    ap.add_argument('--lz-bohr', type=float, default=LZ_BOHR)
    ap.add_argument('--area-bohr2', type=float, default=AREA_BOHR2)
    ap.add_argument('--temp-k', type=float, default=300.0, help='temperature the Rana rates use (sbe_eph_temperature_k)')
    ap.add_argument('--vf-ms', type=float, default=1.0e6, help='v_F of the Rana constants [m/s]')
    ap.add_argument('--plot', action='store_true')
    args = ap.parse_args(argv)
    res = [analyze(p, args.lz_bohr, args.area_bohr2, args.temp_k, args.vf_ms) for p in args.prefixes]
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print('# matplotlib not available -- no plot')
            return 0
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
        for r in res:
            lab = os.path.basename(os.path.dirname(r['prefix'])) or r['prefix']
            ax[0].plot(r['t'], r['nH'], label=f'{lab}: Houston')
            if 'ndref' in r:
                ax[0].plot(r['t2'], r['ndref'], '--', label=f'{lab}: ring-visible (dressed ref.)')
            if 'dn_rana' in r:
                ax[1].plot(r['tc'], r['dn_rana'], label=lab)
        for T, ls in ((300, ':'), (1000, '-.')):
            ax[0].axhline(n_intrinsic_cm2(T, args.vf_ms), ls=ls, c='gray', lw=0.8, label=f'n_i({T} K)')
        ax[0].set_yscale('log'); ax[0].set_xlabel('t [fs]'); ax[0].set_ylabel('n_2d [cm$^{-2}$]'); ax[0].legend(fontsize=7)
        ax[1].axhline(0, c='k', lw=0.5); ax[1].set_xlabel('t [fs]'); ax[1].set_ylabel('cumulative Rana dN [pairs/cm$^2$]')
        ax[1].set_title('>0 carrier multiplication, <0 Auger recombination', fontsize=9); ax[1].legend(fontsize=7)
        fig.tight_layout(); fig.savefig('saturation_check.png', dpi=140)
        print('# wrote saturation_check.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
