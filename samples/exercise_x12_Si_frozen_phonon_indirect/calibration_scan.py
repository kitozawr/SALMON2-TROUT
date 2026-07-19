#!/usr/bin/env python3
"""Field-scan calibration of the eph ring's interband (phonon-assisted BTBT)
channel against the Keldysh brackets -- produces the recommended
`sbe_eph_interband_scale` for this grid.

Workflow (per intensity I, same pulse shape as the x12 example):
    run eph-ON  -> on_<I>_nex.data
    run eph-OFF -> off_<I>_nex.data          (coherent baseline)
then:
    python3 calibration_scan.py on_3.0d+10_nex.data off_3.0d+10_nex.data ...
    (pairs in on,off order; intensity parsed from the file name '..._<I>_nex')

For each field the phonon-assisted generation G_SBE = nex_ON(final) -
nex_OFF(final) is compared to the Keldysh rate integral for the indirect gap
treated as direct (the UPPER bound of a phonon-assisted channel: the true
second-order matrix element is below the direct one). If the ratio
G_SBE / n_Keldysh(1.7) is roughly field-independent, one constant rescale of
the gap-straddling eph rates reproduces the bracket:
    sbe_eph_interband_scale ~ 1 / median(ratio).
Run this on YOUR production grid (nk changes the dressing per k-point and the
sigma_E matching, so the factor is grid-specific) -- 4^3 for the quick look,
9^3 for the final tune.

Caveat: at A0 >~ 0.3 a.u. the parabolic Keldysh model is only a bracket
(Bloch-excursion regime); trust the calibration best at the LOWER fields of
the scan, and check the printed per-field ratios for flatness.
"""
import sys, re
import numpy as np
from scipy.special import ellipk, ellipe, dawsn

OMEGA_AU = 0.0046; TW_AU = 800.0            # the x12 pulse
EG_EV = 1.70; MSTAR = 0.30                   # indirect gap as-if-direct bracket
HBAR = 1.0545718e-34; ME = 9.109e-31; QE = 1.602e-19; EPS0 = 8.854e-12
AU_T = 2.4188843e-17

def keldysh_W(E_vm, Eg_eV, hw_eV, mstar):
    if E_vm < 1e5: return 0.0
    m = mstar * ME; Eg = Eg_eV * QE; w = hw_eV * QE / HBAR
    gam = w * np.sqrt(m * Eg) / (QE * E_vm)
    g1 = gam**2 / (1 + gam**2); g2 = 1 / (1 + gam**2)
    K1, E1 = ellipk(g1), ellipe(g1); K2, E2 = ellipk(g2), ellipe(g2)
    x = Eg * (2/np.pi) * np.sqrt(1 + gam**2) / gam * E2 / (HBAR * w)
    nu = np.floor(x + 1)
    s = sum(np.exp(-n * np.pi * (K1 - E1) / E2)
            * dawsn(np.sqrt(np.pi**2 * (2*nu - 2*x + 2*n) / (4 * K2 * E2)))
            for n in range(80))
    W = (2*w/(9*np.pi)) * ((np.sqrt(1+gam**2)/gam) * (m*w/HBAR))**1.5 \
        * np.sqrt(np.pi/(2*K2)) * s * np.exp(-np.pi * nu * (K1 - E1) / E2)
    return W / 1e6                                             # cm^-3 s^-1

def n_keldysh(I_wcm2):
    E0 = np.sqrt(2 * I_wcm2 * 1e4 / (2.998e8 * EPS0))
    ts = np.linspace(1e-20, TW_AU, 2001) * AU_T
    env = np.cos(np.pi * (ts/AU_T - TW_AU/2) / TW_AU)**2
    W = np.array([keldysh_W(E0 * e_, EG_EV, OMEGA_AU * 27.2114, MSTAR) for e_ in env])
    return np.trapezoid(W, ts)

def final_nex(fn):
    return np.loadtxt(fn, comments='#')[-1, 1]

pairs = list(zip(sys.argv[1::2], sys.argv[2::2]))
if not pairs:
    sys.exit(__doc__)
print(f'{"I [W/cm2]":>12} {"G_SBE":>10} {"n_Keldysh(1.7)":>14} {"ratio":>10}')
ratios = []
for on_fn, off_fn in pairs:
    m = re.search(r'([0-9.]+)[de]\+?([0-9]+)', on_fn)
    I = float(m.group(1)) * 10**int(m.group(2))
    G = final_nex(on_fn) - final_nex(off_fn)
    nk_ = n_keldysh(I)
    r = G / max(nk_, 1e-30)
    ratios.append(r)
    print(f'{I:12.1e} {G:10.2e} {nk_:14.2e} {r:10.1e}')
med = np.median(ratios)
print(f'\nmedian ratio = {med:.1e}  (flat across fields => one factor suffices)')
print(f'recommended  sbe_eph_interband_scale = {1.0/med:.2e}  (for THIS grid/pulse)')
