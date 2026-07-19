#!/usr/bin/env python3
"""Chefonov-style rate-equation benchmark for exercise x12 (Si, I=1e11 W/cm2).

Integrates the model rate equation
    dn/dt = W_Keldysh(E(t))                      [photo/tunnel generation]
on the SAME Acos2 pulse the example drives, and compares against the SBE
nex(t) (on_nex.data / off_nex.data from the README workflow).

Two Keldysh curves bracket the physics:
  * Eg = 3.39 eV  -- the EPM DIRECT Gamma gap: what a parabolic two-band
    direct-transition model predicts with NO phonon assistance;
  * Eg = 1.70 eV  -- the INDIRECT gap treated as if it were direct: an
    UPPER BOUND for any phonon-assisted channel (the true second-order
    matrix element is smaller than the direct one it is compared to).

Interpretation guide (see README "Rate-equation benchmark"):
  - SBE OFF above Keldysh(3.39): expected here -- A0 = 0.37 a.u. sweeps ~60%
    of the half-BZ (eEa ~ 0.5 eV >> hbar*w), far outside the parabolic
    two-band regime Keldysh assumes; the coherent SBE part cannot be judged
    by this comparison at this field.
  - SBE ON above Keldysh(1.7): the collision-assisted conversion exceeds the
    as-if-direct upper bound => the ABSOLUTE generation rate is an upper
    estimate (transient Houston population real-ified at the saturated
    intervalley nu_sat, energy matching width sigma_E). Distributions
    (which k / which band) are unaffected by this scale factor.
  - With impact ionization on, add   dn/dt += nu_II(eps(t)) * n   using the
    SAME registry constants the run prints (Keldysh II fit P, Eth) to
    benchmark the avalanche stage separately.

Usage:  python3 rate_benchmark.py     (after the README on/off workflow)
"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe, dawsn

# ---- pulse of the example (a.u.) -------------------------------------------
I_WCM2 = 1.0e11; OMEGA_AU = 0.0046; TW_AU = 800.0
HBAR = 1.0545718e-34; ME = 9.109e-31; QE = 1.602e-19; EPS0 = 8.854e-12
AU_T = 2.4188843e-17; AU_E_VM = 5.14220675e11
E0_VM = np.sqrt(2.0 * I_WCM2 * 1e4 / (2.998e8 * EPS0))

def keldysh_W(E_vm, Eg_eV, hw_eV, mstar):
    """Keldysh (1965) interband photoionization rate [cm^-3 s^-1]."""
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
    Q = np.sqrt(np.pi / (2 * K2)) * s
    W = (2*w/(9*np.pi)) * ((np.sqrt(1+gam**2)/gam) * (m*w/HBAR))**1.5 \
        * Q * np.exp(-np.pi * nu * (K1 - E1) / E2)
    return W / 1e6                                    # m^-3 -> cm^-3

# envelope of |E(t)| for the Acos2 vector potential (peak per half-cycle)
t_au = np.linspace(0, TW_AU, 4001); dt_s = (t_au[1] - t_au[0]) * AU_T
env = np.cos(np.pi * (t_au - TW_AU/2) / TW_AU)**2
curves = {}
for Eg, mst, lbl in [(3.39, 0.18, 'Keldysh, direct 3.39 eV'),
                     (1.70, 0.30, 'Keldysh, 1.7 eV as-if-direct (ph-assisted UPPER bound)')]:
    W = np.array([keldysh_W(E0_VM * e_, Eg, OMEGA_AU * 27.2114, mst) for e_ in env])
    curves[lbl] = (t_au * AU_T * 1e15, np.cumsum(W) * dt_s)   # fs, cm^-3

fig, ax = plt.subplots(figsize=(7.2, 4.8))
for lbl, (tf, n) in curves.items():
    ax.semilogy(tf, np.maximum(n, 1e6), '--', lw=1.8, label=lbl)
for lbl, fn, c in [('SBE eph OFF (coherent)', 'off_nex.data', '#c44'),
                   ('SBE eph ON (phonon ring)', 'on_nex.data', '#28a')]:
    try:
        d = np.loadtxt(fn, comments='#')
        ax.semilogy(d[:, 0] * 0.0242, np.maximum(d[:, 1], 1e6), c, lw=2.2, label=lbl)
    except OSError:
        print(f'({fn} not found -- run the README on/off workflow first)')
ax.set(xlabel='time [fs]', ylabel=r'$n$ [cm$^{-3}$]', ylim=(1e8, 1e23))
ax.set_title('x12 vs Keldysh rate equation -- ABSOLUTE-rate benchmark\n'
             r'(A$_0$=0.37 a.u. sweeps ~60% of the half-BZ: parabolic Keldysh is a bracket, not truth)',
             fontsize=9)
ax.legend(fontsize=8); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('rate_benchmark.png', dpi=115)

for lbl, (tf, n) in curves.items():
    print(f'{lbl}: n_final = {n[-1]:.2e} cm^-3')
for lbl, fn in [('SBE OFF', 'off_nex.data'), ('SBE ON ', 'on_nex.data')]:
    try:
        d = np.loadtxt(fn, comments='#')
        print(f'{lbl}: n_final = {d[-1, 1]:.2e} cm^-3')
    except OSError:
        pass
print('saved rate_benchmark.png')
