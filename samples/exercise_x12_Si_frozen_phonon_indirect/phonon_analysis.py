#!/usr/bin/env python3
"""Phonon-assisted indirect generation & cooling in Si (exercise x12).

Run the SBE step TWICE -- once with yn_sbe_eph='y' (ON) and once with 'n' (OFF)
-- copying the outputs aside between runs, e.g.:

    ../../build/salmon < Si_frozen_phonon_rt.inp > on.log          # eph = 'y'
    cp Si_prim_sbe_nex.data on_nex.data
    cp Si_prim_sbe_nex_k_real.data on_nkr.data
    sed "s/yn_sbe_eph               = 'y'/yn_sbe_eph               = 'n'/;\\
         s/yn_sbe_eph_acoustic      = 'y'/yn_sbe_eph_acoustic      = 'n'/" \\
        Si_frozen_phonon_rt.inp > off.inp
    ../../build/salmon < off.inp > off.log                          # eph = 'n'
    cp Si_prim_sbe_nex.data off_nex.data
    cp Si_prim_sbe_nex_k_real.data off_nkr.data
    python3 phonon_analysis.py

Produces phonon_assisted_demo.png (nex(t), energy distribution, cooling curve)
and prints the headline numbers. Electrons must read 8.000 in both logs.
"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
AU = 27.2114

def eigen_lcb(fn='Si_prim_eigen.data'):
    E = {}; ik = 0
    for line in open(fn):
        if line.startswith('# ik'): ik = int(line.split('=')[1]); E[ik] = []
        elif not line.startswith('#') and line.strip(): E[ik].append(float(line.split()[1]))
    ks = sorted(E)
    return np.array([E[k][4] for k in ks])            # band 5 = lowest conduction band

def nex_series(fn):
    d = np.loadtxt(fn, comments='#'); return d[:, 0] * 0.0242, d[:, 1]   # fs, cm^-3

def kblocks(fn):
    out = []
    for ch in open(fn).read().split('# t =')[1:]:
        L = ch.splitlines(); t = float(L[0].split()[0])
        r = [[float(x) for x in l.split()] for l in L[1:] if l.strip() and not l.startswith('#')]
        out.append((t, np.array(r)))
    return out

b5 = eigen_lcb(); cbm = b5.min()
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

for lbl, fn, c in [('eph OFF (field only)', 'off_nex.data', '#c44'),
                   ('eph ON (phonon)', 'on_nex.data', '#28a')]:
    t, n = nex_series(fn); ax[0].semilogy(t, n, c, lw=2.2, label=lbl)
ax[0].set(xlabel='time [fs]', ylabel=r'$n_{ex}$ (lowest CB) [cm$^{-3}$]')
ax[0].set_title('Phonon-assisted carrier generation\n(indirect Si: band-edge absorption IS phonon-assisted)', fontsize=9)
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

for lbl, fn, c in [('eph OFF', 'off_nkr.data', '#c44'), ('eph ON', 'on_nkr.data', '#28a')]:
    _, d = kblocks(fn)[-1]; nex = d[:, 4]; e = (b5 - cbm) * AU
    bins = np.linspace(0, 2.2, 12); h = np.zeros(len(bins) - 1)
    for ee, p in zip(e, nex):
        j = np.searchsorted(bins, ee) - 1
        if 0 <= j < len(h): h[j] += p
    h /= max(h.sum(), 1e-30)
    ax[1].step((bins[:-1] + bins[1:]) / 2, h, where='mid', color=c, lw=2.2, label=lbl)
ax[1].axvline(0, color='k', ls=':', lw=1)
ax[1].set(xlabel=r'$E - E_{CBM}$ [eV]', ylabel='fraction of lowest-CB carriers')
ax[1].set_title('Where the carriers end up\nphonon COOLS them to the X-valley', fontsize=9)
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)

for lbl, fn, c in [('eph OFF', 'off_nkr.data', '#c44'), ('eph ON', 'on_nkr.data', '#28a')]:
    ts, em = [], []
    for t, d in kblocks(fn):
        nex = d[:, 4]; tot = nex.sum()
        if tot > 1e-12:
            ts.append(t * 0.0242); em.append(((b5 - cbm) * AU * nex).sum() / tot)
    ax[2].plot(ts, em, c, lw=2.2, marker='o', ms=3, label=lbl)
ax[2].set(xlabel='time [fs]', ylabel=r'$\langle E - E_{CBM}\rangle$ [eV]')
ax[2].set_title('Carrier cooling by phonon emission', fontsize=9)
ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)

plt.suptitle('Phonon-assisted indirect (top-valence -> bottom-conduction) generation & cooling in Si '
             '| 2V+2C frozen Houston window, I=1e11 W/cm2', fontsize=9.5)
plt.tight_layout(rect=[0, 0, 1, 0.94]); plt.savefig('phonon_assisted_demo.png', dpi=115)

for lbl, fn in [('eph OFF', 'off_nkr.data'), ('eph ON ', 'on_nkr.data')]:
    _, d = kblocks(fn)[-1]; nex = d[:, 4]; tot = nex.sum()
    em = (b5 * nex).sum() / max(tot, 1e-30)
    print(f'{lbl}: sum n_ex(k)={tot:.3e}  <E-CBM>={(em - cbm) * AU:6.3f} eV')
print('saved phonon_assisted_demo.png')
