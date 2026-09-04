#!/usr/bin/env python3
"""The driving transient: waveform, vector potential and spectrum.

    python3 plot_drive.py DAST_E100kVcm.txt [more files ...] [--out drive.png]
                          [--ef-ev 0.2] [--vf-au 0.439]

The input files are the ones `make_inputs.py` writes for `ae_shape1 = 'input'`:
columns t [fs], A_x, A_y, A_z [fs V/Angstrom]. Everything else is derived:
E(t) = -dA/dt, the spectrum |E(omega)|^2, and the two quantities the physics of a
doped sheet is measured against --

    A_0 = max |A(t)|                the displacement of the Fermi sea in reciprocal
                                    space (velocity gauge), against k_F = E_F/hbar v_F
    hbar*omega_0                    the intensity centroid of the spectrum,
                                    which sets the photon energy the sheet responds to

A single-cycle transient has no single frequency: the panel marks the centroid and the
FWHM band that `transmission.py` integrates over, and prints the DC content, which is
what makes A_0 (not E/omega of a nominal carrier) the right measure of the drive.
"""
import argparse
import os
import sys

import numpy as np

AU_T_FS = 0.02418884326505
AU_E_VM = 5.14220675e11
AU_EV = 27.211386245988
VF_DEF = 0.439


def load(path):
    d = np.loadtxt(path, comments='#')
    t = d[:, 0]                                   # fs
    A = d[:, 1:4]                                 # fs V/Angstrom
    ax = int(np.argmax(np.abs(A).max(axis=0)))
    a = A[:, ax] * 1e10 / AU_E_VM / AU_T_FS       # a.u.
    tau = t / AU_T_FS                             # a.u.
    e = -np.gradient(a, tau)                      # a.u.
    return t, a, e, 'xyz'[ax]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+')
    ap.add_argument('--ef-ev', type=float, default=0.2, help='doping to mark k_F against A_0')
    ap.add_argument('--vf-au', type=float, default=VF_DEF)
    ap.add_argument('--out', default='drive.png')
    args = ap.parse_args(argv)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('# matplotlib not available'); return 0
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.0))
    kF = (args.ef_ev / AU_EV) / args.vf_au
    for path in args.files:
        t, a, e, comp = load(path)
        lab = os.path.basename(path).replace('.txt', '')
        E0 = np.abs(e).max() * AU_E_VM / 1e5
        A0 = np.abs(a).max()
        ax[0].plot(t, e * AU_E_VM / 1e5, lw=1.3, label=f'{lab}: $E_0$ = {E0:.0f} kV/cm')
        ax[1].plot(t, a, lw=1.3, label=f'$A_0$ = {A0:.4f} a.u. = {A0 / kF:.2f} $k_F$')
        # zero-pad: the record is one cycle long, so the bare FFT bin (2 pi / T) is as
        # wide as the spectrum itself; padding interpolates the transform, it adds no
        # information but makes the centroid and the FWHM band readable.
        dt = (t[1] - t[0]) / AU_T_FS
        npad = 16 * len(e)
        w = np.fft.rfftfreq(npad, d=dt) * 2 * np.pi
        P = np.abs(np.fft.rfft(e, n=npad))**2
        m = w > 0
        w0 = np.sum(w[m] * P[m]) / np.sum(P[m])
        wpk = w[m][int(np.argmax(P[m]))]
        ax[2].plot(w[m] * AU_EV * 1e3, P[m] / P[m].max(), lw=1.3,
                   label=f'{lab.split("_")[-1]}: peak {wpk * AU_EV * 1e3:.1f} meV'
                         f' = {wpk * AU_EV / 4.1357e-3:.2f} THz')
        half = P[m] >= 0.5 * P[m].max()
        wl, wh = w[m][half][0] * AU_EV * 1e3, w[m][half][-1] * AU_EV * 1e3
        ax[2].axvspan(wl, wh, color='#f1c40f', alpha=0.18)
        ax[2].axvline(wpk * AU_EV * 1e3, ls=':', c='#7f8c8d', lw=1)
        print(f'# {lab}: E_0 = {E0:.2f} kV/cm, A_0 = {A0:.5f} a.u. = {A0 / kF:.2f} k_F (E_F = {args.ef_ev} eV),'
              f' spectral peak {wpk * AU_EV * 1e3:.2f} meV ({wpk * AU_EV / 4.1357e-3:.3f} THz),'
              f' centroid {w0 * AU_EV * 1e3:.2f} meV, FWHM {wl:.1f}-{wh:.1f} meV, polarisation {comp}')
    ax[0].set_xlabel('t [fs]'); ax[0].set_ylabel('$E(t)$ [kV/cm]'); ax[0].grid(alpha=0.25)
    ax[0].set_title('DAST single-cycle transient (the drive)', fontsize=9); ax[0].legend(fontsize=8)
    ax[1].axhline(kF, ls='--', c='#c0392b', lw=1)
    ax[1].text(t[len(t) // 6], kF * 1.05, f'$k_F$ at $E_F$ = {args.ef_ev:g} eV', color='#c0392b', fontsize=8)
    ax[1].set_xlabel('t [fs]'); ax[1].set_ylabel('$A(t)$ [a.u.] = displacement of the Fermi sea')
    ax[1].set_title('Vector potential: what displaces the sea in $k$-space', fontsize=9)
    ax[1].grid(alpha=0.25); ax[1].legend(fontsize=8)
    ax[2].set_xlim(0, 60); ax[2].set_xlabel(r'$\hbar\omega$ [meV]'); ax[2].set_ylabel(r'$|E(\omega)|^2$ (normalised)')
    ax[2].set_title('Spectrum; shaded = FWHM band used for Re $\\sigma$', fontsize=9)
    ax[2].grid(alpha=0.25); ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(args.out, dpi=150)
    print(f'# wrote {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
