#!/usr/bin/env python3
"""
hhg_spectrum.py -- high-harmonic-generation (HHG) spectrum from a SALMON SBE run.

Reads SYSNAME_sbe_rt.data (the gauge-invariant TOTAL matter current Jm and the
external field E_ext) and emits the HHG intensity on a LOG scale versus harmonic
order. The radiated intensity is the acceleration form

    S(omega) = | omega * FFT[ w(t) J(t) ] |^2         (~ |FFT[dJ/dt]|^2)

with a Hann window w(t) to suppress spectral leakage. The driven axis (the field
polarization) and the fundamental omega_0 are auto-detected from E_ext; the
current is projected on that axis (parallel HHG). For a centrosymmetric crystal
(Si) driven along a low-index axis only ODD harmonics appear -- a good sanity
check against your TDDFT reference.

Usage:
    python3 hhg_spectrum.py [SYSNAME_sbe_rt.data] [--nmax 25] [--out hhg.png]
    (defaults: the single *_sbe_rt.data in the cwd, nmax=25, hhg_spectrum.png)

Only the CLEAN (velocity-gauge basis NOT overflowing) part of a run is physical:
if the SBE printed 'VG basis edge reached' the current beyond that time is
contaminated -- weaken the field or raise nstate and re-run (see wiki/06).
"""
import sys
import glob
import argparse
import numpy as np


def load_rt(path):
    """Return t[fs], E_ext[N,3], Jm[N,3] from a SALMON _sbe_rt.data file."""
    t, E, J = [], [], []
    for ln in open(path):
        if ln.startswith('#') or not ln.strip():
            continue
        p = ln.split()
        # columns: 1:t 5:6:7 E_ext 14:15:16 Jm  (1-indexed -> 0-indexed 0, 4:7, 13:16)
        t.append(float(p[0]))
        E.append([float(p[4]), float(p[5]), float(p[6])])
        J.append([float(p[13]), float(p[14]), float(p[15])])
    return np.asarray(t), np.asarray(E), np.asarray(J)


def dominant_axis(E):
    """Index of the field polarization axis (largest peak |E| component)."""
    return int(np.argmax(np.max(np.abs(E), axis=0)))


def fundamental_omega(t_fs, e_ax):
    """Fundamental angular frequency [rad/fs] from the E-field FFT peak."""
    dt = np.mean(np.diff(t_fs))
    f = np.fft.rfftfreq(len(t_fs), d=dt)                 # cycles/fs
    A = np.abs(np.fft.rfft(e_ax * np.hanning(len(e_ax))))
    A[0] = 0.0
    f0 = f[np.argmax(A)]                                  # cycles/fs
    return 2.0 * np.pi * f0, f0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rt', nargs='?', default=None, help='SYSNAME_sbe_rt.data')
    ap.add_argument('--nmax', type=float, default=25.0, help='max harmonic order')
    ap.add_argument('--out', default='hhg_spectrum.png')
    ap.add_argument('--tmax', type=float, default=None,
                    help='truncate at t<=tmax fs (drop VG-overflow tail)')
    a = ap.parse_args()

    rt = a.rt or (sorted(glob.glob('*_sbe_rt.data')) + [None])[0]
    if rt is None:
        sys.exit('no *_sbe_rt.data found (pass the path explicitly)')
    t, E, J = load_rt(rt)
    if a.tmax is not None:
        m = t <= a.tmax
        t, E, J = t[m], E[m], J[m]
    if len(t) < 16:
        sys.exit(f'too few time samples in {rt}')

    ax = dominant_axis(E)
    axis_name = 'xyz'[ax]
    w0, f0 = fundamental_omega(t, E[:, ax])
    dt = np.mean(np.diff(t))

    # Hann-windowed acceleration-form HHG intensity S(omega)=|omega FFT[w J]|^2
    win = np.hanning(len(t))
    Jw = np.fft.rfft(J[:, ax] * win)
    f = np.fft.rfftfreq(len(t), d=dt)                    # cycles/fs
    w = 2.0 * np.pi * f
    S = (w ** 2) * np.abs(Jw) ** 2
    order = f / max(f0, 1e-30)                           # harmonic order n = omega/omega0

    # normalize so the fundamental peak = 1
    ref = S[(order > 0.5) & (order < 1.5)]
    S = S / (ref.max() if ref.size and ref.max() > 0 else S.max())

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axp = plt.subplots(figsize=(9, 5))
    keep = order <= a.nmax
    axp.semilogy(order[keep], np.maximum(S[keep], 1e-16), lw=1.1, color='navy')
    for n in range(1, int(a.nmax) + 1):
        axp.axvline(n, color='0.85', lw=0.6, zorder=0)
        if n % 2 == 1:                                   # odd harmonics (centrosymmetric)
            axp.axvline(n, color='tab:orange', lw=0.8, ls='--', alpha=0.6, zorder=0)
    axp.set_xlim(0, a.nmax)
    axp.set_ylim(1e-12, 5.0)
    axp.set_xlabel(f'harmonic order  n = $\\omega/\\omega_0$   '
                   f'($\\omega_0$ = {f0:.4f} cyc/fs = {f0*4.1357:.3f} eV, drive $\\parallel$ {axis_name})')
    axp.set_ylabel(r'HHG intensity $|\omega\,\tilde J(\omega)|^2$  (norm., log)')
    axp.set_title(f'HHG spectrum from {rt}\n(dashed orange = odd orders; '
                  f'even orders suppressed for a centrosymmetric crystal)')
    axp.grid(True, which='both', alpha=0.15)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f'# drive axis = {axis_name}, omega_0 = {f0:.5f} cyc/fs '
          f'({f0*4.1357:.4f} eV), samples = {len(t)}, dt = {dt:.5f} fs')
    print(f'# wrote {a.out}')


if __name__ == '__main__':
    main()
