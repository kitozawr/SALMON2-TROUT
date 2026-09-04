#!/usr/bin/env python3
"""Doped graphene sheet: the initial Fermi-Dirac occupation, the Dirac-cone Drude
weight, the mesh a uniform k-grid needs to carry a Fermi surface, and the two
field scales that control an intraband (Drude) response.

Pins, against closed-form results:

  1. FD occupation on a half-shifted MP mesh of the Dirac cone reproduces the
     analytic sheet density n = k_F^2/pi (g = 4) once the mesh resolves the Fermi
     circle, and the error shrinks as the mesh is refined. This is what
     `gs_info_ssbe` does when sbe_ef_ev / sbe_temp_init_k are set.
  2. The mesh criterion itself: the number of partially occupied k-points scales
     as the Fermi perimeter over the spacing; the solver warns below 20.
  3. Drude weight D(mu, T) = 2 kT ln[2 cosh(mu/2kT)]: degenerate limit |mu|,
     Dirac-point limit 2 kT ln 2, and -- the physically decisive statement --
     at FIXED carrier density heating changes D by at most ~10 % before turning
     round, so carrier heating alone cannot bleach a doped sheet by tens of per cent.
  4. Sheet conductance from a measured, substrate-included transmission (the
     inverse of the Fresnel + sheet formula used to compare with experiment).
  5. The two saturation scales of a doped Dirac sheet at frequency omega:
     the ballistic excursion A_0 = E/omega and the collisional one E*tau, each
     against k_F -- the onset of current (velocity) saturation.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'samples', 'exercise_x14_graphene_self_induced_transparency'))

AU_EV = 27.211386245988
KB_AU = 3.166811563e-6
BOHR_CM = 0.52917721067e-8
A_BOHR = 4.648726
VF = 0.439                       # 43-PW Ramanujam pi-EPM cone slope (0.96e6 m/s)
C_AU = 137.035999084
Z0 = 4.0 * np.pi / C_AU
SIGMA_UNIV = 0.25

fails = []


def check(name, got, want, tol, rel=True):
    err = abs(got - want) / abs(want) if rel and want != 0 else abs(got - want)
    ok = err <= tol
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: got {got:.6g}, want {want:.6g} ({'rel' if rel else 'abs'} err {err:.2e} <= {tol:.0e})")
    if not ok:
        fails.append(name)


def dirac_mesh(n, ef_au, kT_au):
    """Half-shifted MP mesh of the hexagonal cell; FD occupation of the two Dirac
    bands around K; returns (sheet density [cm^-2], partially occupied points)."""
    b = 4.0 * np.pi / (np.sqrt(3.0) * A_BOHR)
    b1 = np.array([b, 0.0]); b2 = np.array([b * 0.5, b * np.sqrt(3.0) / 2.0])
    i = (2 * np.arange(1, n + 1) - n - 1) / (2.0 * n)
    I, J = np.meshgrid(i, i, indexing='ij')
    kx = I * b1[0] + J * b2[0]; ky = I * b1[1] + J * b2[1]
    # the two inequivalent corners of the first BZ
    K1 = (2.0 * b1 + b2) / 3.0
    K2 = -K1
    dk = np.minimum(np.hypot(kx - K1[0], ky - K1[1]), np.hypot(kx - K2[0], ky - K2[1]))
    # umklapp: also measure through the neighbouring cells so points near the zone
    # edge find their nearest corner
    for g in (b1, b2, b1 + b2, -b1, -b2, -b1 - b2):
        for Kc in (K1, K2):
            dk = np.minimum(dk, np.hypot(kx - Kc[0] - g[0], ky - Kc[1] - g[1]))
    e_c = VF * dk          # conduction band from the Dirac point
    e_v = -VF * dk
    f_c = 2.0 / (1.0 + np.exp(np.clip((e_c - ef_au) / kT_au, -60, 60)))
    f_v = 2.0 / (1.0 + np.exp(np.clip((e_v - ef_au) / kT_au, -60, 60)))
    dne = (np.sum(f_c) + np.sum(f_v) - 2.0 * f_c.size) / f_c.size      # per cell
    area = A_BOHR**2 * np.sqrt(3.0) / 2.0
    nfs = int(np.sum((f_c > 0.1) & (f_c < 1.9)) + np.sum((f_v > 0.1) & (f_v < 1.9)))
    return dne / (area * BOHR_CM**2), nfs


def drude_weight(mu, kT):
    x = abs(mu) / (2.0 * kT)
    return 2.0 * kT * (x + np.log1p(np.exp(-2.0 * x)))


def n_of(mu, T):
    """n - p of the Dirac cone via the complete Fermi-Dirac integral F_1 = -Li_2(-e^eta)."""
    kk = np.arange(1, 60)
    def F1(eta):
        a = abs(eta); li = np.sum(((-1.0)**(kk + 1)) * np.exp(-a) ** kk / kk**2)
        return li if eta < 0 else 0.5 * a * a + np.pi**2 / 6.0 - li
    kT_ = KB_AU * T
    return (2.0 / np.pi) * (kT_ / VF)**2 * (F1(mu / kT_) - F1(-mu / kT_)) / BOHR_CM**2


def mu_of(n_t, T):
    lo, hi = 0.0, 2.0
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if n_of(mid, T) < n_t else (lo, mid)
    return 0.5 * (lo + hi)


print("== 1/2. FD occupation and Fermi-surface resolution on a uniform mesh ==")
ef = 0.2 / AU_EV
kT = KB_AU * 300.0
kF = ef / VF
prev_err = None
print(f"  E_F = 0.2 eV -> k_F = {kF:.5f} a.u.; n(T=0) = {kF**2 / np.pi / BOHR_CM**2:.4e},"
      f" n(300 K) = {n_of(ef, 300.0):.4e} cm^-2 (the mesh must reproduce the finite-T value)")
n_exact = n_of(ef, 300.0)
for N, tol in ((147, 0.05), (300, 0.01), (600, 0.005)):
    n_mesh, nfs = dirac_mesh(N, ef, kT)
    err = abs(n_mesh - n_exact) / n_exact
    print(f"  N = {N:4d}: n = {n_mesh:.4e} cm^-2 ({100 * err:5.2f} % off), partially occupied points = {nfs}")
    check(f"density at N={N}", n_mesh, n_exact, tol)
    if prev_err is not None and err > prev_err + 1e-12:
        fails.append(f"mesh convergence N={N}")
        print(f"  FAIL mesh convergence: error grew from {prev_err:.3f} to {err:.3f}")
    prev_err = err
n24, nfs24 = dirac_mesh(24, ef, kT)
print(f"  N = 24 (the smoke mesh): n = {n24:.3e} cm^-2, partially occupied points = {nfs24}"
      f" -> the solver's < 20 warning fires")
if nfs24 >= 20:
    fails.append("under-resolved mesh must give < 20 partially occupied points")
if n24 > 0.5 * n_exact:
    fails.append("24^2 mesh must badly under-count the doped density")

print("== 3. Dirac-cone Drude weight ==")
check("degenerate limit D -> |mu|", drude_weight(0.5 / AU_EV, KB_AU * 10.0), 0.5 / AU_EV, 1e-6)
check("Dirac point D -> 2 kT ln2", drude_weight(0.0, KB_AU * 300.0), 2.0 * KB_AU * 300.0 * np.log(2.0), 1e-12)


print("  D(T_e)/D(300 K) at fixed carrier density:")
worst_drop = 1.0
for n_t in (1e12, 3e12, 1e13):
    d0 = drude_weight(mu_of(n_t, 300.0), KB_AU * 300.0)
    row = []
    for T in (500, 1000, 1500, 2000, 3000):
        r = drude_weight(mu_of(n_t, T), KB_AU * T) / d0
        row.append(r); worst_drop = min(worst_drop, r)
    print(f"    n = {n_t:.0e} cm^-2: " + " ".join(f"{r:5.3f}" for r in row) + "   (500/1000/1500/2000/3000 K)")
print(f"  deepest relative Drude weight at any of those densities/temperatures: {worst_drop:.3f}")
if worst_drop < 0.80:
    fails.append("heating at fixed density must not drop D below ~0.8 (else the mechanism claim is wrong)")
    print("  FAIL: D drops more than 20 % -- revisit the 'heating alone cannot bleach' statement")
else:
    print("  ok   heating at fixed density changes D by <= 20 %: a large bleaching needs tau, not D")

print("== 4. measured transmission -> sheet conductance ==")


def sheet_conductance(T_meas, n_sub, faces=2):
    t_bare = (4.0 * n_sub / (1.0 + n_sub)**2)**faces
    return (1.0 + n_sub) * (np.sqrt(t_bare / T_meas) - 1.0)


zs = sheet_conductance(0.97746, 1.0)          # free-standing universal sheet
check("universal sheet round-trip", zs / Z0 / SIGMA_UNIV, 1.0, 2e-3)
for T_meas, want in ((0.60, 24.7), (0.70, 14.3)):
    got = sheet_conductance(T_meas, 1.65) / Z0 / SIGMA_UNIV
    check(f"PET sample T = {T_meas}", got, want, 0.02)

print("== 5. saturation scales of a doped Dirac sheet (3.36 THz) ==")
# A_0 is the MEASURED peak vector potential of the scaled DAST transient
# (make_inputs.py manifest): 6.213e-2 a.u. at 100 kV/cm -- a single-cycle pulse,
# so E/omega of the nominal 3.36 THz would understate it.
A0_PER_KVCM = 6.213e-4
for E_kv in (1.0, 10.0, 100.0):
    E = E_kv * 1e5 / 5.14220675e11
    A0 = A0_PER_KVCM * E_kv
    print(f"  E0 = {E_kv:6.1f} kV/cm: ballistic A_0 = {A0:.4f} a.u. = {A0 / kF:6.2f} k_F"
          f" ; collisional eE tau/hbar at tau = 63 fs = {E * 63 / 0.0241888 / kF:6.2f} k_F")
onset = kF / A0_PER_KVCM
print(f"  A_0 = k_F at E_0 = {onset:.1f} kV/cm -> above it the Fermi sea is displaced by more than its own")
print("  radius each half-cycle: the drift velocity saturates at v_F and the DIFFERENTIAL conductivity falls")
check("current-saturation onset (E_F = 0.2 eV, DAST)", onset, 26.9, 0.05)

print()
if fails:
    print(f"FAIL ({len(fails)}): " + "; ".join(fails))
    sys.exit(1)
print("PASS  (doped Dirac sheet: FD occupation and its mesh requirement, Drude weight limits and its weak"
      " temperature dependence at fixed density, transmission -> sheet conductance, current-saturation scales)")
