#!/usr/bin/env python3
"""
test_sheet_transmission.py - the x14 "field before / field after" estimator
(samples/exercise_x14_graphene_self_induced_transparency/transmission.py).

Sheet boundary condition at normal incidence (Hartree a.u., Z0 = 4 pi / c):
    E_t = (2 E_inc - Z0 J_s)/(1 + n_sub),  E_r = E_t - E_inc.
Checks:
  1) no current -> E_t = E_inc, T = 1 (free-standing); with a substrate the
     Fresnel value T = 4 n/(1+n)^2 (n = 1.5 -> 0.96, R = 0.04).
  2) LINEAR universal sheet sigma = e^2/(4 hbar) = 1/4 a.u., solved
     self-consistently (J = sigma E_t): the estimator returns the same E_t,
     T = 1/(1 + pi/2c)^2 = 0.97746, A = pi*alpha/(1+pi*alpha/2)^2 = 0.02241
     (pi*alpha = 0.022925), and T + R + A = 1 to 1e-12.
  3) POINTWISE energy identity for an ARBITRARY (non-linear, delayed) sheet
     current: (c/4pi)(E_inc^2 - E_t^2 - E_r^2) = E_t J_s.
  4) spectral (carrier) T/R/A agree with the fluence-integrated ones for the
     quasi-monochromatic linear sheet (1e-6).
  5) stack_from_one_layer: N sheets predicted from a one-layer run. For a LINEAR
     sheet the prediction is exact against the closed form |2/(2+N Z0 sigma)|^2 at
     any N, and it must be ORDERED (T_1 > T_2 > T_3) and equal to T_1 at N = 1.
  6) check_complete refuses a record shorter than the run's own nt.
Run:  python3 tests/test_sheet_transmission.py   (pure numpy)
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'samples', 'exercise_x14_graphene_self_induced_transparency'))
from transmission import (sheet_fields, fluence_tra, spectral_tra, linear_reference,  # noqa: E402
                          shell_resolution, stack_from_one_layer, check_complete,
                          TruncatedRun, C_AU, SIGMA_UNIV)

nfail = 0


def check(name, cond):
    global nfail
    if not cond:
        print(f"  FAIL: {name}")
        nfail += 1


# a quasi-monochromatic pulse: hw = 0.8 eV, 40-cycle Gaussian, dt = 2 a.u.
w0 = 0.8 / 27.211386245988
t = np.arange(0.0, 2.0e5, 2.0)
tc, tau = 1.0e5, 40.0 * 2 * np.pi / w0 / 2.355
E_inc = 1e-5 * np.exp(-0.5 * ((t - tc) / tau)**2) * np.cos(w0 * (t - tc))
Z0 = 4.0 * np.pi / C_AU

# --- (1) no current -----------------------------------------------------------
E_t, E_r = sheet_fields(E_inc, 0.0 * E_inc, 1.0)
check("J=0 free-standing: E_t == E_inc", np.allclose(E_t, E_inc, atol=0, rtol=1e-14))
check("J=0 free-standing: E_r == 0", np.max(np.abs(E_r)) == 0.0)
T, R, A = fluence_tra(t, E_inc, E_t, E_r, 1.0)
check("J=0 free-standing: T = 1", abs(T - 1) < 1e-14 and abs(R) < 1e-14 and abs(A) < 1e-14)
E_t, E_r = sheet_fields(E_inc, 0.0 * E_inc, 1.5)
T, R, A = fluence_tra(t, E_inc, E_t, E_r, 1.5)
check("J=0 substrate n=1.5: Fresnel T = 0.96", abs(T - 0.96) < 1e-12)
check("J=0 substrate n=1.5: Fresnel R = 0.04", abs(R - 0.04) < 1e-12)
check("J=0 substrate n=1.5: A = 0", abs(A) < 1e-12)

# --- (2) linear universal sheet, self-consistent ----------------------------------
tt = 2.0 / (2.0 + Z0 * SIGMA_UNIV)          # E_t/E_inc of the linear sheet, n_sub = 1
E_t_exact = tt * E_inc
J_s = SIGMA_UNIV * E_t_exact
E_t, E_r = sheet_fields(E_inc, J_s, 1.0)
check("linear sheet: estimator reproduces the self-consistent E_t",
      np.allclose(E_t, E_t_exact, rtol=1e-12, atol=0))
T, R, A = fluence_tra(t, E_inc, E_t, E_r, 1.0)
T_ref = 1.0 / (1.0 + np.pi / (2 * C_AU))**2
pa = np.pi / C_AU
A_ref = pa / (1.0 + pa / 2)**2
print(f"  linear universal sheet: T = {T:.6f} (ref {T_ref:.6f})  R = {R:.3e}  A = {A:.6f} (ref {A_ref:.6f}, pi*alpha = {pa:.6f})")
check("linear sheet: T = 1/(1+pi/2c)^2", abs(T - T_ref) < 1e-10)
check("linear sheet: A = pi*alpha/(1+pi*alpha/2)^2", abs(A - A_ref) < 1e-10)
check("linear sheet: A within 3% of pi*alpha (the graphene universal absorption)", abs(A - pa) / pa < 0.03)
check("linear sheet: T + R + A = 1", abs(T + R + A - 1) < 1e-12)
Tl, Rl, Al = linear_reference(SIGMA_UNIV, 1.0)
check("linear_reference() matches the fluence result", abs(Tl - T) < 1e-10 and abs(Al - A) < 1e-10)

# --- (3) pointwise energy identity for an arbitrary delayed nonlinear current -----
J_arb = 3e-7 * np.roll(E_inc, 150) + 2e-2 * E_inc**3 / max(np.max(np.abs(E_inc))**2, 1e-30)
E_t, E_r = sheet_fields(E_inc, J_arb, 1.0)
lhs = (C_AU / (4 * np.pi)) * (E_inc**2 - E_t**2 - E_r**2)
rhs = E_t * J_arb
check("pointwise identity (c/4pi)(E_i^2-E_t^2-E_r^2) = E_t J_s",
      np.max(np.abs(lhs - rhs)) < 1e-12 * max(np.max(np.abs(rhs)), 1e-30))

# --- (4) spectral vs fluence for the quasi-monochromatic linear sheet ------------
E_t, E_r = sheet_fields(E_inc, J_s, 1.0)
Ts, Rs, As, wc = spectral_tra(t, E_inc, E_t, E_r, 1.0)
check("spectral carrier found at hw = 0.8 eV", abs(wc - w0) / w0 < 2e-2)
check("spectral T equals fluence T for the linear sheet", abs(Ts - T) < 1e-6)
check("spectral A equals fluence A for the linear sheet", abs(As - A) < 1e-6)

# --- resolution advisory helper (documented numbers) ----------------------------
sp12 = shell_resolution(w0, 12)
sp150 = shell_resolution(w0, 150)
print(f"  resonance shell at 0.8 eV: {sp12:.2f} points/radius on 12x12, {sp150:.2f} on 150x150")
check("12x12 mesh cannot resolve the 0.8 eV shell (< 1 point/radius)", sp12 < 1.0)
check("150x150 mesh resolves the 0.8 eV shell (>= 3 points/radius)", sp150 >= 3.0)

# --- (5) N-layer stack predicted from the one-layer run -------------------------
# For a linear sheet the frequency-resolved prediction must reproduce the closed
# form at every N: sigma(omega) is constant, so the two are the same algebra.
E_t1, _ = sheet_fields(E_inc, J_s, 1.0)          # J_s = sigma E_t of section (2)
tprev = 1.1
for N in (1, 2, 3, 5):
    tn = stack_from_one_layer(t, E_inc, E_t1, J_s, N, 1.0)
    ref = abs(2.0 / (2.0 + N * Z0 * SIGMA_UNIV))**2
    check(f"{N}-layer stack from a one-layer linear run == |2/(2+{N}z)|^2", abs(tn - ref) < 1e-9)
    check(f"{N}-layer stack is darker than the {N - 1}-layer one", tn < tprev)
    tprev = tn
check("1-layer stack prediction reproduces the run itself",
      abs(stack_from_one_layer(t, E_inc, E_t1, J_s, 1, 1.0) - T) < 1e-9)

# --- (6) the truncated-record guard ---------------------------------------------
import tempfile  # noqa: E402

with tempfile.TemporaryDirectory() as td:
    rt = os.path.join(td, 'x_sbe_rt.data')
    open(rt, 'w').close()
    with open(os.path.join(td, 'variables.log'), 'w') as fh:
        fh.write('#    nt=  3844\n#    dt= 4.13414E+00\n')
    check("check_complete accepts a full record", check_complete(rt, 3844) == 3844)
    check("check_complete accepts a longer record", check_complete(rt, 3845) == 3844)
    raised = False
    try:
        check_complete(rt, 1764)
    except TruncatedRun:
        raised = True
    check("check_complete rejects a 46 % record", raised)
    check("allow_partial overrides it", check_complete(rt, 1764, allow_partial=True) == 3844)
    sub = os.path.join(td, 'elsewhere')
    os.makedirs(sub)
    check("no variables.log -> nothing to check",
          check_complete(os.path.join(sub, 'x_sbe_rt.data'), 1) is None)

if nfail == 0:
    print("PASS  (sheet BC: Fresnel limit, universal pi*alpha sheet, energy identity, "
          "spectral/fluence, N-layer stack from one run, truncated-record guard)")
    sys.exit(0)
print(f"FAIL ({nfail} checks)")
sys.exit(1)
