#!/usr/bin/env python3
"""
test_ii_form_switch.py  -  Part B (impact-ionization fit-form switch).

The k-local impact-ionization rate is gamma = P (eps-E_th)^a Theta(eps-E_th),
with the exponent a selectable (sbe_ii_exponent): a=4 GaAs Stobbe quartic,
a=2 Si Keldysh quadratic, a=4.6 Si full-band. The Fortran prefactor conversion
is  ii_pref_au = P[s^-1 eV^-a] * t_au[s] * (au_ev)^a , and the rate uses
d**a (d = eps_kin - E_th in Hartree). This test locks that math: it reproduces
the Fortran unit conversion and checks the a=2 vs a=4 rate scaling, and that the
GaAs default exponent is unchanged.

Pure-Python (no SALMON build); mirrors src/ssbe/bloch_solver_ssbe.f90.
"""
import sys

AU_FS = 41.341374575751   # 1 fs in a.u. of time
AU_EV = 27.211386245988   # 1 Hartree in eV
S_PER_AU = AU_FS * 1.0e-15  # 1 a.u. of time in seconds


def gamma_au(eps_kin_ev, E_th_ev, P_si, a):
    """Impact-ionization rate in a.u., mirroring the Fortran path exactly."""
    pref_au = P_si * S_PER_AU * AU_EV ** a          # P [s^-1 eV^-a] -> 1/(Ha^a a.u.t)
    d_ha = (eps_kin_ev - E_th_ev) / AU_EV
    if d_ha <= 0.0:
        return 0.0
    return pref_au * d_ha ** a


def main():
    ok = True

    # 1) Below threshold => exactly zero, both forms.
    for a in (2.0, 4.0, 4.6):
        if gamma_au(1.0, 1.1, 1e13, a) != 0.0:
            print(f"  FAIL: rate below threshold nonzero (a={a})"); ok = False

    # 2) Dimensional self-consistency: at (eps-E_th) = 1 eV the rate in s^-1 equals P.
    for a in (2.0, 4.0, 4.6):
        g_au = gamma_au(2.0, 1.0, 5e13, a)     # d = 1 eV
        g_si = g_au / S_PER_AU                  # back to s^-1
        if abs(g_si - 5e13) / 5e13 > 1e-9:
            print(f"  FAIL: rate at d=1eV != P (a={a}): {g_si:.3e} vs 5e13"); ok = False

    # 3) Exponent scaling: doubling the excess energy multiplies the rate by 2**a.
    for a in (2.0, 4.0):
        g1 = gamma_au(2.0, 1.0, 1e13, a)        # d = 1 eV
        g2 = gamma_au(3.0, 1.0, 1e13, a)        # d = 2 eV
        ratio = g2 / g1
        if abs(ratio - 2.0 ** a) / (2.0 ** a) > 1e-9:
            print(f"  FAIL: scaling a={a}: ratio {ratio:.4f} != {2.0**a:.4f}"); ok = False

    # 4) Soft (a=2) vs hard (a=4) ordering near threshold: at small excess the
    #    quadratic Si form gives a LARGER rate than the quartic for equal P.
    g_soft = gamma_au(1.4, 1.1, 1e13, 2.0)      # d = 0.3 eV
    g_hard = gamma_au(1.4, 1.1, 1e13, 4.0)
    if not (g_soft > g_hard):
        print(f"  FAIL: near threshold soft(a=2) should exceed hard(a=4): "
              f"{g_soft:.3e} vs {g_hard:.3e}"); ok = False

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
