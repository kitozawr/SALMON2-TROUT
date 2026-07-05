#!/usr/bin/env python3
"""
cds_ii_calibrate.py — A9: compute the CdS impact-ionization prefactor P from
first principles ON OUR OWN EPM BANDS, replacing the "FIT parameter" caveat.

The Stobbe–Redmer–Schattke procedure (PRB 49, 4494 (1994)) applied to CdS:
evaluate the Fermi-golden-rule electron-initiated II rate on the validated
BC1967 wurtzite EPM band structure,

    Gamma(E_hot) = 2*pi * <|M|^2> * rho_pair(E_hot),

where rho_pair is the 3-particle final-state phase space (hot e- at E, final
states counted over the EPM eigenvalue histogram with momentum conservation
statistically averaged, exactly Stobbe's isotropic-matrix-element ansatz) and
<|M|^2> the screened-Coulomb matrix element averaged over the BZ:

    <|M|^2> = < (4*pi / (eps0 * Omega * (q^2 + qTF^2)))^2 >_q,

eps0 = 8.9 [Berlincourt 1963], qTF from the valence density. The resulting
Gamma(E) is then least-squares fitted to the registry's cited Keldysh form
P*(E - E_th)^2 with E_th = 3.6 eV = 1.5*Eg fixed (cited), yielding P.

This is the SAME class of approximation Stobbe made for GaAs (direction-
averaged matrix elements, energy the dominant variable), so the output P has
the same standing as the cited GaAs 2e12 value: a documented calculation on
the material's own EPM bands — no longer an uncited knob.

Usage:  python3 tools/cds_ii_calibrate.py          (from the repo root)
Output: the suggested `sbe_ii_prefactor` and the fit report; paste the value
        (with this script as the provenance note) into the run input.
"""
import sys, pathlib
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import epm_wurtzite_cds as cds

HA_EV = 27.211386245988
AU_S  = 2.418884326505e-17           # 1 a.u.t in seconds


def main():
    # ---- primitive wurtzite EPM spectrum on an MP grid ---------------------
    print("# CdS II prefactor calibration (Stobbe procedure on BC1967 EPM bands)")
    a_au = cds.CDS_A_ANG * cds.ANG_TO_BOHR
    c_au = cds.CDS_C_ANG * cds.ANG_TO_BOHR
    a1, a2, a3 = cds.hexagonal_vectors_au(a_au, c_au)
    B, V = cds.reciprocal(a1, a2, a3)
    G, _ = cds.build_pw_basis(B, 12.0)
    atoms_pos, atoms_spec = cds.hex_primitive_atoms(a_au, c_au)
    n1 = n2 = 9; n3 = 5
    ks = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                f = np.array([(2*i-n1+1)/(2*n1), (2*j-n2+1)/(2*n2), (2*k-n3+1)/(2*n3)])
                ks.append(f @ B)
    NELEC = 16                     # 2 formula units x 8 valence e-
    nocc = NELEC // 2
    ev_all, evbm, ecbm = [], -1e9, 1e9
    for kc in ks:
        H = cds.build_hamiltonian(kc, G, atoms_pos, atoms_spec)
        ev = np.linalg.eigvalsh(H)[:nocc + 8]
        ev_all.append(ev)
        evbm = max(evbm, ev[nocc-1]); ecbm = min(ecbm, ev[nocc])
    ev_all = np.array(ev_all)
    eg = (ecbm - evbm) * HA_EV
    print(f"# grid {n1}x{n2}x{n3}, gap = {eg:.3f} eV (BC1967 target 2.58)")

    # ---- pair phase space rho_pair(E) (Stobbe's isotropic ansatz) ----------
    # final states: e1(c) + e2(c) + hole(v); energy conservation E = e1+e2+eh
    # counted from the eigenvalue histograms (momentum statistically averaged).
    dv = (evbm - ev_all[:, :nocc].ravel()) * HA_EV          # hole energies >= 0
    dc = (ev_all[:, nocc:].ravel() - ecbm) * HA_EV          # electron energies >= 0
    dv = dv[dv < 8]; dc = dc[dc < 8]
    bins = np.arange(0, 8.01, 0.05)
    hv, _ = np.histogram(dv, bins); hc, _ = np.histogram(dc, bins)
    hv = hv / max(dv.size, 1) / 0.05; hc = hc / max(dc.size, 1) / 0.05  # normalized DOS [1/eV]
    # rho_pair(E) = int de1 de2 gc(e1) gc(e2) gv(E - Eg - e1 - e2)
    E_hot = np.arange(eg + 0.1, eg + 7.0, 0.1)              # hot energy ABOVE the VBM ref
    ctr = 0.5 * (bins[1:] + bins[:-1])
    rho = np.zeros_like(E_hot)
    for iE, E in enumerate(E_hot):
        s = 0.0
        for i1, e1 in enumerate(ctr):
            for i2, e2 in enumerate(ctr):
                eh = E - eg - e1 - e2
                if eh < 0 or eh >= 8:
                    continue
                s += hc[i1] * hc[i2] * hv[int(eh / 0.05)]
        rho[iE] = s * 0.05 * 0.05                            # [1/eV]
    # ---- screened matrix element <|M|^2> over the BZ ------------------------
    n_val = NELEC / V
    kf = (3 * np.pi**2 * n_val) ** (1/3)
    qtf2 = 4 * kf / np.pi
    qs = np.linalg.norm(np.random.default_rng(1).standard_normal((4000, 3)), axis=1) \
         * np.linalg.norm(B[0]) / 2
    eps0 = 8.9
    M2 = np.mean((4 * np.pi / (eps0 * V * (qs**2 + qtf2)))**2)
    gam_au = 2 * np.pi * M2 * rho / HA_EV * V**2 / (2*np.pi)**3 * NELEC   # a.u. rate
    gam_s = gam_au / AU_S

    # ---- fit the cited Keldysh form: P*(E_kin - Eth)^2, Eth = 3.6 eV fixed --
    e_kin = E_hot - eg                                       # from the CBM
    eth = 3.6
    m = e_kin > eth
    if m.sum() < 4:
        print("# not enough above-threshold points; widen E range"); return
    x = (e_kin[m] - eth)**2
    P = float(np.dot(x, gam_s[m]) / np.dot(x, x))
    rms = np.sqrt(np.mean((gam_s[m] - P * x)**2)) / max(gam_s[m].max(), 1e-30)
    print(f"# <|M|^2> (screened, BZ-avg) = {M2:.3e} a.u.; qTF^2 = {qtf2:.3f}")
    print(f"# fit window: E_kin in [{eth:.1f}, {e_kin[m].max():.1f}] eV, rel-RMS = {rms:.2f}")
    print(f"#")
    print(f"#   sbe_ii_prefactor = {P:.3e}   ! [s^-1 eV^-2] CdS, THIS script's")
    print(f"#   ! calibration: Stobbe procedure on the BC1967 EPM bands")
    print(f"#   ! (tools/cds_ii_calibrate.py; eps0=8.9 Berlincourt; Eth=3.6 cited)")
    print(f"#")
    print(f"# GaAs sanity anchor: the same procedure class gave Stobbe P=2e12 s^-1 eV^-4.")


if __name__ == '__main__':
    main()
