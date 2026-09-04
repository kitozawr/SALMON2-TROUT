#!/usr/bin/env python3
"""Generate the x14 input set: graphene EPM ground state + the field scan.

    python3 make_inputs.py [--nk 147] [--field dast|acos2] [--fields 1,3,10,30,100]
                           [--variants coh,diss,mem] [--dt-fs 0.1] [--tail-fs 100]
                           [--no-sheet] [--temp-k 300] [--eps-r 10] [--sigma-ev 0.1]
                           [--dast-file PATH] [--hw-ev 0.8] [--cycles 8] [--outdir DIR]

DRIVE (default `dast`): the maintainer's DAST optical-rectification THz transient,
single-cycle proxy `DAST_singlecycle_100kV.txt` (exercise x12: A(t) Gaussian,
E = -dA/dt one clean cycle at 3.36 THz, 285 fs), x-polarised = in-plane. It is
rescaled to each requested PEAK FIELD (the file's own peak |dA/dt| is measured,
so "100 kV/cm" means exactly 100 kV/cm), the initial offset is removed and the
ends are cos^2-windowed (10 fs) so A starts/ends at exactly 0 (the SALMON file
reader returns 0 outside the file: an offset would be a delta-spike in E).
`acos2`: the analytic near-IR pulse of the first x14 version (hbar*omega, cycles).

UNITS: everything in A_eV_fs (the field files are A[fs*V/Ang] vs t[fs]).

VARIANTS
    coh   coherent SBE only (no dissipation)
    diss  ring: e-ph (E2g/A1' + acoustic) + 2D Rana Auger/CM at the LATTICE T (Markovian; as x11)
    mem   diss + the graphene 2D collisional-memory analog (yn_sbe_colmem, yn_sbe_colmem_pop)
          + Option-A dressed reference + the TWO-TEMPERATURE Coulomb sector (yn_sbe_rana_te)
All variants carry the 2D-sheet SELF-CONSISTENT field (yn_sbe_sheet_field) unless
--no-sheet: then E_tot/Ac_tot in *_sbe_rt.data are the transmitted field.
`rt_dark_*` are zero-field controls.

GRID: default nk = 147 -- dense (cluster) AND an odd multiple of 3, so the Dirac
point K = (2/3,1/3) is ON the half-shifted MP mesh (the THz physics lives at K);
nk = 24 is the pipeline smoke test.
"""
import argparse
import os
import sys

import numpy as np

AU_EV = 27.211386245988
AU_T_FS = 0.02418884326505
AU_E_VM = 5.14220675e11
AU_E_VANG = AU_E_VM * 1e-10          # 1 a.u. field = 51.422 V/Ang
Z0_SI = 376.730313668
VF_AU = 0.44
SYS = 'graphene_sit'
HERE = os.path.dirname(os.path.abspath(__file__))
DAST_DEFAULT = os.path.join(HERE, '..', 'exercise_x12_Si_frozen_phonon_indirect', 'DAST_singlecycle_100kV.txt')

GS_TEMPLATE = """!########################################################################################!
! x14 graphene PRIMITIVE-CELL EPM ground state (in-SALMON Fortran EPM, theory='epm').        !
!----------------------------------------------------------------------------------------!
! Writes {sys}_k/_eigen/_tm.data read by the rt_*.inp files (generated together by      !
! make_inputs.py -- keep nk consistent). nk = {nk}: dense (cluster) and an ODD multiple of 3 !
! so the Dirac point K = (2/3,1/3) is on the half-shifted Monkhorst-Pack mesh.             !
!                                                                                          !
! BASIS (x14 finding, tests/test_graphene_dirac_levels.py): the 7-plane-wave basis of the  !
! old x11 input (epm_pw_cutoff_ry = 2.94) is NOT closed under the little group of K and    !
! opens a SPURIOUS 0.21 eV gap at the Dirac point. |G|^2 <= 29.4 a.u. (43 PW, shell-        !
! complete) restores the symmetry-protected degeneracy (gap 0) and meets the thesis        !
! acceptance numbers (v_F = 0.96e6 m/s, Gamma bottom -7.8 eV, M dip -2.7 eV).               !
!########################################################################################!

&calculation
    theory = 'epm'
/

&control
    sysname = '{sys}'
/

&units
    unit_system = 'A_eV_fs'
/

&system
    yn_periodic = 'y'
    al(1:3) = 2.46000197d0, 2.46000197d0, 19.99999867d0   ! informational: the EPM fixes a = 2.46 Ang
    nelec  = 2                                     ! 2 atoms x 1 pi e- -> 1 filled valence pi band
    nstate = {nstate}                                    ! bands in the GS/SBE basis (2 = pi/pi* only; see README on the velocity-gauge f-sum rule)
/

&kgrid
    num_kgrid(1:3) = {nk}, {nk}, 1
/

&epm
    epm_material     = 'graphene'
    epm_pw_cutoff_ry = 29.4d0                      ! |G|^2 <= 29.4 a.u.: 43 PW, gapless K (NOT 2.94)
    epm_cell         = 'primitive'
/
"""

RT_TEMPLATE = """!########################################################################################!
! x14 graphene self-induced transparency -- {tag:4s} variant, peak E = {ekv:g} kV/cm
!----------------------------------------------------------------------------------------!
! {drive_line}
! I_peak = {iw:.3e} W/cm^2, A0 = {a0:.2e} a.u. (peak k excursion) -- {regime}
! Field BEFORE = E_ext column; field AFTER = E_tot column (2D-sheet self-consistent field,
! yn_sbe_sheet_field) -> transmission.py.  {comment}
!########################################################################################!

&calculation
  theory = 'sbe'
/

&control
  sysname = '{sys}'
/

&units
  unit_system = 'A_eV_fs'
/

&system
  yn_periodic  = 'y'
  al_vec1(1:3) = 2.46000197d0, 0.00000000d0,  0.00000000d0
  al_vec2(1:3) = 1.22999834d0, 2.13042512d0,  0.00000000d0
  al_vec3(1:3) = 0.00000000d0, 0.00000000d0, 19.99999867d0     ! 20 Ang vacuum: L_z of the sheet cell
  nelec  = 2
  nstate = {nstate}
/

&kgrid
  num_kgrid(1:3) = {nk}, {nk}, 1
/

&tgrid
  dt = {dt:.4f}d0                     ! fs
  nt = {nt}
/

&emfield
{emfield}/

&analysis
  out_projection_k_step = {proj}   ! k-resolved level-population snapshots every {snap_fs:g} fs (plot_levels.py)
/

&epm
  epm_material = 'graphene'
/

&sbe
{sbe}/
"""

EM_DAST = """  ae_shape1   = 'input'
  file_input1 = '{fname}'             ! DAST single-cycle 3.36 THz proxy, scaled to {ekv:g} kV/cm peak
  epdir_re1(1:3) = {epdir}     ! in-plane polarisation ({pol}); the THz beam is at normal incidence (k || z)
"""
EM_ACOS2 = """  ae_shape1      = 'Acos2'
  epdir_re1(1:3) = {epdir}     ! in-plane polarisation ({pol}); the THz beam is at normal incidence (k || z)
  E_amplitude1   = {evang:.6e}          ! V/Ang  (= {ekv:g} kV/cm)
  tw1            = {tw:.3f}d0            ! fs ({cycles:g} cycles)
  omega1         = {hw:.6f}d0            ! eV
"""
EM_DARK = """  ae_shape1 = 'none'                  ! zero-field control
"""

SBE_COMMON_RING = (
    "  yn_sbe_superres       = 'y'    ! ring (mandatory for both graphene channels)\n"
    "  yn_sbe_eph            = 'y'    ! E2g(196 meV) + A1'(160 meV) Kohn-anomaly modes\n"
    "  yn_sbe_eph_acoustic   = 'y'    ! D = 16 eV [Hwang-Das Sarma], TF-screened\n"
    "  yn_sbe_auger          = 'y'    ! 2D Rana Auger / carrier multiplication [R07]\n"
    "  sbe_eph_temperature_k = {temp:.1f}d0  ! lattice / phonon-bath temperature\n"
    "  sbe_coulomb_epsilon   = {eps:.1f}d0   ! substrate eps_r of the Rana rates (10 = R07 benchmark)\n"
    "  sbe_search_sigma_e_ev = {sigma:.3f}d0 ! ring energy-matching width\n"
)
SBE_MEM = (
    "  ! ---- graphene 2D collisional-memory analog (wiki/10 sec. 8.11 / wiki/12) ----\n"
    "  yn_sbe_colmem         = 'y'    ! e-ph gout memory: graphene phonon lines\n"
    "  yn_sbe_colmem_pop     = 'y'    ! e-ph source: phonon lines; Rana source: 2D Dirac-plasmon line\n"
    "  yn_sbe_dressed_ref    = 'y'    ! Option A: removes the Dirac-point rotation background\n"
    "  yn_sbe_rana_te        = 'y'    ! two-temperature Coulomb sector: R-G, Q_TF, plasmon at T_e\n"
    "                                 ! (T_e from the cone moments; lattice cools via e-ph) -> *_sbe_te.data\n"
)
SBE_SHEET = "  yn_sbe_sheet_field    = 'y'    ! 2D-sheet self-consistent (radiation-reaction) field\n"
SBE_SUMRULE = ("  yn_sbe_vg_sumrule     = 'y'    ! velocity-gauge pure-gauge restoration: J(t) -= adiabatic ground-\n"
               "                                 ! state current of the same truncated H_k(A(t)) (parameter-free;\n"
               "                                 ! = eta N_e A/V at linear order, exact beyond; wiki/12 sec. 6a)\n")
SBE_NLAYERS = "  sbe_sheet_nlayers     = {nl}      ! {nl} identical decoupled sheets in the same local field (incoherent bilayer)\n"
SBE_WINDOW = ("  frozen_core_threshold_ev = -15.0d0 ! active (dissipated) window around E_F(Gamma): pi, pi* (+ the\n"
              "  frozen_free_threshold_ev =  14.0d0 ! state degenerate with pi* at Gamma); the unitary runs on all nstate\n")


def load_dast(path):
    d = np.loadtxt(path, comments='#')
    return d[:, 0], d[:, 1:4]


def scaled_dast(t, A, e_target_kvcm, window_fs=10.0):
    """Rescale the file so that peak |E| = e_target [kV/cm]; remove the initial
    offset; cos^2-window the first/last `window_fs` so A(0) = A(end) = 0."""
    A = A - A[0]
    E = -np.gradient(A, t, axis=0)                     # V/Ang (file units: A in fs*V/Ang)
    e_peak_file = np.max(np.abs(E)) * 1e10 / 1e5       # -> kV/cm
    A = A * (e_target_kvcm / e_peak_file)
    w = np.ones_like(t)
    m0 = t < t[0] + window_fs
    m1 = t > t[-1] - window_fs
    w[m0] = np.sin(0.5 * np.pi * (t[m0] - t[0]) / window_fs)**2
    w[m1] = np.sin(0.5 * np.pi * (t[-1] - t[m1]) / window_fs)**2
    A = A * w[:, None]
    E = -np.gradient(A, t, axis=0)
    return A, np.max(np.abs(E)) * 1e10 / 1e5, np.max(np.abs(A)) / (AU_E_VANG * AU_T_FS)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--nk', type=int, default=147)
    ap.add_argument('--nstate', type=int, default=2, help='bands in the basis (2 = production with the pure-gauge restoration; static eta 2/4/8/16 -> 0.30/0.10/0.036/0.030 is removed exactly at any nstate; nstate >= 8 needs --dt-fs 0.05: the 34-90 eV bands leak population in the S4 step at 0.1 fs; see README sec. 7)')
    ap.add_argument('--no-sumrule', action='store_true', help='switch the VG pure-gauge restoration of the current off')
    ap.add_argument('--n-layers', type=int, default=1, help='identical electronically decoupled sheets in the same local field (incoherent/twisted bilayer = 2)')
    ap.add_argument('--snap-fs', type=float, default=50.0, help='interval of the k-resolved level-population snapshots [fs] (0 = final only)')
    ap.add_argument('--field', choices=('dast', 'acos2'), default='dast')
    ap.add_argument('--dast-file', default=DAST_DEFAULT)
    ap.add_argument('--fields', default='1,3,10,30,100', help='peak fields [kV/cm]')
    ap.add_argument('--variants', default='coh,diss,mem')
    ap.add_argument('--dt-fs', type=float, default=0.1)
    ap.add_argument('--tail-fs', type=float, default=100.0, help='field-free evolution after the pulse')
    ap.add_argument('--no-sheet', action='store_true', help='switch the self-consistent sheet field off')
    ap.add_argument('--temp-k', type=float, default=300.0)
    ap.add_argument('--eps-r', type=float, default=10.0)
    ap.add_argument('--sigma-ev', type=float, default=0.1)
    ap.add_argument('--hw-ev', type=float, default=0.8, help='acos2 photon energy')
    ap.add_argument('--cycles', type=float, default=8.0, help='acos2 tw in cycles')
    ap.add_argument('--pol', choices=('x', 'y'), default='x', help='in-plane polarisation: x = along a1 (zigzag), y = perpendicular (armchair)')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()
    if args.nk % 3:
        print(f'WARNING: nk = {args.nk} is not a multiple of 3', file=sys.stderr)
    if args.nk % 2 == 0:
        print(f'NOTE: nk = {args.nk} is even -> K is NOT on the half-shifted MP mesh (odd multiples of 3 put it on)',
              file=sys.stderr)
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, f'{SYS}_epm_gs.inp'), 'w') as fh:
        fh.write(GS_TEMPLATE.format(sys=SYS, nk=args.nk, nstate=args.nstate))

    fields = [float(x) for x in args.fields.split(',')]
    epdir = '1.0d0, 0.0d0, 0.0d0' if args.pol == 'x' else '0.0d0, 1.0d0, 0.0d0'
    variants = args.variants.split(',')
    sheet = ('' if args.no_sheet else SBE_SHEET + (SBE_NLAYERS.format(nl=args.n_layers) if args.n_layers > 1 else '')) \
        + ('' if args.no_sumrule else SBE_SUMRULE) + (SBE_WINDOW if args.nstate > 2 else '')
    man = []
    if args.field == 'dast':
        t, A = load_dast(args.dast_file)
        t_pulse = t[-1] - t[0]
        nt = int(round((t_pulse + args.tail_fs) / args.dt_fs))
        man.append(f'# x14 scan: DAST single-cycle 3.36 THz proxy ({os.path.basename(args.dast_file)}, {t_pulse:.1f} fs)'
                   f' + {args.tail_fs:g} fs tail; nk={args.nk}, dt={args.dt_fs} fs, nt={nt}; sheet field {"off" if args.no_sheet else "ON"}')
    else:
        w_au = args.hw_ev / AU_EV
        tw_fs = args.cycles * 2 * np.pi / w_au * AU_T_FS
        nt = int(round((tw_fs + args.tail_fs) / args.dt_fs))
        man.append(f'# x14 scan: Acos2 hw={args.hw_ev} eV, {args.cycles:g} cycles (tw={tw_fs:.1f} fs) + {args.tail_fs:g} fs tail;'
                   f' nk={args.nk}, dt={args.dt_fs} fs, nt={nt}; sheet field {"off" if args.no_sheet else "ON"}')
    man.append(f'{"E[kV/cm]":>9} {"I[W/cm2]":>10} {"A0[a.u.]":>10}  file/drive')

    for ekv in fields + [0.0]:
        e_vm = ekv * 1e5
        iw = e_vm**2 / (2 * Z0_SI) / 1e4
        if ekv > 0:
            if args.field == 'dast':
                As, e_chk, a0 = scaled_dast(t, A, ekv)
                fname = f'DAST_E{ekv:g}kVcm.txt'
                # the 'input' field file carries the three Cartesian components itself
                # (epdir is not applied to it): put the waveform into the --pol column
                As_pol = As.copy()
                if args.pol == 'y':
                    As_pol[:, 1] = As[:, 0]; As_pol[:, 0] = 0.0
                np.savetxt(os.path.join(args.outdir, fname), np.column_stack([t, As_pol]), fmt='%.6f %.10e %.10e %.10e',
                           header=f'DAST single-cycle 3.36 THz proxy rescaled to peak |E| = {e_chk:.3f} kV/cm; '
                                  f't[fs] Ax Ay Az [fs*V/Ang]; offset removed, 10 fs cos^2 end windows')
                emfield = EM_DAST.format(fname=fname, ekv=ekv, epdir=epdir, pol=args.pol)
                drive_line = f'DAST single-cycle 3.36 THz proxy (x12), rescaled to peak |E| = {e_chk:.2f} kV/cm ({fname}).'
                man.append(f'{ekv:9.3f} {iw:10.3e} {a0:10.3e}  {fname} (peak {e_chk:.2f} kV/cm)')
            else:
                a0 = (e_vm / AU_E_VM) / (args.hw_ev / AU_EV)
                emfield = EM_ACOS2.format(evang=e_vm * 1e-10, ekv=ekv, tw=tw_fs, cycles=args.cycles, hw=args.hw_ev, epdir=epdir, pol=args.pol)
                drive_line = f'Acos2, hbar*omega = {args.hw_ev} eV, {args.cycles:g} cycles.'
                man.append(f'{ekv:9.3f} {iw:10.3e} {a0:10.3e}  Acos2 {args.hw_ev} eV')
        for tag in variants:
            if ekv == 0 and tag == 'coh':
                continue
            if tag == 'coh':
                sbe = "  ! coherent SBE only (no dissipation)\n" + sheet
            elif tag == 'diss':
                sbe = SBE_COMMON_RING.format(temp=args.temp_k, eps=args.eps_r, sigma=args.sigma_ev) + sheet
            elif tag == 'mem':
                sbe = SBE_COMMON_RING.format(temp=args.temp_k, eps=args.eps_r, sigma=args.sigma_ev) + SBE_MEM + sheet
            else:
                sys.exit(f'unknown variant {tag}')
            if ekv == 0:
                emfield = EM_DARK
                drive_line = 'zero-field CONTROL (Rana drift, T_e = bath).'
                a0 = 0.0
            regime = ('THz: A0 >> mesh spacing -> Zener/Landau-Zener pair creation at the Dirac point'
                      if args.field == 'dast' and a0 > 3.0 / args.nk else
                      'perturbative interband (A0 << mesh spacing)')
            name = f'rt_E{ekv:g}kVcm_{tag}.inp' if ekv > 0 else f'rt_dark_{tag}.inp'
            comment = ('Compare coh / diss / mem at the same field.' if ekv > 0 else
                       'Subtract this drift from the driven runs.')
            proj = nt if args.snap_fs <= 0 else max(1, int(round(args.snap_fs / args.dt_fs)))
            with open(os.path.join(args.outdir, name), 'w') as fh:
                fh.write(RT_TEMPLATE.format(sys=SYS, tag=tag, ekv=ekv, drive_line=drive_line, iw=iw, a0=a0, proj=proj, snap_fs=args.snap_fs,
                                            regime=regime, comment=comment, nk=args.nk, nstate=args.nstate,
                                            dt=args.dt_fs, nt=nt, emfield=emfield, sbe=sbe))
    with open(os.path.join(args.outdir, 'scan_manifest.txt'), 'w') as fh:
        fh.write('\n'.join(man) + '\n')
    print('\n'.join(man))
    print(f'wrote {SYS}_epm_gs.inp + rt_*.inp (+ field files) in {args.outdir}')


if __name__ == '__main__':
    main()
