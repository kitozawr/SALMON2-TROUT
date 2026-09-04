#!/usr/bin/env python3
"""Generate the x14 input set: one graphene EPM ground state + the field scan.

    python3 make_inputs.py [--nk 24] [--hw-ev 0.8] [--fields 1,3,10,30,100]
                           [--cycles 8] [--tail-fs 60] [--dt-au 4] [--sigma-ev 0.1]
                           [--temp-k 300] [--eps-r 10]

Writes
    graphene_sit_epm_gs.inp          in-SALMON EPM ground state (43-PW basis, nk x nk x 1)
    rt_E<E>kVcm_<tag>.inp            one RT input per field and per variant:
        coh      coherent SBE only (no dissipation)              -> the Pauli/Rabi bleaching baseline
        diss     ring: e-ph (E2g/A1' + acoustic) + 2D Rana Auger/CM (Markovian, as in x11)
        mem      diss + the graphene 2D collisional-memory analog + dressed reference:
                 yn_sbe_colmem, yn_sbe_colmem_pop (phonon lines for the e-ph sectors, the
                 2D Dirac-plasmon line for the Rana source), yn_sbe_dressed_ref
    rt_dark_<tag>.inp                zero-field controls (the Rana drift toward n_i(T_bath))
    scan_manifest.txt                the field table (kV/cm -> a.u., W/cm^2, A0, Rabi area)

ONE source of truth: the GS and RT inputs share nk / cell / sysname, so run
this script rather than hand-editing copies. Grid guidance (README): nk = 24
is the pipeline smoke test (seconds); the 0.8 eV resonance shell needs
nk >= 150 (>= 3 mesh points per shell radius) -- production on the cluster.
"""
import argparse
import os

AU_EV = 27.211386245988
AU_T_FS = 0.02418884326505
AU_E_VM = 5.14220675e11
Z0_SI = 376.730313668
A_BOHR = 4.648726            # a = 2.46 Ang
LZ_BOHR = 37.794523          # 20 Ang vacuum
VF_AU = 0.44                 # 43-PW EPM slope (~0.96e6 m/s)
SYS = 'graphene_sit'

GS_TEMPLATE = """!########################################################################################!
! x14 graphene PRIMITIVE-CELL EPM ground state (in-SALMON Fortran EPM, theory='epm').        !
!----------------------------------------------------------------------------------------!
! Writes {sys}_k/_eigen/_tm.data read by the rt_*.inp files (generated together by      !
! make_inputs.py -- keep nk consistent).                                                  !
!                                                                                          !
! BASIS (x14 finding, tests/test_graphene_dirac_levels.py): the 7-plane-wave basis of the  !
! x11 input (epm_pw_cutoff_ry = 2.94) is NOT closed under the little group of K and opens  !
! a SPURIOUS 0.21 eV gap at the Dirac point. |G|^2 <= 29.4 a.u. (43 PW, shell-complete)    !
! restores the symmetry-protected degeneracy (gap ~ 1e-5 eV) and meets the thesis          !
! acceptance numbers (v_F = 0.96e6 m/s, Gamma bottom -7.8 eV, M dip -2.7 eV). The SBE       !
! still uses nstate = 2 (the pi/pi* pair), so the larger GS basis costs nothing.           !
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
    nstate = 2                                     ! pi (valence) + pi* (conduction): the Dirac pair
/

&kgrid
    num_kgrid(1:3) = {nk}, {nk}, 1                 ! multiple of 3 -> K = (2/3,1/3) on the mesh
/

&epm
    epm_material     = 'graphene'
    epm_pw_cutoff_ry = 29.4d0                      ! |G|^2 <= 29.4 a.u.: 43 PW, gapless K (NOT 2.94)
    epm_cell         = 'primitive'
/
"""

RT_TEMPLATE = """!########################################################################################!
! x14 graphene self-induced transparency -- {tag:4s} variant, peak E = {ekv:g} kV/cm                 !
!----------------------------------------------------------------------------------------!
! In-plane (x) {shape} pulse, hbar*omega = {hw:.3f} eV, {cycles} cycles (tw) + {tail} fs tail.        !
! I = {iw:.3e} W/cm^2, A0 = E0/omega = {a0:.2e} a.u. (k excursion), Rabi area ~ {area:.2f} rad.   !
! Field before = E_ext column, field after = sheet BC on the Jm column: transmission.py.  !
!{comment}
!########################################################################################!

&calculation
  theory = 'sbe'
/

&control
  sysname = '{sys}'
/

&units
  unit_system = 'au'
/

&system
  yn_periodic  = 'y'
  al_vec1(1:3) = {a:.6f}d0, 0.0d0, 0.0d0
  al_vec2(1:3) = {a2:.6f}d0, {a3:.6f}d0, 0.0d0
  al_vec3(1:3) = 0.0d0, 0.0d0, {lz:.6f}d0          ! 20 Ang vacuum: L_z of the sheet cell
  nelec  = 2
  nstate = 2
/

&kgrid
  num_kgrid(1:3) = {nk}, {nk}, 1
/

&tgrid
  dt = {dt:.4f}d0
  nt = {nt}
/

&emfield
  ae_shape1      = '{shape}'
  epdir_re1(1:3) = 1.0d0, 0.0d0, 0.0d0
  E_amplitude1   = {eau:.6e}
  tw1            = {tw:.3f}d0
  omega1         = {w:.8f}d0
/

&analysis
  out_projection_k_step = {nt}
/

&epm
  epm_material = 'graphene'
/

&sbe
{sbe}/
"""

SBE_BLOCKS = {
    'coh': "  ! coherent SBE only: Pauli-blocking / Rabi bleaching baseline (no ring)\n",
    'diss': (
        "  yn_sbe_superres       = 'y'    ! ring (mandatory for both graphene channels)\n"
        "  yn_sbe_eph            = 'y'    ! E2g(196 meV) + A1'(160 meV) Kohn-anomaly modes\n"
        "  yn_sbe_eph_acoustic   = 'y'    ! D = 16 eV [Hwang-Das Sarma], TF-screened\n"
        "  yn_sbe_auger          = 'y'    ! 2D Rana Auger / carrier multiplication [R07]\n"
        "  sbe_eph_temperature_k = {temp:.1f}d0\n"
        "  sbe_coulomb_epsilon   = {eps:.1f}d0   ! substrate eps_r of the Rana rates (10 = R07 benchmark)\n"
        "  sbe_search_sigma_e_ev = {sigma:.3f}d0 ! ring energy-matching width (mesh-matched)\n"
    ),
    'mem': (
        "  yn_sbe_superres       = 'y'\n"
        "  yn_sbe_eph            = 'y'\n"
        "  yn_sbe_eph_acoustic   = 'y'\n"
        "  yn_sbe_auger          = 'y'\n"
        "  sbe_eph_temperature_k = {temp:.1f}d0\n"
        "  sbe_coulomb_epsilon   = {eps:.1f}d0\n"
        "  sbe_search_sigma_e_ev = {sigma:.3f}d0\n"
        "  ! ---- graphene 2D collisional-memory analog (wiki/10 sec. 8.11) ----\n"
        "  yn_sbe_colmem         = 'y'    ! e-ph gout memory: graphene phonon lines\n"
        "  yn_sbe_colmem_pop     = 'y'    ! e-ph source filtered with the phonon lines;\n"
        "                                 ! Rana (Coulomb) source with the 2D Dirac-plasmon line\n"
        "  yn_sbe_dressed_ref    = 'y'    ! Option A: removes the Dirac-point rotation background\n"
    ),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--nk', type=int, default=24)
    ap.add_argument('--hw-ev', type=float, default=0.8)
    ap.add_argument('--fields', default='1,3,10,30,100', help='peak fields [kV/cm]')
    ap.add_argument('--cycles', type=float, default=8.0, help='tw in optical cycles (Acos2 envelope)')
    ap.add_argument('--tail-fs', type=float, default=60.0)
    ap.add_argument('--dt-au', type=float, default=4.0)
    ap.add_argument('--sigma-ev', type=float, default=0.1)
    ap.add_argument('--temp-k', type=float, default=300.0)
    ap.add_argument('--eps-r', type=float, default=10.0)
    ap.add_argument('--variants', default='coh,diss,mem')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()
    if args.nk % 3:
        print(f'WARNING: nk = {args.nk} is not a multiple of 3 -> the K point is off the mesh')

    w = args.hw_ev / AU_EV
    period = 2 * 3.141592653589793 / w
    tw = args.cycles * period
    nt = int(round((tw + args.tail_fs / AU_T_FS) / args.dt_au))
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, f'{SYS}_epm_gs.inp'), 'w') as fh:
        fh.write(GS_TEMPLATE.format(sys=SYS, nk=args.nk))

    fields = [float(x) for x in args.fields.split(',')]
    variants = args.variants.split(',')
    man = [f'# x14 field scan: nk={args.nk}, hw={args.hw_ev} eV, tw={tw * AU_T_FS:.1f} fs ({args.cycles} cycles), '
           f'tail={args.tail_fs} fs, dt={args.dt_au} a.u., nt={nt}',
           f'{"E[kV/cm]":>9} {"E[a.u.]":>11} {"I[W/cm2]":>10} {"A0[a.u.]":>10} {"Rabi area":>9}']
    for ekv in fields + [0.0]:
        e_vm = ekv * 1e5
        eau = e_vm / AU_E_VM
        iw = e_vm**2 / (2 * Z0_SI) / 1e4
        a0 = eau / w
        area = a0 * VF_AU * tw / 2          # Omega_R ~ A0 v_F (Dirac dipole), half-envelope
        if ekv > 0:
            man.append(f'{ekv:9.3f} {eau:11.4e} {iw:10.3e} {a0:10.3e} {area:9.3f}')
        for tag in variants:
            sbe = SBE_BLOCKS[tag].format(temp=args.temp_k, eps=args.eps_r, sigma=args.sigma_ev)
            name = f'rt_E{ekv:g}kVcm_{tag}.inp' if ekv > 0 else f'rt_dark_{tag}.inp'
            if ekv == 0 and tag == 'coh':
                continue          # a dark coherent run is trivially static
            comment = (' zero-field CONTROL: the Rana channel relaxes the T=0 initial state toward the'
                       ' n_i(T_bath) thermal density -- subtract this drift from the driven runs.  '
                       if ekv == 0 else
                       ' Compare coh / diss / mem at the same field: the memory+dressed-reference'
                       ' variant is the one whose ring sees only REAL carriers (wiki/10 sec. 8.11).   ')
            with open(os.path.join(args.outdir, name), 'w') as fh:
                fh.write(RT_TEMPLATE.format(
                    sys=SYS, tag=tag, ekv=ekv, shape='Acos2', hw=args.hw_ev, cycles=args.cycles,
                    tail=args.tail_fs, iw=iw, a0=a0, area=area, comment=comment,
                    a=A_BOHR, a2=A_BOHR / 2, a3=A_BOHR * 3**0.5 / 2, lz=LZ_BOHR, nk=args.nk,
                    dt=args.dt_au, nt=nt, eau=eau, tw=tw, w=w, sbe=sbe))
    with open(os.path.join(args.outdir, 'scan_manifest.txt'), 'w') as fh:
        fh.write('\n'.join(man) + '\n')
    print('\n'.join(man))
    print(f'wrote {SYS}_epm_gs.inp + {len(fields) * len(variants) + len(variants) - 1} rt_*.inp in {args.outdir}')


if __name__ == '__main__':
    main()
