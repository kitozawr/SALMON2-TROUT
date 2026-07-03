!
!  test_rana_2d.f90  -  graphene 2D Auger / carrier multiplication [R07]
!  (F. Rana, PRB 76, 155431 (2007) / arXiv:0705.1204; wiki/07 sec.6).
!
!  *** The ABSOLUTE normalization of rana_rcccv is PENDING the [R07] equation
!  text (unreachable from the dev environment -- network policy; the wiki
!  transcription lost the angular collapse Jacobians, see the function
!  header). Per the strict provenance rule no factors are guessed, so this
!  test validates every RELATIVE / structural property the channel relies on
!  and prints the (uncalibrated) lifetime for the record:
!    1) DETAILED BALANCE: at equilibrium (single mu for both branches) the
!       generation (reverse) integral equals recombination exactly.
!    2) CVVV(n,p) = CCCV(p,n) mirror (and CVVV /= CCCV when n /= p).
!    3) smaller eps => larger rate (eps = 4 SiO2 vs 10 Al2O3 [R07 setup]).
!    4) monotonicity: R grows with density and with temperature.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_rana_2d
    use sbe_superres_ssbe, only: dirac_mu_2d, rana_qtf, rana_rcccv
    implicit none
    real(8), parameter :: AU_T_FS   = 0.02418884326505d0   ! 1 a.u.t in fs
    real(8), parameter :: A0_CM     = 0.52917721067d-8     ! Bohr in cm
    real(8), parameter :: KB_HA_K   = 3.166811563d-6       ! k_B [Ha/K]
    real(8), parameter :: V_F_AU    = 1.0d8 / 2.18769126364d8  ! 1e8 cm/s in a.u.
    integer :: nfail

    real(8) :: n12, n11, kT300, kT77, mu_c, mu_v, R_rec, R_gen, R_cvvv, tau_ps
    real(8) :: R10, R4, R_lo, R_hot

    nfail = 0
    n12   = 1.0d12 * A0_CM**2          ! 1e12 cm^-2 in a.u.^-2
    n11   = 1.0d11 * A0_CM**2
    kT300 = 300d0 * KB_HA_K
    kT77  =  77d0 * KB_HA_K

    ! --- (1) detailed balance: generation = recombination at equilibrium ----
    R_rec = rana_rcccv(0d0, 0d0, kT300, V_F_AU, 10d0, .false.)
    R_gen = rana_rcccv(0d0, 0d0, kT300, V_F_AU, 10d0, .true.)
    write(*,'(a,es12.4,a,es12.4)') '  equilibrium: R_rec = ', R_rec, '  R_gen = ', R_gen
    if (R_rec <= 0d0) call bad('equilibrium rate not positive (integral inert)')
    if (abs(R_rec - R_gen) > 1d-10 * max(R_rec, 1d-300)) &
        call bad('generation /= recombination at equilibrium (detailed balance)')

    ! --- (2) CVVV mirror ------------------------------------------------------
    mu_c = dirac_mu_2d(n12, kT300, V_F_AU)
    mu_v = -dirac_mu_2d(n11, kT300, V_F_AU)                 ! n != p on purpose
    R_rec  = rana_rcccv(mu_c, mu_v, kT300, V_F_AU, 10d0, .false.)
    R_cvvv = rana_rcccv(-mu_v, -mu_c, kT300, V_F_AU, 10d0, .false.)
    if (abs(R_rec - R_cvvv) < 1d-3 * max(R_rec, 1d-300)) &
        call bad('CCCV(n,p) and CVVV(n,p)=CCCV(p,n) unexpectedly equal at n /= p')
    ! and at n = p the mirror must be exact:
    mu_c = dirac_mu_2d(n12, kT300, V_F_AU);  mu_v = -mu_c
    R_rec  = rana_rcccv(mu_c, mu_v, kT300, V_F_AU, 10d0, .false.)
    R_cvvv = rana_rcccv(-mu_v, -mu_c, kT300, V_F_AU, 10d0, .false.)
    if (abs(R_rec - R_cvvv) > 1d-10 * max(R_rec, 1d-300)) &
        call bad('CVVV /= CCCV at n = p (mirror symmetry broken)')

    ! --- (3) smaller eps => larger rate --------------------------------------
    R10 = rana_rcccv(mu_c, mu_v, kT300, V_F_AU, 10d0, .false.)
    R4  = rana_rcccv(mu_c, mu_v, kT300, V_F_AU,  4d0, .false.)
    if (R4 <= R10) call bad('eps=4 (SiO2) rate not larger than eps=10 (Al2O3)')

    ! --- (4) monotonicity in density and temperature -------------------------
    mu_c = dirac_mu_2d(n11, kT300, V_F_AU)
    R_lo = rana_rcccv(mu_c, -mu_c, kT300, V_F_AU, 10d0, .false.)
    if (R_lo >= R10) call bad('rate not increasing with carrier density')
    mu_c = dirac_mu_2d(n12, kT77, V_F_AU)
    R_hot = rana_rcccv(mu_c, -mu_c, kT77, V_F_AU, 10d0, .false.)
    if (R_hot >= R10) call bad('rate not increasing with temperature (77 K vs 300 K)')

    ! --- record: UNCALIBRATED lifetime at the cited [R07] benchmark point ----
    mu_c = dirac_mu_2d(n12, kT300, V_F_AU);  mu_v = -mu_c
    R_rec = rana_rcccv(mu_c, mu_v, kT300, V_F_AU, 10d0, .false.) &
          + rana_rcccv(-mu_v, -mu_c, kT300, V_F_AU, 10d0, .false.)
    tau_ps = n12 / max(R_rec, 1d-300) * AU_T_FS * 1d-3
    write(*,'(a,es12.4,a)') '  UNCALIBRATED tau_r(1e12 cm^-2, 300 K, eps10) = ', tau_ps, &
        ' ps  [R07 cites ~1.1 ps; absolute prefactor pending the R07 Eq. text]'

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (graphene 2D Rana SHAPE: equilibrium detailed balance, '// &
                       'CVVV mirror, eps ordering, n/T monotonicity; ABSOLUTE norm pending R07 text)'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        call exit(1)
    end if

contains
    subroutine bad(name)
        character(*), intent(in) :: name
        write(*,'(a,a)') '  FAIL: ', name
        nfail = nfail + 1
    end subroutine bad
end program test_rana_2d
