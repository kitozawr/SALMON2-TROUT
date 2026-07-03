!
!  test_rana_2d.f90  -  graphene 2D Auger / carrier multiplication [R07]
!  (F. Rana, PRB 76, 155431 (2007) / arXiv:0705.1204; wiki/07 sec.6).
!
!  Implements [R07 Eqs. (13),(14),(17)] verbatim (verified against the journal
!  text). Validates BOTH the cited absolute lifetime benchmarks and the
!  structural properties:
!    1) tau_r(n=p=1e12 cm^-2, 300 K, eps_r=10) ~ 1-2 ps [R07 Fig.4; the
!       minority-carrier lifetime at p=1e12 is 1.1 ps] -- window 0.5..3 ps.
!    2) tau_r > 5 ps at 1e11 cm^-2 (300 K) and > 1 ps at 1e12 (77 K) [R07].
!    3) DETAILED BALANCE: at equilibrium (single mu for both branches) the
!       generation (reverse) integral equals recombination exactly
!       [R07: G = R in thermal equilibrium].
!    4) CVVV(n,p) = CCCV(p,n) mirror (and CVVV /= CCCV when n /= p).
!    5) smaller eps => larger rate (eps = 4 SiO2 vs 10 Al2O3 [R07 Fig.5]).
!    6) R grows with density; the GENERATION rate grows with temperature
!       [R07: "the generation rate increases with temperature"].
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

    ! --- (4) monotonicity: R with density; GENERATION with temperature -------
    mu_c = dirac_mu_2d(n11, kT300, V_F_AU)
    R_lo = rana_rcccv(mu_c, -mu_c, kT300, V_F_AU, 10d0, .false.)
    if (R_lo >= R10) call bad('rate not increasing with carrier density')
    R_hot = rana_rcccv(0d0, 0d0, kT300, V_F_AU, 10d0, .true.)
    R_lo  = rana_rcccv(0d0, 0d0, kT77,  V_F_AU, 10d0, .true.)
    if (R_hot <= R_lo) call bad('generation rate not increasing with temperature')

    ! --- (5) CITED absolute lifetime benchmarks [R07] ------------------------
    mu_c = dirac_mu_2d(n12, kT300, V_F_AU);  mu_v = -mu_c
    R_rec = rana_rcccv(mu_c, mu_v, kT300, V_F_AU, 10d0, .false.) &
          + rana_rcccv(-mu_v, -mu_c, kT300, V_F_AU, 10d0, .false.)
    tau_ps = n12 / max(R_rec, 1d-300) * AU_T_FS * 1d-3
    write(*,'(a,f8.3,a)') '  tau_r(1e12 cm^-2, 300 K, eps10) = ', tau_ps, &
        ' ps  [R07 Fig.4 ~1-2 ps, minority lifetime 1.1 ps]'
    if (tau_ps < 0.5d0 .or. tau_ps > 3.0d0) &
        call bad('tau_r(1e12, 300 K) outside the R07 benchmark window 0.5..3 ps')
    mu_c = dirac_mu_2d(n11, kT300, V_F_AU);  mu_v = -mu_c
    R_rec = 2d0 * rana_rcccv(mu_c, mu_v, kT300, V_F_AU, 10d0, .false.)
    tau_ps = n11 / max(R_rec, 1d-300) * AU_T_FS * 1d-3
    write(*,'(a,f8.3,a)') '  tau_r(1e11 cm^-2, 300 K, eps10) = ', tau_ps, ' ps  [cited > 5 ps]'
    if (tau_ps <= 5.0d0) call bad('tau_r(1e11, 300 K) not > 5 ps')
    mu_c = dirac_mu_2d(n12, kT77, V_F_AU);  mu_v = -mu_c
    R_rec = 2d0 * rana_rcccv(mu_c, mu_v, kT77, V_F_AU, 10d0, .false.)
    tau_ps = n12 / max(R_rec, 1d-300) * AU_T_FS * 1d-3
    write(*,'(a,f8.3,a)') '  tau_r(1e12 cm^-2,  77 K, eps10) = ', tau_ps, ' ps  [cited > 1 ps]'
    if (tau_ps <= 1.0d0) call bad('tau_r(1e12, 77 K) not > 1 ps')

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (graphene 2D Rana [R07 Eqs.13/14/17]: cited lifetime benchmarks, '// &
                       'equilibrium detailed balance, CVVV mirror, eps ordering, monotonicity)'
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
