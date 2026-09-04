!
!  test_dirac_te_fit.f90 - two-temperature model of the graphene Coulomb sector
!  (wiki/12 sec. 3): the carrier temperature T_e and the quasi-Fermi levels are
!  READ from the moments (n, p, energy) of the Dirac-cone populations
!  (dirac_fit_te); the lattice stays at the phonon-bath T and the e-ph kinetics
!  cool the carriers. This test pins the moment integrals and the inversion.
!
!  Checks:
!    1) closed form: eps(mu=0, T) = (2/pi)(3 zeta(3)/2)(kT)^3/v^2 to 1e-3.
!    2) EXPLICIT 2D k-mesh moments (g = 4, Cartesian mesh around K) of thermal
!       electron/hole populations at (mu_c, mu_h, T) agree with dirac_n_2d /
!       dirac_e_2d to 0.5 %, and dirac_fit_te recovers T within 1 % and the mu's
!       within 1 % of kT (a non-degenerate 1000 K plasma; an intrinsic 2000 K one).
!    3) degenerate, cold: E_F = 0.3 eV at 300 K (kT/E_F = 0.09, the thermal part of
!       the energy is ~2 %): T recovered within 5 %, mu within 1 %.
!    4) no carriers -> kT_e = kT_min (bath); energy at/below the bath value -> bath.
!    5) at fixed n the fitted T_e is monotone in the energy (hotter -> larger T_e).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_dirac_te_fit
    use sbe_superres_ssbe, only: dirac_n_2d, dirac_e_2d, dirac_mu_2d, dirac_fit_te
    implicit none
    real(8), parameter :: PI      = 3.14159265358979323846d0
    real(8), parameter :: KB_HA_K = 3.166811563d-6
    real(8), parameter :: AU_EV   = 27.211386245988d0
    real(8), parameter :: V_F_AU  = 1.0d8 / 2.18769126364d8
    real(8), parameter :: ZETA3   = 1.2020569031595943d0
    integer :: nfail
    real(8) :: kT, e_num, e_ref, kTe, muc, muh, n, p, eps, kTmin, kTmax, e1, e2, kTe1, kTe2

    nfail = 0
    kTmin = 300d0 * KB_HA_K;  kTmax = 2d0

    ! --- (1) intrinsic closed form -------------------------------------------
    kT = 1000d0 * KB_HA_K
    e_num = dirac_e_2d(0d0, kT, V_F_AU)
    e_ref = (2d0 / PI) * (1.5d0 * ZETA3) * kT**3 / V_F_AU**2
    if (abs(e_num - e_ref) > 1d-3 * e_ref) call bad('dirac_e_2d(mu=0) /= (2/pi)(3 zeta3/2)(kT)^3/v^2')

    ! --- (2) explicit mesh moments + fit ---------------------------------------
    call mesh_case(1000d0, 0.05d0 / AU_EV, 0.03d0 / AU_EV, 'non-degenerate 1000 K')
    call mesh_case(2000d0, 0d0, 0d0, 'intrinsic 2000 K')

    ! --- (3) degenerate cold (analytic moments) ---------------------------------
    kT = 300d0 * KB_HA_K
    n = dirac_n_2d(0.3d0 / AU_EV, kT, V_F_AU);  p = dirac_n_2d(0d0, kT, V_F_AU)
    eps = dirac_e_2d(0.3d0 / AU_EV, kT, V_F_AU) + dirac_e_2d(0d0, kT, V_F_AU)
    call dirac_fit_te(n, p, eps, V_F_AU, 0.5d0 * kTmin, kTmax, kTe, muc, muh)
    write(*,'(a,f8.1,a,f8.4,a)') '  degenerate 300 K, E_F = 0.3 eV: T_e = ', kTe / KB_HA_K, &
        ' K, mu_c = ', muc * AU_EV, ' eV'
    if (abs(kTe - kT) > 5d-2 * kT) call bad('degenerate case: T_e off by > 5 %')
    if (abs(muc - 0.3d0 / AU_EV) > 1d-2 * kT) call bad('degenerate case: mu_c off by > 1 % of kT')

    ! --- (4) fallbacks ------------------------------------------------------------
    call dirac_fit_te(0d0, 0d0, 0d0, V_F_AU, kTmin, kTmax, kTe, muc, muh)
    if (abs(kTe - kTmin) > 1d-14) call bad('empty cone must return the bath temperature')
    n = dirac_n_2d(0d0, kTmin, V_F_AU)
    eps = 0.5d0 * (dirac_e_2d(0d0, kTmin, V_F_AU) * 2d0)          ! colder than the bath
    call dirac_fit_te(n, n, eps, V_F_AU, kTmin, kTmax, kTe, muc, muh)
    if (abs(kTe - kTmin) > 1d-14) call bad('energy at/below the bath value must clamp to the bath T')

    ! --- (5) monotone in the energy at fixed n ------------------------------------
    n = dirac_n_2d(0.02d0 / AU_EV, 1000d0 * KB_HA_K, V_F_AU)
    e1 = 2d0 * dirac_e_2d(dirac_mu_2d(n, 1000d0 * KB_HA_K, V_F_AU), 1000d0 * KB_HA_K, V_F_AU)
    e2 = 2d0 * dirac_e_2d(dirac_mu_2d(n, 3000d0 * KB_HA_K, V_F_AU), 3000d0 * KB_HA_K, V_F_AU)
    call dirac_fit_te(n, n, e1, V_F_AU, kTmin, kTmax, kTe1, muc, muh)
    call dirac_fit_te(n, n, e2, V_F_AU, kTmin, kTmax, kTe2, muc, muh)
    write(*,'(a,f8.1,a,f8.1,a)') '  fixed n: T_e(e1) = ', kTe1 / KB_HA_K, ' K (expect 1000), T_e(e2) = ', &
        kTe2 / KB_HA_K, ' K (expect 3000)'
    if (.not. (kTe2 > kTe1)) call bad('T_e not monotone in the energy')
    if (abs(kTe1 / KB_HA_K - 1000d0) > 10d0 .or. abs(kTe2 / KB_HA_K - 3000d0) > 30d0) &
        call bad('T_e inversion off by > 1 % at fixed n')

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (two-temperature fit: closed form, explicit-mesh moments, T_e/mu recovery, degenerate limit, fallbacks, monotone)'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if

contains

    subroutine mesh_case(T_K, mu_c, mu_h, label)
        real(8), intent(in) :: T_K, mu_c, mu_h
        character(*), intent(in) :: label
        integer, parameter :: NM = 700
        real(8) :: kT, kmax, dk, kx, ky, kk, fe, fh, n_m, p_m, e_m, n_a, p_a, e_a
        real(8) :: kTe, muc, muh
        integer :: i, j
        kT = T_K * KB_HA_K
        kmax = (max(mu_c, mu_h, 0d0) + 20d0 * kT) / V_F_AU
        dk = 2d0 * kmax / NM
        n_m = 0d0;  p_m = 0d0;  e_m = 0d0
        do i = 1, NM
            kx = -kmax + (i - 0.5d0) * dk
            do j = 1, NM
                ky = -kmax + (j - 0.5d0) * dk
                kk = sqrt(kx * kx + ky * ky)
                fe = 1d0 / (exp(min((V_F_AU * kk - mu_c) / kT, 60d0)) + 1d0)
                fh = 1d0 / (exp(min((V_F_AU * kk - mu_h) / kT, 60d0)) + 1d0)
                n_m = n_m + fe
                p_m = p_m + fh
                e_m = e_m + V_F_AU * kk * (fe + fh)
            end do
        end do
        n_m = n_m * 4d0 * dk * dk / (2d0 * PI)**2       ! g = 4 (spin x valley)
        p_m = p_m * 4d0 * dk * dk / (2d0 * PI)**2
        e_m = e_m * 4d0 * dk * dk / (2d0 * PI)**2
        n_a = dirac_n_2d(mu_c, kT, V_F_AU);  p_a = dirac_n_2d(mu_h, kT, V_F_AU)
        e_a = dirac_e_2d(mu_c, kT, V_F_AU) + dirac_e_2d(mu_h, kT, V_F_AU)
        if (abs(n_m - n_a) > 5d-3 * n_a) call bad('mesh electron density /= dirac_n_2d (' // label // ')')
        if (abs(p_m - p_a) > 5d-3 * p_a) call bad('mesh hole density /= dirac_n_2d (' // label // ')')
        if (abs(e_m - e_a) > 5d-3 * e_a) call bad('mesh energy density /= dirac_e_2d (' // label // ')')
        call dirac_fit_te(n_m, p_m, e_m, V_F_AU, 0.5d0 * kTmin, kTmax, kTe, muc, muh)
        write(*,'(a,a,a,f8.1,a,f8.1,a,f8.4,a,f8.4,a)') '  mesh moments (', label, '): T_e = ', kTe / KB_HA_K, &
            ' K (true ', T_K, '); mu_c = ', muc * AU_EV, ' eV, mu_h = ', muh * AU_EV, ' eV'
        if (abs(kTe - kT) > 1d-2 * kT) call bad('T_e recovery off by > 1 % (' // label // ')')
        if (abs(muc - mu_c) > 1d-2 * kT .or. abs(muh - mu_h) > 1d-2 * kT) &
            call bad('mu recovery off by > 1 % of kT (' // label // ')')
    end subroutine mesh_case

    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_dirac_te_fit
