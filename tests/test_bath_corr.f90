!
!  test_bath_corr.f90  -  SFSB bath correlation function C(t) [B25 Eq. (5)]
!  (Boroumand et al., Rep. Prog. Phys. 88, 070501 (2025); wiki/10 sec. 6).
!
!  Validates the numeric C(t) table against every closed form the letter and
!  the W(w) = jo g(|w|)/w normalization admit:
!    1) C(0) = 0 exactly; Re C <= 0 for all tau (decay, never gain).
!    2) OHMIC Im C(t) = 2 jo atan(wc t)          (exact, any T).
!    3) OHMIC Re C(t; T=0) = -jo ln(1+wc^2 t^2)  (exact).
!    4) DEBYE Im C(t) = pi jo (1-exp(-wc t))     (exact, any T).
!    5) HIGH-T ANCHOR [B25 sec. 2, printed]: d(Re C)/dt -> -2 pi kB T jo / hbar,
!       i.e. exp[C] -> exp[-t/T2] with T2 = hbar/(2 pi kB T jo) = bath_t2_high_t
!       -- for BOTH ohmic and debye profiles (the anchor pins the W(w) ~ jo/w
!       normalization independently of the cutoff shape).
!    6) RTA model: C(t) = -t/T2 exactly (given T2, and the derived default).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_bath_corr
    use sbe_superres_ssbe, only: bath_corr_table, bath_t2_high_t
    implicit none
    real(8), parameter :: PI      = 3.14159265358979323846d0
    real(8), parameter :: KB_HA_K = 3.166811563d-6   ! k_B [Ha/K]
    integer, parameter :: NT = 2000
    real(8), parameter :: DTAU = 0.5d0               ! a.u.t; tau up to 1000
    real(8), parameter :: JO = 0.13d0, WC = 0.03d0   ! dimensionless / Ha
    complex(8) :: ctab(0:NT)
    real(8) :: kT, tau, ref, slope, t2
    integer :: m, nfail

    nfail = 0

    ! --- OHMIC at T = 300 K: C(0) = 0, Re <= 0, exact Im ------------------
    kT = 300d0 * KB_HA_K
    call bath_corr_table(NT, DTAU, kT, JO, WC, 'ohmic', ctab)
    if (abs(ctab(0)) > 1d-12) call bad('ohmic C(0) /= 0')
    do m = 0, NT
        if (real(ctab(m)) > 1d-10) call bad('ohmic Re C > 0 (bath gain)')
    end do
    do m = 200, NT, 600
        tau = m * DTAU
        ref = 2d0 * JO * atan(WC * tau)
        if (abs(aimag(ctab(m)) - ref) > 1d-4 * max(abs(ref), 1d-3)) then
            write(*,'(a,f9.1,2es14.6)') '  ohmic Im mismatch at tau=', tau, aimag(ctab(m)), ref
            call bad('ohmic Im C /= 2 jo atan(wc t)')
        end if
    end do

    ! --- OHMIC at T = 0: exact Re ------------------------------------------
    call bath_corr_table(NT, DTAU, 0d0, JO, WC, 'ohmic', ctab)
    do m = 200, NT, 600
        tau = m * DTAU
        ref = -JO * log(1d0 + (WC * tau)**2)
        if (abs(real(ctab(m)) - ref) > 1d-4 * max(abs(ref), 1d-3)) then
            write(*,'(a,f9.1,2es14.6)') '  ohmic T=0 Re mismatch at tau=', tau, real(ctab(m)), ref
            call bad('ohmic Re C(T=0) /= -jo ln(1+wc^2 t^2)')
        end if
    end do

    ! --- DEBYE at T = 300 K: exact Im ---------------------------------------
    kT = 300d0 * KB_HA_K
    call bath_corr_table(NT, DTAU, kT, JO, WC, 'debye', ctab)
    if (abs(ctab(0)) > 1d-12) call bad('debye C(0) /= 0')
    do m = 200, NT, 600
        tau = m * DTAU
        ref = PI * JO * (1d0 - exp(-WC * tau))
        if (abs(aimag(ctab(m)) - ref) > 1d-3 * max(abs(ref), 1d-3)) then
            write(*,'(a,f9.1,2es14.6)') '  debye Im mismatch at tau=', tau, aimag(ctab(m)), ref
            call bad('debye Im C /= pi jo (1-exp(-wc t))')
        end if
    end do

    ! --- HIGH-T anchor: slope of Re C -> -2 pi kT jo (both profiles) --------
    ! At kT = 3e4 K (kT >> wc) and wc*tau in [18, 30] the exact ohmic high-T
    ! slope is -4 kT jo atan(wc tau) -> -2 pi kT jo with O(1/(wc tau)) + O(wc/kT)
    ! corrections; 5% tolerance covers both.
    kT = 3d4 * KB_HA_K
    call bath_corr_table(NT, DTAU, kT, JO, WC, 'ohmic', ctab)
    slope = (real(ctab(2000)) - real(ctab(1200))) / (800d0 * DTAU)
    ref = -2d0 * PI * kT * JO
    write(*,'(a,2es14.6)') '  ohmic high-T slope vs -2 pi kT jo: ', slope, ref
    if (abs(slope - ref) > 0.05d0 * abs(ref)) call bad('ohmic high-T slope anchor')
    t2 = bath_t2_high_t(kT, JO)
    if (abs(t2 - 1d0 / (2d0 * PI * kT * JO)) > 1d-12) call bad('bath_t2_high_t formula')
    call bath_corr_table(NT, DTAU, kT, JO, WC, 'debye', ctab)
    slope = (real(ctab(2000)) - real(ctab(1200))) / (800d0 * DTAU)
    write(*,'(a,2es14.6)') '  debye high-T slope vs -2 pi kT jo: ', slope, ref
    if (abs(slope - ref) > 0.05d0 * abs(ref)) call bad('debye high-T slope anchor')

    ! --- RTA: exact linear decay --------------------------------------------
    call bath_corr_table(NT, DTAU, kT, JO, WC, 'rta', ctab, t2_rta=120d0)
    if (abs(real(ctab(1000)) + 1000d0 * DTAU / 120d0) > 1d-12) call bad('rta C /= -t/T2 (given T2)')
    if (abs(aimag(ctab(1000))) > 0d0) call bad('rta Im C /= 0')
    call bath_corr_table(NT, DTAU, kT, JO, WC, 'rta', ctab)
    if (abs(real(ctab(1000)) + 1000d0 * DTAU / bath_t2_high_t(kT, JO)) > 1d-10) &
        call bad('rta C /= -t/T2 (derived T2)')

    ! --- jo <= 0: bath off ---------------------------------------------------
    call bath_corr_table(NT, DTAU, kT, 0d0, WC, 'ohmic', ctab)
    if (maxval(abs(ctab)) > 0d0) call bad('jo = 0 must give C = 0')

    if (nfail == 0) then
        write(*,'(a)') 'PASS'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if

contains
    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_bath_corr
