!
!  test_rana_saturation.f90 - graphene 2D Rana Auger/CM: population SATURATION
!  at the generation/recombination balance (exercise x14; wiki/07 sec.6).
!
!  Physics pinned here: on the gapless Dirac cone impact ionization (carrier
!  multiplication, CVCC) and Auger recombination (CCCV) are exact time-reverses
!  [R07 Eqs. 14/17]. Their balance R = G holds iff the two quasi-Fermi levels
!  coincide, mu_c = mu_v; for a symmetric pair population n = p that is mu = 0,
!  i.e. the INTRINSIC thermal density of the cone
!
!      n_i(T) = (pi/6) (k_B T / hbar v_F)^2        (g = 4: spin x valley)
!             = 8.1e10 cm^-2 at 300 K.
!
!  Above n_i the pair population net-recombines, below it net-generates:
!  n_i(T) is the density at which the population SATURATES ("Auger and
!  ionization balance each other above some density"). NOTE: T is the
!  temperature the rates are evaluated at -- in the solver that is the e-ph
!  BATH temperature (rana_kt_au), not a carrier temperature T_e (not extracted
!  yet; see the x14 README) -- so the solver's plateau is n_i(T_bath).
!
!  Checks:
!    1) dirac_n_2d(mu=0, kT, v) == (pi/6)(kT/v)^2 to 1e-3 (the FD integral).
!    2) sign structure of R - G at n = p: < 0 (generation) at 0.5 n_i, > 0
!       (recombination) at 2 n_i; the bisection root lies within 1 % of n_i
!       at 300 K and 1000 K.
!    3) n_0(1000 K)/n_0(300 K) = (1000/300)^2 to 3 % (the T^2 law).
!    4) TWO-SIDED SATURATION through the CPTP channel rana_auger_dpop:
!       iterating from n = 0.1 n_i (carrier-multiplication growth) and from
!       n = 10 n_i (Auger decay) both converge MONOTONICALLY to n_i within 2 %.
!    5) trace conserved to machine precision at every iteration.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_rana_saturation
    use sbe_superres_ssbe, only: dirac_n_2d, dirac_mu_2d, rana_rcccv, rana_auger_dpop
    implicit none
    real(8), parameter :: PI      = 3.14159265358979323846d0
    real(8), parameter :: A0_CM   = 0.52917721067d-8         ! Bohr in cm
    real(8), parameter :: KB_HA_K = 3.166811563d-6           ! k_B [Ha/K]
    real(8), parameter :: V_F_AU  = 1.0d8 / 2.18769126364d8  ! 1e8 cm/s (the R07 constant)
    real(8), parameter :: EPS_R   = 10d0                     ! R07 benchmark substrate
    integer :: nfail
    real(8) :: kT300, kT1000, ni300, ni1000, n0_300, n0_1000, nnum, rn

    nfail = 0
    kT300  = 300d0  * KB_HA_K
    kT1000 = 1000d0 * KB_HA_K
    ni300  = PI / 6d0 * (kT300  / V_F_AU)**2
    ni1000 = PI / 6d0 * (kT1000 / V_F_AU)**2
    write(*,'(a)') '  balance (saturation) density n_i(T) = (pi/6)(k_B T/hbar v_F)^2:'
    write(*,'(a,es11.3,a,es11.3,a,es11.3,a)') '    300 K: ', ni300 / A0_CM**2, &
        '   1000 K: ', ni1000 / A0_CM**2, '   3000 K: ', 9d0 * ni1000 / A0_CM**2, '  cm^-2'

    ! --- (1) the FD integral reproduces the closed form --------------------
    nnum = dirac_n_2d(0d0, kT300, V_F_AU)
    if (abs(nnum - ni300) > 1d-3 * ni300) call bad('dirac_n_2d(mu=0) /= (pi/6)(kT/v)^2')

    ! --- (2) sign structure + root of R - G at n = p ------------------------
    rn = rnet_sym(0.5d0 * ni300, kT300)
    if (rn >= 0d0) call bad('net rate not negative (generation) at 0.5 n_i')
    rn = rnet_sym(2d0 * ni300, kT300)
    if (rn <= 0d0) call bad('net rate not positive (recombination) at 2 n_i')
    n0_300  = root_sym(kT300,  ni300)
    n0_1000 = root_sym(kT1000, ni1000)
    write(*,'(a,es11.3,a,f7.4)') '  R=G root  300 K: n_0 = ', n0_300  / A0_CM**2, ' cm^-2, n_0/n_i = ', n0_300  / ni300
    write(*,'(a,es11.3,a,f7.4)') '  R=G root 1000 K: n_0 = ', n0_1000 / A0_CM**2, ' cm^-2, n_0/n_i = ', n0_1000 / ni1000
    if (abs(n0_300  - ni300)  > 1d-2 * ni300)  call bad('R=G root off n_i(300 K) by > 1 %')
    if (abs(n0_1000 - ni1000) > 1d-2 * ni1000) call bad('R=G root off n_i(1000 K) by > 1 %')

    ! --- (3) T^2 law ----------------------------------------------------------
    if (abs(n0_1000 / n0_300 - (1000d0 / 300d0)**2) > 3d-2 * (1000d0 / 300d0)**2) &
        call bad('n_0(T) does not scale as T^2')

    ! --- (4)+(5) two-sided saturation through the CPTP channel ----------------
    call relax_to_fixed_point(0.1d0 * ni300, kT300, ni300, 'from below (carrier multiplication)')
    call relax_to_fixed_point(10d0  * ni300, kT300, ni300, 'from above (Auger recombination)  ')

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (graphene Rana balance: n_0 = n_i(T) = (pi/6)(k_B T/hbar v_F)^2, '// &
                       'T^2 law, two-sided monotone CPTP saturation)'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if

contains

    ! net rate R - G [a.u.^-2 a.u.t^-1] for a symmetric pair population n = p
    function rnet_sym(n, kT) result(r)
        real(8), intent(in) :: n, kT
        real(8) :: r, mu_c, mu_v
        mu_c =  dirac_mu_2d(n, kT, V_F_AU)
        mu_v = -mu_c
        r = rana_rcccv( mu_c,  mu_v, kT, V_F_AU, EPS_R, .false.) &
          + rana_rcccv(-mu_v, -mu_c, kT, V_F_AU, EPS_R, .false.) &
          - rana_rcccv( mu_c,  mu_v, kT, V_F_AU, EPS_R, .true.)  &
          - rana_rcccv(-mu_v, -mu_c, kT, V_F_AU, EPS_R, .true.)
    end function rnet_sym

    ! bisection root of R - G on [0.2 n_i, 5 n_i]
    function root_sym(kT, ni) result(n0)
        real(8), intent(in) :: kT, ni
        real(8) :: n0, lo, hi, mid
        integer :: it
        lo = 0.2d0 * ni
        hi = 5d0 * ni
        do it = 1, 40
            mid = 0.5d0 * (lo + hi)
            if (rnet_sym(mid, kT) < 0d0) then
                lo = mid
            else
                hi = mid
            end if
        end do
        n0 = 0.5d0 * (lo + hi)
    end function root_sym

    ! iterate the CPTP channel on a two-band cone from n_start until stationary
    subroutine relax_to_fixed_point(n_start, kT, ni, label)
        real(8), intent(in) :: n_start, kT, ni
        character(*), intent(in) :: label
        integer, parameter :: NK = 48, NBA = 2, IV = 1, IC = 2
        real(8), parameter :: OCC = 2d0, F0 = 0.05d0
        real(8) :: f(NBA, NK), dpop(NBA, NK), eval(NBA, NK)
        real(8) :: area, n, nprev, tau, rnet, tr0, tr
        integer :: a, it
        logical :: monotone
        ! Dirac-like 2-band spectrum (VB below 0, CB above): eval only feeds the
        ! trace-neutral A6 energy shuffle; the density balance is spectrum-free.
        do a = 1, NK
            eval(IV, a) = -0.20d0 * dble(a) / dble(NK)
            eval(IC, a) = +0.40d0 * dble(a) / dble(NK)
        end do
        ! uniform per-k populations: n2d = sum_k f_cb/(nk*area) = F0/area
        f(IC, :) = F0
        f(IV, :) = OCC - F0
        area = F0 / n_start
        tau  = 2000d0                          ! ~48 fs per channel application
        tr0  = sum(f)
        nprev = n_of(f, area)
        monotone = .true.
        do it = 1, 6000
            call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, EPS_R, &
                                 tau, dpop, rnet)
            f = f + dpop
            tr = sum(f)
            if (abs(tr - tr0) > 1d-10 * tr0) then
                call bad('trace not conserved during relaxation ' // label)
                exit
            end if
            n = n_of(f, area)
            if (n_start < ni) then
                if (n < nprev - 1d-12 * ni) monotone = .false.
            else
                if (n > nprev + 1d-12 * ni) monotone = .false.
            end if
            if (abs(n - nprev) < 1d-6 * ni .and. it > 10) exit
            nprev = n
        end do
        write(*,'(a,a,a,es11.3,a,f7.4,a,i0,a)') '  relaxation ', label, ': n_end = ', &
            n / A0_CM**2, ' cm^-2 (n_end/n_i = ', n / ni, ') after ', it, ' channel steps'
        if (.not. monotone) call bad('non-monotone approach to the balance density ' // label)
        if (abs(n - ni) > 2d-2 * ni) call bad('population did not saturate at n_i within 2 % ' // label)
    end subroutine relax_to_fixed_point

    function n_of(f, area) result(n)
        real(8), intent(in) :: f(:, :), area
        real(8) :: n
        n = sum(f(2, :)) / (dble(size(f, 2)) * area)
    end function n_of

    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_rana_saturation
