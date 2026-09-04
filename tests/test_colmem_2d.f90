!
!  test_colmem_2d.f90 - graphene 2D collisional-memory analog (wiki/10 sec. 8.11):
!  the Coulomb (Rana) sector's memory line = the 2D Dirac plasmon of the
!  instantaneous e-h plasma at the collision's screening momentum Q_TF, and the
!  plasmon-line population filter that feeds rana_auger_dpop.
!
!  Checks:
!    1) dirac_plasmon_2d: degenerate single-branch limit (T -> 0, holes
!       intrinsic) reproduces w_pl^2 = 2 E_F Q_TF/eps_r [Hwang-Das Sarma 2007]
!       with E_F = mu_c to 1e-6; the intrinsic plasma (mu = 0) gives the
!       thermal Drude weight 2 kT ln2 per branch; w_pl grows monotonically with
!       density and sits at 10..500 meV for n = 1e11..1e12 cm^-2 at 300 K.
!    2) plasmon-line filter (colmem_lines + colmem_pop_filter, 2 lines):
!       a CONSTANT density is a machine-exact fixed point (calibrated R07 rates
!       untouched); a 2*w_laser modulation (hw = 0.8 eV) is suppressed to the
!       line set's steady-state response |R(2w)| (within 15 %), i.e. the
!       virtual breathing filters out of the collision SOURCE.
!    3) rana_auger_dpop with the n2d_in/p2d_in overrides: equal to the internal
!       densities -> identical dpop (1e-14); a filtered (smaller) source ->
!       weaker net rate, trace still exactly 0, bounds respected.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_colmem_2d
    use sbe_superres_ssbe, only: dirac_plasmon_2d, rana_qtf, dirac_mu_2d, colmem_lines, &
                                 colmem_pop_filter, colmem_pop_init, bose_factor, rana_auger_dpop
    implicit none
    real(8), parameter :: PI      = 3.14159265358979323846d0
    real(8), parameter :: A0_CM   = 0.52917721067d-8
    real(8), parameter :: KB_HA_K = 3.166811563d-6
    real(8), parameter :: AU_EV   = 27.211386245988d0
    real(8), parameter :: V_F_AU  = 1.0d8 / 2.18769126364d8
    real(8), parameter :: EPS_R   = 10d0
    integer :: nfail
    real(8) :: kT, kT0, mu, wpl, wref, qtf, w11, w12, w0
    complex(8) :: cl(2), mul(2), z(2), resp
    integer :: nl, it, j
    real(8) :: tauc, tau, n0, n_t, ftil, wlas, amp, fmax, fmin, r_amp, tt

    nfail = 0
    kT = 300d0 * KB_HA_K

    ! --- (1) plasmon line ------------------------------------------------------
    kT0 = 1d-5                                        ! ~3 K: degenerate limit
    mu  = 0.01d0                                      ! E_F = 0.27 eV electrons; holes intrinsic
    wpl = dirac_plasmon_2d(mu, 0d0, kT0, V_F_AU, EPS_R)
    qtf = rana_qtf(mu, 0d0, kT0, V_F_AU, EPS_R)
    wref = sqrt(2d0 * (mu + 2d0 * kT0 * log(2d0)) * qtf / EPS_R)
    if (abs(wpl - wref) > 1d-6 * wref) call bad('degenerate limit: w_pl^2 /= 2 E_F Q_TF/eps_r')
    w0 = dirac_plasmon_2d(0d0, 0d0, kT, V_F_AU, EPS_R)
    qtf = rana_qtf(0d0, 0d0, kT, V_F_AU, EPS_R)
    wref = sqrt(2d0 * (4d0 * kT * log(2d0)) * qtf / EPS_R)
    if (abs(w0 - wref) > 1d-10 * wref) call bad('intrinsic plasma: Drude weight /= 2 kT ln2 per branch')
    mu  = dirac_mu_2d(1d11 * A0_CM**2, kT, V_F_AU)
    w11 = dirac_plasmon_2d(mu, -mu, kT, V_F_AU, EPS_R)
    mu  = dirac_mu_2d(1d12 * A0_CM**2, kT, V_F_AU)
    w12 = dirac_plasmon_2d(mu, -mu, kT, V_F_AU, EPS_R)
    write(*,'(a,f7.2,a,f7.2,a,f7.2,a)') '  w_pl(300 K): intrinsic ', w0 * AU_EV * 1d3, &
        ' meV,  n=p=1e11: ', w11 * AU_EV * 1d3, ' meV,  n=p=1e12: ', w12 * AU_EV * 1d3, ' meV'
    if (.not. (w0 < w11 .and. w11 < w12)) call bad('w_pl not monotone in density')
    if (w11 * AU_EV * 1d3 < 10d0 .or. w12 * AU_EV * 1d3 > 500d0) &
        call bad('w_pl outside the expected 10..500 meV window')

    ! --- (2) plasmon-line filter -----------------------------------------------
    tauc = 1d0 / (0.1d0 / AU_EV)                      ! tau_c = hbar/sigma_E, sigma_E = 0.1 eV
    call colmem_lines(1, (/ w12 /), (/ 1d0 /), (/ bose_factor(w12, kT) /), tauc, nl, cl, mul)
    if (nl /= 2) call bad('plasmon line set should have 2 lines (emission + absorption)')
    tau = 4d0
    n0  = 1d12 * A0_CM**2
    call colmem_pop_init(nl, mul(1:nl), tau, n0, z(1:nl))
    do it = 1, 500
        call colmem_pop_filter(nl, cl(1:nl), mul(1:nl), tau, n0, z(1:nl), ftil)
    end do
    if (abs(ftil - n0) > 1d-12 * n0) call bad('constant density is not a fixed point of the plasmon filter')
    wlas = 0.8d0 / AU_EV
    resp = (0d0, 0d0)
    do j = 1, nl
        resp = resp + cl(j) / (mul(j) + cmplx(0d0, 2d0 * wlas, 8))
    end do
    r_amp = abs(resp)
    call colmem_pop_init(nl, mul(1:nl), tau, n0, z(1:nl))
    fmax = -huge(1d0);  fmin = huge(1d0)
    do it = 1, 4000
        tt  = it * tau
        n_t = n0 * (1d0 + 0.5d0 * cos(2d0 * wlas * tt))
        call colmem_pop_filter(nl, cl(1:nl), mul(1:nl), tau, n_t, z(1:nl), ftil)
        if (it > 3000) then
            fmax = max(fmax, ftil);  fmin = min(fmin, ftil)
        end if
    end do
    amp = 0.5d0 * (fmax - fmin) / (0.5d0 * n0)
    write(*,'(a,f7.4,a,f7.4,a)') '  2*w_laser breathing transmitted: ', amp, '  (line-set |R(2w)| = ', r_amp, ')'
    if (abs(amp - r_amp) > 0.15d0 * r_amp) call bad('filtered modulation amplitude /= |R(2w)|')
    if (amp > 0.2d0) call bad('the 2*w_laser breathing is not suppressed (amp > 0.2)')
    if (abs(0.5d0 * (fmax + fmin) - n0) > 0.05d0 * n0) call bad('the slow (mean) density does not pass (R(0) = 1)')

    ! --- (3) rana_auger_dpop source overrides ----------------------------------
    call check_overrides()

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (2D colmem analog: Dirac-plasmon line limits, plasmon-line source filter, Rana overrides)'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if

contains

    subroutine check_overrides()
        integer, parameter :: NK = 64, NBA = 2, IV = 1, IC = 2
        real(8), parameter :: OCC = 2d0, F0 = 0.02d0
        real(8) :: f(NBA, NK), eval(NBA, NK), dpop1(NBA, NK), dpop2(NBA, NK), dpop3(NBA, NK)
        real(8) :: area, r1, r2, r3, n2d, p2d
        integer :: a
        do a = 1, NK
            eval(IV, a) = -0.10d0 * dble(a) / dble(NK)
            eval(IC, a) = +0.30d0 * dble(a) / dble(NK)
        end do
        f(IC, :) = F0
        f(IV, :) = OCC - F0
        area = F0 / (1d12 * A0_CM**2)                       ! n = p = 1e12 cm^-2 (recombination regime)
        n2d = sum(f(IC:NBA, :)) / (dble(NK) * area)
        p2d = (dble(IV) * OCC - sum(f(1:IV, :)) / dble(NK)) / area
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, EPS_R, 41.34d0, dpop1, r1)
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, EPS_R, 41.34d0, dpop2, r2, &
                             n2d_in=n2d, p2d_in=p2d)
        if (maxval(abs(dpop1 - dpop2)) > 1d-14 .or. abs(r1 - r2) > 1d-14 * abs(r1)) &
            call bad('override equal to the internal densities must give identical dpop')
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, EPS_R, 41.34d0, dpop3, r3, &
                             n2d_in=0.3d0 * n2d, p2d_in=0.3d0 * p2d)
        write(*,'(a,es11.3,a,es11.3)') '  net rate: raw source ', r1, '   filtered (x0.3) source ', r3
        if (.not. (r1 > 0d0 .and. r3 > 0d0 .and. r3 < r1)) &
            call bad('filtered (smaller) source must give a weaker net recombination rate')
        if (abs(sum(dpop3)) > 1d-13) call bad('trace not conserved with the source override')
        if (any(f + dpop3 < -1d-14) .or. any(f + dpop3 > OCC + 1d-14)) call bad('bounds violated with the override')
    end subroutine check_overrides

    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_colmem_2d
