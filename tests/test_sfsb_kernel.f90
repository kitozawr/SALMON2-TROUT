!
!  test_sfsb_kernel.f90  -  SFSB memory-integral ionization [B25 Eq. (3)]
!  (Boroumand et al., Rep. Prog. Phys. 88, 070501 (2025); wiki/10 sec. 6).
!
!  Validates the Volterra stepper sfsb_nc_series on a driven two-level system
!  against INDEPENDENT integrators and the letter's qualitative physics:
!    1) C = 0, weak resonant pulse: nc(t_end) matches an RK4 solution of the
!       exact two-level TDSE (all orders) to a few % -- pins the Om = 2 d E
!       convention, the Stark-shifted phase and the double-integral factors.
!    2) C = -tau/T2 (RTA): the exponential kernel is EXACTLY equivalent to the
!       Markovian ODE dI/dt = phi - I/T2; RK4 on that ODE must reproduce the
!       history sum (the memory machinery reduces correctly to Markov).
!    3) Dephasing ionization [B25 Fig. 1(a)]: for a BELOW-gap (3-photon) drive
!       the RTA bath enhances ionization by orders of magnitude over C = 0.
!    4) Dephasing-SUPPRESSED ionization [B25 Fig. 3]: low-T strong-coupling
!       Ohmic bath suppresses the same 3-photon ionization (Im C acts as a
!       dynamic bandgap ADDITION).
!    5) Killing Im C flips suppression toward enhancement [B25 Fig. 3(c)].
!    6) Window truncation: with a short-T2 RTA kernel, nwin = nt/2 equals the
!       full history (the truncation used in production is safe).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_sfsb_kernel
    use sbe_superres_ssbe, only: sfsb_nc_series, bath_corr_table
    implicit none
    integer, parameter :: NT = 8000
    real(8), parameter :: DT = 0.25d0            ! a.u.t; pulse length 2000 a.u.
    real(8), parameter :: TP = NT * DT
    real(8), parameter :: DE = 0.3d0             ! two-level gap [Ha]
    real(8), parameter :: DIP = 0.2d0            ! transition dipole [a.u.]
    complex(8) :: om(0:NT), ctab(0:NT), ctab2(0:NT)
    real(8) :: es(0:NT), nc(0:NT), nc0(0:NT), ncw(0:NT)
    real(8) :: w0, e0, t, ref, nc_tdse, nc_rta, nc_ohm, nc_imoff
    integer :: i, nfail

    nfail = 0

    ! =============== (1) C = 0 vs exact TDSE (resonant, weak) ===============
    w0 = DE
    e0 = 1.0d-4
    call fill_pulse(w0, e0)
    ctab(:) = (0d0, 0d0)
    call sfsb_nc_series(NT, DT, om, es, ctab, 0, nc)
    nc_tdse = tdse_nc(w0, e0)
    write(*,'(a,2es14.6)') '  resonant nc: sfsb vs TDSE(RK4) = ', nc(NT), nc_tdse
    if (nc(NT) <= 0d0) call bad('resonant nc not positive')
    if (abs(nc(NT) - nc_tdse) > 0.05d0 * nc_tdse) &
        call bad('C=0 sfsb disagrees with the exact TDSE at 2nd order')

    ! =============== (2) RTA kernel == Markovian ODE ========================
    call bath_corr_table(NT, DT, 1d-2, 0.1d0, 0.03d0, 'rta', ctab, t2_rta=150d0)
    call sfsb_nc_series(NT, DT, om, es, ctab, 0, nc)
    ref = rta_ode_nc(150d0)
    write(*,'(a,2es14.6)') '  RTA nc: history-sum vs Markov ODE(RK4) = ', nc(NT), ref
    if (abs(nc(NT) - ref) > 0.02d0 * max(ref, 1d-300)) &
        call bad('RTA memory kernel /= equivalent Markovian ODE')

    ! =============== (6) window truncation (kernel dead past T2) ============
    call sfsb_nc_series(NT, DT, om, es, ctab, NT / 2, ncw)
    if (abs(ncw(NT) - nc(NT)) > 1d-3 * max(nc(NT), 1d-300)) &
        call bad('nwin = nt/2 differs from full history for a short-T2 kernel')

    ! =============== (3) dephasing ionization (below-gap drive) =============
    w0 = DE / 3d0                                 ! 3-photon
    e0 = 0.05d0                                   ! Om_peak = 0.02 << DE
    call fill_pulse(w0, e0)
    ctab(:) = (0d0, 0d0)
    call sfsb_nc_series(NT, DT, om, es, ctab, 0, nc0)
    call bath_corr_table(NT, DT, 1d-2, 0.1d0, 0.03d0, 'rta', ctab, t2_rta=100d0)
    call sfsb_nc_series(NT, DT, om, es, ctab, 0, nc)
    nc_rta = nc(NT)
    write(*,'(a,2es14.6)') '  3-photon nc: no bath vs RTA = ', nc0(NT), nc_rta
    ! the virtual excursion during the pulse must RETURN without a bath
    write(*,'(a,2es14.6)') '  no-bath peak (virtual) vs final (real) = ', maxval(nc0), nc0(NT)
    if (nc0(NT) <= 0d0) call bad('no-bath 3-photon nc not positive (noise floor?)')
    if (nc_rta < 10d0 * nc0(NT)) &
        call bad('RTA must enhance below-gap ionization by >= 10x (dephasing ionization)')
    if (maxval(nc0) < 30d0 * nc0(NT)) &
        call bad('no-bath virtual excursion must exceed the surviving real part')

    ! =============== (4) low-T strong-coupling suppression ==================
    call bath_corr_table(NT, DT, 0d0, 1.5d0, 2.5d0 * w0, 'ohmic', ctab)
    call sfsb_nc_series(NT, DT, om, es, ctab, 0, nc)
    nc_ohm = nc(NT)
    write(*,'(a,2es14.6)') '  3-photon nc: no bath vs ohmic(T=0, jo=1.5) = ', nc0(NT), nc_ohm
    if (nc_ohm >= nc0(NT)) &
        call bad('low-T strong-coupling ohmic bath must SUPPRESS ionization')

    ! =============== (5) Im C -> 0 flips toward enhancement =================
    do i = 0, NT
        ctab2(i) = cmplx(real(ctab(i)), 0d0, 8)
    end do
    call sfsb_nc_series(NT, DT, om, es, ctab2, 0, nc)
    nc_imoff = nc(NT)
    write(*,'(a,2es14.6)') '  ohmic nc: Im C on vs off = ', nc_ohm, nc_imoff
    if (nc_imoff <= nc_ohm) &
        call bad('killing Im C must increase nc (the phase is the suppressor)')

    if (nfail == 0) then
        write(*,'(a)') 'PASS'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if

contains

    ! sin^2-envelope pulse E(t) = e0 sin^2(pi t/TP) cos(w(t-TP/2)) on [0, TP];
    ! fills om = 2 d E and the Stark-shifted gap es = sqrt(DE^2 + |om|^2).
    subroutine fill_pulse(w, e0in)
        real(8), intent(in) :: w, e0in
        integer :: ii
        real(8) :: tt, ef
        do ii = 0, NT
            tt = ii * DT
            ef = e0in * sin(3.14159265358979324d0 * tt / TP)**2 * cos(w * (tt - 0.5d0 * TP))
            om(ii) = cmplx(2d0 * DIP * ef, 0d0, 8)
            es(ii) = sqrt(DE**2 + abs(om(ii))**2)
        end do
    end subroutine fill_pulse

    ! Exact two-level TDSE (length gauge, H = [[0, V],[V, DE]], V = d E(t)),
    ! RK4 at dt/5 -- independent of the Volterra machinery.
    function tdse_nc(w, e0in) result(pc)
        real(8), intent(in) :: w, e0in
        real(8) :: pc
        complex(8) :: psi(2), k1(2), k2(2), k3(2), k4(2)
        real(8) :: h, tt
        integer :: ii, nfine
        nfine = 5 * NT
        h = DT / 5d0
        psi = (/ (1d0, 0d0), (0d0, 0d0) /)
        do ii = 0, nfine - 1
            tt = ii * h
            k1 = rhs(psi,                tt,            w, e0in)
            k2 = rhs(psi + 0.5d0*h*k1,   tt + 0.5d0*h,  w, e0in)
            k3 = rhs(psi + 0.5d0*h*k2,   tt + 0.5d0*h,  w, e0in)
            k4 = rhs(psi + h*k3,         tt + h,        w, e0in)
            psi = psi + (h / 6d0) * (k1 + 2d0*k2 + 2d0*k3 + k4)
        end do
        pc = abs(psi(2))**2
    end function tdse_nc

    function rhs(psi, tt, w, e0in) result(dp)
        complex(8), intent(in) :: psi(2)
        real(8), intent(in) :: tt, w, e0in
        complex(8) :: dp(2)
        real(8) :: v
        v = DIP * e0in * sin(3.14159265358979324d0 * min(max(tt, 0d0), TP) / TP)**2 &
            * cos(w * (tt - 0.5d0 * TP))
        dp(1) = (0d0, -1d0) * (v * psi(2))
        dp(2) = (0d0, -1d0) * (v * psi(1) + DE * psi(2))
    end function rhs

    ! Markovian equivalent of the exponential kernel: dI/dt = phi(t) - I/T2,
    ! nc' = Re[conj(om) e^{i th} I]/2, RK4 at dt/5 with linear interpolation
    ! of phi -- independent of the history sum.
    function rta_ode_nc(t2) result(pc)
        real(8), intent(in) :: t2
        real(8) :: pc
        real(8) :: th(0:NT)
        complex(8) :: phi(0:NT), Iacc, kk1, kk2, kk3, kk4
        real(8) :: h, gcur, gprev
        integer :: ii, s
        th(0) = 0d0
        do ii = 1, NT
            th(ii) = th(ii - 1) + 0.5d0 * (es(ii - 1) + es(ii)) * DT
        end do
        do ii = 0, NT
            phi(ii) = om(ii) * exp(cmplx(0d0, -th(ii), 8))
        end do
        h = DT / 5d0
        Iacc = (0d0, 0d0)
        pc = 0d0
        gprev = 0d0
        do ii = 0, NT - 1
            do s = 0, 4
                kk1 = phin(phi, ii + s/5d0)        - Iacc / t2
                kk2 = phin(phi, ii + (s+0.5d0)/5d0) - (Iacc + 0.5d0*h*kk1) / t2
                kk3 = phin(phi, ii + (s+0.5d0)/5d0) - (Iacc + 0.5d0*h*kk2) / t2
                kk4 = phin(phi, ii + (s+1d0)/5d0)   - (Iacc + h*kk3) / t2
                Iacc = Iacc + (h / 6d0) * (kk1 + 2d0*kk2 + 2d0*kk3 + kk4)
            end do
            gcur = 0.5d0 * real(conjg(om(ii + 1)) * exp(cmplx(0d0, th(ii + 1), 8)) * Iacc)
            pc = pc + 0.5d0 * (gprev + gcur) * DT
            gprev = gcur
        end do
    end function rta_ode_nc

    function phin(phi, x) result(p)
        complex(8), intent(in) :: phi(0:NT)
        real(8), intent(in) :: x
        complex(8) :: p
        integer :: i0
        real(8) :: f
        i0 = min(int(x), NT - 1)
        f = x - i0
        p = (1d0 - f) * phi(i0) + f * phi(i0 + 1)
    end function phin

    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_sfsb_kernel
