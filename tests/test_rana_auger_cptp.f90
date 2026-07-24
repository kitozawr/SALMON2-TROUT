!
!  test_rana_auger_cptp.f90  -  graphene 2D Rana Auger/CM as a CPTP channel
!  (rana_auger_dpop: the wiki/00 TODO-1 wiring of the validated [R07]
!  primitives into a population channel).
!
!  Checks:
!    1) TRACE: sum(dpop) = 0 to machine precision (number conservation).
!    2) BOUNDS: 0 <= f + dpop <= occ_max everywhere (CPTP safety).
!    3) EQUILIBRIUM FIXED POINT: thermal populations (same mu=0 for both
!       branches) => R = G (detailed balance) => dpop == 0 identically.
!    4) NET RECOMBINATION: an excess e-h pair population relaxes CB -> VB
!       (sum over CB of dpop < 0, sum over VB > 0), and the implied lifetime
!       matches the R07 scale (~ps at 1e12 cm^-2, 300 K, eps 10).
!    5) SATURATION: tau -> huge stays within min(avail, room) (no overshoot).
!    6) EMPTY: no carriers -> no-op.
!
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_rana_auger_cptp
    use sbe_superres_ssbe, only: rana_auger_dpop, dirac_mu_2d
    implicit none
    real(8), parameter :: A0_CM   = 0.52917721067d-8      ! Bohr in cm
    real(8), parameter :: KB_HA_K = 3.166811563d-6        ! k_B [Ha/K]
    real(8), parameter :: V_F_AU  = 1.0d8 / 2.18769126364d8
    real(8), parameter :: AU_T_FS = 0.02418884326505d0
    integer, parameter :: NK = 64, NBA = 2, IV = 1, IC = 2
    real(8), parameter :: OCC = 2d0                       ! scalar occupation

    real(8) :: f(NBA, NK), dpop(NBA, NK), rnet, eval(NBA, NK), ebal
    real(8) :: kT, area, n_target, ftot_cb, dsum, tau, tau_ps
    real(8) :: avail, room
    integer :: ik, nfail, a

    nfail = 0
    kT   = 300d0 * KB_HA_K
    ! Dirac-like two-band spectrum: VB mirrored below zero, CB spread above --
    ! gives the A6 energy shuffle genuine phase space above the CB mean.
    do a = 1, NK
        eval(IV, a) = -0.10d0 * dble(a) / dble(NK)
        eval(IC, a) = +0.30d0 * dble(a) / dble(NK)   ! wide CB: room for the A6 shuffle
    end do
    ! Pick the 2D cell area so that a convenient per-k population gives
    ! n = 1e12 cm^-2: n_au = ftot_cb/(nk*area).
    n_target = 1.0d12 * A0_CM**2                          ! a.u.^-2
    ftot_cb  = 0.02d0 * NK                                ! 0.02 e-/k in the CB
    area     = ftot_cb / (dble(NK) * n_target)
    tau      = 41.34d0                                    ! ~1 fs in a.u.t

    ! ---- (6) EMPTY: full valence, empty conduction -> no-op --------------
    f(IV, :) = OCC
    f(IC, :) = 0d0
    call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                         tau, dpop, rnet)
    if (maxval(abs(dpop)) > 1d-14) call bad('empty CB is not a no-op')

    ! ---- (4) NET RECOMBINATION of an excess e-h pair population ----------
    ! n = p = 1e12 cm^-2 spread uniformly (equal hole count in the VB).
    ! populate the LOW quarter of the CB (a cooled distribution): the A6
    ! third-carrier target Ec_bar + E_rel then lies INSIDE the band window.
    f(IC, :) = 0d0
    do ik = 1, NK/4
        f(IC, ik) = ftot_cb / dble(NK/4)
    end do
    f(IV, :) = OCC - ftot_cb / dble(NK)
    call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                         tau, dpop, rnet)
    ! (1) trace conservation
    dsum = sum(dpop)
    if (abs(dsum) > 1d-12 * max(sum(abs(dpop)), 1d-300)) &
        call bad('trace not conserved (sum dpop /= 0)')
    ! (2) bounds
    if (minval(f + dpop) < -1d-14) call bad('negative population produced')
    if (maxval(f + dpop) > OCC + 1d-14) call bad('population above occ_max')
    ! direction: excess pairs must NET-recombine (CB loses, VB gains)
    if (rnet <= 0d0) call bad('excess pairs did not give net recombination (rnet<=0)')
    if (sum(dpop(IC, :)) >= 0d0) call bad('CB did not lose population')
    if (sum(dpop(IV, :)) <= 0d0) call bad('VB did not gain population')
    ! R07 lifetime scale: tau_r = n/(R-G) ~ 0.5..3 ps at 1e12 / 300 K / eps10
    ! (the net rate at n = p = 1e12 is close to the gross R: G is tiny there).
    tau_ps = (ftot_cb / (dble(NK) * area)) / rnet * AU_T_FS * 1d-3
    write(*,'(a,f8.3,a)') '  implied tau_r(1e12 cm^-2, 300 K) = ', tau_ps, &
        ' ps  [R07 window 0.5..3 ps]'
    if (tau_ps < 0.5d0 .or. tau_ps > 3.0d0) &
        call bad('implied lifetime outside the R07 benchmark window')
    ! A6: energy bookkeeping -- the third-carrier shuffle must cancel most of
    ! the electronic energy the naive pair transfer would have destroyed.
    ! Without the shuffle, sum(E*dpop) = -(Ec_bar - Ev_bar)*dn < 0 strictly.
    block
        real(8) :: erel_t, dn_t
        dn_t = -sum(dpop(IV, :))                       ! pairs moved (VB gain)
        erel_t = sum(eval(IC,:)*f(IC,:))/max(sum(f(IC,:)),1d-30) &
               - sum(eval(IV,:)*(OCC-f(IV,:)))/max(sum(OCC-f(IV,:)),1d-30)
        ebal = sum(eval * dpop)
        ! naive pair transfer alone would give ebal ~ -dn*erel_t; the shuffle
        ! must cancel most of it (Gaussian-width residual allowed)
        if (abs(ebal) > 0.35d0 * abs(dn_t * erel_t)) &
            call bad('A6: energy not balanced by the third-carrier shuffle')
    end block
    ! trace must STILL be exact with the shuffle on
    if (abs(sum(dpop)) > 1d-12 * max(sum(abs(dpop)), 1d-300)) &
        call bad('A6: shuffle broke trace conservation')

    ! ---- (3) EQUILIBRIUM FIXED POINT: mu = 0 thermal populations ----------
    ! Same THERMAL pair density on both branches (electrons in CB = holes in
    ! VB at mu=0): the kernel inverts both to mu ~ 0 => R = G => dpop = 0.
    block
        real(8) :: n_th, fth
        n_th = thermal_density(kT, V_F_AU)     ! mu = 0 sheet density [a.u.^-2]
        fth  = n_th * dble(NK) * area / dble(NK)
        do ik = 1, NK
            f(IC, ik) = fth
            f(IV, ik) = OCC - fth
        end do
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                             tau, dpop, rnet)
        ! dpop must be tiny relative to what the same tau does to the excess
        ! case (rnet is the residual of the R = G cancellation; the mu
        ! inversion is bisection-exact to ~1e-12, so the residual is numerical)
        if (sum(abs(dpop)) > 1d-6 * ftot_cb) &
            call bad('equilibrium (mu=0) populations are not a fixed point')
    end block

    ! ---- (5) SATURATION: tau -> huge never overshoots ---------------------
    do ik = 1, NK
        f(IC, ik) = ftot_cb / dble(NK)
        f(IV, ik) = OCC - ftot_cb / dble(NK)
    end do
    call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                         1d30, dpop, rnet)
    avail = sum(f(IC, :))
    room  = dble(IV) * OCC * dble(NK) - sum(f(1:IV, :))
    if (-sum(dpop(IC, :)) > min(avail, room) + 1d-10) &
        call bad('tau->inf overshoots min(avail, room)')
    if (minval(f + dpop) < -1d-12) call bad('tau->inf negative population')

    ! ---- (7) cost-preserving energy-mask (interface uniformity) -----------
    ! Rana is mean-field (no cost win), so the window is threaded only for a
    ! uniform ring interface. It must: (a) reproduce the default when wide,
    ! (b) be a no-op when it excludes every band, (c) stay trace-exact under a
    ! partial window (drain/fill only in-window states).
    block
        real(8) :: dref(NBA, NK), dwide(NBA, NK), dwin(NBA, NK), rn
        f(IC, :) = 0d0
        do ik = 1, NK/4
            f(IC, ik) = ftot_cb / dble(NK/4)
        end do
        f(IV, :) = OCC - ftot_cb / dble(NK)
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                             tau, dref, rn)
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                             tau, dwide, rn, e_src_lo=-1d30, e_src_hi=1d30)
        if (maxval(abs(dwide - dref)) > 1d-12) call bad('mask wide window /= default')
        ! window above every band (eval spans [-0.1, 0.3]) -> exact no-op
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                             tau, dwin, rn, e_src_lo=100d0, e_src_hi=200d0)
        if (maxval(abs(dwin)) > 1d-14) call bad('mask excluding all states is not a no-op')
        ! partial window (keep CB, drop the deepest VB) must stay trace-exact
        call rana_auger_dpop(NK, NBA, eval, f, OCC, IV, IC, area, kT, V_F_AU, 10d0, &
                             tau, dwin, rn, e_src_lo=-0.06d0, e_src_hi=1d30)
        if (abs(sum(dwin)) > 1d-12 * max(sum(abs(dwin)), 1d-300)) &
            call bad('mask: partial window broke trace conservation')
    end block

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (2D Rana Auger CPTP channel: trace, bounds, equilibrium '// &
                       'fixed point, R07 lifetime scale, saturation, empty no-op)'
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

    ! mu = 0 thermal sheet density of one Dirac branch (same integral as the
    ! module's dirac_n_2d, reproduced here for the fixed-point test input).
    function thermal_density(kTl, v) result(n)
        real(8), intent(in) :: kTl, v
        real(8) :: n, k, dk, kmax
        integer :: i
        integer, parameter :: NG = 400
        kmax = 20d0 * kTl / v
        dk = kmax / NG
        n = 0d0
        do i = 1, NG
            k = (i - 0.5d0) * dk
            n = n + k / (exp(min(v * k / kTl, 60d0)) + 1d0)
        end do
        n = n * dk * 4d0 / (2d0 * acos(-1d0))
    end function thermal_density
end program test_rana_auger_cptp
