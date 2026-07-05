!
!  test_ii_interk_cptp.f90  -  momentum-conserving inter-k impact ionization.
!
!  Tests ii_interk_dpop: the 2-particle event hot-e(k1,c) + valence-e(k2,v) ->
!  e(k1',c) + e(k2',c) with crystal momentum k1+k2=k1'+k2' (mod G, via the MP
!  index map) and a broadened Fermi-golden-rule rate. The net dpop must:
!    1) conserve the total trace EXACTLY (sum(dpop)=0) -- CPTP,
!    2) MULTIPLY carriers: conduction population grows, valence drops (a hole),
!    3) keep populations in [0,occ_max],
!    4) be a no-op below threshold.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_ii_interk_cptp
    use sbe_superres_ssbe, only: ii_interk_dpop, build_vq_table, t_ring_opts
    implicit none
    integer, parameter :: nk = 4, nba = 2, iv = 1, ic = 2
    integer :: kn(3), kidx(3, nk), klut(0:nk-1), nfail, a, k
    real(8) :: eval(nba, nk), f(nba, nk), dpop(nba, nk)
    real(8) :: occ_max, a2half, ecbm, eth, pref, expo, sigma, tau
    ! CDRB-screened Cartesian umklapp weight parameters (cubic test cell,
    ! Si-like eps_inf; q2reg regularises the q=0 term of the discrete sum)
    real(8) :: bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
    real(8) :: dcb, dvb, total
    real(8), parameter :: TOL = 1d-11

    nfail = 0
    kn = (/ 2, 2, 1 /)
    ! triples (m1,m2,0) and the flattened lookup lidx = m1 + 2*(m2 + 2*m3)
    kidx(:,1) = (/0,0,0/); kidx(:,2) = (/1,0,0/)
    kidx(:,3) = (/0,1,0/); kidx(:,4) = (/1,1,0/)
    do k = 1, nk
        klut(kidx(1,k) + kn(1)*(kidx(2,k) + kn(2)*kidx(3,k))) = k
    end do

    occ_max = 2d0; a2half = 0d0; pref = 1d0; expo = 2d0
    sigma = 0.3d0; tau = 0.1d0
    bmat = 0d0; bmat(1,1) = 1d0; bmat(2,2) = 1d0; bmat(3,3) = 1d0
    eps_inf = 12d0; qtf2 = 1d0; wp2 = 1d0; lambda2 = 0d0; q2reg = 0.1d0
    eval(iv, :) = 0d0                                   ! flat valence
    eval(ic, :) = (/ 1.0d0, 0.5d0, 0.5d0, 0.5d0 /)      ! k1 hot, others at the edge
    ecbm = 0.5d0; eth = 0.2d0
    f(iv, :) = 2d0                                      ! valence full
    f(ic, :) = (/ 1.0d0, 0d0, 0d0, 0d0 /)               ! one hot conduction carrier @ k1

    call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                        pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dpop)

    ! (1) trace conserved
    total = sum(dpop)
    call chk("trace conserved (sum dpop = 0)", total, 0d0, TOL)

    ! (2) carrier multiplication: conduction up, valence down
    dcb = sum(dpop(ic, :)); dvb = sum(dpop(iv, :))
    if (dcb <= 0d0) call bad("conduction population did not grow (no multiplication)")
    if (dvb >= 0d0) call bad("valence population did not drop (no hole created)")
    call chk("conduction gain == valence loss", dcb, -dvb, TOL)
    ! the hot primary (ic,k1) must LOSE population (it relaxes / ionizes)
    if (dpop(ic, 1) >= 0d0) call bad("hot primary (ic,k1) did not lose population")

    ! (3) populations stay in [0, occ_max]
    do k = 1, nk
        do a = 1, nba
            if (f(a,k)+dpop(a,k) < -1d-9 .or. f(a,k)+dpop(a,k) > occ_max+1d-9) &
                call bad("population left [0, occ_max]")
        end do
    end do

    ! (4) below threshold -> exact no-op
    block
        real(8) :: dz(nba, nk)
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, 10d0, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dz)
        call chk("below threshold -> no-op (max|dpop|)", maxval(abs(dz)), 0d0, 1d-12)
    end block

    ! (5) dynamic free-carrier screening lambda^2 > 0 (GaAs dyn_lambda path)
    ! must REDUCE the ionization rate monotonically (|V(q)|^2 ~ 1/(q^2+lambda^2))
    ! while staying exactly trace-conserving; lambda^2 -> inf kills the channel.
    block
        real(8) :: dlam(nba, nk), dinf(nba, nk)
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, 5d0, q2reg, sigma, tau, dlam)
        call chk("lambda2>0: trace conserved", sum(dlam), 0d0, TOL)
        if (sum(dlam(ic, :)) >= sum(dpop(ic, :))) &
            call bad("lambda2>0 did not reduce the ionization rate")
        if (sum(dlam(ic, :)) <= 0d0) &
            call bad("moderate lambda2 should not kill the channel entirely")
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, 1d30, q2reg, sigma, tau, dinf)
        call chk("lambda2->inf: channel off", maxval(abs(dinf)), 0d0, 1d-12)
    end block

    ! (6) i1_lo/i1_hi subrange additivity (the MPI rank split): the sum of the
    ! two half-range dpops must equal the full-range dpop to round-off.
    block
        real(8) :: dlo(nba, nk), dhi(nba, nk)
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dlo, &
                        i1_lo=1, i1_hi=2)
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dhi, &
                        i1_lo=3, i1_hi=4)
        call chk("rank-split additivity (max|full - lo - hi|)", &
                 maxval(abs(dpop - dlo - dhi)), 0d0, 1d-13)
    end block

    ! (7) B1 vq table: BIT-IDENTICAL to the direct evaluation
    block
        type(t_ring_opts) :: op
        real(8) :: dt2(nba, nk)
        allocate(op%vq_tab((2*kn(1)-1)*(2*kn(2)-1)*(2*kn(3)-1)))
        call build_vq_table(kn, bmat, eps_inf, qtf2, wp2, lambda2, q2reg, op%vq_tab)
        op%use_tab = .true.
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dt2, opts=op)
        call chk("B1 vq table bit-identical (max|diff|)", maxval(abs(dt2 - dpop)), 0d0, 1d-300)
        ! B3: a floor above the max kills the channel; tiny floor is a no-op
        op%vq_floor = 2d0 * maxval(op%vq_tab)
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dt2, opts=op)
        call chk("B3 floor > max -> channel off", maxval(abs(dt2)), 0d0, 1d-14)
    end block

    ! (8) A5 Franz-Keldysh softening: sub-threshold state ionizes, trace exact
    block
        type(t_ring_opts) :: op
        real(8) :: dt2(nba, nk)
        op%fk_theta = 0.3d0
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, 0.6d0, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dt2, opts=op)
        if (sum(dt2(ic, :)) <= 0d0) call bad("A5: FK tail did not activate the sub-threshold state")
        call chk("A5 trace conserved", sum(dt2), 0d0, TOL)
        ! hard threshold (no opts) at eth=0.6 must stay dark (ekin=0.5)
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, 0.6d0, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dt2)
        call chk("A5 control: hard threshold dark", maxval(abs(dt2)), 0d0, 1d-14)
    end block

    ! (9) A1 phonon sidebands: an off-shell quadruple becomes reachable
    block
        type(t_ring_opts) :: op
        real(8) :: dt2(nba, nk), dt3(nba, nk)
        ! narrow sigma so the direct delta misses everything
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, 0.02d0, tau, dt2)
        op%phassist = 1d0;  op%nph = 1
        allocate(op%hw(1), op%nbb(1), op%wrel(1))
        op%hw = 0.5d0;  op%nbb = 0.1d0;  op%wrel = 1d0   ! sideband at the 0.5 mismatch
        call ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, 0.02d0, tau, dt3, opts=op)
        if (sum(dt3(ic, :)) <= sum(dt2(ic, :)) + 1d-14) &
            call bad("A1: phonon sideband did not open the off-shell channel")
        call chk("A1 trace conserved", sum(dt3), 0d0, TOL)
    end block

    ! (10) A2 hole-initiated channel: deep hole -> pair created, trace exact
    block
        type(t_ring_opts) :: op
        real(8) :: dt2(nba, nk), ev2(nba, nk), f2(nba, nk)
        ev2(iv, :) = 0d0;  ev2(iv, 1) = -1.0d0        ! k1 hosts a DEEP valence state
        ev2(ic, :) = 0.5d0
        f2(iv, :) = 2d0;  f2(iv, 1) = 1.0d0           ! deep hole at k1
        f2(ic, :) = 0d0
        op%pref_h = 1d0;  op%evbm = 0d0
        call ii_interk_dpop(nk, nba, ev2, f2, occ_max, a2half, ecbm, 0.2d0, &
                            pref, expo, iv, ic, kidx, kn, klut, &
                        bmat, eps_inf, qtf2, wp2, lambda2, q2reg, 0.3d0, tau, dt2, opts=op)
        ! note: the ELECTRON channel is dark here (empty CB), only holes act
        call chk("A2 trace conserved", sum(dt2), 0d0, TOL)
        if (dt2(iv, 1) <= 0d0) call bad("A2: deep hole was not filled")
        if (sum(dt2(ic, :)) <= 0d0) call bad("A2: no pair electron created (hhe)")
    end block

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (inter-k impact ionization: trace-conserving, '// &
                       'carrier-multiplying, bounded, no-op below threshold; B1 table bit-identical, '// &
                       'B3 floor, A5 FK tail, A1 sidebands, A2 hole channel)'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains
    subroutine chk(name, got, want, tol)
        character(*), intent(in) :: name
        real(8), intent(in) :: got, want, tol
        if (abs(got-want) > tol*max(1d0,abs(want))) then
            write(*,'(a,a,a,es16.8,a,es16.8)') '  FAIL: ',name,' got=',got,' want=',want
            nfail = nfail + 1
        end if
    end subroutine chk
    subroutine bad(name)
        character(*), intent(in) :: name
        write(*,'(a,a)') '  FAIL: ', name
        nfail = nfail + 1
    end subroutine bad
end program test_ii_interk_cptp
