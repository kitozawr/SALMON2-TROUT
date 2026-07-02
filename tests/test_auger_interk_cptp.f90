!
!  test_auger_interk_cptp.f90  -  momentum-conserving inter-k Auger recombination
!  = the exact time-reverse of the nonlocal impact ionization (detailed balance).
!
!  Tests auger_interk_dpop: two conduction e- (k1',c),(k2',c) + a hole (k2,v) ->
!  one recombines into the hole, the other is promoted to the hot state (k1,c),
!  with the SAME momentum map, |V(q)|^2 and energy broadening as ii_interk_dpop
!  and REVERSED occupation factors. The net dpop must:
!    1) conserve the total trace EXACTLY (sum(dpop)=0) -- CPTP,
!    2) RECOMBINE: conduction population drops, valence gains (hole filled),
!       and the hot destination (ic,k1) gains (promotion),
!    3) keep populations in [0, occ_max],
!    4) be a no-op with no conduction pair (empty CB) or below threshold,
!    5) DETAILED BALANCE: for Fermi-Dirac occupations (same T, mu) on an
!       energy-conserving quadruple set, the linear-regime net
!       dpop_II + dpop_Auger must vanish identically (equilibrium fixed point)
!       -- the defining property of the time-reversed pair [Rana 2007 Sec 2].
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_auger_interk_cptp
    use sbe_superres_ssbe, only: ii_interk_dpop, auger_interk_dpop
    implicit none
    integer, parameter :: nk = 4, nba = 2, iv = 1, ic = 2
    integer :: kn(3), kidx(3, nk), klut(0:nk-1), nfail, a, k
    real(8) :: eval(nba, nk), f(nba, nk), dpop(nba, nk)
    real(8) :: occ_max, a2half, ecbm, eth, pref, expo, kappa2, sigma, tau
    real(8) :: dcb, dvb, total
    real(8), parameter :: TOL = 1d-11

    nfail = 0
    kn = (/ 2, 2, 1 /)
    kidx(:,1) = (/0,0,0/); kidx(:,2) = (/1,0,0/)
    kidx(:,3) = (/0,1,0/); kidx(:,4) = (/1,1,0/)
    do k = 1, nk
        klut(kidx(1,k) + kn(1)*(kidx(2,k) + kn(2)*kidx(3,k))) = k
    end do

    occ_max = 2d0; a2half = 0d0; pref = 1d0; expo = 2d0
    kappa2 = 0.1d0; sigma = 0.3d0; tau = 0.1d0
    eval(iv, :) = 0d0                                   ! flat valence
    eval(ic, :) = (/ 1.0d0, 0.5d0, 0.5d0, 0.5d0 /)      ! k1 hot, others at the edge
    ecbm = 0.5d0; eth = 0.2d0

    ! --- forward Auger event: CB pair + a hole -> recombination + promotion ----
    f(iv, :) = (/ 2d0, 1.0d0, 2d0, 2d0 /)               ! a hole at k2
    f(ic, :) = (/ 0d0, 1.0d0, 1.0d0, 1.0d0 /)           ! CB pair at the edge; hot empty
    call auger_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                           pref, expo, iv, ic, kidx, kn, klut, kappa2, sigma, tau, dpop)

    ! (1) trace conserved
    total = sum(dpop)
    call chk("trace conserved (sum dpop = 0)", total, 0d0, TOL)

    ! (2) recombination: conduction down, valence up, hot destination gains
    dcb = sum(dpop(ic, :)); dvb = sum(dpop(iv, :))
    if (dcb >= 0d0) call bad("conduction population did not drop (no recombination)")
    if (dvb <= 0d0) call bad("valence population did not gain (hole not filled)")
    call chk("valence gain == conduction loss", dvb, -dcb, TOL)
    if (dpop(ic, 1) <= 0d0) call bad("hot destination (ic,k1) did not gain (no promotion)")

    ! (3) populations stay in [0, occ_max]
    do k = 1, nk
        do a = 1, nba
            if (f(a,k)+dpop(a,k) < -1d-9 .or. f(a,k)+dpop(a,k) > occ_max+1d-9) &
                call bad("population left [0, occ_max]")
        end do
    end do

    ! (4a) empty conduction band -> exact no-op (no pair to recombine)
    block
        real(8) :: dz(nba, nk)
        real(8) :: f0(nba, nk)
        f0(iv, :) = 2d0;  f0(ic, :) = 0d0
        call auger_interk_dpop(nk, nba, eval, f0, occ_max, a2half, ecbm, eth, &
                               pref, expo, iv, ic, kidx, kn, klut, kappa2, sigma, tau, dz)
        call chk("empty CB -> no-op (max|dpop|)", maxval(abs(dz)), 0d0, 1d-12)
    end block
    ! (4b) below threshold (huge eth) -> exact no-op
    block
        real(8) :: dz(nba, nk)
        call auger_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, 10d0, &
                               pref, expo, iv, ic, kidx, kn, klut, kappa2, sigma, tau, dz)
        call chk("below threshold -> no-op (max|dpop|)", maxval(abs(dz)), 0d0, 1d-12)
    end block

    ! (5) DETAILED BALANCE: Fermi-Dirac occupations at (kT=0.25, mu=0.25) on the
    ! energy-conserving quadruple set (1.0 + 0.0 = 0.5 + 0.5); a NARROW sigma
    ! suppresses the off-shell (broadening) terms, and a tiny tau puts both
    ! kernels in the linear regime where the caps f*(Gamma tau) / room*(Gamma tau)
    ! realize the exact f1 f2 (1-f3)(1-f4) = (1-f1)(1-f2) f3 f4 FD identity.
    block
        real(8) :: dpi(nba, nk), dpa(nba, nk), fd(nba, nk), scale
        real(8) :: kT, mu
        integer :: kk, aa
        kT = 0.25d0;  mu = 0.25d0
        do kk = 1, nk
            do aa = 1, nba
                fd(aa, kk) = occ_max / (exp((eval(aa, kk) - mu) / kT) + 1d0)
            end do
        end do
        call ii_interk_dpop(nk, nba, eval, fd, occ_max, a2half, ecbm, eth, &
                            pref, expo, iv, ic, kidx, kn, klut, kappa2, 0.02d0, 1d-6, dpi)
        call auger_interk_dpop(nk, nba, eval, fd, occ_max, a2half, ecbm, eth, &
                               pref, expo, iv, ic, kidx, kn, klut, kappa2, 0.02d0, 1d-6, dpa)
        scale = max(maxval(abs(dpi)), 1d-300)
        if (scale < 1d-30) call bad("detailed balance: II kernel inert (test not exercising)")
        call chk("detailed balance: max|dpop_II + dpop_Auger| / max|dpop_II|", &
                 maxval(abs(dpi + dpa)) / scale, 0d0, 1d-6)
    end block

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (inter-k Auger: trace-conserving, recombining+promoting, '// &
                       'bounded, no-pair/sub-threshold no-op, FD detailed balance with the II kernel)'
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
end program test_auger_interk_cptp
