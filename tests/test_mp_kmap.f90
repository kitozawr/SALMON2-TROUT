!
!  test_mp_kmap.f90  -  Monkhorst-Pack momentum-conservation index map.
!
!  Tests mp_grid_triple / mp_partner_triple used by the inter-k / nonlocal
!  collision channels (impact ionization, e-e). On a regular MP mesh:
!    * mp_grid_triple recovers the integer index triple from the reduced
!      coordinate exactly (residual ~ 0 for a true MP point),
!    * mp_partner_triple gives k2' = k1 + k2 - k1' (mod G) by integer arithmetic,
!      so crystal momentum is conserved: k1 + k2 == k1' + k2' (mod 1).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_mp_kmap
    use sbe_superres_ssbe, only: mp_grid_triple, mp_partner_triple
    implicit none
    integer, parameter :: n(3) = (/ 4, 4, 2 /)
    integer :: nk, i, j, l, ik, nfail, a, b, c, ma(3), mb(3), mc(3), m2p(3)
    integer, allocatable :: idx(:,:)
    real(8), allocatable :: kred(:,:)
    real(8) :: resid, kc(3), ksum(3), dd
    real(8), parameter :: TOL = 1d-12

    nfail = 0
    nk = n(1) * n(2) * n(3)
    allocate(kred(3, nk), idx(3, nk))

    ! Build the standard MP mesh k_m = (2m - n + 1)/(2n) and verify the round-trip.
    ik = 0
    do i = 0, n(1)-1
      do j = 0, n(2)-1
        do l = 0, n(3)-1
            ik = ik + 1
            kred(1, ik) = (2d0*i - n(1) + 1d0)/(2d0*n(1))
            kred(2, ik) = (2d0*j - n(2) + 1d0)/(2d0*n(2))
            kred(3, ik) = (2d0*l - n(3) + 1d0)/(2d0*n(3))
            call mp_grid_triple(kred(:, ik), n, idx(:, ik), resid)
            if (resid > TOL) call bad("mp_grid_triple residual nonzero on an MP point")
            if (idx(1,ik) /= i .or. idx(2,ik) /= j .or. idx(3,ik) /= l) &
                call bad("mp_grid_triple recovered the wrong index triple")
        end do
      end do
    end do

    ! Momentum conservation for a spread of quadruples (k1=a, k2=b, k1'=c).
    do a = 1, nk, 3
      do b = 1, nk, 5
        do c = 1, nk, 7
            ma = idx(:, a); mb = idx(:, b); mc = idx(:, c)
            call mp_partner_triple(ma, mb, mc, n, m2p)
            ! reconstruct k2' from its triple and check k1 + k2 == k1' + k2' (mod 1)
            kc(1) = (2d0*m2p(1) - n(1) + 1d0)/(2d0*n(1))
            kc(2) = (2d0*m2p(2) - n(2) + 1d0)/(2d0*n(2))
            kc(3) = (2d0*m2p(3) - n(3) + 1d0)/(2d0*n(3))
            ksum = kred(:,a) + kred(:,b) - kred(:,c) - kc
            dd = maxval(abs(ksum - anint(ksum)))      ! 0 mod 1 ?
            if (dd > 1d-12) call bad("momentum not conserved: k1+k2 /= k1'+k2' (mod G)")
            ! partner index must be a valid grid point in range
            if (any(m2p < 0) .or. m2p(1) >= n(1) .or. m2p(2) >= n(2) .or. m2p(3) >= n(3)) &
                call bad("partner triple out of range")
        end do
      end do
    end do

    ! A non-MP point (random offset) must report a nonzero residual (gate signal).
    block
        integer :: mm(3)
        call mp_grid_triple((/ 0.13d0, 0.0d0, 0.0d0 /), n, mm, resid)
        if (resid <= 1d-6) call bad("non-MP point reported as exact (residual ~ 0)")
    end block

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (MP momentum map: exact index round-trip, '// &
                       'k1+k2=k1''+k2'' mod G, non-MP residual flagged)'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains
    subroutine bad(name)
        character(*), intent(in) :: name
        write(*,'(a,a)') '  FAIL: ', name
        nfail = nfail + 1
    end subroutine bad
end program test_mp_kmap
