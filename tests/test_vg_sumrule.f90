!
!  test_vg_sumrule.f90 - velocity-gauge f-sum-rule completeness of a truncated
!  basis (wiki/12 sec. 6): eta_a = 1 - <sum_{m/=n} 2|p_nm^a|^2/(eps_m-eps_n)>
!  over the occupied states = the fraction of the diamagnetic current A N_e/V
!  that the nb-band basis cannot cancel (a spurious reactive current ~ E/omega).
!
!  Checks (synthetic 3-band model, k-independent, one filled band occ = 2):
!    1) p12x = 0.2, de12 = 0.2 -> 0.40 ; p13x = 0.1, de13 = 0.5 -> 0.04 :
!       S_x = 0.44, eta_x = 0.56 ; nothing along y -> eta_y = 1 ; z likewise.
!    2) unequal k-weights and an unoccupied band do not change the result
!       (the sum is occupation- and weight-normalized).
!    3) a degenerate pair (|de| < de_min) is skipped, not divided by zero.
!    4) a "complete" toy basis (S = 1 by construction) gives eta = 0.
!    5) vg_gs_current_k (the pure-gauge restoration actually applied in
!       calc_current_bloch): on an analytic 2-level H(A) = eps + A p_x,
!       J_gs(A) = occ [dE_0/dA + A] (Hellmann-Feynman) at A = 1e-5, 0.05, 0.3
!       against central finite differences of the exact eigenvalue, the linear
!       limit J_gs/(occ A) -> 1 - S = eta, and the nonlinear departure from it.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_vg_sumrule
    use sbe_superres_ssbe, only: vg_sumrule_eta, vg_gs_current_k
    implicit none
    integer, parameter :: NK = 3, NB = 3
    real(8) :: eigen(NB, NK), occ(NB, NK), w(NK), eta(3), scap(3)
    complex(8) :: p(NB, NB, 3, NK)
    integer :: nfail, ik
    ! (5) two-level pure-gauge reference
    real(8), parameter :: d0 = 0.04d0, pvv = 0.10d0, pcc = -0.05d0, h = 1d-6
    complex(8), parameter :: q = (0.06d0, 0.04d0)          ! 2|q|^2/d0 = 0.26 -> eta_lin = 0.74
    complex(8) :: W2(2, 2), p2(2, 2, 3)
    real(8) :: A3(3), occ2(2), jgs(3), e0p, e0m, jhf, eta_lin, slope_a
    real(8), parameter :: a_list(3) = (/ 1d-5, 0.05d0, 0.3d0 /)
    integer :: ia
    nfail = 0

    ! --- (1)+(2) ---------------------------------------------------------------
    eigen(1, :) = 0d0;  eigen(2, :) = 0.2d0;  eigen(3, :) = 0.5d0
    occ = 0d0;  occ(1, :) = 2d0
    w = (/ 0.2d0, 0.5d0, 0.3d0 /)
    p = (0d0, 0d0)
    do ik = 1, NK
        p(1, 2, 1, ik) = (0.2d0, 0d0);  p(2, 1, 1, ik) = (0.2d0, 0d0)
        p(1, 3, 1, ik) = (0d0, 0.1d0);  p(3, 1, 1, ik) = (0d0, -0.1d0)
        p(2, 3, 1, ik) = (0.7d0, 0d0);  p(3, 2, 1, ik) = (0.7d0, 0d0)   ! unoccupied-unoccupied: must not count
    end do
    call vg_sumrule_eta(NK, NB, eigen, p, occ, w, 1d-8, eta, scap)
    write(*,'(a,3f9.5)') '  S   = ', scap
    write(*,'(a,3f9.5)') '  eta = ', eta
    if (abs(scap(1) - 0.44d0) > 1d-12) call bad('S_x /= 0.44')
    if (abs(eta(1) - 0.56d0) > 1d-12) call bad('eta_x /= 0.56')
    if (abs(eta(2) - 1d0) > 1d-12 .or. abs(eta(3) - 1d0) > 1d-12) call bad('eta_y/z /= 1 with no coupling')

    ! --- (3) degenerate pair skipped -------------------------------------------
    eigen(3, 2) = eigen(1, 2)                      ! band 3 degenerate with band 1 at k = 2
    call vg_sumrule_eta(NK, NB, eigen, p, occ, w, 1d-8, eta, scap)
    if (.not. (scap(1) == scap(1))) call bad('NaN from a degenerate pair')          ! NaN check
    if (abs(scap(1) - (0.40d0 + 0.04d0 * (1d0 - w(2)))) > 1d-12) call bad('degenerate pair not skipped correctly')

    ! --- (4) complete toy basis: S = 1 ------------------------------------------
    eigen(1, :) = 0d0;  eigen(2, :) = 0.5d0;  eigen(3, :) = 1d0
    p = (0d0, 0d0)
    do ik = 1, NK
        p(1, 2, 1, ik) = cmplx(sqrt(0.5d0 * 0.5d0 / 2d0), 0d0, 8)   ! 2|p|^2/0.5 = 0.5
        p(1, 3, 1, ik) = cmplx(sqrt(0.5d0 * 1.0d0 / 2d0), 0d0, 8)   ! 2|p|^2/1.0 = 0.5
    end do
    call vg_sumrule_eta(NK, NB, eigen, p, occ, w, 1d-8, eta, scap)
    if (abs(eta(1)) > 1d-12) call bad('complete toy basis must give eta_x = 0')

    ! --- (5) pure-gauge restoration: adiabatic ground-state current ------------
    ! H(A) = [[A pvv, A q], [A q*, d0 + A pcc]] ; exact lowest eigenvector from the 2x2 closed form.
    p2 = (0d0, 0d0)
    p2(1, 1, 1) = cmplx(pvv, 0d0, 8);  p2(2, 2, 1) = cmplx(pcc, 0d0, 8)
    p2(1, 2, 1) = q;                   p2(2, 1, 1) = conjg(q)
    occ2 = (/ 2d0, 0d0 /)
    eta_lin = 1d0 - 2d0 * abs(q)**2 / d0                      ! = 1 - S for this pair
    do ia = 1, 3
        A3 = 0d0
        A3(1) = a_list(ia)
        call eig2(A3(1), W2)
        call vg_gs_current_k(2, W2, p2, A3, occ2, jgs)
        e0p = e0_of(A3(1) + h);  e0m = e0_of(A3(1) - h)
        jhf = occ2(1) * ((e0p - e0m) / (2d0 * h) + A3(1))    ! Hellmann-Feynman: occ [dE_0/dA + A]
        slope_a = (jgs(1) - occ2(1) * pvv) / (occ2(1) * A3(1))  ! subtract the A = 0 band velocity
        write(*,'(a,f7.4,a,es14.6,a,es14.6,a,f9.5,a,f9.5)') '  A =', A3(1), '  J_gs =', jgs(1), &
            '  occ[dE0/dA + A] =', jhf, '   (J_gs - occ pvv)/(occ A) =', slope_a, '   eta_lin =', eta_lin
        if (abs(jgs(1) - jhf) > 1d-7) call bad('J_gs /= Hellmann-Feynman derivative')
        if (abs(jgs(2)) > 1d-14 .or. abs(jgs(3)) > 1d-14) call bad('J_gs must vanish along uncoupled axes')
        if (ia == 1 .and. abs(slope_a - eta_lin) > 2d-3) call bad('linear limit of J_gs is not eta N A')
        if (ia == 3 .and. abs(slope_a - eta_lin) < 1d-2) call bad('J_gs at large A must depart from the linear eta form')
    end do

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (VG f-sum rule: captured strength, eta, weights/occupation normalization, degenerate skip,'// &
                       ' completeness; pure-gauge J_gs = Hellmann-Feynman, linear limit eta, nonlinear beyond)'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if
contains
    pure function e0_of(a) result(e0)      ! exact lowest eigenvalue of the 2x2 H(A)
        real(8), intent(in) :: a
        real(8) :: e0, hm, hd
        hm = 0.5d0 * (a * pvv + d0 + a * pcc);  hd = 0.5d0 * (a * pvv - d0 - a * pcc)
        e0 = hm - sqrt(hd**2 + a**2 * abs(q)**2)
    end function e0_of
    subroutine eig2(a, W)                  ! eigenvectors (columns, ascending) of the 2x2 H(A)
        real(8), intent(in) :: a
        complex(8), intent(out) :: W(2, 2)
        real(8) :: hd, r, nrm
        complex(8) :: off, v1, v2
        hd = 0.5d0 * (a * pvv - d0 - a * pcc);  off = a * q
        r = sqrt(hd**2 + abs(off)**2)
        ! lowest: (H - e0) v = 0 with H11 - e0 = hd + r  ->  v = (-off, hd + r) up to normalisation
        v1 = -off;  v2 = cmplx(hd + r, 0d0, 8)
        nrm = sqrt(abs(v1)**2 + abs(v2)**2);  W(1, 1) = v1 / nrm;  W(2, 1) = v2 / nrm
        ! highest: orthogonal complement
        W(1, 2) = -conjg(W(2, 1));  W(2, 2) = conjg(W(1, 1))
    end subroutine eig2
    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_vg_sumrule
