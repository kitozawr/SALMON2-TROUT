!
!  test_vg_basis_nb.f90  -  VG basis sufficiency & N_b convergence primitives.
!
!  Encodes the three convergence criteria from the wiki page "VG Basis
!  Sufficiency & N_b Convergence":
!    (a) vg_ptop_exceeds  - top-band occupation threshold logic;
!    (b) vg_conv_error    - relative L2 convergence metric between N_b and N_b+D;
!    (c) vg_eta_admixture - dimensionless admixture of the first discarded band;
!        vg_trunc_shift2  - 2nd-order truncation shift the Houston basis inherits.
!
!  The load-bearing physics checks (no LAPACK needed -- a tiny self-contained
!  Jacobi eigensolver is included):
!    * Hylleraas-Undheim-MacDonald interlacing / upper-bound: eigenvalues of the
!      truncated N_b principal submatrix of a Hermitian H are >= the true ones,
!      and decrease monotonically as N_b grows. This is WHY truncation can only
!      raise levels, and why the Houston basis cannot cure an insufficient N_b.
!    * vg_trunc_shift2 reproduces the actual lowest-eigenvalue shift of the
!      projected H_VG = diag(eps) + A*pi in the perturbative (eta << 1) regime,
!      and the residual grows ~A^4 (confirming it is the leading-order error).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_vg_basis_nb
    use sbe_superres_ssbe
    implicit none
    integer :: nfail
    real(8), parameter :: TOL = 1d-9

    nfail = 0

    call test_eta()
    call test_conv_error()
    call test_ptop()
    call test_interlacing()
    call test_shift_formula()

    if (nfail == 0) then
        write(*,'(a)') 'PASS'; call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains

    ! ---- criterion (c): admixture parameter ---------------------------------
    subroutine test_eta()
        real(8) :: eta
        ! eta = A_max |pi| / |gap|
        eta = vg_eta_admixture(0.85d0, 0.97d0, 2.0d0)
        call chk("eta value", eta, 0.85d0*0.97d0/2.0d0, TOL)
        ! larger gap -> smaller eta (safer to discard)
        if (.not. (vg_eta_admixture(0.85d0,0.97d0,20d0) < &
                   vg_eta_admixture(0.85d0,0.97d0,2d0))) call bad("eta not monotone in gap")
        ! degenerate gap -> nonperturbative (huge)
        if (vg_eta_admixture(1d0,1d0,0d0) < 1d30) call bad("eta(zero gap) not huge")
    end subroutine test_eta

    ! ---- criterion (b): relative L2 convergence metric ----------------------
    subroutine test_conv_error()
        real(8) :: a(4), b(4), e
        a = (/ 1d0, 2d0, 3d0, 4d0 /)
        ! identical observables -> zero error
        e = vg_conv_error(a, a, 4)
        call chk("conv_error self = 0", e, 0d0, TOL)
        ! known case: b = a scaled, error = ||b-a||/||b||
        b = 1.1d0 * a
        call chk("conv_error scaled", vg_conv_error(a, b, 4), &
                 sqrt(sum((b-a)**2)/sum(b**2)), TOL)
        ! a larger discrepancy gives a larger metric
        if (.not. (vg_conv_error(a, 2d0*a, 4) > vg_conv_error(a, 1.1d0*a, 4))) &
            call bad("conv_error not monotone")
    end subroutine test_conv_error

    ! ---- criterion (a): top-band threshold ----------------------------------
    subroutine test_ptop()
        if (.not. vg_ptop_exceeds(2d-3, 1d-3)) call bad("ptop should exceed")
        if (vg_ptop_exceeds(1d-4, 1d-3))       call bad("ptop should not exceed")
    end subroutine test_ptop

    ! ---- Hylleraas-Undheim-MacDonald interlacing / upper bound --------------
    subroutine test_interlacing()
        integer, parameter :: n = 6
        real(8) :: H(n,n), evfull(n), evtrunc(5), evtrunc4(4)
        integer :: i, j
        ! a fixed symmetric "H_VG-like" matrix (diagonal trend + off-diagonal coupling)
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    H(i,j) = dble(i)              ! ascending diagonal (eps_n)
                else
                    H(i,j) = 0.30d0 / dble(abs(i-j)+1)   ! interband coupling (pi-like)
                end if
            end do
        end do
        call sym_eig(H, n, evfull)
        ! truncate to the leading 5x5 and 4x4 principal submatrices
        call sym_eig(H(1:5,1:5), 5, evtrunc)
        call sym_eig(H(1:4,1:4), 4, evtrunc4)
        ! (i) upper bound: each truncated level >= the corresponding true level
        do i = 1, 5
            if (evtrunc(i) < evfull(i) - TOL) call bad("truncated level below true (5)")
        end do
        do i = 1, 4
            if (evtrunc4(i) < evfull(i) - TOL) call bad("truncated level below true (4)")
        end do
        ! (ii) interlacing/MacDonald: a smaller basis can only raise a level,
        ! so eps^(4)_i >= eps^(5)_i >= eps^(6)_i for the retained lowest levels.
        do i = 1, 4
            if (evtrunc4(i) < evtrunc(i) - TOL) call bad("4-basis below 5-basis")
            if (evtrunc(i)  < evfull(i)  - TOL) call bad("5-basis below 6-basis")
        end do
        ! the lowest level is the most accurate; the top retained one drifts most
        if (.not. ((evtrunc4(4)-evfull(4)) >= (evtrunc4(1)-evfull(1)) - TOL)) &
            call bad("top retained level not the most corrupted")
    end subroutine test_interlacing

    ! ---- 2nd-order shift formula vs the true projected-eigenvalue shift ------
    subroutine test_shift_formula()
        ! Level a = lowest band; couple it to two discarded bands via A*pi.
        ! H = diag(eps) + A*pi (real symmetric).  Compare the EXACT shift of the
        ! lowest eigenvalue (truncate-to-1x1 reference = eps_a) to vg_trunc_shift2
        ! in the perturbative regime, and confirm the residual scales ~A^4.
        integer, parameter :: m = 3
        real(8) :: eps(m), pivec(2), eps_disc(2)
        real(8) :: A, sh_pred1, sh_pred2, sh_true1, sh_true2, res1, res2, ratio
        eps      = (/ 0.0d0, 5.0d0, 8.0d0 /)   ! a + 2 discarded bands
        pivec    = (/ 0.7d0, 0.4d0 /)          ! pi_{a,c}
        eps_disc = (/ 5.0d0, 8.0d0 /)
        ! small field
        A = 0.02d0
        sh_pred1 = vg_trunc_shift2(A, pivec, eps(1), eps_disc, 2)
        sh_true1 = lowest_shift(eps, pivec, A)
        res1 = abs(sh_true1 - sh_pred1)
        ! the 2nd-order formula must match to better than 1% of the shift here
        if (res1 > 1d-2 * abs(sh_true1)) call bad("shift formula off in perturbative regime")
        ! halve the field: leading error ~A^4 -> residual drops by ~16x
        A = 0.01d0
        sh_pred2 = vg_trunc_shift2(A, pivec, eps(1), eps_disc, 2)
        sh_true2 = lowest_shift(eps, pivec, A)
        res2 = abs(sh_true2 - sh_pred2)
        ratio = res1 / max(res2, 1d-300)
        if (ratio < 8d0) call bad("truncation residual does not scale ~A^4")
        ! sign: coupling to HIGHER bands pushes the level DOWN (negative shift)
        if (sh_pred1 >= 0d0) call bad("shift to higher bands not negative")
    end subroutine test_shift_formula

    ! exact shift of the lowest eigenvalue of diag(eps)+A*pi relative to eps(1)
    function lowest_shift(eps, pivec, A) result(sh)
        real(8), intent(in) :: eps(3), pivec(2), A
        real(8) :: sh, H(3,3), ev(3)
        H = 0d0
        H(1,1) = eps(1); H(2,2) = eps(2); H(3,3) = eps(3)
        H(1,2) = A*pivec(1); H(2,1) = H(1,2)
        H(1,3) = A*pivec(2); H(3,1) = H(1,3)
        call sym_eig(H, 3, ev)
        sh = ev(1) - eps(1)
    end function lowest_shift

    ! --- minimal cyclic Jacobi eigensolver for a real symmetric matrix -------
    ! Returns eigenvalues in ascending order. Self-contained (no LAPACK).
    subroutine sym_eig(Ain, n, ev)
        integer, intent(in) :: n
        real(8), intent(in) :: Ain(n,n)
        real(8), intent(out) :: ev(n)
        real(8) :: A(n,n), theta, t, c, s, tau, off, aii, ajj, aij, tmp
        integer :: p, q, i, sweep
        A = Ain
        do sweep = 1, 100
            off = 0d0
            do p = 1, n-1
                do q = p+1, n
                    off = off + A(p,q)**2
                end do
            end do
            if (off < 1d-30) exit
            do p = 1, n-1
                do q = p+1, n
                    if (abs(A(p,q)) < 1d-300) cycle
                    theta = (A(q,q) - A(p,p)) / (2d0*A(p,q))
                    t = sign(1d0, theta) / (abs(theta) + sqrt(theta*theta + 1d0))
                    c = 1d0 / sqrt(t*t + 1d0)
                    s = t * c
                    tau = s / (1d0 + c)
                    aij = A(p,q)
                    aii = A(p,p); ajj = A(q,q)
                    A(p,p) = aii - t*aij
                    A(q,q) = ajj + t*aij
                    A(p,q) = 0d0; A(q,p) = 0d0
                    do i = 1, n
                        if (i /= p .and. i /= q) then
                            tmp    = A(i,p)
                            A(i,p) = tmp - s*(A(i,q) + tau*tmp)
                            A(p,i) = A(i,p)
                            A(i,q) = A(i,q) + s*(tmp - tau*A(i,q))
                            A(q,i) = A(i,q)
                        end if
                    end do
                end do
            end do
        end do
        do i = 1, n
            ev(i) = A(i,i)
        end do
        call sort_asc(ev, n)
    end subroutine sym_eig

    subroutine sort_asc(x, n)
        integer, intent(in) :: n
        real(8), intent(inout) :: x(n)
        integer :: i, j
        real(8) :: tmp
        do i = 1, n-1
            do j = i+1, n
                if (x(j) < x(i)) then
                    tmp = x(i); x(i) = x(j); x(j) = tmp
                end if
            end do
        end do
    end subroutine sort_asc

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

end program test_vg_basis_nb
