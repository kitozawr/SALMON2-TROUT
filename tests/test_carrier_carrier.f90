!
!  test_carrier_carrier.f90  -  Part F carrier-carrier thermalization (Part F).
!
!  Tests carrier_carrier_relax (and the fit_fermi_dirac it uses): the CPTP map
!  must conserve total NUMBER and ENERGY of the adiabatic populations exactly,
!  damp the coherences by (1-alpha), keep populations in [0,occ], drive the
!  distribution toward the fitted Fermi-Dirac, and be a no-op when the level set
!  is empty/full or not Fermi-Dirac representable (population inversion).
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_carrier_carrier
    use sbe_superres_ssbe
    implicit none
    integer, parameter :: n = 6
    integer :: nfail, a, b, it
    real(8) :: eps(n), occ, nu, tau, N0, E0, N1, E1, beta, mu, ftgt(n)
    complex(8) :: rho(n,n)
    logical :: ok
    real(8), parameter :: TOL = 1d-9

    nfail = 0
    occ = 2d0; nu = 5d-3; tau = 0.5d0
    eps = (/ -0.20d0, -0.10d0, -0.02d0, 0.02d0, 0.10d0, 0.20d0 /)

    ! a non-thermal Hermitian rho with coherences, populations in [0,occ]
    rho = (0d0,0d0)
    rho(1,1)=(1.8d0,0d0); rho(2,2)=(1.2d0,0d0); rho(3,3)=(1.5d0,0d0)
    rho(4,4)=(0.3d0,0d0); rho(5,5)=(0.6d0,0d0); rho(6,6)=(0.1d0,0d0)
    rho(1,3)=(0.2d0,0.1d0); rho(3,1)=conjg(rho(1,3))
    rho(2,5)=(0.15d0,-0.05d0); rho(5,2)=conjg(rho(2,5))

    call moments(rho, eps, n, N0, E0)
    call carrier_carrier_relax(n, rho, eps, occ, nu, tau)
    call moments(rho, eps, n, N1, E1)

    call chk("number conserved", N1, N0, 1d-9)
    call chk("energy conserved", E1, E0, 1d-9)
    ! coherence damped by (1-alpha)
    block
        real(8) :: alpha
        alpha = 1d0 - exp(-nu*tau)
        call chk("coherence damped (1-alpha)", abs(rho(1,3)), &
                 (1d0-alpha)*abs(cmplx(0.2d0,0.1d0, 8)), 1d-9)
    end block
    ! populations stay in [0,occ] and Hermiticity preserved
    do a = 1, n
        if (real(rho(a,a)) < -TOL .or. real(rho(a,a)) > occ+TOL) call bad("population left [0,occ]")
    end do
    if (abs(rho(1,3)-conjg(rho(3,1))) > TOL) call bad("Hermiticity broken")

    ! repeated application drives toward the fixed FD point (number/energy fixed).
    ! Use a fast rate so the geometric relaxation (1-alpha)^it -> 0 converges.
    do it = 1, 200
        call carrier_carrier_relax(n, rho, eps, occ, 2d0, tau)
    end do
    call moments(rho, eps, n, N1, E1)
    call chk("number conserved (converged)", N1, N0, 1d-7)
    call chk("energy conserved (converged)", E1, E0, 1d-7)
    ! at the fixed point the diagonal equals occ*f_FD fitted to (N0,E0)
    call fit_fermi_dirac(n, eps, N0/occ, E0/occ, beta, mu, ftgt, ok)
    if (.not. ok) call bad("FD fit failed for the test moments")
    do a = 1, n
        call chk("converged to FD", real(rho(a,a)), occ*ftgt(a), 1d-5)
    end do

    ! FD fit reproduces its own moments
    block
        real(8) :: Nf, Ef
        Nf = 0d0; Ef = 0d0
        do a = 1, n
            Nf = Nf + ftgt(a); Ef = Ef + eps(a)*ftgt(a)
        end do
        call chk("FD fit matches N", Nf, N0/occ, 1d-9)
        call chk("FD fit matches E", Ef, E0/occ, 1d-9)
    end block

    ! no-op when populations are not FD-representable (inversion: all weight high)
    rho = (0d0,0d0)
    rho(6,6)=(2d0,0d0); rho(5,5)=(2d0,0d0)   ! only the highest levels filled
    call moments(rho, eps, n, N0, E0)
    call carrier_carrier_relax(n, rho, eps, occ, nu, tau)
    call moments(rho, eps, n, N1, E1)
    call chk("inversion -> no-op (number)", N1, N0, TOL)
    call chk("inversion -> no-op (energy)", E1, E0, TOL)

    if (nfail == 0) then
        write(*,'(a)') 'PASS'; call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains
    subroutine moments(r, e, m, Nout, Eout)
        integer, intent(in) :: m
        complex(8), intent(in) :: r(m,m)
        real(8), intent(in) :: e(m)
        real(8), intent(out) :: Nout, Eout
        integer :: a
        Nout = 0d0; Eout = 0d0
        do a = 1, m
            Nout = Nout + real(r(a,a)); Eout = Eout + e(a)*real(r(a,a))
        end do
    end subroutine moments
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
end program test_carrier_carrier
