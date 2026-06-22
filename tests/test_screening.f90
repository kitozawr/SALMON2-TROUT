!
!  test_screening.f90  -  Part G screening primitives (sbe_superres_ssbe).
!
!  Unit tests of the dielectric screening models for the carrier-carrier
!  channel: Thomas-Fermi/Debye, static Lindhard/RPA (with the 2k_F kink), bulk
!  plasmon frequency, and the coupled LO-phonon-plasmon (LOPC) branches (Vieta
!  invariants + anticrossing). Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_screening
    use sbe_superres_ssbe
    implicit none
    integer :: nfail
    real(8) :: F, eps, wp2, wLp2, wLm2, wLO2, wTO2
    real(8), parameter :: TOL = 1d-9
    real(8), parameter :: PI = 3.14159265358979323846d0

    nfail = 0

    ! --- Thomas-Fermi eps(q) = 1 + kappa^2/q^2 -------------------------------
    call chk("eps_TF(q=kappa)=2", eps_thomas_fermi(2d0, 4d0), 2d0, TOL)   ! q^2=4, k2=4
    if (eps_thomas_fermi(1d-6, 1d0) < 1d6) call bad("eps_TF(q->0) not large")
    call chk("eps_TF(q->inf)->1", eps_thomas_fermi(1d6, 1d0), 1d0, 1d-6)

    ! --- TF / Debye kappa^2 limits -------------------------------------------
    call chk("kappa_TF^2 formula", tf_kappa2_degenerate(1d0, 1d0), &
             4d0*(3d0/PI)**(1d0/3d0), 1d-9)
    call chk("kappa_D^2 formula", debye_kappa2(2d0, 4d0, 0.5d0), &
             4d0*PI*2d0/(4d0*0.5d0), 1d-9)
    call chk("kappa_D^2(T=0 guard)=0", debye_kappa2(1d0, 1d0, 0d0), 0d0, TOL)

    ! --- Lindhard F(x): F(0)=1, F(1)=1/2, F(0.5) known value -----------------
    call chk("F(0)=1", lindhard_F(0d0), 1d0, TOL)
    call chk("F(1)=1/2", lindhard_F(1d0), 0.5d0, TOL)
    F = 0.5d0 + (1d0-0.25d0)/(4d0*0.5d0)*log(abs(1.5d0/0.5d0))
    call chk("F(0.5) closed form", lindhard_F(0.5d0), F, TOL)
    if (.not. (lindhard_F(0.5d0) > lindhard_F(1d0) .and. &
               lindhard_F(1d0) > lindhard_F(2d0))) call bad("F(x) not decreasing")

    ! --- static Lindhard reduces to Thomas-Fermi at small q (x->0) -----------
    eps = eps_lindhard_static(1d-3, 1d0, 4d0)
    call chk("eps_Lindhard(small q)->eps_TF", eps, eps_thomas_fermi(1d-3, 4d0), 1d-3)
    call chk("eps_Lindhard(no carriers)=eps_TF", &
             eps_lindhard_static(0.7d0, 0d0, 4d0), eps_thomas_fermi(0.7d0, 4d0), TOL)

    ! --- plasmon frequency wp^2 = 4 pi n/(eps_inf m*) ------------------------
    call chk("plasmon wp^2 formula", plasmon_freq2(1d0, 2d0, 0.5d0), &
             4d0*PI/(2d0*0.5d0), 1d-9)

    ! --- LOPC branches: Vieta invariants + anticrossing ----------------------
    wLO2 = 36d0**2; wTO2 = 33.6d0**2
    wp2  = 36d0**2                ! plasmon == LO -> anticrossing
    call lopc_branches(wp2, wLO2, wTO2, wLp2, wLm2)
    call chk("LOPC sum = wp^2+wLO^2", wLp2+wLm2, wp2+wLO2, 1d-6)
    call chk("LOPC product = wp^2 wTO^2", wLp2*wLm2, wp2*wTO2, 1d-3)
    if (.not. (wLp2 > wLO2 .and. wLm2 < wTO2)) call bad("LOPC branches not split around modes")

    if (nfail == 0) then
        write(*,'(a)') 'PASS'; call exit(0)
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
end program test_screening
