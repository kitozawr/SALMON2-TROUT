!
!  test_eph_cptp.f90  -  Part C5/C6 CPTP gate for the e-ph relaxation map.
!
!  Tests amp_damp_channel (the population-relaxing GKLS jump used by the
!  electron-phonon channel): exact transfer formulas, trace preservation,
!  coherence damping e^{-g/2}, Hermiticity, populations staying in [0,1], and
!  genuine complete positivity for the qubit case (det >= 0 after the map).
!  Standalone gfortran (uses src/ssbe/sbe_superres_ssbe.f90).
!
program test_eph_cptp
    use sbe_superres_ssbe, only: amp_damp_channel
    implicit none
    integer :: nfail, it
    complex(8) :: rho(3,3), rho2(2,2)
    real(8) :: g, gh, tr0, tr1, det, p1, p2, gamma, tau
    real(8), parameter :: TOL = 1d-12

    nfail = 0

    ! --- exact transfer formulas + trace + coherence damping (3-level) -------
    rho = (0d0, 0d0)
    rho(1,1) = (0.6d0, 0d0); rho(2,2) = (0.3d0, 0d0); rho(3,3) = (0.1d0, 0d0)
    rho(1,3) = (0.2d0, 0.1d0); rho(3,1) = conjg(rho(1,3))
    tr0 = real(rho(1,1)+rho(2,2)+rho(3,3))
    gamma = 3.0d0; tau = 0.25d0
    g = exp(-gamma*tau); gh = sqrt(g)
    call amp_damp_channel(3, rho, 1, 2, gamma, tau)   ! transfer 1 -> 2
    tr1 = real(rho(1,1)+rho(2,2)+rho(3,3))
    call chk("trace preserved", tr1, tr0, TOL)
    call chk("rho11 -> g*0.6", real(rho(1,1)), g*0.6d0, TOL)
    call chk("rho22 -> 0.3+(1-g)*0.6", real(rho(2,2)), 0.3d0+(1d0-g)*0.6d0, TOL)
    call chk("coherence rho13 damped e^{-g/2}", abs(rho(1,3)), gh*abs(dcmplx(0.2d0,0.1d0)), TOL)
    ! Hermiticity
    if (abs(rho(1,3)-conjg(rho(3,1))) > TOL) call bad("Hermiticity broken")

    ! --- populations stay in [0,1] under repeated random transfers -----------
    rho = (0d0,0d0); rho(1,1)=(0.5d0,0d0); rho(2,2)=(0.4d0,0d0); rho(3,3)=(0.1d0,0d0)
    do it = 1, 50
        call amp_damp_channel(3, rho, mod(it,3)+1, mod(it+1,3)+1, 1.0d0, 0.1d0)
    end do
    tr1 = 0d0
    do it = 1, 3
        p1 = real(rho(it,it)); tr1 = tr1 + p1
        if (p1 < -TOL .or. p1 > 1d0+TOL) call bad("population left [0,1]")
    end do
    call chk("trace preserved (repeated)", tr1, 1d0, 1d-10)

    ! --- complete positivity (qubit): det(rho') >= 0 -------------------------
    ! start from a valid 2x2 density matrix (PSD, trace 1) with strong coherence
    rho2(1,1)=(0.5d0,0d0); rho2(2,2)=(0.5d0,0d0)
    rho2(1,2)=(0.45d0,0d0); rho2(2,1)=(0.45d0,0d0)     ! det = 0.25-0.2025 = 0.0475 >= 0
    do it = 1, 20
        call amp_damp_channel(2, rho2, 1, 2, 0.7d0, 0.2d0)
        det = real(rho2(1,1))*real(rho2(2,2)) - abs(rho2(1,2))**2
        if (det < -1d-12) call bad("qubit det < 0 (positivity violated)")
        if (abs(real(rho2(1,1)+rho2(2,2)) - 1d0) > 1d-10) call bad("qubit trace drift")
    end do

    ! --- gamma=0 is the identity --------------------------------------------
    rho = (0d0,0d0); rho(1,1)=(0.7d0,0d0); rho(2,2)=(0.3d0,0d0); rho(1,2)=(0.1d0,0.2d0)
    rho(2,1)=conjg(rho(1,2))
    rho2(1,1)=rho(1,1)  ! stash
    call amp_damp_channel(3, rho, 1, 2, 0.0d0, 0.5d0)
    if (abs(real(rho(1,1))-0.7d0) > TOL) call bad("gamma=0 not identity")

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
end program test_eph_cptp
