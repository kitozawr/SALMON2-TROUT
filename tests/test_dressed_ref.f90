!
!  test_dressed_ref.f90 - Option A dressed-reference carrier measure
!  (wiki/10 sec. 3A / 8.10). The ring channels measure carriers against the
!  field-rotated ground state instead of the static {occ, 0} reference.
!
!  Checks on a synthetic nba=3 rotation W(theta) (valence 1..nv mixed with
!  the conduction states):
!    1) sum_a delta0_a = 0 EXACTLY (unitarity => trace-neutral measure).
!    2) theta = 0 (no field): delta0 = 0 identically.
!    3) FROZEN (non-reacting) state f = f0[W]: f_eff = f - delta0 = {occ,0}
!       exactly -- the channels see equilibrium, zero fabrication.
!    4) ADIABATIC state f = {occ,0}: upper excess is NEGATIVE (clamped to 0
!       by the caller) -- also protected.
!    5) A real carrier added on top of the frozen background survives in
!       f_eff exactly.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_dressed_ref
    use sbe_superres_ssbe, only: dressed_ref_delta
    implicit none
    integer, parameter :: NBA = 3, NV = 1
    real(8), parameter :: OCC = 2.0d0, TH = 0.4d0
    complex(8) :: W(NBA, NBA)
    real(8) :: delta(NBA), f(NBA), feff(NBA)
    integer :: a, nfail

    nfail = 0

    ! --- (2) identity rotation => zero delta --------------------------------
    W = (0d0, 0d0)
    do a = 1, NBA
        W(a, a) = (1d0, 0d0)
    end do
    call dressed_ref_delta(NBA, NV, OCC, W, delta)
    if (maxval(abs(delta)) > 1d-14) call bad('delta0 /= 0 at zero field')

    ! --- rotation mixing valence(1) with conduction(2); band 3 untouched ----
    W = (0d0, 0d0)
    W(1, 1) = cmplx(cos(TH), 0d0, 8);  W(2, 1) = cmplx(sin(TH), 0d0, 8)
    W(1, 2) = cmplx(-sin(TH), 0d0, 8); W(2, 2) = cmplx(cos(TH), 0d0, 8)
    W(3, 3) = (1d0, 0d0)
    call dressed_ref_delta(NBA, NV, OCC, W, delta)

    ! --- (1) trace-neutral ---------------------------------------------------
    if (abs(sum(delta)) > 1d-14) call bad('sum delta0 /= 0 (measure not trace-neutral)')
    write(*,'(a,3es14.6)') '  delta0 = ', delta

    ! --- (3) frozen state: f = f0[W] => f_eff = {occ, 0} exactly ------------
    do a = 1, NBA
        f(a) = OCC * abs(W(1, a))**2       ! the rotated GS populations
    end do
    feff = f - delta
    if (abs(feff(1) - OCC) > 1d-14 .or. abs(feff(2)) > 1d-14 .or. abs(feff(3)) > 1d-14) &
        call bad('frozen (non-reacting) state must map to equilibrium exactly')
    write(*,'(a,3es14.6)') '  frozen f_eff = ', feff

    ! --- (4) adiabatic state: f = {occ,0} => upper excess negative ----------
    f = 0d0; f(1) = OCC
    feff = f - delta
    if (feff(2) > 1d-14) call bad('adiabatic state must not show positive upper excess')

    ! --- (5) real carrier on top of the frozen background survives ----------
    do a = 1, NBA
        f(a) = OCC * abs(W(1, a))**2
    end do
    f(3) = f(3) + 0.123d0                  ! genuine carrier in band 3
    feff = f - delta
    if (abs(feff(3) - 0.123d0) > 1d-14) call bad('real carrier must survive the reference subtraction')

    if (nfail == 0) then
        write(*,'(a)') 'PASS'
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'
        stop 1
    end if

contains
    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(2a)') '  FAILED: ', msg
        nfail = nfail + 1
    end subroutine bad
end program test_dressed_ref
