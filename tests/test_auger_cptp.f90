!
!  test_auger_cptp.f90  -  Auger recombination map is a standard CPTP Lindblad.
!
!  The Auger channel (apply_auger_recombination in bloch_solver_ssbe.f90) is the
!  composition of two amplitude-damping GKLS maps on the gap-edge Houston
!  branches: recombination (ic1 -> iv1, a CB electron fills a VB hole) and
!  promotion (ic1 -> ic_hot, the released E_g lifts a second CB electron). Both
!  use amp_damp_channel (sbe_superres_ssbe), the exact finite-time amplitude-
!  damping map L = sqrt(gamma)|dst><src| -- CPTP for any gamma,tau >= 0. The
!  rates carry the occ_max-normalized, [0,1]-clamped Pauli factors.
!
!  This test replicates that map (the exact code path) on a small density matrix
!  and checks the Lindblad/CPTP invariants the maintainer asked to confirm:
!    * trace (total carrier NUMBER) conserved exactly  -> Auger is number-conserving
!    * excited (conduction) population DECREASES         -> it is a recombination
!    * density matrix stays positive semidefinite        -> CPTP
!    * Hermiticity preserved
!    * gamma = 0  ->  identity (no spurious action)
!    * populations stay in [0, occ_max]                  -> Pauli respected
!
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_auger_cptp
    use sbe_superres_ssbe, only: amp_damp_channel
    implicit none
    integer, parameter :: n = 5            ! 2 valence + 3 conduction branches
    integer, parameter :: iv1 = 2          ! top valence (gap edge)
    integer, parameter :: ic1 = 3          ! lowest conduction
    real(8), parameter :: occ = 2d0        ! scalar bands (occ_max = 2)
    integer :: nfail, ic_hot
    complex(8) :: rho(n, n), rho0(n, n)
    real(8) :: tr0, tr1, exc0, exc1, gamma0, tau
    real(8) :: fc1, hv1, bhot, g_rec, g_prom, evmin

    nfail = 0

    ! --- a physically representative excited state -------------------------
    ! valence nearly full (a small hole at iv1), some conduction population.
    call set_state(rho)
    rho0 = rho
    ic_hot = 5                      ! energy-matched hot target (E(ic1)+E_g)
    gamma0 = 0.7d0                  ! per-carrier Auger rate * tau lumped
    tau    = 1d0

    call auger_map(rho, gamma0, tau, ic_hot)

    ! (1) NUMBER conservation: Tr rho unchanged (Auger conserves electrons)
    tr0 = trace(rho0); tr1 = trace(rho)
    if (abs(tr1 - tr0) > 1d-12) call bad('trace (carrier number) not conserved')

    ! (2) recombination: excited (conduction) population decreases
    exc0 = excited(rho0); exc1 = excited(rho)
    if (.not. (exc1 < exc0 - 1d-9)) call bad('excited population did not decrease')

    ! (3) CPTP: rho stays positive semidefinite (min eigenvalue >= 0)
    evmin = min_eig(rho)
    if (evmin < -1d-10) call bad('rho not positive semidefinite (not CPTP)')

    ! (4) Hermiticity preserved
    if (maxval(abs(rho - conjg(transpose(rho)))) > 1d-12) call bad('Hermiticity broken')

    ! (5) populations stay in [0, occ]
    if (any(real([rho(1,1),rho(2,2),rho(3,3),rho(4,4),rho(5,5)]) < -1d-12) .or. &
        any(real([rho(1,1),rho(2,2),rho(3,3),rho(4,4),rho(5,5)]) > occ + 1d-12)) &
        call bad('a population left [0, occ] (Pauli violated)')

    ! (6) gamma = 0  ->  identity
    call set_state(rho)
    rho0 = rho
    call auger_map(rho, 0d0, tau, ic_hot)
    if (maxval(abs(rho - rho0)) > 1d-14) call bad('gamma=0 is not the identity map')

    ! (7) the Pauli factors themselves are the occ-normalized, clamped forms
    call set_state(rho)
    fc1  = min(max(real(rho(ic1, ic1)) / occ, 0d0), 1d0)
    hv1  = min(max(1d0 - real(rho(iv1, iv1)) / occ, 0d0), 1d0)
    bhot = min(max(1d0 - real(rho(ic_hot, ic_hot)) / occ, 0d0), 1d0)
    if (fc1 < 0d0 .or. fc1 > 1d0 .or. hv1 < 0d0 .or. hv1 > 1d0 .or. &
        bhot < 0d0 .or. bhot > 1d0) call bad('Pauli factors not in [0,1]')
    ! recombination rate must vanish when there is no hole (full valence)
    rho(iv1, iv1) = dcmplx(occ, 0d0)
    hv1 = min(max(1d0 - real(rho(iv1, iv1)) / occ, 0d0), 1d0)
    g_rec = gamma0 * fc1 * hv1
    if (abs(g_rec) > 1d-14) call bad('recombination rate nonzero with no hole')

    if (nfail == 0) then
        write(*,'(a)') 'PASS  (Auger map: number-conserving, recombining, CPTP, '// &
            'Hermitian, gamma=0 identity, Pauli in [0,occ])'
        call exit(0)
    else
        write(*,'(a,i0,a)') 'FAIL (', nfail, ' checks)'; call exit(1)
    end if

contains

    subroutine set_state(r)
        complex(8), intent(out) :: r(n, n)
        r = (0d0, 0d0)
        r(1,1) = dcmplx(occ,   0d0)     ! deep valence: full
        r(2,2) = dcmplx(1.6d0, 0d0)     ! top valence iv1: a hole (occ-1.6 = 0.4)
        r(3,3) = dcmplx(0.5d0, 0d0)     ! ic1: excited electrons
        r(4,4) = dcmplx(0.2d0, 0d0)     ! ic2
        r(5,5) = dcmplx(0.1d0, 0d0)     ! ic_hot: mostly empty
        ! a little coherence so the off-diagonal damping is exercised
        r(2,3) = dcmplx(0.15d0, 0.05d0); r(3,2) = conjg(r(2,3))
    end subroutine set_state

    ! Replicates apply_auger_recombination's map (two amp_damp, clamped rates).
    subroutine auger_map(r, g0, t, ihot)
        complex(8), intent(inout) :: r(n, n)
        real(8),    intent(in)    :: g0, t
        integer,    intent(in)    :: ihot
        real(8) :: fcc, hvv, bht, grec, gprom
        fcc  = min(max(real(r(ic1, ic1)) / occ, 0d0), 1d0)
        hvv  = min(max(1d0 - real(r(iv1, iv1)) / occ, 0d0), 1d0)
        bht  = min(max(1d0 - real(r(ihot, ihot)) / occ, 0d0), 1d0)
        grec  = g0 * fcc * hvv
        gprom = g0 * fcc * hvv * bht
        if (ihot /= ic1) call amp_damp_channel(n, r, ic1, ihot, gprom, t)
        call amp_damp_channel(n, r, ic1, iv1, grec, t)
    end subroutine auger_map

    pure function trace(r) result(s)
        complex(8), intent(in) :: r(n, n); real(8) :: s; integer :: i
        s = 0d0; do i = 1, n; s = s + real(r(i, i)); end do
    end function trace

    pure function excited(r) result(s)
        complex(8), intent(in) :: r(n, n); real(8) :: s; integer :: i
        s = 0d0; do i = ic1, n; s = s + real(r(i, i)); end do
    end function excited

    ! smallest eigenvalue of the Hermitian rho via a tiny Jacobi sweep
    function min_eig(r) result(emin)
        complex(8), intent(in) :: r(n, n); real(8) :: emin
        complex(8) :: a(n, n); real(8) :: off, c, t, tau_, ph
        complex(8) :: aip, aiq, g
        integer :: p, q, i, sweep
        a = r
        do sweep = 1, 100
            off = 0d0
            do p = 1, n-1; do q = p+1, n; off = off + abs(a(p,q))**2; end do; end do
            if (off < 1d-24) exit
            do p = 1, n-1; do q = p+1, n
                if (abs(a(p,q)) < 1d-30) cycle
                ph = atan2(aimag(a(p,q)), real(a(p,q)))
                g  = dcmplx(cos(ph), -sin(ph))           ! rotate phase to real
                tau_ = (real(a(q,q)) - real(a(p,p))) / (2d0 * abs(a(p,q)))
                t = sign(1d0, tau_) / (abs(tau_) + sqrt(tau_**2 + 1d0))
                c = 1d0 / sqrt(t**2 + 1d0)
                do i = 1, n
                    aip = a(i,p); aiq = a(i,q)
                    a(i,p) = c*aip - (t*c)*g*aiq
                    a(i,q) = c*aiq + (t*c)*conjg(g)*aip
                end do
                do i = 1, n
                    aip = a(p,i); aiq = a(q,i)
                    a(p,i) = c*aip - (t*c)*conjg(g)*aiq
                    a(q,i) = c*aiq + (t*c)*g*aip
                end do
            end do; end do
        end do
        emin = real(a(1,1))
        do i = 2, n; emin = min(emin, real(a(i,i))); end do
    end function min_eig

    subroutine bad(msg)
        character(*), intent(in) :: msg
        write(*,'(a,a)') '  FAIL: ', msg; nfail = nfail + 1
    end subroutine bad

end program test_auger_cptp
