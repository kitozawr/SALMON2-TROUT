!
!  test_colmem_lines.f90 - collisional-memory dephasing lines (wiki/10 sec. 8.6)
!  [MT99] Meier-Tannor JCP 111, 3365 (1999) Lorentzian lines built from the
!  cited phonon table; the Markovian gout damping is the zero-frequency anchor.
!
!  Checks:
!    1) Markov anchor: R(0) = Re sum c/mu = 1 EXACTLY after normalization.
!    2) Detailed balance: per mode, the emission/absorption line weights are
!       (N+1)/N (before the common normalization).
!    3) Markov limit: tau_c -> 0 makes the response flat, R(w) -> R(0) = 1
!       for |w| far beyond the mode energies.
!    4) Non-Markovian suppression: with a finite tau_c, modulation far above
!       the lines is barely damped, R(w >> w_p) << 1 -- the bath cannot
!       follow sub-correlation-time dynamics [B25 Fig 5(b)].
!    5) Resonance structure: R peaks near the phonon lines.
!    6) Toy propagation: z-updates (decay + source, the ring discretization)
!       damp a STATIC coherence at exactly the Markovian rate g (vs the
!       analytic exponential), and damp a fast-modulated coherence LESS.
!  Standalone gfortran (uses sbe_superres_ssbe.f90).
!
program test_colmem_lines
    use sbe_superres_ssbe, only: colmem_lines, colmem_response, bose_factor, mev_to_ha
    implicit none
    integer, parameter :: NPH = 3
    real(8), parameter :: KB_HA_K = 3.166811563d-6
    real(8) :: hw(NPH), wrel(NPH), nb(NPH), kT, tauc, r0, rfar, rres, ghost
    complex(8) :: cl(2*NPH), mul(2*NPH)
    integer :: nl, p, nfail
    ! toy propagation
    integer, parameter :: NT = 200000
    real(8), parameter :: DT = 0.5d0, GOUT = 2.0d-4
    complex(8), allocatable :: z(:)
    complex(8) :: rho, mem
    real(8) :: amp_static, amp_fast, amp_ref, wmod
    integer :: it, j

    nfail = 0
    kT = 300d0 * KB_HA_K
    hw(1) = mev_to_ha(19d0); hw(2) = mev_to_ha(51d0); hw(3) = mev_to_ha(63d0)
    wrel = (/ 0.2d0, 0.3d0, 0.5d0 /)
    do p = 1, NPH
        nb(p) = bose_factor(hw(p), kT)
    end do
    tauc = 1d0 / (0.2d0 / 27.211386245988d0)      ! 1/sigma_E, sigma_E = 0.2 eV

    call colmem_lines(NPH, hw, wrel, nb, tauc, nl, cl, mul)
    if (nl < NPH) call bad('missing lines')

    ! --- (1) anchor ---------------------------------------------------------
    r0 = colmem_response(nl, cl, mul, 0d0)
    write(*,'(a,i3,a,es14.6)') '  lines = ', nl, '   R(0) = ', r0
    if (abs(r0 - 1d0) > 1d-12) call bad('Markov anchor R(0) /= 1')

    ! --- (2) detailed balance (line pair ratio = (N+1)/N per mode) ----------
    ! lines come in (emission, absorption) pairs per mode in order
    if (abs(real(cl(1))/real(cl(2)) - (nb(1)+1d0)/nb(1)) > 1d-10) &
        call bad('emission/absorption weight ratio /= (N+1)/N')

    ! --- (3) Markov limit: tiny tau_c => flat response ----------------------
    call colmem_lines(NPH, hw, wrel, nb, 1d-2, nl, cl, mul)
    rfar = colmem_response(nl, cl, mul, 50d0 * hw(3))
    write(*,'(a,es14.6)') '  Markov-limit R(50*w_max) = ', rfar
    if (abs(rfar - 1d0) > 0.05d0) call bad('tau_c -> 0 must give a flat (Markovian) response')

    ! --- (4)+(5) finite memory: suppression far above the lines, peak near --
    call colmem_lines(NPH, hw, wrel, nb, tauc, nl, cl, mul)
    rfar = colmem_response(nl, cl, mul, 30d0 * hw(3))
    rres = colmem_response(nl, cl, mul, hw(3))
    ghost = colmem_response(nl, cl, mul, 3d0 * hw(3))
    write(*,'(a,3es14.6)') '  R(w3), R(3 w3), R(30 w3) = ', rres, ghost, rfar
    if (rfar > 0.05d0) call bad('fast modulation must be nearly undamped (memory)')
    if (rres < ghost) call bad('response must peak near the phonon lines')

    ! --- (6) toy propagation of the ring discretization ---------------------
    ! static coherence, damped by the memory term at rate GOUT*R(0) = GOUT
    allocate(z(nl))
    rho = (1d0, 0d0); z = (0d0, 0d0)
    do it = 1, NT
        mem = (0d0, 0d0)
        do j = 1, nl
            mem = mem + cl(j) * z(j)
        end do
        do j = 1, nl
            z(j) = z(j) * exp(-mul(j) * DT) + rho * DT
        end do
        rho = rho - GOUT * DT * mem
    end do
    ! compare decay RATES (a small rate error compounds over Gamma*t; and the
    ! coherence's own decay g shifts the response by O(g*tau_c) -- keep g << 1/tau_c)
    amp_static = -log(max(abs(rho), 1d-300)) / (NT * DT)
    amp_ref = GOUT
    write(*,'(a,2es14.6)') '  static decay rate: memory vs g = ', amp_static, amp_ref
    if (abs(amp_static - amp_ref) > 0.05d0 * amp_ref) &
        call bad('static coherence must decay at the Markovian rate (anchor)')

    ! fast-modulated coherence: apply the same updates to rho(t) = e^{-i w t}
    ! (modulation far above the lines) -- damping must be far weaker
    wmod = 20d0 * hw(3)
    rho = (1d0, 0d0); z = (0d0, 0d0)
    do it = 1, NT
        mem = (0d0, 0d0)
        do j = 1, nl
            mem = mem + cl(j) * z(j)
        end do
        do j = 1, nl
            z(j) = z(j) * exp(-mul(j) * DT) + rho * DT
        end do
        rho = (rho - GOUT * DT * mem) * exp(cmplx(0d0, -wmod * DT, 8))
    end do
    amp_fast = -log(max(abs(rho), 1d-300)) / (NT * DT)
    write(*,'(a,2es14.6)') '  fast-modulated decay rate vs g = ', amp_fast, GOUT
    if (amp_fast > 0.2d0 * GOUT) &
        call bad('fast-modulated coherence must be much less damped than Markov')

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
end program test_colmem_lines
