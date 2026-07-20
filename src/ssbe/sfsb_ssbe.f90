! =============================================================================
! SFSB (strong-field spin-boson) non-Markovian ionization mode.
! [B25] Boroumand, Thorpe, Bart, Parks, Toutounji, Vampa, Brabec, Wang,
! "Strong field physics in open quantum systems",
! Rep. Prog. Phys. 88, 070501 (2025) -- transcribed in wiki/10 sec. 6.
!
! Second-order (Dyson) conduction population per k and (v,c) pair,
!
!   nc(K,t) = (1/2) Re int^t dt1 int^t1 dt2 Om*(K_t1,t1) Om(K_t2,t2)
!                                     exp[ i S(t1,t2) + C(t1-t2) ]     (Eq. 3)
!
! with Om = 2 d(K_t).E(t), S = int Es/hbar, Es = sqrt(E^2+|hbar Om|^2) the
! Stark-shifted gap, K_t = K + A(t) (a.u.), and ALL environment influence in
! the bath correlation function C(t1-t2) (Eq. 5; sbe_superres bath_corr_table).
! The exp[C] factor does not factorize across time steps -- the inner integral
! is a true history (memory) sum, evaluated by sfsb_nc_series.
!
! Scope / validity (documented, honest):
!  - perturbative in the drive (2nd order): valid for small ionized fraction
!    (multiphoton, gamma > 1); NOT for tunneling depletion (gamma << 1).
!  - the k-line is 1D: requires num_kgrid = (N,1,1) and E || b1, so the
!    K + A(t) trajectory stays on the sampled line ([B25] uses the same 1D
!    Gamma-M reduction and verifies it against full 3D).
!  - two-band reduction per conduction band: gap = band-edge E_c - E_vtop
!    (continuous by eigenvalue sorting), coupling = bright-state quadrature
!    sum over the top-nv valence manifold (the gauge-invariant object; the
!    individual |d_vc| of degenerate members are gauge-random per k). Berry
!    connections / channel phases neglected -- the [B25] idealization.
!  - populations do not feed back (no depletion) -- channels sum independently.
!
! Outputs: SYSNAME_sfsb_nex.data (nc per cell volume vs t, same normalization
! as _sbe_nex.data) and SYSNAME_sfsb_nck.data (k-resolved blocks, _sbe_nex_k
! format) -- both plottable with the existing tooling.
! =============================================================================
module sfsb_ssbe
    implicit none
contains

subroutine main_sfsb_ssbe(icomm)
    use salmon_global
    use communication
    use gs_info_ssbe
    use em_field
    use util_ssbe, only: split_range
    use filesystem, only: get_filehandle
    use phys_constants, only: au_fs, au_ev, kB_au
    use sbe_superres_ssbe, only: bath_corr_table, bath_t2_high_t, sfsb_nc_series
    implicit none
    integer, intent(in) :: icomm

    real(8), parameter :: pi_c = 3.14159265358979323846d0
    type(s_sbe_gs_info) :: gs
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: nk, nb_vb, irank, nproc
    integer :: fh_nex, fh_nck
    real(8) :: occ_max, kT, wc, t2r, b1hat(3), a1(3), b1(3)
    integer :: it, ik, jk, ipair, npair, m
    integer :: stride, nts, nwin
    real(8) :: dts, esmax, dmax, emax_f, dq(3), t

    ! sorted 1D k-line
    integer, allocatable :: kord(:)          ! kord(i) = ik of i-th line point
    real(8), allocatable :: q1line(:)        ! sorted reduced q1
    real(8) :: h1, q10

    ! per-pair line data (gauge-smoothed real dipole + gap) and their periodic
    ! cubic-spline second derivatives. Piecewise-LINEAR sampling of de/d along
    ! the K_t trajectory puts slope kinks into the phase Theta(t) at every
    ! grid crossing; that broadband kink noise can exceed the true multiphoton
    ! signal by orders of magnitude (verified on the CdS line: 96 vs 192
    ! points changed nex 6x). The C2 spline removes it.
    real(8), allocatable :: de_line(:, :), dre_line(:, :)   ! (nline, npair)
    real(8), allocatable :: de_m2(:, :), dre_m2(:, :)       ! spline d2f/dx2
    integer, allocatable :: pc(:)

    ! per-step sampled series
    real(8), allocatable :: efld(:), aq1(:)                 ! strided E(t).b1hat, dq1(t)
    complex(8), allocatable :: ctab(:)
    real(8), allocatable :: nck(:, :), nck_l(:, :)
    real(8), allocatable :: nex_t(:), nex_l(:)
    integer, allocatable :: ik_lo(:), ik_hi(:)
    real(8) :: wsum

    call comm_get_groupinfo(icomm, irank, nproc)

    ! ------------------------------------------------------------------ GS --
    nk = num_kgrid(1) * num_kgrid(2) * num_kgrid(3)
    call init_sbe_gs_info(gs, sysname, base_directory, &
        & nk, nstate, nelec, al_vec1, al_vec2, al_vec3, .false., icomm)

    if (num_kgrid(2) /= 1 .or. num_kgrid(3) /= 1) then
        if (irank == 0) write(*, '(a)') &
            '# ERROR: yn_sbe_sfsb requires a 1D k-line: num_kgrid = (N,1,1).'
        error stop 'sfsb: num_kgrid must be (N,1,1)'
    end if

    select case (trim(sbe_bath_model))
    case ('none', 'rta', 'ohmic', 'debye')
    case default
        error stop 'sfsb: sbe_bath_model must be none|ohmic|debye|rta'
    end select
    wc = sbe_bath_wc_ev / au_ev
    if ((trim(sbe_bath_model) == 'ohmic' .or. trim(sbe_bath_model) == 'debye') &
        .and. wc <= 0d0) then
        error stop 'sfsb: sbe_bath_wc_ev must be > 0 for the ohmic/debye bath'
    end if
    kT = max(sbe_bath_temperature_k, 0d0) * kB_au

    if (yn_sbe_spinor == 'y') then
        nb_vb = nelec
        occ_max = 1d0
    else
        nb_vb = nelec / 2
        occ_max = 2d0
    end if
    ! One line per conduction band; the valence partner is selected POINTWISE
    ! as the BRIGHT member of the top-sbe_sfsb_nv near-degenerate manifold
    ! (argmax |d.e_hat|). Energy-sorted band indices swap character at
    ! (avoided) crossings inside a degenerate manifold -- a fixed (v,c) pair
    ! line is then a bright/dark STEP function whose interpolation noise can
    ! exceed the true multiphoton signal by many orders (verified on the CdS
    ! Gamma-M line: v7/v8 swap, |d| steps 2.3 -> 1e-9; nex was grid-divergent).
    ! Bright-line tracking follows the allowed transition smoothly; the dark
    ! partner's own channel (|d| ~ 1e-9) is negligible and dropped. [B25]'s
    ! model is exactly this two-band (one bright valence) reduction.
    if (sbe_sfsb_nv < 1 .or. sbe_sfsb_nv > nb_vb) error stop 'sfsb: bad sbe_sfsb_nv'
    if (sbe_sfsb_nc < 1 .or. nb_vb + sbe_sfsb_nc > gs%nb) error stop 'sfsb: bad sbe_sfsb_nc'
    npair = sbe_sfsb_nc
    allocate(pc(npair))
    do ipair = 1, npair
        pc(ipair) = nb_vb + ipair
    end do

    ! ------------------------------------------------------------- field ----
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, -1, nt+1, Ac_ext_t)

    a1(1:3) = gs%a_matrix(1:3, 1)
    b1(1:3) = gs%b_matrix(1, 1:3)
    b1hat = b1 / sqrt(sum(b1**2))

    ! K_t = K + A(t): reduced displacement dq_i = A . a_i / (2 pi). The 1D
    ! line samples only q1, so the field must not move q2/q3: A || b1
    ! (b1 is orthogonal to a2 and a3 by construction).
    dq = 0d0
    do it = 0, nt
        dq(1) = max(dq(1), abs(dot_product(Ac_ext_t(1:3, it), gs%a_matrix(1:3, 1))))
        dq(2) = max(dq(2), abs(dot_product(Ac_ext_t(1:3, it), gs%a_matrix(1:3, 2))))
        dq(3) = max(dq(3), abs(dot_product(Ac_ext_t(1:3, it), gs%a_matrix(1:3, 3))))
    end do
    dq = dq / (2d0 * pi_c)
    ! tolerance 1e-4 of the full reduced BZ: al_vec inputs carry ~6 digits, so
    ! an algebraically-exact E || b1 leaves a roundoff-level transverse drift;
    ! 1e-4 of the BZ is far below the k-line resolution and band variation.
    if (dq(2) > 1d-4 .or. dq(3) > 1d-4) then
        if (irank == 0) write(*, '(a,3es12.3)') &
            '# ERROR: sfsb 1D mode requires E || b1; reduced |dq| per axis = ', dq
        error stop 'sfsb: field polarization must be along b1 (epdir || b1)'
    end if

    ! --------------------------------------------------- 1D line ordering ---
    allocate(kord(nk), q1line(nk))
    do ik = 1, nk
        kord(ik) = ik
    end do
    ! insertion sort by reduced q1 (nk is small)
    do ik = 2, nk
        m = kord(ik)
        jk = ik - 1
        do while (jk >= 1)
            if (gs%kpoint(1, kord(jk)) <= gs%kpoint(1, m)) exit
            kord(jk + 1) = kord(jk)
            jk = jk - 1
        end do
        kord(jk + 1) = m
    end do
    do ik = 1, nk
        q1line(ik) = gs%kpoint(1, kord(ik))
    end do
    h1 = 1d0 / dble(nk)     ! MP line spacing in reduced coords
    q10 = q1line(1)
    do ik = 2, nk
        if (abs((q1line(ik) - q1line(ik-1)) - h1) > 1d-6) &
            error stop 'sfsb: k-line is not a uniform Monkhorst-Pack 1D grid'
    end do

    ! ------------------------- pair line data (gap + smoothed dipole) -------
    allocate(de_line(nk, npair), dre_line(nk, npair))
    allocate(de_m2(nk, npair), dre_m2(nk, npair))
    call build_pair_lines()
    do ipair = 1, npair
        call spline_periodic_m2(nk, de_line(:, ipair), de_m2(:, ipair))
        call spline_periodic_m2(nk, dre_line(:, ipair), dre_m2(:, ipair))
    end do

    ! ------------------------------------------------------- stride/time ----
    ! resolve the fastest phase: Stark-shifted gap at the strongest field
    emax_f = 0d0
    do it = 1, nt - 1
        emax_f = max(emax_f, sqrt(sum(((Ac_ext_t(1:3, it+1) - Ac_ext_t(1:3, it-1)) / (2d0 * dt))**2)))
    end do
    dmax = maxval(abs(dre_line))
    esmax = sqrt(maxval(de_line)**2 + (2d0 * dmax * emax_f)**2)
    stride = sbe_sfsb_stride
    if (stride <= 0) stride = max(1, int((2d0 * pi_c / (24d0 * max(esmax, 1d-6))) / dt))
    nts = nt / stride
    dts = stride * dt

    nwin = 0
    if (sbe_bath_memory_fs > 0d0) nwin = max(1, int((sbe_bath_memory_fs / au_fs) / dts))

    ! --------------------------------------------------------- bath table ---
    allocate(ctab(0:nts))
    t2r = -1d0
    if (sbe_bath_rta_t2_fs > 0d0) t2r = sbe_bath_rta_t2_fs / au_fs
    if (trim(sbe_bath_model) == 'none') then
        ctab(:) = (0d0, 0d0)
    else
        call bath_corr_table(nts, dts, kT, sbe_bath_jo, wc, trim(sbe_bath_model), ctab, t2r)
    end if
    if (yn_sbe_bath_imc == 'n') then
        do m = 0, nts
            ctab(m) = cmplx(real(ctab(m)), 0d0, 8)
        end do
    end if

    if (irank == 0) then
        write(*, '(a)')        '# SFSB non-Markovian ionization mode [B25 = RPP 88, 070501 (2025)]'
        write(*, '(a,a)')      '#   bath model        = ', trim(sbe_bath_model)
        write(*, '(a,f10.4)')  '#   jo                = ', sbe_bath_jo
        write(*, '(a,f10.4)')  '#   hbar*wc [eV]      = ', sbe_bath_wc_ev
        write(*, '(a,f10.1)')  '#   T [K]             = ', sbe_bath_temperature_k
        if (trim(sbe_bath_model) == 'rta') then
            if (t2r > 0d0) then
                write(*, '(a,f10.3)') '#   T2 [fs] (input)   = ', t2r * au_fs
            else
                write(*, '(a,f10.3)') '#   T2 [fs] (derived) = ', bath_t2_high_t(kT, sbe_bath_jo) * au_fs
            end if
        end if
        write(*, '(a,i8,a,f8.3,a)') '#   stride            = ', stride, '  (dt_sfsb = ', dts, ' a.u.t)'
        write(*, '(a,i8)')     '#   time samples      = ', nts
        if (nwin > 0) then
            write(*, '(a,i8)') '#   memory window     = ', nwin
        else
            write(*, '(a)')    '#   memory window     = full history'
        end if
        write(*, '(a,i4,a,i4,a,i4)') '#   pairs (v x c)     = ', sbe_sfsb_nv, &
            ' x', sbe_sfsb_nc, '  bands up to', nb_vb + sbe_sfsb_nc
        write(*, '(a,es12.4)') '#   |exp(C)| at t_end = ', abs(exp(ctab(nts)))
    end if

    ! ------------------------------------------------ strided field series --
    allocate(efld(0:nts), aq1(0:nts))
    do it = 0, nts
        m = it * stride
        efld(it) = dot_product(b1hat, -(Ac_ext_t(1:3, m+1) - Ac_ext_t(1:3, m-1)) / (2d0 * dt))
        aq1(it)  = dot_product(Ac_ext_t(1:3, m), a1) / (2d0 * pi_c)
    end do

    ! ------------------------------------------------------------ main loop -
    ! MPI over the line points (split_range) x OpenMP inside each rank; each
    ! (jk, pair) writes only its own nck_l row, so the jk loop is conflict-free.
    allocate(nck(nk, 0:nts), nck_l(nk, 0:nts), nex_t(0:nts), nex_l(0:nts))
    allocate(ik_lo(nproc), ik_hi(nproc))
    call split_range(1, nk, nproc, ik_lo, ik_hi)
    nck_l = 0d0
    !$omp parallel do schedule(dynamic)
    do jk = ik_lo(irank + 1), ik_hi(irank + 1)     ! jk indexes the sorted line
        call sfsb_line_point(jk)
    end do
    !$omp end parallel do
    call comm_summation(nck_l, nck, nk * (nts + 1), icomm)

    ! nex(t) = kweight-average per cell / volume (same convention as _sbe_nex)
    wsum = sum(gs%kweight)
    do it = 0, nts
        nex_l(it) = 0d0
        do ik = 1, nk
            nex_l(it) = nex_l(it) + gs%kweight(ik) * nck(ik, it)
        end do
        nex_t(it) = nex_l(it) / (wsum * gs%volume)
    end do

    ! --------------------------------------------------------------- output -
    if (irank == 0) then
        fh_nex = get_filehandle()
        open(unit=fh_nex, file=trim(base_directory)//trim(sysname)//"_sfsb_nex.data", action="write")
        write(fh_nex, '(a)') '# SFSB memory-integral ionization [B25 RPP 88, 070501 (2025)]'
        write(fh_nex, '(a)') '# 1:Time[a.u.] 2:nex[a.u.] (conduction population per cell volume)'
        do it = 0, nts
            write(fh_nex, '(2es25.15e3)') it * dts, nex_t(it)
        end do
        close(fh_nex)

        fh_nck = get_filehandle()
        open(unit=fh_nck, file=trim(base_directory)//trim(sysname)//"_sfsb_nck.data", action="write")
        write(fh_nck, '(a)') '# SFSB k-resolved conduction population nc(K,t) [B25 Eq. (3)]'
        write(fh_nck, '(a)') '# blocks: t[a.u.]; rows: ik kx ky kz nc'
        do it = 0, nts
            t = it * dts
            if (it == nts .or. mod(it * stride, max(out_projection_k_step, 1)) < stride) then
                write(fh_nck, '(a,es25.15e3)') '# time = ', t
                do jk = 1, nk
                    ik = kord(jk)
                    write(fh_nck, '(i8,4es25.15e3)') ik, gs%kpoint(1:3, ik), nck(ik, it)
                end do
            end if
        end do
        close(fh_nck)

        write(*, '(a,es14.6)') '# SFSB done. nex(t_end) [a.u./cell volume] = ', nex_t(nts)
        write(*, '(a,es14.6,a)') '#            nex(t_end) = ', &
            nex_t(nts) * 1d24 / (0.529177210903d0**3), ' cm^-3'
    end if

    call comm_sync_all(icomm)

contains

    ! One k-line point: sample Om(K_t,t)/Es(K_t,t) along the A(t) trajectory
    ! (periodic linear interpolation on the sorted line) for every (v,c) pair
    ! and run the memory-integral stepper. All locals are per-call (thread-safe);
    ! the only host write is this jk's own nck_l row.
    subroutine sfsb_line_point(jj)
        integer, intent(in) :: jj
        complex(8) :: oml(0:nts)
        real(8) :: esl(0:nts), ncl(0:nts)
        integer :: p, itt, ii0, ii1
        real(8) :: xx, ff, dee, ddr
        real(8) :: c0, c1
        do p = 1, npair
            do itt = 0, nts
                ! K_t position on the periodic line (fractional index)
                xx = (q1line(jj) + aq1(itt) - q10) / h1
                xx = xx - dble(nk) * floor(xx / dble(nk))
                ii0 = int(xx)
                ff  = xx - ii0
                ii0 = 1 + mod(ii0, nk)
                ii1 = 1 + mod(ii0, nk)
                ! periodic cubic spline (in the fractional coordinate, h = 1)
                c0 = ((1d0 - ff)**3 - (1d0 - ff)) / 6d0
                c1 = (ff**3 - ff) / 6d0
                dee = (1d0 - ff) * de_line(ii0, p) + ff * de_line(ii1, p) &
                    + c0 * de_m2(ii0, p) + c1 * de_m2(ii1, p)
                ddr = (1d0 - ff) * dre_line(ii0, p) + ff * dre_line(ii1, p) &
                    + c0 * dre_m2(ii0, p) + c1 * dre_m2(ii1, p)
                oml(itt) = cmplx(2d0 * ddr * efld(itt), 0d0, 8)   ! Om = 2 d.E [B25]
                esl(itt) = sqrt(dee**2 + abs(oml(itt))**2)         ! Stark-shifted gap
            end do
            call sfsb_nc_series(nts, dts, oml, esl, ctab, nwin, ncl)
            nck_l(kord(jj), :) = nck_l(kord(jj), :) + occ_max * ncl(:)
        end do
    end subroutine sfsb_line_point

    ! Two-band-reduction lines per conduction band p [B25's E(K), d(K)]:
    !   gap      de(K) = E_c(K) - E_vtop(K) -- the band-edge gap. Sorted
    !            eigenvalues are CONTINUOUS in k (only the character swaps at
    !            crossings), so this line has at most slope kinks, never jumps.
    !   coupling d(K)  = sqrt( sum_{v in top-nv manifold} |d_vc . e_hat|^2 )
    !            -- the bright-state coupling of the valence manifold. The
    !            individual |d_vc| of (near-)degenerate members are NOT gauge
    !            invariant (the GS solver returns arbitrary per-k mixtures);
    !            only this quadrature sum is, so only it is smooth in k.
    ! Berry connections / residual channel phases are neglected -- the [B25]
    ! two-band idealization (their smooth parameterized E(K), d(K)).
    ! A selection-rule WALL (channel exactly forbidden over part of the line,
    ! e.g. CdS Gamma-M with E perp c: allowed pocket |q1| < 0.1 then S = 0)
    ! makes any two-band reduction non-smooth; detect and warn.
    subroutine build_pair_lines()
        integer :: jjk, iik, p, vv
        complex(8) :: dproj
        real(8) :: s2, step_max, rng
        do p = 1, npair
            do jjk = 1, nk
                iik = kord(jjk)
                de_line(jjk, p) = gs%eigen(pc(p), iik) - gs%eigen(nb_vb, iik)
                s2 = 0d0
                do vv = nb_vb - sbe_sfsb_nv + 1, nb_vb
                    dproj = gs%d_matrix(pc(p), vv, 1, iik) * b1hat(1) &
                          + gs%d_matrix(pc(p), vv, 2, iik) * b1hat(2) &
                          + gs%d_matrix(pc(p), vv, 3, iik) * b1hat(3)
                    s2 = s2 + abs(dproj)**2
                end do
                dre_line(jjk, p) = sqrt(s2)
            end do
            ! wall detector: an adjacent-point coupling jump comparable to the
            ! full range means the channel is not a smooth two-band line here
            step_max = 0d0
            do jjk = 1, nk
                step_max = max(step_max, abs(dre_line(1 + mod(jjk, nk), p) - dre_line(jjk, p)))
            end do
            rng = maxval(dre_line(:, p)) - minval(dre_line(:, p))
            if (irank == 0 .and. step_max > 0.5d0 * max(rng, 1d-12)) then
                write(*, '(a,i0,a)') '# WARNING: sfsb channel c=', pc(p), &
                    ': the coupling line jumps by > 50% of its range between'
                write(*, '(a)') '#   adjacent k-points -- a selection-rule wall or an unresolved'
                write(*, '(a)') '#   character swap. The two-band reduction is NOT smooth on this'
                write(*, '(a)') '#   line/polarization; refine the k-line or change polarization.'
            end if
        end do
    end subroutine build_pair_lines

end subroutine main_sfsb_ssbe

! Periodic (cyclic) cubic-spline second derivatives m for uniform unit
! spacing: m(i-1) + 4 m(i) + m(i+1) = 6 (f(i-1) - 2 f(i) + f(i+1)),
! indices mod n. Cyclic tridiagonal via Sherman-Morrison on two Thomas solves.
subroutine spline_periodic_m2(n, f, m)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: f(n)
    real(8), intent(out) :: m(n)
    real(8), allocatable :: d(:), y(:), z(:), u(:)
    real(8) :: gamma, fact
    integer :: i

    if (n < 3) then
        m = 0d0
        return
    end if
    allocate(d(n), y(n), z(n), u(n))
    do i = 1, n
        d(i) = 6d0 * (f(1 + mod(i - 2 + n, n)) - 2d0 * f(i) + f(1 + mod(i, n)))
    end do
    gamma = -4d0
    u = 0d0
    u(1) = gamma
    u(n) = 1d0
    call thomas_mod(n, gamma, d, y)
    call thomas_mod(n, gamma, u, z)
    fact = (y(1) + y(n) / gamma) / (1d0 + z(1) + z(n) / gamma)
    m = y - fact * z
    deallocate(d, y, z, u)
end subroutine spline_periodic_m2

! Thomas solve of the MODIFIED tridiagonal (1, 4, 1) system with
! b(1) = 4 - gamma and b(n) = 4 - 1/gamma (the cyclic correction rows).
subroutine thomas_mod(n, gamma, d, x)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: gamma, d(n)
    real(8), intent(out) :: x(n)
    real(8), allocatable :: b(:), cp(:)
    integer :: i

    allocate(b(n), cp(n))
    b(1) = 4d0 - gamma
    do i = 2, n - 1
        b(i) = 4d0
    end do
    b(n) = 4d0 - 1d0 / gamma
    x(1) = d(1) / b(1)
    cp(1) = 1d0 / b(1)
    do i = 2, n
        b(i) = b(i) - cp(i - 1)
        cp(i) = 1d0 / b(i)
        x(i) = (d(i) - x(i - 1)) / b(i)
    end do
    do i = n - 1, 1, -1
        x(i) = x(i) - cp(i) * x(i + 1)
    end do
    deallocate(b, cp)
end subroutine thomas_mod

end module sfsb_ssbe
