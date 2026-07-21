module realtime_ssbe
    implicit none
contains

subroutine main_realtime_ssbe(icomm)
    use salmon_global
    use communication
    use gs_info_ssbe
    use bloch_solver_ssbe
    use em_field
    use datafile_ssbe
    use input_checker_sbe
    use filesystem, only: get_filehandle
    use phys_constants, only: au_fs
    use sbe_superres_ssbe, only: vg_ptop_exceeds
    use iso_fortran_env, only: error_unit
    implicit none
    integer, intent(in) :: icomm

    type(s_sbe_bloch_solver) :: sbe
    type(s_sbe_gs_info) :: gs
    real(8) :: t, E(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it
    real(8) :: energy, tr_all, tr_vb
    integer :: nproc, irank, ierr
    integer :: fh_sbe_rt, fh_sbe_rt_energy, fh_sbe_nex, fh_sbe_nex_k
    integer :: fh_sbe_channels, it0, ck_unit
    logical :: ck_exists
    character(256) :: ck_file
    integer :: fh_sbe_nex_k_unfold
    integer :: fh_sbe_nex_k_real, fh_sbe_nex_k_unfold_real
    integer :: fh_sbe_nex_k_lev_real
    integer :: fh_sbe_intra_current
    real(8) :: jmat_intra(3)
    integer :: nk
    integer :: ib_lcb, nb_vb
    integer :: ib_top
    real(8) :: ptop
    real(8), parameter :: PTOP_TOL = 1.0d-3   ! VG basis-edge occupation tolerance
    real(8), allocatable :: pop_k(:)
    real(8), allocatable :: pop_top_k(:)
    real(8), allocatable :: pop_lev_k(:, :, :)
    real(8), allocatable :: pop_k_real(:)
    real(8), allocatable :: pop_lev_k_real(:, :, :)
    real(8), allocatable :: pop4(:, :)            ! 4 gap-edge diabatic levels (primitive)

    call comm_get_groupinfo(icomm, irank, nproc)

    if (.not. check_input_variables_sbe()) return

    ! Read ground state electronic system:
    nk = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
    call init_sbe_gs_info(gs, sysname, base_directory, &
        & nk, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., icomm)        
    
    ! Initialization of SBE solver and density matrix:
    call init_sbe_bloch_solver(sbe, gs, nstate_sbe(1), icomm)
    sbe%flag_vnl_correction = (yn_vnl_correction == 'y')

    ! Prepare external pulse (allocated from -1 to nt+1 for safe central differences).
    ! Pass the FULL index range: calc_Ac_ext_t's dummy is explicit-shape (is:ie),
    ! so passing the whole (-1:nt+1) array with is=0 sequence-associates it one
    ! slot early -- the field was shifted by one dt and slot nt+1 was left
    ! uninitialized (a dA/dt spike in the last step). i=-1 (t=-dt) fills as 0.
    allocate(Ac_ext_t(1:3, -1:nt+1))
    Ac_ext_t(:, :) = 0.0d0
    ! MPI-safe external-field construction. Only the root rank builds the pulse
    ! (which, for ae_shape='input', OPENS AND READS file_input1); the whole
    ! trajectory is then broadcast to every rank. This removes the N-way race in
    ! which all ranks simultaneously read the same text file from a shared
    ! filesystem: on a distributed machine a lagging or partially-staged file
    ! hands some ranks a short/empty read (n_dat small -> out-of-bounds in the
    ! interpolator) while others succeed, so the driven field silently diverges
    ! across ranks or one rank segfaults. Reading once and broadcasting makes the
    ! field bit-identical everywhere and confines any file error to a single,
    ! diagnosable rank. The barrier pins all ranks before the read so a
    ! non-synchronized start cannot race ahead of the root's open.
    call comm_sync_all(icomm)
    if (irank == 0) call calc_Ac_ext_t(0.0d0, dt, -1, nt+1, Ac_ext_t)
    call comm_bcast(Ac_ext_t, icomm, 0)

    ! Field-scale audit line (one-time): the peak |A|, the peak |E| = -dA/dt,
    ! and the peak crystal-momentum excursion. file_input1 pulses MUST be in the
    ! run's unit_system (A_eV_fs: A in fs*V/Angstrom, 1 a.u. = 1.2442 fs*V/A) --
    ! a file generated in the wrong units drives a silently wrong field; this
    ! line makes the actually-driven field visible for any pulse source.
    if (comm_is_root(irank)) then
        block
            real(8) :: amx, emx, an, en
            integer :: itf
            amx = 0d0; emx = 0d0
            do itf = 0, nt
                an = sqrt(sum(Ac_ext_t(:, itf)**2))
                en = sqrt(sum(((Ac_ext_t(:, itf+1) - Ac_ext_t(:, itf-1)) / (2d0*dt))**2))
                if (an > amx) amx = an
                if (en > emx) emx = en
            end do
            write(*, '(a,es11.4,a,es11.4,a,f10.3,a)') &
                '# field audit: peak |A| = ', amx, ' a.u. (= peak dk excursion), '// &
                'peak |E| = ', emx, ' a.u. = ', emx * 5.14220675d3, ' MV/cm'
        end block
    end if

    ! Initial energy and fields
    energy = 0.0d0
    E(:) = 0.0d0
    Jmat(:) = 0.0d0

    ! Number of valence bands: nelec spinor bands (occupation 1 each) when the
    ! GS input comes from a spin-orbit split system, nelec/2 scalar bands otherwise.
    if (yn_sbe_spinor == 'y') then
        nb_vb = nelec
    else
        nb_vb = nelec / 2
    end if

    ! Lowest conduction band index (Houston-basis population output)
    ib_lcb = nb_vb + 1
    allocate(pop_k(1:nk), pop_k_real(1:nk))
    if (gs%have_unfold) allocate(pop_lev_k(1:4, 1:4, 1:nk), pop_lev_k_real(1:4, 1:4, 1:nk))
    if (.not. gs%have_unfold) allocate(pop4(1:4, 1:nk))   ! primitive 4-level output

    ! Highest band carried into the dynamics (top of the active subspace) -- the
    ! VG basis edge. We monitor its peak adiabatic occupation P_top as the cheap
    ! necessary condition for N_b sufficiency (wiki: "VG Basis Sufficiency &
    ! N_b Convergence", criterion (a)): population reaching the basis edge means
    ! bands above the cutoff would have mattered too -> enlarge N_b.
    ib_top = 0
    if (sbe%n_active_bands > 0) ib_top = sbe%active_idx(sbe%n_active_bands)
    allocate(pop_top_k(1:nk))

    if (irank == 0) then
        ! SYSNAME_sbe_rt.data
        fh_sbe_rt = get_filehandle()
        open(unit=fh_sbe_rt, file=trim(base_directory)//trim(sysname)//"_sbe_rt.data", action="write")
        call write_sbe_rt_header(fh_sbe_rt)
        ! SYSNAME_sbe_rt_energy.data
        fh_sbe_rt_energy = get_filehandle()
        open(unit=fh_sbe_rt_energy, file=trim(base_directory)//trim(sysname)//"_sbe_rt_energy.data", action="write")
        call write_sbe_rt_energy_header(fh_sbe_rt_energy)
        ! SYSNAME_sbe_nex.data
        fh_sbe_nex = get_filehandle()
        open(unit=fh_sbe_nex, file=trim(base_directory)//trim(sysname)//"_sbe_nex.data", action="write")
        call write_sbe_nex_header(fh_sbe_nex)
        ! SYSNAME_sbe_nex_k.data (instantaneous Houston-basis LCB population)
        fh_sbe_nex_k = get_filehandle()
        open(unit=fh_sbe_nex_k, file=trim(base_directory)//trim(sysname)//"_sbe_nex_k.data", action="write")
        call write_sbe_nex_k_header(fh_sbe_nex_k, nk)
        ! SYSNAME_sbe_nex_k_real.data: REAL carriers only (fixed-basis diabatic
        ! LCB occupation, k-resolved n_ex) -- no reversible A^2(t) virtual breathing
        fh_sbe_nex_k_real = get_filehandle()
        open(unit=fh_sbe_nex_k_real, &
            & file=trim(base_directory)//trim(sysname)//"_sbe_nex_k_real.data", action="write")
        call write_sbe_nex_k_real_header(fh_sbe_nex_k_real, nk)
        ! C1: per-channel dissipation ledger (ring channels; cumulative)
        fh_sbe_channels = get_filehandle()
        open(unit=fh_sbe_channels, file=trim(base_directory)//trim(sysname)//"_sbe_channels.data", action="write")
        write(fh_sbe_channels, '(a)') '# C1 per-channel ledger (CUMULATIVE, per cell): ring channels only'
        write(fh_sbe_channels, '(a)') '# dN = conduction-population change (pairs created > 0), dE [Ha]'
        write(fh_sbe_channels, '(a)') '# t[au]  dN_eph dE_eph  dN_ii dE_ii  dN_auger dE_auger  dN_rana dE_rana'
        ! SYSNAME_sbe_nex_k_unfold.data: populations of PHYSICAL primitive
        ! bands at the unfolded primitive k-points (only with an unfold map)
        if (gs%have_unfold) then
            fh_sbe_nex_k_unfold = get_filehandle()
            open(unit=fh_sbe_nex_k_unfold, &
                & file=trim(base_directory)//trim(sysname)//"_sbe_nex_k_unfold.data", action="write")
            call write_sbe_nex_k_unfold_header(fh_sbe_nex_k_unfold, nk)
            ! REAL-carrier unfolded twin
            fh_sbe_nex_k_unfold_real = get_filehandle()
            open(unit=fh_sbe_nex_k_unfold_real, &
                & file=trim(base_directory)//trim(sysname)//"_sbe_nex_k_unfold_real.data", action="write")
            call write_sbe_nex_k_unfold_header(fh_sbe_nex_k_unfold_real, nk)
        else
            ! Primitive cell (no unfold map): four gap-edge diabatic populations
            ! (VB-1, VB, CB1, CB2) per k, so --spectral can colour all four bands.
            fh_sbe_nex_k_lev_real = get_filehandle()
            open(unit=fh_sbe_nex_k_lev_real, &
                & file=trim(base_directory)//trim(sysname)//"_sbe_nex_k_lev_real.data", action="write")
            call write_sbe_nex_k_lev_header(fh_sbe_nex_k_lev_real, nk)
        end if
        ! SYSNAME_sbe_intra_current.data: intra-band (Houston) current
        if (yn_out_intraband_current == 'y') then
            fh_sbe_intra_current = get_filehandle()
            open(unit=fh_sbe_intra_current, &
                & file=trim(base_directory)//trim(sysname)//"_sbe_intra_current.data", action="write")
            call write_sbe_intra_current_header(fh_sbe_intra_current)
        end if
        ! Stdout logs:
        write(*, "(a)") " time-step time[fs] Current(xyz)[a.u.]                     electrons   Total energy[au]"
        write(*, "(a)") "---------------------------------------------------------------------------------------"
    end if

    ! Write initial (t=0) nex_k block before the time loop. The LCB/CB1 outputs
    ! are zero at equilibrium, but the 4-level file must carry the REAL diabatic
    ! occupations (valence FULL, conduction empty) so the carrier (hole/electron)
    ! colour scale is not poisoned by a spurious t=0 "hole" -- compute them
    ! (collective) before the irank-0 write.
    if (.not. gs%have_unfold) then
        call calc_diabatic_population_k(sbe, ib_lcb-2, pop_k_real, icomm); pop4(1,:) = pop_k_real
        call calc_diabatic_population_k(sbe, ib_lcb-1, pop_k_real, icomm); pop4(2,:) = pop_k_real
        call calc_diabatic_population_k(sbe, ib_lcb,   pop_k_real, icomm); pop4(3,:) = pop_k_real
        call calc_diabatic_population_k(sbe, ib_lcb+1, pop_k_real, icomm); pop4(4,:) = pop_k_real
    end if
    if (irank == 0) then
        pop_k = 0.0d0
        call write_sbe_nex_k_block(fh_sbe_nex_k, 0.0d0, nk, gs%kpoint, pop_k)
        call write_sbe_nex_k_block(fh_sbe_nex_k_real, 0.0d0, nk, gs%kpoint, pop_k)
        flush(fh_sbe_nex_k)
        flush(fh_sbe_nex_k_real)
        if (.not. gs%have_unfold) then
            call write_sbe_nex_k_lev_block(fh_sbe_nex_k_lev_real, 0.0d0, nk, gs%kpoint, pop4)
            flush(fh_sbe_nex_k_lev_real)
        end if
        if (gs%have_unfold) then
            pop_lev_k = 0.0d0
            call write_sbe_nex_k_unfold_block(fh_sbe_nex_k_unfold, 0.0d0, nk, &
                & gs%kpoint, gs%unfold_offset, pop_lev_k, gs%n_coset)
            call write_sbe_nex_k_unfold_block(fh_sbe_nex_k_unfold_real, 0.0d0, nk, &
                & gs%kpoint, gs%unfold_offset, pop_lev_k, gs%n_coset)
            flush(fh_sbe_nex_k_unfold)
            flush(fh_sbe_nex_k_unfold_real)
        end if
    end if

    call comm_sync_all(icomm)

    ! B4: checkpoint restart -- resume rho / X_branch / step index from the
    ! per-rank stream file (same nproc required; the field is recomputed
    ! deterministically from the input, so nothing else needs saving).
    it0 = 1
    if (yn_sbe_checkpoint_restart == 'y') then
        write(ck_file, '(a,a,a,i5.5,a)') trim(base_directory), trim(sysname), '_sbe_ckpt_r', irank, '.bin'
        inquire(file=trim(ck_file), exist=ck_exists)
        if (.not. ck_exists) then
            if (irank == 0) write(*,'(a)') '# ERROR: yn_sbe_checkpoint_restart but no checkpoint file'
            error stop 'B4: checkpoint file missing'
        end if
        open(newunit=ck_unit, file=trim(ck_file), form='unformatted', access='stream', action='read')
        read(ck_unit) it0, energy, sbe%led_dn, sbe%led_de
        read(ck_unit) sbe%rho(:, :, sbe%ik_min:sbe%ik_max)
        if (allocated(sbe%X_branch)) read(ck_unit) sbe%X_branch(:, sbe%ik_min:sbe%ik_max)
        close(ck_unit)
        it0 = it0 + 1
        if (irank == 0) write(*,'(a,i8)') '# B4: resumed from checkpoint, continuing at step ', it0
    end if

    ! Realtime calculation
    do it = it0, nt
        t = dt * it
        
        !---------------------------------------------------------------
        ! CF4(Gauss-Legendre)+Yoshida unitary step combined with strictly
        ! CPTP Strang/Hadamard Kuhn-Zurek dephasing: evolve rho from
        ! t=(it-1)*dt to t=it*dt. The propagator interpolates A(t) at the
        ! internal Gauss-Legendre/Yoshida sub-nodes from the field values
        ! at the step endpoints supplied here.
        !---------------------------------------------------------------
        call dt_evolve_bloch_cf4(sbe, gs, t - dt, dt, Ac_ext_t(:, it - 1), Ac_ext_t(:, it))
        
        !---------------------------------------------------------------
        ! Calculate Current J(t) at t=it*dt (after evolution)
        !---------------------------------------------------------------
        call calc_current_bloch(sbe, gs, Ac_ext_t(:, it), Jmat, icomm)
        
        !---------------------------------------------------------------
        ! Calculate E-field at t=it*dt
        ! E = -dA/dt. Use central difference: E(t) = -(A(t+dt) - A(t-dt)) / (2*dt)
        ! Ac_ext_t is allocated from -1 to nt+1, safe for all it >= 1
        !---------------------------------------------------------------
        E(:) = -(Ac_ext_t(:, it+1) - Ac_ext_t(:, it-1)) / (2.0d0 * dt)
        
        ! Energy update: dW = -E·J·V·dt (work done by field on electrons)
        energy = energy + dot_product(E(1:3), -Jmat(1:3)) * gs%volume * dt
        
        if (irank == 0) then
            call write_sbe_rt_line(fh_sbe_rt, &
                & t, Ac_ext_t(1:3, it), E(1:3), Ac_ext_t(1:3, it), E(1:3), Jmat(1:3))
        end if

        ! Intra-band (Houston-basis) current -- physical intra/inter split in
        ! the velocity gauge (the total J above is the gauge-invariant sum).
        ! Written every step like the total current, for direct comparison.
        if (yn_out_intraband_current == 'y') then
            call calc_intraband_current_houston(sbe, gs, Ac_ext_t(:, it), jmat_intra, icomm)
            if (irank == 0) &
                call write_sbe_intra_current_line(fh_sbe_intra_current, t, jmat_intra(1:3))
        end if

        if (mod(it, out_rt_energy_step) == 0) then
            tr_all = calc_trace(sbe, gs, nstate_sbe(1), icomm)
            if (irank == 0) then
                call write_sbe_rt_energy_line(fh_sbe_rt_energy, t, energy, energy)
                write(*, "(i6,f12.3,3es12.3,2f12.3)") it, t * au_fs, Jmat(1:3), tr_all, energy
                ! C1: cumulative per-channel ledger (identical on all ranks)
                write(fh_sbe_channels, '(f14.4,8es14.5)') t, &
                    sbe%led_dn(1), sbe%led_de(1), sbe%led_dn(2), sbe%led_de(2), &
                    sbe%led_dn(3), sbe%led_de(3), sbe%led_dn(4), sbe%led_de(4)
            end if
        end if

        ! B4: periodic checkpoint (per-rank stream file, overwritten in place)
        if (sbe_checkpoint_step > 0) then
            if (mod(it, sbe_checkpoint_step) == 0) then
                write(ck_file, '(a,a,a,i5.5,a)') trim(base_directory), trim(sysname), '_sbe_ckpt_r', irank, '.bin'
                open(newunit=ck_unit, file=trim(ck_file), form='unformatted', access='stream', action='write')
                write(ck_unit) it, energy, sbe%led_dn, sbe%led_de
                write(ck_unit) sbe%rho(:, :, sbe%ik_min:sbe%ik_max)
                if (allocated(sbe%X_branch)) write(ck_unit) sbe%X_branch(:, sbe%ik_min:sbe%ik_max)
                close(ck_unit)
            end if
        end if

        if (mod(it, out_projection_step) == 0) then
            tr_all = calc_trace(sbe, gs, nstate_sbe(1), icomm)
            tr_vb = calc_trace(sbe, gs, nb_vb, icomm)
            if (irank == 0) then
                call write_sbe_nex_line(fh_sbe_nex, t, (tr_all - tr_vb) / gs%volume, (nelec - tr_vb) / gs%volume)
            end if
        end if

        ! Houston-basis population of the lowest conduction band, per k-point.
        ! Written far less often than _sbe_nex.data (out_projection_k_step,
        ! default 10x out_projection_step) since this output scales with nk.
        if (mod(it, out_projection_k_step) == 0) then
            call calc_bloch_population_k(sbe, gs, Ac_ext_t(:, it), ib_lcb, pop_k, icomm)
            ! Real carriers only (diabatic / fixed-basis LCB occupation)
            call calc_diabatic_population_k(sbe, ib_lcb, pop_k_real, icomm)
            if (irank == 0) then
                call write_sbe_nex_k_block(fh_sbe_nex_k, t, nk, gs%kpoint, pop_k)
                call write_sbe_nex_k_block(fh_sbe_nex_k_real, t, nk, gs%kpoint, pop_k_real)
            end if
            ! Primitive cell: 4 gap-edge diabatic populations (VB-1, VB, CB1, CB2)
            ! for the --spectral 4-band colouring (uses pop_k_real as scratch).
            if (.not. gs%have_unfold) then
                pop4 = 0d0
                call calc_diabatic_population_k(sbe, ib_lcb-2, pop_k_real, icomm); pop4(1,:) = pop_k_real
                call calc_diabatic_population_k(sbe, ib_lcb-1, pop_k_real, icomm); pop4(2,:) = pop_k_real
                call calc_diabatic_population_k(sbe, ib_lcb,   pop_k_real, icomm); pop4(3,:) = pop_k_real
                call calc_diabatic_population_k(sbe, ib_lcb+1, pop_k_real, icomm); pop4(4,:) = pop_k_real
                if (irank == 0) call write_sbe_nex_k_lev_block(fh_sbe_nex_k_lev_real, t, nk, gs%kpoint, pop4)
            end if
            ! Physical (unfolded) CB1 populations per primitive BZ point
            if (gs%have_unfold) then
                call calc_unfolded_population_k(sbe, gs, Ac_ext_t(:, it), pop_lev_k, icomm)
                call calc_diabatic_unfolded_population_k(sbe, gs, pop_lev_k_real, icomm)
                if (irank == 0) then
                    call write_sbe_nex_k_unfold_block(fh_sbe_nex_k_unfold, t, nk, &
                        & gs%kpoint, gs%unfold_offset, pop_lev_k, gs%n_coset)
                    call write_sbe_nex_k_unfold_block(fh_sbe_nex_k_unfold_real, t, nk, &
                        & gs%kpoint, gs%unfold_offset, pop_lev_k_real, gs%n_coset)
                end if
            end if

            ! VG basis-sufficiency monitor (criterion (a)): peak adiabatic
            ! occupation of the top retained band over the BZ. If it exceeds the
            ! tolerance the field is pushing population to the basis edge -- the
            ! band budget N_b is too small (a separate axis from the PW cutoff,
            ! and NOT cured by the Houston basis). Warn on the error channel and
            ! CONTINUE; the user re-runs with more bands and an N_b convergence
            ! study (criterion (b)) to confirm.
            if (ib_top > 0) then
                call calc_bloch_population_k(sbe, gs, Ac_ext_t(:, it), ib_top, pop_top_k, icomm)
                ptop = maxval(pop_top_k)
                if (irank == 0 .and. vg_ptop_exceeds(ptop, PTOP_TOL)) then
                    write(error_unit, '(a,es10.3,a,es10.3,a,f10.3,a,i0,a)') &
                        ' WARNING: VG basis edge reached -- P_top = ', ptop, &
                        ' > ', PTOP_TOL, ' at t = ', t * au_fs, &
                        ' fs (top band ', ib_top, '). Increase N_b (nstate) and re-check convergence.'
                end if
            end if
        end if

        if (mod(it, 500) == 0) then
            if (irank == 0) then
                flush(fh_sbe_rt)
                flush(fh_sbe_rt_energy)
                flush(fh_sbe_nex)
                flush(fh_sbe_nex_k)
                flush(fh_sbe_nex_k_real)
                if (gs%have_unfold) then
                    flush(fh_sbe_nex_k_unfold)
                    flush(fh_sbe_nex_k_unfold_real)
                else
                    flush(fh_sbe_nex_k_lev_real)
                end if
                if (yn_out_intraband_current == 'y') flush(fh_sbe_intra_current)
            end if
        end if
    end do

    call comm_sync_all(icomm)

    if (irank == 0) then
        close(fh_sbe_rt)
        close(fh_sbe_rt_energy)
        close(fh_sbe_nex)
        close(fh_sbe_nex_k)
        close(fh_sbe_nex_k_real)
        if (gs%have_unfold) then
            close(fh_sbe_nex_k_unfold)
            close(fh_sbe_nex_k_unfold_real)
        else
            close(fh_sbe_nex_k_lev_real)
        end if
        if (yn_out_intraband_current == 'y') close(fh_sbe_intra_current)
    end if

    deallocate(pop_k, pop_k_real)
    deallocate(pop_top_k)
    if (allocated(pop_lev_k)) deallocate(pop_lev_k)
    if (allocated(pop_lev_k_real)) deallocate(pop_lev_k_real)
    if (allocated(pop4)) deallocate(pop4)

    return
end subroutine main_realtime_ssbe

end module realtime_ssbe
