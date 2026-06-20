module bloch_solver_ssbe
    use math_constants, only: pi, zi
    use phys_constants, only: au_ev
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use gs_info_ssbe
    use util_ssbe, only: split_range
    implicit none

    private
    public :: s_sbe_bloch_solver, init_sbe_bloch_solver, calc_current_bloch, &
              dt_evolve_bloch_cf4, calc_trace, calc_energy, calc_bloch_population_k, &
              calc_unfolded_population_k

    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
        logical :: flag_vnl_correction

        ! Frozen core handling
        logical, allocatable :: is_active(:) ! .true. if band is active, .false. if frozen
        integer :: n_active_bands = 0
        integer, allocatable :: active_idx(:)  ! Mapping: 1..n_active -> global band index

        ! Houston-basis branch (wave-packet) positions X_a(k,t) used by the
        ! Kuhn-Zurek/Caldeira-Leggett dephasing kernel. By the origin-shift
        ! invariance of the dephasing map (only differences X_a-X_b enter),
        ! the choice X_a(0)=0 is physically irrelevant; it merely fixes a
        ! reproducible reference for restarts.
        real(8), allocatable :: X_branch(:, :)  ! (1:nb, ik_min:ik_max)

        ! Kuhn-Zurek/Caldeira-Leggett decoherence: lambda = kB*T / tau_m
        real(8) :: lambda_decoh = 0d0
        logical :: flag_decoh   = .false.

        ! k-local impact-ionization Lindblad channel (Stobbe rate fit):
        !   gamma(e_kin) = P (e_kin - E_th)^4 Theta(e_kin - E_th),
        ! e_kin measured from the global field-free CBM. Each event is split
        ! (Hartree-Fock / two-particle closure) into two frozen-rate
        ! amplitude-damping channels: primary relaxation h -> h' and cold-pair
        ! creation v1 -> c1, applied in the same Houston basis as the
        ! Kuhn-Zurek dephasing (no extra ZHEEV).
        logical :: flag_impact  = .false.
        real(8) :: ii_pref_au   = 0d0   ! P in 1/(Ha^a a.u.time)
        real(8) :: ii_exponent  = 4d0   ! fit exponent a (4 GaAs Stobbe, 2 Si Keldysh)
        real(8) :: ii_eth_au    = 0d0   ! threshold E_th [Ha]
        real(8) :: ii_ramp_au   = 0d0   ! linear Theta-smoothing width [Ha]
        real(8) :: ii_ecbm_au   = 0d0   ! global CBM of the field-free bands [Ha]
        real(8) :: ii_eg_au     = 0d0   ! band gap E_g [Ha] (primary loses E_g)
        integer :: nv_act       = 0     ! valence branches inside the active subspace
        real(8) :: occ_max      = 2d0   ! 2 (scalar bands) / 1 (spinor bands)

        ! Sublattice (band-unfolding) awareness of the impact-ionization
        ! channel. When the unfolding weights gs%unfold_w are available
        ! (non-zero), the k-local two-particle event is resolved per FCC
        ! sublattice so that the secondary pair is created in the SAME
        ! primitive-cell sector (same primitive crystal momentum) as the
        ! primary -- restoring primitive momentum conservation that the
        ! cubic folding hides, and removing the spurious cross-sublattice
        ! ("false generation") events. Reduces identically to the folded
        ! channel when the weights are absent.
        logical :: flag_unfold_ii = .false.

        ! Coulomb (time-dependent Hartree-Fock / exchange) renormalization,
        ! Golde-Kira-Meier-Koch SBE (Phys. Status Solidi B 248, 863 (2011),
        ! Eqs. 4-5). The exchange (Fock) self-energy
        !   Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q)
        ! is a NON-k-local mean field (couples all k). Added to H_VG it
        ! reproduces both the renormalized single-particle energies
        ! (diagonal: eps~ = eps - sum_q V_{k-q} f_q) and the renormalized
        ! Rabi frequency (off-diagonal: Omega = d.E + sum_q V_{k-q} p_q),
        ! with the (1-f_e-f_h) Pauli factor coming for free from the
        ! von Neumann commutator. Frozen over a dt step (mean-field predictor).
        logical :: flag_coulomb = .false.
        real(8) :: coul_pref     = 0d0  ! strength * 4 pi / (eps * Omega_cell * Nk)
        real(8) :: coul_screen2  = 0d0  ! kappa^2 [1/Bohr^2] (Yukawa regulariser)
        integer :: icomm         = 0    ! MPI communicator (for the non-local exchange sum)
        complex(8), allocatable :: sigma_hf(:, :, :) ! (nba, nba, ik_min:ik_max) exchange Sigma
    end type

    !=========================================================================
    ! CF4 (commutator-free Magnus 4) + Suzuki-Yoshida composition constants
    !=========================================================================
    ! Two-point Gauss-Legendre nodes on [0,1]: c = 1/2 -+ sqrt(3)/6
    real(8), parameter :: cf4_c1 = 0.21132486540518713d0
    real(8), parameter :: cf4_c2 = 0.78867513459481287d0
    ! CF4 combination weights: alpha = 1/4 -+ sqrt(3)/6
    real(8), parameter :: cf4_alpha1 =  0.53867513459481287d0
    real(8), parameter :: cf4_alpha2 = -0.03867513459481287d0
    ! Suzuki-Yoshida triple-jump constants (4th-order composition of a
    ! 2nd-order base scheme): p1 + p2 + p1 = 1
    real(8), parameter :: yoshida_p1 =  1.35120719196d0
    real(8), parameter :: yoshida_p2 = -1.70241438392d0

contains

subroutine init_sbe_bloch_solver(sbe, gs, nb_sbe, icomm)
    use util_ssbe
    use communication, only: comm_get_groupinfo, comm_summation, comm_bcast
    use salmon_global, only: frozen_core_threshold_ev, frozen_free_threshold_ev, &
                             sbe_decoh_temperature_k, sbe_decoh_tau_m_fs, yn_sbe_spinor, &
                             yn_sbe_impact_ionization, sbe_ii_prefactor, &
                             sbe_ii_threshold_ev, sbe_ii_ramp_ev, &
                             sbe_ii_form, sbe_ii_exponent, &
                             yn_sbe_coulomb, sbe_coulomb_epsilon, &
                             sbe_coulomb_strength, sbe_coulomb_screen_au
    use math_constants, only: pi
    use phys_constants, only: au_fs, kB_au, au_ev
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr, count_active
    integer, allocatable :: itbl_min(:), itbl_max(:)
    real(8) :: eigen_ev, fermi_energy_ev
    integer, allocatable :: is_active_buf(:)
    integer :: homo_idx, lumo_idx

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    ! =========================================================================
    ! ИНИЦИАЛИЗАЦИЯ rho: ИСПОЛЬЗУЕМ gs%occup (КАК В ОРИГИНАЛЕ!)
    ! =========================================================================
    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)  ! ← КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ
        end do
    end do
    
    sbe%flag_vnl_correction = .false.

    ! =========================================================================
    ! Houston-basis branch positions X_a(k,t): explicit zero-init.
    ! By the invariance theorem of the dephasing kernel under a global shift
    ! of the X_a origin (only X_a-X_b enters exp[-lambda(X_a-X_b)^2 tau]),
    ! the choice X_a(t=0) = 0 carries no physical content; it only fixes a
    ! reproducible convention for restarts.
    ! =========================================================================
    allocate(sbe%X_branch(1:sbe%nb, sbe%ik_min:sbe%ik_max))
    sbe%X_branch = 0d0

    ! =========================================================================
    ! Kuhn-Zurek/Caldeira-Leggett decoherence strength: lambda = kB*T / tau_m
    ! Both temperature and relaxation time must be positive to enable it;
    ! otherwise the dephasing map reduces identically to the identity (D=0),
    ! which is trivially CPTP.
    ! =========================================================================
    if (sbe_decoh_temperature_k > 0d0 .and. sbe_decoh_tau_m_fs > 0d0) then
        sbe%lambda_decoh = kB_au * sbe_decoh_temperature_k / (sbe_decoh_tau_m_fs / au_fs)
        sbe%flag_decoh   = .true.
    else
        sbe%lambda_decoh = 0d0
        sbe%flag_decoh   = .false.
    end if

    ! =========================================================================
    ! ЛОГИКА is_active (Frozen Core)
    ! =========================================================================
    
    ! 1. Calculate Fermi Energy (Assumes closed-shell / even number of electrons)

    if (yn_sbe_spinor == 'y') then
        ! Spinor bands: one electron per band, ne occupied bands
        homo_idx = gs%ne
    else
        homo_idx = gs%ne / 2
        if (mod(gs%ne, 2) /= 0 .and. irank == 0) then
            write(*, '(a)') 'WARNING: Odd number of electrons. Fermi energy assumes closed-shell.'
        end if
    end if
    lumo_idx = homo_idx + 1

    ! gs%eigen is stored in atomic units (Hartree). Convert to eV here so the
    ! frozen-core window thresholds (frozen_core/free_threshold_ev, genuine eV
    ! inputs) are compared in eV as named, and the diagnostic labels are honest.
    ! This only affects which bands are flagged active and the printout; the
    ! dynamics always use gs%eigen in a.u. directly.
    fermi_energy_ev = ((gs%eigen(homo_idx, 1) + gs%eigen(lumo_idx, 1)) * 0.5d0) * au_ev

    ! 2. Initialize active bands array
    allocate(sbe%is_active(1:sbe%nb))
    sbe%is_active = .false.  
    sbe%n_active_bands = 0   
    
    ! 3. Determine active bands on root rank
    if (irank == 0) then
        do ib = 1, sbe%nb
            eigen_ev = gs%eigen(ib, 1) * au_ev
            ! Note: Ensure frozen_core_threshold_ev is negative if it represents a window below E_F
            if (eigen_ev > fermi_energy_ev + frozen_core_threshold_ev .and. &
                eigen_ev < fermi_energy_ev + frozen_free_threshold_ev) then
                sbe%is_active(ib) = .true.
                sbe%n_active_bands = sbe%n_active_bands + 1
            end if
        end do
    end if
    
    ! 4. Broadcast n_active_bands
    call comm_bcast(sbe%n_active_bands, icomm, 0)
    
    ! 5. Broadcast is_active logical array
    if (sbe%nb > 0) then
        allocate(is_active_buf(1:sbe%nb))
        
        if (irank == 0) then
            ! Modern Fortran: use merge() instead of a verbose do-loop
            is_active_buf = merge(1, 0, sbe%is_active)
        end if
        
        call comm_bcast(is_active_buf, icomm, 0)
        
        ! Element-wise comparison replaces the verbose do-loop
        sbe%is_active = (is_active_buf == 1)
        
        deallocate(is_active_buf)
    end if

    ! 6. Build active_idx array
    if (sbe%n_active_bands > 0) then
        allocate(sbe%active_idx(sbe%n_active_bands))
        count_active = 0
        do ib = 1, sbe%nb
            if (sbe%is_active(ib)) then
                count_active = count_active + 1
                sbe%active_idx(count_active) = ib
            end if
        end do
    else
        ! Modern Fortran handles zero-sized arrays natively. 
        ! If downstream legacy code crashes on size 0, revert to allocate(sbe%active_idx(1))
        allocate(sbe%active_idx(0))  
    end if

    ! =========================================================================
    ! k-local impact-ionization channel (optional, yn_sbe_impact_ionization)
    ! =========================================================================
    sbe%occ_max = merge(1d0, 2d0, yn_sbe_spinor == 'y')
    sbe%flag_impact = (yn_sbe_impact_ionization == 'y')
    if (sbe%flag_impact) then
        ! Fit-form exponent a: GaAs Stobbe quartic (a=4, hard threshold)
        ! [Stobbe-Redmer-Schattke, PRB 49, 4494 (1994)]; Si Keldysh quadratic
        ! (a=2, soft near-gap threshold) [Keldysh, JETP 21, 1135 (1965);
        ! Cartier et al., APL 62, 3339 (1993)]; Si full-band option a=4.6
        ! [Kamakura et al., JAP 75, 3500 (1994)].
        sbe%ii_exponent = sbe_ii_exponent
        ! Prefactor P [s^-1 eV^-a] -> [1/(Ha^a a.u.time)]:
        ! rate_au = P * t_au[s] * (dE[Ha] * au_ev)^a
        sbe%ii_pref_au = sbe_ii_prefactor * (au_fs * 1d-15) * au_ev**sbe_ii_exponent
        sbe%ii_eth_au  = sbe_ii_threshold_ev / au_ev
        sbe%ii_ramp_au = sbe_ii_ramp_ev / au_ev
        ! Global CBM of the field-free band structure (kinetic-energy zero of
        ! the Stobbe fit) and the gap lost by the primary electron per event.
        if (homo_idx + 1 > gs%nb) stop "impact ionization: no conduction bands"
        sbe%ii_ecbm_au = minval(gs%eigen(homo_idx + 1, :))
        sbe%ii_eg_au   = gs%eg_au
        ! Valence branches inside the active subspace: v1 = sbe%nv_act,
        ! c1 = sbe%nv_act + 1 in active (energy-ordered Houston) indexing.
        sbe%nv_act = 0
        do ib = 1, sbe%n_active_bands
            if (sbe%active_idx(ib) <= homo_idx) sbe%nv_act = sbe%nv_act + 1
        end do
        ! Sublattice resolution is enabled iff the unfolding weights were
        ! loaded (gs%unfold_w not all zero). When absent, the channel falls
        ! back to the original folded (single-pool) treatment.
        sbe%flag_unfold_ii = .false.
        if (allocated(gs%unfold_w)) then
            if (maxval(abs(gs%unfold_w)) > 1d-12) sbe%flag_unfold_ii = .true.
        end if
        if (irank == 0) then
            write(*, '(a)') '# impact ionization (k-local Lindblad, Stobbe fit) enabled:'
            write(*, '(a,ES12.5,a)') '#   P     = ', sbe%ii_pref_au, ' 1/(Ha^4 a.u.t)'
            write(*, '(a,ES12.5,a,ES12.5,a)') '#   E_th  = ', sbe%ii_eth_au, ' Ha, ramp = ', sbe%ii_ramp_au, ' Ha'
            write(*, '(a,ES12.5,a,ES12.5,a)') '#   E_CBM = ', sbe%ii_ecbm_au, ' Ha, E_g  = ', sbe%ii_eg_au, ' Ha'
            write(*, '(a,i4,a,i4)') '#   active valence branches = ', sbe%nv_act, ' / ', sbe%n_active_bands
            if (sbe%flag_unfold_ii) then
                write(*, '(a)') '#   sublattice-resolved (unfolding weights loaded): primitive'
                write(*, '(a)') '#   momentum conservation enforced over 4 FCC sublattices'
            else
                write(*, '(a)') '#   folded (single-pool) treatment: no unfolding weights found'
            end if
            if (sbe%nv_act < 1 .or. sbe%nv_act >= sbe%n_active_bands) &
                write(*, '(a)') '#   WARNING: no valence or no conduction branches active -- channel is inert'
        end if
    end if

    ! =========================================================================
    ! Coulomb (time-dependent Hartree-Fock / exchange) renormalization
    ! [Golde-Kira-Meier-Koch, Phys. Status Solidi B 248, 863 (2011)].
    ! Exchange self-energy Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q) with the
    ! statically-screened Coulomb V(q) = strength * 4 pi / (eps Omega Nk (q^2+kappa^2)).
    ! Stored in the velocity-gauge Bloch basis and added to H_VG; the Houston
    ! basis (eigenbasis of H_VG + Sigma) the dissipators use becomes the
    ! Coulomb-renormalized one. Off in the no-Coulomb limit (yn_sbe_coulomb='n').
    ! =========================================================================
    sbe%flag_coulomb = (yn_sbe_coulomb == 'y')
    if (sbe%flag_coulomb) then
        ! Discrete exchange-sum prefactor: continuum (1/(2pi)^3) int d^3q  ->
        ! (1/(Omega_cell Nk)) sum_grid; V(q) = 4 pi / (eps q^2) in a.u.
        sbe%coul_pref = sbe_coulomb_strength * 4d0 * pi &
                      / (max(sbe_coulomb_epsilon, 1d-12) * gs%volume * dble(gs%nk))
        sbe%coul_screen2 = sbe_coulomb_screen_au**2
        sbe%icomm = icomm
        if (sbe%n_active_bands > 0) &
            allocate(sbe%sigma_hf(1:sbe%n_active_bands, 1:sbe%n_active_bands, &
                                  sbe%ik_min:sbe%ik_max))
        if (irank == 0) then
            write(*, '(a)') '# Coulomb HF (Golde-Kira-Meier-Koch SBE) enabled:'
            write(*, '(a,ES12.5)') '#   exchange prefactor 4pi*str/(eps*Omega*Nk) = ', sbe%coul_pref
            write(*, '(a,f8.3,a,ES12.5,a)') '#   eps = ', sbe_coulomb_epsilon, &
                ', screening kappa^2 = ', sbe%coul_screen2, ' 1/Bohr^2'
            write(*, '(a)') '#   NOTE: non-k-local mean field, O(Nk^2) per step (frozen over dt)'
        end if
    end if

    ! 7. Diagnostic Print
    if (irank == 0) then
        write(*, '(a)') '=========================================='
        write(*, '(a)') 'DIAGNOSTIC: Frozen Core Check'
        write(*, '(a, f8.2, a)') '  frozen_core_threshold_ev = ', frozen_core_threshold_ev, ' eV'
        write(*, '(a, f8.2, a)') '  frozen_free_threshold_ev = ', frozen_free_threshold_ev, ' eV'
        write(*, '(a, f12.4, a)') '  Fermi energy (eV)      = ', fermi_energy_ev, ' eV'
        write(*, '(a, i4, a, i4)') '  n_active_bands         = ', sbe%n_active_bands, ' / ', sbe%nb
        write(*, '(a)') '----------------------------------------'
        write(*, '(a)') '  Band energies relative to Fermi level:'
        
        do ib = 1, min(sbe%nb, 100)  ! Print first 100 bands
            eigen_ev = gs%eigen(ib, 1) * au_ev
            write(*, '(a, i3, a, f10.4, a, f8.2, a, l1)') &
                '    Band ', ib, ': E = ', eigen_ev, ' eV, E-E_F = ', &
                (eigen_ev - fermi_energy_ev), ' eV, active = ', sbe%is_active(ib)
        end do
        
        if (sbe%nb > 100) write(*, '(a)') '    ... (more bands)'
        write(*, '(a)') '=========================================='
    end if

end subroutine

subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, nb
    complex(8) :: tmp1(1:3), tmp(1:3), v_mat(1:sbe%nb, 1:sbe%nb)
    complex(8) :: trace_val

    nb = sbe%nb
    tmp1 = 0d0

    !$omp parallel do default(shared) private(ik, idir, ib, jb, v_mat, trace_val) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            v_mat = gs%p_tm_matrix(:, :, idir, ik)
            do ib = 1, nb
                v_mat(ib, ib) = v_mat(ib, ib) + Ac(idir)
            end do
            if (sbe%flag_vnl_correction) then
                v_mat = v_mat + gs%rvnl_tm_matrix(:, :, idir, ik)
            endif

            trace_val = 0d0
            do ib = 1, nb
                do jb = 1, nb
                    trace_val = trace_val + v_mat(ib, jb) * sbe%rho(jb, ib, ik)
                end do
            end do
            tmp1(idir) = tmp1(idir) + gs%kweight(ik) * trace_val
        end do
    end do
    !$omp end parallel do

    call comm_summation(tmp1, tmp, 3, icomm)
    jmat(:) = real(tmp(:)) / (sum(gs%kweight) * gs%volume)
end subroutine calc_current_bloch


function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication, only: comm_get_groupinfo, comm_summation
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1

    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, nb_max
            tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    tr = tmp / sum(gs%kweight)

    return
end function calc_trace


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs_info), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir
    real(8) :: tmp1, tmp, energy
    
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_mod_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    energy = tmp / sum(gs%kweight)

    return
end function calc_energy


!=============================================================================
! CF4 (commutator-free Magnus, 4th order) propagator on Gauss-Legendre nodes,
! composed via the Suzuki-Yoshida triple-jump for the unitary part, combined
! with a strictly CPTP Kuhn-Zurek/Caldeira-Leggett dephasing map applied
! through Strang splitting with an exact Hadamard/Gaussian kernel in the
! instantaneous (Houston) eigenbasis:
!
!   rho(t+h) = D(h/2) o [ S2(p1 h) o S2(p2 h) o S2(p1 h) ] o D(h/2) [rho(t)]
!
! IMPORTANT: the Suzuki-Yoshida composition wraps ONLY the unitary CF4
! sub-steps S2(.), never the dephasing map D(.). A unitary step run with a
! negative sub-step (p2 h < 0) is simply a unitary rotation run backwards in
! time -- always valid. A dephasing step run with tau < 0 would replace the
! Hadamard/Gaussian kernel exp[-lambda (X_a-X_b)^2 tau] by its reciprocal,
! exp[+lambda (X_a-X_b)^2 |tau|], which is not positive semi-definite (it
! fails the Schoenberg/Bochner criterion for an RBF kernel and the Schur
! product theorem for the Hadamard map), and would break completely positive
! trace preservation. Hence D(h/2) is applied twice, with tau=+h/2>0 each
! time, by Strang splitting around the (always-safe) unitary composition.
!=============================================================================

subroutine dt_evolve_bloch_cf4(sbe, gs, t_start, dt, Ac_begin, Ac_end)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info), intent(inout) :: gs
    real(8), intent(in) :: t_start  ! time at the beginning of the step, rho(t_start) -> rho(t_start+dt)
    real(8), intent(in) :: dt
    real(8), intent(in) :: Ac_begin(1:3)  ! external vector potential A(t_start)
    real(8), intent(in) :: Ac_end(1:3)    ! external vector potential A(t_start+dt)

    real(8) :: tau_sub(3), t_sub(3)
    real(8) :: t_node(2, 3), s_node
    real(8) :: Ac_node(1:3, 2, 3)
    integer :: isub

    integer :: ik, nb, nba, i, j, idir, in, im

    complex(8), allocatable :: p_active(:, :, :)
    complex(8), allocatable :: rho_a(:, :)
    complex(8), allocatable :: H1(:, :), H2(:, :), HVG(:, :)
    real(8),    allocatable :: eigen_active(:)
    real(8),    allocatable :: V_begin(:), V_end(:), X_a(:)
    real(8),    allocatable :: w_act_sub(:, :)   ! (4, nba) field-free sublattice weights of active bands
    complex(8), allocatable :: p_k_full(:, :, :)
    complex(8), allocatable :: rho_n_full(:, :)
    integer :: s

    nb  = sbe%nb
    nba = sbe%n_active_bands

    !-------------------------------------------------------------------------
    ! The external field is known only at the step endpoints (the analytic
    ! pulse in realtime_ssbe, or the macroscopic Maxwell field in
    ! multiscale_ssbe -- both callers supply Ac(t_start) and Ac(t_start+dt)).
    ! CF4(Gauss-Legendre)+Yoshida needs A at several intermediate sub-nodes;
    ! we obtain them by linear interpolation in time,
    !   A(t_start + s*dt) = (1-s) Ac_begin + s Ac_end,   s in [0,1],
    ! which is consistent with the existing multiscale convention (compare
    ! the "linear interpolation for A(t+dt/2)" used by the previous ETDRK4
    ! step) and strictly more accurate than the old approach of treating A as
    ! constant over the whole step. These nodes are identical for every
    ! k-point, so they are evaluated once before the OpenMP/k-point loop.
    !-------------------------------------------------------------------------
    tau_sub(1) = yoshida_p1 * dt
    tau_sub(2) = yoshida_p2 * dt
    tau_sub(3) = yoshida_p1 * dt

    t_sub(1) = t_start
    t_sub(2) = t_sub(1) + tau_sub(1)
    t_sub(3) = t_sub(2) + tau_sub(2)
    ! t_sub(3) + tau_sub(3) = t_start + dt, since p1 + p2 + p1 = 1

    do isub = 1, 3
        t_node(1, isub) = t_sub(isub) + cf4_c1 * tau_sub(isub)
        t_node(2, isub) = t_sub(isub) + cf4_c2 * tau_sub(isub)

        s_node = (t_node(1, isub) - t_start) / dt
        Ac_node(:, 1, isub) = (1d0 - s_node) * Ac_begin + s_node * Ac_end

        s_node = (t_node(2, isub) - t_start) / dt
        Ac_node(:, 2, isub) = (1d0 - s_node) * Ac_begin + s_node * Ac_end
    end do

    ! Coulomb HF mean field, frozen over this dt step: form the non-k-local
    ! exchange self-energy Sigma(k) from rho(t_start) BEFORE the (OpenMP) k-loop
    ! (it is an MPI all-gather + convolution, must run collectively once).
    if (sbe%flag_coulomb .and. nba > 0) call compute_coulomb_selfenergy(sbe, gs)

    !$omp parallel default(shared) &
    !$omp    private(ik, i, j, idir, in, im, isub, s) &
    !$omp    private(p_active, rho_a, H1, H2, HVG, eigen_active, V_begin, V_end, X_a, w_act_sub) &
    !$omp    private(p_k_full, rho_n_full)

    if (nba > 0) then
        allocate(p_active(nba, nba, 3), rho_a(nba, nba))
        allocate(H1(nba, nba), H2(nba, nba), HVG(nba, nba))
        allocate(eigen_active(nba), V_begin(nba), V_end(nba), X_a(nba))
        allocate(w_act_sub(4, nba))
    end if
    allocate(p_k_full(nb, nb, 1:3), rho_n_full(nb, nb))

    !$omp do
    do ik = sbe%ik_min, sbe%ik_max

        if (nba > 0) then
            p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
            if (sbe%flag_vnl_correction) p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)
            rho_n_full(:, :) = sbe%rho(:, :, ik)

            ! Restrict to the active subspace (frozen core/free bands excluded)
            do idir = 1, 3
                do j = 1, nba
                    im = sbe%active_idx(j)
                    do i = 1, nba
                        in = sbe%active_idx(i)
                        p_active(i, j, idir) = p_k_full(in, im, idir)
                    end do
                end do
            end do
            do i = 1, nba
                eigen_active(i) = gs%eigen(sbe%active_idx(i), ik)
            end do
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                rho_a(i, j) = rho_n_full(in, im)
            end do; end do
            do i = 1, nba
                X_a(i) = sbe%X_branch(sbe%active_idx(i), ik)
            end do

            ! Field-free sublattice weights of the active bands at this k
            ! (projected onto the instantaneous Houston branches inside
            ! houston_dissipate). Only consumed by the unfolding-aware
            ! impact-ionization channel.
            if (sbe%flag_impact .and. sbe%flag_unfold_ii) then
                do i = 1, nba
                    do s = 1, 4
                        w_act_sub(s, i) = gs%unfold_w(s, sbe%active_idx(i), ik)
                    end do
                end do
            else if (sbe%flag_impact) then
                w_act_sub = 0d0
            end if

            !-----------------------------------------------------------------
            ! Step 1: D(h/2) -- Strang dissipative half-step: Hadamard
            ! Kuhn-Zurek dephasing and/or impact-ionization channels, both in
            ! the same Houston basis (one shared ZHEEV), tau = +h/2 > 0
            !-----------------------------------------------------------------
            if (sbe%flag_decoh .or. sbe%flag_impact) then
                call build_HVG(nba, eigen_active, p_active, Ac_begin, HVG)
                if (sbe%flag_coulomb) HVG = HVG + sbe%sigma_hf(:, :, ik)
                call houston_dissipate(sbe, nba, rho_a, HVG, p_active, Ac_begin, X_a, &
                                       0.5d0 * dt, V_begin, w_act_sub)
            else
                V_begin = 0d0
            end if

            !-----------------------------------------------------------------
            ! Step 2: S4_unitary = S2(p1 h) o S2(p2 h) o S2(p1 h)
            ! Each S2(tau) is a CF4 (two-exponential) commutator-free Magnus
            ! step on the two Gauss-Legendre nodes spanning that sub-interval.
            ! A negative tau (the middle Yoshida jump) is just a backward-time
            ! unitary rotation -- exact and unconditionally safe.
            !-----------------------------------------------------------------
            do isub = 1, 3
                call build_HVG(nba, eigen_active, p_active, Ac_node(:, 1, isub), H1)
                call build_HVG(nba, eigen_active, p_active, Ac_node(:, 2, isub), H2)
                if (sbe%flag_coulomb) then
                    H1 = H1 + sbe%sigma_hf(:, :, ik)
                    H2 = H2 + sbe%sigma_hf(:, :, ik)
                end if
                call cf4_unitary_step(nba, rho_a, H1, H2, tau_sub(isub))
            end do

            !-----------------------------------------------------------------
            ! Step 3: D(h/2) -- Strang dissipative half-step (see Step 1)
            !-----------------------------------------------------------------
            if (sbe%flag_decoh .or. sbe%flag_impact) then
                call build_HVG(nba, eigen_active, p_active, Ac_end, HVG)
                if (sbe%flag_coulomb) HVG = HVG + sbe%sigma_hf(:, :, ik)
                call houston_dissipate(sbe, nba, rho_a, HVG, p_active, Ac_end, X_a, &
                                       0.5d0 * dt, V_end, w_act_sub)
            else
                V_end = 0d0
            end if

            ! Branch-position update via the midpoint (average of endpoint)
            ! velocities -- consistent with the overall 4th-order accuracy of
            ! CF4 (a forward-Euler X_a += V_a(t_start)*h would degrade the
            ! scheme to 1st order in the branch coordinates).
            do i = 1, nba
                sbe%X_branch(sbe%active_idx(i), ik) = X_a(i) + 0.5d0 * (V_begin(i) + V_end(i)) * dt
            end do

            ! Scatter the evolved active block back into the full matrix
            do j = 1, nba; do i = 1, nba
                in = sbe%active_idx(i); im = sbe%active_idx(j)
                rho_n_full(in, im) = rho_a(i, j)
            end do; end do
            sbe%rho(:, :, ik) = rho_n_full(:, :)
        end if

        ! Hermiticity (numerical safeguard)
        do j = 1, nb; do i = 1, nb
            sbe%rho(i, j, ik) = 0.5d0 * (sbe%rho(i, j, ik) + conjg(sbe%rho(j, i, ik)))
        end do; end do

        ! Freeze deep core/high-energy zones
        do j = 1, nb; do i = 1, nb
            if (.not. (sbe%is_active(i) .and. sbe%is_active(j))) then
                if (i == j) then
                    if (gs%occup(i, ik) > 0.5d0) then
                        ! Ground-state occupation: 2 (scalar bands) or 1 (spinor bands)
                        sbe%rho(i, j, ik) = dcmplx(gs%occup(i, ik), 0.0d0)
                    else
                        sbe%rho(i, j, ik) = dcmplx(0.0d0, 0.0d0)
                    end if
                else
                    sbe%rho(i, j, ik) = dcmplx(0.0d0, 0.0d0)
                end if
            end if
        end do; end do

    end do
    !$omp end do

    if (nba > 0) then
        deallocate(p_active, rho_a, H1, H2, HVG, eigen_active, V_begin, V_end, X_a, w_act_sub)
    end if
    deallocate(p_k_full, rho_n_full)
    !$omp end parallel

end subroutine dt_evolve_bloch_cf4


! Build the instantaneous velocity-gauge Hamiltonian in the active subspace:
!   H_VG(t) = diag(eigen) + A(t) . pi
subroutine build_HVG(nba, eigen_active, p_active, Ac, H)
    implicit none
    integer,    intent(in)  :: nba
    real(8),    intent(in)  :: eigen_active(nba)
    complex(8), intent(in)  :: p_active(nba, nba, 3)
    real(8),    intent(in)  :: Ac(3)
    complex(8), intent(out) :: H(nba, nba)
    integer :: i, idir

    H = Ac(1) * p_active(:, :, 1) + Ac(2) * p_active(:, :, 2) + Ac(3) * p_active(:, :, 3)
    do i = 1, nba
        H(i, i) = H(i, i) + eigen_active(i)
    end do
end subroutine build_HVG


!=============================================================================
! Time-dependent Hartree-Fock (exchange/Fock) Coulomb self-energy, in the
! velocity-gauge stationary-Bloch active basis [Golde-Kira-Meier-Koch, Phys.
! Status Solidi B 248, 863 (2011), Eqs. 4-5]:
!
!   Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q),
!   V(p) = strength * 4 pi / ( eps Omega_cell Nk (|p|^2 + kappa^2) )   [a.u.]
!
! This is a NON-k-local mean field: it convolves the active-band density
! matrix with the screened Coulomb kernel over the whole BZ. The convolution
! is gauge-covariant under the uniform velocity-gauge shift k -> k - A(t)
! (A cancels in k - q), so it is evaluated directly on the grid-k rho. Adding
! Sigma to H_VG renormalizes the diagonal energies (eps~ = eps - sum_q V f_q)
! and the off-diagonal Rabi term (Omega = d.E + sum_q V p_q); the (1-f_e-f_h)
! Pauli factor then follows from the von Neumann commutator [H+Sigma, rho].
! Sigma is Hermitian (rho Hermitian, V real), so the propagation stays unitary.
!
! Equilibrium subtraction: the EPM/DFT bands carry NO explicit Coulomb
! exchange, so Sigma is built from the DEVIATION from the ground state,
! drho = rho - rho_0 (rho_0 = diag(gs%occup)). All of drho's pieces vanish at
! t=0 (no excited electrons/holes, no polarization), so Sigma(t=0)=0 and the
! equilibrium gap stays the EPM gap; only the dynamical (carrier-induced)
! renormalization is added -- exactly the f^e, f^h, p of Eqs. 4-5.
!
! MPI: rho is k-partitioned (ik_min:ik_max). We zero-pad the local block and
! comm_summation to obtain the full-BZ drho on every rank (an all-gather), then
! each rank forms Sigma for its local k. Cost O(Nk^2 nba^2) per call; evaluated
! once per dt step (the mean field is frozen over the Strang/CF4 sub-steps).
!=============================================================================
subroutine compute_coulomb_selfenergy(sbe, gs)
    use communication, only: comm_summation
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs_info),      intent(in)    :: gs
    integer :: nba, ik, iq, i, j, in, im
    real(8) :: dr(3), dkx, dky, dkz, q2, vkq
    complex(8), allocatable :: rho_loc(:, :, :), rho_all(:, :, :)

    nba = sbe%n_active_bands
    if (nba < 1) return

    allocate(rho_loc(nba, nba, sbe%nk), rho_all(nba, nba, sbe%nk))
    rho_loc = (0d0, 0d0)
    do ik = sbe%ik_min, sbe%ik_max
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                rho_loc(i, j, ik) = sbe%rho(in, im, ik)
            end do
            ! subtract the ground-state occupation (rho_0) on the diagonal
            rho_loc(j, j, ik) = rho_loc(j, j, ik) - dcmplx(gs%occup(im, ik), 0d0)
        end do
    end do
    call comm_summation(rho_loc, rho_all, nba * nba * sbe%nk, sbe%icomm)

    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp    private(ik, iq, i, j, dr, dkx, dky, dkz, q2, vkq)
    do ik = sbe%ik_min, sbe%ik_max
        sbe%sigma_hf(:, :, ik) = (0d0, 0d0)
        do iq = 1, sbe%nk
            ! minimum image of (k - q) in reduced coords, then -> Cartesian
            dr(1) = gs%kpoint(1, ik) - gs%kpoint(1, iq)
            dr(2) = gs%kpoint(2, ik) - gs%kpoint(2, iq)
            dr(3) = gs%kpoint(3, ik) - gs%kpoint(3, iq)
            dr(1) = dr(1) - dnint(dr(1))
            dr(2) = dr(2) - dnint(dr(2))
            dr(3) = dr(3) - dnint(dr(3))
            dkx = dr(1)*gs%b_matrix(1,1) + dr(2)*gs%b_matrix(2,1) + dr(3)*gs%b_matrix(3,1)
            dky = dr(1)*gs%b_matrix(1,2) + dr(2)*gs%b_matrix(2,2) + dr(3)*gs%b_matrix(3,2)
            dkz = dr(1)*gs%b_matrix(1,3) + dr(2)*gs%b_matrix(2,3) + dr(3)*gs%b_matrix(3,3)
            q2 = dkx*dkx + dky*dky + dkz*dkz
            if (q2 < 1d-12) cycle          ! exclude q = k (the V(0) self term)
            vkq = sbe%coul_pref / (q2 + sbe%coul_screen2)
            do j = 1, nba
                do i = 1, nba
                    sbe%sigma_hf(i, j, ik) = sbe%sigma_hf(i, j, ik) - vkq * rho_all(i, j, iq)
                end do
            end do
        end do
    end do
    !$omp end parallel do

    deallocate(rho_loc, rho_all)
end subroutine compute_coulomb_selfenergy


! Single CF4 (commutator-free Magnus, 4th order) sub-step of length tau,
! evaluated on the two Gauss-Legendre Hamiltonians H1=H(t+c1*tau), H2=H(t+c2*tau):
!   Omega1 = tau (alpha1 H1 + alpha2 H2),  Omega2 = tau (alpha2 H1 + alpha1 H2)
!   rho <- exp(-i Omega2) exp(-i Omega1) rho exp(+i Omega1) exp(+i Omega2)
! Implemented as two successive exact unitary rotations (each built from an
! eigendecomposition of the Hermitian generator, so no Pade/Krylov truncation
! error is introduced -- the propagator is exactly unitary to machine precision).
subroutine cf4_unitary_step(nba, rho, H1, H2, tau)
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: H1(nba, nba), H2(nba, nba)
    real(8),    intent(in)    :: tau
    complex(8) :: Omega(nba, nba)

    Omega = tau * (cf4_alpha1 * H1 + cf4_alpha2 * H2)
    call apply_unitary_rotation(nba, rho, Omega)

    Omega = tau * (cf4_alpha2 * H1 + cf4_alpha1 * H2)
    call apply_unitary_rotation(nba, rho, Omega)
end subroutine cf4_unitary_step


! Apply rho -> U rho U^dagger with U = exp(-i*Omega) for Hermitian Omega,
! computed exactly via eigendecomposition Omega = W diag(lambda) W^dagger:
!   U rho U^dagger = W [ exp(-i lambda_i) (W^dagger rho W)_ij exp(+i lambda_j) ] W^dagger
subroutine apply_unitary_rotation(nba, rho, Omega)
    use eigen_lapack, only: eigen_zheev
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: Omega(nba, nba)

    real(8)    :: evals(nba)
    complex(8) :: W(nba, nba), t1(nba, nba), t2(nba, nba)
    integer :: i, j

    call eigen_zheev(Omega, evals, W)

    call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, rho, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,   nba, dcmplx(0d0, 0d0), t2, nba)

    do j = 1, nba
        do i = 1, nba
            t2(i, j) = t2(i, j) * exp(dcmplx(0d0, -(evals(i) - evals(j))))
        end do
    end do

    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, t2, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'C', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,  nba, dcmplx(0d0, 0d0), rho, nba)
end subroutine apply_unitary_rotation


! Strang dissipative half-step in the instantaneous Houston (adiabatic) basis.
! One ZHEEV serves both dissipative channels:
!   1) diagonalize the instantaneous H_VG(t) -> Houston basis U, {eps_a}
!   2) rotate rho~ = U^dagger rho U
!   3) [flag_decoh ] rho~_ab <- exp[-lambda (X_a - X_b)^2 * tau] * rho~_ab
!      (Hadamard product with a PSD Gram/RBF matrix => exactly CPTP, by
!      Schoenberg/Bochner positive-definiteness + the Schur product theorem)
!   4) [flag_impact] k-local impact-ionization jump channels (see
!      apply_impact_ionization below) -- each an exact amplitude-damping map,
!      CPTP for any tau >= 0
!   5) rotate back rho = U rho~ U^dagger
! Also returns the instantaneous branch (group) velocities in the field
! polarization direction, V_a = [(U^dagger pi U)_aa . e_hat] + (A . e_hat)
! (computed only when the Kuhn-Zurek dephasing needs the branch positions).
subroutine houston_dissipate(sbe, nba, rho, H, p_active, Ac, X, tau, V, w_act_sub)
    use eigen_lapack, only: eigen_zheev
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    complex(8), intent(in)    :: H(nba, nba)
    complex(8), intent(in)    :: p_active(nba, nba, 3)
    real(8),    intent(in)    :: Ac(3)
    real(8),    intent(in)    :: X(nba)
    real(8),    intent(in)    :: tau
    real(8),    intent(out)   :: V(nba)
    real(8),    intent(in)    :: w_act_sub(4, nba)  ! field-free sublattice weights of the active bands

    real(8)    :: evals(nba), ehat(3), Ac_norm
    complex(8) :: W(nba, nba), t1(nba, nba), t2(nba, nba)
    real(8)    :: wsub_branch(4, nba), wcoef
    integer :: i, j, idir, a, s

    call eigen_zheev(H, evals, W)

    ! rho~ = U^dagger rho U
    call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, rho, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,   nba, dcmplx(0d0, 0d0), t2, nba)

    ! Exact Hadamard/Gaussian dephasing kernel (PSD for tau >= 0)
    if (sbe%flag_decoh) then
        do j = 1, nba
            do i = 1, nba
                t2(i, j) = t2(i, j) * exp(-sbe%lambda_decoh * (X(i) - X(j))**2 * tau)
            end do
        end do
    end if

    ! k-local impact-ionization channels (threshold-gated, frozen rates)
    if (sbe%flag_impact) then
        ! Project the field-free sublattice weights of the active Bloch bands
        ! onto the instantaneous Houston branches a:
        !   w_s(a) = sum_i |U_ia|^2 w_s(band_i),   sum_s w_s(a) = 1.
        ! This makes the sublattice character field-aware (correct under the
        ! Houston rotation U) and degenerates to the bare band weights at A=0.
        if (sbe%flag_unfold_ii) then
            wsub_branch = 0d0
            do a = 1, nba
                do i = 1, nba
                    wcoef = real(W(i, a))**2 + aimag(W(i, a))**2
                    do s = 1, 4
                        wsub_branch(s, a) = wsub_branch(s, a) + wcoef * w_act_sub(s, i)
                    end do
                end do
            end do
        else
            wsub_branch = 0d0
        end if
        call apply_impact_ionization(sbe, nba, t2, evals, Ac, tau, wsub_branch, sbe%flag_unfold_ii)
    end if

    ! rho = U rho~ U^dagger
    call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, t2, nba, dcmplx(0d0, 0d0), t1, nba)
    call ZGEMM('N', 'C', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W,  nba, dcmplx(0d0, 0d0), rho, nba)

    ! Branch velocities, projected on the polarization direction of A(t);
    ! only the Kuhn-Zurek branch positions consume them.
    V = 0d0
    if (sbe%flag_decoh) then
        Ac_norm = sqrt(dot_product(Ac, Ac))
        if (Ac_norm > 1.0d-12) then
            ehat = Ac / Ac_norm
        else
            ehat = (/ 1d0, 0d0, 0d0 /)
        end if

        do idir = 1, 3
            call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0, 0d0), W,  nba, p_active(:, :, idir), nba, &
                       dcmplx(0d0, 0d0), t1, nba)
            call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0, 0d0), t1, nba, W, nba, dcmplx(0d0, 0d0), t2, nba)
            do i = 1, nba
                V(i) = V(i) + ehat(idir) * (real(t2(i, i)) + Ac(idir))
            end do
        end do
    end if
end subroutine houston_dissipate


!=============================================================================
! k-local impact ionization (Stobbe-Redmer-Schattke rate fit, PRB 49, 4494)
! in the Houston basis. The quartic two-particle event
!   A_h = sqrt(gamma_St) c+_h' c+_c1 c_v1 c_h
! is closed k-locally (no momentum transfer) and Hartree-Fock-factorized
! (two-particle closure, Rosati-Iotti-Dolcini-Rossi PRB 90, 125140) into two
! effective single-particle amplitude-damping channels with FROZEN scalar
! rates (partner population and Pauli blockers enter as factors clamped to
! [0,1], guaranteeing Gamma >= 0 => each map is exactly CPTP):
!   rel : branch h  -> h'  with Gamma_rel  = gamma_St * f_v1 (1-f_c1)(1-f_h')
!   pair: branch v1 -> c1  with Gamma_pair = gamma_St * f_h  (1-f_c1)(1-f_h')
! "Cold pair": secondaries are born at the band-edge branches v1/c1 with no
! kinetic energy; the primary drops to the conduction branch closest to
! eps_h - E_g. The hot set is gated on the kinetic energy from the global CBM,
!   eps_kin = eps_h + |A|^2/2 - E_CBM
! (the A^2/2 scalar dropped in H_VG is restored here, where the comparison
! against the field-free constant E_CBM requires it; by the Houston identity
! this equals E_h(k+A) - E_CBM -- the scale on which the Stobbe fit is
! defined). For most k-points and times the gate is empty and the channel
! costs O(N_C) comparisons -- the "rare impact events" mechanism.
!
! SUBLATTICE (BAND-UNFOLDING) RESOLUTION (use_unfold = .true.)
! ------------------------------------------------------------
! Under cubic folding, one grid k hides 4 distinct primitive crystal
! momenta -- the 4 FCC sublattices. Genuine impact ionization conserves
! the PRIMITIVE crystal momentum, so the strict k-local event must close
! WITHIN a single sublattice: the primary (h -> h') and the secondary pair
! (v1 -> c1) must all sit in the same primitive-cell sector. The folded
! treatment above ignores this and lets a primary on sublattice s create
! its secondary pair at the GLOBAL band edges v1/c1 (which may belong to a
! different sublattice) -- a momentum-non-conserving "false generation"
! event that inflates the carrier multiplication rate.
!
! With the spectral weights w_s(a) (sum_s w_s(a) = 1, projected onto the
! Houston branches) we resolve the event per sublattice s = 1..4:
!   * the secondary pair is created at the s-RESOLVED band edges
!     v1(s), c1(s) (the valence/conduction branches whose weight is
!     dominant on s, taken at the gap edge),
!   * the primary's receiving branch h'(s) is the conduction branch with
!     weight on s closest to eps_h - E_g,
!   * both channel rates are scaled by w_s(h) -- the primary's fraction on
!     s -- so that summing over s reproduces the primary's total
!     relaxation rate. In the no-folding limit (one sublattice, w = 1)
!     v1(s)=v1, c1(s)=c1, h'(s)=h' and the scheme reduces EXACTLY to the
!     folded channel above.
! Every sub-channel is still an exact CPTP amplitude-damping map (rates
! >= 0 by construction), so the composition remains CPTP.
!=============================================================================
subroutine apply_impact_ionization(sbe, nba, rho_ad, evals, Ac, tau, wsub, use_unfold)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho_ad(nba, nba)   ! adiabatic-basis rho~
    real(8),    intent(in)    :: evals(nba)         ! Houston branch energies
    real(8),    intent(in)    :: Ac(3)
    real(8),    intent(in)    :: tau
    real(8),    intent(in)    :: wsub(4, nba)       ! per-branch sublattice weights (sum_s = 1)
    logical,    intent(in)    :: use_unfold

    integer :: iv1, ic1, ih, ihp, a, s
    integer :: iv1_sub(4), ic1_sub(4), ihp_s
    real(8) :: a2half, ekin, d, gam, etgt, f, ws, wbest
    real(8) :: pv1, phh, bc1, bhp, g_rel, g_pair
    real(8), parameter :: occ_eps = 1d-12
    real(8), parameter :: w_min   = 0.05d0   ! ignore sublattices the primary barely touches
    real(8), parameter :: w_tie   = 1d-12

    if (sbe%nv_act < 1 .or. sbe%nv_act >= nba) return
    iv1 = sbe%nv_act        ! topmost valence branch (energy-ordered)
    ic1 = sbe%nv_act + 1    ! lowest conduction branch
    f = sbe%occ_max
    a2half = 0.5d0 * dot_product(Ac, Ac)

    !-------------------------------------------------------------------------
    ! Folded (single-pool) treatment: no unfolding weights available
    !-------------------------------------------------------------------------
    if (.not. use_unfold) then
        do ih = ic1, nba
            ekin = evals(ih) + a2half - sbe%ii_ecbm_au
            d = ekin - sbe%ii_eth_au
            if (d <= 0d0) cycle
            if (real(rho_ad(ih, ih)) < occ_eps) cycle
            gam = sbe%ii_pref_au * d**sbe%ii_exponent
            if (sbe%ii_ramp_au > 0d0 .and. d < sbe%ii_ramp_au) gam = gam * d / sbe%ii_ramp_au
            etgt = evals(ih) - sbe%ii_eg_au
            ihp = ic1
            do a = ic1, nba
                if (abs(evals(a) - etgt) < abs(evals(ihp) - etgt)) ihp = a
            end do
            if (ihp == ih) cycle
            pv1 = min(max(real(rho_ad(iv1, iv1)) / f, 0d0), 1d0)
            phh = min(max(real(rho_ad(ih,  ih )) / f, 0d0), 1d0)
            bc1 = min(max(1d0 - real(rho_ad(ic1, ic1)) / f, 0d0), 1d0)
            bhp = min(max(1d0 - real(rho_ad(ihp, ihp)) / f, 0d0), 1d0)
            g_rel  = gam * pv1 * bc1 * bhp
            g_pair = gam * phh * bc1 * bhp
            call apply_damping_channel(nba, rho_ad, ih,  ihp, g_rel,  tau)
            call apply_damping_channel(nba, rho_ad, iv1, ic1, g_pair, tau)
        end do
        return
    end if

    !-------------------------------------------------------------------------
    ! Sublattice-resolved (primitive-momentum-conserving) treatment.
    ! Precompute the s-resolved band edges once per k:
    !   v1(s) = valence branch dominant on s, taken at the TOP of the valence
    !   c1(s) = conduction branch dominant on s, taken at the BOTTOM of the
    !           conduction. The energy tie-break recovers the global edges
    !           (nv_act, nv_act+1) when all weights are equal (no folding).
    !-------------------------------------------------------------------------
    do s = 1, 4
        iv1_sub(s) = iv1
        wbest = -1d0
        do a = 1, sbe%nv_act
            if (wsub(s, a) > wbest + w_tie .or. &
                (abs(wsub(s, a) - wbest) <= w_tie .and. evals(a) > evals(iv1_sub(s)))) then
                wbest = wsub(s, a); iv1_sub(s) = a
            end if
        end do
        ic1_sub(s) = ic1
        wbest = -1d0
        do a = ic1, nba
            if (wsub(s, a) > wbest + w_tie .or. &
                (abs(wsub(s, a) - wbest) <= w_tie .and. evals(a) < evals(ic1_sub(s)))) then
                wbest = wsub(s, a); ic1_sub(s) = a
            end if
        end do
    end do

    do ih = ic1, nba
        ! Threshold gate on the kinetic energy from the field-free CBM
        ekin = evals(ih) + a2half - sbe%ii_ecbm_au
        d = ekin - sbe%ii_eth_au
        if (d <= 0d0) cycle
        if (real(rho_ad(ih, ih)) < occ_eps) cycle

        ! Stobbe rate, with optional linear ramp over the fit resolution
        gam = sbe%ii_pref_au * d**sbe%ii_exponent
        if (sbe%ii_ramp_au > 0d0 .and. d < sbe%ii_ramp_au) gam = gam * d / sbe%ii_ramp_au

        etgt = evals(ih) - sbe%ii_eg_au

        ! Distribute the event across the sublattices the primary occupies
        do s = 1, 4
            ws = wsub(s, ih)
            if (ws < w_min) cycle

            ! Receiving branch of the primary on sublattice s: conduction
            ! branch with weight on s closest to eps_h - E_g
            ihp_s = ic1_sub(s)
            do a = ic1, nba
                if (wsub(s, a) < w_min) cycle
                if (abs(evals(a) - etgt) < abs(evals(ihp_s) - etgt)) ihp_s = a
            end do
            if (ihp_s == ih) cycle   ! no receiving state below: event suppressed

            ! Frozen scalar factors at the s-resolved branches, clamped to
            ! [0,1]; the leading w_s(h) splits the primary across sublattices
            ! (sum_s w_s(h) = 1 => total primary relaxation preserved)
            pv1 = min(max(real(rho_ad(iv1_sub(s), iv1_sub(s))) / f, 0d0), 1d0)
            phh = min(max(real(rho_ad(ih,         ih        )) / f, 0d0), 1d0)
            bc1 = min(max(1d0 - real(rho_ad(ic1_sub(s), ic1_sub(s))) / f, 0d0), 1d0)
            bhp = min(max(1d0 - real(rho_ad(ihp_s,      ihp_s     )) / f, 0d0), 1d0)
            g_rel  = gam * ws * pv1 * bc1 * bhp
            g_pair = gam * ws * phh * bc1 * bhp

            ! Sequential exact amplitude-damping maps (each CPTP for tau>=0)
            call apply_damping_channel(nba, rho_ad, ih,         ihp_s,      g_rel,  tau)
            call apply_damping_channel(nba, rho_ad, iv1_sub(s), ic1_sub(s), g_pair, tau)
        end do
    end do
end subroutine apply_impact_ionization


! Exact finite-time map of a single Lindblad jump channel L = sqrt(Gamma)|f><i|
! (amplitude damping i -> f) in the basis where i, f are basis states; O(N):
!   rho_ii -> e^{-G tau} rho_ii
!   rho_ff -> rho_ff + (1 - e^{-G tau}) rho_ii
!   rho_ib -> e^{-G tau/2} rho_ib  (and rho_bi), for all b /= i
! Trace-preserving and completely positive for any Gamma, tau >= 0. The jump
! L rho L^dagger = Gamma rho_ii |f><f| feeds only the f-diagonal: no new
! coherences are created, while the populations and coherences of the source
! branch are damped -- ionization is itself a decoherence channel.
subroutine apply_damping_channel(nba, rho, i_src, i_dst, gamma_ch, tau)
    implicit none
    integer,    intent(in)    :: nba
    complex(8), intent(inout) :: rho(nba, nba)
    integer,    intent(in)    :: i_src, i_dst
    real(8),    intent(in)    :: gamma_ch, tau

    real(8) :: g, gh, transfer
    integer :: b

    if (gamma_ch * tau < 1d-14) return
    g  = exp(-gamma_ch * tau)
    gh = sqrt(g)
    transfer = (1d0 - g) * real(rho(i_src, i_src))

    do b = 1, nba
        if (b == i_src) cycle
        rho(i_src, b) = gh * rho(i_src, b)
        rho(b, i_src) = gh * rho(b, i_src)
    end do
    rho(i_src, i_src) = g * rho(i_src, i_src)
    rho(i_dst, i_dst) = rho(i_dst, i_dst) + transfer
end subroutine apply_damping_channel


! Population of band `ib_target` resolved per k-point, in the stationary
! Bloch (crystal-gauge) eigenbasis.
!
! In the Velocity Gauge (VG) the SBE propagates rho(k,t) at the fixed grid
! crystal momentum k, while the physical electrons are displaced to
!   k'(t) = k - A(t)    (A in a.u., e/hbar = 1)
! The crystal-gauge population at k requires projecting onto the eigenstates
! of H_0(k'), the field-free Hamiltonian evaluated at the shifted momentum.
! To first order in A(t) (valid when |A(t)| << BZ size):
!   H_0(k - A(t)) ≈ H_0(k) + (k'-k)·∂H_0/∂k = H_0(k) - A(t)·p(k)
!                 = diag(eigen) - A·p   (note: MINUS sign, unlike Houston + sign)
! We diagonalise H_crystal to get U_shifted, apply a greedy bipartite match
! on |U_shifted_ij| to correct the energy-sort ambiguity of ZHEEV at near-
! degeneracies, then form rho_crystal = U_sorted^dagger rho_VG U_sorted and
! return Re(rho_crystal[ia_target, ia_target]).
subroutine calc_bloch_population_k(sbe, gs, Ac, ib_target, pop_k, icomm)
    use eigen_lapack, only: eigen_zheev
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    integer,                  intent(in)  :: ib_target
    real(8),                  intent(out) :: pop_k(1:sbe%nk)
    integer,                  intent(in)  :: icomm

    integer :: nba, ia_target, ik, i, j, in, im, n_done, best_i, best_j
    real(8)  :: curr_max
    integer,  allocatable :: zone_map(:)
    logical,  allocatable :: row_used(:), col_used(:)
    real(8),    allocatable :: pop_local(:), evals(:), p_k_full(:,:,:), eigen_a(:)
    complex(8), allocatable :: H(:,:), W(:,:), W_sorted(:,:), t1(:,:), t2(:,:), rho_a(:,:)

    nba = sbe%n_active_bands

    ! Locate ib_target in the active subspace; return all-zero if it is frozen.
    ia_target = 0
    do i = 1, nba
        if (sbe%active_idx(i) == ib_target) then
            ia_target = i
            exit
        end if
    end do
    if (ia_target == 0) then
        pop_k = 0d0
        return
    end if

    allocate(pop_local(1:sbe%nk))
    allocate(evals(nba), H(nba,nba), W(nba,nba), W_sorted(nba,nba), t1(nba,nba), t2(nba,nba))
    allocate(rho_a(nba,nba), p_k_full(sbe%nb, sbe%nb, 3), eigen_a(nba))
    allocate(zone_map(nba), row_used(nba), col_used(nba))
    pop_local = 0d0

    do ik = sbe%ik_min, sbe%ik_max
        ! Build the full p matrix exactly as the propagator does.
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)

        ! Restrict H_VG and rho to the active subspace (mirrors dt_evolve_bloch_cf4).
        do i = 1, nba
            eigen_a(i) = gs%eigen(sbe%active_idx(i), ik)
        end do
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                ! H_0(k') = H_0(k) - A·p  (crystal-gauge shift k' = k - A)
                H(i, j) = -Ac(1)*p_k_full(in,im,1) &
                          - Ac(2)*p_k_full(in,im,2) &
                          - Ac(3)*p_k_full(in,im,3)
                rho_a(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        do i = 1, nba
            H(i, i) = H(i, i) + eigen_a(i)
        end do

        ! Diagonalize H_crystal = H_0(k-A): H = W Lambda W^dagger  (LAPACK ZHEEV)
        call eigen_zheev(H, evals, W)

        ! Overlap-tracking permutation (greedy bipartite match on |W_ij|).
        ! H_0(k) = diag(eigen) => reference U_0 = I => overlap S = W directly.
        row_used = .false.
        col_used = .false.
        zone_map  = 0
        do n_done = 1, nba
            curr_max = -1d0
            best_i = 1; best_j = 1
            do j = 1, nba
                if (col_used(j)) cycle
                do i = 1, nba
                    if (row_used(i)) cycle
                    if (abs(W(i, j)) > curr_max) then
                        curr_max = abs(W(i, j))
                        best_i   = i
                        best_j   = j
                    end if
                end do
            end do
            zone_map(best_j) = best_i
            row_used(best_i) = .true.
            col_used(best_j) = .true.
        end do

        ! Reorder columns: W_sorted(:, zone_map(j)) = W(:, j)
        do j = 1, nba
            W_sorted(:, zone_map(j)) = W(:, j)
        end do

        ! rho_crystal = U_sorted^dagger rho_VG U_sorted
        call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0,0d0), W_sorted, nba, rho_a,    nba, dcmplx(0d0,0d0), t1, nba)
        call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0,0d0), t1,       nba, W_sorted, nba, dcmplx(0d0,0d0), t2, nba)
        pop_local(ik) = real(t2(ia_target, ia_target))
    end do

    call comm_summation(pop_local, pop_k, sbe%nk, icomm)
    deallocate(pop_local, evals, H, W, W_sorted, t1, t2, rho_a, p_k_full, eigen_a)
    deallocate(zone_map, row_used, col_used)
end subroutine calc_bloch_population_k


! Population of the PHYSICAL lowest conduction band (CB1) of each folded
! primitive BZ point, resolved per cubic k-point.
!
! In the cubic-supercell EPM dataset every band at k_sc belongs to exactly
! one FCC sublattice s (= primitive point k_prim = k_sc + G0(s); the folding
! is exact, see the unfold map written by `epm_gaas_reference.py unfoldmap`).
! Energy-ordered supercell branch indices therefore mix DIFFERENT physical
! primitive bands from k to k. This routine projects rho into the crystal
! gauge exactly like calc_bloch_population_k and then accumulates, for each
! sublattice s, the population of its physical CB1 level:
!   spinor input: primitive ranks nv_prim+1 and nv_prim+2 (Kramers pair,
!                 spins summed); scalar input: rank nv_prim+1.
! pop_lev(L, s, ik) is the crystal-gauge population of the L-th physical primitive
! level (spins summed) at the folded primitive point k_prim = k_sc(ik) + G0(s):
!   L = 1 -> VB-1 (second valence from top), 2 -> VB (top valence),
!   L = 3 -> CB1 (lowest conduction),         4 -> CB2 (second conduction).
! Band assignment comes entirely from the EPM `unfoldmap` (gs%unfold_sub/prim).
! Energies for the spectral plot are taken from *_bandpath.data, so none are
! emitted here.
subroutine calc_unfolded_population_k(sbe, gs, Ac, pop_lev, icomm)
    use eigen_lapack, only: eigen_zheev
    use salmon_global, only: yn_sbe_spinor
    implicit none
    type(s_sbe_bloch_solver), intent(in)  :: sbe
    type(s_sbe_gs_info),      intent(in)  :: gs
    real(8),                  intent(in)  :: Ac(3)
    real(8),                  intent(out) :: pop_lev(1:4, 1:4, 1:sbe%nk)
    integer,                  intent(in)  :: icomm

    integer :: nba, ik, i, j, in, im, n_done, best_i, best_j
    integer :: isub, irank_prim, n_spin, nv_phys, pphys, off, islot, s
    real(8)  :: curr_max, popi
    integer,  allocatable :: zone_map(:)
    logical,  allocatable :: row_used(:), col_used(:)
    real(8),    allocatable :: pop_local(:, :, :), evals(:), p_k_full(:,:,:), eigen_a(:)
    complex(8), allocatable :: H(:,:), W(:,:), W_sorted(:,:), t1(:,:), t2(:,:), rho_a(:,:)

    nba = sbe%n_active_bands
    pop_lev = 0d0
    if (.not. gs%have_unfold .or. nba == 0) return

    ! Spin states per physical level in the primitive rank ordering, and the
    ! number of physical valence bands (Kramers doublets for spinor input).
    n_spin  = merge(2, 1, yn_sbe_spinor == 'y')
    nv_phys = gs%nv_prim / n_spin

    allocate(pop_local(1:4, 1:4, 1:sbe%nk))
    allocate(evals(nba), H(nba,nba), W(nba,nba), W_sorted(nba,nba), t1(nba,nba), t2(nba,nba))
    allocate(rho_a(nba,nba), p_k_full(sbe%nb, sbe%nb, 3), eigen_a(nba))
    allocate(zone_map(nba), row_used(nba), col_used(nba))
    pop_local = 0d0

    do ik = sbe%ik_min, sbe%ik_max
        ! Crystal-gauge projection: identical to calc_bloch_population_k
        p_k_full(:, :, :) = gs%p_tm_matrix(:, :, :, ik)
        if (sbe%flag_vnl_correction) &
            p_k_full(:, :, :) = p_k_full(:, :, :) + gs%rvnl_tm_matrix(:, :, :, ik)

        do i = 1, nba
            eigen_a(i) = gs%eigen(sbe%active_idx(i), ik)
        end do
        do j = 1, nba
            im = sbe%active_idx(j)
            do i = 1, nba
                in = sbe%active_idx(i)
                H(i, j) = -Ac(1)*p_k_full(in,im,1) &
                          - Ac(2)*p_k_full(in,im,2) &
                          - Ac(3)*p_k_full(in,im,3)
                rho_a(i, j) = sbe%rho(in, im, ik)
            end do
        end do
        do i = 1, nba
            H(i, i) = H(i, i) + eigen_a(i)
        end do

        call eigen_zheev(H, evals, W)

        row_used = .false.
        col_used = .false.
        zone_map  = 0
        do n_done = 1, nba
            curr_max = -1d0
            best_i = 1; best_j = 1
            do j = 1, nba
                if (col_used(j)) cycle
                do i = 1, nba
                    if (row_used(i)) cycle
                    if (abs(W(i, j)) > curr_max) then
                        curr_max = abs(W(i, j))
                        best_i   = i
                        best_j   = j
                    end if
                end do
            end do
            zone_map(best_j) = best_i
            row_used(best_i) = .true.
            col_used(best_j) = .true.
        end do
        do j = 1, nba
            W_sorted(:, zone_map(j)) = W(:, j)
        end do

        call ZGEMM('C', 'N', nba, nba, nba, dcmplx(1d0,0d0), W_sorted, nba, rho_a,    nba, dcmplx(0d0,0d0), t1, nba)
        call ZGEMM('N', 'N', nba, nba, nba, dcmplx(1d0,0d0), t1,       nba, W_sorted, nba, dcmplx(0d0,0d0), t2, nba)

        ! Accumulate the spin-summed population of the top two valence and the
        ! bottom two conduction physical primitive bands of every sublattice.
        ! The cubic band's primitive rank (from the unfold map) fixes the physical
        ! slot {VB-1, VB, CB1, CB2}; its population is then DISTRIBUTED over the 4
        ! sublattices by the spectral weights gs%unfold_w (sum_s = 1). At a
        ! symmetry degeneracy this splits the population equally among the
        ! equivalent primitive points, instead of an argmax dumping it on one.
        do i = 1, nba
            in = sbe%active_idx(i)
            isub = gs%unfold_sub(in, ik)
            if (isub < 1 .or. isub > 4) cycle
            irank_prim = gs%unfold_prim(in, ik)
            if (irank_prim < 1) cycle
            if (irank_prim <= gs%nv_prim) then
                ! Valence: physical index counted from the band bottom; offset
                ! from the top valence band (off = 0 top, -1 next).
                pphys = (irank_prim + n_spin - 1) / n_spin
                off   = pphys - nv_phys
                if (off == 0) then
                    islot = 2          ! VB (top valence)
                else if (off == -1) then
                    islot = 1          ! VB-1
                else
                    cycle
                end if
            else
                ! Conduction: physical index counted from the band bottom.
                pphys = (irank_prim - gs%nv_prim + n_spin - 1) / n_spin
                if (pphys == 1) then
                    islot = 3          ! CB1
                else if (pphys == 2) then
                    islot = 4          ! CB2
                else
                    cycle
                end if
            end if
            popi = real(t2(i, i))
            do s = 1, 4
                pop_local(islot, s, ik) = pop_local(islot, s, ik) &
                    & + gs%unfold_w(s, in, ik) * popi
            end do
        end do
    end do

    call comm_summation(pop_local, pop_lev, 16 * sbe%nk, icomm)
    deallocate(pop_local, evals, H, W, W_sorted, t1, t2, rho_a, p_k_full, eigen_a)
    deallocate(zone_map, row_used, col_used)
end subroutine calc_unfolded_population_k


end module
