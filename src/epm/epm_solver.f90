!
!  Local Empirical Pseudopotential Method (EPM) ground-state solver
!
!  Builds the zincblende/diamond (GaAs / Si / Si_cb) band structure from the
!  local Cohen-Bergstresser pseudopotential in a plane-wave basis, diagonalizes
!  H(k) at each k-point of the requested Monkhorst-Pack grid, and prepares the
!  data structures that are written (by main_epm) into SYSNAME_k.data /
!  SYSNAME_eigen.data / SYSNAME_tm.data -- i.e. exactly the files that
!  gs_info_ssbe::init_sbe_gs_info reads to start an SBE real-time calculation.
!  This closes the EPM -> SBE chain without any external pre-/post-processing.
!
!  CONVENTION (must match the Python reference epm_gaas_reference.py, which is
!  the source of truth): the SIMPLE-CUBIC 8-atom supercell is used, NOT the FCC
!  primitive cell. The plane-wave cutoff epm_pw_cutoff_ry bounds |G|^2 in
!  (2*pi/a)^2 units (integer shells h^2+k^2+l^2), and the FCC-in-cubic parity
!  selection rule on the form factors folds the 4 primitive BZs into the cubic
!  BZ exactly (32 folded bands / 32 valence electrons for the 8-atom cell).
!  k-points are written in REDUCED coordinates. With these conventions the
!  Fortran and Python solvers produce IDENTICAL k-points, band energies (to
!  ~5e-11 Ha) and occupations, and momentum matrix elements that agree to
!  machine precision up to the unavoidable degenerate-subspace basis freedom of
!  the top (truncation-boundary) conduction band. Verified by building SALMON
!  and diffing against epm_gaas_reference.py for GaAs and Si (scalar).
!
!  Because the pseudopotential is purely LOCAL, the velocity operator reduces to
!  v = p + A(t) (no nonlocal correction, see e.g. Yue & Gaarde, J. Opt. Soc. Am.
!  B 39, 535 (2022), Sec. on local EPM): the matrix rvnl_tm = -i[r,Vnl] is
!  identically zero and is written as such.
!
!  Parallelization: k-points are distributed over MPI ranks (split_range, the
!  same helper used by sbe), and the per-k Hamiltonian build / diagonalization /
!  momentum-matrix evaluation is additionally parallelized over OpenMP threads.
!  All LAPACK work arrays used inside the parallel region are subroutine-local
!  allocatables (eigen_zheev), hence thread-safe without extra bookkeeping.
!
module epm_solver
    use math_constants, only: pi, zi
    use phys_constants, only: au_ev
    use communication
    use epm_cohen_bergstresser
    implicit none

    private
    public :: s_epm_info, init_epm_info, run_epm_calculation, write_epm_files, finalize_epm_info

    type s_epm_info
        character(32) :: material
        character(16) :: cell = 'primitive'  ! 'primitive' (non-orthogonal, no folding) | 'folded' (cubic SC)

        ! Lattice / reciprocal lattice / cell volume
        real(8) :: a_lattice                 ! lattice constant [a.u.]
        real(8) :: a_matrix(3,3)             ! columns: a1, a2, a3
        real(8) :: b_matrix(3,3)             ! rows:    b1, b2, b3 (2*pi convention)
        real(8) :: volume
        real(8) :: tau(3)                    ! zincblende two-atom basis displacement (cubic path)

        ! Non-cubic (primitive-cell, no-folding) path: graphene, ... The cubic
        ! GaAs/Si path leaves these at their defaults and is byte-unchanged.
        logical :: noncubic = .false.
        integer :: natom = 0                 ! atoms in the basis (non-cubic)
        real(8) :: a_used = 0d0              ! in-plane lattice constant actually used [a.u.]
        real(8), allocatable :: tau_atoms(:,:)  ! (3,natom) basis positions [a.u.]

        ! Plane-wave basis (k-independent set of reciprocal lattice vectors G)
        integer :: npw
        integer,  allocatable :: Gindex(:,:) ! (3,npw) integer Miller indices m1,m2,m3
        real(8),  allocatable :: Gcart(:,:)  ! (3,npw) Cartesian G vectors [a.u.]
        integer,  allocatable :: G2(:)       ! (npw)   |G|^2 in units of (2*pi/a)^2 (rounded)

        ! k-point grid (Monkhorst-Pack, uniform weights -- no symmetry reduction)
        integer :: nk
        real(8), allocatable :: kpoint(:,:)  ! (3,nk) Cartesian [a.u.]
        real(8), allocatable :: kfrac(:,:)   ! (3,nk) reduced fractional coords (MP f_i)
        real(8), allocatable :: kweight(:)   ! (nk)

        ! Output data (filled by run_epm_calculation; full arrays on every rank
        ! after the cross-rank reduction so that rank 0 alone can write files)
        integer :: nb, ne
        real(8),    allocatable :: eigen(:,:)      ! (nb,nk)  [Ha]
        real(8),    allocatable :: occup(:,:)      ! (nb,nk)
        complex(8), allocatable :: p_tm(:,:,:,:)   ! (nb,nb,3,nk) = <u_m|p|u_n>
    end type s_epm_info

contains

    !=========================================================================
    ! Build lattice, plane-wave basis and k-grid (identical on every rank)
    !=========================================================================
    subroutine init_epm_info(epm, icomm)
        use salmon_global, only: epm_material, epm_lattice_constant_au, epm_pw_cutoff_ry, &
                                 epm_cell, num_kgrid, nstate, nelec
        use epm_noncubic, only: nc_is_noncubic, nc_lattice_and_basis
        implicit none
        type(s_epm_info), intent(out) :: epm
        integer, intent(in) :: icomm
        integer :: irank, nproc
        logical :: zincblende

        call comm_get_groupinfo(icomm, irank, nproc)

        epm%material   = epm_material
        epm%a_lattice  = epm_lattice_constant_au
        epm%cell       = epm_cell

        ! Material class: zincblende/diamond (GaAs/Si/Si_cb -- the Cohen-Bergstresser
        ! cubic form factors) vs the genuinely non-orthogonal cells (graphene; CdS
        ! pending) handled by the epm_noncubic module.
        zincblende = (trim(epm%material) == 'GaAs' .or. trim(epm%material) == 'Si' &
                      .or. trim(epm%material) == 'Si_kunikiyo' .or. trim(epm%material) == 'Si_cb')

        if (zincblende .and. trim(epm%cell) == 'folded') then
            ! ---- Legacy FOLDED path: simple-cubic 8-atom supercell + FCC-in-cubic
            ! parity folding (verified byte-equal to epm_gaas_reference.py; feeds
            ! the SBE folding/unfold pipeline). Opt-in via epm_cell='folded'. ----
            epm%noncubic = .false.
            call cb_lattice_vectors_sc(epm%a_lattice, epm%a_matrix(1:3,1), epm%a_matrix(1:3,2), epm%a_matrix(1:3,3))
            epm%tau(1:3) = cb_tau_zincblende(epm%a_lattice)
            call calc_reciprocal_lattice(epm%a_matrix, epm%b_matrix, epm%volume)
            call build_plane_wave_basis(epm, epm_pw_cutoff_ry)
        else if (zincblende) then
            ! ---- PRIMITIVE (default): FCC 2-atom non-orthogonal cell, NO folding.
            ! The FCC reciprocal (BCC) is the cubic (2pi/a)(h,k,l) with h,k,l all
            ! same parity, so the primitive basis is the cubic plane-wave basis
            ! restricted to ONE parity class -- reusing the SAME Cohen-Bergstresser
            ! Hamiltonian (build_hamiltonian) and momentum machinery verbatim.
            ! Matches epm_gaas_primitive.py (interchangeable output). ----
            epm%noncubic = .false.                 ! reuses the cubic CB Hamiltonian
            call fcc_primitive_vectors(epm%a_lattice, epm%a_matrix)
            epm%tau(1:3) = cb_tau_zincblende(epm%a_lattice)
            call calc_reciprocal_lattice(epm%a_matrix, epm%b_matrix, epm%volume)
            call build_plane_wave_basis_primitive(epm, epm_pw_cutoff_ry)
        else if (nc_is_noncubic(epm%material)) then
            ! ---- Non-orthogonal primitive cell with a GENERAL structure factor
            ! (graphene honeycomb; CdS wurtzite pending): epm_noncubic supplies the
            ! lattice + basis atoms + form factors. ----
            epm%noncubic = .true.
            call nc_lattice_and_basis(epm%material, epm%a_lattice, epm%a_matrix, &
                                      epm%natom, epm%tau_atoms, epm%a_used)
            call calc_reciprocal_lattice(epm%a_matrix, epm%b_matrix, epm%volume)
            call build_plane_wave_basis_noncubic(epm, epm_pw_cutoff_ry)
        else
            if (irank == 0) write(*,'(a,a)') '# ERROR: Fortran EPM has no path for material ', trim(epm%material)
            error stop 'epm_solver: unsupported epm_material (CdS primitive not yet ported -- use the Python reference)'
        end if

        epm%nk = num_kgrid(1) * num_kgrid(2) * num_kgrid(3)
        allocate(epm%kpoint(1:3, 1:epm%nk), epm%kfrac(1:3, 1:epm%nk), epm%kweight(1:epm%nk))
        call build_monkhorst_pack_grid(epm, num_kgrid)

        epm%nb = nstate
        epm%ne = nelec

        allocate(epm%eigen(1:epm%nb, 1:epm%nk))
        allocate(epm%occup(1:epm%nb, 1:epm%nk))
        allocate(epm%p_tm(1:epm%nb, 1:epm%nb, 1:3, 1:epm%nk))
        epm%eigen = 0d0
        epm%occup = 0d0
        epm%p_tm  = (0d0, 0d0)

        if (irank == 0) then
            write(*,'(a)')          '# EPM (local empirical pseudopotential)'
            write(*,'(a,a)')        '#   material           = ', trim(epm%material)
            if (trim(epm%cell) == 'folded') then
                write(*,'(a)')      '#   cell               = folded (cubic supercell, FCC-in-cubic folding)'
            else
                write(*,'(a)')      '#   cell               = primitive (non-orthogonal, no folding)'
            end if
            write(*,'(a,es12.5,a)') '#   lattice constant a = ', epm%a_lattice, ' a.u.'
            write(*,'(a,i8)')       '#   plane waves        = ', epm%npw
            write(*,'(a,i8)')       '#   k-points           = ', epm%nk
            write(*,'(a,i6,a,i6)')  '#   bands requested    = ', epm%nb, ' / valence electrons = ', epm%ne
        end if
    end subroutine init_epm_info


    subroutine finalize_epm_info(epm)
        implicit none
        type(s_epm_info), intent(inout) :: epm
        if (allocated(epm%Gindex))  deallocate(epm%Gindex)
        if (allocated(epm%Gcart))   deallocate(epm%Gcart)
        if (allocated(epm%G2))      deallocate(epm%G2)
        if (allocated(epm%kpoint))  deallocate(epm%kpoint)
        if (allocated(epm%kfrac))   deallocate(epm%kfrac)
        if (allocated(epm%kweight)) deallocate(epm%kweight)
        if (allocated(epm%tau_atoms)) deallocate(epm%tau_atoms)
        if (allocated(epm%eigen))   deallocate(epm%eigen)
        if (allocated(epm%occup))   deallocate(epm%occup)
        if (allocated(epm%p_tm))    deallocate(epm%p_tm)
    end subroutine finalize_epm_info


    ! b_i = 2*pi * (a_j x a_k) / V,  V = a1.(a2 x a3)
    subroutine calc_reciprocal_lattice(a_matrix, b_matrix, volume)
        implicit none
        real(8), intent(in)  :: a_matrix(3,3)
        real(8), intent(out) :: b_matrix(3,3)
        real(8), intent(out) :: volume
        real(8) :: a1(3), a2(3), a3(3), a23(3), a31(3), a12(3)

        a1 = a_matrix(1:3,1); a2 = a_matrix(1:3,2); a3 = a_matrix(1:3,3)
        a23 = cross(a2, a3)
        a31 = cross(a3, a1)
        a12 = cross(a1, a2)
        volume = dot_product(a1, a23)

        b_matrix(1,1:3) = (2d0*pi/volume) * a23(1:3)
        b_matrix(2,1:3) = (2d0*pi/volume) * a31(1:3)
        b_matrix(3,1:3) = (2d0*pi/volume) * a12(1:3)
    end subroutine calc_reciprocal_lattice


    function cross(u, v) result(w)
        implicit none
        real(8), intent(in) :: u(3), v(3)
        real(8) :: w(3)
        w(1) = u(2)*v(3) - u(3)*v(2)
        w(2) = u(3)*v(1) - u(1)*v(3)
        w(3) = u(1)*v(2) - u(2)*v(1)
    end function cross


    ! FCC PRIMITIVE real-space vectors (columns a1,a2,a3) [a.u.]:
    !   a1 = a/2 (0,1,1), a2 = a/2 (1,0,1), a3 = a/2 (1,1,0),  V = a^3/4.
    ! The 2-atom zincblende/diamond primitive cell (no folding). Matches
    ! epm_gaas_primitive.py::fcc_primitive_vectors_au.
    subroutine fcc_primitive_vectors(a_lattice, a_matrix)
        implicit none
        real(8), intent(in)  :: a_lattice
        real(8), intent(out) :: a_matrix(3,3)
        a_matrix(1:3,1) = 0.5d0 * a_lattice * (/ 0d0, 1d0, 1d0 /)
        a_matrix(1:3,2) = 0.5d0 * a_lattice * (/ 1d0, 0d0, 1d0 /)
        a_matrix(1:3,3) = 0.5d0 * a_lattice * (/ 1d0, 1d0, 0d0 /)
    end subroutine fcc_primitive_vectors


    !=========================================================================
    ! FCC-primitive plane-wave basis = cubic (2*pi/a)(h,k,l) with h,k,l ALL the
    ! same parity (the BCC reciprocal lattice = FCC reciprocal), |G|^2 <= cutoff
    ! in (2*pi/a)^2 units (integer shells h^2+k^2+l^2). This is the cubic basis
    ! restricted to ONE parity class (coset 0) -- NO folding. Gindex carries the
    ! cubic Miller indices (h,k,l) so build_hamiltonian's parity rule + form
    ! factors apply verbatim (dG is always all-even here, so V is fully active).
    ! Matches epm_gaas_primitive.py::build_pw_basis_fcc.
    !=========================================================================
    subroutine build_plane_wave_basis_primitive(epm, cutoff_ry)
        implicit none
        type(s_epm_info), intent(inout) :: epm
        real(8), intent(in) :: cutoff_ry
        real(8) :: twopi_a
        integer :: nmax, h, k, l, g2, n, npw_max
        integer, allocatable :: idx_tmp(:,:), g2_tmp(:)
        real(8), allocatable :: gcart_tmp(:,:)

        twopi_a = 2d0 * pi / epm%a_lattice
        nmax = ceiling(sqrt(cutoff_ry)) + 1
        npw_max = (2*nmax+1)**3
        allocate(idx_tmp(3, npw_max), gcart_tmp(3, npw_max), g2_tmp(npw_max))

        n = 0
        do h = -nmax, nmax
            do k = -nmax, nmax
                do l = -nmax, nmax
                    ! keep only one parity class (h,k,l all even or all odd)
                    if (.not. (mod(h-k,2) == 0 .and. mod(k-l,2) == 0)) cycle
                    g2 = h*h + k*k + l*l
                    if (dble(g2) <= cutoff_ry + 1d-8) then
                        n = n + 1
                        idx_tmp(1:3, n)   = (/ h, k, l /)
                        gcart_tmp(1:3, n) = twopi_a * (/ dble(h), dble(k), dble(l) /)
                        g2_tmp(n)         = g2
                    end if
                end do
            end do
        end do

        epm%npw = n
        allocate(epm%Gindex(3, n), epm%Gcart(3, n), epm%G2(n))
        epm%Gindex(1:3, 1:n) = idx_tmp(1:3, 1:n)
        epm%Gcart(1:3, 1:n)  = gcart_tmp(1:3, 1:n)
        epm%G2(1:n)          = g2_tmp(1:n)
        deallocate(idx_tmp, gcart_tmp, g2_tmp)
    end subroutine build_plane_wave_basis_primitive


    !=========================================================================
    ! General non-orthogonal plane-wave basis (graphene; CdS pending): keep all
    ! G = m1 b1 + m2 b2 + m3 b3 with the KINETIC energy 0.5|G|^2 [Ha] <= cutoff
    ! (epm_pw_cutoff_ry interpreted as a Ry kinetic cutoff: |G_au|^2 <= cutoff_ry,
    ! since E_kin[Ry] = |G_au|^2). For a 2D sheet embedded in a large vacuum cell
    ! the b3 (vacuum) direction is restricted to m3=0 (no kz dispersion), matching
    ! the strictly-2D Python reference. Gindex is unused on this path (the general
    ! structure factor uses Cartesian dG); G2 stores the rounded Cartesian |G|^2.
    !=========================================================================
    subroutine build_plane_wave_basis_noncubic(epm, cutoff_ry)
        implicit none
        type(s_epm_info), intent(inout) :: epm
        real(8), intent(in) :: cutoff_ry
        real(8) :: Gtmp(3), g2au
        integer :: nmax(3), m1, m2, m3, m3lo, m3hi, n, npw_max
        integer, allocatable :: idx_tmp(:,:), g2_tmp(:)
        real(8), allocatable :: gcart_tmp(:,:)
        logical :: twod

        ! 2D-in-vacuum materials: no dispersion along the vacuum axis (m3 = 0).
        twod = (trim(epm%material) == 'graphene')
        nmax(1) = ceiling(sqrt(cutoff_ry) / sqrt(dot_product(epm%b_matrix(1,1:3),epm%b_matrix(1,1:3)))) + 2
        nmax(2) = ceiling(sqrt(cutoff_ry) / sqrt(dot_product(epm%b_matrix(2,1:3),epm%b_matrix(2,1:3)))) + 2
        nmax(3) = ceiling(sqrt(cutoff_ry) / sqrt(dot_product(epm%b_matrix(3,1:3),epm%b_matrix(3,1:3)))) + 2
        if (twod) nmax(3) = 0
        m3lo = -nmax(3);  m3hi = nmax(3)

        npw_max = (2*nmax(1)+1) * (2*nmax(2)+1) * (2*nmax(3)+1)
        allocate(idx_tmp(3, npw_max), gcart_tmp(3, npw_max), g2_tmp(npw_max))
        n = 0
        do m1 = -nmax(1), nmax(1)
            do m2 = -nmax(2), nmax(2)
                do m3 = m3lo, m3hi
                    Gtmp(1:3) = m1*epm%b_matrix(1,1:3) + m2*epm%b_matrix(2,1:3) + m3*epm%b_matrix(3,1:3)
                    g2au = dot_product(Gtmp, Gtmp)
                    if (g2au <= cutoff_ry + 1d-8) then
                        n = n + 1
                        idx_tmp(1:3, n)   = (/ m1, m2, m3 /)
                        gcart_tmp(1:3, n) = Gtmp(1:3)
                        g2_tmp(n)         = nint(g2au)
                    end if
                end do
            end do
        end do

        epm%npw = n
        allocate(epm%Gindex(3, n), epm%Gcart(3, n), epm%G2(n))
        epm%Gindex(1:3, 1:n) = idx_tmp(1:3, 1:n)
        epm%Gcart(1:3, 1:n)  = gcart_tmp(1:3, 1:n)
        epm%G2(1:n)          = g2_tmp(1:n)
        deallocate(idx_tmp, gcart_tmp, g2_tmp)
    end subroutine build_plane_wave_basis_noncubic


    !=========================================================================
    ! General non-orthogonal Hamiltonian (graphene; CdS pending):
    !   H_ij = 0.5|k+G_i|^2 delta_ij + V_S(|dG|^2) * S(dG),
    !   S(dG) = sum_a exp(-i dG . tau_a)   (general multi-atom structure factor;
    !   V^A = 0 for the centrosymmetric materials handled here).
    ! Mirrors epm_graphene.py::build_hamiltonian (struct_norm=1, primitive).
    !=========================================================================
    subroutine build_hamiltonian_noncubic(epm, kvec, H)
        use epm_noncubic, only: nc_form_factor
        implicit none
        type(s_epm_info), intent(in) :: epm
        real(8), intent(in) :: kvec(3)
        complex(8), intent(out) :: H(epm%npw, epm%npw)
        integer :: i, j, ia
        real(8) :: kpg(3), dG(3), VS, dg2, ph
        complex(8) :: S

        do j = 1, epm%npw
            do i = 1, epm%npw
                if (i == j) then
                    kpg(1:3) = kvec(1:3) + epm%Gcart(1:3, i)
                    H(i, j) = dcmplx(0.5d0 * dot_product(kpg, kpg), 0d0)
                else
                    dG(1:3) = epm%Gcart(1:3, i) - epm%Gcart(1:3, j)
                    dg2 = dot_product(dG, dG)
                    VS = nc_form_factor(epm%material, dg2, epm%a_used)
                    if (VS == 0d0) then
                        H(i, j) = (0d0, 0d0)
                    else
                        S = (0d0, 0d0)
                        do ia = 1, epm%natom
                            ph = dot_product(dG, epm%tau_atoms(1:3, ia))
                            S = S + dcmplx(cos(ph), -sin(ph))   ! exp(-i dG.tau_a)
                        end do
                        H(i, j) = dcmplx(VS, 0d0) * S
                    end if
                end if
            end do
        end do
    end subroutine build_hamiltonian_noncubic


    !=========================================================================
    ! Plane-wave basis: fixed (k-independent) set of reciprocal lattice
    ! vectors G = m1*b1 + m2*b2 + m3*b3 with |G|^2 <= cutoff in units of
    ! (2*pi/a)^2 (i.e. h^2+k^2+l^2 <= epm_pw_cutoff_ry for the simple-cubic
    ! basis). This is the SAME convention as the Python reference
    ! epm_gaas_reference.py (cutoff on the integer shell index), so the two
    ! solvers select an IDENTICAL basis for the same epm_pw_cutoff_ry.
    !=========================================================================
    subroutine build_plane_wave_basis(epm, cutoff_ry)
        implicit none
        type(s_epm_info), intent(inout) :: epm
        real(8), intent(in) :: cutoff_ry
        real(8) :: gcut2
        real(8) :: g2_units
        real(8) :: Gtmp(3)
        integer :: nmax(3), m1, m2, m3, n, npw_max
        integer, allocatable :: idx_tmp(:,:)
        real(8), allocatable :: gcart_tmp(:,:)
        integer, allocatable :: g2_tmp(:)
        real(8) :: g2cart_to_units

        ! cutoff is on |G|^2 in (2*pi/a)^2 units (integer shells); convert a
        ! Cartesian |G|^2 [a.u.^2] to those units with (a/2*pi)^2.
        gcut2 = cutoff_ry
        g2cart_to_units = (epm%a_lattice / (2d0*pi))**2
        nmax(1:3) = ceiling(sqrt(gcut2)) + 1

        npw_max = (2*nmax(1)+1) * (2*nmax(2)+1) * (2*nmax(3)+1)
        allocate(idx_tmp(3, npw_max), gcart_tmp(3, npw_max), g2_tmp(npw_max))

        n = 0
        do m1 = -nmax(1), nmax(1)
            do m2 = -nmax(2), nmax(2)
                do m3 = -nmax(3), nmax(3)
                    Gtmp(1:3) = m1*epm%b_matrix(1,1:3) + m2*epm%b_matrix(2,1:3) + m3*epm%b_matrix(3,1:3)
                    g2_units = dot_product(Gtmp, Gtmp) * g2cart_to_units
                    if (g2_units <= gcut2 + 1.0d-8) then
                        n = n + 1
                        idx_tmp(1:3, n) = (/ m1, m2, m3 /)
                        gcart_tmp(1:3, n) = Gtmp(1:3)
                        g2_tmp(n) = nint(g2_units)
                    end if
                end do
            end do
        end do

        epm%npw = n
        allocate(epm%Gindex(3, n), epm%Gcart(3, n), epm%G2(n))
        epm%Gindex(1:3, 1:n) = idx_tmp(1:3, 1:n)
        epm%Gcart(1:3, 1:n)  = gcart_tmp(1:3, 1:n)
        epm%G2(1:n)          = g2_tmp(1:n)

        deallocate(idx_tmp, gcart_tmp, g2_tmp)
    end subroutine build_plane_wave_basis


    !=========================================================================
    ! Uniform Monkhorst-Pack grid with nk = n1*n2*n3 points and equal weights
    ! 1/nk (no symmetry reduction). gs_info_ssbe places no requirement on the
    ! ordering/weights of k-points beyond "they sum to a consistent total", so
    ! any convention is acceptable as long as nk matches num_kgrid in the SBE
    ! input -- which it does by construction here.
    !=========================================================================
    subroutine build_monkhorst_pack_grid(epm, num_kgrid)
        implicit none
        type(s_epm_info), intent(inout) :: epm
        integer, intent(in) :: num_kgrid(3)
        integer :: n1, n2, n3, i1, i2, i3, ik
        real(8) :: f1, f2, f3

        n1 = num_kgrid(1); n2 = num_kgrid(2); n3 = num_kgrid(3)
        ik = 0
        do i1 = 1, n1
            do i2 = 1, n2
                do i3 = 1, n3
                    ik = ik + 1
                    f1 = dble(2*i1 - n1 - 1) / dble(2*n1)
                    f2 = dble(2*i2 - n2 - 1) / dble(2*n2)
                    f3 = dble(2*i3 - n3 - 1) / dble(2*n3)
                    epm%kpoint(1:3, ik) = f1*epm%b_matrix(1,1:3) + f2*epm%b_matrix(2,1:3) + f3*epm%b_matrix(3,1:3)
                    ! reduced fractional coords in the reciprocal basis are exactly
                    ! (f1,f2,f3) (k = f1 b1 + f2 b2 + f3 b3) -- written verbatim for
                    ! non-orthogonal cells (matches kfrac = k @ inv(b) of the ref).
                    epm%kfrac(1:3, ik) = (/ f1, f2, f3 /)
                    epm%kweight(ik) = 1.0d0 / dble(epm%nk)
                end do
            end do
        end do
    end subroutine build_monkhorst_pack_grid


    !=========================================================================
    ! H_{G,G'}(k) = (1/2)|k+G|^2 delta_{G,G'}
    !             + [V^S(|dG|^2) cos(dG."tau") + i V^A(|dG|^2) sin(dG."tau")]
    !               * [parity(dG)]
    !
    ! with dG = G-G' and the FCC-in-cubic PARITY SELECTION RULE: the 8-atom
    ! supercell structure factor vanishes unless the integer indices
    ! (dh,dk,dl) of dG all share the same parity (all even or all odd). This is
    ! the band-folding trick of the Python reference (epm_gaas_reference.py):
    ! it makes H block-diagonal over the 4 FCC reciprocal cosets to machine
    ! precision, folding the 4 primitive BZs into the cubic BZ. (For the
    ! Cohen-Bergstresser shells |dG|^2 in {3,4,8,11} the rule is automatically
    ! satisfied, but it is enforced explicitly for exactness and generality.)
    !
    ! (kinetic term carries the standard Hartree-atomic-unit factor 1/2; the
    !  form factors V^S, V^A returned by cb_get_form_factors are already
    !  converted Ry -> Ha, i.e. divided by 2.)
    !=========================================================================
    subroutine build_hamiltonian(epm, kvec, H)
        implicit none
        type(s_epm_info), intent(in) :: epm
        real(8), intent(in) :: kvec(3)
        complex(8), intent(out) :: H(epm%npw, epm%npw)
        integer :: i, j, dG2, dh, dk, dl
        real(8) :: kpg(3), dG(3), VS, VA, phase

        do j = 1, epm%npw
            do i = 1, epm%npw
                if (i == j) then
                    kpg(1:3) = kvec(1:3) + epm%Gcart(1:3, i)
                    H(i, j) = dcmplx(0.5d0 * dot_product(kpg, kpg), 0d0)
                else
                    dh = epm%Gindex(1, i) - epm%Gindex(1, j)
                    dk = epm%Gindex(2, i) - epm%Gindex(2, j)
                    dl = epm%Gindex(3, i) - epm%Gindex(3, j)
                    ! parity selection rule: dh,dk,dl all same parity
                    if (mod(dh - dk, 2) == 0 .and. mod(dk - dl, 2) == 0) then
                        dG2 = dh*dh + dk*dk + dl*dl
                        call cb_get_form_factors(epm%material, dG2, VS, VA)
                        if (VS == 0d0 .and. VA == 0d0) then
                            H(i, j) = (0d0, 0d0)
                        else
                            dG(1:3) = epm%Gcart(1:3, i) - epm%Gcart(1:3, j)
                            phase = dot_product(dG, epm%tau)
                            H(i, j) = dcmplx(VS * cos(phase), VA * sin(phase))
                        end if
                    else
                        H(i, j) = (0d0, 0d0)
                    end if
                end if
            end do
        end do
    end subroutine build_hamiltonian


    ! p_{mn}(k) = sum_G conjg(c_m(G)) * (k+G) * c_n(G)   (diagonal in the plane-wave basis)
    subroutine calc_momentum_matrix(epm, kvec, evec, p_mn)
        implicit none
        type(s_epm_info), intent(in) :: epm
        real(8), intent(in) :: kvec(3)
        complex(8), intent(in)  :: evec(epm%npw, epm%nb)   ! columns = eigenvectors (lowest epm%nb states)
        complex(8), intent(out) :: p_mn(epm%nb, epm%nb, 3)
        complex(8), allocatable :: Dc(:,:)
        integer :: idir, ig
        real(8) :: kpg

        allocate(Dc(epm%npw, epm%nb))
        do idir = 1, 3
            do ig = 1, epm%npw
                kpg = kvec(idir) + epm%Gcart(idir, ig)
                Dc(ig, 1:epm%nb) = dcmplx(kpg, 0d0) * evec(ig, 1:epm%nb)
            end do
            call ZGEMM('C', 'N', epm%nb, epm%nb, epm%npw, dcmplx(1d0,0d0), &
                       evec, epm%npw, Dc, epm%npw, dcmplx(0d0,0d0), p_mn(:,:,idir), epm%nb)
        end do
        deallocate(Dc)
    end subroutine calc_momentum_matrix


    !=========================================================================
    ! Main driver: distribute k-points over MPI ranks (split_range, as in
    ! init_sbe_bloch_solver), diagonalize H(k) and build p_mn for each local k
    ! (OpenMP-parallel over k within a rank), then reduce the disjoint
    ! per-rank contributions onto every rank with comm_summation(..., dest=0)
    ! -- a zero-padding all-to-one sum is exact here because the k-ranges of
    ! different ranks never overlap. Rank 0 ends up holding the full dataset
    ! and is the only one that writes the output files (write_epm_files).
    !=========================================================================
    subroutine run_epm_calculation(epm, icomm)
        use util_ssbe, only: split_range
        implicit none
        type(s_epm_info), intent(inout) :: epm
        integer, intent(in) :: icomm
        integer :: irank, nproc, ik_min, ik_max, ik, ib
        integer, allocatable :: itbl_min(:), itbl_max(:)
        real(8),    allocatable :: eigen_local(:,:), occup_local(:,:)
        complex(8), allocatable :: p_tm_local(:,:,:,:)
        complex(8), allocatable :: H(:,:), evec(:,:), p_mn(:,:,:)
        real(8),    allocatable :: eval(:)
        integer :: nb, npw, nocc

        call comm_get_groupinfo(icomm, irank, nproc)

        nb  = epm%nb
        npw = epm%npw

        allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))
        call split_range(1, epm%nk, nproc, itbl_min, itbl_max)
        ik_min = itbl_min(irank)
        ik_max = itbl_max(irank)

        allocate(eigen_local(nb, epm%nk), occup_local(nb, epm%nk), p_tm_local(nb, nb, 3, epm%nk))
        eigen_local = 0d0
        occup_local = 0d0
        p_tm_local  = (0d0, 0d0)

        nocc = epm%ne / 2  ! doubly-occupied valence bands (closed-shell GaAs: 8 valence e- / 2 atoms)

        !$omp parallel default(shared) private(ik, ib, H, evec, p_mn, eval)
        allocate(H(npw, npw), evec(npw, npw), eval(npw), p_mn(nb, nb, 3))
        !$omp do schedule(dynamic)
        do ik = ik_min, ik_max
            if (epm%noncubic) then
                call build_hamiltonian_noncubic(epm, epm%kpoint(1:3, ik), H)
            else
                call build_hamiltonian(epm, epm%kpoint(1:3, ik), H)
            end if
            call eigen_zheev_wrap(H, eval, evec)

            do ib = 1, nb
                eigen_local(ib, ik) = eval(ib)
                if (ib <= nocc) then
                    occup_local(ib, ik) = 2.0d0
                else
                    occup_local(ib, ik) = 0.0d0
                end if
            end do

            call calc_momentum_matrix(epm, epm%kpoint(1:3, ik), evec(1:npw, 1:nb), p_mn)
            p_tm_local(1:nb, 1:nb, 1:3, ik) = p_mn(1:nb, 1:nb, 1:3)
        end do
        !$omp end do
        deallocate(H, evec, eval, p_mn)
        !$omp end parallel

        ! Disjoint-range reduction: sum over ranks reproduces the full array
        ! because each rank contributes zeros outside its own [ik_min,ik_max].
        call comm_summation(eigen_local, epm%eigen, size(epm%eigen), icomm)
        call comm_summation(occup_local, epm%occup, size(epm%occup), icomm)
        call comm_summation(p_tm_local,  epm%p_tm,  size(epm%p_tm),  icomm)

        deallocate(eigen_local, occup_local, p_tm_local)
        deallocate(itbl_min, itbl_max)
    end subroutine run_epm_calculation


    ! Thin wrapper so that eigen_zheev (module eigen_lapack, src/gs) is the
    ! single source of truth for the LAPACK ZHEEV call; all of its work arrays
    ! are subroutine-local allocatables, hence safe to call concurrently from
    ! different OpenMP threads.
    subroutine eigen_zheev_wrap(h, e, v)
        use eigen_lapack, only: eigen_zheev
        implicit none
        complex(8), intent(in)  :: h(:,:)
        real(8),    intent(out) :: e(:)
        complex(8), intent(out) :: v(:,:)
        call eigen_zheev(h, e, v)
    end subroutine eigen_zheev_wrap


    !=========================================================================
    ! Write SYSNAME_k.data / SYSNAME_eigen.data / SYSNAME_tm.data in exactly
    ! the format expected by gs_info_ssbe::read_k_data / read_eigen_data /
    ! read_tm_data (src/ssbe/gs_info_ssbe.f90). Those routines parse with
    ! free-format read(fh,*) but rely on a *fixed number of header lines* and
    ! a *fixed number/order of numeric fields* per data line -- both are
    ! reproduced verbatim below. Must be called on rank 0 only (the caller,
    ! main_epm, guards this with comm_is_root).
    !=========================================================================
    subroutine write_epm_files(epm, sysname, gs_directory)
        use filesystem, only: get_filehandle
        implicit none
        type(s_epm_info), intent(in) :: epm
        character(*), intent(in) :: sysname
        character(*), intent(in) :: gs_directory
        integer :: fh, ik, ib, jb, idir
        real(8) :: b_diag(3), kred(3)

        ! --- SYSNAME_k.data ---------------------------------------------------
        ! read_k_data skips any number of '#'/blank header lines then reads
        ! "ik, kx,ky,kz, weight". k-points are written in REDUCED (fractional)
        ! coordinates -- the convention the SBE uses and the Python reference emits.
        ! For a NON-orthogonal primitive cell the reciprocal vectors are written
        ! into the header (# b1/# b2/# b3) so the plotter can un-shear the triclinic
        ! grid; the fractional coords are the well-defined MP f_i (epm%kfrac). For
        ! the legacy folded cubic cell (diagonal b) kfrac == kpoint/b_diag exactly.
        fh = get_filehandle()
        open(unit=fh, file=trim(gs_directory)//trim(sysname)//'_k.data', action='write', status='replace')
        write(fh, '(A)') '# k-point data'
        write(fh, '(A)') '# generated by EPM ('//trim(epm%cell)//' cell)'
        write(fh, '(A,A,A,I8)') '# material = ', trim(epm%material), ', nk = ', epm%nk
        if (trim(epm%cell) /= 'folded') then
            do idir = 1, 3
                write(fh, '(A,I1,A,3E18.10,A)') '# b', idir, ' = ', epm%b_matrix(idir,1:3), '  [a.u.]'
            end do
        end if
        write(fh, '(A)') '# units: kx,ky,kz [reduced fractional of the reciprocal lattice], weight'
        write(fh, '(A)') '# ik, kx, ky, kz, weight'
        if (trim(epm%cell) == 'folded') then
            b_diag(1) = epm%b_matrix(1,1)
            b_diag(2) = epm%b_matrix(2,2)
            b_diag(3) = epm%b_matrix(3,3)
            do ik = 1, epm%nk
                kred(1:3) = epm%kpoint(1:3, ik) / b_diag(1:3)
                write(fh, '(I6, 4E18.10)') ik, kred(1), kred(2), kred(3), epm%kweight(ik)
            end do
        else
            do ik = 1, epm%nk
                write(fh, '(I6, 4E18.10)') ik, epm%kfrac(1,ik), epm%kfrac(2,ik), epm%kfrac(3,ik), epm%kweight(ik)
            end do
        end if
        close(fh)

        ! --- SYSNAME_eigen.data ------------------------------------------------
        ! read_eigen_data consumes exactly 3 header lines, then for each ik:
        ! one header line followed by nb lines "ib, energy[Ha], occupancy"
        fh = get_filehandle()
        open(unit=fh, file=trim(gs_directory)//trim(sysname)//'_eigen.data', action='write', status='replace')
        write(fh, '(A)') '# eigenvalue data'
        write(fh, '(A)') '# generated by EPM (Cohen-Bergstresser local pseudopotential)'
        write(fh, '(A,I6,A,I6)') '# nk = ', epm%nk, ', nb = ', epm%nb
        do ik = 1, epm%nk
            write(fh, '(A,I6)') '# ik = ', ik
            do ib = 1, epm%nb
                write(fh, '(I6, 2E18.10)') ib, epm%eigen(ib, ik), epm%occup(ib, ik)
            end do
        end do
        close(fh)

        ! --- SYSNAME_tm.data ----------------------------------------------------
        ! read_tm_data consumes exactly 3 header lines, then nk*nb*nb lines of
        ! "ik, ib, jb, Re px, Im px, Re py, Im py, Re pz, Im pz" (block 1: p_tm),
        ! then ONE more header line, then the same nk*nb*nb lines for rvnl_tm
        ! (block 2). Local pseudopotential => rvnl_tm is identically zero
        ! (v = p + A, no nonlocal velocity correction; see module header).
        fh = get_filehandle()
        open(unit=fh, file=trim(gs_directory)//trim(sysname)//'_tm.data', action='write', status='replace')
        write(fh, '(A)') '# transition matrix data'
        write(fh, '(A)') '# generated by EPM (Cohen-Bergstresser local pseudopotential)'
        write(fh, '(A)') '# block 1: p_tm = <u_m|p|u_n>  (ik, ib, jb, Re px, Im px, Re py, Im py, Re pz, Im pz)'
        do ik = 1, epm%nk
            do ib = 1, epm%nb
                do jb = 1, epm%nb
                    write(fh, '(3I6, 6E18.10)') ik, ib, jb, &
                        & (real(epm%p_tm(ib,jb,idir,ik)), aimag(epm%p_tm(ib,jb,idir,ik)), idir = 1, 3)
                end do
            end do
        end do
        write(fh, '(A)') '# block 2: rvnl_tm = -i[r,Vnl]  (all zero: local pseudopotential, no nonlocal correction)'
        do ik = 1, epm%nk
            do ib = 1, epm%nb
                do jb = 1, epm%nb
                    write(fh, '(3I6, 6E18.10)') ik, ib, jb, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0
                end do
            end do
        end do
        close(fh)

        write(*,'(a,a)') '# EPM: wrote ground-state data files for sysname = ', trim(sysname)
    end subroutine write_epm_files

end module epm_solver
