! Ground State Date Storage Module:

module gs_info_ssbe
    use math_constants, only: pi, zI
    use phys_constants, only: au_ev, kB_au
    implicit none

    type s_sbe_gs_info
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kpoint(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: occup_ref(:, :)   ! UNDOPED, T -> 0 filling: the pure-gauge reference of the f-sum-rule restoration (wiki/12 sec. 6a); never the doped/thermal occupation
        real(8), allocatable :: delta_omega(:, :, :)
        complex(8), allocatable :: p_mod_matrix(:, :, :, :)
        ! p_tm_matrix = <u|p|u>
        complex(8), allocatable :: p_tm_matrix(:, :, :, :)
        ! rvnl_tm_matrix = <u|-i[r, Vnl]|u>
        complex(8), allocatable :: rvnl_tm_matrix(:, :, :, :)
        complex(8), allocatable :: d_matrix(:, :, :, :)

        !k-space grid and geometry information
        !NOTE: prepred for uniformally distributed k-grid....
        !integer :: num_kgrid(1:3)

        ! Minimum band gap in atomic units (for gauge-covariant decoherence)
        real(8) :: eg_au

        ! Optional unfold map (SYSNAME_unfold.data from the cubic-supercell
        ! EPM): assigns every supercell band to its FCC sublattice (= folded
        ! primitive BZ point k_prim = k_sc + G0) and primitive band index.
        ! Used to output populations of PHYSICAL primitive bands instead of
        ! energy-ordered supercell branches.
        logical :: have_unfold = .false.
        integer :: nv_prim = 0
        ! Number of folding cosets/sublattices: 4 for the cubic FCC fold
        ! (GaAs/Si), 2 for the wurtzite CdS / rectangular graphene folds. The
        ! arrays below are dimensioned for the maximum (4); slots > n_coset are
        ! left zero (so the 4-coset path is byte-unchanged) and the output/
        ! population loops run only over 1..n_coset.
        integer :: n_coset = 4
        integer, allocatable :: unfold_sub(:, :)     ! (nb, nk) dominant coset 1..n_coset
        integer, allocatable :: unfold_prim(:, :)    ! (nb, nk) primitive band rank
        real(8), allocatable :: unfold_w(:, :, :)    ! (4, nb, nk) spectral weights, sum_s = 1
        real(8) :: unfold_offset(1:3, 1:4)           ! G0 in sc reduced coords (1..n_coset used)
    end type


contains


subroutine init_sbe_gs_info(gs, sysname, gs_directory, nk, nb, ne, a1, a2, a3, read_bin, icomm)
    use communication
    use filesystem, only: open_filehandle, get_filehandle
    use salmon_global, only: yn_sbe_spinor, sbe_ef_ev, sbe_temp_init_k
    implicit none
    type(s_sbe_gs_info), intent(inout) :: gs
    character(*), intent(in) :: sysname
    character(*), intent(in) :: gs_directory
    integer, intent(in) :: nk
    integer, intent(in) :: nb
    integer, intent(in) :: ne
    real(8), intent(in) :: a1(1:3), a2(1:3), a3(1:3)
    logical, intent(in) :: read_bin
    integer, intent(in) :: icomm
    integer :: irank, nproc

    call comm_get_groupinfo(icomm, irank, nproc)

    gs%nk = nk
    gs%nb = nb
    gs%ne = ne
    !gs%num_kgrid(1:3) = num_kgrid(1:3)

    !Calculate b_matrix, volume_cell and volume_bz from a1..a3 vector.
    call calc_lattice_info()

    allocate(gs%kpoint(1:3, 1:nk))
    allocate(gs%kweight(1:nk))
    allocate(gs%eigen(1:nb, 1:nk))
    allocate(gs%occup(1:nb, 1:nk))
    allocate(gs%occup_ref(1:nb, 1:nk))
    allocate(gs%delta_omega(1:nb, 1:nb, 1:nk))
    allocate(gs%p_mod_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%d_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%p_tm_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%rvnl_tm_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%unfold_sub(1:nb, 1:nk), gs%unfold_prim(1:nb, 1:nk))
    allocate(gs%unfold_w(1:4, 1:nb, 1:nk))
    gs%unfold_sub = 0
    gs%unfold_prim = 0
    gs%unfold_w = 0d0
    gs%unfold_offset = 0d0

    if (irank == 0) then
        if (read_bin) then
            !Retrieve all data from binray
            write(*,*) "# read_sbe_gs_bin"
            call read_sbe_gs_bin()
        else
            !Retrieve eigenenergies from 'SYSNAME_eigen.data':
            write(*, '(a)') "# read_eigen_data"
            call read_eigen_data()
            !Retrieve k-points from 'SYSNAME_k.data':
            write(*, '(a)') "# read_k_data"
            call read_k_data()
            !Retrieve transition matrix from 'SYSNAME_tm.data':
            write(*, '(a)') "# read_tm_data"
            call read_tm_data()
            !Optional: band -> sublattice unfold map 'SYSNAME_unfold.data'
            call read_unfold_data()
            !Export all data from binray
            write(*, '(a)') "# save_sbe_gs_bin"
            call save_sbe_gs_bin()
        end if
    end if

    call comm_bcast(gs%kpoint, icomm, 0)
    call comm_bcast(gs%kweight, icomm, 0)
    call comm_bcast(gs%eigen, icomm, 0)
    call comm_bcast(gs%occup, icomm, 0)
    call comm_bcast(gs%p_tm_matrix, icomm, 0)
    call comm_bcast(gs%rvnl_tm_matrix, icomm, 0)
    call comm_bcast(gs%have_unfold, icomm, 0)
    if (gs%have_unfold) then
        call comm_bcast(gs%nv_prim, icomm, 0)
        call comm_bcast(gs%n_coset, icomm, 0)
        call comm_bcast(gs%unfold_sub, icomm, 0)
        call comm_bcast(gs%unfold_prim, icomm, 0)
        call comm_bcast(gs%unfold_w, icomm, 0)
        call comm_bcast(gs%unfold_offset, icomm, 0)
    end if

    !Calculate omega and d_matrix (neglecting diagonal part):
    if (irank == 0) write(*,"(a)") "# prepare_matrix"

    call prepare_matrix()
    call comm_bcast(gs%p_mod_matrix, icomm, 0)
    call comm_bcast(gs%delta_omega, icomm, 0)
    call comm_bcast(gs%d_matrix, icomm, 0) ! Experimental

    !Initial Occupation Number
    gs%occup(:,:) = 0d0 !!Experimental!!
    if (yn_sbe_spinor == 'y') then
        ! Spinor (spin-orbit split) input: one electron per spinor band
        gs%occup(1:ne,:) = 1d0
    else
        gs%occup(1:(ne/2),:) = 2d0 !!Experimental!!
    end if

    ! T -> 0 limit at an exactly degenerate, partially filled level group (a Dirac
    ! point sitting on the k-mesh, a metal's Fermi crossing on a mesh point): the
    ! integer filling above picks LAPACK's arbitrary basis inside the degenerate
    ! subspace -- a broken-symmetry state that carries a spurious current of order
    ! v_F/N_k per such k-point (decisive for a 2D sheet at low THz fields: with the
    ! self-consistent sheet field it acts as a relay that pins the local field to
    ! zero). The correct zero-temperature density matrix is the group average
    ! (identity on the degenerate block): equal fractional occupation, no current,
    ! invariant under any rotation of the block -- and the pure-gauge reference of
    ! wiki/12 sec. 6a becomes continuous in A. Inert for gapped materials.
    block
        integer :: ik, i1, i2, nsym
        real(8) :: osum
        real(8), parameter :: deg_tol = 1d-6   ! Ha
        nsym = 0
        do ik = 1, nk
            i1 = 1
            do while (i1 <= nb)
                i2 = i1
                do while (i2 < nb)
                    if (abs(gs%eigen(i2 + 1, ik) - gs%eigen(i2, ik)) > deg_tol) exit
                    i2 = i2 + 1
                end do
                if (i2 > i1) then
                    if (maxval(gs%occup(i1:i2, ik)) - minval(gs%occup(i1:i2, ik)) > 1d-12) then
                        osum = sum(gs%occup(i1:i2, ik))
                        gs%occup(i1:i2, ik) = osum / dble(i2 - i1 + 1)
                        nsym = nsym + 1
                    end if
                end if
                i1 = i2 + 1
            end do
        end do
        if (irank == 0 .and. nsym > 0) write(*, '(a,i0,a)') &
            '# occupations: degenerate partially filled level groups averaged at ', nsym, &
            ' k-point(s) (T -> 0 density matrix; e.g. the Dirac point on the mesh)'
    end block

    ! =====================================================================
    ! The UNDOPED, T -> 0 filling is kept as gs%occup_ref: it is the reference
    ! of the velocity-gauge pure-gauge restoration (wiki/12 sec. 6a). That
    ! subtraction must remove only the truncation artifact of the FILLED sea;
    ! evaluated on a doped/thermal occupation it would also remove the
    ! physical intraband (Drude) current of the doping carriers, which is the
    ! very quantity a doped sheet is run for.
    ! =====================================================================
    gs%occup_ref(:, :) = gs%occup(:, :)

    ! =====================================================================
    ! Doped / finite-temperature INITIAL occupation (wiki/12 sec. 4a):
    !   f_n(k) = occ_max * f_FD(eps_n(k); mu, T_init),   mu = E_F^undoped + sbe_ef_ev
    ! E_F^undoped = midpoint of the HOMO ceiling and the LUMO floor (the Dirac
    ! point for graphene, mid-gap for a semiconductor). This is what gives a
    ! partially filled band, hence a Drude weight and an intraband conductivity;
    ! with sbe_ef_ev = 0 and T_init = 0 nothing changes (integer filling).
    ! The added charge is NOT compensated (a gated/adsorbate-doped sheet); the
    ! electron count per cell changes accordingly and is reported.
    ! =====================================================================
    if (abs(sbe_ef_ev) > 0d0 .or. sbe_temp_init_k > 0d0) then
        block
            integer :: ik, ib, nb_vb2, nfs
            real(8) :: e_ref, mu_au, beta, dne, area, n2d_cm2, x
            real(8), parameter :: BOHR_CM = 0.52917721067d-8
            if (yn_sbe_spinor == 'y') then
                nb_vb2 = gs%ne
            else
                nb_vb2 = gs%ne / 2
            end if
            if (nb_vb2 >= 1 .and. nb_vb2 < nb) then
                e_ref = 0.5d0 * (maxval(gs%eigen(nb_vb2, :)) + minval(gs%eigen(nb_vb2 + 1, :)))
            else
                e_ref = 0d0
            end if
            mu_au = e_ref + sbe_ef_ev / au_ev
            if (sbe_temp_init_k > 0d0) then
                beta = 1d0 / (kB_au * sbe_temp_init_k)
                do ik = 1, nk
                    do ib = 1, nb
                        x = beta * (gs%eigen(ib, ik) - mu_au)
                        if (x > 60d0) then
                            gs%occup(ib, ik) = 0d0
                        else if (x < -60d0) then
                            gs%occup(ib, ik) = occ_max_gs()
                        else
                            gs%occup(ib, ik) = occ_max_gs() / (1d0 + exp(x))
                        end if
                    end do
                end do
            else
                do ik = 1, nk
                    do ib = 1, nb
                        if (gs%eigen(ib, ik) < mu_au - 1d-9) then
                            gs%occup(ib, ik) = occ_max_gs()
                        else if (gs%eigen(ib, ik) > mu_au + 1d-9) then
                            gs%occup(ib, ik) = 0d0
                        else
                            gs%occup(ib, ik) = 0.5d0 * occ_max_gs()
                        end if
                    end do
                end do
            end if
            ! Fermi-surface resolution: a doped metal is representable on a
            ! uniform mesh only if enough k-points fall in the FD transition
            ! ring. Count the partially occupied points (5..95 % filling); for
            ! a Dirac cone at E_F this is ~ the perimeter 2 pi k_F / dk per
            ! valley. Too few -> the density and the Drude weight are set by a
            ! handful of points and the intraband response is meaningless.
            nfs = count(gs%occup > 0.05d0 * occ_max_gs() .and. gs%occup < 0.95d0 * occ_max_gs())
            dne = sum(gs%occup) / dble(nk) - sum(gs%occup_ref) / dble(nk)
            area = sqrt(max(dot_product(cross3(gs%a_matrix(1:3,1), gs%a_matrix(1:3,2)), &
                                        cross3(gs%a_matrix(1:3,1), gs%a_matrix(1:3,2))), 1d-300))
            n2d_cm2 = dne / (area * BOHR_CM**2)
            if (irank == 0) then
                write(*, '(a,f8.4,a,f8.1,a)') '# doped/thermal initial occupation: E_F = ', sbe_ef_ev, &
                    ' eV from the undoped level, T_init = ', sbe_temp_init_k, ' K'
                write(*, '(a,es12.5,a,es12.5,a)') '#   added charge = ', dne, &
                    ' e/cell  ->  sheet density = ', n2d_cm2, ' cm^-2 (2D cell area used)'
                write(*, '(a)') '#   pure-gauge f-sum-rule reference stays the UNDOPED filling (wiki/12 sec. 6a)'
                write(*, '(a,i0,a)') '#   partially occupied k-points (Fermi surface on the mesh): ', nfs, &
                    ' -- the intraband/Drude response is carried by these'
                if (nfs == 0 .and. abs(dne) > 0d0) then
                    ! Qualitatively worse than "few": the doping charge has NO partially
                    ! occupied state to occupy, so it lands entirely on fully filled or
                    ! empty levels and the initial state is not a Fermi-Dirac metal at
                    ! all. Observed consequence (33^2, E_F = 0.2 eV, ring on): population
                    ! inversion and a sheet with NEGATIVE absorption, A = -0.24, the
                    ! electrons pumping 1.2e-3 eV/cell into the field.
                    write(*, '(a)') &
                        '#   *** ERROR-LEVEL WARNING: ZERO partially occupied k-points. The Fermi disc'// &
                        ' contains no mesh point, so the requested doping cannot be represented:'
                    write(*, '(a)') &
                        '#       the added charge lands on fully occupied / empty levels, the initial'// &
                        ' state is not a metal, and the run can develop population inversion and GAIN'
                    write(*, '(a)') &
                        '#       (negative absorption). Raise num_kgrid until k_F = E_F/hbar v_F spans'// &
                        ' at least a few mesh spacings, or raise |sbe_ef_ev|. Do not use this run.'
                else if (nfs < 20 .and. abs(dne) > 0d0) then
                    write(*, '(a)') &
                        '#   WARNING: the Fermi surface is UNDER-RESOLVED (< 20 partially occupied k-points):'// &
                        ' the carrier density and the Drude weight are set by a few mesh points.'// &
                        ' Raise num_kgrid (k_F = E_F/hbar v_F must exceed a few mesh spacings) or |sbe_ef_ev|.'
                end if
            end if
        end block
    end if

    ! Calculate minimum band gap in atomic units (for gauge-covariant decoherence)
    call calc_eg_au()

contains

    pure function occ_max_gs() result(o)
        real(8) :: o
        o = merge(1d0, 2d0, yn_sbe_spinor == 'y')
    end function occ_max_gs

    pure function cross3(u, v) result(w)
        real(8), intent(in) :: u(3), v(3)
        real(8) :: w(3)
        w(1) = u(2) * v(3) - u(3) * v(2)
        w(2) = u(3) * v(1) - u(1) * v(3)
        w(3) = u(1) * v(2) - u(2) * v(1)
    end function cross3

    ! Calculate lattice and reciprocal vectors
    subroutine calc_lattice_info()
        implicit none
        real(8) :: a12(1:3), a23(1:3), a31(1:3), volume
        real(8) :: b1(1:3), b2(1:3), b3(1:3)

        a12(1) = a1(2) * a2(3) - a1(3) * a2(2)
        a12(2) = a1(3) * a2(1) - a1(1) * a2(3)
        a12(3) = a1(1) * a2(2) - a1(2) * a2(1)
        a23(1) = a2(2) * a3(3) - a2(3) * a3(2)
        a23(2) = a2(3) * a3(1) - a2(1) * a3(3)
        a23(3) = a2(1) * a3(2) - a2(2) * a3(1)
        a31(1) = a3(2) * a1(3) - a3(3) * a1(2)
        a31(2) = a3(3) * a1(1) - a3(1) * a1(3)
        a31(3) = a3(1) * a1(2) - a3(2) * a1(1)
        volume = dot_product(a12, a3)
        b1(1:3) = (2d0 * pi / volume) * a23(1:3)
        b2(1:3) = (2d0 * pi / volume) * a31(1:3)
        b3(1:3) = (2d0 * pi / volume) * a12(1:3)

        gs%a_matrix(1:3, 1) = a1(1:3)
        gs%a_matrix(1:3, 2) = a2(1:3)
        gs%a_matrix(1:3, 3) = a3(1:3)
        gs%b_matrix(1, 1:3) = b1(1:3)
        gs%b_matrix(2, 1:3) = b2(1:3)
        gs%b_matrix(3, 1:3) = b3(1:3)
        gs%volume = volume
    end subroutine calc_lattice_info


    ! Read k-point coordinates from SALMON's output file. Robustly skips ANY
    ! number of leading comment ('#') / blank lines, so the EPM may add header
    ! metadata (e.g. the non-orthogonal reciprocal vectors '# b1/# b2/# b3')
    ! without breaking the reader -- the data rows start with the integer ik.
    subroutine read_k_data()
        implicit none
        character(256) :: line
        integer :: fh, ik, iik, ios
        real(8) :: tmp(4)
        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_k.data', 'old')
        ik = 0
        do
            read(fh, "(a)", iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            ik = ik + 1
            read(line, *) iik, tmp(1:4)
            if (ik .ne. iik) stop "ik mismatch"
            gs%kpoint(1:3, ik) = tmp(1:3)
            gs%kweight(ik) = tmp(4)
            if (ik == nk) exit
        end do
        close(fh)
    end subroutine read_k_data


    ! Read eigenvalue data from SALMON's output file.
    ! Unit-aware: the SBE works in Hartree, but a DFT ground state run with
    ! unit_system='A_eV_fs' writes esp in eV -- write_eigen states its unit in
    ! the 3rd header line ("# 1:io, 2:esp[eV], 3:occ"). Detect that tag and
    ! convert; files without it (EPM writers, unit_system='au') are already a.u.
    subroutine read_eigen_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, iik, iib, ib
        real(8) :: tmp(2)
        real(8) :: e_conv

        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_eigen.data', 'old')
        e_conv = 1d0
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        if (index(dummy, 'esp[eV]') > 0) then
            e_conv = 1d0 / au_ev
            write(*, '(a)') "# read_eigen_data: esp[eV] header detected -> converting to Hartree"
        end if
        do ik = 1, nk
            read(fh, "(a)") dummy
            do ib = 1, nb
                read(fh, *) iib, tmp(1:2)
                if (ib .ne. iib) stop "ib mismatch"
                gs%eigen(ib, ik) = tmp(1) * e_conv
                ! gs%occup(ib, ik) = ctmp(2)
            end do
        end do
        close(fh)
    end subroutine read_eigen_data




    ! Read transition dipole moment from SALMON's output file
    subroutine read_tm_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, ib, jb, iik, iib, jjb
        real(8) :: tmp(1:6)


        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_tm.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    gs%p_tm_matrix(ib, jb, 1, ik) = cmplx(tmp(1), tmp(2), 8)
                    gs%p_tm_matrix(ib, jb, 2, ik) = cmplx(tmp(3), tmp(4), 8)
                    gs%p_tm_matrix(ib, jb, 3, ik) = cmplx(tmp(5), tmp(6), 8)
                end do
            end do
        end do
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    gs%rvnl_tm_matrix(ib, jb, 1, ik) = cmplx(tmp(1), tmp(2), 8)
                    gs%rvnl_tm_matrix(ib, jb, 2, ik) = cmplx(tmp(3), tmp(4), 8)
                    gs%rvnl_tm_matrix(ib, jb, 3, ik) = cmplx(tmp(5), tmp(6), 8)
                end do
            end do
        end do


        close(fh)
    end subroutine read_tm_data


    ! Optional: read the band -> (FCC sublattice, primitive band) unfold map
    ! written by `epm_gaas_reference.py unfoldmap`. Absence is not an error:
    ! the unfolded population output is simply disabled.
    subroutine read_unfold_data()
        implicit none
        character(256) :: dummy
        character(512) :: fpath
        logical :: exists
        integer :: fh, ik, ib, iik, iib, isub, ibprim, i, nnk, nnb, ioff(3), ios
        real(8) :: w(4)

        fpath = trim(gs_directory) // trim(sysname) // '_unfold.data'
        inquire(file=trim(fpath), exist=exists)
        gs%have_unfold = .false.
        if (.not. exists) then
            ! Absence is not an error: the physical (unfolded) population output
            ! is simply disabled. Logged so the user can tell whether the file
            ! was looked for in the directory they expect (gs_directory).
            write(*, '(a)') "# read_unfold_data: no unfold map found, " // &
                & "physical-band population output disabled"
            write(*, '(a)') "#   (searched: " // trim(fpath) // ")"
            return
        end if

        write(*, '(a)') "# read_unfold_data: " // trim(fpath)
        fh = open_filehandle(trim(fpath), 'old')
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        ! data header: "nk nb nv_prim [n_coset]". n_coset=4 (FCC cubic, GaAs/Si)
        ! or 2 (wurtzite CdS / rectangular graphene). Legacy files omit it -> 4.
        read(fh, "(a)") dummy
        read(dummy, *, iostat=ios) nnk, nnb, gs%nv_prim, gs%n_coset
        if (ios .ne. 0) then
            read(dummy, *) nnk, nnb, gs%nv_prim
            gs%n_coset = 4
        end if
        if (gs%n_coset < 1 .or. gs%n_coset > 4) then
            write(*, '(a,i0)') "# read_unfold_data: unsupported n_coset = ", gs%n_coset
            stop "unfold map: n_coset must be 1..4"
        end if
        if (nnk .ne. nk) then
            write(*, '(a,i0,a,i0)') "# read_unfold_data: nk mismatch -- file has ", &
                & nnk, ", SBE run expects ", nk
            stop "unfold map: nk mismatch"
        end if
        if (nnb .ne. nb) then
            write(*, '(a,i0,a,i0,a)') "# read_unfold_data: nb mismatch -- file has ", &
                & nnb, ", SBE run expects ", nb, &
                & " (check nstate and yn_sbe_spinor vs. the EPM dataset)"
            stop "unfold map: nb mismatch"
        end if
        read(fh, "(a)") dummy
        do i = 1, gs%n_coset
            read(fh, *) isub, ioff(1:3)
            if (isub .ne. i) stop "unfold map: offset index mismatch"
            gs%unfold_offset(1:3, i) = dble(ioff(1:3))
        end do
        read(fh, "(a)") dummy
        do ik = 1, nk
            do ib = 1, nb
                w = 0d0
                read(fh, *) iik, iib, isub, ibprim, w(1:gs%n_coset)
                if (ik .ne. iik) stop "unfold map: ik mismatch"
                if (ib .ne. iib) stop "unfold map: ib mismatch"
                gs%unfold_sub(ib, ik) = isub
                gs%unfold_prim(ib, ik) = ibprim
                gs%unfold_w(1:4, ib, ik) = w(1:4)
            end do
        end do
        close(fh)
        gs%have_unfold = .true.
        write(*, '(a,i0)') "# unfold map loaded: nv_prim = ", gs%nv_prim
    end subroutine read_unfold_data


    subroutine read_sbe_gs_bin()
        implicit none
        integer :: fh
        ! fh = get_filehandle()
        ! open(fh, file=trim(gs_directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='old')
        ! read(fh) gs%kpoint
        ! read(fh) gs%kweight
        ! read(fh) gs%eigen
        ! read(fh) gs%p_mod_matrix
        ! read(fh) gs%rvnl_tm_matrix
        ! ! read(fh) gs%prod_dk
        ! close(fh)
        ! return
    end subroutine read_sbe_gs_bin


    subroutine save_sbe_gs_bin()
        implicit none
        integer :: fh
        ! fh = get_filehandle()
        ! open(fh, file=trim(gs_directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='replace')
        ! write(fh) gs%kpoint
        ! write(fh) gs%kweight
        ! write(fh) gs%eigen
        ! write(fh) gs%p_mod_matrix
        ! write(fh) gs%rvnl_tm_matrix
        ! ! write(fh) gs%prod_dk
        ! close(fh)
        ! return
    end subroutine save_sbe_gs_bin




    subroutine prepare_matrix()
        implicit none
        integer :: ik, ib, jb
        real(8), parameter :: omega_eps = 1d-9

        gs%p_mod_matrix = gs%p_tm_matrix + gs%rvnl_tm_matrix

        do ik=1, nk
            do ib=1, nb
                do jb=1, nb
                    gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    if (omega_eps < abs(gs%delta_omega(ib, jb, ik))) then
                        ! gs%d_matrix(ib, jb, 1:3, ik) = &
                        !     & (zi * gs%p_mod_matrix(ib, jb, 1:3, ik) - gs%rvnl_tm_matrix(ib, jb, 1:3, ik)) &
                        !     & / gs%delta_omega(ib, jb, ik)
                        ! gs%d_matrix(ib, jb, 1:3, ik) = &
                        !     & zi * (gs%p_mod_matrix(ib, jb, 1:3, ik) +  gs%rvnl_tm_matrix(ib, jb, 1:3, ik)) &
                        !     & / gs%delta_omega(ib, jb, ik)
                        gs%d_matrix(ib, jb, 1:3, ik) = &
                            & zi * (gs%p_mod_matrix(ib, jb, 1:3, ik)) &
                            & / gs%delta_omega(ib, jb, ik)
                    else
                        gs%d_matrix(ib, jb, 1:3, ik) = 0d0
                    end if
                end do
            end do
        end do
    end subroutine prepare_matrix

    ! Calculate minimum band gap in atomic units
    subroutine calc_eg_au()
        use salmon_global, only: eg_ev
        implicit none
        integer :: ik, ib_cb, ib_vb, nb_vb
        real(8) :: eg_tmp

        ! Highest occupied band: ne spinor bands (occupation 1) or ne/2 scalar bands (occupation 2)
        if (yn_sbe_spinor == 'y') then
            nb_vb = gs%ne
        else
            nb_vb = gs%ne / 2
        end if

        ! Check if user specified eg_ev = -1 (automatic calculation from band structure)
        if (eg_ev < 0.0d0) then
            ! Automatic calculation: find minimum band gap across all k-points
            gs%eg_au = 1.0d99  ! Initialize with large value

            do ik = 1, gs%nk
                ! For each k-point, find the minimum gap between conduction and valence bands
                do ib_cb = nb_vb + 1, gs%nb
                    do ib_vb = 1, nb_vb
                        eg_tmp = gs%eigen(ib_cb, ik) - gs%eigen(ib_vb, ik)
                        if (eg_tmp > 0.0d0 .and. eg_tmp < gs%eg_au) then
                            gs%eg_au = eg_tmp
                        end if
                    end do
                end do
            end do
            
            ! Ensure eg_au is positive and reasonable
            if (gs%eg_au < 1.0d-6) gs%eg_au = 1.0d-6
            
            if (irank == 0) then
                write(*,'("# info: auto-calculated minimum band gap =",ES12.5," au (=",ES12.5," eV)")') &
                    gs%eg_au, gs%eg_au * au_ev
            end if
            return
        endif
        
        ! Use user-specified value (default: 1.5 eV), convert from eV to atomic units
        gs%eg_au = eg_ev / au_ev
        
        if (irank == 0) then
            write(*,'("# info: using eg_ev =",ES12.5," eV (= ",ES12.5," au)")') eg_ev, gs%eg_au
        end if
    end subroutine calc_eg_au


end subroutine init_sbe_gs_info

end module gs_info_ssbe

