! Ground State Date Storage Module:

module gs_info_ssbe
    use math_constants, only: pi, zI
    use phys_constants, only: au_ev
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
        integer, allocatable :: unfold_sub(:, :)     ! (nb, nk) sublattice 1..4
        integer, allocatable :: unfold_prim(:, :)    ! (nb, nk) primitive band rank
        real(8) :: unfold_offset(1:3, 1:4)           ! G0 in sc reduced coords
    end type


contains


subroutine init_sbe_gs_info(gs, sysname, gs_directory, nk, nb, ne, a1, a2, a3, read_bin, icomm)
    use communication
    use filesystem, only: open_filehandle, get_filehandle
    use salmon_global, only: yn_sbe_spinor
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
    allocate(gs%delta_omega(1:nb, 1:nb, 1:nk))
    allocate(gs%p_mod_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%d_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%p_tm_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%rvnl_tm_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%unfold_sub(1:nb, 1:nk), gs%unfold_prim(1:nb, 1:nk))
    gs%unfold_sub = 0
    gs%unfold_prim = 0
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
        call comm_bcast(gs%unfold_sub, icomm, 0)
        call comm_bcast(gs%unfold_prim, icomm, 0)
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

    ! Calculate minimum band gap in atomic units (for gauge-covariant decoherence)
    call calc_eg_au()

contains

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


    ! Read k-point coordinates from SALMON's output file
    subroutine read_k_data()
        implicit none
        character(256) :: dummy
        integer :: fh, ik, iik
        real(8) :: tmp(4)
        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_k.data', 'old')
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        do ik = 1, nk
            read(fh, *) iik, tmp(1:4)
            if (ik .ne. iik) stop "ik mismatch"
            gs%kpoint(1:3, ik) = tmp(1:3)
            gs%kweight(ik) = tmp(4)
        end do
        close(fh)
    end subroutine read_k_data


    ! Read eigenvalue data from SALMON's output file
    subroutine read_eigen_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, iik, iib, ib
        real(8) :: tmp(2)

        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_eigen.data', 'old')
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        read(fh, "(a)") dummy
        do ik = 1, nk
            read(fh, "(a)") dummy
            do ib = 1, nb
                read(fh, *) iib, tmp(1:2)
                if (ib .ne. iib) stop "ib mismatch"
                gs%eigen(ib, ik) = tmp(1)
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
                    gs%p_tm_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%p_tm_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%p_tm_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
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
                    gs%rvnl_tm_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%rvnl_tm_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%rvnl_tm_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
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
        integer :: fh, ik, ib, iik, iib, isub, ibprim, i, nnk, nnb, ioff(3)
        real(8) :: w

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
        read(fh, *) nnk, nnb, gs%nv_prim
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
        do i = 1, 4
            read(fh, *) isub, ioff(1:3)
            if (isub .ne. i) stop "unfold map: offset index mismatch"
            gs%unfold_offset(1:3, i) = dble(ioff(1:3))
        end do
        read(fh, "(a)") dummy
        do ik = 1, nk
            do ib = 1, nb
                read(fh, *) iik, iib, isub, ibprim, w
                if (ik .ne. iik) stop "unfold map: ik mismatch"
                if (ib .ne. iib) stop "unfold map: ib mismatch"
                gs%unfold_sub(ib, ik) = isub
                gs%unfold_prim(ib, ik) = ibprim
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

