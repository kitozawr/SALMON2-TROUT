!
!  Empirical Pseudopotential Method (EPM): Cohen-Bergstresser local form factors
!
!  Reference: M. L. Cohen and T. K. Bergstresser, Phys. Rev. 141, 789 (1966).
!  Local pseudopotential form factors for GaAs (zincblende), tabulated at the
!  squared reciprocal-lattice-vector magnitudes G^2 = 3, 4, 8, 11 in units of
!  (2*pi/a)^2. Values are quoted by the original authors in Rydberg; this
!  module returns them already converted to Hartree atomic units (the
!  convention used throughout the rest of SALMON), i.e. divided by 2.
!
module epm_cohen_bergstresser
    implicit none

    private
    public :: cb_get_form_factors, cb_tau_zincblende, cb_lattice_vectors_fcc
    public :: cb_load_formfactor_file

    ! --- DFT-fitted form-factor table (epm_material=='file') -----------------
    ! Loaded once (on rank 0, then used everywhere after the standard init-time
    ! broadcast of the lattice/basis is irrelevant here because every rank calls
    ! cb_load_formfactor_file independently from init_epm_info). Shells absent
    ! from the table return (0,0), exactly like the built-in tables.
    integer, parameter :: cb_max_shells = 64
    logical, save      :: cb_file_loaded = .false.
    integer, save      :: cb_file_nff = 0
    integer, save      :: cb_file_g2(cb_max_shells) = 0
    real(8), save      :: cb_file_vs_ry(cb_max_shells) = 0d0
    real(8), save      :: cb_file_va_ry(cb_max_shells) = 0d0

contains

    ! Load a DFT-fitted local form-factor table written by
    ! tools/dft_to_epm/dft_to_epm.py. Lines beginning with '#' are comments; data
    ! lines are "G2  VS_ry  VA_ry" (G2 = |G|^2 in (2*pi/a)^2 units, integer; the
    ! form factors are in Rydberg, converted to Hartree by cb_get_form_factors
    ! exactly like the built-in tables). Called once from init_epm_info when
    ! epm_material=='file'.
    subroutine cb_load_formfactor_file(path)
        use filesystem, only: get_filehandle
        implicit none
        character(*), intent(in) :: path
        integer :: fh, ios, g2
        real(8) :: vs_ry, va_ry
        character(256) :: line

        if (cb_file_loaded) return
        if (len_trim(path) == 0) then
            stop 'epm_cohen_bergstresser: epm_material=="file" requires epm_formfactor_file'
        end if

        cb_file_nff = 0
        fh = get_filehandle()
        open(unit=fh, file=trim(path), action='read', status='old', iostat=ios)
        if (ios /= 0) then
            write(*,'(a,a)') '# EPM: cannot open form-factor file: ', trim(path)
            stop 'epm_cohen_bergstresser: failed to open epm_formfactor_file'
        end if

        do
            read(fh, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#' .or. line(1:1) == '!') cycle
            read(line, *, iostat=ios) g2, vs_ry, va_ry
            if (ios /= 0) cycle
            if (cb_file_nff >= cb_max_shells) then
                stop 'epm_cohen_bergstresser: too many form-factor shells (raise cb_max_shells)'
            end if
            cb_file_nff = cb_file_nff + 1
            cb_file_g2(cb_file_nff)    = g2
            cb_file_vs_ry(cb_file_nff) = vs_ry
            cb_file_va_ry(cb_file_nff) = va_ry
        end do
        close(fh)

        if (cb_file_nff == 0) then
            stop 'epm_cohen_bergstresser: form-factor file contained no usable shells'
        end if
        cb_file_loaded = .true.

        write(*,'(a,i0,a,a)') '# EPM: loaded ', cb_file_nff, &
            & ' DFT-fitted form-factor shells from ', trim(path)
    end subroutine cb_load_formfactor_file

    ! Symmetric/antisymmetric local pseudopotential form factors V^S(G^2), V^A(G^2)
    ! in Hartree atomic units. G2 is the squared length of the reciprocal lattice
    ! vector in units of (2*pi/a)^2 (an integer for the fcc/zincblende lattice).
    subroutine cb_get_form_factors(material, G2, VS_ha, VA_ha)
        implicit none
        character(*), intent(in)  :: material
        integer,      intent(in)  :: G2
        real(8),      intent(out) :: VS_ha, VA_ha
        real(8), parameter :: ry_to_ha = 0.5d0
        real(8) :: VS_ry, VA_ry

        VS_ry = 0d0
        VA_ry = 0d0

        select case (trim(material))
        case ('GaAs')
            ! Cohen-Bergstresser (1966), Table 2:
            !   V^S(3)  = -0.23 Ry,  V^A(3)  = +0.07 Ry
            !   V^S(4)  =  0.00 Ry,  V^A(4)  = +0.05 Ry   (V^S(4)=0: structure factor vanishes)
            !   V^S(8)  = +0.01 Ry,  V^A(8)  =  0.00 Ry
            !   V^S(11) = +0.06 Ry,  V^A(11) = +0.01 Ry
            select case (G2)
            case (3)
                VS_ry = -0.23d0;  VA_ry =  0.07d0
            case (4)
                VS_ry =  0.00d0;  VA_ry =  0.05d0
            case (8)
                VS_ry =  0.01d0;  VA_ry =  0.00d0
            case (11)
                VS_ry =  0.06d0;  VA_ry =  0.01d0
            end select
        case ('Si', 'Si_kunikiyo')
            ! Silicon (diamond): two IDENTICAL atoms per primitive cell, so the
            ! antisymmetric structure factor vanishes -> V^A == 0 for all shells.
            ! Default (production) symmetric form factors, Kunikiyo Table I:
            !   V^S(3) = -0.2258 Ry, V^S(8) = +0.05698 Ry, V^S(11) = +0.070709 Ry
            ! [T. Kunikiyo et al., J. Appl. Phys. 75, 297 (1994), Table I]
            ! Validated: converged indirect gap ~1.06 eV (Kunikiyo's own calc
            ! 1.068 eV; exp 1.12 eV), CBM at 0.86*(2pi/a) along <100>.
            select case (G2)
            case (3)
                VS_ry = -0.2258d0;   VA_ry = 0.0d0
            case (8)
                VS_ry =  0.05698d0;  VA_ry = 0.0d0
            case (11)
                VS_ry =  0.070709d0; VA_ry = 0.0d0
            end select
        case ('Si_cb')
            ! Silicon (diamond), Cohen-Bergstresser (1966) symmetric set
            ! (validation alternative): V^S(3)=-0.21, V^S(8)=+0.04, V^S(11)=+0.08 Ry.
            ! V^A == 0 (diamond). [Cohen & Bergstresser, Phys. Rev. 141, 789 (1966)]
            select case (G2)
            case (3)
                VS_ry = -0.21d0;  VA_ry = 0.0d0
            case (8)
                VS_ry =  0.04d0;  VA_ry = 0.0d0
            case (11)
                VS_ry =  0.08d0;  VA_ry = 0.0d0
            end select
        case ('file')
            ! DFT-fitted local form factors loaded by cb_load_formfactor_file
            ! (tools/dft_to_epm). Shells absent from the table give (0,0), as for
            ! the built-in materials.
            block
                integer :: ish
                if (.not. cb_file_loaded) then
                    stop 'epm_cohen_bergstresser: epm_material=="file" but no table loaded'
                end if
                do ish = 1, cb_file_nff
                    if (cb_file_g2(ish) == G2) then
                        VS_ry = cb_file_vs_ry(ish)
                        VA_ry = cb_file_va_ry(ish)
                        exit
                    end if
                end do
            end block
        case default
            stop 'epm_cohen_bergstresser: unsupported epm_material (use "GaAs", "Si", "Si_cb" or "file")'
        end select

        VS_ha = VS_ry * ry_to_ha
        VA_ha = VA_ry * ry_to_ha
    end subroutine cb_get_form_factors


    ! Zincblende internal displacement between the two-atom basis, with the
    ! origin placed midway between the cation and the anion (Cohen-Bergstresser
    ! convention): tau = (a/8)*(1,1,1). This fixes the sign of V^A consistently
    ! with the tabulated form factors above.
    function cb_tau_zincblende(a_lattice) result(tau)
        implicit none
        real(8), intent(in) :: a_lattice
        real(8) :: tau(3)
        tau(1:3) = a_lattice / 8.0d0
    end function cb_tau_zincblende


    ! Conventional fcc primitive lattice vectors for the zincblende structure
    ! (Cohen-Bergstresser convention):
    !   a1 = (a/2)(0,1,1),  a2 = (a/2)(1,0,1),  a3 = (a/2)(1,1,0)
    subroutine cb_lattice_vectors_fcc(a_lattice, a1, a2, a3)
        implicit none
        real(8), intent(in)  :: a_lattice
        real(8), intent(out) :: a1(3), a2(3), a3(3)
        real(8) :: h
        h = 0.5d0 * a_lattice
        a1 = (/ 0d0, h,   h   /)
        a2 = (/ h,   0d0, h   /)
        a3 = (/ h,   h,   0d0 /)
    end subroutine cb_lattice_vectors_fcc

end module epm_cohen_bergstresser
