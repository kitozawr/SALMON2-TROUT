!
!  Non-cubic (primitive-cell, no-folding) EPM geometry + local form factors.
!
!  The cubic GaAs/Si path (epm_cohen_bergstresser + epm_solver's cubic builder,
!  reused on the single-parity basis for the FCC primitive) is left untouched.
!  This module supplies the genuinely NON-orthogonal primitive cells with a
!  GENERAL multi-atom structure factor:
!
!    * GRAPHENE -- 2-atom hexagonal honeycomb embedded in vacuum (V^A = 0,
!      Ramanujam minimal pi-EPM, 3 symmetric form factors). Mirrors
!      epm_graphene.py / epm_graphene_primitive.py.
!      [Ramanujam, M.S. thesis, Arizona State University (2015).]
!    * CdS WURTZITE -- 4-atom hexagonal primitive (2 Cd + 2 S), polar
!      (inversion broken => V^A != 0). Symmetric/antisymmetric form factors
!      interpolated from the cited BC1967 Table II anchors. Mirrors
!      epm_wurtzite_cds.py / epm_cds_primitive.py.
!      [Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967).]
!
!  General Hamiltonian (Hartree):
!    H_ij = 0.5|k+G_i|^2 delta_ij + S^S(dG) V^S(|dG|^2) + S^A(dG) V^A(|dG|^2)
!    S^S(dG) = (1/norm) sum_a       exp(-i dG.tau_a)
!    S^A(dG) = (1/norm) sum_a P_a   exp(-i dG.tau_a)     (P_a = +1 cation, -1 anion)
!  with the per-material structure-factor normalisation `norm` (graphene 1 --
!  the eV form factors bake in the 2-atom sum; CdS = total atoms n = 4, the
!  BC1967 per-atom volume normalisation). For graphene V^A = 0 and all P_a = +1.
!
module epm_noncubic
    use math_constants, only: pi
    implicit none

    private
    public :: nc_is_noncubic, nc_lattice_and_basis, nc_form_factors

    real(8), parameter :: HA_TO_EV    = 27.211386245988d0
    real(8), parameter :: RY_TO_HA    = 0.5d0
    real(8), parameter :: ANG_TO_BOHR = 1.0d0 / 0.52917721067d0

    ! --- graphene (Ramanujam) --------------------------------------------------
    real(8), parameter :: GR_A_ANG      = 2.46d0          ! lattice constant [Ang]
    real(8), parameter :: GR_VACUUM_ANG = 20.0d0          ! c-axis vacuum [Ang]
    integer, parameter :: GR_NSHELL = 3
    integer, parameter :: GR_SHELL(GR_NSHELL) = (/ 4, 12, 16 /)
    real(8), parameter :: GR_VFF_EV(GR_NSHELL) = (/ -8.23d0, 1.5d0, 0.05d0 /)

    ! --- CdS wurtzite (BC1967 Table I/II) --------------------------------------
    real(8), parameter :: CDS_A_ANG = 4.136d0
    real(8), parameter :: CDS_C_ANG = 4.136d0 * 1.623d0   ! c/a = 1.623
    real(8), parameter :: CDS_U     = 3.0d0 / 8.0d0       ! internal parameter
    ! V^S / V^A anchors [Ry] keyed by the reduced shell number n (physical |G|^2
    ! = n*(2*pi/a_ZB)^2, a_ZB = sqrt2*a_W). Linearly interpolated; 0 beyond n=16.
    integer, parameter :: CDS_NVS = 8
    real(8), parameter :: CDS_VS_N(CDS_NVS)  = (/ 0d0, 3.04d0, 3.43d0, 5.70d0, 9.50d0, 10.67d0, 13.30d0, 16.0d0 /)
    real(8), parameter :: CDS_VS_V(CDS_NVS)  = (/ 0d0, -0.26d0, -0.24d0, -0.20d0, 0.04d0, 0.04d0, 0.02d0, 0d0 /)
    integer, parameter :: CDS_NVA = 9
    real(8), parameter :: CDS_VA_N(CDS_NVA)  = (/ 0d0, 3.04d0, 3.43d0, 5.70d0, 9.50d0, 11.40d0, 12.15d0, 13.30d0, 16.0d0 /)
    real(8), parameter :: CDS_VA_V(CDS_NVA)  = (/ 0d0, 0.23d0, 0.18d0, 0.08d0, 0.05d0, 0.05d0, 0.05d0, 0.03d0, 0d0 /)

contains

    pure logical function nc_is_noncubic(material)
        character(*), intent(in) :: material
        select case (trim(material))
        case ('graphene', 'CdS')
            nc_is_noncubic = .true.
        case default
            nc_is_noncubic = .false.
        end select
    end function nc_is_noncubic

    ! Lattice vectors a_matrix (columns a1,a2,a3) [Bohr], basis atom positions
    ! tau_atoms(3,natom) [Bohr] and species P_a (+1 cation / -1 anion), the
    ! structure-factor normalisation `snorm`, and the in-plane lattice constant
    ! a_used [Bohr] used by the form factors. a_in (epm_lattice_constant_au) is
    ! IGNORED -- both materials have fixed cited lattice constants.
    subroutine nc_lattice_and_basis(material, a_in, a_matrix, natom, tau_atoms, &
                                    spec, snorm, a_used)
        character(*), intent(in)  :: material
        real(8),      intent(in)  :: a_in
        real(8),      intent(out) :: a_matrix(3,3)
        integer,      intent(out) :: natom
        real(8), allocatable, intent(out) :: tau_atoms(:,:)
        real(8), allocatable, intent(out) :: spec(:)
        real(8),      intent(out) :: snorm
        real(8),      intent(out) :: a_used
        real(8) :: a, c, vac, delta(3), frac(3,4)
        integer :: ia

        select case (trim(material))
        case ('graphene')
            ! fixed cited a = 2.46 Ang (the namelist default is the GaAs value).
            a   = GR_A_ANG * ANG_TO_BOHR
            vac = GR_VACUUM_ANG * ANG_TO_BOHR
            a_matrix(1:3,1) = (/ a,        0d0,             0d0 /)
            a_matrix(1:3,2) = (/ 0.5d0*a,  sqrt(3d0)/2d0*a, 0d0 /)
            a_matrix(1:3,3) = (/ 0d0,      0d0,             vac /)
            natom = 2
            allocate(tau_atoms(3, natom), spec(natom))
            delta(1:3) = (a_matrix(1:3,1) + a_matrix(1:3,2)) / 3d0   ! A->B bond
            tau_atoms(1:3,1) = -0.5d0 * delta(1:3)
            tau_atoms(1:3,2) = +0.5d0 * delta(1:3)
            spec(1:2) = (/ 1d0, 1d0 /)       ! single species -> V^A inert
            snorm  = 1d0                      ! eV form factors bake in the 2-atom sum
            a_used = a
        case ('CdS')
            ! wurtzite hexagonal primitive: a1=a(1,0,0), a2=a(-1/2,sqrt3/2,0),
            ! a3=c(0,0,1). 4 atoms: Cd(0,0,0),(1/3,2/3,1/2); S +(0,0,u).
            a = CDS_A_ANG * ANG_TO_BOHR
            c = CDS_C_ANG * ANG_TO_BOHR
            a_matrix(1:3,1) = (/ a,        0d0,             0d0 /)
            a_matrix(1:3,2) = (/ -0.5d0*a, sqrt(3d0)/2d0*a, 0d0 /)
            a_matrix(1:3,3) = (/ 0d0,      0d0,             c   /)
            natom = 4
            allocate(tau_atoms(3, natom), spec(natom))
            frac(1:3,1) = (/ 0d0,       0d0,       0d0           /)   ! Cd
            frac(1:3,2) = (/ 1d0/3d0,   2d0/3d0,   0.5d0         /)   ! Cd
            frac(1:3,3) = (/ 0d0,       0d0,       CDS_U         /)   ! S
            frac(1:3,4) = (/ 1d0/3d0,   2d0/3d0,   0.5d0 + CDS_U /)   ! S
            spec(1:4) = (/ 1d0, 1d0, -1d0, -1d0 /)
            do ia = 1, natom
                tau_atoms(1:3,ia) = frac(1,ia)*a_matrix(1:3,1) &
                                  + frac(2,ia)*a_matrix(1:3,2) + frac(3,ia)*a_matrix(1:3,3)
            end do
            snorm  = dble(natom)             ! BC1967 per-atom (1/n) normalisation
            a_used = a
        case default
            stop 'epm_noncubic: unsupported non-cubic material'
        end select
    end subroutine nc_lattice_and_basis

    ! Symmetric/antisymmetric local form factors V^S, V^A [Ha] given the Cartesian
    ! |dG|^2 [Bohr^-2] and the in-plane lattice constant a_used [Bohr].
    pure subroutine nc_form_factors(material, dg2_au, a_used, vs_ha, va_ha)
        character(*), intent(in)  :: material
        real(8),      intent(in)  :: dg2_au, a_used
        real(8),      intent(out) :: vs_ha, va_ha
        real(8) :: unit, xn, a_zb
        integer :: nr, s
        real(8), parameter :: TOL = 0.2d0

        vs_ha = 0d0;  va_ha = 0d0
        select case (trim(material))
        case ('graphene')
            ! thesis shell unit (2*pi/(sqrt3*a))^2; nonzero only on n=4,12,16.
            unit = (2d0 * pi / (sqrt(3d0) * a_used)) ** 2
            xn = dg2_au / unit
            nr = nint(xn)
            if (abs(xn - dble(nr)) < TOL) then
                do s = 1, GR_NSHELL
                    if (nr == GR_SHELL(s)) then
                        vs_ha = GR_VFF_EV(s) / HA_TO_EV
                        return
                    end if
                end do
            end if
        case ('CdS')
            ! BC1967 anchors keyed by n*(2*pi/a_ZB)^2, a_ZB = sqrt2*a_W. Linear
            ! interpolation in physical |G|^2 (np.interp; 0 outside [0, 16*unit]).
            a_zb = sqrt(2d0) * a_used
            unit = (2d0 * pi / a_zb) ** 2
            vs_ha = lin_interp(dg2_au, CDS_VS_N, CDS_VS_V, CDS_NVS, unit) * RY_TO_HA
            va_ha = lin_interp(dg2_au, CDS_VA_N, CDS_VA_V, CDS_NVA, unit) * RY_TO_HA
        end select
    end subroutine nc_form_factors

    ! np.interp(x, xp*unit, fp, left=0, right=0): xp increasing (reduced shells),
    ! the anchor x-values are xp(i)*unit [Bohr^-2].
    pure function lin_interp(x, xp, fp, n, unit) result(y)
        real(8), intent(in) :: x, xp(:), fp(:), unit
        integer, intent(in) :: n
        real(8) :: y, x0, x1
        integer :: i
        y = 0d0
        if (x < xp(1)*unit .or. x > xp(n)*unit) return
        do i = 1, n-1
            x0 = xp(i)*unit;  x1 = xp(i+1)*unit
            if (x >= x0 .and. x <= x1) then
                if (x1 > x0) then
                    y = fp(i) + (fp(i+1)-fp(i)) * (x - x0) / (x1 - x0)
                else
                    y = fp(i)
                end if
                return
            end if
        end do
    end function lin_interp

end module epm_noncubic
