!
!  Non-cubic (primitive-cell, no-folding) EPM geometry + local form factors.
!
!  The cubic solver (epm_cohen_bergstresser + epm_solver's cubic path) builds the
!  zincblende/diamond GaAs/Si band structure on the SIMPLE-CUBIC 8-atom supercell
!  with the FCC-in-cubic parity selection rule. That convention is hardwired and
!  verified byte-equivalent to the Python reference; it is left untouched.
!
!  This module adds the genuinely NON-orthogonal primitive cells (no folding),
!  starting with monolayer GRAPHENE -- the 2-atom hexagonal honeycomb embedded in
!  a 3D vacuum box. It mirrors the validated Python reference epm_graphene.py /
!  epm_graphene_primitive.py:
!    * lattice a1 = a(1,0,0), a2 = a(1/2, sqrt3/2, 0), a3 = (0,0,vacuum);
!    * 2 carbon atoms at +/- delta/2, delta = (a1+a2)/3 (bond center origin, so
!      the structure factor is real -- centrosymmetric, V^A = 0);
!    * minimal Ramanujam pi-model: 3 symmetric form factors V_S(|G|^2) [eV] on
!      the n = 4, 12, 16 shells, n = |G|^2 / (2*pi/(sqrt3*a))^2 [thesis unit].
!  [Ramanujam, M.S. thesis, Arizona State University (2015); the Dirac cone
!   (gapless at K, v_F ~ 9.6e5 m/s) is symmetry-protected by the honeycomb.]
!
!  Unlike the cubic path, the structure factor is the GENERAL multi-atom sum
!  S(dG) = sum_a exp(-i dG . tau_a) (no parity rule), and the plane-wave basis is
!  built directly from the Cartesian |G|^2 cutoff -- both work for any cell.
!
module epm_noncubic
    use math_constants, only: pi
    implicit none

    private
    public :: nc_is_noncubic, nc_lattice_and_basis, nc_form_factor

    ! eV <-> Ha and Angstrom <-> Bohr (kept local so the module is self-contained)
    real(8), parameter :: HA_TO_EV   = 27.211386245988d0
    real(8), parameter :: ANG_TO_BOHR = 1.0d0 / 0.52917721067d0

    ! --- graphene constants (Python reference epm_graphene.py) -----------------
    real(8), parameter :: GR_A_ANG     = 2.46d0          ! lattice constant [Ang]
    real(8), parameter :: GR_VACUUM_ANG = 20.0d0         ! c-axis vacuum [Ang]
    ! Ramanujam in-plane monolayer form factors [eV] on the n=4,12,16 shells
    ! (n = |G|^2 in units of (2*pi/(sqrt3*a))^2; the thesis sqrt3*a unit).
    integer, parameter :: GR_NSHELL = 3
    integer, parameter :: GR_SHELL(GR_NSHELL) = (/ 4, 12, 16 /)
    real(8), parameter :: GR_VFF_EV(GR_NSHELL) = (/ -8.23d0, 1.5d0, 0.05d0 /)

contains

    ! Is this material handled by the non-cubic primitive path?
    pure logical function nc_is_noncubic(material)
        character(*), intent(in) :: material
        select case (trim(material))
        case ('graphene')
            nc_is_noncubic = .true.
        case default
            nc_is_noncubic = .false.
        end select
    end function nc_is_noncubic

    ! Lattice vectors a1,a2,a3 [Bohr] (columns) and the basis atom positions
    ! tau_atoms(1:3,1:natom) [Bohr]. a_in is the requested in-plane lattice
    ! constant [Bohr]; if <= 0 the cited material default is used.
    subroutine nc_lattice_and_basis(material, a_in, a_matrix, natom, tau_atoms, a_used)
        character(*), intent(in)  :: material
        real(8),      intent(in)  :: a_in
        real(8),      intent(out) :: a_matrix(3,3)        ! columns a1,a2,a3
        integer,      intent(out) :: natom
        real(8), allocatable, intent(out) :: tau_atoms(:,:)
        real(8),      intent(out) :: a_used               ! in-plane a actually used [Bohr]
        real(8) :: a, vac, delta(3)

        select case (trim(material))
        case ('graphene')
            ! graphene's in-plane lattice constant is a FIXED cited material
            ! constant (2.46 Ang); always use it (the namelist epm_lattice_constant_au
            ! defaults to the GaAs value, so honouring it here would be a trap).
            ! Matches epm_graphene.py A_LATT = 2.46 Ang exactly (interchangeable).
            a   = GR_A_ANG * ANG_TO_BOHR
            vac = GR_VACUUM_ANG * ANG_TO_BOHR
            a_matrix(1:3,1) = (/ a,         0d0,              0d0 /)
            a_matrix(1:3,2) = (/ 0.5d0*a,   sqrt(3d0)/2d0*a,  0d0 /)
            a_matrix(1:3,3) = (/ 0d0,       0d0,              vac /)
            natom = 2
            allocate(tau_atoms(3, natom))
            ! delta = (a1+a2)/3 (A->B bond), atoms at +/- delta/2
            delta(1:3) = (a_matrix(1:3,1) + a_matrix(1:3,2)) / 3d0
            tau_atoms(1:3,1) = -0.5d0 * delta(1:3)
            tau_atoms(1:3,2) = +0.5d0 * delta(1:3)
            a_used = a
        case default
            stop 'epm_noncubic: unsupported non-cubic material'
        end select
    end subroutine nc_lattice_and_basis

    ! Symmetric local form factor V_S(|dG|^2) [Ha] for the non-cubic material,
    ! given the Cartesian |dG|^2 [Bohr^-2] and the in-plane lattice constant
    ! a_used [Bohr]. Nonzero only on the cited shells (rounded shell index within
    ! tol). V^A = 0 (centrosymmetric) for graphene.
    pure function nc_form_factor(material, dg2_au, a_used) result(vs_ha)
        character(*), intent(in) :: material
        real(8),      intent(in) :: dg2_au, a_used
        real(8) :: vs_ha
        real(8) :: unit, xn
        integer :: nr, s
        real(8), parameter :: TOL = 0.2d0

        vs_ha = 0d0
        select case (trim(material))
        case ('graphene')
            ! thesis shell unit = (2*pi/(sqrt3*a))^2
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
        end select
    end function nc_form_factor

end module epm_noncubic
