!
!  Nonlocal "super-compute" mode (Part C) -- rate and final-state-search
!  PRIMITIVES. This module contains only PURE functions and cited data tables;
!  it has no SALMON dependencies and no state, so it can be unit-tested on its
!  own (tests/test_superres_rates.f90) and reused by the dynamics integration in
!  later increments. Nothing here is called from the time loop yet -- the whole
!  nonlocal mode is gated behind yn_sbe_superres='n' (default).
!
!  Contents:
!    - nu_saturation        : collision-rate saturation nu(eps) (smooth, no min())
!    - bose_factor          : phonon occupation N_B(hw, T)
!    - gaussian_bin/rect_bin: normalized broadened delta(Delta E) (energy bins, C3)
!    - frohlich_hi_factor   : GaAs polar-optical high-energy asinh factor
!    - ii_rate_general      : general impact-ionization fit P (eps-E_th)^a
!    - bgr_gap_shift_ev     : bandgap-renormalization cube-root law (C7)
!    - Si / GaAs intervalley deformation-potential data tables (C5)
!
!  All energies in a single consistent unit chosen by the caller (the dynamics
!  uses Hartree atomic units); the data tables carry their own stated units.
!
module sbe_superres_ssbe
    implicit none
    private

    real(8), parameter :: PI = 3.14159265358979323846d0

    public :: nu_saturation, bose_factor, gaussian_bin, rect_bin, &
              frohlich_hi_factor, ii_rate_general, bgr_gap_shift_ev

    ! =====================================================================
    ! Silicon intervalley deformation potentials -- Pop "new" set (default).
    ! Order: g-TA, g-LA, g-LO, f-TA, f-LA, f-TO. E [meV], D [1e8 eV/cm].
    ! g-type couple same-<100>-axis valleys; f-type orthogonal-axis valleys.
    ! [Jacoboni & Reggiani, Rev. Mod. Phys. 55, 645 (1983); Pop set compiled in
    !  Jacoboni-Lugli; anchor Canali et al., Phys. Rev. B 15, 3994 (1977)]
    ! =====================================================================
    integer, parameter, public :: SI_N_PHONON = 6
    real(8), parameter, public :: SI_PHONON_E_MEV(SI_N_PHONON) = &
        (/ 10.0d0, 19.0d0, 63.0d0, 19.0d0, 51.0d0, 57.0d0 /)
    real(8), parameter, public :: SI_PHONON_D_1E8EVCM(SI_N_PHONON) = &
        (/  0.3d0,  1.5d0,  6.0d0,  0.5d0,  3.5d0,  1.5d0 /)
    ! Silicon acoustic/transport constants [Jacoboni-Reggiani 1983; Canali 1977]
    real(8), parameter, public :: SI_XI_D_EV     = 9.0d0      ! acoustic def. potential [eV]
    real(8), parameter, public :: SI_RHO_GCM3    = 2.33d0     ! mass density [g/cm^3]
    real(8), parameter, public :: SI_VLA_CMS     = 9.01d5     ! v_LA [cm/s]
    real(8), parameter, public :: SI_VTA_CMS     = 5.23d5     ! v_TA [cm/s]
    real(8), parameter, public :: SI_ALPHA_NP    = 0.5d0      ! non-parabolicity [1/eV]

    ! =====================================================================
    ! GaAs intervalley deformation potentials. Order: Gamma-L, Gamma-X, L-L,
    ! L-X, X-X. E [meV], D [eV/Angstrom]. Plus Frohlich polar-optical constants.
    ! [Fischetti & Laux, Phys. Rev. B 38, 9721 (1988); IEEE TED 38, 634 (1991);
    !  Frohlich: Fawcett-Boardman-Swain, J. Phys. Chem. Solids 31, 1963 (1970)]
    ! =====================================================================
    integer, parameter, public :: GAAS_N_IV = 5
    real(8), parameter, public :: GAAS_IV_E_MEV(GAAS_N_IV) = &
        (/ 27.8d0, 29.9d0, 29.0d0, 29.3d0, 29.9d0 /)
    real(8), parameter, public :: GAAS_IV_D_EVANG(GAAS_N_IV) = &
        (/ 10.0d0, 10.0d0, 10.0d0,  5.0d0,  7.0d0 /)
    real(8), parameter, public :: GAAS_HW_LO_MEV = 36.0d0    ! Frohlich LO phonon [meV]
    real(8), parameter, public :: GAAS_ALPHA_FR  = 0.068d0   ! Frohlich coupling
    real(8), parameter, public :: GAAS_EPS0      = 12.9d0    ! static dielectric
    real(8), parameter, public :: GAAS_EPS_INF   = 10.89d0   ! high-freq dielectric
    real(8), parameter, public :: GAAS_M_GAMMA   = 0.067d0   ! Gamma effective mass [m_e]

contains

    ! Collision-rate saturation: nu(eps) = nu_sat [1 - exp(-(eps/eps0)^n)].
    ! Smooth (no hard min() cutoff -- the derivative discontinuity destabilizes
    ! the stiff solver). Rises through ~0.5-1.5 eV and saturates at nu_sat.
    ! [Meng et al., PRB 91, 075201 (2015); Fischetti & Laux, PRB 38, 9721 (1988)]
    pure function nu_saturation(eps, nu_sat, eps0, nexp) result(nu)
        real(8), intent(in) :: eps, nu_sat, eps0, nexp
        real(8) :: nu, x
        if (eps <= 0d0 .or. eps0 <= 0d0) then
            nu = 0d0
        else
            x = (eps / eps0) ** nexp
            nu = nu_sat * (1d0 - exp(-x))
        end if
    end function nu_saturation

    ! Bose-Einstein phonon occupation N_B = 1/(exp(hw/kT) - 1). hw and kT in the
    ! same units. Guards the hw/kT -> 0 (classical) and large-argument limits.
    pure function bose_factor(hw, kT) result(nb)
        real(8), intent(in) :: hw, kT
        real(8) :: nb, x
        if (kT <= 0d0) then
            nb = 0d0                       ! T=0: no absorption
        else
            x = hw / kT
            if (x < 1d-8) then
                nb = kT / max(hw, 1d-300)  ! classical limit N_B ~ kT/hw
            else if (x > 7.0d2) then
                nb = 0d0                   ! exp overflow guard
            else
                nb = 1d0 / (exp(x) - 1d0)
            end if
        end if
    end function bose_factor

    ! Normalized Gaussian energy bin (broadened delta): area 1, width sigma.
    ! delta(Delta E) -> (1/(sqrt(2 pi) sigma)) exp(-dE^2 / (2 sigma^2)).
    pure function gaussian_bin(dE, sigma) result(w)
        real(8), intent(in) :: dE, sigma
        real(8) :: w, z
        if (sigma <= 0d0) then
            w = 0d0
        else
            z = dE / sigma
            w = exp(-0.5d0 * z * z) / (sqrt(2d0 * PI) * sigma)
        end if
    end function gaussian_bin

    ! Normalized unit-area rectangle energy bin: 1/width for |dE| < width/2, else 0.
    pure function rect_bin(dE, width) result(w)
        real(8), intent(in) :: dE, width
        real(8) :: w
        if (width <= 0d0) then
            w = 0d0
        else if (abs(dE) < 0.5d0 * width) then
            w = 1d0 / width
        else
            w = 0d0
        end if
    end function rect_bin

    ! GaAs polar-optical (Frohlich) EMISSION high-energy factor:
    !   f(E) = (1/sqrt(E)) asinh(sqrt(E/hw0 - 1))   for E > hw0, else 0.
    ! Decays slower than 1/sqrt(E); omitting the asinh under-counts ~20-30% at
    ! 2-3 eV. asinh(y) = log(y + sqrt(y^2+1)). [Fawcett-Boardman-Swain 1970]
    pure function frohlich_hi_factor(E, hw0) result(f)
        real(8), intent(in) :: E, hw0
        real(8) :: f, y
        if (E <= hw0 .or. hw0 <= 0d0) then
            f = 0d0
        else
            y = sqrt(E / hw0 - 1d0)
            f = log(y + sqrt(y * y + 1d0)) / sqrt(E)
        end if
    end function frohlich_hi_factor

    ! General impact-ionization fit rate gamma = P (eps_kin - E_th)^a, gated at
    ! threshold (Theta). a=4 GaAs Stobbe (hard), a=2 Si Keldysh (soft), a=4.6 Si
    ! full-band. All consistent units (caller supplies P in 1/(energy^a time)).
    ! [Stobbe-Redmer-Schattke PRB 49, 4494; Keldysh JETP 21, 1135; Kamakura JAP 75, 3500]
    pure function ii_rate_general(eps_kin, E_th, P, a) result(g)
        real(8), intent(in) :: eps_kin, E_th, P, a
        real(8) :: g, d
        d = eps_kin - E_th
        if (d <= 0d0) then
            g = 0d0
        else
            g = P * d ** a
        end if
    end function ii_rate_general

    ! Electron-hole-plasma bandgap-renormalization cube-root law (C7):
    !   Delta E_gap [eV] = -K * (n [cm^-3])^(1/3),  K ~ 1.9e-8 eV.cm.
    ! Returns a NEGATIVE shift (gap shrinks). [Vashishta & Kalia, PRB 25, 6492 (1982)]
    pure function bgr_gap_shift_ev(n_cm3, K) result(dE)
        real(8), intent(in) :: n_cm3, K
        real(8) :: dE
        if (n_cm3 <= 0d0) then
            dE = 0d0
        else
            dE = -K * n_cm3 ** (1d0 / 3d0)
        end if
    end function bgr_gap_shift_ev

end module sbe_superres_ssbe
