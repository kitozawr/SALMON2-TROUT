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
    ! Unit conversions to Hartree atomic units (a.u.). [CODATA 2018]
    real(8), parameter :: HA_EV    = 27.211386245988d0   ! 1 Hartree in eV
    real(8), parameter :: BOHR_CM  = 0.529177210903d-8   ! 1 Bohr in cm
    real(8), parameter :: BOHR_ANG = 0.529177210903d0    ! 1 Bohr in Angstrom
    real(8), parameter :: ME_G     = 9.1093837015d-28    ! electron mass in g

    public :: nu_saturation, bose_factor, gaussian_bin, rect_bin, &
              frohlich_hi_factor, ii_rate_general, bgr_gap_shift_ev, &
              gaussian_shape, amp_damp_channel, &
              mev_to_ha, d_evcm_to_au, d_evang_to_au, rho_gcm3_to_au, &
              golden_rule_prefactor, eph_thermal_split, &
              eps_thomas_fermi, tf_kappa2_degenerate, debye_kappa2, &
              lindhard_F, eps_lindhard_static, plasmon_freq2, lopc_branches

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
    real(8), parameter, public :: GAAS_HW_TO_MEV = 33.6d0    ! TO phonon [meV] (LOPC)
    real(8), parameter, public :: GAAS_ALPHA_FR  = 0.068d0   ! Frohlich coupling
    real(8), parameter, public :: GAAS_EPS0      = 12.9d0    ! static dielectric
    real(8), parameter, public :: GAAS_EPS_INF   = 10.89d0   ! high-freq dielectric
    real(8), parameter, public :: GAAS_M_GAMMA   = 0.067d0   ! Gamma effective mass [m_e]
    real(8), parameter, public :: GAAS_M_HH      = 0.51d0    ! heavy-hole mass [m_e]
    real(8), parameter, public :: GAAS_M_LH      = 0.082d0   ! light-hole mass [m_e]

    ! =====================================================================
    ! Silicon screening constants (non-polar: no LO-phonon-plasmon coupling).
    ! [eps: std Si; masses: Si effective-mass tables (m_l/m_t, hh/lh)]
    ! =====================================================================
    real(8), parameter, public :: SI_EPS    = 11.7d0    ! static dielectric (non-polar)
    real(8), parameter, public :: SI_M_L    = 0.98d0    ! longitudinal mass [m_e]
    real(8), parameter, public :: SI_M_T    = 0.19d0    ! transverse mass [m_e]
    real(8), parameter, public :: SI_M_HH   = 0.49d0    ! heavy-hole mass [m_e]
    real(8), parameter, public :: SI_M_LH   = 0.16d0    ! light-hole mass [m_e]

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

    ! Peak-normalized Gaussian energy-conservation weight (dimensionless, in
    ! [0,1]): exp(-dE^2 / (2 sigma^2)). Used to weight a discrete partner
    ! transition by how well it conserves energy (vs gaussian_bin, which is
    ! area-normalized for converting a state sum into a rate).
    pure function gaussian_shape(dE, sigma) result(w)
        real(8), intent(in) :: dE, sigma
        real(8) :: w, z
        if (sigma <= 0d0) then
            w = merge(1d0, 0d0, dE == 0d0)
        else
            z = dE / sigma
            w = exp(-0.5d0 * z * z)
        end if
    end function gaussian_shape

    ! Exact finite-time amplitude-damping (population-transfer) CPTP map for a
    ! single GKLS jump L = sqrt(gamma) |dst><src| in the basis where src,dst are
    ! basis states (here the Houston/adiabatic basis):
    !   rho_ss -> e^{-g} rho_ss ;  rho_dd -> rho_dd + (1 - e^{-g}) rho_ss
    !   rho_sb -> e^{-g/2} rho_sb  and  rho_bs -> e^{-g/2} rho_bs  (b /= s)
    ! with g = gamma*tau. Trace-preserving and completely positive for any
    ! gamma,tau >= 0 (a genuine GKLS map). This is the building block of the
    ! population-relaxing electron-phonon channel; the impact-ionization channel
    ! uses an identical (separately validated) copy in bloch_solver_ssbe.
    ! [Lindblad, Commun. Math. Phys. 48, 119 (1976)]
    subroutine amp_damp_channel(nba, rho, i_src, i_dst, gamma, tau)
        integer,    intent(in)    :: nba, i_src, i_dst
        complex(8), intent(inout) :: rho(nba, nba)
        real(8),    intent(in)    :: gamma, tau
        real(8) :: g, gh, transfer
        integer :: b
        if (i_src == i_dst) return
        if (gamma * tau < 1d-14) return
        g  = exp(-gamma * tau)
        gh = sqrt(g)
        transfer = (1d0 - g) * real(rho(i_src, i_src))
        do b = 1, nba
            if (b == i_src) cycle
            rho(i_src, b) = gh * rho(i_src, b)
            rho(b, i_src) = gh * rho(b, i_src)
        end do
        rho(i_src, i_src) = g * rho(i_src, i_src)
        rho(i_dst, i_dst) = rho(i_dst, i_dst) + transfer
    end subroutine amp_damp_channel

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

    ! --- unit conversions to a.u. (for golden-rule e-ph rates, C5) -----------
    pure function mev_to_ha(E_mev) result(E)
        real(8), intent(in) :: E_mev
        real(8) :: E
        E = E_mev * 1d-3 / HA_EV
    end function mev_to_ha

    ! Deformation potential D [eV/cm] -> [Ha/Bohr]. (Si tables are in 1e8 eV/cm,
    ! so pass D_table*1e8.)
    pure function d_evcm_to_au(D_evcm) result(D)
        real(8), intent(in) :: D_evcm
        real(8) :: D
        D = D_evcm * BOHR_CM / HA_EV
    end function d_evcm_to_au

    ! Deformation potential D [eV/Angstrom] -> [Ha/Bohr]. (GaAs intervalley tables.)
    pure function d_evang_to_au(D_eva) result(D)
        real(8), intent(in) :: D_eva
        real(8) :: D
        D = D_eva * BOHR_ANG / HA_EV
    end function d_evang_to_au

    ! Mass density rho [g/cm^3] -> [m_e/Bohr^3].
    pure function rho_gcm3_to_au(rho) result(r)
        real(8), intent(in) :: rho
        real(8) :: r
        r = rho * BOHR_CM**3 / ME_G
    end function rho_gcm3_to_au

    ! Golden-rule deformation-potential rate PREFACTOR pi D^2/(rho omega), all
    ! arguments in a.u. (rate per unit final-state density-of-states; the DOS /
    ! energy-bin factor is applied by the caller). [Jacoboni-Reggiani RMP 55, 645]
    pure function golden_rule_prefactor(D_au, rho_au, omega_au) result(g)
        real(8), intent(in) :: D_au, rho_au, omega_au
        real(8) :: g
        if (rho_au <= 0d0 .or. omega_au <= 0d0) then
            g = 0d0
        else
            g = PI * D_au * D_au / (rho_au * omega_au)
        end if
    end function golden_rule_prefactor

    ! Normalized thermal emission/absorption split for one phonon mode:
    !   fe = (N_B+1)/(2 N_B+1),  fa = N_B/(2 N_B+1).
    ! fe+fa = 1 (so the per-mode total is the mode weight) and fe/fa = (N_B+1)/N_B
    ! (detailed balance). At N_B=0: fe=1, fa=0 (spontaneous emission only).
    pure subroutine eph_thermal_split(Nb, fe, fa)
        real(8), intent(in)  :: Nb
        real(8), intent(out) :: fe, fa
        real(8) :: denom
        denom = 2d0 * Nb + 1d0
        fe = (Nb + 1d0) / denom
        fa = Nb / denom
    end subroutine eph_thermal_split

    ! =====================================================================
    ! Part G -- dielectric screening models for the carrier-carrier channel.
    ! W(q) = V(q)/eps(q[,omega]). Pure functions; the caller supplies the
    ! carrier density, Fermi wavevector and phonon frequencies in a.u.
    ! =====================================================================

    ! Static dielectric, model (a) Thomas-Fermi/Debye: eps(q) = 1 + kappa^2/q^2.
    ! [Ashcroft & Mermin; long-wavelength limit of Lindhard theory]
    pure function eps_thomas_fermi(q, kappa2) result(eps)
        real(8), intent(in) :: q, kappa2
        real(8) :: eps
        if (q <= 0d0) then
            eps = huge(1d0)            ! q->0: perfect screening
        else
            eps = 1d0 + kappa2 / (q * q)
        end if
    end function eps_thomas_fermi

    ! Degenerate (T->0) Thomas-Fermi screening wavevector squared [a.u.]:
    !   kappa_TF^2 = 4 (3 n / pi)^(1/3) / eps  ( = 4 k_F/(pi eps) ).
    ! [arXiv:2312.13059; Ashcroft & Mermin]
    pure function tf_kappa2_degenerate(n_au, eps_bg) result(k2)
        real(8), intent(in) :: n_au, eps_bg
        real(8) :: k2
        if (n_au <= 0d0) then
            k2 = 0d0
        else
            k2 = 4d0 * (3d0 * n_au / PI) ** (1d0 / 3d0) / max(eps_bg, 1d-12)
        end if
    end function tf_kappa2_degenerate

    ! Nondegenerate (Debye-Huckel) screening wavevector squared [a.u.]:
    !   kappa_D^2 = 4 pi n / (eps k_B T).  [Ashcroft & Mermin]
    pure function debye_kappa2(n_au, eps_bg, kT_au) result(k2)
        real(8), intent(in) :: n_au, eps_bg, kT_au
        real(8) :: k2
        if (n_au <= 0d0 .or. kT_au <= 0d0) then
            k2 = 0d0
        else
            k2 = 4d0 * PI * n_au / (max(eps_bg, 1d-12) * kT_au)
        end if
    end function debye_kappa2

    ! 3D static Lindhard function F(x), x = q/(2 k_F):
    !   F(x) = 1/2 + (1-x^2)/(4x) ln|(1+x)/(1-x)|.
    ! F(0)=1 (recovers Thomas-Fermi), F(1)=1/2 (the 2k_F kink -> Friedel
    ! oscillations). [J. Lindhard, Mat.-Fys. Medd. 28, 8 (1954); Ashcroft & Mermin]
    pure function lindhard_F(x) result(F)
        real(8), intent(in) :: x
        real(8) :: F
        if (abs(x) < 1d-8) then
            F = 1d0
        else if (abs(x - 1d0) < 1d-8) then
            F = 0.5d0
        else
            F = 0.5d0 + (1d0 - x * x) / (4d0 * x) * log(abs((1d0 + x) / (1d0 - x)))
        end if
    end function lindhard_F

    ! Static Lindhard/RPA dielectric, model (b, recommended default):
    !   eps(q,0) = 1 + (kappa_TF^2/q^2) F(q/2k_F).
    ! Correct across the q range that dominates the collision integral (Thomas-
    ! Fermi over-screens for q>2k_F). [Lindhard 1954; arXiv:1206.2003]
    pure function eps_lindhard_static(q, kF, kappa2) result(eps)
        real(8), intent(in) :: q, kF, kappa2
        real(8) :: eps
        if (q <= 0d0) then
            eps = huge(1d0)
        else if (kF <= 0d0) then
            eps = 1d0 + kappa2 / (q * q)        ! no carriers: bare TF form
        else
            eps = 1d0 + kappa2 / (q * q) * lindhard_F(q / (2d0 * kF))
        end if
    end function eps_lindhard_static

    ! Bulk plasmon frequency squared [a.u.]: omega_p^2 = 4 pi n/(eps_inf m*).
    ! (the high-frequency plasmon screens with eps_inf, not eps_0). [Mahan]
    pure function plasmon_freq2(n_au, eps_inf, mstar) result(wp2)
        real(8), intent(in) :: n_au, eps_inf, mstar
        real(8) :: wp2
        if (n_au <= 0d0) then
            wp2 = 0d0
        else
            wp2 = 4d0 * PI * n_au / (max(eps_inf, 1d-12) * max(mstar, 1d-12))
        end if
    end function plasmon_freq2

    ! Coupled LO-phonon-plasmon (LOPC) branches, model (c, GaAs only):
    !   omega_{L+/-}^2 = 1/2 [ (wp^2+wLO^2) +/- sqrt((wp^2+wLO^2)^2 - 4 wp^2 wTO^2) ].
    ! Vieta: wLp2+wLm2 = wp^2+wLO^2, wLp2*wLm2 = wp^2 wTO^2. The branches
    ! anticross at wp = wLO. [Varga PR 137, A1896; Mooradian-McWhorter PR 177, 1231]
    pure subroutine lopc_branches(wp2, wLO2, wTO2, wLp2, wLm2)
        real(8), intent(in)  :: wp2, wLO2, wTO2
        real(8), intent(out) :: wLp2, wLm2
        real(8) :: b, disc
        b = wp2 + wLO2
        disc = sqrt(max(b * b - 4d0 * wp2 * wTO2, 0d0))
        wLp2 = 0.5d0 * (b + disc)
        wLm2 = 0.5d0 * (b - disc)
    end subroutine lopc_branches

end module sbe_superres_ssbe
