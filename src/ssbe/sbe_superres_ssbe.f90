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
              lindhard_F, eps_lindhard_static, plasmon_freq2, lopc_branches, &
              eps_cdrb, interk_vq, &
              energy_partner_weights, fermi_dirac, fit_fermi_dirac, &
              carrier_carrier_relax, eph_interk_dpop, ii_interk_dpop, &
              auger_interk_dpop, mp_grid_triple, mp_partner_triple, &
              vg_eta_admixture, vg_trunc_shift2, vg_conv_error, vg_ptop_exceeds, &
              get_material_params

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

    ! =====================================================================
    ! Wurtzite CdS (P6_3mc) -- polar, NON-centrosymmetric II-VI.
    ! Structure: Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967).
    ! DISSIPATION CHANNELS (cited in the CdS physics-methods spec / md):
    !   * Frohlich polar-optical e-ph -- the PRIMARY room-T channel
    !     (hw_LO = 38 meV [Raman]; alpha ~ 0.5 [cyclotron, arXiv:cond-mat/
    !     0107481; path-integral arXiv:2205.11780]). nu_sat is DERIVED from the
    !     cited coupling: the Frohlich rate scale ~ alpha * omega_LO.
    !   * Coulomb HF (static dielectric eps0 ~ 9.0 isotropic avg, eps_inf ~ 5.3
    !     [md; anisotropic ||/_|_ c]).
    !   * Impact ionization: Keldysh soft threshold E_th ~ 1.5 E_g = 3.6 eV
    !     [md, (3/2)E_g rule]; the PREFACTOR is a fit parameter (CdS-specific
    !     value scarce -- md), so the user must set sbe_ii_prefactor explicitly.
    ! Carrier-carrier (e-e/e-h) and Auger recombination are FORBIDDEN for CdS:
    ! the CdS literature only fixes an e-e *timescale* (sub-100fs thermalization
    ! at n >= 1e18 cm^-3 [Shah et al., IEEE JQE 22, 1728 (1986); Elsaesser PRL 66,
    ! 1757 (1991)]), NOT a carrier-carrier *rate*. There is NO verified CdS Auger
    ! coefficient: the previously-cited "C = 2.0e-30 cm^6/s [Haury, PRB 57, 11513
    ! (1998)]" was a FABRICATED reference -- the real Haury et al. paper is PRL 79,
    ! 511 (1997) on ferromagnetism in CdMnTe quantum wells, unrelated to Auger in
    ! CdS, and the coefficient is unconfirmed anywhere. It has been REMOVED; CdS
    ! Auger is gated off (auger_ok=.false.). A user with a verified rate can still
    ! opt in explicitly via sbe_auger_c_cm6s. (Physically CdS is wide-gap, Eg~2.4
    ! eV, so direct Auger is exponentially suppressed and any future channel would
    ! be the phonon/defect-assisted one -- see wiki/07, the nonlocal-Auger task.)
    ! Piezoelectric acoustic (e33/e31/e15 [Berlincourt PR 129,1009]) and
    ! deformation-potential acoustic are cited but NOT yet SBE Lindblad channels.
    ! =====================================================================
    real(8), parameter, public :: CDS_A_BOHR    = 7.8159d0  ! a = 4.136 Ang  [BC1967 Table I]
    real(8), parameter, public :: CDS_ASQ3_BOHR = 13.5376d0 ! b = a*sqrt(3) (orthorhombic)
    real(8), parameter, public :: CDS_C_BOHR    = 12.6877d0 ! c = 6.714 Ang (c/a=1.623) [BC1967]
    real(8), parameter, public :: CDS_U_INT     = 0.375d0   ! internal parameter u=3/8 [BC1967]
    real(8), parameter, public :: CDS_EG_EV     = 2.58d0    ! direct gap at Gamma (low-T) [BC1967 Table I]
    real(8), parameter, public :: CDS_EPS0      = 8.9d0     ! static dielectric (isotropic avg) [Berlincourt 1963]
    real(8), parameter, public :: CDS_EPS_INF   = 5.3d0     ! high-frequency dielectric [Berlincourt 1963]
    ! (NO verified CdS Auger coefficient exists -- the former CDS_AUGER_C with its
    !  "Haury PRB 57, 11513 (1998)" citation was fabricated and has been removed.)
    real(8), parameter, public :: CDS_EE_ACT_N  = 1.0d18    ! e-e activation density [cm^-3] [Shah JQE 22, 1728]
    real(8), parameter, public :: CDS_HW_LO_MEV = 38.0d0    ! Frohlich LO ~305 cm^-1 [Raman; md]
    real(8), parameter, public :: CDS_ALPHA_FR  = 0.5d0     ! Frohlich coupling alpha [md]
    ! Frohlich e-ph rate scale nu_sat = alpha * omega_LO (omega_LO = hw_LO/hbar):
    ! 0.5 * (38e-3 eV / 6.582e-16 eV.s) = 2.89e13 s^-1 (~35 fs; md: sub-100 fs).
    real(8), parameter, public :: CDS_NU_SAT_SI = 2.89d13   ! derived from cited alpha, hw_LO [md]
    real(8), parameter, public :: CDS_II_ETH_EV = 3.6d0     ! II threshold ~1.5 E_g [md, (3/2)E_g]
    real(8), parameter, public :: CDS_RHO_GCM3  = 4.82d0    ! mass density [g/cm^3]
    real(8), parameter, public :: CDS_E33       = 0.385d0   ! piezo [C/m^2] [Berlincourt PR 129,1009]
    real(8), parameter, public :: CDS_E31       = -0.262d0  ! piezo [C/m^2] [Berlincourt 1963]
    real(8), parameter, public :: CDS_E15       = -0.183d0  ! piezo [C/m^2] [Berlincourt 1963]

    ! =====================================================================
    ! Monolayer graphene -- GAPLESS Dirac semimetal, honeycomb (D6h), V^A=0
    ! (centrosymmetric). a = 2.46 Ang = 4.6511 Bohr.
    ! ELECTRON-PHONON (the cited, enabled channel): the two strongly-coupled
    ! optical modes responsible for the Kohn anomalies and dominant hot-carrier
    ! energy relaxation:
    !   * E2g at Gamma (intra-valley LO/TO, q~0): hw = 196 meV (~1580 cm^-1),
    !     EPC <g^2_Gamma> = 0.0405 eV^2.
    !   * A1' at K   (INTER-valley K<->K', the Kohn-anomaly mode): hw = 160 meV,
    !     EPC <g^2_K> = 0.0994 eV^2 (DFT); enhanced ~x2 by GW.
    ! [Piscanec, Lazzeri, Mauri, Ferrari, Robertson, PRL 93, 185503 (2004);
    !  GW enhancement: Lazzeri, Attaccalite, Wirtz, Mauri, PRB 78, 081406 (2008)].
    ! Both modes RELAX carriers DOWN THE CONE / across valleys, i.e. to a DIFFERENT
    ! k -- so graphene e-ph is physical only through the INTER-K ring path
    ! (yn_sbe_superres='y'); the k-local same-k search is inert on the 2-band cone.
    ! nu_sat = optical-phonon emission rate cap; OP emission is the dominant hot-
    ! carrier cooling channel, time scale ~10-100 fs [hot-carrier cooling lit.].
    ! Conservative cap ~5e13 s^-1 (~20 fs); tunable via sbe_eph_nu_sat.
    real(8), parameter, public :: GRAPH_A_BOHR    = 4.6511d0  ! a = 2.46 Ang
    real(8), parameter, public :: GRAPH_HW_E2G_MEV = 196.0d0  ! Gamma E2g optical [Piscanec 2004]
    real(8), parameter, public :: GRAPH_HW_A1P_MEV = 160.0d0  ! K A1' optical (Kohn anomaly) [Piscanec 2004]
    real(8), parameter, public :: GRAPH_G2_E2G    = 0.0405d0  ! <g^2_Gamma> [eV^2] [Piscanec 2004]
    real(8), parameter, public :: GRAPH_G2_A1P    = 0.0994d0  ! <g^2_K> [eV^2] DFT [Piscanec 2004]
    real(8), parameter, public :: GRAPH_GW_K      = 2.0d0     ! GW enhancement of the K coupling [Lazzeri 2008]
    real(8), parameter, public :: GRAPH_NU_SAT_SI = 5.0d13    ! OP-emission rate cap [~20 fs; hot-carrier cooling]

    ! =====================================================================
    ! Material registry -- the SINGLE place that maps a material name to all
    ! the per-material constants the SBE dissipation channels need (dielectric,
    ! impact-ionization fit, electron-phonon table). Adding a material is one
    ! `case` block in get_material_params() plus its name in MAT_SUPPORTED;
    ! every channel then auto-selects through the same struct. All numbers are
    ! the cited constants declared above -- the registry only assembles them.
    ! =====================================================================
    integer, parameter, public :: MAT_MAXPH = 8         ! capacity of a phonon table
    character(*), parameter, public :: MAT_SUPPORTED = 'GaAs, Si, Si_cb, CdS, graphene'

    type, public :: s_material_params
        logical       :: found        = .false.
        character(16) :: name         = ''
        ! crystal structure
        real(8)       :: a_lattice_au = 0d0       ! in-plane lattice constant a [Bohr]
        real(8)       :: cell_au(3)   = 0d0       ! expected &system al(1:3) box [Bohr]
        logical       :: is_diamond   = .false.   ! diamond (V^A=0) vs zincblende/wurtzite
        ! -- PROVENANCE GATES --------------------------------------------------
        ! A channel may be enabled for a material ONLY if its constants are
        ! backed by a cited source for THAT material. No source => the flag is
        ! .false. => the SBE init aborts (error stop) if the channel is enabled.
        ! Never fall back to another material's constants.
        logical       :: ii_ok        = .false.   ! impact-ionization fit cited?
        logical       :: eph_ok       = .false.   ! e-ph rate (nu_sat) cited?
        logical       :: eeh_ok       = .false.   ! carrier-carrier rate cited?
        logical       :: coulomb_ok   = .false.   ! dielectric for Coulomb cited?
        logical       :: auger_ok     = .false.   ! Auger coeff C cited?
        ! Auger recombination (Sec 13): R = C n^3, density-gated above n_gate.
        real(8)       :: auger_c_cm6s    = 0d0    ! Auger coefficient C [cm^6/s]
        real(8)       :: auger_n_gate_cm3 = 0d0   ! activation density [cm^-3]
        ! dielectric (Coulomb HF exchange / carrier screening)
        real(8)       :: eps0         = 1d0       ! static dielectric
        real(8)       :: eps_inf      = 1d0       ! high-frequency dielectric
        ! impact-ionization fit defaults
        character(20) :: ii_form      = 'stobbe_quartic'
        real(8)       :: ii_exponent  = 4d0       ! exponent a in P*(eps-Eth)^a
        real(8)       :: ii_prefactor = 2d12      ! P [s^-1 eV^-a]
        real(8)       :: ii_threshold_ev = 2.1d0  ! Eth above the CBM [eV]
        ! electron-phonon population-relaxing Lindblad table
        real(8)       :: eph_nu_sat_si = 1d14     ! saturation rate [s^-1]
        logical       :: eph_polar    = .false.   ! has a Frohlich polar-LO branch
        integer       :: eph_nph      = 0
        real(8)       :: eph_hw_mev(MAT_MAXPH) = 0d0   ! phonon energies [meV]
        real(8)       :: eph_wraw(MAT_MAXPH)   = 0d0   ! raw (un-normalized) D^2/hw weights
    end type s_material_params

contains

    ! Look up all per-material SBE constants by name. mp%found = .false. for an
    ! unknown material (callers that need a registry value must check it and stop
    ! with a helpful message). This is the ONLY function to extend when adding a
    ! material -- add a `case` and update MAT_SUPPORTED above.
    pure function get_material_params(name) result(mp)
        character(*), intent(in) :: name
        type(s_material_params)  :: mp
        integer :: p
        mp%name = name
        select case (trim(name))
        case ('GaAs')
            mp%found = .true.
            mp%a_lattice_au = 10.68d0
            mp%cell_au = (/ 10.68d0, 10.68d0, 10.68d0 /)   ! cubic
            mp%is_diamond   = .false.
            mp%ii_ok = .true.; mp%eph_ok = .true.; mp%eeh_ok = .true.; mp%coulomb_ok = .true.
            mp%eps0 = GAAS_EPS0;  mp%eps_inf = GAAS_EPS_INF
            mp%ii_form = 'stobbe_quartic'; mp%ii_exponent = 4d0
            mp%ii_prefactor = 2d12;        mp%ii_threshold_ev = 2.1d0
            mp%eph_nu_sat_si = 1.0d14      ! [Fischetti IEEE TED 38, 634 (1991)]
            mp%eph_polar = .true.
            ! Frohlich polar-LO (mode 1) + 5 intervalley modes
            mp%eph_nph = GAAS_N_IV + 1
            mp%eph_hw_mev(1) = GAAS_HW_LO_MEV
            do p = 1, GAAS_N_IV
                mp%eph_hw_mev(p + 1) = GAAS_IV_E_MEV(p)
                mp%eph_wraw(p + 1)   = GAAS_IV_D_EVANG(p)**2 / GAAS_IV_E_MEV(p)
            end do
            mp%eph_wraw(1) = sum(mp%eph_wraw(2:GAAS_N_IV + 1))   ! polar LO ~ dominant
        case ('Si', 'Si_cb')
            mp%found = .true.
            mp%a_lattice_au = 10.26d0
            mp%cell_au = (/ 10.26d0, 10.26d0, 10.26d0 /)   ! cubic
            mp%is_diamond   = .true.
            mp%ii_ok = .true.; mp%eph_ok = .true.; mp%eeh_ok = .true.; mp%coulomb_ok = .true.
            mp%eps0 = SI_EPS;  mp%eps_inf = SI_EPS     ! non-polar: eps_inf = eps0
            mp%ii_form = 'keldysh_quadratic'; mp%ii_exponent = 2d0
            mp%ii_prefactor = 2d12;           mp%ii_threshold_ev = 1.1d0
            mp%eph_nu_sat_si = 1.3d14      ! [Meng PRB 91, 075201 (2015)]
            mp%eph_polar = .false.
            mp%eph_nph = SI_N_PHONON       ! 6 intervalley g/f modes
            do p = 1, SI_N_PHONON
                mp%eph_hw_mev(p) = SI_PHONON_E_MEV(p)
                mp%eph_wraw(p)   = SI_PHONON_D_1E8EVCM(p)**2 / SI_PHONON_E_MEV(p)
            end do
        case ('CdS')
            ! Wurtzite (P6_3mc), polar, NON-centrosymmetric. Orthorhombic cell
            ! (a, a*sqrt3, c) via &system al(1:3) -- NOT cubic.
            ! [structure: Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967)]
            !
            ! CITED, ENABLED channels (the constants below are CdS-specific):
            !   * e-ph   : Frohlich polar-optical LO (PRIMARY room-T channel);
            !              hw_LO=38 meV [Raman], nu_sat=alpha*omega_LO=2.89e13
            !              from alpha=0.5 [cyclotron].
            !   * Coulomb: static dielectric eps0=8.9 [Berlincourt PR 129,1009 (1963)].
            !   * impact : Keldysh soft threshold, E_th=1.5*E_g=3.6 eV; the
            !              PREFACTOR is an uncited fit parameter -> sentinel (-1),
            !              so the run aborts unless the user sets sbe_ii_prefactor.
            !
            ! FORBIDDEN channel (no cited CdS *rate* -> eeh_ok = .false.):
            !   * carrier-carrier (yn_sbe_eeh). The implemented channel is an
            !     FD-thermalization parametrised by a single rate scale nu_cc.
            !     CdS has NO cited carrier-carrier rate, so enabling it with the
            !     default would borrow the generic 1e14 s^-1 scale that is cited
            !     ONLY for GaAs/Si (Goodnick-Lugli PRB 37, 2578; Fischetti-Laux
            !     PRB 38, 9721) -- forbidden by the strict provenance rule. The
            !     CdS citations (sub-100fs thermalization @ n>=1e18 [Shah 1986;
            !     Elsaesser 1991]) describe a *timescale*, not a rate. A user who
            !     supplies their own rate can still opt in via sbe_eeh_nu_sat (the
            !     same explicit-input escape hatch as the II prefactor).
            !
            ! FORBIDDEN channel (NO verified CdS Auger coefficient -> auger_ok=.false.):
            !   * Auger recombination (yn_sbe_auger). The former default
            !     "C = 2.0e-30 cm^6/s [Haury, PRB 57, 11513 (1998)]" was a
            !     FABRICATED reference (the real Haury 1997 PRL is on CdMnTe
            !     ferromagnetism, not CdS Auger) and was removed. CdS Auger is
            !     gated off; a user with a verified C may still opt in via
            !     sbe_auger_c_cm6s.
            mp%found = .true.
            mp%a_lattice_au = CDS_A_BOHR
            mp%cell_au = (/ CDS_A_BOHR, CDS_ASQ3_BOHR, CDS_C_BOHR /)  ! orthorhombic
            mp%is_diamond   = .false.                  ! V^A != 0 (broken inversion)
            mp%coulomb_ok = .true.; mp%eph_ok = .true.; mp%ii_ok = .true.
            mp%eeh_ok = .false.   ! no cited CdS carrier-carrier rate (see above)
            ! Auger recombination FORBIDDEN for CdS: no verified coefficient (the
            ! former "Haury PRB 57, 11513" C was fabricated, removed). A density
            ! gate is kept (Shah 1986) only for a user who supplies their own C
            ! via sbe_auger_c_cm6s; with auger_ok=.false. the default C is 0, so a
            ! plain yn_sbe_auger='y' run aborts (provenance) unless C is given.
            mp%auger_ok = .false.
            mp%auger_c_cm6s = 0d0;  mp%auger_n_gate_cm3 = CDS_EE_ACT_N
            mp%eps0 = CDS_EPS0;  mp%eps_inf = CDS_EPS_INF
            ! Frohlich polar-optical: a single dominant LO mode at 38 meV; the
            ! rate scale nu_sat = alpha*omega_LO is the cited Frohlich coupling.
            mp%eph_polar = .true.
            mp%eph_nph = 1
            mp%eph_hw_mev(1) = CDS_HW_LO_MEV
            mp%eph_wraw(1)   = 1d0                     ! single mode -> weight 1
            mp%eph_nu_sat_si = CDS_NU_SAT_SI
            ! Impact ionization: Keldysh soft threshold; E_th cited, prefactor is
            ! a fit parameter (left to the user via sbe_ii_prefactor).
            mp%ii_form = 'keldysh_quadratic'; mp%ii_exponent = 2d0
            mp%ii_threshold_ev = CDS_II_ETH_EV
            mp%ii_prefactor = -1d0     ! sentinel: NO cited CdS prefactor; the
            ! init requires the user to set sbe_ii_prefactor explicitly (fit param).
        case ('graphene')
            ! GAPLESS Dirac semimetal (honeycomb D6h, V^A=0, centrosymmetric).
            ! a minimal pi-model (Ramanujam 3-form-factor EPM); the Python ref
            ! emits the GS, so this registry entry only supplies the DISSIPATION
            ! constants. Hexagonal 2-atom primitive (or the rectangular 2-fold
            ! folded cell): cell_au below is informational only.
            !
            ! CITED, ENABLED channel:
            !   * e-ph: the two Kohn-anomaly optical modes E2g(Gamma,196 meV) and
            !     A1'(K,160 meV), EPC <g^2_Gamma>=0.0405, <g^2_K>=0.0994 eV^2 (x2 GW)
            !     [Piscanec PRL 93,185503 (2004); Lazzeri PRB 78,081406 (2008)].
            !     ** Both modes are NON-LOCAL on the cone (carriers relax to a
            !        DIFFERENT k), so graphene e-ph is valid ONLY through the
            !        inter-k RING (yn_sbe_superres='y'); the k-local same-k search
            !        is inert on the 2-band Dirac cone. **
            !
            ! FORBIDDEN channels (GAPLESS physics -- the gap-based maps don't apply):
            !   * impact ionization (ii_ok=.false.): the Stobbe/Keldysh fit is a
            !     (eps_kin - E_th)^a law with E_th a GAP threshold -- meaningless for
            !     a gapless cone. Graphene carrier multiplication is the (nearly
            !     thresholdless) 2D collinear process -- the Rana branch of the
            !     nonlocal-Auger task (wiki/07), NOT this gap-threshold channel.
            !   * Auger recombination (auger_ok=.false.): recombination across a gap
            !     with C*n^3 -- no gap, no single cm^6/s C. Again the 2D Rana branch.
            !   * carrier-carrier (eeh_ok=.false.): no cited graphene FD-thermalization
            !     RATE (the graphene e-e coupling alpha=2.2/eps_eff is a coupling, not
            !     a rate); the gapless e-e/CM physics is the Rana branch.
            !   * Coulomb HF (coulomb_ok=.false.): needs the 2D kernel V(q)=2*pi*e^2/
            !     (eps*q), not the 3D 4*pi/(eps*q^2) used by the current Sigma^HF.
            !
            ! Also: a Kuhn-Zurek (single-particle wave-packet) dephasing is
            ! UNPHYSICAL for gapless Dirac carriers (coherence loss is many-body) --
            ! the SBE init aborts on graphene + sbe_decoh_* (see bloch_solver init).
            mp%found = .true.
            mp%a_lattice_au = GRAPH_A_BOHR
            mp%cell_au = (/ GRAPH_A_BOHR, GRAPH_A_BOHR*sqrt(3d0), 0d0 /)  ! informational
            mp%is_diamond = .true.    ! V^A = 0 (centrosymmetric)
            mp%eph_ok = .true.
            mp%ii_ok = .false.; mp%auger_ok = .false.
            mp%eeh_ok = .false.; mp%coulomb_ok = .false.
            ! non-polar optical deformation-potential modes (no Frohlich LO branch)
            mp%eph_polar = .false.
            mp%eph_nph = 2
            mp%eph_hw_mev(1) = GRAPH_HW_E2G_MEV
            mp%eph_hw_mev(2) = GRAPH_HW_A1P_MEV
            ! raw per-mode weight = <g^2>/hw (analogue of D^2/hw); the A1' K-mode
            ! carries the x2 GW enhancement and dominates (Kohn anomaly).
            mp%eph_wraw(1) = GRAPH_G2_E2G            / GRAPH_HW_E2G_MEV
            mp%eph_wraw(2) = GRAPH_G2_A1P * GRAPH_GW_K / GRAPH_HW_A1P_MEV
            mp%eph_nu_sat_si = GRAPH_NU_SAT_SI
        case default
            mp%found = .false.
        end select
    end function get_material_params



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

    ! =====================================================================
    ! CDRB model dielectric function for the inter-k Coulomb channels
    ! [Cappellini, Del Sole, Reining, Bechstedt, PRB 47, 9892 (1993); used for
    !  Auger by Kioupakis et al. (K15), PRB 92, 035207 (2015), Eq. (8)]:
    !   eps(q) = 1 + [ (eps_inf-1)^-1 + alpha (q/q_TF)^2 + q^4/(4 wp^2) ]^-1
    ! in a.u. (hbar = m = e = 1), alpha = 1.563. Limits: eps(0) = eps_inf,
    ! eps(q->inf) -> 1 -- the large-q transparency that makes umklapp
    ! transitions strong (static eps0 over-screens them by ~10x [L90]).
    ! q_TF and wp are those of the VALENCE electron gas (n = nelec/volume).
    ! =====================================================================
    pure function eps_cdrb(q2, eps_inf, qtf2, wp2) result(eps)
        implicit none
        real(8), intent(in) :: q2, eps_inf, qtf2, wp2
        real(8) :: eps, denom
        real(8), parameter :: CDRB_ALPHA = 1.563d0
        if (eps_inf <= 1d0 + 1d-12 .or. qtf2 <= 0d0 .or. wp2 <= 0d0) then
            eps = 1d0
            return
        end if
        denom = 1d0 / (eps_inf - 1d0) + CDRB_ALPHA * q2 / qtf2 + q2 * q2 / (4d0 * wp2)
        eps = 1d0 + 1d0 / denom
    end function eps_cdrb

    ! =====================================================================
    ! Screened inter-k Coulomb weight WITH the umklapp G-sum, in the CARTESIAN
    ! metric of the actual (possibly non-orthogonal) cell:
    !   vq = sum_G 1 / [ eps(|q+G|) (|q+G|^2 + lambda^2 + q2reg) ]
    ! over the 27 neighbouring reciprocal images G = n1 b1 + n2 b2 + n3 b3,
    ! n_i in {-1,0,1} (bmat rows = b1..b3 [a.u.]; df = fractional coordinate
    ! difference of the two grid points, in (-1,1) per component so the 27
    ! images bracket the minimum image). This implements the two [L90]
    ! "must-not-drop" pieces at the overlap-free level (I(G) -> 1): the
    ! umklapp G-sum and the q-dependent eps(q); the Bloch-overlap factors
    ! I_{13}(G) I_{24}(G'-G) need the plane-wave coefficients the SBE does not
    ! carry and remain a refinement. lambda2 = free-carrier screening (Si: 0
    ! by Burt's dynamical argument [L90]); q2reg = grid-scale q->0 regulariser
    ! of the discrete BZ sum (numerical, refines with the k-grid -- NOT a
    ! physics constant, unlike the old fixed reduced-space kappa2 = 0.05).
    ! =====================================================================
    pure function interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg) result(vq)
        implicit none
        real(8), intent(in) :: df(3), bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8) :: vq, qc(3), q2
        integer :: g1, g2, g3
        vq = 0d0
        do g3 = -1, 1
            do g2 = -1, 1
                do g1 = -1, 1
                    qc(1:3) = (df(1) + g1) * bmat(1, 1:3) &
                            + (df(2) + g2) * bmat(2, 1:3) &
                            + (df(3) + g3) * bmat(3, 1:3)
                    q2 = dot_product(qc, qc)
                    vq = vq + 1d0 / (eps_cdrb(q2, eps_inf, qtf2, wp2) * (q2 + lambda2 + q2reg))
                end do
            end do
        end do
    end function interk_vq

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

    ! Part C3: energy-conserving final-state weights (broadened-delta energy
    ! bins; NO tetrahedra, NO Monte-Carlo). For each candidate energy e_cand(i),
    !   w(i) = gaussian_bin(e_cand(i) - e_target, sigma)   if |dE| <= max_window
    !        = 0                                           otherwise,
    ! and wsum = sum(w). max_window realizes the energy-windowed expanding-radius
    ! search (the caller grows it with the hot-carrier excess energy). The caller
    ! normalizes the chosen transitions by wsum. Deterministic (CPTP-safe).
    ! [Stobbe PRB 49, 4494 (1994); Kunikiyo JAP 75, 297 (1994)]
    pure subroutine energy_partner_weights(e_target, e_cand, n_cand, sigma, &
                                           max_window, w, wsum)
        integer, intent(in)  :: n_cand
        real(8), intent(in)  :: e_target, e_cand(n_cand), sigma, max_window
        real(8), intent(out) :: w(n_cand), wsum
        integer :: i
        real(8) :: dE
        wsum = 0d0
        do i = 1, n_cand
            dE = e_cand(i) - e_target
            if (abs(dE) <= max_window) then
                w(i) = gaussian_bin(dE, sigma)
            else
                w(i) = 0d0
            end if
            wsum = wsum + w(i)
        end do
    end subroutine energy_partner_weights

    ! =====================================================================
    ! Part F -- carrier-carrier (e-e/e-h) thermalization map and the
    ! number-and-energy-conserving Fermi-Dirac fit it relaxes toward.
    ! =====================================================================

    pure function fermi_dirac(beta, e, mu) result(f)
        real(8), intent(in) :: beta, e, mu
        real(8) :: f, x
        x = beta * (e - mu)
        if (x > 7d2) then
            f = 0d0
        else if (x < -7d2) then
            f = 1d0
        else
            f = 1d0 / (exp(x) + 1d0)
        end if
    end function fermi_dirac

    ! For fixed beta, bisect mu so sum_a f_FD(eps_a) = ntot; return E and ftgt.
    subroutine e_at_beta(nlev, eps, ntot, beta, Eout, mu, ftgt)
        integer, intent(in)  :: nlev
        real(8), intent(in)  :: eps(nlev), ntot, beta
        real(8), intent(out) :: Eout, mu, ftgt(nlev)
        real(8) :: mlo, mhi, mmid, nsum
        integer :: it, a
        mlo = minval(eps) - 1d0; mhi = maxval(eps) + 1d0
        do it = 1, 100
            mmid = 0.5d0 * (mlo + mhi)
            nsum = 0d0
            do a = 1, nlev
                nsum = nsum + fermi_dirac(beta, eps(a), mmid)
            end do
            if (nsum < ntot) then
                mlo = mmid
            else
                mhi = mmid
            end if
        end do
        mu = 0.5d0 * (mlo + mhi)
        Eout = 0d0
        do a = 1, nlev
            ftgt(a) = fermi_dirac(beta, eps(a), mu)
            Eout = Eout + eps(a) * ftgt(a)
        end do
    end subroutine e_at_beta

    ! Fit (beta, mu) so sum_a f_FD = ntot and sum_a eps_a f_FD = etot, by nested
    ! bisection (outer beta log-scale, inner mu). ok=.false. if etot is outside
    ! [E(T=0), E(T=inf)] (not Fermi-Dirac representable, e.g. population inversion).
    subroutine fit_fermi_dirac(nlev, eps, ntot, etot, beta, mu, ftgt, ok)
        integer, intent(in)  :: nlev
        real(8), intent(in)  :: eps(nlev), ntot, etot
        real(8), intent(out) :: beta, mu, ftgt(nlev)
        logical, intent(out) :: ok
        real(8) :: blo, bhi, bmid, Elo, Ehi, Emid, mtmp
        integer :: it
        ok = .false.
        blo = 1d-4; bhi = 1d4
        call e_at_beta(nlev, eps, ntot, blo, Elo, mtmp, ftgt)
        call e_at_beta(nlev, eps, ntot, bhi, Ehi, mtmp, ftgt)
        if (etot > Elo + 1d-12 .or. etot < Ehi - 1d-12) return
        mu = mtmp; bmid = sqrt(blo * bhi)
        do it = 1, 100
            bmid = sqrt(blo * bhi)
            call e_at_beta(nlev, eps, ntot, bmid, Emid, mu, ftgt)
            if (Emid > etot) then
                blo = bmid
            else
                bhi = bmid
            end if
            if (abs(Emid - etot) < 1d-13 * max(1d0, abs(etot))) exit
        end do
        beta = bmid
        ok = .true.
    end subroutine fit_fermi_dirac

    ! Carrier-carrier CPTP relaxation map (Part F): relax the adiabatic
    ! populations toward the number- and energy-matched Fermi-Dirac and damp the
    ! coherences (EID):  rho -> (1-a) rho + a diag(occ f_FD), a = 1-exp(-nu tau),
    ! a convex combination of identity and a constant-state channel -> CPTP.
    ! Conserves Tr rho (number) and sum_a eps_a rho_aa (energy) exactly. A no-op
    ! when the population set is empty/full or not Fermi-Dirac representable.
    subroutine carrier_carrier_relax(nlev, rho, eps, occ, nu, tau)
        integer,    intent(in)    :: nlev
        complex(8), intent(inout) :: rho(nlev, nlev)
        real(8),    intent(in)    :: eps(nlev), occ, nu, tau
        real(8) :: f(nlev), ftgt(nlev), ntot, etot, alpha, beta, mu
        integer :: a, b
        logical :: ok
        ntot = 0d0; etot = 0d0
        do a = 1, nlev
            f(a) = min(max(real(rho(a, a)) / occ, 0d0), 1d0)
            ntot = ntot + f(a)
            etot = etot + eps(a) * f(a)
        end do
        if (ntot < 1d-9 .or. (dble(nlev) - ntot) < 1d-9) return
        call fit_fermi_dirac(nlev, eps, ntot, etot, beta, mu, ftgt, ok)
        if (.not. ok) return
        alpha = 1d0 - exp(-nu * tau)
        do b = 1, nlev
            do a = 1, nlev
                if (a == b) then
                    rho(a, a) = (1d0 - alpha) * rho(a, a) + dcmplx(alpha * occ * ftgt(a), 0d0)
                else
                    rho(a, b) = (1d0 - alpha) * rho(a, b)
                end if
            end do
        end do
    end subroutine carrier_carrier_relax

    ! =====================================================================
    ! INTER-K e-ph (the "super-mode ring" intervalley channel, Part C5/D).
    ! ---------------------------------------------------------------------
    ! The k-LOCAL apply_eph_relaxation relaxes a carrier to the energy-matched
    ! partner BAND at the SAME k -- correct only when the final valley folds onto
    ! a same-k band (folded cell) or for intra-valley polar-optical. On the
    ! PRIMITIVE cell the intervalley final states live at DIFFERENT k, so the
    ! search must run over all (k,band). This pure routine takes the GATHERED
    ! Houston spectrum and returns the net diagonal population change; the caller
    ! gathers eval/f (all-gather / ring) and applies dpop. Enabled when the ring
    ! (yn_sbe_superres) is on -- "if the ring is on, inter-k goes through it".
    !
    ! EXACTLY trace-conserving (CPTP-safe): each source (ik,a) loses
    !   out_tot = f*(1-exp(-Gamma_out*tau)) <= f   (no negativity),
    ! distributed to destinations proportional to their partial rate, so
    ! sum(dpop) = 0 by construction. Detailed balance (emission fe / absorption
    ! fa from the Bose factor) + Pauli blocking (1 - f_dest/occ_max) at the
    ! destination. Intra-k (jq==ik) is included automatically -- this SUBSUMES
    ! the intra-k channel when the ring is on. O(nk^2 nba^2 nph) all-pairs (same
    ! order as the Coulomb all-pairs sum it rides alongside).
    subroutine eph_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, evbm, &
                               nph, hw, wrel, nb_bose, nu_sat, nu_eps0, nu_n, &
                               sigma, tau, dpop, gout)
        implicit none
        integer, intent(in)  :: nk, nba, nph
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, a2half
        real(8), intent(in)  :: ecbm, evbm, hw(nph), wrel(nph), nb_bose(nph)
        real(8), intent(in)  :: nu_sat, nu_eps0, nu_n, sigma, tau   ! nu_n = saturation exponent
        real(8), intent(out) :: dpop(nba, nk)
        real(8), intent(out), optional :: gout(nba, nk)  ! total out-rate Gamma_out per source (coherence damping)
        integer :: ik, jq, a, b, ip
        real(8) :: eps_kin, nu_a, fe, fa, dE, shp, th, blk, gam, gamtot, out_tot
        real(8) :: gpart(nba, nk)
        real(8), parameter :: occ_eps = 1d-12

        dpop = 0d0
        if (present(gout)) gout = 0d0
        do ik = 1, nk
            do a = 1, nba
                if (f(a, ik) < occ_eps) cycle
                ! carrier kinetic energy from the nearest band edge (restore A^2/2,
                ! the k-independent Houston offset; it cancels in energy MATCHING).
                eps_kin = max(eval(a, ik) + a2half - ecbm, &
                              evbm - (eval(a, ik) + a2half), 0d0)
                nu_a = nu_saturation(eps_kin, nu_sat, nu_eps0, nu_n)
                if (nu_a * tau < 1d-14) cycle

                gpart = 0d0
                gamtot = 0d0
                do ip = 1, nph
                    call eph_thermal_split(nb_bose(ip), fe, fa)
                    do jq = 1, nk
                        do b = 1, nba
                            if (jq == ik .and. b == a) cycle
                            if (eval(b, jq) < eval(a, ik)) then        ! emission (down)
                                dE = abs((eval(a, ik) - eval(b, jq)) - hw(ip)); th = fe
                            else                                       ! absorption (up)
                                dE = abs((eval(b, jq) - eval(a, ik)) - hw(ip)); th = fa
                            end if
                            if (th <= 0d0) cycle
                            shp = gaussian_shape(dE, sigma)
                            if (shp <= 0d0) cycle
                            blk = min(max(1d0 - f(b, jq) / occ_max, 0d0), 1d0)
                            gam = nu_a * wrel(ip) * th * shp * blk
                            gpart(b, jq) = gpart(b, jq) + gam
                            gamtot = gamtot + gam
                        end do
                    end do
                end do

                if (present(gout)) gout(a, ik) = gamtot
                if (gamtot * tau < 1d-14) cycle
                out_tot = f(a, ik) * (1d0 - exp(-gamtot * tau))
                dpop(a, ik) = dpop(a, ik) - out_tot
                do jq = 1, nk
                    do b = 1, nba
                        if (gpart(b, jq) > 0d0) &
                            dpop(b, jq) = dpop(b, jq) + out_tot * gpart(b, jq) / gamtot
                    end do
                end do
            end do
        end do
    end subroutine eph_interk_dpop

    ! =====================================================================
    ! INTER-K (momentum-conserving) impact ionization through the ring.
    ! ---------------------------------------------------------------------
    ! The TRUE 2-particle event: a hot conduction e- (k1, band ih) and a valence
    ! e- (k2, iv) -> two conduction e- at (k1', ic) and (k2', ic) leaving a hole,
    ! with crystal momentum k1 + k2 = k1' + k2' (mod G) [k2' from mp_partner_triple]
    ! and energy E(k1,ih)+E(k2,iv) = E(k1',ic)+E(k2',ic) [broadened Fermi golden
    ! rule]. The threshold magnitude g0 = pref*(eps_kin - E_th)^expo is the cited
    ! Stobbe-fit rate; the screened Coulomb |V(q)|^2 (q = k1-k1') and the energy
    ! delta SHAPE which momentum-conserving final config it goes to. EXACTLY
    ! trace-conserving: each event writes -amt,-amt,+amt,+amt (primary out +
    ! valence out = the two conduction gains), so sum(dpop)=0 by construction.
    ! Primary out capped at f(ih,k1)*(1-exp(-Gamma*tau)) <= f (no negativity).
    ! klut(0:nk-1): flattened MP lookup, triple (m1,m2,m3) -> ik via
    ! lidx = m1 + n1*(m2 + n2*m3). The caller precomputes it once.
    subroutine ii_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                              pref, expo, iv, ic, kidx, kn, klut, &
                              bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dpop)
        implicit none
        integer, intent(in)  :: nk, nba, iv, ic, kidx(3, nk), kn(3), klut(0:nk-1)
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, a2half
        real(8), intent(in)  :: ecbm, eth, pref, expo, sigma, tau
        real(8), intent(in)  :: bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8), intent(out) :: dpop(nba, nk)
        integer :: i1, i1p, i2, ih, ipass, jj, m2p(3), d
        real(8) :: ekin, dd, g0, etgt, vq, df(3), shp, pauli, gpart, gamtot, out_tot, amt
        real(8), parameter :: occ_eps = 1d-12

        dpop = 0d0
        if (iv < 1 .or. ic > nba .or. ic <= iv) return

        do i1 = 1, nk
            do ih = ic, nba
                if (f(ih, i1) < occ_eps) cycle
                ekin = eval(ih, i1) + a2half - ecbm
                dd = ekin - eth
                if (dd <= 0d0) cycle
                g0 = pref * dd ** expo            ! cited Stobbe-fit total magnitude

                ! Pass 1 (ipass=1): accumulate gamtot; Pass 2: distribute out_tot.
                gamtot = 0d0
                out_tot = 0d0
                do ipass = 1, 2
                    if (ipass == 2) then
                        if (gamtot * tau < 1d-14) exit
                        out_tot = f(ih, i1) * (1d0 - exp(-gamtot * tau))
                    end if
                    do i1p = 1, nk
                        ! transferred momentum q = k1 - k1' (+G): CDRB-screened
                        ! Cartesian-metric weight with the 27-image umklapp sum
                        do d = 1, 3
                            df(d) = dble(kidx(d, i1) - kidx(d, i1p)) / dble(max(kn(d), 1))
                        end do
                        vq = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                        do i2 = 1, nk
                            call mp_partner_triple(kidx(:,i1), kidx(:,i2), kidx(:,i1p), kn, m2p)
                            jj = klut(m2p(1) + kn(1) * (m2p(2) + kn(2) * m2p(3))) ! O(1) k2'
                            if (jj < 1) cycle
                            etgt = eval(ih,i1) + eval(iv,i2) - eval(ic,i1p) - eval(ic,jj)
                            shp = gaussian_shape(etgt, sigma)
                            if (shp <= 0d0) cycle
                            pauli = (f(iv,i2) / occ_max) &
                                  * min(max(1d0 - f(ic,i1p)/occ_max, 0d0), 1d0) &
                                  * min(max(1d0 - f(ic,jj )/occ_max, 0d0), 1d0)
                            gpart = g0 * vq * shp * pauli
                            if (ipass == 1) then
                                gamtot = gamtot + gpart
                            else
                                amt = out_tot * gpart / gamtot
                                dpop(ih, i1)  = dpop(ih, i1)  - amt   ! hot primary leaves
                                dpop(ic, i1p) = dpop(ic, i1p) + amt   ! primary relaxed
                                dpop(iv, i2)  = dpop(iv, i2)  - amt   ! valence e- -> hole
                                dpop(ic, jj)  = dpop(ic, jj)  + amt   ! promoted to conduction
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine ii_interk_dpop

    ! =====================================================================
    ! INTER-K (momentum-conserving) AUGER RECOMBINATION through the ring --
    ! the EXACT TIME-REVERSE of ii_interk_dpop (detailed balance).
    ! ---------------------------------------------------------------------
    ! The reverse 2-particle event: two conduction e- at (k1',ic) and (k2',ic)
    ! + a hole at (k2,iv) -> one e- recombines into the hole and the released
    ! energy promotes the other to the hot state (k1,ih). SAME quadruples,
    ! SAME weights g0*|V(q)|^2*delta_sigma as the impact-ionization kernel (no
    ! new constant: the rate scale IS the cited Stobbe/Keldysh II magnitude --
    ! Auger and II share |M|^2, only the occupation factors swap [Rana 2007;
    ! Kioupakis 2015]), REVERSED occupation product
    !   pauli_rev = (1 - f_v/occ) * (f_c1'/occ) * (f_c2'/occ)
    ! and NEGATED dpop signs: +hot, -c1', +valence, -c2'. Trace-conserving by
    ! construction (sum(dpop)=0). The hot-state gain is capped at
    ! (occ - f(ih,k1))*(1-exp(-Gamma*tau)) (no over-filling), mirroring the
    ! II primary-out cap. DETAILED BALANCE: for Fermi-Dirac occupations and an
    ! energy-conserving quadruple, f1 f2 (1-f3)(1-f4) = (1-f1)(1-f2) f3 f4, so
    ! in the linear (Gamma*tau -> 0) regime the net II + Auger dpop vanishes
    ! identically -- the equilibrium-fixed-point unit test.
    ! =====================================================================
    subroutine auger_interk_dpop(nk, nba, eval, f, occ_max, a2half, ecbm, eth, &
                                 pref, expo, iv, ic, kidx, kn, klut, &
                                 bmat, eps_inf, qtf2, wp2, lambda2, q2reg, sigma, tau, dpop)
        implicit none
        integer, intent(in)  :: nk, nba, iv, ic, kidx(3, nk), kn(3), klut(0:nk-1)
        real(8), intent(in)  :: eval(nba, nk), f(nba, nk), occ_max, a2half
        real(8), intent(in)  :: ecbm, eth, pref, expo, sigma, tau
        real(8), intent(in)  :: bmat(3,3), eps_inf, qtf2, wp2, lambda2, q2reg
        real(8), intent(out) :: dpop(nba, nk)
        integer :: i1, i1p, i2, ih, ipass, jj, m2p(3), d
        real(8) :: ekin, dd, g0, etgt, vq, df(3), shp, pauli, gpart, gamtot, in_tot, amt, room
        real(8), parameter :: occ_eps = 1d-12

        dpop = 0d0
        if (iv < 1 .or. ic > nba .or. ic <= iv) return

        do i1 = 1, nk
            do ih = ic, nba
                room = occ_max - f(ih, i1)          ! empty hot-state phase space
                if (room < occ_eps) cycle
                ekin = eval(ih, i1) + a2half - ecbm
                dd = ekin - eth
                if (dd <= 0d0) cycle
                g0 = pref * dd ** expo              ! same II magnitude (shared |M|^2)

                gamtot = 0d0
                in_tot = 0d0
                do ipass = 1, 2
                    if (ipass == 2) then
                        if (gamtot * tau < 1d-14) exit
                        in_tot = room * (1d0 - exp(-gamtot * tau))
                    end if
                    do i1p = 1, nk
                        ! SAME CDRB-screened umklapp weight as the II kernel
                        ! (shared |M|^2 -> detailed balance preserved exactly)
                        do d = 1, 3
                            df(d) = dble(kidx(d, i1) - kidx(d, i1p)) / dble(max(kn(d), 1))
                        end do
                        vq = interk_vq(df, bmat, eps_inf, qtf2, wp2, lambda2, q2reg)
                        do i2 = 1, nk
                            call mp_partner_triple(kidx(:,i1), kidx(:,i2), kidx(:,i1p), kn, m2p)
                            jj = klut(m2p(1) + kn(1) * (m2p(2) + kn(2) * m2p(3)))
                            if (jj < 1) cycle
                            etgt = eval(ih,i1) + eval(iv,i2) - eval(ic,i1p) - eval(ic,jj)
                            shp = gaussian_shape(etgt, sigma)
                            if (shp <= 0d0) cycle
                            ! REVERSED occupations: hole present at (iv,i2), both
                            ! conduction sources occupied at (ic,i1p) and (ic,jj).
                            pauli = min(max(1d0 - f(iv,i2)/occ_max, 0d0), 1d0) &
                                  * (f(ic,i1p) / occ_max) &
                                  * (f(ic,jj ) / occ_max)
                            gpart = g0 * vq * shp * pauli
                            if (ipass == 1) then
                                gamtot = gamtot + gpart
                            else
                                amt = in_tot * gpart / gamtot
                                dpop(ih, i1)  = dpop(ih, i1)  + amt   ! promoted to hot
                                dpop(ic, i1p) = dpop(ic, i1p) - amt   ! promoted-source leaves
                                dpop(iv, i2)  = dpop(iv, i2)  + amt   ! hole filled (recombination)
                                dpop(ic, jj)  = dpop(ic, jj)  - amt   ! recombining e- leaves
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine auger_interk_dpop

    ! =====================================================================
    ! Monkhorst-Pack momentum-conservation index map (for the inter-k / nonlocal
    ! collision channels: impact ionization, e-e). A regular MP mesh has points
    ! k_m = (2m - n + 1)/(2n) per dimension (m = 0..n-1), so crystal-momentum
    ! conservation k1 + k2 - k1' (mod G) is EXACT integer arithmetic on the index
    ! triples: m2' = mod(m1 + m2 - m1', n). mp_grid_triple recovers the triple
    ! from the reduced coordinate and reports the residual (0 for a true MP point,
    ! so the caller can gate the channel off on a non-MP / symmetry-reduced grid).
    ! =====================================================================
    pure subroutine mp_grid_triple(kred, n, m, resid)
        real(8), intent(in)  :: kred(3)
        integer, intent(in)  :: n(3)
        integer, intent(out) :: m(3)
        real(8), intent(out) :: resid     ! max wrapped |kred - reconstructed MP point|
        integer :: d
        real(8) :: km
        resid = 0d0
        do d = 1, 3
            if (n(d) <= 0) then
                m(d) = 0
                cycle
            end if
            m(d) = modulo(nint(kred(d) * n(d) + (n(d) - 1) * 0.5d0), n(d))
            km = (2d0 * m(d) - n(d) + 1d0) / (2d0 * n(d))
            resid = max(resid, abs(modulo(kred(d) - km + 0.5d0, 1d0) - 0.5d0))
        end do
    end subroutine mp_grid_triple

    ! Momentum-conserving partner index triple: m2' = mod(m1 + m2 - m1', n).
    pure subroutine mp_partner_triple(m1, m2, m1p, n, m2p)
        integer, intent(in)  :: m1(3), m2(3), m1p(3), n(3)
        integer, intent(out) :: m2p(3)
        integer :: d
        do d = 1, 3
            m2p(d) = modulo(m1(d) + m2(d) - m1p(d), max(n(d), 1))
        end do
    end subroutine mp_partner_triple

    ! =====================================================================
    ! Velocity-gauge basis sufficiency / N_b convergence primitives.
    ! These quantify whether the band count N_b carried into the dynamics is
    ! large enough -- a correctness axis SEPARATE from the plane-wave cutoff,
    ! and NOT fixed by rotating into the Houston basis (which only diagonalizes
    ! the already-truncated H_VG^(N_b) = P_Nb H_VG P_Nb). See the wiki page
    ! "VG Basis Sufficiency & N_b Convergence".
    ! [Hylleraas-Undheim Z.Phys.65,759(1930); MacDonald PR 43,830(1933);
    !  Wismer-Yakovlev PRB 97,144302(2018)]
    ! =====================================================================

    ! Criterion (c): dimensionless adiabatic admixture of band c into band a
    !   eta_ac = A_max |pi_ac| / |eps_a - eps_c|.
    ! The basis is safe for level a when eta to the FIRST DISCARDED band << 1;
    ! eta >~ 1 means the admixture is nonperturbative and c cannot be dropped.
    pure function vg_eta_admixture(A_max, pi_ac, gap_ac) result(eta)
        real(8), intent(in) :: A_max, pi_ac, gap_ac
        real(8) :: eta
        if (abs(gap_ac) <= 0d0) then
            eta = huge(1d0)
        else
            eta = A_max * abs(pi_ac) / abs(gap_ac)
        end if
    end function vg_eta_admixture

    ! Second-order (Rayleigh-Schrodinger) truncation shift of a retained level a
    ! from the n DISCARDED bands c > N_b:
    !   delta eps_a = sum_c A^2 |pi_ac|^2 / (eps_a - eps_c).
    ! This is exactly the error the Houston basis inherits (it diagonalizes the
    ! projected operator, so this shift is NOT recovered). Valid for eta << 1.
    pure function vg_trunc_shift2(A, pi_ac, eps_a, eps_disc, n) result(dshift)
        integer, intent(in) :: n
        real(8), intent(in) :: A, pi_ac(n), eps_a, eps_disc(n)
        real(8) :: dshift
        integer :: c
        real(8) :: denom
        dshift = 0d0
        do c = 1, n
            denom = eps_a - eps_disc(c)
            if (abs(denom) > 0d0) dshift = dshift + A*A*pi_ac(c)*pi_ac(c) / denom
        end do
    end function vg_trunc_shift2

    ! Criterion (b): relative L2 convergence error between an observable computed
    ! with N_b and with N_b+Delta bands (current, HHG spectrum, carrier number):
    !   eps_conv = ||O_large - O_small|| / ||O_large||.
    pure function vg_conv_error(o_small, o_large, n) result(eps_conv)
        integer, intent(in) :: n
        real(8), intent(in) :: o_small(n), o_large(n)
        real(8) :: eps_conv, num, den
        integer :: i
        num = 0d0; den = 0d0
        do i = 1, n
            num = num + (o_large(i) - o_small(i))**2
            den = den + o_large(i)**2
        end do
        if (den <= 0d0) then
            eps_conv = 0d0
        else
            eps_conv = sqrt(num / den)
        end if
    end function vg_conv_error

    ! Criterion (a): does the top retained adiabatic band's peak occupation
    ! P_top = max_k rho~_{Nb,Nb}(k) exceed the tolerance (default 1e-3)? A .true.
    ! means population reached the basis edge -> enlarge N_b.
    pure function vg_ptop_exceeds(ptop, thr) result(over)
        real(8), intent(in) :: ptop, thr
        logical :: over
        over = (ptop > thr)
    end function vg_ptop_exceeds

end module sbe_superres_ssbe
