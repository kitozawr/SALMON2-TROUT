!
!  Copyright 2017-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module salmon_global
  implicit none

!Parameters for pseudo-potential
  integer, parameter :: maxmki=10
   !shinohara
  integer :: ipsfileform(maxmki)   ! file format for pseudo potential
  character(16)  :: ps_format(maxmki)
! List of pseudopotential file formats
  integer,parameter :: n_Yabana_Bertsch_psformat = 1 !.rps
  integer,parameter :: n_ABINIT_psformat = 2 ! .pspnc
  integer,parameter :: n_FHI_psformat = 3 ! .cpi
  integer,parameter :: n_ABINITFHI_psformat = 4 ! .fhi

! Flag for atomic coordinate type
  integer :: iflag_atom_coor
  integer,parameter :: ntype_atom_coor_none      = 0
  integer,parameter :: ntype_atom_coor_cartesian = 1
  integer,parameter :: ntype_atom_coor_reduced   = 2

! Flag for suppress standard outputs 
  logical :: quiet

! For band
  integer,parameter :: max_num_of_segments = 10

  character(16)  :: calc_mode      !old input variable, but used as a flag; move later

!Input variables
!! &calculation
  character(32)  :: theory
  character(1)   :: yn_md
  character(1)   :: yn_opt
  character(1)   :: yn_dc
  character(1)   :: yn_conventional_from_dcdft
!! &control
  character(256) :: sysname
  character(256) :: base_directory
  character(1)   :: yn_restart
  character(256) :: directory_read_data
  character(1)   :: yn_self_checkpoint
  integer        :: checkpoint_interval
  character(1)   :: yn_reset_step_restart
  character(256) :: read_gs_restart_data
  character(256) :: write_gs_restart_data
  real(8)        :: time_shutdown
  character(20)  :: method_wf_distributor
  integer        :: nblock_wf_distribute
  !remove later
  character(1)   :: write_gs_wfn_k
  character(1)   :: write_rt_wfn_k

!! &units
  character(16)  :: unit_system
  character(16)  :: unit_time
  character(16)  :: unit_length
  character(16)  :: unit_energy
  character(16)  :: unit_charge

!! &parallel
  integer        :: nproc_k
  integer        :: nproc_ob
  integer        :: nproc_rgrid(3)
  character(1)   :: yn_ffte
  character(1)   :: yn_fftw
  character(1)   :: yn_scalapack
  character(1)   :: yn_gramschmidt_blas
  character(1)   :: yn_eigenexa
  character(1)   :: yn_diagonalization_red_mem
  character(32)  :: process_allocation

!! &system
  integer        :: iperiodic  !this is old keyword but still defined here
  character(1)   :: yn_periodic
  character(16)  :: spin
  real(8)        :: al(3)
  real(8)        :: al_vec1(3),al_vec2(3),al_vec3(3)
  integer        :: nstate
  integer        :: nelec
  integer        :: nelec_spin(2)
  real(8)        :: temperature
  real(8)        :: temperature_k
  integer        :: nelem
  integer        :: natom
  character(256) :: file_atom_coor
  character(256) :: file_atom_red_coor
  character(1)   :: yn_spinorbit
  character(3)   :: yn_symmetry
  character(16)  :: absorbing_boundary
  real(8)        :: imagnary_potential_w0
  real(8)        :: imagnary_potential_dr

!! &pseudo
  character(256) :: file_pseudo(maxmki)
  integer        :: lmax_ps(maxmki)
  integer        :: lloc_ps(maxmki)
  integer        :: izatom(maxmki)
  character(1)   :: yn_psmask
  real(8)        :: alpha_mask
  real(8)        :: gamma_mask
  real(8)        :: eta_mask

!! &functional
  character(64)  :: xc !, xcname
  character(64)  :: xname
  character(64)  :: cname
  character(64)  :: alibx
  character(64)  :: alibc
  character(64)  :: alibxc
  real(8)        :: cval

!! &rgrid
  real(8)        :: dl(3)
  integer        :: num_rgrid(3)

!! &kgrid
  integer        :: num_kgrid(3)
  character(256) :: file_kw
  real(8)        :: dk_shift(3)

!! &tgrid
  integer        :: nt
  real(8)        :: dt
  integer        :: gram_schmidt_interval

!! &propagation
  integer        :: n_hamil
  character(16)  :: propagator
  character(1)   :: yn_fix_func
  character(1)   :: yn_predictor_corrector

!! &scf
  character(8)   :: method_init_wf
  integer        :: iseed_number_change
  character(8)   :: method_min
  integer        :: ncg,ncg_init
  character(16)  :: method_mixing
  real(8)        :: mixrate
  integer        :: nmemory_mb
  real(8)        :: alpha_mb
  integer        :: nmemory_p
  real(8)        :: beta_p
  character(1)   :: yn_auto_mixing
  real(8)        :: update_mixing_ratio
  integer        :: nscf
  character(1)   :: yn_subspace_diagonalization
  character(16)  :: convergence
  real(8)        :: threshold
  integer        :: nscf_init_redistribution
  integer        :: nscf_init_no_diagonal
  integer        :: nscf_init_mix_zero
  real(8)        :: conv_gap_mix_zero
  character(16)  :: method_init_density
  real(8)        :: magdir_atom(100)
  character(1)   :: yn_preconditioning
  real(8)        :: alpha_pre

!! &emfield
  character(2)   :: trans_longi
  character(16)  :: ae_shape1
  character(256) :: file_input1
  real(8)        :: e_impulse
  real(8)        :: E_amplitude1
  real(8)        :: I_wcm2_1
  real(8)        :: tw1
  real(8)        :: omega1
  real(8)        :: epdir_re1(3)
  real(8)        :: epdir_im1(3)
  real(8)        :: phi_cep1
  character(16)  :: ae_shape2
  real(8)        :: E_amplitude2
  real(8)        :: I_wcm2_2
  real(8)        :: tw2
  real(8)        :: omega2
  real(8)        :: epdir_re2(3)
  real(8)        :: epdir_im2(3)
  real(8)        :: phi_cep2
  real(8)        :: t1_t2
  real(8)        :: t1_start
  integer        :: num_dipole_source
  real(8)        :: vec_dipole_source(3,2)
  real(8)        :: cood_dipole_source(3,2)
  real(8)        :: rad_dipole_source

!! &singlescale
  character(32)  :: method_singlescale
  real(8)        :: cutoff_G2_emfield
  character(1)   :: yn_symmetrized_stencil
  character(1)   :: yn_put_wall_z_boundary
  real(8)        :: wall_height
  real(8)        :: wall_width

!! &multiscale
  character(16)  :: fdtddim
  character(16)  :: twod_shape
  integer        :: nx_m
  integer        :: ny_m
  integer        :: nz_m
  real(8)        :: hx_m
  real(8)        :: hy_m
  real(8)        :: hz_m
  integer        :: nksplit ! this parameter will be deprecated in a future release
  integer        :: nxysplit ! this parameter will be deprecated in a future release
  integer        :: nxvacl_m ! this parameter will be deprecated in a future release
  integer        :: nxvacr_m ! this parameter will be deprecated in a future release
  integer        :: nxvac_m(1:2)
  integer        :: nyvac_m(1:2)
  integer        :: nzvac_m(1:2)
  integer        :: nx_origin_m
  integer        :: ny_origin_m
  integer        :: nz_origin_m
  integer        :: out_ms_region_ix_m(1:2)
  integer        :: out_ms_region_iy_m(1:2)
  integer        :: out_ms_region_iz_m(1:2)
  character(256) :: file_macropoint
  character(1)   :: set_ini_coor_vel
  integer        :: nmacro_write_group
  integer        :: nmacro_chunk
  real(8)        :: rmat_ms(3, 3)

!! &maxwell
  real(8)        :: al_em(3)
  real(8)        :: dl_em(3)
  integer        :: num_rgrid_em(3)
  real(8)        :: at_em
  real(8)        :: dt_em
  integer        :: nt_em
  character(8)   :: boundary_em(3,2)
  character(256) :: shape_file
  integer        :: media_num
  character(16)  :: media_type(0:200)
  real(8)        :: epsilon_em(0:200)
  real(8)        :: mu_em(0:200)
  real(8)        :: sigma_em(0:200)
  integer        :: pole_num_ld(0:200)
  real(8)        :: omega_p_ld(0:200)
  real(8)        :: f_ld(0:200,1:100)
  real(8)        :: gamma_ld(0:200,1:100)
  real(8)        :: omega_ld(0:200,1:100)
  character(16)  :: wave_input
  real(8)        :: ek_dir1(3)
  real(8)        :: source_loc1(3)
  real(8)        :: gbeam_sigma_plane1(3)
  real(8)        :: gbeam_sigma_line1(3)
  real(8)        :: ek_dir2(3)
  real(8)        :: source_loc2(3)
  real(8)        :: gbeam_sigma_plane2(3)
  real(8)        :: gbeam_sigma_line2(3)
  integer        :: obs_num_em
  integer        :: obs_samp_em
  real(8)        :: obs_loc_em(200,3)
  real(8)        :: obs_plane_ene_em(200,100)
  character(1)   :: yn_obs_plane_em(200)
  character(1)   :: yn_obs_plane_integral_em(200)
  character(1)   :: yn_wf_em
  real(8)        :: film_thickness
  integer        :: media_id_pml(3,2)
  integer        :: media_id_source1
  integer        :: media_id_source2
  real(8)        :: bloch_k_em(3)
  character(4)   :: bloch_real_imag_em(3)
  integer        :: ase_num_em
  real(8)        :: ase_ene_min_em
  real(8)        :: ase_ene_max_em
  real(8)        :: ase_wav_min_em
  real(8)        :: ase_wav_max_em
  integer        :: ase_smedia_id_em
  real(8)        :: ase_box_cent_em(3)
  real(8)        :: ase_box_size_em(3)
  integer        :: art_num_em
  real(8)        :: art_ene_min_em
  real(8)        :: art_ene_max_em
  real(8)        :: art_wav_min_em
  real(8)        :: art_wav_max_em
  integer        :: art_smedia_id_em
  real(8)        :: art_plane_bot_em(3)
  real(8)        :: art_plane_top_em(3)
  character(1)   :: yn_make_shape
  character(1)   :: yn_output_shape
  character(1)   :: yn_copy_x
  character(1)   :: yn_copy_y
  character(1)   :: yn_copy_z
  character(6)   :: rot_type
  integer        :: n_s
  character(32)  :: typ_s(1000)
  integer        :: id_s(1000)
  real(8)        :: inf_s(1000,10)
  real(8)        :: ori_s(1000,3)
  real(8)        :: rot_s(1000,3)
  
!! &analysis
  character(2)   :: projection_option
  real(8)        :: threshold_projection
  integer        :: nenergy
  real(8)        :: de
  integer        :: out_rt_energy_step
  character(1)   :: yn_out_psi
  character(1)   :: yn_out_dos
  character(1)   :: yn_out_dos_set_fe_origin
  real(8)        :: out_dos_start
  real(8)        :: out_dos_end
  integer        :: out_dos_nenergy
  real(8)        :: out_dos_width
  character(16)  :: out_dos_function
  character(1)   :: yn_out_pdos
  character(1)   :: yn_out_dns
  character(1)   :: yn_out_dns_rt
  character(1)   :: yn_out_dns_ac_je
  character(1)   :: yn_out_micro_je
  integer        :: out_dns_rt_step
  integer        :: out_dns_ac_je_step
  integer        :: out_micro_je_step
  character(1)   :: out_old_dns
  character(1)   :: yn_out_dns_trans
  real(8)        :: out_dns_trans_energy
  character(1)   :: yn_out_elf
  character(1)   :: yn_out_elf_rt
  integer        :: out_elf_rt_step
  character(1)   :: yn_out_estatic_rt
  integer        :: out_estatic_rt_step
  character(1)   :: yn_out_rvf_rt
  integer        :: out_rvf_rt_step
  character(1)   :: yn_out_tm
  character(1)   :: yn_out_gs_sgm_eps
  integer        :: out_gs_sgm_eps_mu_nu(2)
  real(8)        :: out_gs_sgm_eps_width
  integer        :: out_projection_step
  integer        :: out_projection_k_step
  integer        :: out_ms_step
  character(16)  :: format_voxel_data
  integer        :: nsplit_voxel_data
  character(1)   :: yn_lr_w0_correction
  character(1)   :: yn_out_intraband_current
  character(1)   :: yn_out_current_decomposed
  integer        :: out_current_decomposed_step
  integer        :: out_rt_spin_step
  character(1)   :: yn_out_mag_decomposed_rt
  character(1)   :: yn_out_mag_micro_rt
  character(1)   :: yn_out_spin_current_decomposed
  character(1)   :: yn_out_spin_current_micro
  character(1)   :: yn_out_rt_energy_components
  character(1)   :: yn_out_perflog
  character(6)   :: format_perflog ! 'stdout','text','csv'
  
!! &poisson
  integer        :: layout_multipole
  integer        :: num_multipole_xyz(3)
  integer        :: lmax_multipole
  real(8)        :: threshold_cg
  character(9)   :: method_poisson ! 'cg','ft','dirichlet'

!! &ewald
  integer        :: newald
  real(8)        :: aewald
  real(8)        :: cutoff_r
  real(8)        :: cutoff_r_buff
  real(8)        :: cutoff_g

!! &opt
  integer        :: nopt
  real(8)        :: max_step_len_adjust
  real(8)        :: convrg_opt_fmax
  character(5)   :: method_opt ! 'bfgs','steep','fire'
  real(8)        :: step_steep
  real(8)        :: step_fire

!! &md
  character(10)  :: ensemble
  character(20)  :: thermostat
  integer        :: step_velocity_scaling
  integer        :: step_update_ps
  real(8)        :: temperature0_ion_k
  character(1)   :: yn_set_ini_velocity
  character(256) :: file_ini_velocity
  real(8)        :: thermostat_tau
  character(1)   :: yn_stop_system_momt

!! &jellium
  character(1)   :: yn_jm
  character(1)   :: yn_charge_neutral_jm
  character(1)   :: yn_output_dns_jm
  character(256) :: shape_file_jm
  integer        :: num_jm
  real(8)        :: rs_bohr_jm(200)
  integer        :: sphere_nion_jm(200)
  real(8)        :: sphere_loc_jm(200,3)

!! &atomic_coor
!! &atomic_red_coor
integer,allocatable :: kion(:)
real(8),allocatable :: Rion(:,:)
real(8),allocatable :: Rion_red(:,:)
character(1),allocatable :: flag_opt_atom(:)
character(256),allocatable :: atom_name(:)

!! &code
  character(1) :: yn_want_stencil_hand_vectorization
  character(1) :: yn_want_communication_overlapping
  character(10) :: stencil_openmp_mode  ! 'auto', 'orbital', 'rgrid'
  character(10) :: current_openmp_mode  ! 'auto', 'orbital', 'rgrid'
  character(10) :: force_openmp_mode    ! 'auto', 'orbital', 'rgrid'

!! &band
  character(3) :: lattice
  integer :: nref_band
  real(8) :: tol_esp_diff
  integer :: num_of_segments
  integer :: ndiv_segment(max_num_of_segments)
  real(8) :: kpt(3,max_num_of_segments+1)
  character(1) :: kpt_label(max_num_of_segments+1)
  
  !! &sbe 
  character(1)   :: yn_vnl_correction
  integer        :: num_sbe
  character(256) :: sysname_sbe(1:200)
  integer        :: nk_sbe(1:200)
  integer        :: nstate_sbe(1:200)
  integer        :: nelec_sbe(1:200)
  real(8)        :: al_sbe(3,200)
  real(8)        :: al_vec1_sbe(3,200),al_vec2_sbe(3,200),al_vec3_sbe(3,200)
  integer        :: norder_correction
  real(8)        :: t2_sbe_fs
  real(8)        :: eg_ev
  real(8)        :: frozen_core_threshold_ev
  real(8)        :: frozen_free_threshold_ev
  ! Kuhn-Zurek/Caldeira-Leggett decoherence (Houston-basis wave-packet separation):
  ! lambda_au = kB_au * sbe_decoh_temperature_k / (hartree_kelvin_relationship * tau_m_au)
  real(8)        :: sbe_decoh_temperature_k
  real(8)        :: sbe_decoh_tau_m_fs
  ! 'y': ground-state input files come from a spinor (spin-orbit split) system:
  !      2*Nb spinor bands, occupation 1 per band, nelec valence bands (not nelec/2)
  character(1)   :: yn_sbe_spinor
  ! k-local impact-ionization Lindblad channel (Stobbe-Redmer-Schattke fit, GaAs):
  ! gamma(e_kin) = P * (e_kin - E_th)^4 * Theta(e_kin - E_th), e_kin from the CBM
  character(1)   :: yn_sbe_impact_ionization
  real(8)        :: sbe_ii_prefactor      ! P [s^-1 eV^-a]
  real(8)        :: sbe_ii_threshold_ev   ! E_th [eV]
  real(8)        :: sbe_ii_ramp_ev        ! linear Theta-smoothing width delta_E [eV]
  ! Impact-ionization fit-form switch: gamma = P (eps-E_th)^a Theta(eps-E_th).
  ! GaAs: Stobbe quartic (a=4, hard). Si: Keldysh quadratic (a=2, soft).
  character(32)  :: sbe_ii_form           ! 'stobbe_quartic' | 'keldysh_quadratic'
  real(8)        :: sbe_ii_exponent       ! fit exponent a (operative value)
  ! Coulomb (time-dependent Hartree-Fock / exchange) renormalization
  ! [Golde-Kira-Meier-Koch, Phys. Status Solidi B 248, 863 (2011)]:
  ! Sigma_nm(k) = - sum_{q/=k} V(k-q) rho_nm(q), V(p)=str*4pi/(eps*Omega*Nk*(p^2+kappa^2))
  character(1)   :: yn_sbe_coulomb
  real(8)        :: sbe_coulomb_epsilon    ! background dielectric constant eps
  real(8)        :: sbe_coulomb_strength   ! overall scaling of the exchange kernel
  real(8)        :: sbe_coulomb_screen_au  ! Yukawa screening kappa [1/Bohr]
  ! Hartree-Fock folding fix: project Sigma^HF block-diagonally onto the 4 FCC
  ! sublattice sectors (zero the spurious inter-sublattice exchange coupling
  ! that the cubic-cell band folding creates). [Popescu-Zunger PRB 85, 085201]
  character(1)   :: yn_sbe_hf_sublattice_proj

  ! Coset block-diagonal projection of the FIELD coupling (momentum matrix p,
  ! hence H_VG and the velocity). A translationally invariant perturbation
  ! conserves primitive crystal momentum, so its matrix elements between bands
  ! of DIFFERENT cosets (different primitive k folded to the same supercell k)
  ! are exactly zero; the band-folding/degeneracy mixing in the EPM eigenvectors
  ! creates spurious inter-coset coupling that artificially hybridizes the
  ! valleys (dense avoided crossings -> cascade tunneling Gamma~L~X). Projecting
  ! p block-diagonal over the cosets keeps rho block-diagonal and restores the
  ! correct per-valley (primitive) Zener physics. Same mechanism as the HF
  ! sublattice projection, applied to the core propagator. [Popescu-Zunger 2012]
  character(1)   :: yn_sbe_coset_proj

  ! Nonlocal "super-compute" mode (Part C): genuine momentum-exchange impact
  ! ionization + population-relaxing electron-phonon Lindblad. All OFF by
  ! default; the k-local channels remain the fast default. [scaffolding only at
  ! this stage -- the dynamics integration follows in later increments]
  character(1)   :: yn_sbe_superres        ! master switch for the nonlocal mode
  character(1)   :: yn_sbe_eph             ! population-relaxing e-ph Lindblad
  real(8)        :: sbe_eph_temperature_k  ! phonon bath T_ph [K] (Bose factor)
  real(8)        :: sbe_eph_nu_sat         ! collision-rate saturation nu_sat [s^-1]
  real(8)        :: sbe_eph_eps0_ev        ! nu(eps) saturation onset eps_0 [eV]
  real(8)        :: sbe_eph_n              ! nu(eps) shape exponent n
  character(1)   :: yn_sbe_bgr_threshold   ! density-dependent II threshold (BGR)
  real(8)        :: sbe_bgr_n_gate         ! apply BGR shift only above n [cm^-3]
  real(8)        :: sbe_bgr_coeff          ! BGR coefficient K [eV cm] (cube-root law)
  real(8)        :: sbe_search_sigma_e_ev  ! energy-bin width sigma_E [eV] (<=0: grid-matched)
  real(8)        :: sbe_ring_gate_fs       ! ring virtual-transient gate time [fs] (0: off = default; <0: auto 2*pi/Egap; >0 manual). EXPERIMENTAL: in sub-cycle fields the dressing follows the envelope and outlives any lifetime gate.
  real(8)        :: sbe_eph_interband_scale ! calibration factor on the eph ring's gap-straddling (phonon-assisted BTBT) partial rates; 1 = cited nu_sat upper estimate
  ! SFSB non-Markovian heat bath (Boroumand et al., Rep. Prog. Phys. 88, 070501
  ! (2025); wiki/10 sec. 6): second-order memory-integral ionization mode with
  ! the bath correlation kernel C(t) instead of the Markovian T2 dephasing.
  character(1)   :: yn_sbe_sfsb            ! 'y': run the SFSB memory-integral ionization mode (1D k-line, E || b1)
  character(16)  :: sbe_bath_model         ! bath spectral model: 'none'|'ohmic'|'debye'|'rta' [B25 Eq. (5)]
  real(8)        :: sbe_bath_jo            ! dimensionless bath coupling jo [B25]
  real(8)        :: sbe_bath_wc_ev         ! bath cutoff hbar*wc [eV]; required > 0 for ohmic/debye
  real(8)        :: sbe_bath_temperature_k ! bath (electron/lattice) temperature [K]; 0 = T=0 allowed
  real(8)        :: sbe_bath_rta_t2_fs     ! RTA T2 [fs]; <=0: derived T2 = hbar/(2 pi kB T jo) [B25 sec. 2]
  real(8)        :: sbe_bath_memory_fs     ! memory-kernel truncation window [fs]; <=0: full history
  character(1)   :: yn_sbe_bath_imc        ! 'n': zero Im C(t) -- diagnostic only [B25 Fig. 3(c)]
  integer        :: sbe_sfsb_nv            ! # top valence bands in the SFSB pair sum
  integer        :: sbe_sfsb_nc            ! # bottom conduction bands in the SFSB pair sum
  integer        :: sbe_sfsb_stride        ! integrate every Nth field sample; 0: auto from the Stark-shifted gap
  ! Carrier-carrier (e-e/e-h) scattering (Part F): CPTP relaxation of the
  ! adiabatic populations toward a Fermi-Dirac with the same number AND energy
  ! (intraband thermalization + EID), at a screened-Coulomb rate.
  character(1)   :: yn_sbe_eeh             ! enable carrier-carrier channel
  real(8)        :: sbe_eeh_nu_sat         ! carrier-carrier rate scale [s^-1]
  ! Auger recombination (Sec 13): density-gated, number-conserving CPTP channel.
  ! Per-carrier rate gamma = C n^2 (so R = C n^3); C cited per material.
  character(1)   :: yn_sbe_auger           ! enable Auger recombination channel
  real(8)        :: sbe_auger_c_cm6s       ! Auger coeff C [cm^6/s]; <=0: material default
  real(8)        :: sbe_auger_n_gate_cm3   ! activation density [cm^-3]; <=0: material default
  ! ring II/Auger refinements (wiki/00 proposed block, maintainer-approved 2026-07-04)
  real(8)        :: sbe_ring_vq_floor      ! B3: skip quadruples with vq < floor*max(vq); 0 = off
  character(1)   :: yn_sbe_ii_fk_soften    ! A5: Franz-Keldysh field-softened II threshold
  real(8)        :: sbe_ii_fk_mu           ! A5: reduced mass [m_e] of hbar*theta=(F^2/2mu)^(1/3); required >0 when on
  real(8)        :: sbe_ii_phassist        ! A1: phonon-assisted II/Auger sideband strength; 0 = off
  character(1)   :: yn_sbe_ii_holes        ! A2: hole-initiated II + its Auger reverse (Cp/Cn from registry)
  character(1)   :: yn_sbe_eph_acoustic    ! A4: quasi-elastic acoustic deformation e-ph mode
  real(8)        :: sbe_eph_ac_xi_ev       ! A4: acoustic deformation potential OVERRIDE [eV];
                                           !     <=0 -> registry default. Graphene: substrate-
                                           !     dependent (Hwang-Das Sarma D=16 eV is the upper
                                           !     literature scale; on some substrates it is lower).
  ! SBE checkpoint/restart (B4)
  integer        :: sbe_checkpoint_step    ! write rho checkpoint every N steps; 0 = off
  character(1)   :: yn_sbe_checkpoint_restart ! resume from SYSNAME_sbe_ckpt_rank*.bin

  !! &epm (local empirical pseudopotential method, Cohen-Bergstresser)
  character(32)  :: epm_material
  real(8)        :: epm_lattice_constant_au
  real(8)        :: epm_pw_cutoff_ry
  character(16)  :: epm_cell        ! 'primitive' (default) | 'folded' (cubic supercell, GaAs/Si only)

  !! &dc
  integer        :: num_fragment(3)
  integer        :: num_rgrid_buffer(3)
  integer        :: nproc_rgrid_tot(3)
  character(256) :: file_atom_coor_frag
  real(8)        :: xi_dc
  character(1)   :: yn_dc_lcfo
  character(1)   :: yn_dc_lcfo_diag
  integer        :: nstate_frag
  real(8)        :: energy_cut
  real(8)        :: lambda_cut

end module salmon_global
