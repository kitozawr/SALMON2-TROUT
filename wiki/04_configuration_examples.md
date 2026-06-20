# Configuration & Examples

> All namelist parameters and runnable examples. ✅ = implemented & in code; 🚧 = planned. Physics in [Physics Models](01_physics_models.md); defaults in [Constants](02_constants.md).

## &sbe parameters

### Existing (✅ implemented, unchanged)
| Parameter | Units | Default | Description |
|---|---|---|---|
| sbe_decoh_temperature_k | K | −1.0 | Bath T for Kuhn-Zurek λ=k_B T/τ_m. Both this and τ_m > 0 to enable. |
| sbe_decoh_tau_m_fs | fs | −1.0 | Momentum-relaxation time τ_m. |
| frozen_core_threshold_ev | eV | 0.0 | Freeze bands below E_F + threshold (negative, e.g. −15). |
| frozen_free_threshold_ev | eV | 0.0 | Freeze bands above E_F + threshold (positive, e.g. +20). |
| yn_sbe_spinor | — | 'n' | 'y': spinor (spin-orbit) input. Combine with yn_vnl_correction='y'. |
| yn_sbe_impact_ionization | — | 'n' | 'y': k-local impact-ionization Lindblad channel. |
| sbe_ii_prefactor | s⁻¹eV⁻ᵃ | 2.0e12 | Fit prefactor P. |
| sbe_ii_threshold_ev | eV | 2.1 | Threshold E_th above field-free CBM. |
| sbe_ii_ramp_ev | eV | 0.2 | Linear Θ-smoothing width; ≤0 = hard step. |
| yn_sbe_coulomb | — | 'n' | 'y': Coulomb TDHF exchange (Golde-Kira-Meier-Koch). Non-k-local. |
| sbe_coulomb_epsilon | — | 12.9 | Background dielectric (GaAs). |
| sbe_coulomb_strength | — | 1.0 | Exchange-kernel scaling (0 disables; >1 enhances). |
| sbe_coulomb_screen_au | Bohr⁻¹ | 0.0 | Yukawa screening κ; 0 = bare with q=0 excluded. |

### Part B (✅ implemented)
| Parameter | Units | Default | Description |
|---|---|---|---|
| sbe_ii_form | — | 'stobbe_quartic' | II fit form: 'stobbe_quartic' (GaAs, a=4) or 'keldysh_quadratic' (Si, a=2). |
| sbe_ii_exponent | — | 4.0 | Fit exponent a (operative; 2 for Si soft, 4.6 for Si full-band). |

### Super-compute / Part C–E (🚧 planned — flags reserved)
| Parameter | Units | Default | Description |
|---|---|---|---|
| yn_sbe_superres | — | 'n' | 🚧 nonlocal super-compute mode (nonlocal II + e-ph, ring MPI). |
| yn_sbe_eph | — | 'n' | 🚧 electron-phonon Lindblad (relaxes populations). Toggle Zurek off. |
| sbe_eph_temperature_k | K | 300.0 | 🚧 phonon bath T_ph for N_B. |
| sbe_eph_nu_sat | s⁻¹ | material | 🚧 saturation rate (Si 1.3e14, GaAs 1e14). |
| sbe_eph_eps0_ev | eV | 0.8 | 🚧 saturation onset ε₀. |
| sbe_eph_n | — | 2 | 🚧 saturation shape exponent n. |
| yn_sbe_bgr_threshold | — | 'n' | 🚧 density-dependent II threshold. |
| sbe_bgr_n_gate | cm⁻³ | 5.0e18 | 🚧 apply BGR shift only above this density. |
| sbe_bgr_coeff | eV·cm | 1.9e-8 | 🚧 BGR coefficient K (tunable [1.9,3.8]e-8). |
| yn_sbe_hf_sublattice_proj | — | 'y' | 🚧 (Part E) project Σ^HF block-diagonally onto 4 FCC sublattices. |
| sbe_search_sigma_e_ev | eV | grid-matched | 🚧 energy-bin width σ_E for the final-state search. |

## &epm parameters ✅
| Parameter | Units | Default | Description |
|---|---|---|---|
| epm_material | — | 'GaAs' | 'GaAs' (zincblende) or 'Si' / 'Si_cb' (diamond, V^A=0). |
| epm_lattice_constant_au | Bohr | 10.68 | Lattice constant a (Si: 10.26). |
| epm_pw_cutoff_ry | Ry | 11.1 | Plane-wave cutoff |k+G|². |

## &analysis (output cadence) ✅
| Parameter | Default | Description |
|---|---|---|
| out_rt_energy_step | 10 | Stride for SYSNAME_sbe_rt_energy.data + stdout. |
| out_projection_step | 100 | Stride for SYSNAME_sbe_nex.data (excited e/h, summed over k). |
| out_projection_k_step | 1000 | Stride for SYSNAME_sbe_nex_k.data (k-resolved Houston-basis lowest-CB). |

## Example: GaAs EPM → SBE pipeline ✅
```fortran
! Step 1: ground state
&calculation
  theory = 'epm'
/
&epm
  epm_material            = 'GaAs'
  epm_lattice_constant_au = 10.68d0
  epm_pw_cutoff_ry        = 11.1d0
/
```
Then `theory='sbe'` reading the generated files (sysname, lattice, num_kgrid, nstate, nelec must match).

## Example: Silicon EPM → SBE with soft-threshold impact ionization ✅ (Parts A,B)
```fortran
! Step 1: Si ground state (scalar, 4x4x4)
&calculation
  theory = 'epm'
/
&control
  sysname = 'Si_cubic'
/
&system
  yn_periodic = 'y'
  al(1:3) = 10.26d0, 10.26d0, 10.26d0
  nelec  = 32
  nstate = 32
/
&kgrid
  num_kgrid(1:3) = 4, 4, 4
/
&epm
  epm_material            = 'Si'        ! diamond, V^A=0, Kunikiyo form factors
  epm_lattice_constant_au = 10.26d0     ! 5.431 Angstrom
  epm_pw_cutoff_ry        = 11.1d0
/
```
```fortran
! Step 2: SBE with Si soft-threshold (quadratic) impact ionization
&calculation
  theory = 'sbe'
/
&sbe
  yn_sbe_impact_ionization = 'y'
  sbe_ii_form              = 'keldysh_quadratic'   ! Si soft threshold
  sbe_ii_exponent          = 2
  sbe_ii_threshold_ev      = 1.1d0
  sbe_ii_ramp_ev           = 0.2d0
/
```

## Example: Silicon super-compute mode 🚧 (Part C — not yet implemented)
Target (Chefonov THz bleaching): `yn_sbe_superres='y'`, `yn_sbe_eph='y'`, `sbe_eph_nu_sat=1.3d14`, Zurek off, `yn_sbe_bgr_threshold='y'`. Validation staging: (1) ~8.5% bleaching plateau at ~5 MV/cm with e-ph alone; (2) enable II, ~2× transmission drop at >10–15 MV/cm. [Chefonov et al., PRB 98, 165206 (2018)]

## Building ✅
```sh
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"
```
For production (and the ring-pipeline super-mode): `-D USE_MPI=ON` with an MPI Fortran compiler.

## Tests ✅
Self-contained Python tests in [`../tests/`](../tests/) (reuse the standalone EPM machinery; no SALMON build needed for most). Run all: `python3 tests/run_all.py`.
