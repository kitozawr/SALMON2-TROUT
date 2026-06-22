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

### Super-compute / Parts C–G (✅ implemented; all default OFF)
| Parameter | Units | Default | Description |
|---|---|---|---|
| yn_sbe_superres | — | 'n' | nonlocal super-compute master switch: ring-MPI Σ^HF (D) + nonlocal-partner II (C4). |
| yn_sbe_eph | — | 'n' | electron-phonon population-relaxing Lindblad (C5; full phonon table). Toggle Zurek off. |
| sbe_eph_temperature_k | K | 300.0 | phonon bath T_ph for N_B. |
| sbe_eph_nu_sat | s⁻¹ | material | saturation rate (Si 1.3e14, GaAs 1e14). |
| sbe_eph_eps0_ev | eV | 0.8 | saturation onset ε₀ in ν(ε)=ν_sat[1−exp(−(ε/ε₀)^n)]. |
| sbe_eph_n | — | 2 | saturation shape exponent n. |
| yn_sbe_eeh | — | 'n' | carrier-carrier (e-e/e-h) CPTP thermalization to a Fermi-Dirac (F). |
| sbe_eeh_nu_sat | s⁻¹ | 1e14 | carrier-carrier rate scale. |
| yn_sbe_bgr_threshold | — | 'n' | density-dependent II threshold E_th(t)=E_th0−|ΔE_BGR(n)| (C7). |
| sbe_bgr_n_gate | cm⁻³ | 5.0e18 | apply BGR shift only above this density. |
| sbe_bgr_coeff | eV·cm | 1.9e-8 | BGR coefficient K (tunable [1.9,3.8]e-8). |
| yn_sbe_hf_sublattice_proj | — | 'y' | project Σ^HF block-diagonally onto 4 FCC sublattices (folding fix, E). |
| sbe_search_sigma_e_ev | eV | grid-matched | energy-bin width σ_E for the final-state search (C3). |

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

## Example: Silicon super-compute mode ✅ (Parts C–G implemented; physical validation pending)
```fortran
&sbe
  yn_sbe_superres          = 'y'   ! ring-MPI Sigma^HF + nonlocal-partner II
  yn_sbe_eph               = 'y'   ! e-ph cooling (drives bleaching)
  sbe_eph_nu_sat           = 1.3d14
  yn_sbe_eeh               = 'y'   ! carrier-carrier thermalization
  yn_sbe_impact_ionization = 'y'
  sbe_ii_form              = 'keldysh_quadratic'
  sbe_ii_exponent          = 2
  sbe_ii_threshold_ev      = 1.1d0
  yn_sbe_bgr_threshold     = 'y'
  sbe_decoh_temperature_k  = -1.0d0   ! Zurek off (e-ph provides decoherence)
/
```
All channels are CPTP and gated OFF by default. Validation staging (Chefonov Si THz bleaching, not yet run here): (1) ~8.5% bleaching plateau at ~5 MV/cm with e-ph alone; (2) enable II, ~2× transmission drop at >10–15 MV/cm. [Chefonov et al., PRB 98, 165206 (2018)]

## Recipes by material & mode ✅

All recipes assume the matching EPM ground state was generated first (see the
GaAs / Si EPM examples above). Every flag below defaults OFF; turning none of
them on reproduces the byte-for-byte legacy GaAs run. The `&sbe` block can be
combined freely — the recipes just bundle the physically-sensible defaults.
Constants and their primary-source citations: [Constants](02_constants.md).

> **Important:** when `yn_sbe_eph='y'`, the SBE-step input must ALSO carry an
> `&epm` block with the right `epm_material` — the electron-phonon channel
> selects its phonon table (Si: 6 intervalley g/f modes; GaAs: Fröhlich LO + 5
> intervalley) from `epm_material`, which **defaults to `'GaAs'`** if absent.
> A runnable two-step example lives in
> [`../samples/exercise_x3_bulkSi_epm_bloch_superres/`](../samples/exercise_x3_bulkSi_epm_bloch_superres/).

### GaAs — per-channel quick reference
| Mode | Minimal `&sbe` | Material defaults (auto if unset) |
|---|---|---|
| Coulomb TDHF only | `yn_sbe_coulomb='y'` | `sbe_coulomb_epsilon=12.9` (ε∞ GaAs) |
| Impact ionization (quartic) | `yn_sbe_impact_ionization='y'`, `sbe_ii_form='stobbe_quartic'`, `sbe_ii_exponent=4`, `sbe_ii_threshold_ev=2.1` | hard threshold ~1.5·E_g |
| e-ph cooling | `yn_sbe_eph='y'` | Fröhlich LO 36 meV + 5 intervalley; `ν_sat=1e14 s⁻¹` |
| carrier-carrier | `yn_sbe_eeh='y'` | `ν_cc=1e14 s⁻¹`; ε for screening = 12.9 |

### GaAs — full super-compute recipe
```fortran
&epm
  epm_material            = 'GaAs'    ! zincblende, V^A≠0
  epm_lattice_constant_au = 10.68d0
/
&sbe
  yn_sbe_superres          = 'y'   ! ring-MPI Sigma^HF + nonlocal-partner II
  yn_sbe_coulomb           = 'y'   ! TDHF exchange (extreme-THz regime)
  sbe_coulomb_epsilon      = 12.9d0
  yn_sbe_eph               = 'y'   ! Fröhlich-LO-dominated cooling
  sbe_eph_nu_sat           = 1.0d14
  sbe_eph_temperature_k    = 300.0d0
  yn_sbe_eeh               = 'y'   ! carrier-carrier thermalization
  sbe_eeh_nu_sat           = 1.0d14
  yn_sbe_impact_ionization = 'y'
  sbe_ii_form              = 'stobbe_quartic'   ! GaAs hard threshold, a=4
  sbe_ii_exponent          = 4
  sbe_ii_threshold_ev      = 2.1d0
  yn_sbe_bgr_threshold     = 'y'   ! density-dependent gap shrinkage
  sbe_decoh_temperature_k  = -1.0d0   ! Zurek off (e-ph provides decoherence)
/
```

### Si — per-channel quick reference
| Mode | Minimal `&sbe` | Material defaults (auto if unset) |
|---|---|---|
| Impact ionization (soft) | `yn_sbe_impact_ionization='y'`, `sbe_ii_form='keldysh_quadratic'`, `sbe_ii_exponent=2`, `sbe_ii_threshold_ev=1.1` | soft near-gap threshold |
| Impact ionization (full-band) | as above but `sbe_ii_exponent=4.6` | Kamakura full-band fit |
| e-ph cooling | `yn_sbe_eph='y'` | 6 intervalley g/f phonons; `ν_sat=1.3e14 s⁻¹` |
| carrier-carrier | `yn_sbe_eeh='y'` | `ν_cc=1e14 s⁻¹`; ε for screening = 11.7 |

### Si — full super-compute recipe
```fortran
&epm
  epm_material            = 'Si'     ! diamond, V^A=0, Kunikiyo form factors
  epm_lattice_constant_au = 10.26d0
/
&sbe
  yn_sbe_superres          = 'y'
  yn_sbe_eph               = 'y'   ! intervalley-dominated cooling
  sbe_eph_nu_sat           = 1.3d14
  sbe_eph_temperature_k    = 300.0d0
  yn_sbe_eeh               = 'y'
  yn_sbe_impact_ionization = 'y'
  sbe_ii_form              = 'keldysh_quadratic'   ! Si soft threshold, a=2
  sbe_ii_exponent          = 2
  sbe_ii_threshold_ev      = 1.1d0
  yn_sbe_bgr_threshold     = 'y'
  sbe_decoh_temperature_k  = -1.0d0
/
```

### Maxwell-SBE multiscale with the new channels ✅
The multiscale driver (`theory='maxwell_sbe'`, `src/ssbe/multiscale_ssbe.f90`)
runs one independent SBE cell per macropoint, driven by that point's
macroscopic Maxwell (Weyl-gauge FDTD) field. Every `&sbe` channel above works
unchanged — just add the multiscale grid / media blocks and the `&sbe`
material-pointer parameters (`num_sbe`, `sysname_sbe(:)`, `nk_sbe(:)`,
`nstate_sbe(:)`, `nelec_sbe(:)`, `al_vec*_sbe`). The nonlocal reductions
(Coulomb ring, BGR, nonlocal-II) are confined to each macropoint's own k-grid,
so momentum exchange stays local to the cell. Banner diagnostics print once
(global rank 0, first macropoint).
```fortran
&calculation
  theory = 'maxwell_sbe'
/
&control
  sysname = 'Si_ms'
/
&multiscale
  fdtddim     = '1d'
  nx_m        = 2000      ! Maxwell propagation cells
  ny_m        = 1
  nz_m        = 1
  hx_m        = 250.0d0   ! [a.u.] coarse Maxwell grid spacing
  hy_m        = 250.0d0
  hz_m        = 250.0d0
/
&sbe
  num_sbe        = 1               ! one SBE material in this run
  sysname_sbe(1) = 'Si_cubic'      ! reads the Si_cubic EPM output as the cell
  nk_sbe(1)      = 64              ! must match the EPM ground state (4x4x4)
  nstate_sbe(1)  = 32
  nelec_sbe(1)   = 32
  ! --- new dissipation channels, applied per macropoint ---
  yn_sbe_eph               = 'y'   ! e-ph cooling -> Drude-conductivity bleaching
  sbe_eph_nu_sat           = 1.3d14
  yn_sbe_impact_ionization = 'y'   ! enable for the high-field transmission drop
  sbe_ii_form              = 'keldysh_quadratic'
  sbe_ii_exponent          = 2
  sbe_ii_threshold_ev      = 1.1d0
  sbe_decoh_temperature_k  = -1.0d0
/
```
Run on ≥ nmacro MPI ranks for one-macropoint-per-rank; with fewer ranks the
k-points of each macropoint are split across the ranks sharing it (the ring
and reductions then run over that per-macropoint group `icomm_macro`).

### Unit conventions for the new parameters
Inputs are in the named units; the solver converts to atomic units internally
(audited): rates `s⁻¹ → a.u.⁻¹` (`×au_fs·1e-15`), energies `eV → Ha` (`/au_ev`),
densities `cm⁻³` (gate), deformation potentials from the cited tables. You
never pass a.u. directly for these knobs. See [Constants](02_constants.md).

## Building ✅
```sh
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"
```
For production (and the ring-pipeline super-mode): `-D USE_MPI=ON` with an MPI Fortran compiler.

## Tests ✅
Self-contained Python tests in [`../tests/`](../tests/) (reuse the standalone EPM machinery; no SALMON build needed for most). Run all: `python3 tests/run_all.py`.
