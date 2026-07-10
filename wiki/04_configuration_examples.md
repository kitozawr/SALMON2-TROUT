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
| yn_sbe_auger | — | 'n' | Auger recombination, number-conserving CPTP. **With the ring** (`yn_sbe_superres='y'`+`yn_sbe_impact_ionization='y'`) it is the exact **time-reverse of the inter-k impact ionization** — same \|M\|²/ε(q)/umklapp weight, reversed Fermi factors, **no separate C** (rate scale = the cited II magnitude, detailed balance). **Ring off:** the k-local density-gated γ=C·n² (R=C·n³) fallback, which requires an explicit `sbe_auger_c_cm6s` (no verified default). See wiki/07. |
| sbe_auger_c_cm6s | cm⁶/s | −1.0 | Auger coefficient C; ≤0 ⇒ material default, but **no material has a verified default** (the former CdS 2.0e-30 carried a fabricated citation, removed) — must be set explicitly. |
| sbe_auger_n_gate_cm3 | cm⁻³ | −1.0 | activation density; ≤0 ⇒ material default (CdS 1e18 [Shah 1986], used only with a user-supplied C). |
| yn_sbe_bgr_threshold | — | 'n' | density-dependent II threshold E_th(t)=E_th0−|ΔE_BGR(n)| (C7). |
| sbe_bgr_n_gate | cm⁻³ | 5.0e18 | apply BGR shift only above this density. |
| sbe_bgr_coeff | eV·cm | 1.9e-8 | BGR coefficient K (tunable [1.9,3.8]e-8). |
| yn_sbe_hf_sublattice_proj | — | 'y' | project Σ^HF block-diagonally onto 4 FCC sublattices (folding fix, E). |
| sbe_search_sigma_e_ev | eV | grid-matched | energy-bin width σ_E for the final-state search (C3). |

## &epm parameters ✅
| Parameter | Units | Default | Description |
|---|---|---|---|
| epm_material | — | 'GaAs' | **Fortran `theory='epm'` (cubic only):** 'GaAs' (zincblende) / 'Si' / 'Si_cb' (diamond, V^A=0). Also the **dissipation-channel material pointer** in the `&sbe` step — here 'CdS' is additionally valid (registry). |
| epm_lattice_constant_au | Bohr | 10.68 | Lattice constant a (Si: 10.26). |
| epm_pw_cutoff_ry | — | 11.1 | Plane-wave cutoff on \|G\|² in **(2π/a)² units** (integer shells h²+k²+l²); same convention as the Python reference. |

> **Non-cubic materials (CdS, graphene) use the Python EPM references**, not the
> Fortran `theory='epm'` (which is cubic-only). Generate their ground state +
> band path + 2-coset unfold map with `python3 epm_wurtzite_cds.py` /
> `python3 epm_graphene.py` (emits `SYSNAME_k/_eigen/_tm.data`, `_bandpath.data`,
> `_unfold.data`), then run `theory='sbe'` on the result. The cubic Fortran EPM
> is **verified byte-equivalent** to the Python reference for GaAs/Si (scalar).

## &analysis (output cadence) ✅
| Parameter | Default | Description |
|---|---|---|
| out_rt_energy_step | 10 | Stride for SYSNAME_sbe_rt_energy.data + stdout. |
| out_projection_step | 100 | Stride for SYSNAME_sbe_nex.data (excited e/h, summed over k). |
| out_projection_k_step | 1000 | Stride for the k-resolved population maps (below). |
| yn_out_intraband_current | n | If `y`, write SYSNAME_sbe_intra_current.data: the intra-band (drift) current in the Houston basis, every step. |

**Per-k population maps** (written every `out_projection_k_step`):
- **SYSNAME_sbe_nex_k_real.data** — *real carriers only*: the fixed-basis (diabatic) lowest-CB occupation, i.e. the k-resolved excited-electron count n_ex. No reversible A²(t) virtual-polarization breathing; accumulates monotonically and freezes when the field passes. **This is the carrier map to use.** (`_unfold_real` twin per primitive BZ point when a `_unfold.data` map is present; N cosets — 4 cubic / 2 wurtzite-rectangular.) A small residual (~20 %) can remain at the most strongly interband-coupled folds (e.g. the GaAs L-valley); it is 0 % at Γ and 0 % in the BZ total.
- **SYSNAME_sbe_nex_k.data** — instantaneous Houston-basis lowest-CB population (`_unfold` twin). Physical *during* the pulse but carries the reversible virtual breathing (∝A(t)²); equals the real map after the pulse. Kept for diagnostics.

**Intra-band current** (`yn_out_intraband_current='y'`): in the velocity gauge only the total current (SYSNAME_sbe_rt.data) is gauge invariant; its intra/inter split is physical in the Houston basis. J_intra is the Boltzmann drift (vanishes when the field is off); J_inter = J_total − J_intra is the interband polarization. [T. Otobe, PRB 94, 235152 (2016)]

**Plotting:** `python3 plot_sbe_results.py -i <dir> -o <outdir>` renders the RT
observables, conductivity, the intra-band current, the **real-carrier** k–t
population maps (folded + unfolded), and band structure. The instantaneous
(breathing) Houston maps are plotted only with `--instantaneous` (or when no
`_real` file is present). Use `--lattice wurtzite` for CdS; the clean primitive
band path (`SYSNAME_bandpath.data`, material-agnostic) plots via `--only-bands`.

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

## Runnable examples per material (`samples/`) ✅

Self-contained two-step EPM→SBE exercises, one folder per material, each with a
`README.md` and the input file(s). The `x`-series are the TROUT (Bloch/EPM-SBE)
additions. **Naming:** GaAs & Si use the cubic Fortran EPM (`theory='epm'`,
self-contained `.inp`); CdS & graphene use the **Python EPM** for step 1 (the
Fortran EPM is zincblende/diamond only). All run on a sparse 4×4×4 (graphene
4×4×1) grid and demonstrate the real-carrier maps + the Houston intra-band
current.

| Exercise | Material | Step 1 (ground state) | Highlights |
|---|---|---|---|
| [`exercise_x3_bulkSi_epm_bloch_superres`](../samples/exercise_x3_bulkSi_epm_bloch_superres/) | Si (diamond) | `salmon < Si_epm_gs.inp` | super-compute mode, all CPTP dissipation channels |
| [`exercise_x4_GaAs_epm_bloch_realcarrier`](../samples/exercise_x4_GaAs_epm_bloch_realcarrier/) | GaAs (zincblende) | `salmon < GaAs_epm_gs.inp` | real-carrier maps (Γ populated, no breathing) + intra-band current |
| [`exercise_x5_CdS_wurtzite_epm_bloch`](../samples/exercise_x5_CdS_wurtzite_epm_bloch/) | CdS (wurtzite) | `python3 epm_wurtzite_cds.py gs` | 2-coset unfold, direct 2.5 eV gap |
| [`exercise_x6_graphene_epm_bloch`](../samples/exercise_x6_graphene_epm_bloch/) | graphene (π-model) | `python3 epm_graphene.py gs` | 2-coset unfold, Dirac carriers, **in-plane** field |

The Python EPM `gs` mode (`epm_wurtzite_cds.py gs`, `epm_graphene.py gs`) emits
`SYSNAME_k/_eigen/_tm/_unfold/_bandpath.data` into the working directory without
the slow convergence validation (the Hamiltonian build is vectorized: CdS GS
≈ 30 s, graphene ≈ 1 s). Run it where you run `salmon`. Then:

```sh
./build/salmon < <material>_sbe_rt.inp          # step 2
python3 plot_sbe_results.py -i . -o plots --snapshots
```

Each SBE input sets `yn_out_intraband_current='y'` and runs past the pulse so
the real-carrier maps (`*_sbe_nex_k_real.data`, `*_unfold_real.data`) settle to
their field-free residual. All four conserve the trace (GaAs/Si/CdS = 32,
graphene = 4) to machine precision. Verified end-to-end with the current build.

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

### CdS (wurtzite) — Python EPM → SBE ✅
CdS is **not** a Fortran-EPM material; generate its ground state with the Python
reference, which emits the orthorhombic `al(1:3)=(a, a√3, c)` cell dataset:
```sh
python3 epm_wurtzite_cds.py        # -> CdS_k/_eigen/_tm.data, _bandpath.data, _unfold.data
```
The orthorhombic cell vector is `al(1:3) = (7.8159, 13.5375, 12.6852)` Bohr.
Per-channel provenance (registry): e-ph (Fröhlich LO 38 meV), Coulomb (ε₀=8.9),
and impact ionization (E_th=3.6 eV; **you must set `sbe_ii_prefactor`** — no cited
CdS value) are cited/enabled. **Carrier-carrier and Auger are both FORBIDDEN** for
CdS (no cited e-e rate; and no verified Auger C — the former "Haury 1998" value was
a fabricated citation, removed): `yn_sbe_eeh='y'` and `yn_sbe_auger='y'` both abort.
```fortran
&calculation
  theory = 'sbe'
/
&control
  sysname = 'CdS'
/
&system
  yn_periodic = 'y'
  al(1:3) = 7.8159d0, 13.5375d0, 12.6852d0   ! (a, a*sqrt3, c)
  nelec  = 32
  nstate = 32
/
&kgrid
  num_kgrid(1:3) = 4, 4, 4
/
&epm
  epm_material = 'CdS'             ! dissipation-channel pointer (wurtzite registry)
/
&sbe
  yn_sbe_eph    = 'y'             ! Fröhlich polar-optical LO (primary room-T channel)
  yn_sbe_coulomb = 'y'            ! auto eps0 = 8.9
  yn_sbe_auger  = 'y'            ! Auger recombination (auto C=2e-30, n_gate=1e18; rare)
  ! yn_sbe_impact_ionization = 'y'  ! needs sbe_ii_prefactor (uncited) -> set explicitly
  ! yn_sbe_eeh = 'y'               ! FORBIDDEN for CdS (no cited e-e rate) -> aborts
/
```

### graphene (monolayer, π-model) — Python EPM → SBE ✅ (bands / clean dynamics)
graphene uses the Python reference (rectangular 4-atom cell, 2D sheet in a 3D
vacuum box); it is a **minimal π-model** (1 π electron/atom → `nelec=4`,
`nstate=8`, Dirac cone = lowest band pair).
```sh
python3 epm_graphene.py            # -> graphene_k/_eigen/_tm.data, _bandpath.data, _unfold.data
```
`al(1:3) = (a, √3a, vacuum) = (4.6487, 8.0518, 37.7945)` Bohr; `num_kgrid(3)=1`
(no dispersion along the vacuum axis). **No dissipation channels are wired for
graphene yet** (its registry entry — e-ph E2g/A1', gapless-CM Auger, no-Kuhn-Zurek
policy — is a TODO; any material-dependent channel currently aborts). Use it for
the **clean** (no-dissipation) SBE and the band structure:
```fortran
&calculation
  theory = 'sbe'
/
&control
  sysname = 'graphene'
/
&system
  yn_periodic = 'y'
  al(1:3) = 4.6487d0, 8.0518d0, 37.7945d0   ! (a, sqrt3*a, vacuum)
  nelec  = 4
  nstate = 8
/
&kgrid
  num_kgrid(1:3) = 4, 4, 1
/
&sbe
  ! clean run (Dirac-cone dynamics); dissipation channels not yet wired
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


---

# APPENDIX — consolidated parameter reference & pipelines (moved verbatim from the README, 2026-07-04 restructure)

**This appendix is the canonical, most-current parameter reference** — where the
older per-part tables above and this appendix disagree, the appendix wins (it
tracked every merged channel through PR #65). The README now keeps only a short
configuration philosophy + a pointer here.

## Configuration Parameters

The `&sbe` namelist now accepts the following parameters:

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `sbe_decoh_temperature_k` | K | `-1.0d0` | Bath temperature $T$ for the Kuhn-Zurek/Caldeira-Leggett dephasing rate $\lambda=k_B T/\tau_m$. Both this and `sbe_decoh_tau_m_fs` must be `> 0` to enable decoherence. |
| `sbe_decoh_tau_m_fs` | fs | `-1.0d0` | Wave-packet momentum-relaxation time $\tau_m$ entering $\lambda=k_B T/\tau_m$. |
| `frozen_core_threshold_ev` | eV | `0.0d0` | Freeze bands below $E_F + \text{threshold}$. (Use negative values, e.g., `-15.0`). |
| `frozen_free_threshold_ev` | eV | `0.0d0` | Freeze bands above $E_F + \text{threshold}$. (Use positive values, e.g., `+20.0`). |
| `yn_sbe_spinor` | — | `'n'` | `'y'`: ground-state input files come from a **spinor (spin-orbit split)** system — occupation 1 per spinor band, `nelec` valence bands instead of `nelec/2`. Combine with `yn_vnl_correction='y'` when the dataset carries the $\hat v_{SO}=\nabla_k\hat H_{SO}$ correction in `rvnl_tm`. |
| `yn_sbe_impact_ionization` | — | `'n'` | `'y'`: enable the **k-local impact-ionization** Lindblad channel. Fully optional; threshold-gated, so it costs ~nothing while no populated branch exceeds $E_{\rm th}$. The fit (form, exponent, prefactor, threshold) auto-selects from the material registry — see below. |
| `sbe_ii_prefactor` | s⁻¹eV⁻ᵃ | `-1.0d0` | Fit prefactor $P$ in $\gamma=P(\varepsilon^{\rm kin}-E_{\rm th})^a$. `≤ 0` ⇒ material default. |
| `sbe_ii_threshold_ev` | eV | `-1.0d0` | Ionization threshold $E_{\rm th}$ above the field-free CBM. `< 0` ⇒ material default (GaAs `2.1`, Si `1.1`). |
| `sbe_ii_ramp_ev` | eV | `0.2d0` | Linear $\Theta$-smoothing width (the fit's energy resolution); `<= 0` gives a hard step. |
| `yn_sbe_coulomb` | — | `'n'` | `'y'`: enable the **Coulomb (time-dependent Hartree–Fock / exchange) renormalization** (§8, Golde–Kira–Meier–Koch). Non-k-local mean field, $O(N_k^2)$ per step — off by default, best on modest grids. |
| `sbe_coulomb_epsilon` | — | `-1.0d0` | Background dielectric constant $\varepsilon$ screening the exchange kernel. `≤ 0` ⇒ material default (GaAs `12.9`, Si `11.7`). |
| `sbe_coulomb_strength` | — | `1.0d0` | Overall scaling of the exchange kernel (set `0` to disable while leaving the flag on; `>1` to enhance). |
| `sbe_coulomb_screen_au` | Bohr⁻¹ | `0.0d0` | Yukawa screening $\kappa$ regularizing $V(q)\propto1/(q^2+\kappa^2)$; `0` = bare Coulomb with the $q=0$ self-term excluded. |

*Note: Internal conversions to atomic units (Hartree) are handled automatically (`kB_au`, `au_fs`).*

### Super-compute / dissipation channels (Parts A–G)

These `&sbe` parameters add Silicon support and the optional nonlocal CPTP dissipation channels. **Every flag below defaults OFF (or to a value that leaves the block inert), so a run that sets none of them is byte-for-byte identical to the legacy GaAs dynamics.**

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `sbe_ii_form` | — | `'auto'` | Impact-ionization fit form. `'auto'` ⇒ material default; or set `'stobbe_quartic'` (GaAs, $a=4$, hard threshold) / `'keldysh_quadratic'` (Si, $a=2$, soft threshold). |
| `sbe_ii_exponent` | — | `-1.0d0` | Operative fit exponent $a$ in $\gamma=P(\varepsilon^{\rm kin}-E_{\rm th})^a$. `≤ 0` ⇒ material default (GaAs `4`, Si `2`; set `4.6` for Si full-band). |
| `yn_sbe_superres` | — | `'n'` | `'y'`: nonlocal **super-compute** master switch — ring-pipeline-MPI $\Sigma^{\rm HF}$ + nonlocal-partner impact ionization + dissipator sub-cycling. |
| `yn_sbe_eph` | — | `'n'` | `'y'`: **electron–phonon** population-relaxing Lindblad (full phonon table; toggle Kuhn-Zurek off when using it). |
| `sbe_eph_temperature_k` | K | `300.0d0` | Phonon-bath temperature $T_{\rm ph}$ for the Bose factor $N_B$. |
| `sbe_eph_nu_sat` | s⁻¹ | `-1.0d0` | e-ph saturation rate. `≤ 0` ⇒ **material default** (GaAs `1.0e14`, Si `1.3e14`). |
| `sbe_eph_eps0_ev` | eV | `0.8d0` | Onset $\varepsilon_0$ of the rate saturation $\nu(\varepsilon)=\nu_{\rm sat}[1-e^{-(\varepsilon/\varepsilon_0)^n}]$. |
| `sbe_eph_n` | — | `2` | Saturation shape exponent $n$. |
| `yn_sbe_eeh` | — | `'n'` | `'y'`: **carrier–carrier** (e-e/e-h) CPTP thermalization to a Fermi–Dirac (conserves number + energy). |
| `yn_sbe_auger` | — | `'n'` | `'y'`: **Auger recombination**, number-conserving CPTP. **Two paths:** (a) *with the nonlocal ring* (`yn_sbe_superres='y'` + `yn_sbe_impact_ionization='y'`) it is the **exact time-reverse of the inter-k impact ionization** — same \|M\|²/ε(q)/umklapp weight, reversed Fermi factors — so it needs **no separate $C$** (the rate scale *is* the cited II magnitude, by detailed balance; verified G=R at equilibrium). (b) *k-local fallback* (ring off): the density-gated $\gamma=C n^2$ ($R=C n^3$) Lindblad, inert below $n_{\rm gate}$, which **requires** an explicit `sbe_auger_c_cm6s` (no material ships a verified default). See [`wiki/07_nonlocal_auger.md`](wiki/07_nonlocal_auger.md). |
| `sbe_auger_c_cm6s` | cm⁶/s | `-1.0d0` | Auger coefficient $C$. `≤ 0` ⇒ material default, but **no material has a verified default** (the former CdS `2.0e-30` carried a fabricated citation and was removed), so this must be set explicitly to use the channel. |
| `sbe_auger_n_gate_cm3` | cm⁻³ | `-1.0d0` | Activation density $n_{\rm gate}$. `≤ 0` ⇒ material default (CdS `1e18` [Shah 1986], used only with a user-supplied $C$). |
| `sbe_eeh_nu_sat` | s⁻¹ | `-1.0d0` | Carrier-carrier rate scale. `≤ 0` ⇒ `1.0e14` default (both materials). |
| `yn_sbe_bgr_threshold` | — | `'n'` | `'y'`: density-dependent **band-gap-renormalized** II threshold $E_{\rm th}(t)=E_{\rm th0}-\lvert K\,n^{1/3}\rvert$ (needs impact ionization on). |
| `sbe_bgr_n_gate` | cm⁻³ | `5.0d18` | Apply the BGR shift only above this carrier density. |
| `sbe_bgr_coeff` | eV·cm | `1.9d-8` | BGR coefficient $K$ (tunable in $[1.9,3.8]\times10^{-8}$). |
| `yn_sbe_hf_sublattice_proj` | — | `'y'` | Project $\Sigma^{\rm HF}$ block-diagonally onto the 4 FCC sublattices (folding fix). Inert unless Coulomb is on **and** unfold weights are present. |
| `sbe_search_sigma_e_ev` | eV | grid-matched | Energy-bin width $\sigma_E$ for the final-state partner search (C3). |

### 2026-07 approved refinements (A1–A9 / B3–B4; all default OFF/inert)

The physics of each item is specified in `wiki/00_implementation_status.md`
(the ✅-marked A/B/C entries) and `wiki/08_master_equation.md`; this table is
the input-parameter contract. All of these act **inside the ring channels**
(`yn_sbe_superres='y'` + the corresponding channel flag).

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `sbe_ii_phassist` | — | `0.0d0` | **A1 phonon-assisted II/Auger sidebands.** `> 0`: the energy-match $\delta_\sigma(\Delta E)$ of the ring II *and* Auger kernels gains emission/absorption sidebands at $\Delta E = \pm\hbar\omega_p$ for every mode of the cited per-material phonon table, weighted by $(N_B{+}1)/N_B$ (Bose factors **swapped** in the Auger kernel ⇒ detailed balance holds exactly per sideband). The value is the overall sideband strength (K15-style knob, `1.0` = the natural scale); `0` = off. Needs `yn_sbe_eph` material constants (the phonon table). |
| `yn_sbe_ii_holes` | — | `'n'` | **A2 hole-initiated impact ionization (hhe) + its Auger reverse.** Mirrors the electron quadruple with the valence-side stencil; rate scale = `sbe_ii_prefactor` × the **cited** $C_p/C_n$ ratio from the registry (Si `0.354` [L90/Dziewior-Schmid], GaAs `4.8` [S14]; CdS has none ⇒ forbidden). Hole kinetic energy is measured from the instantaneous Houston VBM. |
| `yn_sbe_ii_fk_soften` | — | `'n'` | **A5 Franz–Keldysh-softened II threshold.** Replaces the hard $(\varepsilon-E_{\rm th})_+$ with the softplus $s = \hbar\theta\,\ln(1+e^{(\varepsilon-E_{\rm th})/\hbar\theta})$, where $\hbar\theta(t) = (F(t)^2/2\mu)^{1/3}$ from the **instantaneous field** — the Quade–Schöll–Rossi "no fixed threshold at MV/cm" physics. → hard threshold as $F\to 0$. |
| `sbe_ii_fk_mu` | mₑ | `0.06d0` | Reduced pair mass $\mu$ in the electro-optic energy $\hbar\theta=(F^2/2\mu)^{1/3}$ (A5). |
| `sbe_ring_vq_floor` | — | `0.0d0` | **B3 relative floor** on the screened Coulomb weight $\lvert V(q)\rvert^2$ in the ring II/Auger: quadruples with $v_q <$ `floor`·max$(v_q)$ are skipped (speed knob). `0` = off = bit-identical to the full sum; `1` kills the channel. |
| `yn_sbe_eph_acoustic` | — | `'n'` | **A4 quasi-elastic acoustic deformation-potential cooling.** Appends a grid-resolved acoustic mode ($\hbar\omega_{\rm ac} = c_s\,q_{\rm grid}$, $D_{\rm ac}=\Xi_d\,q$) to the cited phonon table — continues the cooling below the optical $\hbar\omega$ (removes the freeze-out). Constants per material (all cited): Si $\Xi_d=9.0$ eV, GaAs $7.0$ eV [Fischetti-Laux], CdS $E_1=14.5$ eV [Rode 1970], graphene $D=16$ eV [Hwang–Das Sarma]. **Always applied with the Thomas–Fermi screen** $S(q)=[q/(q+q_{\rm TF})]^2$ built from the instantaneous carrier density (mandatory: the bare CdS $E_1$ would be unphysical at $n\gtrsim10^{18}$ cm⁻³; $S\to1$ at low $n$). |
| `sbe_eph_ac_xi_ev` | eV | `-1.0d0` | Override of the acoustic deformation potential (A4). `≤ 0` ⇒ the registry default above. The graphene $D$ is **substrate-dependent** (H&DS 16 eV is the upper literature scale) — set it here per substrate. |
| `sbe_checkpoint_step` | steps | `0` | **B4 crash-safe checkpoints:** every N steps each rank streams its ρ slab (+ Houston branch, step index, energy, channel ledger) to a per-rank `.bin` file. `0` = off. |
| `yn_sbe_checkpoint_restart` | — | `'n'` | `'y'`: resume from the newest checkpoint set (same-nproc restart; the field is recomputed deterministically). |

**CPTP limiter (automatic, not a flag).** The per-step population transfer of
every ring channel is applied through a **global CPTP limiter**: if the summed
`dpop` of a step would overdraw any state (below 0 or above `occ`), the whole
field is scaled by the largest $s\in[0,1]$ that keeps all states in bounds —
the trace stays *exactly* conserved (a per-state clip would silently create
particles). When the limiter engages it prints
`# ring CPTP limiter engaged: dpop scaled by s = ...` **once** — that message
means the time step is too large for the instantaneous collision flux
(reduce `dt`, or accept the rate-limited dynamics).

The EPM material is chosen in `&epm` (`epm_material = 'GaAs'` | `'Si'` | `'Si_cb'` | `'CdS'` | `'graphene'`); see [EPM ground-state solver](#epm-ground-state-solver-epm). `'CdS'` is wurtzite — its cell is the orthorhombic vector `al(1:3) = (a, a√3, c)`, not a single cubic constant — with three cited, enabled dissipation channels (e-ph, Coulomb, impact ionization); carrier-carrier (e-e) and Auger are **forbidden** for CdS (no cited e-e rate; and no verified CdS Auger coefficient — the previously-cited "Haury 1998" value was a fabricated reference and was removed — see the matrix above and `wiki/02`). **The `&epm` block must also appear in the SBE-step input** whenever a channel auto-selects per-material constants (electron-phonon, or impact ionization / Coulomb left at their `'auto'`/sentinel defaults): the solver reads `epm_material` from the registry and **`epm_material` itself defaults to `'GaAs'` if the block is absent**, so a Si run without it would silently use GaAs constants.

### Baseline (clean) run and how to disable each effect

**Clean run (pure unitary SBE — no dissipation, no decoherence).** All channel flags default OFF and the decoherence rate is disabled by its `-1` sentinels, so a `&sbe` block that sets *none* of the channel parameters runs the bare CF4 / Suzuki–Yoshida / Strang propagation — the reference for isolating the effect of any single block:

```fortran
&sbe
  ! nothing here -> clean unitary baseline:
  !   decoherence OFF (sbe_decoh_* = -1), impact ionization OFF,
  !   Coulomb OFF, e-ph OFF, carrier-carrier OFF, BGR OFF, super-res OFF.
  yn_vnl_correction = 'n'
/
```

Enable one block at a time and compare against this baseline. The table gives the **enable** switch, how to **fully disable** it (the clean-baseline value), and the per-material settings (values that are auto-selected on the `epm_material` switch are marked *auto*):

| Effect | Enable | Fully disable (baseline) | GaAs | Si | CdS |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Kuhn-Zurek decoherence** | `sbe_decoh_temperature_k > 0` **and** `sbe_decoh_tau_m_fs > 0` | either `≤ 0` (default `-1.0` ⇒ $D=0$, exactly the clean run) | $T$, $\tau_m$ user-set | same | ✅ (material-independent) |
| **Impact ionization** | `yn_sbe_impact_ionization = 'y'` | `'n'` (default) | *auto:* `stobbe_quartic`, $a=4$, $E_{\rm th}=2.1$ | *auto:* `keldysh_quadratic`, $a=2$ ($4.6$ full-band), $E_{\rm th}=1.1$ | *auto:* `keldysh_quadratic`, $E_{\rm th}=3.6$; **must set `sbe_ii_prefactor`** (fit param) |
| **Electron–phonon** | `yn_sbe_eph = 'y'` (+ `&epm epm_material`) | `'n'` (default) | *auto:* Fröhlich-LO 36 meV + 5 intervalley, `ν_sat=1.0e14` | *auto:* 6 intervalley g/f, `ν_sat=1.3e14` | *auto:* Fröhlich-LO 38 meV, `ν_sat=2.9e13` (=α·ω_LO) |
| **Carrier–carrier (e-e/e-h)** | `yn_sbe_eeh = 'y'` | `'n'` (default) | `ν_cc=1e14` (default), screening *auto* $\varepsilon=12.9$ | `ν_cc=1e14` (default), screening *auto* $\varepsilon=11.7$ | ⛔ **forbidden** (no cited rate) |
| **Auger recombination** | `yn_sbe_auger = 'y'` (needs explicit `sbe_auger_c_cm6s`) | `'n'` (default) | ⛔ (no verified $C$ yet) | ⛔ (no verified $C$ yet) | ⛔ (no verified $C$; the "Haury 1998" value was fabricated, removed) |
| **Coulomb HF** | `yn_sbe_coulomb = 'y'` | `'n'` (default) **or** `sbe_coulomb_strength=0` | *auto:* $\varepsilon=12.9$ | *auto:* $\varepsilon=11.7$ | *auto:* $\varepsilon=8.9$ |
| **BGR-shifted II threshold** | `yn_sbe_bgr_threshold = 'y'` (needs II on) | `'n'` (default) | `sbe_bgr_n_gate=5e18`, `sbe_bgr_coeff=1.9e-8` | same defaults (tune $K$) | ⛔ (needs II) |
| **HF sublattice projection** | `yn_sbe_hf_sublattice_proj = 'y'` (default; needs Coulomb + unfold) | `'n'` | applies to folded cubic cell | applies to folded cubic cell | n/a |
| **Nonlocal super-compute** | `yn_sbe_superres = 'y'` | `'n'` (default) | ring $\Sigma^{\rm HF}$ + nonlocal-partner II | same | (depends on the channels above) |

**`graphene` (gapless Dirac, not a column above):** its only cited dissipation channel is **electron–phonon** — the two Kohn-anomaly optical modes E2g(Γ, 196 meV) and A1′(K, 160 meV), EPC ⟨g²⟩ from [Piscanec PRL 93, 185503 (2004)] (×2 GW, [Lazzeri PRB 78, 081406 (2008)]). Both modes relax carriers to a *different* k on the cone, so graphene e-ph is physical **only with the ring** (`yn_sbe_eph='y'` **and** `yn_sbe_superres='y'`). The 3D gap-based / many-body channels are **forbidden** on the gapless cone (the k-local impact ionization, the C·n³ Auger, carrier–carrier, Coulomb HF all `error stop`) — graphene carrier multiplication is the nearly-thresholdless **2D Rana process**, which has its own branch: the collinear-collapsed CCCV/CVCC integral (`rana_rcccv`/`dirac_mu_2d`/`rana_qtf`, [Rana PRB 76, 155431 (2007)]) is **implemented and validated** against the cited lifetimes (τ_r≈1.1 ps at n=10¹² cm⁻², 300 K; `tests/test_rana_2d.f90`). *Wiring this into a live CPTP SBE channel is the one remaining Auger TODO* — until then graphene `auger_ok` stays `.false.`. Kuhn–Zurek decoherence (`sbe_decoh_*`) is also **forbidden** for graphene (gapless coherence loss is many-body), aborting with a clear message.

**Provenance rule (strict):** a channel is enabled for a material **only** if its constants are backed by a cited source for that material — *no source ⇒ forbidden*, and the run **aborts** (`error stop`) rather than borrow another material's numbers. The per-material constants come from a single **material registry** (`get_material_params` in `src/ssbe/sbe_superres_ssbe.f90`), which carries per-channel provenance gates `ii_ok / eph_ok / eeh_ok / coulomb_ok`. **CdS** (wurtzite): the band structure is validated against [BC1967](#references--theoretical-background), and three dissipation channels use cited CdS-specific constants — Fröhlich polar-optical e-ph (the primary room-T channel), Coulomb (ε=8.9 [Berlincourt 1963]), and impact ionization (E_th=3.6 eV; the prefactor is a *fit parameter* with no cited value, so the run aborts unless you set `sbe_ii_prefactor`). **Carrier-carrier (e-e) is forbidden** (`eeh_ok=.false.`): there is no cited CdS e-e *rate*, so the FD-thermalization channel would have to borrow the GaAs/Si 1e14 scale. The CdS literature only fixes an e-e *timescale* (sub-100 fs at n ≳ 1e18 [Shah 1986; Elsaesser 1991]), not a rate. **There is no verified CdS Auger coefficient**: the value formerly carried as "C=2.0e-30 cm⁶/s [Haury 1998]" was a fabricated reference (the real Haury et al. paper is [PRL 79, 511 (1997)](https://doi.org/10.1103/PhysRevLett.79.511) on CdMnTe ferromagnetism, unrelated) and was removed, so CdS Auger is gated off. The piezoelectric and deformation-acoustic channels are cited but not yet implemented. **Adding a material is one cited `case` block.** The full effect-support matrix with sources is in [`wiki/02_constants.md`](wiki/02_constants.md).

### Real-time output frequency (`&analysis`)

Real-time SBE propagation writes three diagnostic files (`SYSNAME_sbe_rt_energy.data`, `SYSNAME_sbe_nex.data`, `SYSNAME_sbe_nex_k.data`), each on its own cadence selectable in the `&analysis` namelist. The k-resolved file in particular can grow to gigabytes for dense k-grids/long runs, so its default stride is ten times coarser than the band-projection output:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `out_rt_energy_step` | `10` | Stride (in time steps) for `SYSNAME_sbe_rt_energy.data` (total energy / trace) and the stdout progress log. |
| `out_projection_step` | `100` | Stride for `SYSNAME_sbe_nex.data` (number of excited electrons/holes, summed over k). |
| `out_projection_k_step` | `1000` | Stride for `SYSNAME_sbe_nex_k.data` (Houston-basis population of the lowest conduction band, resolved per k-point). Defaults to 10× `out_projection_step` to avoid producing terabyte-scale output on dense k-grids; increase the stride (larger value) further for very large `nk`/`nt`. |

`SYSNAME_sbe_nex_k.data` reports, for every saved time `t`, one block of `nk` lines `ik, kx, ky, kz, population_lcb`, where `population_lcb = (W^\dagger \rho W)_{aa}` is the diagonal element of the density matrix rotated into the instantaneous Houston (adiabatic) eigenbasis $W$ of $H_{VG}(t)$ for the lowest conduction band $a = N_{elec}/2+1$ — i.e. the same gauge-independent basis used internally by the CPTP dephasing step. With `yn_sbe_spinor='y'` the lowest conduction band is $a = N_{elec}+1$ (the lower spin sub-band of the first conduction level).

### EPM ground-state solver (`&epm`)

The `&epm` namelist configures the local-EPM ground-state solver (`theory='epm'`):

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `epm_material` | — | `'GaAs'` | `'GaAs'` \| `'Si'` \| `'Si_cb'` \| `'CdS'` \| `'graphene'` — selects the cited local form factors (see the materials table below). |
| `epm_cell` | — | `'primitive'` | `'primitive'`: non-orthogonal primitive cell, NO folding (GaAs/Si FCC 2-atom; graphene & CdS hexagonal) — verified interchangeable with the Python primitive refs. `'folded'`: legacy cubic 8-atom supercell + FCC-in-cubic folding (GaAs/Si/Si_cb only; needed for the SBE unfold pipeline). For graphene/CdS only `'primitive'` is valid. |
| `epm_lattice_constant_au` | Bohr | `10.68d0` | Lattice constant $a$ [Bohr] (zincblende/diamond; the FCC primitive cell is built from it). Graphene fixes its own cited $a=2.46$ Å and ignores this. |
| `epm_pw_cutoff_ry` | Ry | `11.1d0` | Plane-wave cutoff. For the zincblende cells it bounds $|G|^2$ in $(2\pi/a)^2$ units (integer shells, same as the cubic ref). For the general non-orthogonal path (graphene) it is a **kinetic** cutoff $|G|^2\,[\text{a.u.}]\le$ value (e.g. `2.94` = 40 eV → the 7-PW Ramanujam basis). |

#### Supported materials & EPM references (Python-primary)

**The Python EPM references are the source of truth — each is validated against its cited benchmark.** The fast in-SALMON MPI EPM (`src/epm/`) is *secondary* but is **verified interchangeable** with the Python references. The cell is chosen by `epm_cell` (`'primitive'`, the default, or `'folded'`):
* **`epm_cell='primitive'` (default) — non-orthogonal primitive cell, NO folding.** GaAs/Si/Si_cb use the FCC 2-atom cell (the Cohen-Bergstresser Hamiltonian on the single-parity plane-wave basis); graphene uses the 2-atom hexagonal honeycomb (Ramanujam π-EPM); CdS uses the 4-atom hexagonal wurtzite cell (polar, V^A≠0, BC1967 form factors interpolated over the cited shells) — graphene and CdS via a general multi-atom structure factor `S(dG)=Σₐ Pₐ e^{-i dG·τₐ}`. **Verified interchangeable with the Python primitive references** (`epm_{gaas,si,graphene,cds}_primitive.py`): identical k-points, band energies to **5×10⁻¹¹ Ha** (GaAs/Si/CdS) / **~6×10⁻⁸ eV** (graphene), identical occupations.
* **`epm_cell='folded'` — legacy cubic 8-atom supercell + FCC-in-cubic parity band-folding** (GaAs/Si/Si_cb only; feeds the SBE folding/unfold pipeline). Same `(2π/a)²`-unit cutoff and reduced-coordinate `k.data`; **verified byte-equivalent to the folded Python reference** (`epm_gaas_reference.py`).

The validated Python references:

| Material | Module | Structure / folding | Validation (cited benchmark) | EPM→SBE pipeline |
| :--- | :--- | :--- | :--- | :--- |
| GaAs | `epm_gaas_reference.py` | zincblende (V^A≠0), 4-fold FCC fold | Cohen-Bergstresser PR 141, 789 (1966) | ✅ full (folding + GS files) |
| **Si** (Kunikiyo, default) | `epm_si_reference.py` | diamond (V^A≡0), 4-fold FCC fold | indirect gap **1.059 eV** (Kunikiyo calc 1.068), CBM @ 0.850·2π/a [Kunikiyo JAP 75, 297 (1994)] | ✅ full |
| **Si_cb** (Cohen-Bergstresser) | `epm_si_reference.py --variant Si_cb` | diamond, same machinery as Si | indirect gap **0.818 eV**, CBM @ 0.850 [CB PR 141, 789 (1966)] | ✅ full |
| CdS | `epm_wurtzite_cds.py` | wurtzite (polar), orthorhombic `al(1:3)`, **2-fold fold (verified exact)** | direct gap **2.55 eV** vs Bergstresser-Cohen PR 164, 1069 (1967) 2.58 eV | ✅ Python emits GS files + **band path** + **2-coset unfold map**; ran end-to-end in the SBE (trace conserved) |
| graphene | `epm_graphene.py` | honeycomb (D6h, V_A=0), **rectangular 4-atom cell + 2-fold fold (verified exact)** | **zero gap** at Dirac K, **v_F=9.6×10⁵ m/s**, Γ=−7.78 eV, M=−2.70 eV [Ramanujam thesis, ASU 2015] | ✅ Python emits GS files + **band path** + **2-coset unfold map**; ran end-to-end in the SBE (trace conserved; unfolded carriers localize at the Dirac points) |

> **Pipeline status (honest):** **all four materials emit SBE ground-state files** (`SYSNAME_k/_eigen/_tm.data`) and run in the SBE, in **two interchangeable representations**:
>
> * **PRIMITIVE (the default, `epm_cell='primitive'`) — recommended.** The non-orthogonal primitive cell, NO folding. The **Fortran `theory='epm'`** path now does **all four** materials (GaAs/Si FCC 2-atom; graphene & CdS hexagonal) — each **verified interchangeable** with the Python primitive references (`epm_{gaas,si,graphene,cds}_primitive.py`): band energies to 5×10⁻¹¹ Ha (GaAs/Si/CdS) / ~6×10⁻⁸ eV (graphene). So the whole EPM→SBE primitive pipeline runs in-SALMON for every material (`samples/exercise_x7_primitive_cell_epm/*_epm_gs.inp`), including the **spinor GaAs** path (`yn_spinorbit='y'` — verified interchangeable with the Python `INCLUDE_SPIN_ORBIT` reference to 5×10⁻¹¹ Ha).
> * **FOLDED (`epm_cell='folded'`) — legacy, for the band-unfold map.** The cubic/orthorhombic supercell + an `_unfold.data` coset map. GaAs/Si fold in-SALMON (Fortran, 4-coset FCC) or via the Python ref; **CdS and graphene use the Python references** (`epm_wurtzite_cds.py` / `epm_graphene.py`) on their folded supercells (CdS orthorhombic 2-fold; graphene rectangular 4-atom 2-fold — both verified exact). The unfold infrastructure (`gs_info_ssbe`, `bloch_solver`, `datafile`, the plotter) is generalized from the hardcoded **4 FCC cosets to N** (4 = cubic, 2 = wurtzite/rectangular); legacy 4-coset GaAs maps still read.
>
> See `wiki/00_implementation_status.md` for the full TODO.

> **graphene is a minimal π-model:** the Ramanujam 3-form-factor local EPM represents the π/π* frontier (1 π electron per carbon, Dirac cone = the lowest band pair), not the full 4-electron valence. The emitted dataset uses `nelec = 4` (4-atom cell), `nstate = 8`, occupation 2 per filled π band, Fermi level at the Dirac point — a clean low-energy model for THz/Dirac physics, not a full-band one.

**`Si` vs `Si_cb`:** identical machinery (diamond, V^A≡0, a=10.26 Bohr, 4-fold folding) — the *only* difference is the V^S form-factor triplet. `Si` uses Kunikiyo (−0.2258, +0.05698, +0.070709 Ry → gap 1.059 eV); `Si_cb` uses Cohen-Bergstresser (−0.21, +0.04, +0.08 Ry → gap 0.818 eV). `Si` (Kunikiyo) is the default as it matches the modern Si gap target; `Si_cb` is for cross-validation. Run `python3 tests/run_all.py` to validate all references.

## Examples

### Minimal SBE Input Example

```fortran
&calculation
  theory = 'sbe'
/

&sbe
  ! ... standard SALMON SBE system parameters ...

  ! ---------------------------------------------------------
  ! 1. Kuhn-Zurek/Caldeira-Leggett Decoherence (strictly CPTP)
  ! ---------------------------------------------------------
  ! lambda = kB*T / tau_m;  enabled only when both are > 0
  sbe_decoh_temperature_k = 300.0d0
  sbe_decoh_tau_m_fs      = 10.0d0

  ! ---------------------------------------------------------
  ! 2. Frozen Core / Active Subspace Optimization
  ! ---------------------------------------------------------
  ! Only evolve bands within ±15 eV of the Fermi level non-linearly.
  ! Deep core bands will only undergo exact linear phase oscillation.
  frozen_core_threshold_ev = -15.0d0
  frozen_free_threshold_ev =  15.0d0
/
```

**Reverting to default behavior:**
* Set `sbe_decoh_temperature_k` and/or `sbe_decoh_tau_m_fs` to a non-positive value to recover the original purely-coherent (no dephasing, $D\equiv 0$, trivially CPTP) behavior.
* Set both `frozen_core_threshold_ev` and `frozen_free_threshold_ev` to `0.0d0` to force all bands into the active nonlinear subspace.

### Minimal EPM → SBE Pipeline Example

#### Standalone Python reference (`epm_gaas_reference.py`)

For quick debugging without building/running SALMON, the repository root also contains `epm_gaas_reference.py` -- a monolithic, single-machine NumPy/SciPy reimplementation of the GaAs Cohen-Bergstresser local-EPM solver (no MPI/OpenMP). It builds the same lattice/plane-wave basis/Hamiltonian/momentum matrices as `src/epm`, and writes byte-compatible `SYSNAME_k.data`/`_eigen.data`/`_tm.data` files that `gs_info_ssbe` can read directly -- so its output can be diffed against the Fortran `theory='epm'` run, or fed straight into an SBE real-time calculation. All parameters (lattice constant, plane-wave cutoff, k-grid, number of bands/electrons, sysname) are hardcoded constants at the top of the script -- including the spinor switch `INCLUDE_SPIN_ORBIT` (see the spinor pipeline example below) -- edit them there and run:

```sh
python3 epm_gaas_reference.py
```

This is a debugging aid only -- `theory='epm'` in SALMON remains the primary, MPI/OpenMP-parallel ground-state path.

```fortran
! Step 1: ground state via local EPM (writes GaAs_k/_eigen/_tm.data)
&calculation
  theory = 'epm'
/
&epm
  epm_material            = 'GaAs'
  epm_lattice_constant_au = 10.68d0
  epm_pw_cutoff_ry        = 11.1d0
/
```
```fortran
! Step 2: real-time SBE propagation reading the files generated above
&calculation
  theory = 'sbe'
/
&system
  ! sysname, lattice vectors, num_kgrid, nstate, nelec must match the EPM run
/
```

### Spinor (spin-orbit) EPM → SBE Pipeline Example

Step 1 — generate the spin-orbit split ground state with the Python reference (the spinor switch is a hardcoded constant at the top of the script):

```sh
# epm_gaas_reference.py:  INCLUDE_SPIN_ORBIT = True   (default)
python3 epm_gaas_reference.py
# writes GaAs_cubic_so_k.data / _eigen.data / _tm.data:
#   64 spin-orbit split bands (occupation 1 per band),
#   mu auto-calibrated at Gamma to Delta0 = 0.341 eV (Gamma8-Gamma7),
#   v_SO = grad_k H_SO written analytically into block 2 (rvnl_tm)
```

Step 2 — real-time SBE propagation on the spinor dataset (note `nstate` doubled, `yn_sbe_spinor` and `yn_vnl_correction` both `'y'`):

```fortran
&calculation
  theory = 'sbe'
/
&control
  sysname = 'GaAs_cubic_so'
/
&units
  unit_system = 'au'
/
&system
  yn_periodic = 'y'
  al(1:3) = 10.68d0, 10.68d0, 10.68d0   ! must match the EPM run
  nelec  = 32
  nstate = 64                            ! 2*Nb spinor bands
/
&kgrid
  num_kgrid(1:3) = 4, 4, 4               ! must match the EPM run
/
&tgrid
  dt = 0.05d0
  nt = 20000
/
&emfield
  ae_shape1 = "Acos2"
  epdir_re1(1:3) = 0.0d0, 0.0d0, 1.0d0
  I_wcm2_1 = 1.0d+11
  tw1 = 500.0d0
  omega1 = 0.056d0
/
&sbe
  yn_sbe_spinor     = 'y'   ! spinor input: occupation 1/band, nelec valence bands
  yn_vnl_correction = 'y'   ! use pi = p + v_SO from rvnl_tm everywhere
/
```

Step 3 — plot (spin pairs are summed into levels automatically):

```sh
python3 plot_sbe_results.py
```

Setting `INCLUDE_SPIN_ORBIT = False` in the script restores the scalar pipeline (`GaAs_cubic`, 32 bands, occupation 2 per band) byte-for-byte; the SBE input then keeps `yn_sbe_spinor = 'n'` (default).

### Band-structure calculation (`theory='dft_band'`)

`theory='dft_band'` diagonalizes the **converged** Kohn-Sham Hamiltonian at k-points along a high-symmetry path and writes the eigenvalues to `band.dat` **and** to `SYSNAME_bandpath.data` (the same plotter/SBE-spectral contract the EPM emits). It is a post-processing step: run a normal `theory='dft'` ground state first, then restart from it. A ready-to-run pair lives in `samples/exercise_04_bulkSi_gs/` (`Si_gs.inp` + `Si_band.inp`).

> **Feeding real DFT levels straight into the SBE.** Beyond band structure, a `theory='dft'` run with `yn_out_tm='y'` writes `SYSNAME_k/_eigen/_tm.data` in the exact `gs_info_ssbe` format (reduced k with the `# b1/b2/b3` header for non-orthogonal cells; `esp[eV]` auto-converted to Hartree; both `_tm.data` blocks, the nonlocal one genuinely non-zero for a real pseudopotential) — so a real (even rough) DFT band structure can replace the EPM as the SBE ground state. The SBE's active/frozen-core window (`frozen_core_threshold_ev`) freezes the deep DFT bonding bands the EPM never had. End-to-end example: [`samples/exercise_x9_bulkSi_dft_sbe/`](samples/exercise_x9_bulkSi_dft_sbe/) (Si FCC primitive, explicit `al_vec`, rough DFT GS → `dft_band` path → short SBE, electrons conserved). **Gotcha:** point `dft_band`'s `base_directory` at a sub-directory (e.g. `./band`) so its path-k `_k.data`/`_eigen.data` don't overwrite the MP-grid GS files the SBE reads.

```sh
cd samples/exercise_04_bulkSi_gs

# 1. Ground state (writes the restart directory data_for_restart/)
salmon < Si_gs.inp

# 2. dft_band restarts from ./restart — point it at the GS output
ln -s data_for_restart restart

# 3. Band structure along L-G-X-M-G  ->  band.dat
salmon < Si_band.inp
```

The path is given explicitly in the `&band` namelist (reduced reciprocal coordinates):

```fortran
&calculation
  theory = 'dft_band'
/
&control
  sysname    = 'Si'
  yn_restart = 'y'      ! restart from the ground-state density in ./restart
/
&band
  lattice         = 'non'              ! use the explicit kpt/ndiv_segment path below
  nref_band       = 20                 ! converge eigenvalues up to this band index
  tol_esp_diff    = 1.0d-5             ! per-band convergence tolerance on |dE| (a.u.)
  num_of_segments = 4                  ! L-G-X-M-G : 4 segments, 5 end points
  ndiv_segment(1:4) = 16, 16, 16, 16   ! k-points per segment
  kpt(1:3,1) = 0.5d0, 0.5d0, 0.5d0     ! L
  kpt(1:3,2) = 0.0d0, 0.0d0, 0.0d0     ! G
  kpt(1:3,3) = 0.5d0, 0.0d0, 0.0d0     ! X
  kpt(1:3,4) = 0.5d0, 0.5d0, 0.0d0     ! M
  kpt(1:3,5) = 0.0d0, 0.0d0, 0.0d0     ! G
  kpt_label(1) = 'L'
  kpt_label(2) = 'G'
  kpt_label(3) = 'X'
  kpt_label(4) = 'M'
  kpt_label(5) = 'G'
/
```

`band.dat` starts with a small header (`Number_of_Bands`, `Number_of_kpt_in_each_block`, `Number_of_blocks`), then one `ik  k_red(1:3)  k_cart(1:3)` line per k-point, followed by `ik  ib  energy(spin...)` eigenvalue lines (energies in Hartree). For the sample above the silicon valence-band top sits at $\Gamma$ with the conduction-band minimum near $X$ (indirect gap), as expected for an LDA silicon band structure.

`plot_sbe_results.py` plots `band.dat` directly (it is picked up automatically alongside the other band-structure files): energies are converted to eV and shifted to a valence-band-maximum reference, with vertical guides drawn at the detected path nodes (direction changes). Since `band.dat` carries no occupations, the VBM band index defaults to `nb//2` (half filling); override it with `--band-vbm IDX`.

```sh
cp plot_sbe_results.py /path/to/band_calculation/
cd /path/to/band_calculation/
python3 plot_sbe_results.py --only-bands --energy-range -13 7   # -> sbe_plots/band_dat_band.png
```

| `&band` parameter | Default | Description |
| :--- | :--- | :--- |
| `lattice` | `''` | `'non'`: take the path from `kpt`/`ndiv_segment` below. `'sc'`/`'fcc'`/`'bcc'`/`'hex'`: use a built-in default path for that Bravais lattice. |
| `nref_band` | `0` | Eigenvalues are converged (and convergence is checked) up to this band index. |
| `tol_esp_diff` | `1.0d-5` | Per-band convergence tolerance on the eigenvalue change between iterations (a.u.). |
| `num_of_segments` | `0` | Number of straight segments in the path (a path of `N` segments has `N+1` end points). |
| `ndiv_segment(:)` | `0` | Number of k-points sampled along each segment. |
| `kpt(1:3,:)` | `0` | Segment end points in **reduced reciprocal** coordinates (one more than `num_of_segments`). |
| `kpt_label(:)` | `''` | Optional labels for the end points (`'G'`, `'X'`, ...). |
