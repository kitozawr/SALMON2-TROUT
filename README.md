# SALMON-TROUT

**T**ime-**R**esolved **O**pen-system **U**ltrafast **T**ransport — a velocity-gauge
semiconductor Bloch equation solver with many-body CPTP Lindblad dissipation for
strong-field THz/MIR physics, on EPM band structures (folded supercells **and**
non-orthogonal primitive cells) or real DFT levels.

A fork of [SALMON](http://salmon-tddft.jp/) (Scalable Ab-initio Light-Matter
simulator for Optics and Nanoscience).

---

## The model in one equation

Everything the solver propagates is one GKLS master equation per k-point,

$$
\partial_t \rho(\mathbf k) = -i\big[H_{VG}(\mathbf k,t) + \Sigma^{\rm HF}[\rho],\rho(\mathbf k)\big]
+ \mathcal D_{\rm KZ} + \mathcal D_{\rm e\text{-}ph} + \mathcal D_{\rm II} + \mathcal D_{\rm Aug} + \mathcal D_{\rm ee},
$$

| term | physics | switch |
|---|---|---|
| $H_{VG} = \varepsilon_n\delta_{nm} + \mathbf A(t)\!\cdot\!\boldsymbol\pi_{nm}$ | velocity-gauge band dynamics; **all** Zener/multiphoton injection lives here | always on |
| $\Sigma^{\rm HF}$ | Coulomb exchange: dynamic gap + Rabi renormalization (Golde–Kira–Meier–Koch) | `yn_sbe_coulomb` |
| $\mathcal D_{\rm KZ}$ | Kuhn–Zurek wave-packet dephasing (exactly CPTP Gaussian-kernel Hadamard map) | `sbe_decoh_*` |
| $\mathcal D_{\rm e\text{-}ph}$ | cited per-material phonon tables; **inter-k through the ring** → true intervalley relaxation | `yn_sbe_eph` |
| $\mathcal D_{\rm II}$ | impact ionization; **momentum-conserving nonlocal** 2-particle events with CDRB ε(q) + umklapp G-sum | `yn_sbe_impact_ionization` |
| $\mathcal D_{\rm Aug}$ | Auger: the exact **time-reverse of the II** (detailed balance, no extra coefficient); graphene = 2D Rana Auger/CM | `yn_sbe_auger` |
| $\mathcal D_{\rm ee}$ | carrier–carrier thermalization to a hot Fermi–Dirac (conserves N **and** E) | `yn_sbe_eeh` |

Every dissipator is a genuine Lindblad generator applied in the instantaneous
**Houston basis** — completely positive and trace-preserving for any step, no
positivity clipping anywhere. **The complete mathematical specification — every
term written out, long formulas included, mapped to flags and routines — is
[`wiki/08_master_equation.md`](wiki/08_master_equation.md).** Read that instead
of a textbook.

## Features

- CF4/Magnus 4th-order exponential propagator (exact unitaries, Suzuki–Yoshida
  composition, Strang-split CPTP dissipation) — [`wiki/03`](wiki/03_numerical_methods.md)
- CPTP Lindblad channels: e-ph, e-e, impact ionization, Auger recombination,
  Kuhn–Zurek dephasing — [`wiki/08`](wiki/08_master_equation.md)
- **Momentum-conserving nonlocal ring channels** (inter-k e-ph, impact
  ionization ↔ Auger detailed-balance pair through one systolic-ring gather)
  with **CDRB model ε(q) + umklapp G-sum**, GaAs dynamic free-carrier
  λ²(n(t)), and the graphene **2D Rana Auger/CM** branch — [`wiki/07`](wiki/07_nonlocal_auger.md)
- Materials: **GaAs** (scalar + spin-orbit spinor), **Si**, **CdS** (wurtzite),
  **graphene** (2D Dirac) — strict per-material provenance gates (no uncited
  constants, ever) — [`wiki/02`](wiki/02_constants.md)
- EPM ground states **fully in-SALMON** (`theory='epm'`) for every material and
  spin mode, primitive **non-orthogonal cells by default** (no folding
  artifacts; the folded supercell + exact unfold pipeline kept for
  cross-checks) — [`wiki/05`](wiki/05_folding_unfolding.md)
- **Real DFT levels into the SBE** (`theory='dft'` + `yn_out_tm='y'` →
  `dft_band` → SBE with the frozen-core window) —
  [`samples/exercise_x09_bulkSi_dft_sbe/`](samples/exercise_x09_bulkSi_dft_sbe/)
- HHG with realistic scattering-limited coherence; gauge-invariant current,
  conductivity σ(ω) and STFT maps, 3D BZ population movies —
  [`wiki/09`](wiki/09_plotting_and_analysis.md)
- **Non-Markovian dissipation** (Boroumand 2025 / Meier–Tannor 1999):
  the SFSB memory-integral ionization reference (`yn_sbe_sfsb`) and
  **collisional-memory dephasing** of the e-ph Lindblad (`yn_sbe_colmem` —
  kernel lines from the cited phonon table; adiabatic coherences damp at the
  calibrated Markovian rate, sub-cycle field-driven dressing is protected) —
  [`wiki/10`](wiki/10_open_quantum_systems_literature.md) §6–8

## Quick start

```sh
# build (serial; same as CI). MPI: -D USE_MPI=ON + mpif90.
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
      -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"

# smallest full pipeline: EPM ground state -> SBE dynamics -> pictures
cd samples/exercise_x07_primitive_cell_epm
../../build/salmon < Si_prim_epm_gs.inp   # writes Si_prim_k/_eigen/_tm.data
../../build/salmon < Si_prim_sbe_rt.inp   # clean velocity-gauge dynamics
python3 ../../plot_sbe_results.py -i . -o plots --snapshots --valleys

# validate everything
python3 tests/run_all.py                  # 20/20
```

Then go straight to **the showcase** — every valid channel of every material at
once, watch the carriers roll into the valleys in 3D:
[`samples/exercise_x11_full_dissipation_showcase/`](samples/exercise_x11_full_dissipation_showcase/).

## Examples (each directory is a self-contained lesson)

| sample | what it teaches |
|---|---|
| [`exercise_x07_primitive_cell_epm/`](samples/exercise_x07_primitive_cell_epm/) | the EPM→SBE pipeline on the primitive cell, all 4 materials, clean dynamics + per-channel notes |
| [`exercise_x08_Si_primitive_hhg_basis/`](samples/exercise_x08_Si_primitive_hhg_basis/) | velocity-gauge **basis sufficiency** + reading the output (band maps, BZ maps, HHG); br3d ↔ multiphoton-injection cross-ref |
| [`exercise_x09_bulkSi_dft_sbe/`](samples/exercise_x09_bulkSi_dft_sbe/) | replacing the EPM with a **real DFT** ground state (rough LDA → SBE, frozen core) |
| [`exercise_x10_Si_full_dissipation/`](samples/exercise_x10_Si_full_dissipation/) | Si with all channels — the bleaching-ready copy of x8 (BGR+KZ variant) |
| [`exercise_x11_full_dissipation_showcase/`](samples/exercise_x11_full_dissipation_showcase/) | **the flagship demo**: all 4 materials × every valid channel (Σ^HF instead of BGR, collisional decoherence instead of KZ, nonlocal ring II), valley-rolling 3D movies |
| [`exercise_x13_GaAs_sfsb_nonmarkovian/`](samples/exercise_x13_GaAs_sfsb_nonmarkovian/) | **non-Markovian dissipation**: the SFSB memory kernel vs the relaxation-time approximation (dephasing ionization ×30, low-T suppression ×25, the Im C flip) |
| [`exercise_04_bulkSi_gs/`](samples/exercise_04_bulkSi_gs/) | plain DFT ground state + `dft_band` band structure |

## Configuration philosophy

Three rules cover the whole namelist surface (the complete parameter reference
with per-material recipes is [`wiki/04`](wiki/04_configuration_examples.md)):

1. **Everything defaults OFF.** An `&sbe` block with no flags is the clean
   unitary baseline — bit-identical to legacy runs. Enable one channel at a
   time and diff against it.
2. **Strict provenance.** A channel enables for a material only if its constants
   carry a citation *for that material* (single registry,
   `get_material_params`); otherwise the run **aborts** rather than borrow
   numbers. Escape hatches are explicit (`sbe_ii_prefactor`,
   `sbe_auger_c_cm6s`, `sbe_eeh_nu_sat`).
3. **No double counting — guarded.** Σ^HF and BGR are mutually exclusive (same
   gap physics; `error stop`); the k-local II/e-ph/Auger gate off when the ring
   upgrades them; Kuhn–Zurek is unphysical for graphene (`error stop`) and
   redundant next to collision channels.

Put the `&epm` block (`epm_material=...`) in the **SBE input too** whenever any
channel auto-selects material constants — it defaults to `'GaAs'` if absent.

## Analysis tools

One plotter — [`plot_sbe_results.py`](plot_sbe_results.py) — renders everything:
energies/nex traces, σ(ω) + STFT conductivity, reduced & Cartesian BZ population
snapshots and k–t maps, `--bz3d`/`--bz3d-voxel` 3D BZ movies with valley
markers, per-frame A(k,E) spectral movies, folded/unfolded/spinor band plots.
Three standalone physics probes quantify the *injection* step on the same EPM
bands: `band_field_coupling.py` (matrix elements), `zener_tunneling_estimate.py`
(Kane/Landau–Zener tunnelling, GaAs), `si_three_photon_isosurfaces.py`
(multiphoton I^N + Franz–Keldysh, Si, 3D isosurfaces + MRI slicer). All
documented in [`wiki/09_plotting_and_analysis.md`](wiki/09_plotting_and_analysis.md).

## Documentation & project wiki

The [`wiki/`](wiki/) is the project's **persistent memory** — the single source
of truth that survives context resets. On resume, read
[`wiki/00`](wiki/00_implementation_status.md) first; **all TODO tracking lives
there** and it is updated in the same commit as any code change.

| page | contents |
|---|---|
| [`00_implementation_status.md`](wiki/00_implementation_status.md) | live tracker, decisions log (do-not-relitigate gotchas), test inventory. **Read first.** |
| [`01_physics_models.md`](wiki/01_physics_models.md) | per-channel modelling assumptions & approximations, cited |
| [`02_constants.md`](wiki/02_constants.md) | every default constant with its primary source; effect-support matrix |
| [`03_numerical_methods.md`](wiki/03_numerical_methods.md) | CF4 Magnus, Yoshida, Strang, Houston basis, ring MPI, CPTP proofs |
| [`04_configuration_examples.md`](wiki/04_configuration_examples.md) | **complete parameter reference** + per-material recipes + pipelines (incl. spinor, dft_band, Maxwell-SBE) |
| [`05_folding_unfolding.md`](wiki/05_folding_unfolding.md) | supercell folding theory, exact-folding proof, N-coset unfold |
| [`06_vg_basis_nb_convergence.md`](wiki/06_vg_basis_nb_convergence.md) | the band-count correctness axis; P_top monitor |
| [`07_nonlocal_auger.md`](wiki/07_nonlocal_auger.md) | nonlocal II ↔ Auger theory, ε(q)/umklapp, 2D Rana branch, source-verified coefficients |
| [`08_master_equation.md`](wiki/08_master_equation.md) | **the full master equation, mathematically complete** — every effect as an explicit term |
| [`09_plotting_and_analysis.md`](wiki/09_plotting_and_analysis.md) | the plotter, unfold pipeline, injection probes |
| [`10_open_quantum_systems_literature.md`](wiki/10_open_quantum_systems_literature.md) | open-quantum-system dissipation review; **[B25] SFSB transcription (§6) + the non-Markovian memory-kernel mode (§7)** |

Conventions (units, Houston basis, CPTP, provenance) are at the top of
[`wiki/08`](wiki/08_master_equation.md) and in the pages above; working
conventions with the maintainer (bounded increments, 4³ scalar test grid,
validation by calculation, signed-off commits) in
[`wiki/00`](wiki/00_implementation_status.md).

## References & theoretical background

1. **Commutator-Free Magnus Integrators:** Blanes, S., & Moan, P. C. *J. Comput. Appl. Math.* 142, 313-330 (2002); Alvermann, A., & Fehske, H. *J. Comput. Phys.* 230, 5930-5956 (2011).
2. **Suzuki-Yoshida Composition:** Yoshida, H. *Phys. Lett. A* 150, 262-268 (1990).
3. **CPTP / Lindblad & RBF-kernel positivity:** Schoenberg, I. J. *Ann. Math.* 39, 811-841 (1938) (Bochner/Schoenberg PSD criterion); Schur product theorem.
4. **Caldeira-Leggett / Kuhn-Zurek Decoherence:** Caldeira, A. O., & Leggett, A. J. *Physica A* 121, 587-616 (1983); Zurek, W. H. *Rev. Mod. Phys.* 75, 715 (2003).
5. **Cohen-Bergstresser Local Pseudopotentials:** Cohen, M. L., & Bergstresser, T. K. *Phys. Rev.* 141, 789 (1966). **Wurtzite CdS:** Bergstresser, T. K., & Cohen, M. L. *Phys. Rev.* 164, 1069 (1967). **Si (Kunikiyo):** Kunikiyo, T. et al. *J. Appl. Phys.* 75, 297 (1994).
6. **Spin-Orbit in EPM:** Weisz, G. *Phys. Rev.* 149, 504 (1966); Bloom, S., & Bergstresser, T. K. *Solid State Commun.* 6, 465 (1968); Chelikowsky, J. R., & Cohen, M. L. *Phys. Rev. B* 14, 556 (1976).
7. **Velocity-Gauge SBE / Houston Basis:** Wismer, M. S., & Yakovlev, V. S. *Phys. Rev. B* 97, 144302 (2018).
8. **Original SALMON SBE:** Sato, S. A. et al. *Phys. Rev. B* 92, 115145 (2015).
9. **Coulomb HF renormalization:** Golde, D., Kira, M., Meier, T., & Koch, S. W. *Phys. Status Solidi B* 248, 863 (2011).
10. **Impact ionization fits:** Stobbe, M., Redmer, R., & Schattke, W. *Phys. Rev. B* 49, 4494 (1994); Keldysh, L. V. *JETP* 21, 1135 (1965).
11. **Nonlocal Auger & impact ionization** (`wiki/07`): Laks, D. B., Neumark, G. F., & Pantelides, S. T. *Phys. Rev. B* 42, 5176 (1990) [L90]; Kioupakis, E. et al. *Phys. Rev. B* 92, 035207 (2015) [K15]; Steiauf, D., Kioupakis, E., & Van de Walle, C. G. *ACS Photonics* 1, 643 (2014) [S14]; Rana, F. *Phys. Rev. B* 76, 155431 (2007) [R07].
12. **Non-Markovian strong-field heat bath** (`wiki/10`): Boroumand, N., Thorpe, A., Bart, G., Parks, A. M., Toutounji, M., Vampa, G., Brabec, T., & Wang, L. *Rep. Prog. Phys.* 88, 070501 (2025) [B25]; Meier, C., & Tannor, D. J. *J. Chem. Phys.* 111, 3365 (1999) (exponential bath decomposition).

## License

SALMON-TROUT is a fork of [SALMON](https://salmon-tddft.jp/) and, like the
upstream project, is distributed under the **Apache License, Version 2.0**
(see [`LICENSE`](LICENSE)). Forking an Apache-2.0 project does not change its
license — the original SALMON code and this fork's additions are all Apache-2.0.

    Original SALMON code:
        Copyright 2017-2026 SALMON developers

    Fork modifications and additions (the src/ssbe SBE dissipation engine,
    the src/epm + epm_*.py EPM tooling, the Python analysis tools, the
    samples, and the wiki):
        Copyright 2026 SALMON-TROUT contributors

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

Attribution and the statement of changes required by Apache-2.0 §4(b)/§4(d)
are recorded in [`NOTICE`](NOTICE). Source files inherited unchanged from
SALMON retain their original copyright headers; bundled third-party
components are listed in [`LICENSE.THIRD-PARTY`](LICENSE.THIRD-PARTY) and
[`ACKNOWLEDGEMENTS`](ACKNOWLEDGEMENTS).
