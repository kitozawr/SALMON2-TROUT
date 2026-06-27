# Implementation Status — live progress tracker

**Read this first on every resume.** It records what is done, what is next, key decisions, and the test inventory. Update it in the same commit as any code change.

Roadmap source: the Si + nonlocal-super-compute task (Parts A–F) plus future crystals/effects. The maintainer drives one bounded increment per session ("continue"). Test grid: **4×4×4 (or 5×5×5), scalar (no spinor)**.

---

## Legend
✅ done & tested  · 🚧 in progress · ⬜ not started · 🔭 future / out of current scope

---

## Status by part

| Part | Item | State | Branch / PR | Notes |
|---|---|---|---|---|
| pre | CF4/Yoshida/Strang propagator, Kuhn-Zurek dephasing, frozen core | ✅ | merged | do NOT touch the integrator |
| pre | k-local impact ionization (Stobbe quartic), GaAs | ✅ | merged | Houston basis, HF-factorized amplitude damping |
| pre | Coulomb TDHF exchange Σ^HF (Golde-Kira-Meier-Koch) | ✅ | merged (#41/#43) | non-k-local, all-gather/step; δρ subtraction; Hermitian |
| pre | EPM GaAs + 4-fold FCC folding + unfold pipeline | ✅ | merged | `epm_gaas_reference.py`, `SYSNAME_unfold.data` |
| **A** | Silicon EPM (`epm_material='Si'`, V^A=0, Kunikiyo) | ✅ | #44 | gap 1.059 eV conv (Kunikiyo calc 1.068), CBM 0.86·2π/a |
| **B** | II fit-form switch (`sbe_ii_form`, `sbe_ii_exponent`) | ✅ | #44 | d**a, prefactor au_ev**a; GaAs a=4 unchanged |
| **E** | HF sublattice-block projection (`yn_sbe_hf_sublattice_proj`) | ✅ | #44 | proj_ij=Σ_s w_s(i)w_s(j) off-diag, diag kept; Hermitian; default 'y' |
| **C** | Nonlocal super-compute (`yn_sbe_superres`) | ✅ | #44 | C0-C8 all done; ring (D) + nl-II + nl-eph + e-e all on the super-mode path |
| C0 | Flags scaffolding (10 params, default OFF) | ✅ | #44 | namelist read+bcast+log; default run byte-for-byte unchanged |
| C-prim | Pure rate/search primitives module `sbe_superres_ssbe.f90` | ✅ | #44 | ν(ε), N_B, bins, Fröhlich asinh, II rate, BGR, amp-damp map, golden-rule prefactor + unit conv, thermal split, Si/GaAs tables |
| C2 | Houston-basis adiabatic populations reuse | ✅ | #44 | e-ph runs in the existing houston_dissipate ZHEEV basis (t2 = U†ρU) |
| C3 | Energy-bin final-state SEARCH (energy_partner_weights) | ✅ | #44 | windowed broadened-delta weights over candidates; deterministic, no MC; unit-tested |
| C4 | Nonlocal impact ionization (momentum exchange) | ✅ | #44 | valence partner + Pauli factors from BZ-averaged occupation (gathered/step); CPTP; end-to-end stable. Full momentum-resolved final states = refinement |
| C5 | e-ph population-relaxing Lindblad (k-local, full phonon table) | ✅ | #44 | Si 6 intervalley / GaAs LO+5; golden-rule weights D²/ħω (norm.), detailed-balance emis/abs, ν(ε) cap, Pauli-clamped, CPTP; trace=32 conserved. Nonlocal version pending (C4/D) |
| C6 | CPTP gate (amplitude-damping map test) | ✅ | #44 | trace, qubit positivity det≥0, transfer formulas, Hermiticity, γ=0 identity |
| C7 | BGR-gated II threshold | ✅ | #44 | running n(t) shifts E_th=E_th0−|K n^⅓| above gate; end-to-end stable |
| C8 | Dissipator sub-cycling | ✅ | #44 | m_sub from eph_numax·τ; II+e-ph split into m CPTP sub-steps; trace conserved at ν_sat=1e18. m=1 (unchanged) when e-ph off |
| **D** | Ring/pipeline MPI (systolic ring, one fused pass) | ✅ | #44 | Σ^HF via ring in super-mode; O(Nk/P) mem; ring==all-gather (6e-21), MPI 1==2 ranks (per-k 0.0). nl-II/eph/e-e accumulators slot in |
| **G** | Screening primitives (TF/Debye, Lindhard/RPA, LOPC) | ✅ | #44 | pure functions + GaAs/Si dielectric constants; unit-tested |
| **F** | Carrier-carrier (e-e/e-h) CPTP thermalization (`yn_sbe_eeh`) | ✅ | #44 | intra-k FD relaxation: conserves number AND energy exactly, EID coherence damping, CPTP; unit-tested + end-to-end. Inter-k momentum-resolved on ring = refinement |
| **MS** | Maxwell-SBE multiscale adaptation (`theory='multiscale_experiment'`) | ✅ | #44 | all A–G channels usable per-macropoint via the shared `init_sbe_bloch_solver`/`dt_evolve_bloch_cf4`; `sbe%icomm` set unconditionally (was Coulomb-only → MPI_COMM_NULL for BGR/nl-II-only); banner-print gated to one macropoint (`verbose`) |
| **VG** | Band budget: VG basis sufficiency & N_b convergence | ✅ | #44 | separate axis from PW cutoff, NOT cured by Houston basis; primitives + test (interlacing/shift formula) + runtime P_top monitor (warns on error channel, continues) |
| doc | Wiki pages 01–06 committed as long-term memory | ✅ | #44 | maintained per increment |

---

## Decisions log (gotchas that bit us / must not be re-litigated)
- **Si gap:** a 3-parameter local EPM gives ~1.06 eV converged = Kunikiyo's own calc (1.068), NOT the 1.12 eV experimental value. The spec's "30 meV from 1.12" is against experiment; treat **1.068 eV (Kunikiyo calc)** as the real target. CBM at 0.86·2π/a along ⟨100⟩.
- **Si diamond:** V^A ≡ 0 for all shells → the existing zincblende `VS·cos + i·VA·sin` machinery is reused verbatim with V^A=0 (purely real diamond). τ=(a/8)(1,1,1) for both.
- **Coulomb Σ^HF uses δρ = ρ − ρ₀** (ρ₀ = diag(occup)) so Σ(t=0)=0 and the equilibrium EPM gap is unchanged. Verified Hermitian; 1-rank vs 2-rank populations bitwise-identical.
- **nex density** = (excited e per cell)/V_cell; the /N_k is already in `calc_trace` via /Σ(kweight) (BZ average). Do NOT divide by N_k again.
- **Max ionization ceiling** (32 e per cubic cell) = 1.77e23 cm⁻³ (full valence inversion).
- **GitHub CI is serial-only by request** (compiler/syntax gate). MPI verified locally.
- **Commits are signed-off** (`git commit -s`), trailers Co-Authored-By + Claude-Session.
- **Part G = screening** (TF/Debye, static Lindhard/RPA = default, dynamic LOPC = GaAs-only); **Part F = carrier-carrier (e-e/e-h)** collision channel that USES G.
- **No HF double-counting:** carrier-carrier (F) is the correlation (2nd-Born/GW) self-energy — dissipative only. The static screened-exchange energy shift stays SOLELY in HF (Σ^HF). Do not add it twice.
- **Carrier-carrier invariants:** conserves Σf_k (number) AND ΣE_k f_k (energy) within the carrier subsystem (it thermalizes to a hot Fermi-Dirac, does NOT relax energy to the lattice). Use both as machine-precision validation invariants.
- **e-e/e-h rate scale:** 1e13–1e14 s⁻¹ at n=1e17–1e19 cm⁻³ (thermalization ~10–200 fs). Static screening under-estimates; dynamic (LOPC) needed for sub-100-fs.
- **Maxwell-SBE multiscale:** the multiscale driver runs ONE independent SBE cell per macropoint, driven by that point's macroscopic Maxwell A(t). It calls the SAME `init_sbe_bloch_solver`/`dt_evolve_bloch_cf4` as the single-cell solver, so every A–G channel flows through the `&sbe` namelist automatically — no per-channel wiring in `multiscale_ssbe`. The nonlocal collectives (Coulomb all-gather/ring, BGR + nl-II global reductions) reduce over the **per-macropoint** group `icomm_macro`, so momentum exchange is correctly confined to each cell's own k-grid. `sbe%icomm` MUST be set for every run (not only Coulomb) — it is now set unconditionally in init.
- **Provenance rule (STRICT — agreed with maintainer):** a channel may be enabled for a material ONLY if its constants are backed by a cited source for THAT material. **No source ⇒ invalid ⇒ forbidden**, and the SBE init `error stop`s — constants are NEVER transferred from another material. The registry carries per-channel gates `ii_ok/eph_ok/eeh_ok/coulomb_ok`. (This corrected an earlier mistake where CdS had estimated/copied dissipation constants — those were removed.)
- **CdS (wurtzite P6₃mc):** EPM band structure VALIDATED + structure registered; all dissipation FORBIDDEN. Cell = orthorhombic `al(1:3)=(a,a√3,c)` (not a single cubic constant). Form factors are the REAL cited **wurtzite** values from **Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967), Table II** (LOCAL potential — spherical atomic potentials, no nonlocal term ⇒ `rvnl_tm=0`; no cited CdS nonlocal parameter exists, so none is fabricated). **Validated:** direct gap at Γ converges to **2.55 eV vs BC1967's 2.58 eV** (|Δ|≈0.03 eV); structure factors match Table II (002/101/102). Two bugs that gave 13 eV were fixed: (1) potential normalized by **total atoms/cell** (1/n, the BC1967 "volume per atom" normalization — not per species); (2) use the **wurtzite** Table II form factors, not the zinc-blende anchors. **Folding to the SBE cell is EXACT:** the orthorhombic supercell (8 atoms, the al-vector cell) block-diagonalizes over the 2 cosets (off-coset |H|≈8e-17, analogue of the cubic 4-fold FCC folding) and reproduces the same 2.54 eV gap; coset 0 = Γ_hex (direct gap), coset 1 = zone-edge (6.2 eV). **CdS dissipation channels use the cited constants from the CdS physics-methods spec (md):** Fröhlich polar-optical e-ph (primary; ħω_LO=38 meV [Raman], ν_sat=α·ω_LO=2.89e13 s⁻¹ from α=0.5 [cyclotron]), Coulomb (ε=9.0 [md]), impact ionization (E_th=1.5·E_g=3.6 eV [md]; **prefactor is a fit parameter** with no cited value → the run aborts unless `sbe_ii_prefactor` is set). Carrier-carrier stays forbidden (no cited CdS rate). Piezoelectric/deformation-acoustic cited but not yet SBE channels. `epm_wurtzite_cds.py`, `tests/test_wurtzite_cds_epm.py`.
- **Material registry (single source of per-material constants):** `get_material_params(name)` in `sbe_superres_ssbe.f90` returns one `s_material_params` struct (dielectric ε0/ε∞, II fit form/exponent/prefactor/threshold, e-ph phonon table + ν_sat, lattice, diamond flag) assembled from the cited module constants. Every channel auto-selects through it; the namelist defaults are sentinels (`sbe_ii_form='auto'`, `sbe_ii_exponent/prefactor/threshold ≤ 0`, `sbe_coulomb_epsilon ≤ 0`) that resolve to the material value, and any explicit namelist value overrides. GaAs resolves to the exact legacy numbers (byte-for-byte). Unknown material + a material-dependent channel ⇒ `error stop` with the supported list. **Adding a material = one `case` block + its name in `MAT_SUPPORTED`** — no edits scattered across channels. (EPM form factors are the other per-material table, already case-based in `epm_cohen_bergstresser.f90`.)
- **VG band budget ≠ PW cutoff:** the band count N_b carried into the dynamics is a correctness axis SEPARATE from the EPM plane-wave cutoff, and the Houston/adiabatic basis does NOT fix an insufficient N_b — it diagonalizes the already-truncated H_VG^(N_b), inheriting the truncation error (Hylleraas-Undheim-MacDonald). Monitor `P_top = max_k ρ̃_{N_b,N_b}` (criterion a): the real-time solver warns on the **error channel** and **continues** when P_top > 1e-3 (it does not stop — the user re-runs with more bands + an N_b convergence study). Re-verify N_b for every new material / driver wavelength; the converged value does not transfer. See wiki/06.
- **a.u. conversion audit (verified):** rates s⁻¹→a.u.⁻¹ via `×(au_fs·1e-15)`; energies eV→Ha via `/au_ev`; `mev_to_ha=E·1e-3/au_ev`; deformation potentials `d_evcm_to_au=D·BOHR_CM/au_ev`, `d_evang_to_au=D·BOHR_ANG/au_ev`; mass density `rho_gcm3_to_au=ρ·BOHR_CM³/ME_G`; carrier density a.u.⁻³→cm⁻³ via `1e24/BOHR_ANG³`. The e-ph `eph_wrel=D²/ħω` is normalized to a dimensionless per-mode weight, so GaAs (eV/Å) vs Si (1e8 eV/cm) deformation-potential units only need within-material consistency (the absolute scale is `ν_sat`).

## CPTP invariants (must hold for every new dissipator)
- Each exp(τD) a genuine GKLS map; clamp all (1−ρ) Pauli factors to [0,1].
- Yoshida wraps ONLY the unitary part; dissipator always τ>0.
- Hermitian-only dissipator must conserve populations exactly (Γ_aa=0).

---

## Test inventory (`tests/`)
| Test | Covers | How |
|---|---|---|
| `test_si_epm_gap.py` | Part A | primitive sublattice-block Si gap ≈1.06 eV, CBM≈0.86·2π/a |
| `test_ii_form_switch.py` | Part B | rate γ=P(ε−E_th)^a scaling for a=2 vs 4 |
| `test_hf_sublattice_proj.py` | Part E | proj zeroes inter-sublattice, keeps diag, Hermitian; real GaAs weights |
| `test_superres_rates.f90` | Part C primitives | ν(ε) limits, N_B, energy-bin area/peak, Fröhlich asinh, II 2^a scaling, BGR −19/−41 meV, data tables (standalone gfortran) |
| `test_eph_cptp.f90` | Part C5/C6 | amplitude-damping map: trace, qubit positivity (det≥0), transfer formulas, coherence damping, Hermiticity, γ=0 identity |
| `test_superres_rates.f90` (extended) | Part C5 | + unit conversions (meV/eV·cm/eV·Å/g·cm⁻³→a.u.), golden_rule_prefactor, eph_thermal_split (fe+fa=1, fe/fa=(N+1)/N) |
| `test_screening.f90` | Part G | eps_TF limits, Lindhard F(0)=1/F(1)=½/monotone, eps_Lindhard→TF at small q, plasmon ω_p², LOPC Vieta sum/product + anticrossing |
| `test_carrier_carrier.f90` | Part F | carrier_carrier_relax conserves number+energy, CPTP positivity, EID damping, converges to FD, inversion no-op; fit_fermi_dirac self-consistency |
| `test_vg_basis_nb.f90` | VG basis sufficiency | eta admixture + threshold, conv-error metric, P_top gate; Hylleraas-Undheim-MacDonald interlacing/upper-bound + 2nd-order truncation-shift formula (self-contained Jacobi solver, no LAPACK) |
| `test_material_registry.f90` | material registry + provenance gates | GaAs/Si constants + all gates `.true.`; CdS structure-only + all dissipation gates `.false.` (forbidden); unknown→not-found |
| `test_wurtzite_cds_epm.py` | CdS EPM geometry + cited FFs | orthorhombic cell from al, 8 atoms (4 Cd+4 S), Hermitian H, broken inversion, exact BC1967 Table II anchors; band solve reported NOT-yet-validated (not asserted) |
| _(add per increment)_ | | |

End-to-end smoke (manual): scalar GaAs 4³, `yn_sbe_superres='y'` + `yn_sbe_eph='y'`
→ runs stable, correct diagnostics (ħω=36 meV, N_B=0.33 @300 K), trace = 32.000
conserved over all steps (CPTP at the dynamics level).

End-to-end smoke (manual, not in run_all): scalar GaAs 4³, `yn_sbe_coulomb='y'`
+ `yn_sbe_hf_sublattice_proj='y'`, weak field → runs stable, finite, projection
diagnostic confirmed. Full "zero spurious Γ→X transfer" weak-field validation
is a longer run, deferred to a dedicated validation pass.

Run all: `python3 tests/run_all.py` (each test prints PASS/FAIL and exits nonzero on failure).

---

## Materials — Python EPM references (THE source of truth) ✅
The **Python EPM references are primary** and each is validated against a cited
benchmark (run `python3 tests/run_all.py`, currently 11/11). The fast in-SALMON
**MPI EPM (`src/epm/`) is SECONDARY** — calibrated against the Python refs and
**deprecated for the non-cubic materials (CdS, graphene) until debugged**;
generate those GS with the Python references.

| Material | Module | Band-validated | Folding | SBE GS files (eigen/tm/k + unfold) |
|---|---|---|---|---|
| GaAs | `epm_gaas_reference.py` | CB 1966 | ✅ 4-fold FCC (exact) | ✅ emitted (full EPM→SBE pipeline) |
| Si (Kunikiyo, default) | `epm_si_reference.py` | gap 1.059 eV, CBM 0.850 | ✅ 4-fold FCC (reuses GaAs) | ✅ emitted (full pipeline) |
| Si_cb (Cohen-Bergstresser) | `epm_si_reference.py --variant Si_cb` | gap 0.818 eV, CBM 0.850 | ✅ 4-fold FCC | ✅ emitted |
| CdS | `epm_wurtzite_cds.py` | gap 2.55 vs 2.58 eV (BC1967) | ✅ 2-fold orth←hex, **verified exact** (off-coset \|H\|≈8e-17) | ❌ **NOT wired** (folding verified in the Python module only; no CdS `_unfold.data`/GS emission) |
| graphene | `epm_graphene.py` | zero gap at K, v_F 9.6e5, Γ −7.78, M −2.70 | ❌ **NO folding** — Config-A primitive cell only (the simpler task; G2 orthorhombic 4-atom + 2-fold not done) | ❌ NOT wired |

**HONEST status (do not overclaim):** only GaAs/Si have the **complete EPM→SBE
pipeline** (folding + GS-file emission). CdS has the **folding implemented and
verified exact** but it is NOT yet plumbed into the SBE (no `_unfold.data`, no
eigen/tm/k emission) — it is at the band+folding validation stage. **graphene
did the SIMPLER task**: primitive-cell band validation, **no folding** (PART G2
not implemented), and the prompt's "skip folding if the SBE accepts a
non-orthogonal cell" shortcut was taken WITHOUT verifying the SBE actually
accepts the hexagonal non-orthogonal cell.

**Si vs Si_cb:** identical machinery (diamond, V^A≡0, a=10.26 Bohr, 4-fold fold);
ONLY the V^S triplet differs (Kunikiyo vs CB). See README "Supported materials".

## Next action on resume — STANDING TODOs (context restarted)
Branch `develop-2.0.0` (merged from PR #44). Order of priority per the maintainer:

1. **Python EPM refs for all 4 materials FIRST, MPI EPM second.** Done: all four
   Python refs validated (table above). The MPI Fortran EPM still needs
   debugging/calibration for CdS + graphene (non-cubic folding) — mark deprecated
   in code until it reproduces the Python-ref gaps. *(README documents this.)*

2. **Sublattice-projection-with-Coulomb is the ONLY correct HF mode — verify on
   ALL materials.** `apply_hf_sublattice_projection` zeroes the spurious
   inter-coset Σ^HF created by folding. The cosets DIFFER per material (GaAs/Si
   = 4-fold FCC; CdS/graphene = 2-fold, different coset vectors). TODO: confirm
   `yn_sbe_hf_sublattice_proj='y'` + `yn_sbe_coulomb='y'` uses the correct
   per-material coset/unfold weights (gs%unfold_w) and is block-diagonal to
   machine precision for each. Currently wired for the FCC 4-fold path; the
   CdS/graphene 2-fold unfold maps must feed the same projector.

3. **GRAPHENE dissipators (G4/G5) + no-Kuhn-Zurek (G6)** — registry entry
   `epm_material='graphene'` NOT yet added. Needed (constants in wiki, PROMPT_graphene_full):
   - e-ph: E2g 196 meV (g²=0.0405 eV²), A1' 160 meV (g²=0.0994 ×2 GW), acoustic
     D_ac=16 eV (tunable); total ~1e10 s⁻¹ @300K, 196 meV threshold step.
   - e-e/Auger: α=2.2/ε_eff, static|dynamical RPA switch, gapless CM up to ~2.
   - **assert: graphene + Kuhn-Zurek flag → error** (coherence is many-body only).
   - is_diamond=.true. (V_A=0, centrosymmetric); odd-only HHG (linear), 6m±1 (circular).

4. **AUGER CPTP Lindblad primitive** (wiki Section 13) — write the
   number-conserving Auger/impact-ion jump-operator channel in
   `sbe_superres_ssbe.f90` (reuse amp_damp/carrier_carrier pattern) + unit test;
   wire to graphene (G5) and **reuse for CdS** (C=2.0e-30 cm⁶/s [Haury 1998],
   density-gated at n≥1e18 [Shah 1986]). CdS `eeh_ok` already documented `.true.`
   in wiki — REGISTRY CODE still has CdS eeh_ok=.false.; flip it + add
   CDS_AUGER_C / CDS_EE_ACTIVATION_N constants to match the wiki.

5. **Clean remaining OUTDATED wiki fragments** — several pages still say "Parts
   A–G / 7 tests / PR #44 / forbidden CdS e-e"; sweep wiki/00,04 for stale text
   (the effect-support matrix + section 12/13 in wiki/02 are now current).

6. **GRAPHENE FOLDING (G2) — graphene EPM is the SIMPLER task so far.** Current
   `epm_graphene.py` is Config-A primitive-cell ONLY (Dirac cone validated), NO
   folding. To bring it level with CdS:
   (a) **FIRST check whether the SBE accepts a NON-ORTHOGONAL (hexagonal) cell**
       (G2.0). If YES → use the 2-atom hexagonal primitive cell, skip folding
       entirely (cleaner). If NO → (b).
   (b) implement the orthorhombic 4-atom rectangular cell (a × √3a, zigzag x /
       armchair y), the 2-fold hex→rect fold, the `graphene_unfold.data` map,
       and assert Dirac fold lands at ⅔ Γ–X (k_x=±0.851 Å⁻¹), block-diagonal to
       machine precision (mirror `orth_folding_check` in epm_wurtzite_cds.py).

7. **WIRE CdS + graphene EPM → SBE GS files.** Only GaAs/Si emit the SBE
   ground-state data (eigen/tm/k + `_unfold.data`). CdS folding is verified IN
   THE PYTHON MODULE but NOT plumbed into the SBE: it does not yet write the
   CdS GS files / unfold map. TODO: have epm_wurtzite_cds.py (and graphene) emit
   the eigen/tm/k.data + per-material `_unfold.data` in the SBE's read contract
   (reuse the GaAs writers), so the full EPM→SBE pipeline closes for all 4.

8. Optional refinements (unchanged): inter-k momentum-resolved F/C4 on the ring;
   dynamic LOPC; absolute golden-rule e-ph prefactor; Chefonov Si bleaching run.

Also pending (validation, not code): a longer Si super-mode run to check the
e-ph cooling actually reduces the Drude conductivity (bleaching) -- needs a Si
dataset (epm_material='Si', scalar) and the carrier-density/current diagnostics.
