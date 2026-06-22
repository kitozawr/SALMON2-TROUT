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

## Next action on resume
**ALL roadmap parts A–G are implemented, tested (7/7) and documented.** Every
new channel is CPTP and gated OFF by default; existing GaAs runs unchanged.
Branch `claude/sbe-silicon-superres`, PR #44.

What remains is **validation + refinement** (not new roadmap items):
1. **Physical validation runs** (need a scalar Si dataset, `epm_material='Si'`):
   reproduce the Chefonov Si THz-bleaching staging (wiki/04 Example 4) — e-ph
   cooling reduces the Drude conductivity; then enable II + carrier-carrier.
   Add carrier-density/current diagnostics if needed.
2. **Refinements (documented as such, optional):** the full inter-k
   momentum-resolved versions of carrier-carrier (F) and impact ionization (C4)
   on the ring (current: F = intra-k FD relaxation, C4 = global-partner
   sourcing — both CPTP and conserving, but not yet full final-state-resolved);
   dynamic LOPC screening wired into a carrier-carrier rate; the golden-rule
   deformation-potential PREFACTOR used for absolute e-ph rates (current: ν_sat
   scale + relative weights).
3. **Future crystals/effects** (the maintainer noted "other crystals" beyond F):
   new materials reuse the EPM form-factor table + folding; new effects slot in
   as gated CPTP channels following the same pattern (primitive in
   sbe_superres_ssbe + unit test + channel + end-to-end smoke + wiki + status).

Earlier detailed plan (kept for reference):
- C3 full energy-windowed expanding-radius partner search (enumeration, no MC).
- C4 nonlocal momentum-exchange impact ionization + nonlocal e-ph (genuine q).
- D systolic-ring MPI replacing the all-gather; ONE fused ring pass for Σ^HF +
  nonlocal II + nonlocal e-ph (+ the F e-e hook); active-subspace compression;
  C1 predictor-corrector.
Validate against the Chefonov Si THz-bleaching staging (wiki/04 Example 4)
before claiming physical correctness. Keep all behind yn_sbe_superres OFF.

Also pending (validation, not code): a longer Si super-mode run to check the
e-ph cooling actually reduces the Drude conductivity (bleaching) -- needs a Si
dataset (epm_material='Si', scalar) and the carrier-density/current diagnostics.
