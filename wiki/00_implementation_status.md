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
| **C** | Nonlocal super-compute (`yn_sbe_superres`) | 🚧 | #44 | flags scaffolded (all OFF); primitives done; integration pending |
| C0 | Flags scaffolding (10 params, default OFF) | ✅ | #44 | namelist read+bcast+log; default run byte-for-byte unchanged |
| C-prim | Pure rate/search primitives module `sbe_superres_ssbe.f90` | ✅ | #44 | ν(ε), N_B, bins, Fröhlich asinh, II rate, BGR, amp-damp map, golden-rule prefactor + unit conv, thermal split, Si/GaAs tables |
| C2 | Houston-basis adiabatic populations reuse | ✅ | #44 | e-ph runs in the existing houston_dissipate ZHEEV basis (t2 = U†ρU) |
| C3 | Energy-bin final-state SEARCH (energy_partner_weights) | ✅ | #44 | windowed broadened-delta weights over candidates; deterministic, no MC; unit-tested |
| C4 | Nonlocal impact ionization (momentum exchange) | ⬜ | — | needs C3 expanding search + ring (D) |
| C5 | e-ph population-relaxing Lindblad (k-local, full phonon table) | ✅ | #44 | Si 6 intervalley / GaAs LO+5; golden-rule weights D²/ħω (norm.), detailed-balance emis/abs, ν(ε) cap, Pauli-clamped, CPTP; trace=32 conserved. Nonlocal version pending (C4/D) |
| C6 | CPTP gate (amplitude-damping map test) | ✅ | #44 | trace, qubit positivity det≥0, transfer formulas, Hermiticity, γ=0 identity |
| C7 | BGR-gated II threshold | ✅ | #44 | running n(t) shifts E_th=E_th0−|K n^⅓| above gate; end-to-end stable |
| C8 | Dissipator sub-cycling | ✅ | #44 | m_sub from eph_numax·τ; II+e-ph split into m CPTP sub-steps; trace conserved at ν_sat=1e18. m=1 (unchanged) when e-ph off |
| **D** | Ring/pipeline MPI (systolic ring, one fused pass) | ✅ | #44 | Σ^HF via ring in super-mode; O(Nk/P) mem; ring==all-gather (6e-21), MPI 1==2 ranks (per-k 0.0). nl-II/eph/e-e accumulators slot in |
| **G** | Screening primitives (TF/Debye, Lindhard/RPA, LOPC) | ✅ | #44 | pure functions + GaAs/Si dielectric constants; unit-tested |
| **F** | Carrier-carrier (e-e/e-h) nonlocal CPTP Lindblad channel | 🚧 | #44 | design spec done; uses G screening; collision integral rides the ring (D) — pending |
| doc | Wiki pages 01–05 committed as long-term memory | ✅ | #44 | maintained per increment |

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
k-local super-mode COMPLETE (C0,C-prim,C2,C5,C6,C7,C8 ✅). Screening primitives
(Part G ✅) and the carrier-carrier design spec (Part F) are in. Every new
channel is CPTP, gated OFF by default.

**The remaining big block is the NONLOCAL infrastructure** (needed by nonlocal
II/e-ph AND by the carrier-carrier channel F). Do it in this order:
1. **C3 energy-windowed expanding-radius partner search** — a routine over the
   gathered all-k adiabatic energies returning energy-conserving partner
   weights (broadened bin), deterministic, no MC. Unit-test on a synthetic grid.
2. **D ring/pipeline MPI** — systolic ring (Plimpton) replacing the all-gather;
   ONE fused pass accumulating Σ^HF + nonlocal II + nonlocal e-ph + e-e. Start by
   refactoring the existing Coulomb all-gather into a ring with a single
   accumulator, verify 1-rank vs 2-rank identical (as the current Coulomb test
   does), then add accumulators.
3. **C4 nonlocal momentum-exchange II + nonlocal e-ph** on the ring.
4. **F carrier-carrier (e-e/e-h)** on the same ring: in/out screened-Coulomb
   collision integral with (1−ρ) Pauli factors and direct−exchange |W̃|²,
   Taj-Rossi CP-Markov + Rosati closure; screening ε(q) from Part G computed
   once/step from the gathered ρ; broadened-delta energy bins; Houston basis;
   predictor-corrector. CONSERVE Σf_k and ΣE_k f_k (validation invariants).
   Do NOT add a static screened-exchange shift (that stays in HF) — no double
   counting. Default screening = static Lindhard (G option b); LOPC (option c)
   GaAs-only, n ≳ 5e17.
Keep everything behind yn_sbe_superres OFF; validate against Chefonov Si
bleaching (wiki/04) before claiming physical correctness.

**Then the final major effort — nonlocal (C4) + ring MPI (D):**
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
