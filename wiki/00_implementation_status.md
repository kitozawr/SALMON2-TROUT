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
| **C** | Nonlocal super-compute (`yn_sbe_superres`) | ⬜ | — | big effort; break into C2/C3/C5 units below |
| C2 | Houston-basis adiabatic populations reuse | ⬜ | — | reuse existing ZHEEV of H_VG |
| C3 | Energy-bin final-state search (no tetrahedra/MC) | ⬜ | — | first testable unit of C |
| C4 | Nonlocal impact ionization (momentum exchange) | ⬜ | — | needs C3 + ring (D) |
| C5 | e-ph population-relaxing Lindblad (Si phonons, GaAs Fröhlich) | ⬜ | — | drives THz bleaching |
| C6 | T1/T2 bookkeeping (½ only on T1) | ⬜ | — | unit test: Hermitian dissipator conserves populations |
| C7 | BGR-gated II threshold | ⬜ | — | `yn_sbe_bgr_threshold`, gate 5e18 |
| C8 | Dissipator sub-cycling | ⬜ | — | when ν·(h/2) ≳ 0.2 |
| **D** | Ring/pipeline MPI (replace all-gather in super-mode) | ⬜ | — | one fused pass for Σ^HF + nl-II + nl-eph |
| **F** | e-e scattering — architecture TODO hook only | ⬜ | — | comment + ring hook, no coefficients |
| doc | Wiki pages 01–05 committed as long-term memory | 🚧 | #44+ | this commit establishes them |

---

## Decisions log (gotchas that bit us / must not be re-litigated)
- **Si gap:** a 3-parameter local EPM gives ~1.06 eV converged = Kunikiyo's own calc (1.068), NOT the 1.12 eV experimental value. The spec's "30 meV from 1.12" is against experiment; treat **1.068 eV (Kunikiyo calc)** as the real target. CBM at 0.86·2π/a along ⟨100⟩.
- **Si diamond:** V^A ≡ 0 for all shells → the existing zincblende `VS·cos + i·VA·sin` machinery is reused verbatim with V^A=0 (purely real diamond). τ=(a/8)(1,1,1) for both.
- **Coulomb Σ^HF uses δρ = ρ − ρ₀** (ρ₀ = diag(occup)) so Σ(t=0)=0 and the equilibrium EPM gap is unchanged. Verified Hermitian; 1-rank vs 2-rank populations bitwise-identical.
- **nex density** = (excited e per cell)/V_cell; the /N_k is already in `calc_trace` via /Σ(kweight) (BZ average). Do NOT divide by N_k again.
- **Max ionization ceiling** (32 e per cubic cell) = 1.77e23 cm⁻³ (full valence inversion).
- **GitHub CI is serial-only by request** (compiler/syntax gate). MPI verified locally.
- **Commits are signed-off** (`git commit -s`), trailers Co-Authored-By + Claude-Session.

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
| _(add per increment)_ | | |

End-to-end smoke (manual, not in run_all): scalar GaAs 4³, `yn_sbe_coulomb='y'`
+ `yn_sbe_hf_sublattice_proj='y'`, weak field → runs stable, finite, projection
diagnostic confirmed. Full "zero spurious Γ→X transfer" weak-field validation
is a longer run, deferred to a dedicated validation pass.

Run all: `python3 tests/run_all.py` (each test prints PASS/FAIL and exits nonzero on failure).

---

## Next action on resume
**Start Part C** — the nonlocal super-compute mode. It is large, so do it in
testable sub-units (commit each):
1. **Flags scaffolding** — add `yn_sbe_superres`, `yn_sbe_eph`, `sbe_eph_*`,
   `yn_sbe_bgr_threshold`, `sbe_bgr_*`, `sbe_search_sigma_e_ev` to
   salmon_global + inputoutput (default OFF). No behaviour change yet. Test:
   build + a run with all OFF is byte-for-byte identical to before.
2. **C3 energy-bin final-state search** — a standalone, unit-testable routine:
   given adiabatic energies {ε_a(k)} and a target ΔE, return the normalized
   broadened-bin weights over partner states (Gaussian/rectangle width σ_E ~
   mean level spacing). Deterministic, no MC. Test: weights are normalized,
   peak at energy match, symmetric.
3. **C5 e-ph rate functions** — pure functions ν(ε) (smooth saturation),
   Bose N_B(T), golden-rule deformation-potential rates (Si table) and the
   Fröhlich asinh factor (GaAs). Test each against the cited limits.
Then assemble C5 into a population-relaxing Lindblad half-step (Houston basis,
predictor-corrector C1/C6, sub-cycling C8), with the C6 Hermitian-dissipator
population-conservation unit test as the CPTP gate. Ring MPI (D) comes after a
working serial super-mode. Keep every sub-unit behind `yn_sbe_superres` OFF.
