# Open-Quantum-System Dissipation in Strong Fields — Literature Review

**Status: literature-review task (seeded 2026-07-19). Assigned to Fable to extend/maintain.
2026-07-20: the [B25] prescription is transcribed (§6) and implemented+validated in-repo (§7).**
This page collects the *verified* literature behind the "virtual-vs-real carrier"
problem and ranks the candidate universal fixes. It is a decision document for the
open coding task tracked in `wiki/00` ("universal virtual/real carrier separation").
Every citation below is tagged with how it was verified; **do not treat an
unverified entry as settled** — fetch the paper before relying on its equations.

---

## 1. The problem this review exists to solve

The momentum-resolved Markovian dissipators (the "ring": e-ph, impact ionization,
Auger) act on the **Houston-basis populations** `f_a(k,t)`. During a strong pulse
those populations contain, besides real excited carriers, the **field-induced
virtual polarization** ("dressing") — the reversible admixture of conduction
character that follows the field envelope and returns to the valence band at
switch-off. A Markovian collision operator scatters this dressing as if it were
real population, which manifests as:

- **Unphysical scaling.** Coarse-grid scan of exercise x12 (Si, sub-cycle deep-IR):
  the phonon-assisted generation `G_SBE = nex_ON − nex_OFF` scales as **`G ∝ I^0.94`**
  (near-linear), whereas true band-to-band tunneling is non-perturbative
  (Keldysh/Hurkx exponential `exp(−B/√I)`). Ratio `G_SBE / n_Keldysh` runs from
  `2×10⁷` at `3×10¹⁰` to `68` at `10¹² W/cm²` — a five-order-of-magnitude drift,
  so no constant rescale fixes the *shape* (see `samples/exercise_x12…/calibration_scan.py`).
- **Envelope-following dressing.** Coherent-only diagnostic (dissipation off):
  the Houston CB population correlates with `A²(t)` at **0.99** and **99.9 %**
  returns to the VB at switch-off. Reality is defined only by *post-pulse survival*,
  so a lifetime gate cannot separate dressing from real carriers in a sub-cycle
  field (verified: the `sbe_ring_gate_fs` gate passes ~60 % of the dressing).

The literature question: **in which basis / relative to which reference state
should the dissipator act so that reversible dressing carries zero excess
population and only genuine (non-adiabatic) excitation scatters — without
assuming a periodic (CW) field?**

---

## 2. Verified bibliography

Verification legend: **[V]** bibliographic facts (title/authors/journal/year)
confirmed via search; **[A]** claim confirmed at abstract level only;
**[U]** could not verify — treat as a lead, not a source.

### Non-Markovian machinery (maintainer-supplied PDFs, 2026-07-20)
- **[V/P] C. Meier, D. J. Tannor, "Non-Markovian evolution of the density
  operator in the presence of strong laser fields", J. Chem. Phys. 111, 3365
  (1999)** [MT99]. Nakajima–Zwanzig to 2nd order in the system–bath coupling,
  full memory + initial correlations, **arbitrary strong time-dependent
  fields**; the memory kernel is made time-LOCAL by an exponential
  decomposition of the bath correlation function → coupled primary + N
  auxiliary density matrices. Transcribed §8.1.
- **[V/P] M. Merkli, G. P. Berman, R. Sayre, "Electron transfer reactions:
  generalized spin-boson approach", J. Math. Chem. 51, 890 (2013)** [M13]
  (= [B25] ref [75]). Mathematically rigorous (dynamical-resonance method)
  dynamics of a donor–acceptor system coupled to a thermal bath, with BOTH
  energy-conserving (dephasing, g-type) and energy-exchange (a-type)
  couplings, all orders in the tunneling V, controlled remainders uniform in
  t; **multi-level donor/acceptor reduces exactly to a two-level model in the
  collective-state basis with √N-rescaled couplings**. Transcribed §8.2.

### Primary anchor — strong-field open systems
- **[V/P] N. Boroumand, A. Thorpe, G. Bart, A. M. Parks, M. Toutounji, G. Vampa,
  T. Brabec, L. Wang, "Strong field physics in open quantum systems",
  Rep. Prog. Phys. 88, 070501 (2025)** (= arXiv:2502.10240). Brabec/Vampa-group
  Letter on dephasing/dissipation for a strong-field-driven electron coupled to
  a bath. **[P] = transcribed from the maintainer-supplied journal PDF
  (2026-07-20); the full verified transcription is §6 below and the in-repo
  implementation + validation is §7.** Referred to as **[B25]** throughout.

### Gauge / basis of the SBE dissipator
- **[V] T. G. Jensen, L. B. Madsen, "Structure gauges and laser gauges for the
  semiconductor Bloch equations in high-order harmonic generation in solids",
  Phys. Rev. A 101, 053411 (2020)** (arXiv:2003.02961). Gauge freedom of the
  transition dipoles; gauge-invariant observables.
- **[V/A] "Semiconductor Bloch equations in Wannier gauge with well-behaved
  dephasing", arXiv:2508.07869 (2025).** **[A]** identifies exactly our failure
  mode: *the standard dephasing operator is ill-defined at band
  crossings/anticrossings and produces artifacts in the carrier distribution*
  because the operator's matrix elements change rapidly there. Directly relevant
  to the "conversion at Houston anticrossings" our scan localized.

### Field-dependent scattering (quantum kinetics)
- **[V] A. Di Carlo / W. Quade / F. Rossi / T. Kuhn** line on the **intracollisional
  field effect (ICFE)**: field-dependent, collision-broadened scattering rates
  from a density-matrix / NEGF collision integral. Concrete rate reference:
  **[V] J. A. Madureira, M. Bonitz et al., "Impact-ionization rates of
  semiconductors in an electric field: the effect of collisional broadening",
  J. Appl. Phys. 90, 829 (2001)** — the impact-ionization threshold *softens*
  and the rate is enhanced near threshold in a field. (The code's existing
  Franz–Keldysh softening `sbe_ii_fk_*` is the mean-field shadow of this.)

### Floquet open-system theory (background + caveat)
- **[V] N. Tsuji, T. Oka, H. Aoki, "Correlated electron systems periodically
  driven out of equilibrium: Floquet + DMFT formalism", Phys. Rev. B 78, 235124
  (2008).** Non-perturbative driven steady state (Falicov–Kimball demo).
- **[V] T. Oka, S. Kitamura, "Floquet Engineering of Quantum Materials",
  Annu. Rev. Condens. Matter Phys. 10, 387–408 (2019).** Review.
- **[V] A. Schnell, A. Eckardt, S. Denisov, "Is there a Floquet Lindbladian?",
  Phys. Rev. B 101, 100301(R) (2020).** **Key caveat:** a strict time-independent
  Floquet–Lindbladian *does not exist in general* — only in restricted parameter
  regions. Undercuts treating "Floquet–Lindblad" as a turnkey gold standard.

### Corrections to a mis-cited plan (record)
A prior AI-generated survey (Haiku) proposed four Floquet citations. Audit:
- **[U] "M. Schuler et al., Floquet-Born-Markov master equation, PRB 87, 035111
  (2013)"** — could not be confirmed; likely spurious. Real Floquet–Markov
  origins: **[V]** M. Grifoni & P. Hänggi, Phys. Rep. 304, 229 (1998); S. Kohler,
  J. Lehmann, P. Hänggi, Phys. Rep. 406, 379 (2005).
- **"Sentef et al., Nat. Commun. 6, 7047 (2015)"** — real number, **fabricated
  title**: the actual paper is *"Theory of Floquet band formation and local
  pseudospin textures in pump-probe photoemission of graphene"* (nothing about
  e-ph coupling). The attributed "phonons cool into the lower Floquet band"
  insight was invented.
- Tsuji-Oka-Aoki (2008) and Oka-Kitamura (2019) are real (pages 387–**408**, not
  413; author is **Takashi** Oka) but their attributed "insights" were loosely
  distorted.

---

## 3. Candidate universal fixes, ranked

### A. Dissipate in the instantaneous field-DRESSED (adiabatic) basis — RECOMMENDED
Rediagonalize the full instantaneous light-matter Hamiltonian `H(k,t)` (band
energies **plus** the interband `A·p`/dipole coupling) at each step and let the
dissipator act on the **dressed** populations. A state that adiabatically follows
the field stays in the lower dressed band ⇒ its excess population is **zero** ⇒ the
dressing is *never* scattered; only non-adiabatic (Landau–Zener) transfer to the
upper dressed band — the *real* carriers — is dissipated. Equivalently: relax the
density matrix toward the *instantaneous dressed reference* `ρ₀[A(t)]`, not the
field-free ground state.
- **Why universal:** parameter-free; it is the **instantaneous (non-periodic)
  limit of Floquet–Born–Markov**, so it needs no CW/periodicity assumption and
  survives sub-cycle THz (unlike adiabatic Floquet, which needs `Ω ≫ envelope
  rate`). The Houston unitary already produces the correct real (LZ) excitation;
  this only stops the ring from double-counting the reversible part.
- **Backing:** dressed-state master equations (Boroumand review §; Floquet-Markov
  in the CW limit); the "well-behaved dephasing in the right basis" argument
  (arXiv:2508.07869); ICFE (the field belongs *inside* the collision, i.e. in the
  states the collision connects).
- **Risk / open question (pending the code-basis audit):** whether the ring can be
  moved to the dressed basis cheaply (reuse an existing per-step diagonalization)
  or whether it requires a new per-k-per-step eigendecomposition + reworking the
  momentum-conservation map. Velocity-gauge dressed-state convergence in a
  truncated band window must be checked.

### B. Full Floquet–Born–Markov — REJECTED for this solver
Rigorous CW gold standard, but (i) strict Floquet–Lindbladian may not exist
(Schnell 2020); (ii) the adiabatic-Floquet approximation for pulses breaks down
for sub-cycle/few-cycle THz (giant non-adiabatic couplings) — precisely our regime;
(iii) it discards the Houston basis's exact sub-cycle field tracking. Keep as the
*theoretical explanation* of why Houston + Markov overcounts, not as the method.

### C. Intracollisional field effect / collisional broadening — COMPLEMENTARY
The quantum-kinetic (NEGF/density-matrix) route makes the scattering *rate itself*
field- and energy-dependent with a broadened threshold. It is the rigorous origin
of the correct field scaling but is **non-Markovian** (memory kernels) and costly.
The code already has the mean-field shadow (`sbe_ii_fk_*`, Franz–Keldysh softening).
Use it to *calibrate* option A, not as the primary mechanism.
**Concretized (2026-07-20): the [B25] SFSB memory kernel — the coherence-sector
member of this class — is implemented as the quantitative reference bracket
(§6–7, `yn_sbe_sfsb`).**

### D. Phenomenological BTBT source + interband knob — CURRENT FALLBACK
`sbe_eph_interband_scale` (gap-straddling rate factor) and/or a Hurkx/Kane source
term. Honest but **per-working-point** (the scan shows a single constant cannot
reproduce the field dependence). Retain as an explicit, documented approximation
until option A lands.

---

## 4. Recommendation

Pursue **option A** (instantaneous dressed-basis / adiabatic-reference dissipation)
as the universal, parameter-free fix, explained in the paper/wiki via **B**
(Floquet — why the artifact exists) and calibrated against **C** (ICFE/Keldysh
brackets), with **D** kept as the labelled fallback. Feasibility and the exact
code path depend on the ring-basis audit (in progress) recorded in `wiki/00`.

Suggested manuscript framing (once option A is validated):
> *The Markovian Lindblad dissipator in the Houston basis scatters field-induced
> virtual interband polarization, giving an unphysical near-linear generation rate
> instead of the non-perturbative Keldysh/Hurkx law. Floquet–Born–Markov theory
> resolves this by diagonalizing the light–matter coupling, but a strict
> Floquet–Lindbladian need not exist [Schnell 2020] and its adiabatic form breaks
> for sub-cycle THz transients. We therefore dissipate in the instantaneous
> field-dressed basis — the non-periodic limit of Floquet–Born–Markov — retaining
> exact sub-cycle field tracking while excluding reversible dressing from the
> collision integral.*

---

## 5. What is NOT verified (do before publishing)

- ~~The Boroumand review's **exact** dissipation prescription and equations~~
  ✅ **DONE (2026-07-20)** — transcribed from the journal PDF into §6 and
  implemented/validated (§7).
- Whether arXiv:2508.07869's "well-behaved dephasing" basis is *the* dressed basis
  or a Wannier construction (abstract only).
- The Schuler 2013 citation (treat as non-existent until found).
- The [B25] **supplement** (spectral-density definitions section IV, ZnO
  parameters table) is still NOT available — only the letter itself. The
  Ohmic/Debye profile normalizations in §6.2 are pinned by the letter's own
  printed anchors, not by the supplement.

---

## 6. [B25] — verified transcription (from the journal PDF, 2026-07-20)

### 6.1 The SFSB model

Single-electron two-band system + bosonic bath, linear diagonal coupling
[B25 Eq. (1)]:

$$H = -\tfrac12 E(K_t,t)\,\sigma_z + \tfrac12\hbar\Omega(K_t,t)\,\sigma_x
+ \sum_q \hbar\omega_q b_q^\dagger b_q + \sigma_z \sum_q g_q (b_q + b_q^\dagger)$$

with $K_t = K + eA(t)/\hbar$, generalized Rabi frequency
$\Omega = (2e/\hbar)\,d(K_t)E(t)$, transition dipole $d$ and gap $E$ at the
shifted momentum. The bath couples **only diagonally** (pure dephasing — no
direct bath-driven interband transitions); ionization is still affected
because laser and bath dynamics couple. After a polaron transformation +
interaction picture [B25 Eq. (2)]:

$$H_I = -\tfrac{E}{2}\sigma_z + \tfrac{\hbar\Omega}{2}
\left(\sigma_+ D^\dagger + \sigma_- D\right),\qquad
D = \exp\Big\{-\sum_q g_q\big(b_q^\dagger - b_q\big)/(\hbar\omega_q)\Big\}$$

— note the $1/(\hbar\omega_q)$ in the shift operator: the bath weight that
appears downstream is $(g_q/\hbar\omega_q)^2$.

### 6.2 The result: memory-kernel ionization

Second-order Dyson expansion, bath traced out ⇒ [B25 Eqs. (3)–(4)]:

$$n_c(K,t) = \tfrac12\,\mathrm{Re}\int_{-\infty}^{t}\! dt_1
\int_{-\infty}^{t_1}\! dt_2\; \Omega^*(K_{t_1},t_1)\,\Omega(K_{t_2},t_2)\,
\exp\big[\,i S(t_1,t_2) + C(t_1-t_2)\,\big],\qquad
n_c(t)=\int_{BZ} n_c(K,t)\,dK$$

with the dynamic-Stark-shifted action
$S(t_1,t_2)=\int_{t_2}^{t_1} E_s(K_\tau,\tau)\,d\tau/\hbar$,
$E_s=\sqrt{E^2+|\hbar\Omega|^2}$. **ALL environment influence is in the bath
correlation function** [B25 Eq. (5)]:

$$C(t_1-t_2) \approx \int_{-\infty}^{\infty} J(\omega)\Big[\, i\sin\omega(t_1{-}t_2)
- \big(1-\cos\omega(t_1{-}t_2)\big)\coth\tfrac{\hbar\omega}{2k_BT} \Big]\,d\omega$$

- $e^{C(\tau)}$ does **not** factorize across time steps — this pair-time
  structure IS the non-Markovianity ("the future evolution depends on the past
  history").
- Temperature enters only via the coth.
- $J(\omega)$ = continuum spectral density with two parameters: dimensionless
  coupling $j_o$ and cutoff $\omega_c$. Models named: Debye, Ohmic,
  Under-Damped Brownian, Gaussian, Shifted-Gaussian (**definitions are in the
  unavailable supplement §IV** — see the normalization note below).
- **Printed anchor** [B25 §2]: for the Debye bath at high T,
  $C \to -(t_1-t_2)/T_2$ with $T_2 = \hbar/(2\pi k_B T j_o)$ — the
  relaxation-time approximation is this limit. The other models do NOT become
  linear-in-t at high T.
- Ionization ratio $\eta = n_c(j_o\neq0)/n_c(j_o=0)|_{t=\infty}$ [B25 Eq. (6)].

**Normalization (ours, anchored):** with the letter's $D$-operator convention
the continuum weight is $W(\omega) = \sum_q (g_q/\hbar\omega_q)^2\,
\delta(\omega-\omega_q)$ extended oddly; we implement
$W(\omega)=j_o\,g(|\omega|)/\omega$ with $g_{\rm ohmic}=e^{-\omega/\omega_c}$,
$g_{\rm debye}=\omega_c^2/(\omega^2+\omega_c^2)$. For ANY cutoff with
$g(0)=1$ this reproduces the printed high-T anchor
$T_2=\hbar/(2\pi k_BTj_o)$ **exactly** (test_bath_corr), plus the closed
forms $\mathrm{Im}\,C_{\rm ohmic}=2j_o\arctan(\omega_c t)$,
$\mathrm{Re}\,C_{\rm ohmic}(T{=}0)=-j_o\ln(1+\omega_c^2t^2)$,
$\mathrm{Im}\,C_{\rm debye}=\pi j_o(1-e^{-\omega_c t})$. The supplement's own
Gaussian/Shifted-Gaussian/UDB definitions are NOT implemented (no source).

### 6.3 The letter's results (what a correct implementation must reproduce)

1. **Dephasing ionization** [Fig 1]: RTA (T₂ = 6 fs) overestimates
   under-resonant ionization by ~5 orders of magnitude (ZnO, 9-photon,
   E₀ = 5×10⁸ V/m) — "when the phase relationship is disrupted, virtually
   excited electrons are prevented from returning to the ground state."
2. **Physical baths are benign at realistic T** [Fig 2(a)]: Ohmic/UDB/
   Gaussian/Shifted-Gaussian all give negligible η changes at low T
   ("as detected by experiments"); enhancement appears only at extreme
   temperatures. Debye + RTA rise unphysically (the long ω-tail).
3. **Debye → RTA at very high T** [Fig 2(b)].
4. **Suppression at low T / strong coupling** [Fig 3]: the phase of C(t)
   (Im part) acts as a **dynamic bandgap addition** ⇒ dephasing-SUPPRESSED
   ionization; setting Im C = 0 flips suppression into enhancement [Fig 3(c)].
5. **Distinct k-signature** [Fig 4]: dephasing ionization spreads nc(K) far
   beyond the coherent distribution (ARPES-detectable).
6. **Bath influence shrinks in the tunneling regime** [Fig 5(b)]: at γ < 1
   the heat bath cannot follow the sub-cycle ionization step.
7. Stated limits: 2nd-order Dyson (small n_c); thermal-equilibrium bath
   (no dynamic T(t)); no e-e scattering (needs a fermionic bath / Keldysh);
   two-band (multi-band extensions cited: their refs [74–76]).

### 6.4 Connection to the §1 problem (why this page cares)

[B25]'s "dephasing ionization" IS the §1 virtual/real disease in the
**coherence sector**: the Markovian T₂ real-ifies the returning virtual
polarization exactly like the Markovian ring collision operators real-ify the
Houston dressing. And the paper's fix is the same *kind* of fix option C
anticipated: keep the field-driven dynamics exact, make the environment
non-Markovian (memory kernel), and the pathological conversion disappears at
realistic temperatures. The SFSB model is therefore the quantitative
**reference bracket** for how much generation a physical environment can
add — option A's dressed-basis dissipator should land between the coherent
result and the SFSB-with-physical-bath result, far below RTA.

---

## 7. SFSB in this repo — implementation + validation (2026-07-20)

### 7.1 What was implemented

| piece | where | what |
|---|---|---|
| `bath_corr_table` | `sbe_superres_ssbe.f90` | C(τ) table [B25 Eq. (5)]: ohmic / debye / rta profiles, Simpson ω-integral with the ω→0 limit node, coth series switches; `bath_t2_high_t` = ħ/(2πk_BTj_o) |
| `sfsb_nc_series` | `sbe_superres_ssbe.f90` | pure Volterra stepper for Eq. (3): e^{iS} factorizes across steps, e^{C} does not — true history sum, truncatable window |
| `yn_sbe_sfsb='y'` mode | `sfsb_ssbe.f90` | k-line driver: 1D MP line (num_kgrid=(N,1,1), E∥b1 enforced), K_t=K+A(t) trajectory with periodic **cubic-spline** interpolation, auto stride from the Stark gap, MPI×OpenMP over k, `_sfsb_nex/_sfsb_nck` outputs |
| inputs | `&sbe` | `sbe_bath_model/jo/wc_ev/temperature_k/rta_t2_fs/memory_fs`, `yn_sbe_bath_imc`, `sbe_sfsb_nv/nc/stride` |
| tests | `test_bath_corr.f90`, `test_sfsb_kernel.f90` | kernel closed forms + high-T anchor; stepper vs exact TDSE (RK4, 0.003%), RTA kernel ≡ Markov ODE (RK4), dephasing ionization ×1700 on a toy two-level, T=0 suppression, Im-C flip, window truncation |
| example | `samples/exercise_x13_GaAs_sfsb_nonmarkovian/` | GaAs Γ–L, the paper's parameter points + `sfsb_validation.png` |

### 7.2 The two-band-reduction traps (do not relitigate)

Raw energy-sorted GS lines are NOT a two-band model; two traps found and
solved on the way (both produced grid-DIVERGENT, orders-too-large nex):

1. **Degenerate-manifold gauge randomness.** Individual |d_vc| between
   members of a (near-)degenerate manifold are gauge-random per k (the GS
   solver returns arbitrary per-k mixtures); sorted band indices also swap
   character at avoided crossings (CdS: v7/v8 min splitting 3.6e-5 Ha; the
   |d| line stepped 2.3 → 1e-9 between adjacent points). Only the
   **quadrature (bright-state) sum** over the manifold is gauge-invariant and
   smooth ⇒ coupling line = √Σ_v|d·ê|², gap line = band-edge E_c−E_v,top
   (continuous by sorting). With this construction GaAs Γ–L converges:
   nex = 2.957/2.962/2.9625e18 cm⁻³ at 96/192/384 points.
2. **Selection-rule walls.** On CdS Γ–M (E⊥c) the whole band-edge channel is
   allowed only in a pocket |q₁| ≲ 0.1 and **exactly zero** beyond (the bright
   conduction character migrates diabatically to c10+): no smooth two-band
   line exists on that line/polarization. The mode detects the wall
   (adjacent-point coupling jump > 50 % of range) and WARNS; use a clean line
   (GaAs Γ–L) or another polarization.

### 7.3 Validation (GaAs Γ–L, 192 pts, λ₀ = 3.2 µm, E₀ = 2×10⁸ V/m, γ ≈ 2.15)

All seven [§6.3] claims reproduce — figure
`samples/exercise_x13_GaAs_sfsb_nonmarkovian/sfsb_validation.png`:
RTA(6 fs) η = **30.6**; ohmic j₀=1/2.1ω₀/300 K η = **0.041** (suppressed);
Im C := 0 flips to η = **158**; same bath at 2×10⁴ K η = **101**; ohmic/debye
j₀=0.1/0.1ω₀ give η ≈ 1.08 for T ≤ 3×10³ K rising only above 10⁴ K; nc(K) of
dephasing ionization is a flat BZ-wide background vs the sharp coherent peak
[Fig-4 signature]; at γ ≈ 0.43 the RTA excess collapses 30.6 → 1.4 [Fig-5(b)].
Debye at 3×10⁴ K is enhanced (η 5.3) but still ×29 below the derived-T₂ RTA —
full convergence Debye→RTA needs T beyond our scan and the τ ≫ 1/ω_c regime;
recorded as a quantitative difference from [B25 Fig 2(b)], not a discrepancy
in trend.

### 7.4 Route into the multiband master equation (proposed, NOT yet authorized)

A proposal circulated to bolt the kernel into the Lindblad as
∂ρ/∂t|diss = ∫K(t−t′)ρ(t′)dt′ with K(τ) ∝ Re[e^{C(τ)}]·L_jump ρ L_jump†.
**Rejected as formulated**, three reasons:

1. **Wrong object.** e^{C(τ)} is the *polaron cumulant* attached to the two
   DRIVE vertices Ω(t₁)Ω(t₂) (that is what [B25 Eq. (2)–(3)] says: after the
   polaron transform the bath rides the σ± drive term). A Nakajima–Zwanzig/
   Born memory kernel on a jump operator would carry the *second-order* bath
   correlator ⟨B(τ)B(0)⟩, not the exponentiated cumulant. Mixing the two
   expansions has no derivation.
2. **Re[·] alone is provably wrong physics.** Killing Im C flips suppression
   into enhancement (test_sfsb_kernel: 5.2e-12 vs 0.16, eleven orders;
   validation: η 0.041 → 158). The bath phase (dynamic gap addition) is not
   decoration — it is the mechanism [B25 Fig 3(c)].
3. **CP safety.** A raw time-nonlocal ρ-history kernel is not
   CP-by-construction — the negative-population/trace-drift diseases the
   frozen-window work (PRs #80–#87) just eliminated would return. And storing
   ρ(t′) history for all (k, nb²) is prohibitive.

**The correct Stage-2 route (bounded, CP-safe, no history storage):**
Meier–Tannor exponential decomposition ([B25 ref 45]; Meier & Tannor,
J. Chem. Phys. 111, 3365 (1999)) of the kernel,
$e^{C(\tau)} \approx \sum_j a_j e^{-\gamma_j\tau}$ (complex γ_j, fitted once
from `bath_corr_table`). The inner history integral then becomes a handful of
**auxiliary polarization ODEs per (k, v-c pair)**:
$\dot p_j = -(i\Delta E + \gamma_j)p_j + i a_j (\Omega/2) w$, with the
coherence ρ_cv = Σ_j p_j — memory WITHOUT storing history, k-local, exactly
reduces to the current RTA/Kuhn-Zurek for a single real exponential, and each
sub-step is a damped linear map (CP-compatible stepping). Cost: n_aux (~3–6)
extra complex fields per (k, pair). This slots into the existing coherence
update; the population (ring) channels keep their separate ICFE story (§3C).
**Awaiting maintainer go-ahead.**

---

## 8. Stage-2 PROPOSAL — non-Markovian dephasing inside the multiband master equation (2026-07-20, AWAITING MAINTAINER APPROVAL)

The maintainer supplied the two sources ([MT99], [M13]) and asked for a
thought-through fix. This section is the design; no solver code is written yet.

### 8.1 [MT99] — verified transcription (the machinery)

Model: Caldeira–Leggett bath, coupling H_sb = −f(x)·Σc_i x_i, plus an
**arbitrary strong time-dependent field W(x,t) treated non-perturbatively**.
NZ projection to 2nd order in the coupling λ gives the time-nonlocal master
equation [MT99 Eq. (13′)]

$$\dot\rho_s(t) = L_s^{\rm eff}(t)\rho_s(t) + \int_{-\infty}^{t}\! K(t,t')\rho_s(t')\,dt',$$

$$K(t,t') = \lambda^2 L^-\big(a(t{-}t')\,\mathcal{T}e^{\int_{t'}^{t}L_s}L^-
+ b(t{-}t')\,\mathcal{T}e^{\int_{t'}^{t}L_s}L^+\big),$$

with $L^-=-i[f,\cdot\,]$, $L^+=[f,\cdot\,]_+ -2\chi$ (χ = thermal mean of f),
and the bath correlation function c(t) = a(t) − i b(t) [MT99 Eq. (8)]. Two
essential structural facts: **the field enters the memory kernel** (bath
influence is modified by strong fields), and the kernel is nonlocal — naive
storage of ρ(t′) does not even suffice (the field-dressed backward propagator
is needed).

The fix [MT99 Sec. IV]: write a(t) = Σ_k α_k^r e^{γ_k^r t},
b(t) = Σ_k α_k^i e^{γ_k^i t} (their Lorentzian parametrization of J(ω),
Eq. (15)–(17), incl. Matsubara terms; **"we view Eq. (15) merely as a
numerical decomposition"** — a direct numeric exponential fit is endorsed).
Then auxiliary matrices ρ_k (Eqs. (18)–(21)) convert the problem to the
time-LOCAL coupled system [MT99 Eq. (22)]:

$$\dot\rho_s = L_s^{\rm eff}\rho_s + \lambda\Big(\sum_k \alpha_k^r L^-\rho_k^r
+ \sum_k \alpha_k^i L^-\rho_k^i\Big),$$
$$\dot\rho_k^r = (L_s + \gamma_k^r)\rho_k^r + \lambda L^-\rho_s,\qquad
\dot\rho_k^i = (L_s + \gamma_k^i)\rho_k^i + \lambda L^+\rho_s.$$

Properties (all printed in [MT99]): the external field enters only
block-diagonally (each matrix rides the SAME system propagator); cost =
(N+1)× a Markovian propagation for arbitrarily long memory; genuine
irreversible decay + the correct (to O(λ²)) thermal equilibrium; initial
correlations = "memory from the past" (the t<0 field-free equilibrium
history), naturally included by starting the coupled system in its field-free
steady state.

### 8.2 [M13] — verified transcription (the rigor + multilevel)

Generalized spin-boson [M13 Eq. (1.3)]:

$$H = \begin{pmatrix}E_D & V\\ V & E_A\end{pmatrix} + H_R
+ \lambda\begin{pmatrix}g_D & a\\ a & g_A\end{pmatrix}\otimes\varphi(h),$$

g-type (diagonal) coupling = energy-conserving **dephasing**; a-type
(off-diagonal) = **energy exchange** (bath-driven transport even at V = 0 —
the physics missing from a pure-σ_z model). Dynamical-resonance method: the
full reduced density matrix for ALL t ≥ 0, perturbation theory in λ with
controlled remainders **uniform in t**, all orders in V; validity |λ| ≪
E_D−E_A. Relaxation rates [Eqs. (1.8)/(1.10)] carry the bath through
coth(βΔ/2)J(Δ) at the SYSTEM transition frequency Δ = √((E_D−E_A)²+4V²).
**Multilevel acceptor** [Eqs. (1.12)–(1.14)]: an N_A-fold degenerate acceptor
reduces EXACTLY, in the basis of the donor state and the collective state
σ_A = N_A^{-1/2}Σ|A_i⟩, to the two-level model with rescaled couplings
V → V√N_A, a → a√N_A. Rate grows with N_A (γ ∝ N_A in the small-hopping
regime [Eq. (1.16)]); the transfer separation does NOT [Eq. (1.11)].

**Retroactive grounding of #88:** the SFSB mode's bright-manifold quadrature
coupling √Σ_v|d·ê|² IS [M13]'s collective-state √N rescaling — the two-band
reduction we derived from gauge invariance is the rigorous multilevel
reduction of [M13].

### 8.3 The proposed fix, mapped onto the solver

**What is being fixed.** The production dephasing is Markovian: Kuhn–Zurek
`exp[−λ(X_a−X_b)²τ]` with constant λ = k_BT/τ_m (a white-noise/high-T limit)
plus constant per-channel gout damping. [B25] (§6) shows this class fabricates
ionization ("dephasing ionization") by real-ifying returning virtual
polarization — the coherence-sector twin of the §1 disease. The merged
`yn_sbe_sfsb` mode (#88) brackets the correct answer but is 2nd-order/1D —
not the production propagator.

**The fix: an MT auxiliary-mode non-Markovian dephasing channel** replacing
the constant-rate dephasing, [MT99] structure with [B25]'s bath:

1. **Coupling operator** f = diag(s_a/2), s_a = +1 conduction / −1 valence —
   the multiband σ_z of [B25], rigorously justified for the (near-)degenerate
   band-edge manifolds by [M13]'s collective reduction. **No a-type
   (energy-exchange) bath coupling**: bath-driven interband population
   transfer is already modeled by the collision channels (e-ph/II/Auger) —
   adding it here would double-count (same rule class as BGR↔Σ^HF).
   Optional later: f from the Kuhn–Zurek branch positions X_a (the
   colored-noise upgrade of the exact same channel; shares 100% of the
   machinery since both f are diagonal).
2. **Kernel data**: c(τ) ≡ −C″(τ) of the SAME shipped bath models
   (W(ω)=j₀g(|ω|)/ω; ohmic/debye, T₂ anchor). With f = σ_z/2 the exact
   undriven coherence decay is exp[C(t)] **including the Im C phase** (the
   dynamic gap addition, [B25 Fig 3c]) — since C(t) = −∫₀ᵗ(t−τ)c(τ)dτ
   identically. This is the sign/normalization anchor test.
3. **Exponential decomposition**: numeric Prony / least-squares fit of
   a(τ) = Re c, b(τ) = −Im c to Σα e^{γτ} ([MT99]-endorsed), n_exp ≈ 3–6 +
   Matsubara at finite T; fit-quality unit test against the table.
4. **Wiring** (concrete, from `dt_evolve_bloch_cf4`): auxiliary matrices
   σ_k(nb_full or nba?, per k) stored as `sbe%rho_mt(:,:,:,k)`;
   Strang: [half MT-coupling + e^{γ_k dt/2} decay] → **the existing S4/CF4
   unitary applied to ρ_s AND every σ_k** (refactor `cf4_unitary_step` to
   build the two exponentials once per (k, sub-step) and apply to 1+n
   matrices — the exponentials dominate, so the marginal cost is GEMMs) →
   [half MT-coupling]. Because f is DIAGONAL, the coupling sub-step acts
   ELEMENTWISE: (L⁻ρ)_{ab} = −i(s_a−s_b)/2·ρ_{ab}, (L⁺ρ)_{ab} =
   ((s_a+s_b)/2−2χ)ρ_{ab} — each (a,b) element of (ρ_s, {σ_k}) closes into a
   small (1+n)-dimensional constant-coefficient linear ODE whose exponential
   is precomputed ONCE per element-class (s_a−s_b ∈ {0,±2}) per step — O(1)
   overhead. Populations (s_a=s_b diagonal) are untouched by the coupling —
   pure dephasing, exactly [B25]'s structure (the bath dephases, the LASER
   converts).
5. **Conservation/honesty**: trace conservation EXACT (only L⁻ commutators
   feed ρ_s); Hermiticity kept by conjugate-pair exponentials + the existing
   re-Hermitization. Positivity is NOT guaranteed by construction (TC2/
   Redfield class — the fundamental non-Markovian↔CP tension, Schnell 2020):
   weak-coupling validity documented (|λ| ≪ gap, [M13 Eq. (1.2)]), runtime
   Houston-diagonal monitor (warn/stop per existing discipline).
6. **Initial correlations**: start the coupled (ρ_s, σ_k) system in its
   field-free steady state (pre-relax a few bath correlation times before
   t=0) — [MT99]'s "time flows from −∞" form; cheap and clean.
7. **Frozen window**: f and the auxiliaries live in the active window (same
   scope as every dissipator); the CP-extension bookkeeping of the
   active↔frozen coherences applies unchanged.
8. **Guards**: mutual exclusion with `sbe_decoh_*` (two dephasing models on
   ⇒ error stop); requires `sbe_bath_model` ohmic/debye (rta would be the
   Markovian limit — allowed for A/B testing); graphene: replaces the
   forbidden KZ (many-body coherence argument re-examined — the MT channel
   with a physical bath may be legitimate there; needs its own decision).
9. **Cost**: ×(1+n_exp) memory for ρ (e.g. ×5), runtime ×3–5 (GEMM-bound);
   acceptable as a super-mode option.

### 8.4 Validation plan (after approval)

1. Undriven two-level, weak coupling: coherence decay ≡ exp[C(t)] (Re decay
   AND Im phase) — pins signs/normalization against the merged C(τ) table.
2. High-T Debye: MT channel → constant-T₂ behaviour (the RTA limit).
3. Weak-field driven two-level: full-SBE MT channel vs the `yn_sbe_sfsb`
   reference (#88) — must agree in the 2nd-order regime (this cross-check is
   exactly what the SFSB mode was built to provide).
4. Production demo (the headline): Si/GaAs 4³ multiband, below-gap pulse —
   nex(coherent) vs nex(KZ 300 K) vs nex(MT ohmic 300 K); expected
   MT ≈ coherent ≪ KZ (the fabricated dephasing ionization removed from the
   full solver).
5. [M13] check: rate vs conduction-manifold size (γ ∝ N in the small-V
   regime) on a toy multiband spectrum.

### 8.5 Decision points for the maintainer

- (a) approve the channel as specified (f = band-character σ_z/2, no a-type
  coupling)? (b) also want the X_a-coupled variant as an option?
- accept the TC2 positivity caveat (monitored, weak-coupling domain) — or
  require a CP-enforcing projection step after each dt?
- graphene: allow the MT channel where KZ is forbidden, or keep both off?
- bounded-increment split: (1) kernels c(τ) + Prony fit + tests; (2) solver
  wiring + CF4 refactor; (3) validation + docs — one PR each or one branch?
