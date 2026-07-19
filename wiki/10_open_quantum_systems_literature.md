# Open-Quantum-System Dissipation in Strong Fields — Literature Review

**Status: literature-review task (seeded 2026-07-19). Assigned to Fable to extend/maintain.**
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

### Primary anchor — strong-field open systems
- **[V] N. Boroumand, A. Thorpe, G. Bart, A. M. Parks, M. Toutounji, G. Vampa,
  T. Brabec, L. Wang, "Strong field physics in open quantum systems",
  Rep. Prog. Phys. 88, 070501 (2025)** (= arXiv:2502.10240). Brabec/Vampa-group
  review of dephasing/dissipation for a strong-field-driven electron coupled to a
  bath. **[A]** frames the environment as a heat bath and dephasing as
  environment-induced loss of phase coherence. *This is the authoritative recent
  reference; its dissipation prescription must be transcribed directly from the
  PDF (in-repo upload) before we cite specific equations.*

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

- The Boroumand review's **exact** dissipation prescription and equations
  (fetch blocked; PDF is in the session upload — transcribe directly).
- Whether arXiv:2508.07869's "well-behaved dephasing" basis is *the* dressed basis
  or a Wannier construction (abstract only).
- The Schuler 2013 citation (treat as non-existent until found).
