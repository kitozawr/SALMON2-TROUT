# Numerical Methods

> How the [Physics Models](01_physics_models.md) equations are integrated. §1–4, 10–12 ✅ implemented; §5–9 🚧 planned (super-mode). Constants in [Constants](02_constants.md).

## 1. Operator splitting (master scheme) ✅
Per timestep h:
> ρ(t+h) = D(h/2) ∘ [ S2(p1 h) ∘ S2(p2 h) ∘ S2(p1 h) ] ∘ D(h/2) [ρ(t)]

The inner triple product is the **Suzuki-Yoshida 4th-order composition** of the 2nd-order symmetric unitary step S2; D is the dissipative step. **Yoshida wraps ONLY the unitary part:** the middle sub-step has p2·h < 0; for a unitary this is a harmless backward-time rotation, but applying it to the dissipator would turn exp[−λΔX²τ] into exp[+...], not PSD, breaking CPTP. So D is applied only with τ = +h/2. [Yoshida, PLA 150, 262 (1990); Hatano & Suzuki, LNP 679, 37 (2005)]

## 2. CF4 unitary step S2 (commutator-free Magnus 4th order) ✅
On Gauss-Legendre nodes c1,2 = ½ ∓ √3/6:
- H1 = H_VG(t + c1 τ), H2 = H_VG(t + c2 τ)
- Ω1 = τ(α1 H1 + α2 H2), Ω2 = τ(α2 H1 + α1 H2), α1 = ¼+√3/6, α2 = ¼−√3/6
- ρ → exp(−iΩ2) exp(−iΩ1) ρ exp(+iΩ1) exp(+iΩ2)

Each exponential built **exactly from an eigendecomposition** of the Hermitian generator (unitary to machine precision). Two ZHEEV per S2. [Blanes & Moan, JCAM 142, 313 (2002); Alvermann & Fehske, JCP 230, 5930 (2011)]

**Why CF4 over ETDRK4:** the huge wave-packet excursion ΔX and high VG orbitals force a tiny ETDRK4 step (limited by 1/E_g and A·p). CF4 exponentiates H_VG, so the step is limited only by the field frequency.

## 3. Dissipative step D — Kuhn-Zurek (Hadamard/Gaussian, exactly CPTP) ✅
- Diagonalize H_VG(t) → Houston basis U(t), {ε_a}.
- ρ̃ = U†ρU; multiply ρ̃_ab ← exp[−λ(X_a−X_b)²τ] ρ̃_ab (PSD Gaussian kernel → CPTP for τ≥0).
- Rotate back; advance X_a += ½(V_a(t)+V_a(t+h))h (midpoint velocities, matching CF4 4th order). V_a = (U†πU)_aa + A(t).

## 4. Houston (adiabatic) basis for all dissipative channels ✅
All dissipative channels operate in the **same instantaneous adiabatic basis U(t)** from the ZHEEV of H_VG already computed for CF4 — **no extra diagonalization**. Energies, velocities, Pauli blockers from **adiabatic populations ρ̃ = U†ρU**, not Bloch populations (virtual in VG). [Yue & Gaarde, JOSA B 39, 535 (2022)]

## 5. Nonlocal super-compute dissipator: insertion and CPTP 🚧
In super-mode the nonlocal II + e-ph dissipators are a **single Strang half-step OUTSIDE the k-grid**:
> [ D_nl(h/2) ] ∘ [ existing per-k step (CF4 + k-local D) ] ∘ [ D_nl(h/2) ]
- **CPTP:** each exp(τ D_nl) with a GKLS generator is CPTP for any τ≥0; the whole composition is CPTP regardless of k-coupling. [Lindblad, CMP 48, 119 (1976)]
- **Order:** Strang 2nd-order; caps the dissipative sector at 2nd order while unitary CF4 stays 4th. D is a perturbative correction.

## 6. Predictor-corrector (NOT full freeze) — protects Pauli blocking 🚧
At ν ~ 1e14 s⁻¹ populations and (1−ρ) change within a step; freezing rates can drive (1−ρ) < 0.
1. Evaluate rates and Σ^HF at ρ(t) (predictor).
2. Take the dissipative half-step.
3. Re-evaluate at ρ(t+h/2) (corrector); apply.
4. **Clamp all (1−ρ) Pauli factors to [0,1].**
2nd-order, keeps eigenvalues in [0,1]. [Rosati-Iotti-Dolcini-Rossi, PRB 90, 125140 (2014)]

## 7. Energy-conserving final-state search — ENERGY BINS (no tetrahedra, no Monte-Carlo) 🚧
- Replace δ(ΔE) by a **normalized broadened bin** (Gaussian or unit-area rectangle, width σ_E ~ mean inter-level spacing). [Stobbe 0.2 eV rects; Kunikiyo 5 meV bins]
- **Full deterministic pair enumeration** is the reference/fallback for small grids.
- **No Monte-Carlo** — injects stochastic noise, breaks exact CPTP/trace.
- **Energy-windowed expanding search radius:** small partner shell near threshold; radius grows with excess energy; global sweep at high energy. Partners found **by energy, never by hardcoded valley coordinates** → material-universal (Si X-valleys 0.85·2π/a ⟨100⟩; GaAs L-valleys ⟨111⟩).

## 8. Sub-cycling the dissipator 🚧
When ν_max(h/2) ≳ 0.2, split the dissipative half-step into m CPTP sub-steps [exp((τ/2m)D)]^m. Each CPTP; choose m so (τ/2m)ν_max ≲ 0.1. [Lindblad 1976]

## 9. Ring/pipeline MPI (super-mode; replaces all-gather) 🚧
All-gather stores a full copy of ρ on every rank → RAM overflow (HF near node-memory limit). Use a **systolic ring**:
- Each rank owns its N_k/P block + **one transit buffer**. Blocks circulate via MPI_Sendrecv in P steps; each step computes partial pairwise contributions between resident and transit block, then forwards. After P−1 hops the full nonlocal pairwise sum is complete.
- **Memory O(N_k/P + one block)** — does NOT grow with P. Communication O(P) sends of O(N_k/P).
- **One fused ring pass serves ALL nonlocal sums:** Σ^HF + nonlocal II + nonlocal e-ph in a single circulation, with the §6 predictor-corrector.
- **Active-subspace compression:** circulate only n_act×n_act gap-edge blocks. [Plimpton, JCP 117, 1 (1995); Plimpton & Hendrickson, JCC 17, 326 (1996)]

## 10. Hartree-Fock sublattice-block projection (cubic-cell folding fix — Part E) 🚧
In the 8-atom cubic cell, bands fold 4-fold. An unrestricted Σ^HF spuriously couples states from different primitive-BZ sectors. A translationally invariant Coulomb operator conserves primitive crystal momentum → the inter-sublattice Fock element is **exactly zero**; any nonzero value is a folding artifact.

**Fix:** project Σ^HF block-diagonally onto the 4 FCC sublattice sectors — keep intra-sublattice, zero inter-sublattice. Equivalent: unfold ρ to primitive k via the spectral weights w_s (`SYSNAME_unfold.data`), apply Σ^HF with primitive momentum conservation, fold back. **Validation:** with projection ON, a weak-field (<100 kV/cm) run shows zero spurious Γ→X/Γ→L transfer from exchange. [Popescu & Zunger, PRB 85, 085201 (2012); Ku-Berlijn-Lee, PRL 104, 216401 (2010)]

## 11. Exact current operator ✅
j(t) = −Σ_k Tr[(π + A(t)) ρ^k] (a.u.), with the SAME π (incl. v_SO for spinor) used in the dynamics — no perturbative expansion. Hermiticity (ρ=ρ†) enforced each step. Spectrum S(ω) ∝ ω²|j(ω)|².

## 12. Frozen-core / active-subspace optimization ✅

**Frozen core reduces the Houston/dissipator work — it is NOT a velocity-gauge
cutoff.** The Strang step `D(h/2) U(h) D(h/2)` splits into two subspaces:

* **U(h) — the reversible velocity-gauge unitary — runs on the FULL band
  basis** (`dt_evolve_bloch_cf4`: `build_HVG(nb, …)` + `cf4_unitary_step(nb, …)`
  on all `nb` bands). This is the whole point of the frozen scheme: a strong
  field with a large `A(t)` pushes population **up/down through the "frozen"
  high/deep bands and brings it back** — the frozen bands hold the field-dressed
  virtual population, so **VG basis sufficiency is preserved**. The pure-SBE
  unitary is cheap (O(nk·nb³), linear in nk), and the EPM already lets you pick
  `nstate` large enough for any `A`.
* **D(h/2) — the CPTP dissipators — act only in the active Houston window**
  (bands inside `frozen_core_threshold_ev`/`frozen_free_threshold_ev`). Each half
  step truncates ρ to the active block, transforms to the instantaneous Houston
  basis, applies the channels, and glues the block back into the full matrix.
  This is where the saving lives: fewer ZHEEV, a smaller Σ^HF (its cost/traffic
  scale as n_act², so ALWAYS pair `yn_sbe_coulomb` with a frozen window), and —
  the big one **in super mode** — a smaller ring, whose O(nk³·n_act) collision
  sum shrinks directly with the active-band count.

There is **no freeze-reset**: inactive bands evolve reversibly under U and are
simply skipped by the dissipators. Consequently, with **clean dynamics** the
frozen-core flag has *zero* effect on the current (verified to machine precision:
a 5×10¹² V/cm Si run with 3/8 bands frozen reproduces the all-active current to
1.7×10⁻¹⁵ relative). Population output still reads the four gap-edge bands
regardless of how many are active — so freeze everything except the gap edge,
keep `nstate` large for the basis, and pay only for the window you dissipate.

*(Historical note: earlier builds truncated U to the active block and reset the
frozen bands to ground occupation every step — that WAS a VG cutoff and broke
basis sufficiency at strong fields; fixed 2026-07-12.)*

**Exact current.** J = Tr[(π + A)ρ] (a.u.) — no perturbative splitting; the
velocity-gauge inter/intra-band compensation is exact. Hermiticity stabilization
(ρ = ½(ρ+ρ†)) is enforced each step for real currents and FFT stability.

## 13. Time-step recommendation: 4 active bands (2 VB + 2 CB, all unfrozen in VG) at strong fields ✅

Four constraints set `dt`; with a 2V+2C Houston window the **field-dressed
band-span phase** is almost always the binding one. All formulas in Hartree
atomic units (1 a.u.t = 0.0241889 fs); the showcase inputs are in `A_eV_fs`,
conversions given inline.

**1. The binding constraint — dressed active-window span.** The CF4/Magnus
step must resolve the fastest interband phase of the *dressed* Hamiltonian
$H_{VG}=\varepsilon_n + \mathbf A(t)\cdot\boldsymbol\pi$. Two ingredients:

* the bare window span $W = E_{\max}({\rm CB2}) - E_{\min}({\rm VB{-}1})$ —
  for the Si/GaAs gap-edge quartet $W \approx 11\text{–}12$ eV
  ($0.42\text{–}0.45$ Ha);
* the peak velocity-gauge coupling $\Omega = A_{\rm peak}\,\max|\pi_{vc}|$,
  with $A_{\rm peak} = F_{\rm peak}/\omega$. At $F = 10$ MV/cm
  ($1.94\times10^{-3}$ a.u.) and $\hbar\omega = 0.125$ eV:
  $A_{\rm peak} \approx 0.42$ a.u., $|\pi_{vc}|\approx0.3\text{–}0.5$ a.u.
  $\Rightarrow \Omega \approx 0.13\text{–}0.2$ Ha. (The $A^2/2$ shift is
  band-uniform — it cancels in every energy difference and costs nothing.)

The dressed span $\widetilde W = \sqrt{W^2 + 4\Omega^2} \approx 0.55\text{–}0.6$
Ha at 10 MV/cm. The practical 4th-order phase rule:

$$
\boxed{\;dt \lesssim \frac{0.2}{\widetilde W}\ \ \text{(populations/valley
movies)},\qquad dt \lesssim \frac{0.05}{\widetilde W}\ \ \text{(currents/HHG)}\;}
$$

**For 2V+2C Si/GaAs at fields up to 10 MV/cm this gives:**

| target | dt [a.u.] | dt [fs] |
|---|---|---|
| population dynamics, BZ movies, nex | **0.3–0.4** | **0.008–0.010** |
| gauge-invariant currents, σ(ω), HHG | **0.08–0.1** | **0.002–0.0025** |

(The x11 defaults dt = 0.4 a.u./0.0097 fs sit exactly on the population
criterion at moderate field; at the full 10 MV/cm drop to ~0.3 a.u.)

**2. Laser period (never binding for THz/MIR).** $T = 2\pi/\omega$; even at
800 nm ($\hbar\omega=1.55$ eV, $T \approx 110$ a.u.) the population rule above
already gives ≥ 275 points/cycle.

**3. Houston-basis rotation (subordinate, but real).** The dissipators act in
the instantaneous eigenbasis; its rotation rate peaks at swept anticrossings,
$\dot\theta \sim F\,|\pi_{vc}|/\Delta E_{\rm gap}$ — at 10 MV/cm and a 1.4 eV
gap this is ~0.016 rad/a.u.t, so the span rule (≥25× stricter) covers it.

**4. Collision flux — watch the CPTP limiter.** Each ring channel transfers
$f\,(1-e^{-\Gamma dt})$ against step-START populations; many independent
quadruples can pile onto one sink within a long step. The global CPTP limiter
keeps this exactly trace-conserving by scaling the whole transfer field, and
prints once:

```
# ring CPTP limiter engaged: dpop scaled by s = ...
```

`s < 1` is the runtime dt-diagnostic: the collision rates are being
**capped**, not resolved — multiply `dt` by roughly the reported `s` to
restore faithful rates (e.g. the 0.5 fs step that motivated the limiter
reported `s = 0.118`, i.e. a ~8× too-large step for that flux). For faithful
linearized channel rates keep $\Gamma_{\max}\,dt \lesssim 0.1$.

**Strong-field sanity check (basis, not dt):** with only 2 CB unfrozen, verify
the VG basis-sufficiency monitor stays quiet (`grep "basis edge" *_rt.log`
empty, P_top small — wiki/06); above ~10 MV/cm on 4 bands the top-band
population, not the time step, is usually the first thing to fail.
