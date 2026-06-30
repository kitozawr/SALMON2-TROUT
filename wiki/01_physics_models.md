# Physics Models & Approximations

> Status: master equation + Σ^HF + k-local II + Kuhn-Zurek ✅ implemented; e-ph (§6), nonlocal II, BGR (§9) 🚧 planned (see [status](00_implementation_status.md)). All citations are to primary sources. Constants in [Constants](02_constants.md); schemes in [Numerical Methods](03_numerical_methods.md).

## 1. The master equation

The density matrix ρ(k,t) (size 2Nb×2Nb spinor, or Nb×Nb scalar) at each k-point obeys

> dρ(k,t)/dt = −i[ H_VG(k,t) + Σ^HF[ρ](k,t) , ρ(k,t) ] + D_KZ[ρ] + D_II[ρ] + D_eph[ρ]

- **H_VG = H₀(k) + A(t)·π** — velocity-gauge band Hamiltonian; π = p + A(t) (+ v_SO for spinor; + v_nl only for a nonlocal pseudopotential, which the local EPM lacks). The scalar A²/2 term is dropped from H_VG (commutes with everything) but **restored** in the impact-ionization kinetic-energy gate (§4). [velocity-gauge SBE: Wismer & Yakovlev, PRB 97, 144302 (2018); Houston basis: Yue & Gaarde, JOSA B 39, 535 (2022)]
- **Σ^HF** — optional Coulomb time-dependent Hartree-Fock exchange self-energy (§5). ✅
- **D_KZ** — Kuhn-Zurek pure dephasing (§2). ✅
- **D_II** — impact ionization (§3–4). ✅ (k-local)
- **D_eph** — electron-phonon scattering (§6), super-compute mode. 🚧

## 2. Kuhn-Zurek / Caldeira-Leggett dephasing (strictly CPTP) ✅

Phenomenological −ρ_nm/T2 dephasing is generally not completely positive. This fork uses a wave-packet dephasing model that is exactly CPTP by construction:
- At each step, diagonalize H_VG(t) → Houston/adiabatic basis U(t) and branch positions X_a(t), propagated by group velocities V_a = (U†πU)_aa + A(t).
- Rotate ρ̃ = U†ρU, then dephase via the **Hadamard/Gaussian (RBF) kernel** ρ̃_ab ← exp[−λ(X_a−X_b)²τ] ρ̃_ab.
- The Gaussian kernel is positive-definite (Schoenberg/Bochner); by the Schur product theorem the Hadamard map is CPTP for **any** τ ≥ 0 — no positivity violations, no ad-hoc clipping.
- Rate λ = k_B T / τ_m. **Γ_aa = 0** — decays coherences only; adiabatic populations conserved (pure dephasing).

[Caldeira-Leggett, Physica A 121, 587 (1983); Zurek, RMP 75, 715 (2003); Schoenberg, Ann. Math. 39, 811 (1938); Kuhn et al., PRB 82, 075204 (2010)]

**Limitation:** constant τ_m. The geometric (X_a−X_b)² separation captures decoherence from wave-packet drift but NOT the intrinsic k-dependence of the e-ph coupling strength; the HHG wavelength-scaling exponent will be on the shallow side.

## 3. Impact ionization as a k-local Lindblad (default fast mode) ✅

A hot conduction electron above threshold ionizes a valence electron, creating an e-h pair. Quartic event A_h = √γ c†_h' c†_c1 c_v1 c_h, closed **k-locally** (no momentum transfer), HF-factorized into two **frozen-rate amplitude-damping channels** in the Houston basis: primary relaxation h→h' (conduction branch nearest ε_h−E_g) and cold-pair creation v1→c1. Partner populations and Pauli blockers clamped to [0,1] → every map exactly CPTP. In the Boltzmann limit the diagonals reproduce W_h = γ ρ̃_hh ρ̃_v1v1 (1−ρ̃_c1c1)(1−ρ̃_h'h'); ionization also destroys participating coherences. [Rosati-Iotti-Dolcini-Rossi, PRB 90, 125140 (2014); Taj & Rossi, PRA 78, 052113 (2008)]

## 4. Impact-ionization rate: material-dependent fit form ✅ (form switch is Part B)

General form: **γ_II(ε_kin) = P (ε_kin − E_th)^a Θ(ε_kin − E_th)**, ε_kin from the field-free CBM (the dropped A²/2 restored here via the Houston identity ε_kin = E_h(k+A) − E_CBM). Θ smoothed by a linear ramp of width σ_E.
- **GaAs — hard, quartic (a=4):** P=2e12 s⁻¹eV⁻⁴, E_th=2.1 eV. [Stobbe-Redmer-Schattke, PRB 49, 4494 (1994), Eq. 11]
- **Silicon — soft, quadratic (a=2):** E_th=1.1 eV (near the gap, not 3/2·Eg). [Keldysh, JETP 21, 1135 (1965); Cartier et al., APL 62, 3339 (1993)]
- **Silicon — full-band option:** a=4.6, E_th=1.15 eV. [Kamakura et al., JAP 75, 3500 (1994)]

**Why Si soft, GaAs hard:** threshold behaviour follows the near-threshold DOS; GaAs's steeper DOS gives a harder threshold. [Stobbe 1994; Sano & Yoshii, PRB 45, 4171 (1992)]

**Fit limitations:** direction-averaged; electron-initiated only; no phonon-assisted ionization or field-induced threshold softening (at MV/cm there is strictly no fixed threshold). [Quade-Schöll-Rossi-Jacoboni, PRB 50, 7398 (1994)]

## 5. Coulomb time-dependent Hartree-Fock (optional) ✅

Σ^HF_nm(k) = −Σ_{q≠k} V(k−q) δρ_nm(q), added to H_VG, reproducing both the renormalized single-particle energies (diagonal) and the renormalized Rabi frequency (off-diagonal), with the (1−f_e−f_h) Pauli factor emerging from the commutator. Uses δρ = ρ − ρ₀ so Σ^HF(t=0)=0 (equilibrium gap unchanged). Hermitian → CPTP-safe. [Golde-Kira-Meier-Koch, Phys. Status Solidi B 248, 863 (2011)]

**Folding-artifact bug and fix (Part E, 🚧):** in the 8-atom cubic cell, bands fold 4-fold. An unrestricted Σ^HF spuriously couples states from different primitive-BZ sectors (Γ- to X-derived), ejecting electrons into the wrong valley even at <100 kV/cm. A translationally invariant Coulomb operator conserves primitive crystal momentum → inter-sublattice Fock coupling is exactly zero; project Σ^HF block-diagonally onto the 4 FCC sublattice sectors. See [Folding](05_folding_unfolding.md). [Popescu & Zunger, PRB 85, 085201 (2012); Ku-Berlijn-Lee, PRL 104, 216401 (2010)]

**Cost:** non-k-local, O(N_k²)/step, one MPI all-gather/step. In super-compute mode the exchange folds into the ring pipeline (§7).

## 6. Electron-phonon scattering as a population-relaxing Lindblad (super-compute) ✅

Unlike Kuhn-Zurek, this channel **relaxes populations** (Γ_aa ≠ 0) and reproduces THz bleaching; replaces (user toggles off) Zurek in super-mode. Built as a **sum of explicit channels with smoothed thresholds**:

> **Note (gate fix):** the dissipative half-step that applies this channel and the carrier-carrier channel (§10) is now entered whenever **any** dissipator is on (`flag_decoh ∨ flag_impact ∨ flag_eph ∨ flag_eeh ∨ flag_auger`). A previous gate checked only `decoh ∨ impact`, so a run with **only** `yn_sbe_eph='y'` (or only `yn_sbe_eeh='y'`) silently did nothing. Both now work standalone (verified scalar GaAs); clean/decoh/impact configs are byte-unchanged.

> ν_total(ε) = ν_intra(ε) + Σ_iv ν_iv(ε) Θ_smooth(ε − E_iv)

Jump operators (emission/absorption): L^em_{k,c→c'} = √ν_em c†_{k−q,c'} c_{k,c}, ν_em ∝ (N_B+1); absorption ν_ab ∝ N_B. Deformation-potential matrix elements. Golden-rule rate 1/τ(E) = (πD²)/(ρVω)(N_B + ½ ∓ ½) DOS(E ∓ ħω). [Jacoboni & Reggiani, RMP 55, 645 (1983)]

- **Silicon (non-polar):** intervalley deformation-potential scattering via six phonons (3 g-type same-axis, 3 f-type orthogonal) + intravalley acoustic. No Fröhlich. [Jacoboni-Reggiani 1983; Canali et al., PRB 15, 3994 (1977)]
- **GaAs (polar):** Fröhlich LO + intervalley + acoustic. The **Fröhlich emission rate must include the asinh/log factor**: W ∝ (1/√E)·asinh(√(E/ħω₀−1)); omitting it under-counts ~20–30% at 2–3 eV. [Fawcett-Boardman-Swain, JPCS 31, 1963 (1970)]

**Collision-rate saturation** (drives Si bleaching): ν(ε) = ν_sat[1 − exp(−(ε/ε₀)^n)], n=2, ε₀~0.8 eV. **Smooth form, never a hard min() cutoff.** [Meng et al., PRB 91, 075201 (2015); Fischetti & Laux, PRB 38, 9721 (1988)]

## 7. THz bleaching mechanisms differ by material (same code, different dominant channel) 🚧
- **Silicon:** no Gunn effect; bleaching is pure ν(ε) collision-rate saturation (σ = ne²/(m*ν)). [Chefonov et al., PRB 98, 165206 (2018)]
- **GaAs:** Gunn / Ridley-Watkins-Hilsum intervalley transfer Γ→L (separation 0.29 eV); mobility drops as electrons move to heavy L valleys. ν(ε) saturation secondary. [Gunn, SSC 1, 88 (1963); Ridley & Watkins, PPS 78, 293 (1961); Hilsum, Proc. IRE 50, 185 (1962)]

## 8. T1/T2 coherence-decay bookkeeping (critical correctness point) 🚧
For the population-relaxing e-ph channel, coherence ρ̃_if decays at
> Γ₂(i,f) = [Kuhn-Zurek pure dephasing, FULL strength] + ½(ν_i + ν_f)

The **½ applies ONLY to the population-relaxation (T1) contribution** (from −½{L†L,ρ}). The Hermitian pure-dephasing Kuhn-Zurek channel enters at **FULL strength** (rate from (X_a−X_b)², no extra ½). **Do NOT apply a blanket γ_deph=½ν to all channels.** Unit test: a Hermitian-only dissipator conserves all populations (Γ_aa=0). [Breuer & Petruccione, OUP 2002]

**Occupation-max Pauli normalization (scalar vs spinor).** The diagonal density-matrix entries run in **[0, occ_max]** with occ_max = 2 for scalar bands (`yn_sbe_spinor='n'`, two electrons/band) or 1 for spinor bands. **Every** population-changing channel — e-ph (§6), impact ionization (§3), carrier-carrier (§10), Auger (§12) — therefore builds its Pauli factors as the **fractional** occupation `ρ/occ_max` (presence) and `1 − ρ/occ_max` (blocking), each clamped to [0,1] (`sbe%occ_max = merge(1,2, spinor)`). A bare `1−ρ` would wrongly clamp a half-full scalar band (ρ=1.5) to a zero blocking factor; `1−ρ/2 = 0.25` correctly keeps the room. The **off-diagonal (coherence) damping is occupation-independent** (the GKLS √(e^{−Γτ}) factor / the `(1−α)` convex weight), as it must be — coherence damping is not Pauli-blocked.

## 9. Bandgap renormalization coupling to the II threshold 🚧
E_th is tied to the gap; Σ^HF shrinks the gap with density. Make E_th density-dependent **only above n > 5e18 cm⁻³**:
> E_th(t) = E_th0 − |ΔE_BGR(n(t))|, ΔE_gap[eV] = −1.9e-8 (n[cm⁻³])^(1/3)

−19 meV at 1e18, −41 meV at 1e19 (~2–5% of a ~1.1 eV Si threshold). [Vashishta & Kalia, PRB 25, 6492 (1982)]. Distinct from dopant BGN (√n) [Lanyon & Tuft, IEEE TED ED-26, 1014 (1979)]. K carries a factor-~2 ambiguity; treat [1.9, 3.8]e-8 eV·cm as tunable.

## 10. Carrier-carrier (e-e / e-h) scattering — Part F (channel ✅, screening G ✅)
A **second-Born / GW statically-(or dynamically-)screened-Coulomb collision integral cast into CPTP Lindblad form**. The diagonal (population) part is the in−out screened-Coulomb collision integral with (1−f) Pauli factors and a **direct−exchange** matrix element

> Γ_cc(k₁) = (2π/ħ) Σ_{k₂,q,λ} |W̃(q)|² { f_{k₂}(1−f_{k₁+q})(1−f_{k₂−q}) − (1−f_{k₂}) f_{k₁+q} f_{k₂−q} } δ(E_{k₁}+E_{k₂}−E_{k₁+q}−E_{k₂−q}),  |W̃|² = W(W*_dir − W*_exch),

and the off-diagonal part is the **excitation-induced dephasing (EID)** of the polarization (γ = γ₀ + a·n). It becomes CPTP via density-fluctuation jump operators L_q = Σ_k c†_{k+q}c_k under the **Taj-Rossi** completely-positive Markov limit [PRA 78, 052113 (2008)] with the **Rosati-Iotti-Dolcini-Rossi** nonlinear single-particle closure [PRB 90, 125140 (2014)] (carries the (1−ρ) Pauli factors; preserves Tr ρ=N and 0≤ρ≤1).

**Conservation (validation invariants):** conserves total carrier **number** Σf_k AND **energy** ΣE_k f_k within the carrier subsystem (thermalizes to a hot Fermi-Dirac), but does **not** relax energy to the lattice (that is e-ph, §6).

**No double-counting with HF (§5):** HF is the first-order coherent self-energy (energy/Rabi renormalization, no scattering); carrier-carrier is the correlation (2nd-Born/GW) self-energy (real scattering + dephasing). Add carrier-carrier ONLY as the dissipative channel — do not also add a static screened-exchange shift from the same W. [Baym Φ-derivable; Σ=Σ^HF+Σ^(2B)]

**Computation:** O(N_k²) double sum over (k₂,q) → rides the **same ring pass** (Part D) as Σ^HF + nonlocal II + nonlocal e-ph (one extra accumulator, no new communication); ε(q) computed once/step from the gathered ρ; broadened-delta bins; Houston basis; predictor-corrector. **Rate scale 1e13–1e14 s⁻¹** at n=1e17–1e19 cm⁻³ (~10–200 fs). [Goodnick-Lugli PRB 37, 2578; Fischetti-Laux PRB 38, 9721; Mocatti et al. arXiv:2512.08618 (2025) Eq.55; EID: Honold PRB 40, 6442 (1989), Wang PRL 71, 1261 (1993)]

## 11. Dielectric screening — Part G (✅ primitives implemented)
W(q)=V(q)/ε(q[,ω]); three selectable models (pure functions in `sbe_superres_ssbe`, unit-tested):
- **(a) Thomas-Fermi/Debye:** ε(q)=1+κ²/q², degenerate κ_TF²=4(3n/π)^⅓/ε or nondegenerate κ_D²=4πn/(εk_BT). [Ashcroft-Mermin]
- **(b) static Lindhard/RPA — default:** ε(q,0)=1+(κ_TF²/q²)F(q/2k_F), F(x)=½+(1−x²)/(4x)ln|(1+x)/(1−x)| (2k_F kink → Friedel). TF over-screens for q>2k_F. [Lindhard 1954]
- **(c) dynamic Lindhard / LOPC — GaAs only, n≳5e17:** ω_{L±}²=½[(ω_p²+ω_LO²)±√((ω_p²+ω_LO²)²−4ω_p²ω_TO²)], ω_p²=4πn/(ε_∞m*), single-plasmon-pole. **Disabled for Si** (non-polar). [Varga PR 137, A1896; Mooradian-McWhorter PR 177, 1231 (1969)]

Static screening under-estimates the rate; dynamic (LOPC) is needed for sub-100-fs thermalization (Elsaesser/Shah).

## 12. Auger recombination as a number-conserving CPTP Lindblad ✅ channel (k-local); ⛔ no cited C yet
The inverse of impact ionization (§3): instead of a hot carrier creating an e–h pair, two carriers + a hole give one hot carrier and a destroyed pair. It is a **3-body, density-gated** process whose recombination rate goes as **R = C·n³** (the per-carrier rate is γ = C·n²), with a per-material coefficient C and an activation density n_gate.

> **Provenance note (2026-06-30):** the value formerly used as the CdS default, "C = 2.0×10⁻³⁰ cm⁶/s [Haury et al., PRB 57, 11513 (1998)]", was a **fabricated citation** (the real Haury et al. paper is PRL 79, 511 (1997) on CdMnTe ferromagnetism, unrelated to Auger in CdS; the coefficient is unconfirmed). It has been **removed** and CdS Auger is gated off. **No material currently ships a verified C**, so the channel requires an explicit `sbe_auger_c_cm6s`. Cited per-material coefficients (GaAs / Si / graphene) and the full **nonlocal** rebuild are the subject of `wiki/07_nonlocal_auger.md`.

**Gap-edge mean-field (HF-factorized) closure**, in the Houston/adiabatic basis on the gap-edge branches (top valence iv1, lowest conduction ic1, hot target ic_hot energy-matched to E(ic1)+E_g):
- **recombination** ic1 → iv1: a conduction electron fills a valence hole, destroying an e–h pair;
- **promotion** ic1 → ic_hot: the released gap energy E_g lifts a second conduction electron to the hot state.

Both are realized as the **exact finite-time amplitude-damping GKLS map** `amp_damp_channel` (the same CPTP primitive as §6 e-ph), at the rate γ carrying the occ_max-normalized, [0,1]-clamped Pauli factors (CB electron present `f_c/occ`, VB hole present `1−f_v/occ`, hot target empty `1−f_hot/occ`). Because `amp_damp` is trace-preserving, the **total carrier number is conserved** (Auger rearranges, it does not remove electrons); energy is conserved to the mean-field order (same as §3 impact ionization). Recombination self-limits — as the holes fill, `1−f_v/occ → 0` and the rate vanishes.

**It is a RARE channel.** Auger acts on the **real (Houston/adiabatic) populations**, not the virtual driving polarization, and C is typically tiny — so at ordinary fields the dynamics are essentially unchanged with the channel on vs off; it only becomes visible at very high real carrier density / strong fields. Its role is to be present and exactly CPTP, density-gated, with the n³ law. **Provenance-gated:** only materials with a verified C may enable it — currently **none** ship a default, so `yn_sbe_auger='y'` `error stop`s unless `sbe_auger_c_cm6s` is set. The graphene gapless carrier-multiplication variant (nearly thresholdless) is a separate, model-dependent channel. Both — cited per-material 3D coefficients and the graphene 2D branch — are folded into the **nonlocal-Auger** rebuild (`wiki/07_nonlocal_auger.md`), which reuses the nonlocal impact-ionization kernel by detailed balance. [GKLS: Lindblad 1976; Taj-Rossi PRA 78, 052113 (2008). CPTP invariants: `tests/test_auger_cptp.f90`.]
