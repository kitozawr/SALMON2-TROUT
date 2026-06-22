# Band Folding & Unfolding

> Why the cubic-cell representation folds the primitive bands, why the folding is exact, how it is unfolded, and why it matters for Hartree-Fock. Applies identically to GaAs (zincblende) and Silicon (diamond) — both FCC. Folding/unfold pipeline ✅ implemented; HF sublattice projection (§4) 🚧 Part E.

## 1. The folding ✅
The conventional cubic 8-atom cell is a supercell of 4 primitive FCC cells. The MP-grid band plot shows the primitive bands **folded 4-fold**: every cubic k-point carries the states of 4 primitive BZ points, and the conduction manifold appears as 4 overlaid copies of CB1/CB2/CB3. These dense crossings are an artifact of the supercell representation, not physics. [Popescu & Zunger, PRB 85, 085201 (2012)]

## 2. The folding is exact (parity selection rule) ✅
The cubic reciprocal lattice is an index-4 sublattice of the FCC reciprocal lattice; the 4 primitive k-points folding to one cubic k are distinguished by the **4 cosets** of the cubic reciprocal-lattice subgroup. The Hamiltonian is **block-diagonal over these 4 FCC reciprocal sublattices to machine precision** (a parity/translation selection rule), so the clean primitive picture is recovered exactly. The code asserts this block-diagonality at runtime.

## 3. Unfolding pipeline (does not require regenerating the ground state) ✅
Three stages — EPM → SALMON → plotter — reusing the ground state as-is:
1. **EPM builds the unfold map (once, cheap).** Re-diagonalize the MP-grid cubic Hamiltonians, assign every cubic band to the 4 FCC sublattices, write spectral weights w_s = |⟨ψ|P_s|ψ⟩|² (Σ_s w_s = 1) and the energy-ranked primitive-band index to `SYSNAME_unfold.data`. Must be on the same k-grid as the ground state. The GS dataset is NOT touched.
2. **SALMON runs the dynamics.** With `SYSNAME_unfold.data` present, the SBE writes `SYSNAME_sbe_nex_k_unfold.data`: the crystal-gauge population of the physical levels (VB-1, VB, CB1, CB2; spins summed) at each primitive point k_prim = k_sc + G₀(s), distributed over sublattices by the spectral weights (not a hard argmax), so a symmetry degeneracy splits equally and the result is symmetric by construction.
3. **Plotter visualizes** unfolded (primitive FCC BZ) and folded (cubic) k-maps, plus optional A(k,E) spectral views.

## 4. Why folding matters for Hartree-Fock (the bug and the fix) 🚧 Part E
A Fock exchange Σ^HF that couples all bands at a given cubic k will **spuriously couple states belonging to different primitive-BZ sectors** (e.g. a Γ-derived state to an X-derived state). Physically this ejects an electron into the wrong valley even in weak fields (<100 kV/cm) — an artifact, because:

> A translationally invariant two-body operator (the Coulomb interaction) conserves primitive crystal momentum. Two states folding to the same cubic k but belonging to different primitive k (different cosets) have different primitive crystal momenta, so the exact Coulomb/Fock matrix element between them is ZERO. Any nonzero inter-sublattice Fock coupling is purely a folding artifact.

**Fix (sublattice-block projection):** project Σ^HF block-diagonally onto the 4 FCC sublattice sectors — keep intra-sublattice matrix elements, set inter-sublattice ones to zero. Equivalent: unfold ρ to primitive k (via the w_s weights), apply Σ^HF with primitive momentum conservation, fold back. Consistent with the exact block-diagonality of §2. Controlled by `yn_sbe_hf_sublattice_proj` (default-on whenever HF is on in a folded cubic cell). [Popescu & Zunger, PRB 85, 085201 (2012); Ku-Berlijn-Lee, PRL 104, 216401 (2010)]

**Validation:** with the projection ON, a weak-field run must show zero spurious Γ→X/Γ→L population transfer attributable to exchange. If inter-valley population appears in weak fields with HF on but not off, the projection is mis-applied.

## 5. Silicon vs GaAs ✅
The folding/unfolding is **identical** for both — diamond and zincblende share the FCC Bravais lattice. The GaAs unfold machinery transfers to Si verbatim; only the EPM form-factor table changes (and V^A=0 for Si). High-symmetry valley positions differ (Si: 6 Δ-valleys near X at 0.85·2π/a along ⟨100⟩; GaAs: 4 L-valleys along ⟨111⟩), but the unfolding code finds sublattice character by spectral weight, not by hardcoded valley coordinates — no material-specific change needed.
