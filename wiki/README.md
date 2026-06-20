# SALMON2-SBE Wiki

Long-term reference and **persistent project memory** for the SALMON2-SBE fork. This directory is the single source of truth that survives context resets: on every "continue", read [`00_implementation_status.md`](00_implementation_status.md) first to recover exactly where the work stands.

A fork of [SALMON](http://salmon-tddft.jp/) extending the Semiconductor Bloch Equations (SBE) module with a CF4/Suzuki-Yoshida exponential propagator, strictly-CPTP Lindblad dissipation (Kuhn-Zurek dephasing, impact ionization, electron-phonon scattering), a self-contained local Empirical Pseudopotential (EPM) ground-state solver for GaAs and Silicon, and an optional nonlocal "super-compute" mode with a ring-pipeline MPI dissipator.

## Pages
0. **[Implementation Status](00_implementation_status.md)** — live progress tracker, decisions log, test inventory, next steps. **Read first on resume.**
1. **[Physics Models & Approximations](01_physics_models.md)** — master equation, dissipation channels, assumptions. Cited.
2. **[Constants & Coefficients](02_constants.md)** — every default value with its primary source.
3. **[Numerical Methods](03_numerical_methods.md)** — CF4 Magnus, Yoshida, Strang, Houston basis, energy-bin search, predictor-corrector, ring-pipeline MPI, sublattice-block HF projection, CPTP proofs.
4. **[Configuration & Examples](04_configuration_examples.md)** — all namelist parameters, runnable examples, pipelines.
5. **[Band Folding & Unfolding](05_folding_unfolding.md)** — cubic-cell 4-fold FCC folding, exact-folding proof, unfold/refold pipeline, why it matters for Hartree-Fock.

## Conventions
- **Units:** EPM structural block in Rydberg (kinetic |k+G|², lengths in Bohr); SBE dynamics in Hartree atomic units (ħ=mₑ=|e|=1). 1 Ry = ½ Hartree.
- **Basis:** dissipative channels operate in the instantaneous **Houston/adiabatic basis** U(t) (eigenbasis of H_VG(t)), NOT the field-free Bloch basis. Adiabatic populations ρ̃ = U†ρU are physical; Bloch populations are virtual during the pulse in the velocity gauge.
- **CPTP:** every dissipator is a genuine GKLS Lindblad generator; finite-step maps exp(τD) are CP and trace-preserving for any τ ≥ 0.
- **Provenance:** every default constant carries a primary-source citation in **both** the code comment and the wiki. Do not introduce a number without a source.

## Working conventions for this project (agreed with maintainer)
- Develop **step by step**, one bounded increment per session; never cram.
- **Test grid:** 4×4×4 (or 5×5×5), **scalar (no spinor)** — sufficient for all validation here.
- Every increment: **working code + clean code + detailed docs (this wiki) + a test**. Commit signed-off, push.
- All new capability behind flags defaulting OFF (or Si-only) so existing GaAs runs stay byte-for-byte unchanged.
- Tests live in [`../tests/`](../tests/); each is self-contained and documented.
