# Exercise x3 — Silicon EPM → SBE super-compute mode

A self-contained two-step Silicon run: a local-EPM ground state (no
pseudopotential file, no atomic coordinates) feeding a velocity-gauge SBE
real-time propagation with the optional CPTP dissipation channels
(electron-phonon cooling, carrier-carrier thermalization, impact ionization,
band-gap-renormalized threshold). All on a scalar 4×4×4 grid.

This is the runnable companion to the recipes in
[`../../wiki/04_configuration_examples.md`](../../wiki/04_configuration_examples.md).

## Run

```sh
# build (serial is fine for 4x4x4):
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"

# step 1: EPM ground state -> Si_epm_k.data / _eigen.data / _tm.data (+ unfold map)
./build/salmon < Si_epm_gs.inp

# step 2: SBE real-time on that ground state
./build/salmon < Si_epm_sbe_rt.inp
```

`sysname`, `al`, `num_kgrid`, `nstate`, `nelec` are identical in both inputs —
they must match for step 2 to read step 1's output.

## What to look at

- `Si_epm_sbe_rt.data` — current J(t), vector potential, E-field.
- `Si_epm_sbe_rt_energy.data` — total electronic energy + total trace (should
  stay at 32 to machine precision: the channels are CPTP).
- `Si_epm_sbe_nex.data` — excited carrier density n(t) [cm⁻³].
- `Si_epm_sbe_nex_k.data` / `_unfold.data` — k-resolved conduction populations.
- stdout — the channel banners (printed once) report ħω of each phonon mode,
  N_B, the impact-ionization threshold, etc.

## Tuning

- **Start from the plain run:** set every `yn_sbe_*` to `'n'` to recover the
  bare CF4/Yoshida/Strang dynamics, then switch channels on one at a time.
- **Field strength:** `I_wcm2_1`. Impact ionization and BGR only matter at
  strong fields — push toward 1e13–1e14 W/cm² to see the threshold open.
- **GaAs instead:** rebuild the ground state with `epm_material='GaAs'`,
  `epm_lattice_constant_au=10.68`, `al=10.68`, and use the GaAs recipe
  (`sbe_ii_form='stobbe_quartic'`, `sbe_ii_exponent=4`,
  `sbe_ii_threshold_ev=2.1`, `sbe_eph_nu_sat=1.0d14`).
- **Multiscale (Maxwell-SBE):** see the `theory='maxwell_sbe'` recipe in
  wiki/04; the same `&sbe` channel block applies per macropoint.

All defaults carry primary-source citations in the code and in
[`../../wiki/02_constants.md`](../../wiki/02_constants.md).
