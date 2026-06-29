# Exercise x4 — GaAs EPM → SBE: real-carrier maps & intra-band current

The GaAs reference run (the 100%-validated material). A self-contained
two-step run: a local-EPM zincblende ground state feeding a velocity-gauge SBE
real-time propagation, with **no dissipation** (plain CF4/Yoshida/Strang). The
focus is the per-k carrier visualisation:

- **`*_sbe_nex_k_real.data`** — the *real carriers only* (fixed-basis diabatic
  occupation, the k-resolved n_ex). No reversible A²(t) virtual-polarization
  breathing: it accumulates and **freezes once the field passes**. This is the
  map to read.
- **`*_sbe_nex_k.data`** — the instantaneous Houston-basis population. Physical
  *during* the pulse but carries the virtual breathing (∝A(t)²); equals the
  real map after the pulse. Kept for diagnostics.
- **`*_sbe_intra_current.data`** — the intra-band (drift) current in the
  Houston basis. In the velocity gauge only the total current is gauge
  invariant; the intra/inter split is physical in the Houston basis, and
  J_intra vanishes once the field is off.

## Run

```sh
# build (serial is fine for 4x4x4):
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"

# step 1: EPM ground state -> GaAs_epm_k.data / _eigen.data / _tm.data
./build/salmon < GaAs_epm_gs.inp

# step 2: SBE real-time on that ground state
./build/salmon < GaAs_epm_sbe_rt.inp

# plot: real-carrier maps (default), intra-band current, conductivity, ...
python3 plot_sbe_results.py -i . -o plots --snapshots
```

`sysname`, `al`, `num_kgrid`, `nstate`, `nelec` are identical in both inputs —
they must match for step 2 to read step 1's output.

## What to look at

- `plots/nex_k_real_snap_*.png` — the **real-carrier** snapshots. After the
  pulse (last frames) Γ is populated and the L-valley folds at the cube corner
  — which dominate *mid-pulse as virtual polarization* — are the **least**
  populated. The folded LCB band mixes valleys, so read the field-containing
  planes (kx–kz / ky–kz).
- `plots/sbe_intra_current.png` — J_intra (drift, vanishes after the pulse) vs
  J_total (gauge-invariant) and J_inter = J_total − J_intra (interband
  polarization, keeps oscillating).
- `GaAs_epm_sbe_rt_energy.data` — total trace stays at 32 to machine precision.
- Add `--instantaneous` to `plot_sbe_results.py` to also see the breathing
  Houston maps for comparison.

## Tuning

- **Field strength / colour:** `I_wcm2_1`. Below-gap (ω = 0.06 a.u. ≈ 1.63 eV)
  drives multiphoton/tunnelling excitation spread over k; raise ω toward the
  2.5 eV gap for a more Γ-localised response.
- **Dissipation:** every `yn_sbe_*` defaults OFF. Turn them on one at a time
  (see `exercise_x3` and wiki/04) to add e-ph cooling, carrier-carrier
  thermalization, impact ionization, etc. — the channels are CPTP, trace stays
  at 32.
- **Spinor (spin-orbit) GaAs:** the Python reference `epm_gaas_reference.py`
  emits the spin-orbit GS (`GaAs_cubic_so_*`); set `yn_sbe_spinor='y'` and the
  maps then carry occupation 1 per band.
- **Unfolded (per-coset) maps:** the cubic Fortran EPM emits only the folded
  GS, so this exercise shows the *folded* real-carrier map. For the per-coset
  UNFOLDED primitive-BZ maps, generate the GS + 4-fold FCC unfold map with the
  Python EPM (`python3 epm_gaas_reference.py unfoldmap`) — the 2-coset
  unfolding is demonstrated end-to-end in exercises **x5** (CdS) and **x6**
  (graphene).

See [`../../wiki/04_configuration_examples.md`](../../wiki/04_configuration_examples.md)
for the per-channel reference and [`../../wiki/05_folding_unfolding.md`](../../wiki/05_folding_unfolding.md)
for the real-vs-Houston population discussion.
