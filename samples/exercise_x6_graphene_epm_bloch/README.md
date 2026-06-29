# Exercise x6 — Monolayer graphene EPM → SBE: 2-coset folding, Dirac carriers

Monolayer graphene in the **rectangular 4-atom cell** (zigzag x = a, armchair
y = √3·a, vacuum along z), a **2-fold (2-coset)** supercell of the 2-atom
hexagonal primitive cell. The ground state comes from the **Python EPM** (a
π-band model: `nelec=4`, `nstate=8`); the cubic Fortran EPM does not do this
cell. Graphene is **gapless** — the folded K point lands on the Dirac point.

**The field must be IN-PLANE.** An out-of-plane (z) field does not couple to
the in-plane π bands and produces exactly zero excitation; this input drives
`epdir = x`.

## Run

```sh
# build:
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"

# step 1: Python EPM ground state -> graphene_k/_eigen/_tm/_unfold/_bandpath.data
python3 epm_graphene.py gs

# step 2: SBE real-time on that ground state (in-plane field)
./build/salmon < graphene_sbe_rt.inp

# plot
python3 plot_sbe_results.py -i . -o plots --snapshots
```

Run `epm_graphene.py gs` where you will run `salmon` (it emits into the working
directory). `sysname='graphene'`, `al`, `num_kgrid=4,4,1`, `nelec=4`,
`nstate=8` must match.

## What to look at

- `plots/nex_k_unfold_real_snap_*.png` — the **2-coset** unfolded CB1 on the
  primitive BZ (real carriers). They localize at the **Dirac K points**
  (K = (2/3, 1/3) reduced), zero near Γ — the unfold places carriers in the
  correct primitive-BZ sectors.
- `plots/sbe_intra_current.png` — intra-band (Houston) vs total current.
- `graphene_bandpath.data` (step 1) — the clean primitive bands; the folded K
  is gapless (verified). Plot with `--only-bands`.
- `graphene_sbe_rt_energy.data` — total trace stays at 4 (CPTP).

## Notes / tuning

- **In-plane only:** `epdir_re1 = 1,0,0` (or `0,1,0`). A z-field gives zero
  carriers (correct physics for a 2D sheet).
- **Real vs Houston:** the real-carrier map (default) has no A²(t) breathing;
  for graphene the Houston map overshoots the post-pulse value by ~300 % at the
  field peak (all virtual). Add `--instantaneous` to compare.
- Graphene has no registry entry for the dissipation channels yet (gapless,
  metal-like); keep `yn_sbe_*` off (the plain dynamics is what is validated).

See [`../../wiki/05_folding_unfolding.md`](../../wiki/05_folding_unfolding.md) §6–7
(2-fold folds, N-coset unfold) and [`../../wiki/04_configuration_examples.md`](../../wiki/04_configuration_examples.md).
