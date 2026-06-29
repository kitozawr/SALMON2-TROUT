# Exercise x5 — CdS (wurtzite) EPM → SBE: 2-coset folding, real-carrier maps

CdS is a **wurtzite** crystal: the SBE runs in the orthorhombic √3×1×1
supercell `al = (a, a√3, c)`, a **2-fold (2-coset)** supercell of the hexagonal
primitive cell (the wurtzite analogue of the cubic 4-fold FCC folding). The
ground state comes from the **Python EPM** — the cubic Fortran EPM is
zincblende/diamond only and does not do wurtzite.

Direct gap at Γ ≈ 2.55 eV (Bergstresser–Cohen 1967: 2.58 eV).

## Run

```sh
# build:
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"

# step 1: Python EPM ground state -> CdS_k/_eigen/_tm/_unfold/_bandpath.data
#         (run from the repo root; emits into the current directory)
python3 epm_wurtzite_cds.py gs

# step 2: SBE real-time on that ground state
./build/salmon < CdS_sbe_rt.inp

# plot (use --lattice wurtzite for the primitive band path)
python3 plot_sbe_results.py -i . -o plots --snapshots --lattice wurtzite
```

`epm_wurtzite_cds.py gs` writes the files into the working directory, so run it
where you will run `salmon` (or copy the `CdS_*.data` files next to the input).
`sysname='CdS'`, `al`, `num_kgrid=4,4,4`, `nstate=nelec=32` must match.

## What to look at

- `plots/nex_k_unfold_real_snap_*.png` — the **2-coset** unfolded CB1 on the
  primitive BZ (real carriers). Verifies the wurtzite folding/unfolding places
  carriers in the correct primitive-BZ sectors.
- `plots/nex_k_real_snap_*.png` — folded supercell real-carrier LCB map.
- `plots/sbe_intra_current.png` — intra-band (Houston) vs total current.
- `CdS_bandpath.data` (from step 1) — the clean primitive hexagonal bands;
  `python3 plot_sbe_results.py -i . -o plots --only-bands --lattice wurtzite`.
- `CdS_sbe_rt_energy.data` — total trace stays at 32 (CPTP).

## Notes / tuning

- **Provenance gates:** CdS forbids the electron-hole carrier-carrier channel
  (`yn_sbe_eeh`) — no cited material constant — and `error stop`s if you enable
  it. **Auger** (`yn_sbe_auger`) is available with the cited CdS coefficient and
  is exactly CPTP, but it is a rare, high-field event (see wiki/01 §13).
- **Field:** ω = 0.06 a.u. (1.63 eV) is below the 2.5 eV gap (multiphoton /
  tunnelling). Raise `I_wcm2_1` or ω for stronger excitation.
- The EPM is a local optical model: it reproduces the direct gap and the
  conduction/upper-valence dispersion but not the deep Cd-4d semicore.

See [`../../wiki/05_folding_unfolding.md`](../../wiki/05_folding_unfolding.md) §6–7
(2-fold folds, N-coset unfold) and [`../../wiki/04_configuration_examples.md`](../../wiki/04_configuration_examples.md).
