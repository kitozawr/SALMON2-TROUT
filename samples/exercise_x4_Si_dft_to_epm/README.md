# Exercise x4 — DFT → EPM form factors → EPM ground state (Silicon)

This sample turns a **SALMON DFT** band structure into **EPM-compatible local
form factors** and then runs the local-EPM ground state from them — the route to
use `theory='epm'` (and the downstream SBE pipeline) for crystals that have **no
tabulated Cohen-Bergstresser form factors**.

It reuses the existing low-precision cubic Si DFT input in
`samples/exercise_04_bulkSi_gs/` (LDA, `num_rgrid=12³`, `4×4×4` k-points), so no
new pseudopotential or high-accuracy run is required — accuracy of the *DFT* step
is deliberately modest; the point is the extraction workflow.

## Files

| File | Purpose |
|------|---------|
| `run_dft_to_epm.sh`   | Runs DFT → `dft_band` → `tools/dft_to_epm/dft_to_epm.py`. Writes `work/Si_fromDFT_epm_formfactors.data`. |
| `Si_epm_fromdft.inp`  | `theory='epm'` input that reads those fitted factors via `epm_material='file'`. |

## Run it

```bash
SALMON=/path/to/build/salmon

# 1) DFT -> band.dat -> fitted form factors  (in ./work)
SALMON="$SALMON" bash run_dft_to_epm.sh

# 2) build the EPM ground state from the fitted factors
cp work/Si_fromDFT_epm_formfactors.data .
"$SALMON" < Si_epm_fromdft.inp        # -> Si_epm_fromdft_{k,eigen,tm}.data
```

## What to expect

The fit reproduces the DFT valence + low-conduction bands with a band RMS of
~0.2 eV (this low-precision LDA run), giving local form factors of the right
order and sign, e.g.

```
shell |G|^2   V^S [Ry]
   3        -0.14
   8         0.04
  11         0.15
```

compare to the built-in Kunikiyo Si set (`-0.2258, 0.05698, 0.070709 Ry`). They
differ because (i) the DFT here is coarse and (ii) DFT's *nonlocal* pseudopotential
is being projected onto a purely *local* EPM — exactly the approximation EPM makes.
Improve the agreement by converging the DFT (`num_rgrid`, k-points, `nref_band`),
widening the fit window, or adding shells.

The resulting `Si_epm_fromdft_{k,eigen,tm}.data` are in the standard EPM→SBE
format and can be fed straight into a `theory='sbe'` run (see the main README's
*Minimal EPM → SBE Pipeline*), with `sysname='Si_epm_fromdft'` and matching
`al`, `num_kgrid`, `nstate`, `nelec`.

See `tools/dft_to_epm/README.md` for all options (zincblende `V^A`, primitive vs
cubic cell, choosing shells, etc.).
