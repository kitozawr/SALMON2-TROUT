# exercise_x9 — bulk Si: REAL DFT levels into the SBE (instead of the EPM)

The `_k.data / _eigen.data / _tm.data` files the SBE reads were originally
modeled on SALMON's DFT output — this exercise closes that loop again: a
**real (deliberately rough) DFT calculation** of Si on the **FCC primitive
cell** feeds the velocity-gauge SBE directly, plus a **real DFT band path**
from `theory='dft_band'`. No empirical pseudopotential anywhere.

What it demonstrates:

* the exact file contract DFT → SBE (reduced k + `# b1/b2/b3` header for the
  non-orthogonal cell; `esp[eV]` unit tag auto-converted to Hartree by the SBE
  reader; both `_tm.data` blocks — with DFT the nonlocal `<u|-i[r,Vnl]|u>`
  block is genuinely nonzero, the EPM writes zeros there),
* the **primitive cell spelled with explicit `al_vec1..3`** (a/2 off-diagonal
  rows) — NOT the cubic 8-atom `al(1:3)=5.43` cell of the original
  exercise_04; the SBE cell must be the primitive one,
* the **active/frozen-core window**: DFT drags the deep bonding s-band along
  (≈ −10.5 eV below E_F at the first k-point); `frozen_core_threshold_ev`
  freezes it out of the dynamics (11 of 12 bands active),
* `theory='dft_band'` now emits `SYSNAME_bandpath.data` in the plotter/SBE
  spectral contract (same format as the Python EPM band paths).

**Deliberately rough** (proof of concept): 4×4×4 k-grid, 12³ real-space grid,
SCF threshold 1e-6. The point is the correct *transfer* of k/eigen/tm and the
band path, not converged levels. The rough LDA answer is still recognizably
silicon: indirect gap ≈ 0.6 eV with the CBM on the Δ line near X (LDA
underestimates the experimental 1.17 eV), Γ-conduction ≈ 2.5 eV.

## Run (serial; a few minutes total)

```bash
cd samples/exercise_x9_bulkSi_dft_sbe

# 1/3 DFT ground state -> Si_k.data, Si_eigen.data, Si_tm.data, data_for_restart/
../../build/salmon < Si_dft_gs.inp | tee dft_gs.out

# 2/3 DFT band path (restarts from the GS density; outputs go to ./band/ so the
#     path-k Si_k.data/Si_eigen.data do NOT overwrite the MP-grid files!)
../../build/salmon < Si_dft_band.inp | tee dft_band.out
cp band/Si_bandpath.data .          # co-locate with the SBE outputs for plotting

# 3/3 short SBE run on the DFT levels (60 fs, A||z pulse from 10.txt)
../../build/salmon < Si_sbe_rt.inp | tee sbe_rt.out

# plots: band path, band structure, A(k,E) spectral frames, k-t maps, currents
python3 ../../plot_sbe_results.py -i . -o plots --spectral --valleys
```

> **DFT runs carry no `# material =` header.** That header (and the band-path
> auto-detection it drives, `--lattice auto`) is written only by the EPM. A DFT
> `_k.data` has none, so `--lattice auto` falls back to **fcc** — which is
> correct here (bulk Si is diamond/fcc, path `L-Γ-X-W-K`) and for any cubic
> material. You do **not** add `&epm epm_material=…` for a DFT run (that block is
> not read in `theory='dft'`). For a *non-cubic* DFT material, just pass the path
> to the plotter explicitly, e.g. `--lattice wurtzite` (or `--band-path A Γ M K
> Γ`). The `--spectral` band path above is unaffected either way: it reads the
> exact nodes the DFT `theory='dft_band'` step wrote into `Si_bandpath.data`.

## What to check

* **SBE startup diagnostics** (`sbe_rt.out`):
  `# read_eigen_data: esp[eV] header detected -> converting to Hartree`,
  `cell = PRIMITIVE`, and the Frozen Core Check table — band 1
  (E−E_F ≈ −10.5 eV) must be `active = F`, bands 2–12 active (11/12).
* **Electron number**: the `electrons` column of the run log / `Si_sbe_rt.data`
  stays at 8.000 for all steps (CPTP trace with the frozen core on).
* **Excitation**: `Si_sbe_nex.data` reaches ~2e22 cm⁻³ by 60 fs with `10.txt`
  (real carriers; use `1.txt` for a 10× weaker, near-perturbative drive).
* **Band path** (`plots/bandpath_Si_bandpath.png`): VBM at Γ, CBM on Γ–X near
  X (indirect), rough-LDA gap ≈ 0.6 eV. `plots/band_dat_band.png` is the same
  data via the legacy `band.dat`.
* **Spectral frames** (`plots/spectral_frames/`): the 4 gap-edge bands on the
  DFT skeleton; carriers pool near the Δ/X conduction minimum as the pulse
  ramps.

## Field files

`A/eV/fs` units: column 1 = time [fs] (never rescale it), columns 2–4 = vector
potential Ax,Ay,Az. `10.txt` (A_max = 0.28, along z) and `1.txt` (A_max =
0.028) are the same pulse shape at different amplitudes — see exercise_x8 for
the full amplitude-sweep / basis-sufficiency methodology.

## Notes & pitfalls

* `nelec/nstate/num_kgrid/al_vec` **must match** between the DFT and SBE
  inputs (they size and interpret the GS files).
* The **non-orthogonal cell forbids r-space domain decomposition**: run
  serially or give MPI ranks that fit k/orbital parallelization (the automatic
  distribution now does this by itself; explicit `nproc_rgrid > 1` aborts).
* The band step writes its own `Si_k.data`/`Si_eigen.data` (path k-points,
  weight 0!) — that is why `Si_dft_band.inp` sets `base_directory='./band'`.
  Never let a `dft_band` run share a directory namespace with the GS files the
  SBE reads.
* `Si_rps.dat` is the same norm-conserving Si pseudopotential as exercise_04
  (`izatom=14`, `lloc_ps=2`).
* Tighten `num_rgrid`, `num_kgrid`, `threshold` (and raise `nstate`) for
  production physics; re-verify basis sufficiency (wiki/06, exercise_x8) after
  any change.
