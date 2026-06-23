# `dft_to_epm` — EPM local form factors from a SALMON DFT band structure

The local-EPM solver in this fork (`theory='epm'`, `src/epm`) needs a small table
of Cohen-Bergstresser-style local form factors `V^S(|G|^2)`, `V^A(|G|^2)`. Those
are tabulated only for a few diamond / zincblende semiconductors (built-in:
`GaAs`, `Si`, `Si_cb`). **For any other crystal there is no table.**

`dft_to_epm.py` removes that limitation: run a normal SALMON DFT band structure
for *your* crystal, and this tool **fits** the EPM local form factors so the EPM
band structure reproduces the DFT one (least squares). The result plugs straight
back into `theory='epm'` via the new `epm_material='file'` reader — no
recompilation, closing a **DFT → EPM → SBE** pipeline for arbitrary materials.

This is the standard *semi-empirical* route (band fitting) to EPM form factors,
and is in the same spirit as the machine-learned pseudopotentials of the vendored
[`external/DeePseudopot`](../../external/DeePseudopot) project
(Lin, Coley-O'Rourke & Rabani, *npj Comput. Mater.* **11**, 381 (2025)).

## Two extraction methods (`--method`)

Both methods fit the **same** EPM plane-wave Hamiltonian to the DFT bands; they
differ in how the local potential is parametrised:

| `--method` | Free parameters | When to use |
| :--- | :--- | :--- |
| `lsq` (default) | the form factors `V^S`/`V^A` at each requested shell, independently | fast; **exact** when the reference is itself an EPM band structure; the classic Cohen-Bergstresser picture. |
| `zunger` | the **DeePseudopot analytic Zunger form** `V(q)=a0(q²−a1)/(a2·exp(a3 q²)−1)`, `a0..a3` per species, then sampled at the shells | a *smooth, physically constrained, transferable* local potential — same functional form DeePseudopot fits to ab-initio bands; regularises sparse/noisy data and extrapolates sensibly to shells that were not independently free. |

`--method zunger` **pulls in the vendored DeePseudopot code as a module**: when
`torch` is installed it calls `external/DeePseudopot/utils/pp_func.py::pot_func`
directly (the run reports `zunger backend: deepseudopot(vendored)`); otherwise it
falls back to a NumPy evaluation of the identical closed form. Install the
upstream stack with `pip install -r external/DeePseudopot/requirements.txt` to use
the real module path.

> Note on accuracy: with only the 3–4 classic shells free, both methods reach the
> same band RMS (a 4-parameter Zunger curve can hit 3 shell values exactly). The
> Zunger advantage shows up as a smooth `V(q)` over *all* `|G|` (so you can emit
> extra shells), physical regularisation, and a direct bridge to the full
> DeePseudopot NN model — DeePseudopot itself is more accurate still because it
> additionally fits a continuous NN `v(q)` plus nonlocal / SOC / strain channels
> and richer targets (gaps, masses, deformation potentials).

## What it does

1. Reads a SALMON DFT band structure:
   * `--format band_dat` — `band.dat` from a `theory='dft_band'` run, **or**
   * `--format salmon_eigen` — `SYSNAME_eigen.data` + `SYSNAME_k.data`
     (the EPM/SBE ground-state format).
   * **Energies in both formats are Hartree (a.u.)** — SALMON ignores `&units`
     when writing them — so the default `--dft-energy-unit ha` is correct.
2. Builds the *same* plane-wave EPM Hamiltonian SALMON uses, with the form
   factors as free parameters:
   * `--cell primitive` mirrors `src/epm/epm_solver.f90` (2-atom FCC primitive
     cell; cutoff on `|G|^2` in **a.u.**, matching `epm_pw_cutoff_ry`).
   * `--cell cubic` mirrors `epm_gaas_reference.py` (conventional 8-atom cubic
     cell with the FCC parity selection rule; cutoff on `|G|^2` in `(2π/a)²`).
   * **The cell must match the DFT run**: cubic DFT → `--cell cubic`; a 2-atom
     primitive DFT (or `theory='epm'` output) → `--cell primitive`.
3. Fits `V^S` (and `V^A` for zincblende, via `--shells-a`) at the requested
   shells, eliminating the rigid DFT↔EPM energy-zero offset analytically.
4. Writes:
   * `<prefix>_epm_formfactors.data` — the form-factor table.
   * `<prefix>_fit_report.txt` — fit quality + factors in Ry and Ha.
   * stdout snippets for both EPMs.

### Consuming the result

* **Python EPM — `epm_gaas_reference.py` (primary, maintained):**
  ```python
  import epm_gaas_reference as epm
  epm.load_form_factor_file('Si_fromDFT_epm_formfactors.data')
  # then run with MATERIAL = 'file'
  ```
  (or set `epm.FORM_FACTOR_FILE = '...'` / paste the printed
  `_FITTED_FORM_FACTORS_RY` dict).
* **Fortran `theory='epm'` (deprecated):** `&epm` with `epm_material='file'`,
  `epm_formfactor_file='Si_fromDFT_epm_formfactors.data'`. The tool also prints a
  ready-to-paste Fortran `case` block for `cb_get_form_factors`.

## Quick start (Silicon, cubic cell)

```bash
SALMON=/path/to/build/salmon

# 1. DFT ground state + band structure (uses the existing low-precision sample)
cd samples/exercise_04_bulkSi_gs
$SALMON < Si_gs.inp                       # -> data_for_restart/
ln -sfn data_for_restart restart
$SALMON < Si_band.inp                     # -> band.dat

# 2. Extract EPM form factors from the DFT bands
python3 ../../tools/dft_to_epm/dft_to_epm.py \
    --dft band.dat --format band_dat --cell cubic \
    --a-lattice-au 10.2626 --cutoff-ry 11.1 \
    --material-name Si_fromDFT --shells-s 3,8,11 \
    --nval 16 --nbands-fit 18 --weight-valence 3.0 \
    --out-prefix Si_fromDFT
# -> Si_fromDFT_epm_formfactors.data, Si_fromDFT_fit_report.txt

# 3. Use them in theory='epm' (then feed the GS files to theory='sbe')
#    &epm: epm_material='file', epm_formfactor_file='Si_fromDFT_epm_formfactors.data'
```

A runnable end-to-end script lives in
`samples/exercise_x4_Si_dft_to_epm/run_dft_to_epm.sh`.

## Adapting to a new crystal

* **Diamond / monoatomic-basis** (Si, Ge, C): leave `--shells-a` empty (`V^A≡0`).
* **Zincblende / two-species** (GaAs, InSb, …): pass `--shells-a 3,4,11`
  (the shells where the antisymmetric structure factor is non-zero) so the
  cation/anion asymmetry is fitted too.
* Choose `--shells-s`/`--shells-a` to match which `|G|^2` shells your `&epm`
  basis actually couples (for the classic set: `3, 4, 8, 11`).
* Set `--a-lattice-au` to the **same** lattice constant you will pass as
  `epm_lattice_constant_au`, and `--nval` to (valence electrons)/2 for the cell.

## Tests

```bash
python3 tools/dft_to_epm/tests/test_recovery.py
```
Synthesises EPM bands with known Si factors and checks the fitter recovers them
to ~1e-14 for both forward models (no SALMON build needed).
