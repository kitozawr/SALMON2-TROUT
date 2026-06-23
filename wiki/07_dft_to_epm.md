# DFT → EPM: extracting local form factors for any crystal

> Maintained EPM = the **Python** `epm_gaas_reference.py`. The Fortran
> `theory='epm'` is deprecated (still works for back-compat).
> Tooling: [`tools/dft_to_epm/`](../tools/dft_to_epm); vendored reference:
> [`external/DeePseudopot/`](../external/DeePseudopot).

## 1. The problem

The local EPM Hamiltonian couples plane waves through a handful of **local form
factors** keyed by the reciprocal-shell index `|G|²` (in `(2π/a)²` units):

```
H_{G,G'}(k) = ½|k+G|² δ_{G,G'} + V^S(|G-G'|²) cos(ΔG·τ) + i V^A(|G-G'|²) sin(ΔG·τ)
```

Cohen-Bergstresser / Kunikiyo tabulated `V^S, V^A` only for a few diamond /
zincblende semiconductors (built in: `GaAs`, `Si`, `Si_cb`). For **any other
crystal there is no table** — but SALMON's DFT (`theory='dft'`) still produces a
Kohn-Sham band structure. `tools/dft_to_epm` turns that band structure into
EPM-compatible form factors by **fitting** them so the EPM bands reproduce the
DFT bands (least squares). This is the textbook definition of *empirical* form
factors, and is the same family of method as the vendored **DeePseudopot** project
(which fits its pseudopotential to ab-initio bands read from BerkeleyGW/QE).

## 2. Two fitting methods (`--method`)

Both share the **same** EPM forward model (a NumPy mirror of the SALMON
Hamiltonian — primitive `epm_solver.f90` *or* cubic `epm_gaas_reference.py`) and
both eliminate the rigid DFT↔EPM energy-zero offset analytically. They differ in
how the local potential is parametrised:

| `--method` | Free parameters | Character |
| :--- | :--- | :--- |
| `lsq` (default) | the form factors at each requested shell, independently (`N_shell` numbers) | unconstrained per-shell fit; the classic Cohen-Bergstresser picture |
| `zunger` | the analytic Wang-Zunger local form `V(q)=a0(q²−a1)/(a2·exp(a3 q²)−1)`, `a0..a3` per species, then sampled at the shells | smooth, physically constrained, transferable `V(q)` over **all** `|G|` |

`--method zunger` **pulls the vendored DeePseudopot code in as a module**: with
`torch` installed it calls `external/DeePseudopot/utils/pp_func.py::pot_func`
directly (the run prints `backend: deepseudopot(vendored)`); otherwise a NumPy
fallback evaluates the identical closed form. Install upstream deps with
`pip install -r external/DeePseudopot/requirements.txt`.

## 3. So — is `zunger` *more accurate for the discrete shells*?

Short answer: **not meaningfully for the minimal 3–4-shell case with good data;
it becomes more robust/accurate as you add shells or as the data degrades.** This
was measured directly — generate Si EPM bands from the *known* Kunikiyo factors,
add band noise / restrict k-points, fit both ways, and compare the recovered shell
values to the truth (mean `|V_fit − V_true|`, in milli-Rydberg, averaged over
seeds):

| scenario | shells | lsq mean err | zunger mean err |
| :--- | :--- | ---: | ---: |
| clean, nk=12 | 3, 8, 11 | **0.00 mRy** | **0.00 mRy** |
| 50 meV noise, nk=12 | 3, 8, 11 | 3.98 | 3.96 |
| 100 meV noise, nk=8 | 3, 8, 11 | 7.13 | 6.84 |
| 200 meV noise, nk=6 | 3, 8, 11 | 22.37 | 21.31 |
| 50 meV noise, nk=6, **over-parametrised** | 3, 8, 11, **16** (true V₁₆=0) | 15.86 | **14.15** |

Reading of the table:

* **Clean, sufficient data → identical.** `lsq` is the unconstrained optimum, and
  a 4-parameter Zunger curve has enough freedom to hit 3 shell values exactly, so
  both land on the truth and give the same band RMS. (This is also why on the real
  low-precision cubic Si DFT both modes return the *same* `V^S ≈ −0.14, 0.04,
  0.15 Ry` at RMS ≈ 0.21 eV.)
* **Noisy data → `zunger` marginally better** (a few %): the smooth physical form
  regularises against noise.
* **Over-parametrised → `zunger` clearly better.** Asking for shell 16 (whose true
  value is 0), `lsq` invents a spurious `V₁₆` (35 mRy error) by overfitting noise,
  while the Zunger form keeps it small (26 mRy) — and its low-shell values stay
  sane. The more shells you request relative to the information in the bands, the
  more `lsq` overfits and the more `zunger` wins.

So the honest picture: **the discrete-shell representation itself caps accuracy**,
and within it `lsq` and `zunger` are equivalent when well-determined. `zunger`'s
real value is (i) robustness when the fit is under-determined or noisy, (ii) a
smooth `V(q)` you can sample at *any*/extra shells, and (iii) a direct bridge to
the **full** DeePseudopot model — which is more accurate still because it fits a
continuous NN `v(q)` plus nonlocal / SOC / strain channels against richer targets
(gaps, masses, deformation potentials), not just 3–4 local numbers.

## 4. Worked example — Silicon (conventional cubic cell)

Reuses the existing low-precision DFT sample (`samples/exercise_04_bulkSi_gs`,
LDA, `num_rgrid=12³`, `4×4×4` k). A runnable script is in
`samples/exercise_x4_Si_dft_to_epm/run_dft_to_epm.sh`.

```bash
SALMON=/path/to/build/salmon

# (1) DFT ground state + band structure -> band.dat (energies are Hartree)
cd samples/exercise_04_bulkSi_gs
$SALMON < Si_gs.inp ;  ln -sfn data_for_restart restart ;  $SALMON < Si_band.inp

# (2a) least-squares extraction
python3 ../../tools/dft_to_epm/dft_to_epm.py \
    --dft band.dat --format band_dat --cell cubic \
    --a-lattice-au 10.2626 --cutoff-ry 11.1 --material-name Si_fromDFT \
    --shells-s 3,8,11 --nval 16 --nbands-fit 18 --weight-valence 3.0 \
    --method lsq --out-prefix Si_lsq

# (2b) Zunger extraction via the vendored DeePseudopot module
python3 ../../tools/dft_to_epm/dft_to_epm.py \
    --dft band.dat --format band_dat --cell cubic \
    --a-lattice-au 10.2626 --cutoff-ry 11.1 --material-name Si_fromDFT \
    --shells-s 3,8,11 --nval 16 --nbands-fit 18 --weight-valence 3.0 \
    --method zunger --out-prefix Si_zunger
```

Both print a fit report and write `*_epm_formfactors.data`. Expected (this coarse
DFT): `V^S(3,8,11) ≈ −0.140, 0.042, 0.150 Ry`, band RMS ≈ 0.21 eV; the `zunger`
run reports `zunger backend: deepseudopot(vendored)` and the fitted `a0..a3`.

### Consume the result

* **Python EPM (primary):**
  ```python
  import epm_gaas_reference as epm
  epm.load_form_factor_file('Si_zunger_epm_formfactors.data')
  # then run the script with MATERIAL = 'file'
  ```
* **Fortran `theory='epm'` (deprecated):** in `&epm`, set
  `epm_material='file'` and `epm_formfactor_file='Si_zunger_epm_formfactors.data'`
  (see `samples/exercise_x4_Si_dft_to_epm/Si_epm_fromdft.inp`).

From there the standard `SYSNAME_{k,eigen,tm}.data` flow feeds `theory='sbe'`
unchanged.

## 5. Adapting to a new crystal

* **Diamond / monoatomic basis** (Si, Ge, C): leave `--shells-a` empty (`V^A≡0`),
  one species.
* **Zincblende / two species** (GaAs, InSb, …): pass `--shells-a 3,4,11`; `zunger`
  then fits two species (`a0..a3` each) and forms `V^S=(V_cat+V_an)/2`,
  `V^A=(V_cat−V_an)/2`.
* Match `--cell` to the DFT cell (cubic 8-atom vs 2-atom primitive), set
  `--a-lattice-au` to the same `a` you will give the EPM, and `--nval` to
  (valence electrons)/2 for that cell.
* Prefer `--method zunger` when you fit **many** shells or the DFT bands are
  coarse/noisy; `lsq` is fine (and fastest) for the minimal classic shell set.

## 6. Tests

```bash
python3 tools/dft_to_epm/tests/test_recovery.py
```
Recovers known Si factors to ~1e-14 (`lsq`, both cells) and ~1e-15 (`zunger`,
vendored backend) — no SALMON build required.

## References

* M. L. Cohen & T. K. Bergstresser, *Phys. Rev.* **141**, 789 (1966).
* L.-W. Wang & A. Zunger, *Phys. Rev. B* **51**, 17398 (1995) (analytic local form).
* K. Lin, M. J. Coley-O'Rourke & E. Rabani, *npj Comput. Mater.* **11**, 381
  (2025) (DeePseudopot) — vendored in `external/DeePseudopot/`.
