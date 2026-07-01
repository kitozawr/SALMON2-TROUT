# Exercise x8 — Si primitive-cell SBE: velocity-gauge **basis sufficiency** + output analysis

A focused case study on **monolayer-clean Si** (FCC 2-atom primitive cell, no folding):
drive a strong pulse, learn to **check that the velocity-gauge band basis is
sufficient**, and learn to **read the output** (real-carrier band map, Cartesian
BZ maps, intra-band current, HHG spectrum). Everything runs from `build/salmon`
plus the repo plotters — no hand-editing of data.

All inputs are in **A/eV/fs** units. The field files' first column is time `[fs]`
(never scaled); columns 2–4 are the vector potential `Ax,Ay,Az` — the four files
below are the **same pulse shape (same frequency), only the amplitude scaled**.

| field file | `A_max` | role |
|---|---|---|
| `100.txt` | 2.76 | **over-drives** the basis — use it to *see* the VG overflow |
| `10.txt`  | 0.28 | clean here (the input default) |
| `1.txt`   | 0.028 | weak / perturbative |
| `0.1.txt` | 0.0028 | linear response |

Regenerate any scaled file from `100.txt` (divide only cols 2–4):
```bash
awk '{printf "%s %.15g %.15g %.15g\n", $1, $2/10, $3/10, $4/10}' 100.txt > 10.txt
```

---

## 1. Run the pipeline

```bash
cd samples/exercise_x8_Si_primitive_hhg_basis

# (1) ground state: in-SALMON Fortran EPM (Si FCC primitive, nstate=20)
../../build/salmon < Si_prim_epm_gs.inp > gs.log

# (2) clean primitive band path for the A(k,E) skeleton
#     (the Fortran EPM does not yet emit _bandpath.data -- use the Python ref)
python3 ../../epm_si_primitive.py bandpath          # writes Si_prim_bandpath.data

# (3) real-time SBE (edit file_input1 in the .inp to pick the field amplitude)
../../build/salmon < Si_prim_sbe_rt.inp > rt.log
```

---

## 2. Check the basis is sufficient  ← the point of this exercise

The velocity gauge propagates in a **fixed, finite band set** (`nstate`). A strong
field shifts the crystal momentum `k -> k + A(t)`; if it pushes carriers past the
**top band**, the basis is exhausted and the results are garbage. Two independent
checks (see `wiki/06_vg_basis_nb_convergence.md`):

**(a) `P_top` monitor — the top-band population.** SALMON warns on stderr/log:
```bash
grep "basis edge" rt.log        # 0 lines = basis OK ; any line = overflow (raise nstate or lower A)
```

**(b) the excited-carrier count `nex` must stay physical** (`0 <= nex <= nelec`):
```bash
awk '/^ +[0-9]/{print $2, $7}' rt.log | sort -k2 -g | tail -3   # last col = nex ; must not exceed 8
```
If `nex` blows past `nelec=8` (or goes negative) the basis overflowed.

**Try it:** run with `file_input1='100.txt'` — you will see `basis edge` warnings
around **t≈43 fs** and `nex` shoot up to ~11 (unphysical). Then switch to
`10.txt`: 0 warnings, `nex` stays `<0.3`, `electrons=8.000` conserved throughout.
That contrast **is** the lesson: strengthen the field only as far as the basis
(the `nstate` conduction headroom) can follow — otherwise weaken the field.

> **Two separate knobs.** `P_top`/`nex` overflow ⇒ **field too strong for `nstate`**
> (raise `nstate` or lower `A`). Per-*level* populations leaving `[0, occ]` at a
> *weak* field instead point to **`dt` too large** for the ~14 eV band spread
> (the CF4 needs `dt ≲ 0.05 fs` for the fastest inter-band coherence). `dt=0.1 fs`
> is fast for a first look; drop to `0.02–0.05 fs` to converge the current/HHG
> against your TDDFT reference.

---

## 3. Read the output

```bash
python3 ../../plot_sbe_results.py -i . -o plots --spectral --snapshots --valleys
python3 hhg_spectrum.py Si_prim_sbe_rt.data --out plots/hhg_spectrum.png
```

What each shows:

* **`plots/spectral_frames/nex_k_prim_spectral_path_*.png`** — the A(k,E) **band map
  coloured by OCCUPATION** `f = pop/occ` (valence **full = 1** at t=0, watch it
  deplete; conduction fills). This is the REAL (fixed-basis, diabatic) carrier
  population, not the virtual Houston one. Add `--spectral-excitation` to colour by
  excitation instead (holes in VB, electrons in CB, both 0 at t=0).
* **`plots/*_snap_*` (Cartesian BZ)** — the REAL-carrier LCB population as
  **kx-ky / kx-kz / ky-kz** heatmaps. With `--valleys` the Si high-symmetry points
  are overlaid (Γ, the six X faces, and the Δ-valley minima at 0.85·X along ⟨100⟩),
  so you can check the excitation **hot spots land in the correct Δ-valleys**.
* **`Si_prim_sbe_rt.data` (total current) / `Si_prim_sbe_intra_current.data`** —
  the gauge-invariant matter current `Jm` and its intra-band part (plotted by
  `plot_sbe_results.py` as `Jm_*` / `J_intra`).
* **`plots/hhg_spectrum.png`** — HHG intensity `|ω·J(ω)|²` on a **log scale** vs
  harmonic order; only **odd** orders survive (Si is centrosymmetric) — a clean
  check against your TDDFT harmonics.

---

## Notes
* `nstate=20` (16 conduction bands) is the VG headroom; raise it if a stronger
  field is needed. `nelec/nstate/num_kgrid` must match between the GS and SBE steps.
* The clean SBE reads only band energies + momentum from the GS, so the
  non-orthogonal (triclinic) primitive k-grid needs no special handling; the
  `al_vec` set the metric for the Cartesian BZ un-shearing and the cell volume.
* The total current `Jm` is gauge-invariant (the velocity-gauge A² term is internal);
  the intra/inter split is meaningful only in the Houston basis [Otobe, PRB 94, 235152].
