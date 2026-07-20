# Plotting & analysis tools

Everything that turns the raw SBE/EPM output into pictures and physical
diagnostics. Moved here from the README (2026-07-04 restructure); the README
keeps only the one-line pointers. The real-time output files themselves and
their cadence knobs are described in
[`04_configuration_examples.md`](04_configuration_examples.md) (`&analysis`).

---

## 1. `plot_sbe_results.py` — the one plotter

Self-contained (matplotlib + numpy, `Agg` backend, no interactive windows). Copy
into the calculation directory (or use `-i DIR -o OUTDIR`) and run. It scans for
`SYSNAME_sbe_rt_energy.data`, `SYSNAME_sbe_nex.data`, `SYSNAME_sbe_nex_k.data`
(+ `_lev_real`, `_unfold` variants, `_bandpath.data`, `band.dat`) and produces:

* line plots of total energy and excited-electron/hole counts vs time;
* from `SYSNAME_sbe_rt.data`, the **optical conductivity**
  $\sigma(\omega)=J(\omega)/E(\omega)$ along the driven axis: a global Re/Im
  spectrum over 0–4 THz and a strongly-overlapped **short-time Fourier**
  $\mathrm{Re}\,\sigma(\omega,t)$ map (Hann-windowed, hop = 1 sample by default).
  Tune `--fmax-thz`, `--stft-window-fs`, `--stft-hop`; disable
  `--no-conductivity`. Resolving 0–4 THz needs a ps-scale run (it warns
  otherwise);
* per saved time step, the Houston-basis lowest-conduction-band population as
  three 2D k-grid heatmap slices ($k_x$-$k_y$, $k_x$-$k_z$, $k_y$-$k_z$) — in
  reduced coordinates AND (when the `# b1/b2/b3` header is present, i.e. every
  primitive dataset) un-sheared into the **true Cartesian BZ**
  (`*_cart_snap_*`), plus reduced+Cartesian **k–t maps**;
* `--valleys`: overlay the material's high-symmetry/valley markers (Γ, the six X
  faces, the Δ minima at 0.85·X, L, K/K′...) on the Cartesian maps;
* `--bz3d`: the paper-style **3D Brillouin-zone population plot** per time step —
  the Wigner–Seitz BZ wireframe (Voronoi cell of the reciprocal lattice: FCC →
  truncated octahedron, hexagonal → prism) + the MP k-points coloured & sized by
  population (weak points fade/hidden);
* `--bz3d-voxel`: variant (b) — a **semi-transparent voxel cloud** rendered
  from the un-sheared Cartesian grid (opacity ∝ population): the volumetric
  *shape* of the populated valleys rather than the sampled points. By default
  **unsmoothed** — one cube per populated k-bin, faithful to `--bz3d`, with the
  colorbar on the true peak population. `--voxel-smooth SIGMA` (e.g. `1.0`) adds
  a Gaussian blur for a softer cloud on dense grids;
* `--bz3d-cb-sum`: by default the 3D/BZ views (and the k-maps) show the **lowest**
  conduction band only. This flag instead uses the **sum of BOTH recorded
  conduction bands, CB1 + CB2**, read from the four-level real file
  `SYSNAME_sbe_nex_k_lev_real.data` (the same file that colours the `--spectral`
  band movie — cols `pop_cb1`, `pop_cb2`). Use it when the drive fills the second
  conduction band (strong or high-photon-energy pulses): e.g. on a CdS resonant-Γ
  run at 2×10¹³ W/cm² CB2 already holds ~⅓ of the conduction carriers, invisible
  in the LCB-only view. Outputs get a `_cbsum` tag; falls back to lowest-CB if the
  four-level file is absent;
* `--spectral`: per-frame **A(k,E) band-structure movies** (`spectral_frames/`)
  coloured by occupation (valence full=1 at t=0, depleting; conduction
  filling); `--spectral-excitation` colours by excitation instead. Needs
  `SYSNAME_bandpath.data` (the EPM/`dft_band` emit it);
* band-structure plots from `SYSNAME_k.data`+`_eigen.data` and from `band.dat`
  (`theory='dft_band'`), energies shifted to VBM = 0 (`--band-vbm IDX` to
  override the half-filling default). **The high-symmetry path is auto-selected
  per material** (`--lattice auto`, default): the plotter reads `# material =`
  from the `_k.data` header and picks **fcc** for cubic GaAs/Si (`L-Γ-X-W-K`) or
  **wurtzite** for CdS (`A-Γ-M-K-Γ` — the hexagonal path, leading with the
  **Γ-A** c-axis segment, matching the `_bandpath.data` populations plot). Force
  it with `--lattice fcc|wurtzite`, or set `--band-path` explicitly. (The
  with-population `--spectral` plot always follows the exact path nodes the EPM
  wrote into `_bandpath.data`, so it needs no material knowledge.);
* **Auto-animation (default on).** Every per-frame series — band maps, k-maps,
  BZ snapshots (reduced + Cartesian), `--bz3d`/`--bz3d-voxel`, and the
  `spectral_frames/` — is stitched into `<stem>_anim.gif` (Pillow, the default,
  no external tool) as soon as ≥2 frames are written, so the evolution is
  viewable without a manual `ffmpeg` step. Single images (k–t maps, band
  structures, `*_vs_Time`) are skipped. Knobs: `--anim-format {gif,mp4,both,auto}`
  (default `gif`; `mp4` needs `ffmpeg`, `auto` = mp4-if-ffmpeg-else-gif), `--fps`
  (default 6), `--no-animate` to turn it off;
* `--only-bands`, `--log-cmap`, `--subtract-baseline`, `--snapshots`,
  `--instantaneous` — see `--help`.

```sh
python3 plot_sbe_results.py -i . -o plots --snapshots --valleys --bz3d --bz3d-voxel --spectral
# -> plots/*_anim.mp4 (or .gif) written automatically, one per frame series.
# Manual assembly is still possible if you want a custom frame rate/pattern:
#   ffmpeg -framerate 4 -pattern_type glob -i 'plots/*_bz3d_t*.png' movie.mp4
```

### How the plotter knows which Brillouin-zone cell to draw

The drawn cell is **not chosen by material name** — it is the **Wigner–Seitz
cell of the reciprocal lattice**, reconstructed from three numbers written into
the data itself:

1. **Source of truth — the `_k.data` header.** The non-orthogonal EPM (and the
   DFT export) write the primitive reciprocal vectors into the ground-state
   `SYSNAME_k.data` header as three comment lines

   ```
   # b1 =   <b1x>  <b1y>  <b1z>   [a.u.]
   # b2 =   ...
   # b3 =   ...
   ```

   The plotter reads them with `_read_bmatrix()` (regex on `# b{1,2,3} =`) into a
   3×3 matrix `b_matrix` (rows `b1,b2,b3`). Each `*_sbe_nex_k*.data` frame is
   paired to its sibling `SYSNAME_k.data` by `_bmatrix_for()`, so the cell always
   matches the run that produced the populations.

2. **The cell is the Voronoi cell of the origin.** From `b_matrix` alone:
   * **3D wireframe** (`_bz_wireframe_3d`): take the 27 nearest reciprocal-lattice
     points $\{i\,\mathbf b_1 + j\,\mathbf b_2 + k\,\mathbf b_3 : i,j,k\in\{-1,0,1\}\}$,
     build their Voronoi diagram, and keep the ridges that bound the origin's
     cell — those polygons are the BZ faces, their edges the wireframe.
   * **2D silhouette** (`_bz_outline_2d`): the Bragg planes
     $\mathbf k\cdot\mathbf G = |\mathbf G|^2/2$ of the same lattice, projected
     onto each Cartesian pair.
   * **k-point placement** (`_cartesian_bz_grid`): every k is mapped to its
     nearest reciprocal-lattice image (the same WS wrap) before being binned
     onto the Cartesian voxel/heatmap grid.

   Because it is a pure Voronoi construction, the **shape follows the vectors
   automatically**: an FCC dataset (GaAs, Si) draws a **truncated octahedron**, a
   hexagonal dataset (CdS, graphene) draws a **hexagonal prism** — with nothing
   hardcoded per material. This is the same cell for `--bz3d`, `--bz3d-voxel`,
   the `*_cart_snap_*` slices, and the Cartesian k–t maps.

3. **Fallback when the header is absent.** Legacy cubic/orthogonal datasets do
   not carry `# b1/b2/b3`, so `b_matrix` is `None`: the true-BZ views are skipped
   and the plotter shows the **reduced-coordinate** (fractional-axis) heatmaps
   instead. `--bz3d`/`--bz3d-voxel` therefore require a primitive dataset with
   the header.

4. **Caveat — `--valleys` markers are cubic.** The Γ / X / Δ (0.85·X) / L overlay
   points come from `_cubic_valleys()`, which assumes an **FCC/diamond** cell.
   They are meaningful on GaAs/Si; do **not** trust them on the hexagonal
   CdS/graphene cells (the BZ wireframe there is still correct — only the valley
   markers assume cubic symmetry).

### Spinor (spin-orbit split) datasets

Detected automatically from the occupation column (1 per band). The plotter
**sums the spins of each level**: adjacent Kramers-partner sub-bands merge into
one level (occupations summed 1+1=2, energy = pair mean), so tiny Dresselhaus
splittings don't render as doubled lines while the real spin-orbit splittings
(Γ₈/Γ₇, Δ₀ = 0.341 eV for GaAs) stay visible. `--spin-sum {auto,on,off}`.
For spinor `bandpath` data it additionally renders the **Dresselhaus spin
splitting** Δⱼ(k) panel in meV (zero along [100]/[111] by symmetry, ~10–140 meV
peaks near W/K for GaAs).

---

## 2. Folded vs unfolded band pictures

The cubic 8-atom cell is a supercell of 4 primitive FCC cells: its MP-grid band
plot shows the primitive bands **folded 4-fold** (dense conduction crossings —
a representation artifact, not physics; cf. Popescu–Zunger;
Quan–Rybin–Scheffler–Carbogno PRB 113, 085112 (2026)). The folding is **exact**
here (parity selection rule ⇒ block-diagonal Hamiltonian, asserted at runtime),
so the clean primitive picture is recovered exactly:

```sh
python3 epm_gaas_reference.py bandpath   # SYSNAME_bandpath.data along L-Γ-X-W-K-Γ
python3 plot_sbe_results.py --only-bands # bandpath_*.png (+ spin-splitting panel)
```

*(On the primitive-cell pipeline — `epm_cell='primitive'`, the default — there is
no folding in the first place; `theory='epm'` emits `_bandpath.data` directly.)*

### The unfolded k-resolved population pipeline (folded cells)

The supercell branch index `nelec+1` mixes *different physical primitive bands*
from k to k. The three-stage pipeline — EPM → SALMON → plotter, **no GS
regeneration** — reports instead the population of the physical gap-edge bands
(VB-1, VB, CB1, CB2; spins summed) at every unfolded primitive BZ point:

1. **EPM, once:** `python3 epm_gaas_reference.py unfoldmap` →
   `SYSNAME_unfold.data` — the spectral weights
   $w_s = |\langle\psi|P_s|\psi\rangle|^2$ ($\sum_s w_s = 1$) of every cubic band
   on the 4 FCC sublattices + the energy-ranked primitive band index. Must match
   the GS k-grid (checked). Also `... bandpath` for the A(k,E) skeleton.
2. **SALMON:** with the map present, the SBE writes
   `SYSNAME_sbe_nex_k_unfold.data` — per saved time, the population of the four
   physical levels at each primitive point $k_{\rm prim} = k_{\rm sc}+G_0(s)$,
   distributed over sublattices **by the spectral weights** (not argmax — exact
   at degeneracies). The impact-ionization channel is likewise
   sublattice-resolved when the map is present.
3. **Plotter:** picks `_unfold` up automatically → `nex_k_unfold_*` (primitive-BZ
   maps; X-valley satellites at the zone boundary are *physics*),
   `nex_k_fold_*` (cosets summed back onto the cubic zone), and with
   `--spectral` the per-frame A(k,E) movies.

Theory of the folding itself (parity rule, cosets, N-coset generalization):
[`05_folding_unfolding.md`](05_folding_unfolding.md).

---

## 3. Standalone band-structure / injection-physics probes

All reuse the repo's own EPM machinery (same folding, spin-orbit and momentum
routines) — no external band data.

### `band_field_coupling.py` (GaAs)
Unfolded bands, group velocity, **field-projected interband coupling**
$|\langle cb|p_{\hat e}|v\rangle|^2$ over the Γ₈ manifold, and the direct gap
along [100]/[111]/[110]. Quantifies that the diagonal occupation pattern in
folded maps is band-**folding**, not a matrix-element effect.

### `zener_tunneling_estimate.py` (GaAs — the TUNNELLING probe)
Transverse-k-resolved vertical **Zener/Kane** probability
$P_{\rm Kane}(k_\perp) = \exp(-C\,m_r^{1/2}E_g^{3/2}/F)$ + linear **Landau–Zener**
at genuine off-Γ avoided crossings, vs field; `--map2d` adds the 2-D transverse
birth map $W(k_y,k_z)$ at $k_\parallel = 0$ (the injection "needle" FWHM at Γ,
with the folded conduction geometry overlaid — the proof that diagonal SBE
weight sits at fold positions, not the injection blob).

### `si_three_photon_isosurfaces.py` (Si — the MULTIPHOTON probe)
The power-law $I^N$ counterpart (different prefactor — *not* a tunnelling
exponential): the N-photon direct rate
$W_N(k) \sim I^N |M_N(k)|^2\, g_{\rm FK}(E_{\rm dir}(k)-N\hbar\omega)$ on a 3-D BZ
grid — full LOPT ladder over all intermediate bands, exact polarization
orientation-average (keeps O_h ⇒ IBZ-evaluated), **Franz–Keldysh broadening**
$\hbar\theta = (F^2/2\mu)^{1/3}$ whose sub-edge Airy tail lets the below-edge
3-photon drive (3.0 eV < 3.34 eV) reach the gap along the low-gap **Γ–L ⟨111⟩
valleys**, while 4-photon (4.0 eV) is a genuine shell reaching toward (never at)
X. The $(E_0/2\omega)^{2N}$ prefactor puts orders on one absolute scale
($W_3+W_4$ summable; 4γ share ≈4%→34% from 3→10 MV/cm). Outputs: plotly 3-D
isosurfaces / resonance shells / rate cloud, the summed rate + axis-averaged
projections (PNG+HTML), an MRI-style ortho-slicer (`--slice-quantity`), Γ-L/Γ-X/
Γ-K line scans vs field. The SBE br3d maps show the same Γ–L diagonal pattern
(x8 README note): HHG carriers are seeded along Γ–L by this vertical transition,
then drift to the Δ valleys. **Per-k FK gotcha** (wiki/00 decisions log): do NOT
use the Airy JDOS lineshape per-k — it double-counts the explicit k-sum.

> **Spinor velocity in the probes.** For spinor datasets the physical interband
> velocity is $\langle m|(k+G) + \nabla_k H_{SO}|n\rangle$;
> `band_field_coupling.primitive_bands_momentum` adds the $\mu\nabla_k H_{SO}$
> correction (the Zener script inherits it) — small near Γ, O(1) relative effect
> for the nearly-p-dark heavy-hole transition. **Caveat for all probes:** the
> local EPM undershoots gaps (GaAs Γ: 1.27 vs 1.42/1.52 eV exp.), so quote
> Γ-normalized probabilities / relative geometry unless you scissor-correct.

### `hhg_spectrum.py` (in `samples/exercise_x08...`, `x10`, `x11`)
HHG intensity $|\omega J(\omega)|^2$ vs harmonic order (log scale); only odd
orders for centrosymmetric Si — the clean check against TDDFT harmonics.
