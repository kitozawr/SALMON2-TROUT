# Exercise x11 — the full-dissipation SHOWCASE: all four materials, every valid channel at once

**The demonstration run of this fork.** One EPM→SBE pipeline per material — **GaAs,
Si, CdS, graphene** — each with *every* dissipation/renormalization channel that is
**valid and cited for that material** switched on simultaneously, driven by an
analytic `Acos2` pulse (no field files). The deliverable is not a number but a
*movie*: watch the carriers get born by the vertical (k-conserving) transition and
then **roll into the correct valleys** frame by frame, while impact ionization
multiplies them, Auger recombines them, carrier–carrier thermalizes them and Σ^HF
renormalizes the gap under them.

Where `exercise_x8` (HHG) teaches you to *read* the band-resolved 3D BZ maps, x11
is what the maps are *for*: the same `--bz3d` frames now show population physics —
injection → intervalley transfer → cooling — with every channel contributing.

---

## 1. The channel matrix (what is ON, per material, and WHY the rest is off)

| channel (`&sbe`) | GaAs | Si | CdS | graphene |
|---|---|---|---|---|
| ring backbone `yn_sbe_superres` | ✅ | ✅ | ✅ | ✅ (mandatory) |
| e-ph `yn_sbe_eph` | ✅ polar LO + 5 intervalley | ✅ 6 intervalley g/f | ✅ Fröhlich LO (intra-valley) | ✅ E2g + A1′ (inter-k on the cone) |
| nonlocal II `yn_sbe_impact_ionization` | ✅ Stobbe quartic | ✅ Keldysh quadratic | ✅ E_th cited, **prefactor = FIT** (explicit `sbe_ii_prefactor`) | 🚫 gapless: no threshold law — CM is the Rana channel |
| Auger `yn_sbe_auger` | ✅ ring (II time-reverse, no extra C) | ✅ ring | ✅ ring (inherits the II fit scale) | ✅ **2D Rana [R07]** (quasi-Fermi R−G) |
| carrier-carrier `yn_sbe_eeh` | ✅ | ✅ | 🚫 no cited CdS rate (provenance gate) | 🚫 no cited rate |
| Coulomb Σ^HF `yn_sbe_coulomb` | ✅ | ✅ | ✅ (ε₀ = 8.9) | ✅ **2D sheet kernel** 2π/(ε_r·A·(q+κ)), substrate ε_r |
| dynamic λ²(n(t)) in the ring II/Auger | ✅ auto (registry) | — (λ=0 **correct**, Burt [L90]) | — | — |

**The three deliberate absences (the subtlety of a maximally-loaded *valid* run):**

1. **`yn_sbe_bgr_threshold = 'n'` everywhere.** Σ^HF already renormalizes the
   Houston eigenvalues that the II threshold is measured against; the BGR n^(1/3)
   law is the cheap *stand-in* for exactly that shift when Coulomb is off. Both at
   once would count the same physics twice — the init guard aborts by design.
2. **No Kuhn-Zurek `sbe_decoh_*`.** The collision channels amplitude-damp the
   coherences themselves (that is what CPTP dissipators do); adding single-particle
   wave-packet dephasing on top would double-count decoherence. (For graphene it is
   forbidden outright — gapless Dirac coherence loss is many-body.)
3. **No k-local impact ionization.** With the ring on, the momentum-conserving
   nonlocal II replaces the BZ-averaged k-local imitation automatically (it is
   gated off) — the created e–h pair lands in the *correct* Δ valleys instead of a
   BZ-average.

## 2. Run it (per material; heavy mode)

```bash
cd samples/exercise_x11_full_dissipation_showcase

# --- Si: the flagship valley-rolling movie -----------------------------------
../../build/salmon < Si_epm_gs.inp  > Si_gs.log      # EPM GS (writes Si_prim_*.data)
../../build/salmon < Si_sbe_rt.inp  > Si_rt.log      # the loaded SBE run

# same pattern for the others:
../../build/salmon < GaAs_epm_gs.inp     > GaAs_gs.log
../../build/salmon < GaAs_sbe_rt.inp     > GaAs_rt.log
../../build/salmon < CdS_epm_gs.inp      > CdS_gs.log
../../build/salmon < CdS_sbe_rt.inp      > CdS_rt.log
../../build/salmon < graphene_epm_gs.inp > graphene_gs.log
../../build/salmon < graphene_sbe_rt.inp > graphene_rt.log
```

**Heavy compute.** The ring II/Auger is O(nk³) per step — that is the point (it is
the honest 2-particle momentum-conserving sum). The outer k1 loop is
MPI-distributed (O(nk³/P)), so run the 8³ materials with MPI:

```bash
mpiexec -n 16 ../../build/salmon < Si_sbe_rt.inp > Si_rt.log   # pick P | nk (512)
```

Sanity while it runs: `grep "basis edge" *_rt.log` must stay empty (VG basis
sufficiency, the x8 lesson), and the `electrons` column must stay at
8.000/8.000/16.000/2.000 (CPTP — every channel conserves the trace by
construction).

## 3. Watch the physics (the pictures this exercise exists for)

```bash
python3 ../../plot_sbe_results.py -i . -o plots --snapshots --valleys \
        --bz3d --bz3d-voxel --spectral
```

| material | what the frames show |
|---|---|
| **Si** | cloud born at Γ (+ Γ–L diagonals: the vertical multiphoton channel, cf. `si_three_photon_isosurfaces.py`) → intervalley e-ph **drains it into the six Δ valleys at 0.85·X** — the indirect-gap signature picture |
| **GaAs** | needle-sharp Γ spot (Kane vertical tunnelling) → hot spill into the **eight L half-spots** (Gunn-physics Γ→L transfer); log prints the one-time `dynamic free-carrier screen active, lambda^2 = ...` line as the plasma builds |
| **CdS** | Γ cloud that **stays at Γ** while the Fröhlich ladder cools its energy (polar intra-valley cooling — no intervalley rolling; note the run drives E⊥c because E∥c is Γ9-dipole-forbidden) |
| **graphene** | carriers relax onto the **six K/K′ Dirac corners**; `*_sbe_nex.data` shows the 2D Rana Auger tail — monotone nex decay after the pulse (and *gentle generation* toward the 300 K density before it: detailed balance, not a bug) |

The frame set per run: `plots/*_bz3d_t*.png` (scatter, variant a),
`plots/*_bz3dvox_t*.png` (voxel cloud, variant b — the *shape* of the valleys;
unsmoothed by default, `--voxel-smooth 1.0` for a soft cloud on dense grids),
`plots/*_cart_snap_*.png` (Cartesian slices with valley markers),
`plots/spectral_frames/*.png` (A(k,E) band map: valence depleting, valleys
filling). `out_projection_k_step = 500` gives 12 frames per run — and each frame
series is **auto-assembled into an animation** (`plots/*_bz3d_anim.mp4` /
`_bz3dvox_anim.mp4` / … via ffmpeg, or `.gif` if ffmpeg is absent). No manual
step needed; tune with `--fps` / `--anim-format`, disable with `--no-animate`.
(Manual assembly still works if you want a custom rate:
`ffmpeg -framerate 4 -pattern_type glob -i 'plots/*_bz3d_t*.png' si_rolling.mp4`.)

## 4. Fine print

- **Grids must match** between the GS and SBE steps (`nelec/nstate/num_kgrid`).
  Odd grids (9³, and 7×7×5 for CdS) put Γ on the mesh — sharper injection spots at
  ~1.4× the k-count of 8³.
- **CdS `sbe_ii_prefactor = 1.03e12` is CALIBRATED** on our own BC1967 EPM bands
  by the Stobbe golden-rule procedure (`tools/cds_ii_calibrate.py`, fit RMS 3%) —
  the same standing as Stobbe's GaAs 2e12 (the registry sentinel still aborts
  without an explicit value, by design).
- **CdS acoustic (Rode E₁ = 14.5 eV) is ALWAYS TF-screened** `[q/(q+q_TF)]²` from
  the instantaneous carrier density — the bare DP channel would be unphysical at
  n ≥ 1e18 cm⁻³; screened, it is the fallback cooling below ħω_LO.
- **Graphene ε_r**: the Rana rates default to the R07 benchmark substrate
  (ε_r = 10); on SiO₂ set `sbe_coulomb_epsilon = 4.0d0` (≈2× the Auger rate).
- **Inputs are in `A_eV_fs` units** (lengths Å, energies eV, times fs; `&epm`
  constants stay in a.u. by definition). dt ≈ 0.0097 fs (= 0.4 a.u.) is fine
  for the population movie; for converged currents/HHG on
  top of full dissipation drop to 0.05–0.1 (x8 §2 note).
- Physics spec of every channel, with the full master equation: `wiki/08_master_equation.md`
  (and per-channel details in `wiki/01`, `wiki/07`).
