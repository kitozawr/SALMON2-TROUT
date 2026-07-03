# Exercise x7 — PRIMITIVE-CELL EPM → SBE (non-orthogonal, NO folding)

This is the **decisive-test** workflow: run the SBE on the **primitive cell
directly** — non-orthogonal lattice vectors, no supercell, no band folding, no
cosets, no unfold map, no sublattice projection. The primitive cell **is** the
irreducible problem, so the spurious folded-manifold physics (dense avoided
crossings between cosets) simply does not exist.

**Why this matters.** On the *folded* cubic supercell, GaAs shows an unphysical
"L over-populates Γ" pattern (L/Γ ≈ 760, anti-Zener). On the *primitive* FCC cell
the ordering is the correct Kane one — **Γ ≫ L ≫ X** (L/Γ ≈ 0.01, matching the
Zener tunnelling estimate). Removing folding swings L/Γ by ~10⁵. See
`wiki/05_folding_unfolding.md` §8.

The ground state can come from **either** the in-SALMON **Fortran EPM**
(`theory='epm'`, `epm_cell='primitive'` — now the DEFAULT) **or** the **Python EPM
primitive references**. For **all four** materials the two are **verified
interchangeable** (identical k-points; band energies to 5e-11 Ha for GaAs/Si/CdS,
~6e-8 eV for graphene), so you no longer need Python to generate the GS.
The clean SBE is **k-grid-agnostic** (it propagates band energies + momentum
matrix elements, and the Coulomb kernel is metric-aware), so a non-orthogonal /
triclinic k-grid needs no special handling — the reduced k it stores are labels.

| Material | Cell | Fortran EPM GS step | Python generator | Gap (validated) |
|---|---|---|---|---|
| GaAs (scalar) | FCC 2-atom rhombohedral | `GaAs_prim_epm_gs.inp` ✅ | `epm_gaas_primitive.py` | Γ 1.39 / L 2.68 / X 3.94 eV |
| GaAs (spin-orbit) | same + SO | — (spinor Fortran pending) | `epm_gaas_primitive.py` (`INCLUDE_SPIN_ORBIT`) | Δ₀=0.341 eV, gap 1.27 eV |
| Si | FCC 2-atom (diamond, V^A=0) | `Si_prim_epm_gs.inp` ✅ | `epm_si_primitive.py` | indirect 1.06 eV @ 0.85·X |
| CdS | wurtzite 4-atom **hexagonal** | `CdS_prim_epm_gs.inp` ✅ | `epm_cds_primitive.py` | direct 2.55 eV (BC1967 2.58) |
| graphene | 2-atom **hexagonal** (2D-in-vacuum) | `graphene_prim_epm_gs.inp` ✅ | `epm_graphene_primitive.py` | gapless Dirac at K |

## Build

```sh
cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"
```

## Run (pick a material)

Each EPM generator writes `SYSNAME_k/_eigen/_tm.data` (+ `_bandpath.data`) **into
the current directory** — run it where you will run `salmon`. The `&system`
`al_vec1/2/3`, `nelec`, `nstate`, and `&kgrid num_kgrid` in the `.inp` MUST match
what the generator used.

```sh
# --- GaAs scalar (8x8x8) ---------------------------------------------------
# in-SALMON Fortran EPM (default primitive) -- no Python needed: writes
# GaAs_prim_k/_eigen/_tm.data AND GaAs_prim_bandpath.data (the clean primitive
# band path along L-Gamma-X-W-K-Gamma, verified interchangeable with the
# Python p.main_bandpath() to 5e-11 Ha at every node).
./build/salmon < GaAs_prim_epm_gs.inp
#   (equivalently: python3 -c "import epm_gaas_primitive as p; p.main_gs(); p.main_bandpath()")
./build/salmon < GaAs_prim_sbe_rt.inp
python3 plot_sbe_results.py -i . -o plots --snapshots --spectral

# --- GaAs spin-orbit (4x4x4; 16 spinor bands) ------------------------------
# in-SALMON Fortran spinor EPM (yn_spinorbit='y'; mu calibrated at Gamma to the
# cited Delta0=0.341 eV; verified interchangeable with the Python ref to 5e-11 Ha):
./build/salmon < GaAs_prim_so_epm_gs.inp
#   (equivalently: python3 -c "import epm_gaas_primitive as p; p.INCLUDE_SPIN_ORBIT=True; \
#               p.NUM_KGRID=(4,4,4); p.SYSNAME='GaAs_prim_so'; p.main_gs(); p.main_bandpath()")
./build/salmon < GaAs_prim_so_sbe_rt.inp     # NOTE: needs yn_sbe_spinor='y'
python3 plot_sbe_results.py -i . -o plots --snapshots

# --- Si (8x8x8; cutoff 27 Ry to resolve the Delta valley) ------------------
./build/salmon < Si_prim_epm_gs.inp          # in-SALMON Fortran EPM (cutoff 27 Ry)
python3 -c "import epm_si_primitive as s; s.configure_for_si('Si'); s.prim.main_bandpath()"
./build/salmon < Si_prim_sbe_rt.inp
python3 plot_sbe_results.py -i . -o plots --snapshots --spectral

# --- CdS (7x7x5 hexagonal wurtzite) ----------------------------------------
./build/salmon < CdS_prim_epm_gs.inp         # in-SALMON Fortran EPM (4-atom polar wurtzite)
python3 -c "import epm_cds_primitive as c; c.main_bandpath()"
./build/salmon < CdS_prim_sbe_rt.inp
python3 plot_sbe_results.py -i . -o plots --snapshots

# --- graphene (12x12x1 hexagonal, 2D + vacuum; gapless Dirac) ---------------
./build/salmon < graphene_prim_epm_gs.inp    # in-SALMON Fortran EPM (cutoff 2.94 Ry = 40 eV)
python3 -c "import epm_graphene_primitive as g; g.main_bandpath()"
./build/salmon < graphene_prim_sbe_rt.inp    # in-plane field; trace=2 conserved
python3 plot_sbe_results.py -i . -o plots --snapshots
# -> the kx-ky Cartesian heatmap shows carriers at the six K/K' Dirac corners.
```

## What to look at

- `plots/<sys>_sbe_nex_k_real_cart_snap_*.png` — the **Cartesian-BZ heatmap**
  (real carriers). The plotter reads the reciprocal vectors `# b1/# b2/# b3` from
  the GS `k.data` header and **un-shears the triclinic k-grid into a regular
  Wigner–Seitz Cartesian volume** — a true picture of the Brillouin zone that
  gets smoother the denser the grid. For GaAs/Si carriers sit at Γ (small-gap
  pocket); CdS under a c-axis field is suppressed at Γ (see "selection rule").
- `plots/spectral_frames/<sys>_spectral_*` — **one frame per time step**: the
  conduction band coloured by population (broadened by carrier kinetic energy),
  fixed colour scale. Assemble into a movie:
  `ffmpeg -framerate 8 -pattern_type glob -i 'plots/spectral_frames/*_path_*.png' movie.mp4`.
- `plots/<sys>_sbe_nex_k_real_ktmap_k{x,y,z}.png` — population vs (time, k).
- `plots/band_structure_*` and `plots/bandpath_*` — the clean primitive bands
  (spinor runs also get `bandpath_spin_splitting_*` — zero along ⟨100⟩/⟨111⟩).
- `<sys>_sbe_rt_energy.data` — the total trace is conserved (CPTP).

## KEY GOTCHAS (read before changing the inputs)

1. **Spinor GS needs `yn_sbe_spinor='y'` in `&sbe`.** Without it the solver reads
   the 16-band spin-orbit GS as scalar (`nb_vb = nelec/2`), mis-identifies the
   lowest conduction band, and caps occupation at 2 instead of 1 → garbage
   populations. Scalar runs leave it off (default `'n'`).
2. **Use an ODD k-grid (7³, 9³, 7×7×5…) to sample Γ explicitly.** Even grids
   straddle Γ (nearest point at 0.1875), so the Γ population is never sampled.
   Odd grids put Γ on the mesh and resolve the sharp near-Γ excitation pocket.
3. **Si needs `PW_CUTOFF_RY = 27`** (set by `epm_si_primitive`) — its Δ-valley
   camel-back is unresolved at the GaAs cutoff (11.1) and the CBM mislands at X.
   The cutoff only sizes the one-off GS diagonalization, not the SBE.
4. **CdS with a c-axis (E∥c) field** suppresses excitation AT Γ: the top valence
   Γ9→conduction transition is dipole-forbidden for E∥c (only Γ7 couples), so
   carriers peak off-axis. This is real wurtzite selection-rule physics. Drive
   E⊥c (e.g. `epdir_re1 = 1,0,0`) to populate Γ.

## Dissipators + super-mode (CPTP)

The `&epm` block selects the per-material dissipation tables; `&sbe` enables the
channels (all default OFF). Verified on the primitive cells:
- **Si** + `yn_sbe_superres/eph/eeh/impact_ionization`: trace conserved, energy
  relaxes ~12% post-pulse.
- **CdS** + `yn_sbe_superres/eph/auger`: trace = 16 exact, carriers cooled by the
  Fröhlich polar-optical e-ph + Auger recombination (`eeh` is provenance-forbidden
  for CdS; impact ionization needs an explicit `sbe_ii_prefactor`).

**⚠️ e-ph on the primitive cell — intervalley caveat.** The k-local e-ph search
relaxes band-to-band at the SAME k. That is correct for **CdS** (its only mode is
Fröhlich polar-optical = intra-valley) and for **folded** cells (valleys fold onto
same-k bands), but **Si/GaAs intervalley** modes (Δ↔Δ, Γ–L–X) live at *different*
k in the primitive cell, so true intervalley transfer needs the **inter-k e-ph
that rides the super-mode ring** — gated on `yn_sbe_superres` (design locked in
`wiki/00`, implementation pending). Until then, primitive Si/GaAs e-ph is the
intra-k approximation.

See `wiki/00_implementation_status.md` (the live status + standing TODOs) and
`wiki/05_folding_unfolding.md` §8.
