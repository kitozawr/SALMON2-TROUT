# Exercise x12 — Si frozen-core **phonon-assisted indirect** (Γ→X) generation

Si is an **indirect-gap** semiconductor: the valence-band maximum (Γ) and the
conduction-band minimum (the Δ/X valley) sit at different crystal momenta, so
band-edge absorption **cannot** be a vertical optical transition — it needs a
**phonon** to carry the momentum. This example drives Si with a near-tunneling
IR field in a **frozen-core 2V+2C Houston window** and shows the e-ph ring
opening exactly that indirect channel: turn the phonon off and the real-carrier
density collapses; turn it on and carriers are generated in the **lowest**
conduction band at the **X-valley** and cooled toward the band minimum.

It is the open-system counterpart to the frozen-core discussion in
[`wiki/03` §12](../../wiki/03_numerical_methods.md) and a direct validation that
the dissipators are **exactly electron-number conserving on a frozen window**.

## The active window (2 valence + 2 conduction)

`frozen_core_threshold_ev`/`frozen_free_threshold_ev` select the gap-edge
Houston window; the deep valence and high conduction bands are **frozen from
dissipation** but still evolve under the **full-basis VG unitary** (so a strong
field can push virtual population up through them and back — basis sufficiency).
For the 4³ Si primitive cell at Γ the window is

| active band | character | E−E_F |
|---|---|---|
| 3, 4 | VBM (Γ₂₅′ doublet) — **2 valence** | −1.70 eV |
| 5 | CBM (Γ₂′ singlet) | +1.70 eV |
| 6, 7 | next CB (Γ₁₅ doublet) | +3.51 eV |

so the conduction side is the CBM singlet **+** the degenerate Γ₁₅ pair (the two
lowest conduction *levels*). Bands 1–2 and 8–16 are frozen. `nstate=16` gives
the unitary its headroom; the ring/dissipators only ever see the 5 active bands.

## Run

```bash
cd samples/exercise_x12_Si_frozen_phonon_indirect

# (1) ground state (in-SALMON EPM)
../../build/salmon < Si_prim_epm_gs.inp > gs.log

# (2) phonon ON  (indirect channel open)
../../build/salmon < Si_frozen_phonon_rt.inp > on.log
cp Si_prim_sbe_nex.data on_nex.data ; cp Si_prim_sbe_nex_k_real.data on_nkr.data

# (3) phonon OFF (indirect channel closed) — flip yn_sbe_eph -> 'n'
sed "s/yn_sbe_eph               = 'y'/yn_sbe_eph               = 'n'/;\
     s/yn_sbe_eph_acoustic      = 'y'/yn_sbe_eph_acoustic      = 'n'/" \
    Si_frozen_phonon_rt.inp > off.inp
../../build/salmon < off.inp > off.log
cp Si_prim_sbe_nex.data off_nex.data ; cp Si_prim_sbe_nex_k_real.data off_nkr.data

# (4) plot: nex(t), energy distribution, cooling curve
python3 phonon_analysis.py         # -> phonon_assisted_demo.png
```

## Running the 9×9×9 production benchmark (SLURM / `sbatch`)

The 4³ default above resolves the physics but samples the indirect **X-valley**
shell only coarsely. The Keldysh-bracket verdict (`rate_benchmark.py`,
`calibration_scan.py`) wants the **9³ = 729 k-point** grid. A ready SLURM script
is provided: **[`run_9x9x9_lomonosov.sbatch`](run_9x9x9_lomonosov.sbatch)** —
edit the marked lines and submit with `sbatch run_9x9x9_lomonosov.sbatch`.

### What changes from 4³ → 9³ (only the grid)

Two dedicated input files ship with the grid already set to `9,9,9`:

| step | file | `theory` | `num_kgrid` |
|---|---|---|---|
| 1. ground state | `Si_prim_epm_gs_9x9x9.inp` | `epm` | `9, 9, 9` |
| 2. real time | `Si_frozen_phonon_rt_9x9x9.inp` | `sbe` | `9, 9, 9` |

**`num_kgrid` is the only line that differs from the 4³ pair.** `nelec` (8),
`nstate` (16) and `al` MUST be identical between the GS and RT step — the SBE
reads `Si_prim_{k,eigen,tm}.data` written by the GS, so a mismatch is rejected.
(To scale a *different* exercise, edit `num_kgrid(1:3)` in **both** its GS and
RT `.inp` — nothing else.)

> **Γ-anchored frozen window (9³).** The frozen-core window is built from the
> **Γ point** (the reduced-coordinate k closest to the origin — exactly (0,0,0)
> on the odd 9³ grid, `k-point 365`), identical on every MPI rank, not from
> whichever k happens to be index 1. At Γ, Si shows its clean gap-edge
> degeneracies (the VBM triplet Γ₂₅′ + the CB manifold), so the −3/+6 eV window
> selects **7 active bands** (VBM 2–4 + CB 5–8, Fermi ≈ 11.89 eV) — the correct
> gap-edge set. The run banner prints `Gamma reference: k-point … at (0,0,0)`
> and `n_active_bands = 7 / 16`; a window that selects `< 2` bands aborts with a
> diagnostic (Fermi + window bounds) instead of an opaque crash.

### How SALMON takes the input under `sbatch` — there is **no `-i` flag**

SALMON reads its namelist from **stdin**, so in the batch script you pipe the
`.inp` in:

```bash
srun ./salmon < Si_frozen_phonon_rt_9x9x9.inp > rt.log
```

Only **MPI rank 0** reads stdin — it copies the namelist to `.namelist.tmp`,
which every rank then reads (the parsed values are also broadcast). Practical
consequences:

- **There is no command-line input argument.** `srun ./salmon Si_..._rt.inp`
  does nothing useful; it must be the stdin redirect `< file`.
- **Run from a shared-filesystem directory** (your `scratch`): rank 0 writes
  `.namelist.tmp` and the GS `*.data`, and the other nodes read them back. Your
  crashing path `Z:\...\projects\Si_TROUT_100kV` is fine as long as every node
  mounts it.
- Two SBE runs in the **same directory** overwrite each other's
  `Si_prim_sbe_*.data`; the script `cp`s the outputs to `on_*`/`off_*` between
  steps (and if you launch several fields at once, give each its own subdir).

### How many MPI ranks

729 k-points are split across ranks **evenly (±1 each)** — `nproc` does **not**
have to divide 729. Rules of thumb:

- Use **`nproc ≤ 729`** (more ranks than k-points leaves ranks idle and is
  wasteful — with the fix below it is safe, just pointless).
- For perfect balance pick a **divisor of 729 = 3⁶**: **27** (27 k/rank),
  **81** (9 k/rank), or **243** (3 k/rank). The template uses **81**.
- The inter-k ring does `nproc` communication hops per step, so very large
  `nproc` trades compute for communication; 27–81 is the sweet spot here.
- Hybrid MPI×OpenMP works too (the ring gather is OpenMP-parallel): e.g. 27 MPI
  ranks × `OMP_NUM_THREADS` cores each.

### ⚠️ Build from the fixed code first

The 9³ multi-node run is exactly the one that segfaulted on Lomonosov-2 (rank 29
of 42, in `compute_coulomb_selfenergy_ring`, at the first step) while the same
input ran on the Intel-oneAPI mini-cluster. That was a **non-synchronized
distributed-start** race in the frozen-core broadcast and the field-file read,
now fixed (wiki/00 decisions-log, 2026-07-21). **Rebuild** on Lomonosov-2 from
`develop-2.0.0` (with this fix merged) before the production run — the fix is
compiler-agnostic (plain `MPI_Allreduce`/`MPI_Bcast`/`MPI_Barrier`), so it
behaves the same under OpenMPI+gfortran and Intel MPI. If a genuinely
inconsistent start still occurs you now get a clean, collective `error stop`
with a diagnostic instead of a lone-rank SIGSEGV.

## What you should see (`phonon_assisted_demo.png`)

| | eph OFF (field only) | eph ON (phonon) |
|---|---|---|
| final `nex` | 5.9×10¹⁸ cm⁻³ | 1.1×10²¹ cm⁻³ (**≈180×**) |
| ⟨E−E_CBM⟩ (lowest CB) | 1.0 eV (hot) | 0.42 eV (cooled toward the valley) |
| electrons (every step) | **8.000** | **8.000** |

The phonon opens the indirect Γ→X channel (≈180× more real carriers — most of
Si's band-edge absorption is phonon-assisted), lands them in the **lowest**
conduction band (band-5 / band-6 population ratio ≈ 33×), and cools them from
~1.2 → 0.42 eV above the CBM by phonon emission.

> Reference numbers regenerated **2026-07-19** after the **CP extension**
> (`wiki/00`): active↔frozen coherences now carry the same loss-Kraus factors
> as the active block, so a real scattering event decoheres the reversible
> frozen-band dressing instead of leaving ρ non-positive. The eph-OFF column
> is bit-identical to the pre-fix build.

### What the phonon actually scatters here (real vs virtual)

This pulse is **sub-cycle** (envelope 19 fs < period 2π/ω = 33 fs): the
coherent run shows the Houston conduction population *adiabatically follows
the field envelope* (corr(n_H, A²) = 0.99, peak 0.56 e⁻/supercell) and 99.9 %
of it returns to the valence band at switch-off. The e-ph conversion of this
long-lived **virtual** dressing is golden-rule phonon-assisted tunneling
(Hurkx/Kane BTBT physics — real for indirect Si), but its scale here is set
by ν_sat, the *real-carrier* intervalley rate: treat the absolute yield as an
**upper estimate** (BTBT parametrizations sit orders below it at this field).
A lifetime gate cannot separate this dressing (`sbe_ring_gate_fs`, default
**off**, is kept as an experimental knob: at τ = 2π/E_gap it still passes
~60 % of the envelope-following dressing and only trims sub-cycle fringes,
−41 % here).

## Calibration scan (`calibration_scan.py`) — 4³ worked example, redo on your grid

Field scan at 4³ (this pulse, eph-ON minus eph-OFF vs the Keldysh
1.7-eV-as-direct upper bound):

| I [W/cm²] | G_SBE [cm⁻³] | n_Keldysh(1.7) | ratio |
|---|---|---|---|
| 3×10¹⁰ | 2.9×10²⁰ | 1.3×10¹³ | 2.2×10⁷ |
| 1×10¹¹ | 1.1×10²¹ | 2.7×10¹⁶ | 4.0×10⁴ |
| 3×10¹¹ | 3.1×10²¹ | 2.9×10¹⁸ | 1.1×10³ |
| 1×10¹² | 7.9×10²¹ | 1.2×10²⁰ | 68 |

**The ratio is NOT flat**: G_SBE ∝ I^0.94 (conversion of the *first-order*
virtual polarization, ∝A², through Houston anticrossings within σ_E), while
the bracket is strongly nonlinear (tunneling exponent). So
`sbe_eph_interband_scale` is a **per-working-point** calibration, not a
universal constant: pick it from the ratio at YOUR field on YOUR grid
(9³ production runs), or leave 1.0 and read the yield as an upper estimate.
The σ_E matching width (`sbe_search_sigma_e_ev`) is a *moderate* lever
(4³, I=10¹¹: σ = 0.05/0.10/0.20 eV → nex = 0.44/0.68/1.07×10²¹, G ∝ σ^0.6) —
the conversion at Houston anticrossings survives small σ, so σ alone does not
close the gap to the bracket.

## Absolute-rate benchmark (`rate_benchmark.py`)

`python3 rate_benchmark.py` (after the on/off workflow) integrates the
Chefonov-style rate equation dn/dt = W_Keldysh(E(t)) on the same pulse and
overlays the SBE curves (→ `rate_benchmark.png`). At this field the two
Keldysh curves only *bracket* the physics (A₀ = 0.37 a.u. sweeps ~60 % of the
half-BZ — far outside the parabolic two-band regime), but the comparison makes
the absolute scale visible: treat SBE `nex` magnitudes at γ ≲ 1 as **upper
estimates** until benchmarked on your grid (`nk` convergence, σ_E, ν_sat).
The distributions (which k, which band, cooling) are robust; the absolute
yield is the calibration-sensitive part.

## Adding impact ionization (+ phonon-assisted II)

Impact ionization rides the **same ring** and is verified to stay exactly
electron-number conserving on this frozen window (`sbe_ii_phassist` enables the
phonon-assisted branch). In `Si_frozen_phonon_rt.inp` set both:

```
  yn_sbe_impact_ionization = 'y'
  yn_sbe_ii_holes          = 'y'     ! must be paired with impact_ionization
  sbe_ii_phassist          = 1.0d0
```

Electrons stay **8.000** (checked at I=10¹¹ and 10¹² W/cm², with Auger on, and
with the FULL channel set eph+ac+II+holes+phassist+Auger+eeh+Coulomb at
dt=10 a.u. — 8.000 at every step, no negative nelec/nhole).

## Notes / caveats

- **Basis convergence.** This field (I=10¹¹ W/cm², Keldysh γ≈1) is
  basis-converged: `nex` changes < 4 % from `nstate` 12 → 16. The built-in
  VG basis-edge monitor (`P_top` warning on stderr) stays quiet. **Deeper
  tunneling** (e.g. I=10¹² at this ω) sweeps virtual population across many
  bands — the monitor will warn and you must raise `nstate` (frozen core keeps
  that cheap, since only the active window is dissipated).
- **Grid.** The ring is O(nk³)/step; 4³ is the fast default. The indirect
  X-valley is only coarsely sampled at 4³ — raise `num_kgrid` (and MPI ranks
  dividing nk) to resolve the valley shell.
- **Trace exactness.** The frozen-window ring is trace-exact (`wiki/00`, the
  `ring_apply_dpop` fix): electron number is conserved to machine precision at
  every `dt`, so any residual `nex` is real carrier generation, not a leak.
