# Supercomputer runs — Lomonosov-2 (SLURM, OpenMPI + gfortran)

Shared reference for running SALMON-TROUT on the maintainer's HPC. This is the
page we talk through for every cluster run: the environment, the job recipe, how
to tell a healthy run from a broken one, and the failure catalog we have already
worked through. Keep it current — add each new machine, each new gotcha.

> **Two machines, one code.** The mini-cluster uses **Intel oneAPI (ifort +
> Intel MPI + MKL)**; Lomonosov-2 uses **OpenMPI 4.0.5 + gfortran (gcc-9.1) +
> MKL**. All MPI code in TROUT is compiler-agnostic (plain `MPI_Bcast` /
> `MPI_Allreduce` / `MPI_Barrier`), so a run that is correct on one must be
> correct on the other. When it is not, the cause is the **environment or a
> stale binary**, not the physics — see the failure catalog.

---

## 0. The golden rules (read before every run)

1. **Rebuild clean after every `git pull`/merge.** `rm -rf build/` (or `make
   clean`) then reconfigure + build. Fortran `.mod`/`.o` staleness has already
   cost us a full debugging cycle (a "0 active bands" result that was
   *impossible* from the merged source — it was an un-rebuilt `trout2`).
2. **SALMON reads the input from STDIN**, not a flag: `mpirun ./trout2 <
   input.inp`. Only rank 0 reads stdin (it copies the namelist to
   `.namelist.tmp`, which every rank then reads).
3. **Run from shared scratch.** Rank 0 writes `.namelist.tmp` and the GS
   `*_k/_eigen/_tm.data`; every node reads them back. The working directory must
   be visible to all nodes (Lustre scratch).
4. **Check the run banner** (Section 4) before trusting any output.

---

## 1. Machine profile — Lomonosov-2

| | |
|---|---|
| Scheduler | SLURM (`sbatch`, `srun`/`mpirun`) |
| Compute node | 14 cores (set `--ntasks-per-node=14`) |
| Compiler | gfortran, **gcc-9.1** (`/opt/software/gcc-9.1`) |
| MPI | **OpenMPI 4.0.5** (`/opt/mpi/openmpi-4.0.5-gcc`) |
| BLAS/LAPACK | Intel **MKL** 2019.5 (GNU threading layer) |
| Binary | built as `trout2` in the run dir |

Backtraces are gfortran-mangled (`__bloch_solver_ssbe_MOD_...`) — that is how you
confirm a crash log is from the Lomonosov build, not the Intel one.

---

## 2. Environment (the working `sbatch` preamble)

This is the environment block that runs correctly on Lomonosov-2 (from the
maintainer's job script). Keep it in the sbatch, above the `mpirun` lines:

```bash
export LC_ALL=C
export LANG=C
export MKL_THREADING_LAYER=GNU
export PATH="/opt/software/gcc-9.1/bin:/opt/mpi/openmpi-4.0.5-gcc/bin:$PATH"

export MKLROOT=/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl
export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/software/gcc-9.1/lib64"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/mpi/openmpi-4.0.5-gcc/lib:/opt/mpi/openmpi-4.0.5-gcc/lib64"

export MKL_DYNAMIC=FALSE
export OMP_DYNAMIC=FALSE
export OMPI_MCA_opal_warn_on_missing_libcuda=0   # no GPU on these nodes
export OMPI_MCA_opal_cuda_support=0
export MKL_DEBUG_CPU_TYPE=5                        # AVX2 codepath on the Westmere-labelled build
export MKL_CBWR=COMPATIBLE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=1                           # MPI does the parallelism; keep MKL serial
export MKL_DISABLE_FAST_MM=1
ulimit -s unlimited                                # large automatic arrays live on the stack
```

Notes:
- `MKL_NUM_THREADS=1` + one MPI rank per core is the safe default; hybrid
  MPI×OpenMP is possible (the ring gathers are OpenMP-parallel) via
  `--cpus-per-task > 1` and letting `OMP_NUM_THREADS` follow it.
- `ulimit -s unlimited` matters: several SBE routines use automatic arrays sized
  by `nstate`/`nk`.

---

## 3. Job recipe

The ready script is
[`samples/exercise_x12_Si_frozen_phonon_indirect/run_9x9x9_lomonosov.sbatch`](../samples/exercise_x12_Si_frozen_phonon_indirect/run_9x9x9_lomonosov.sbatch).
Skeleton:

```bash
#SBATCH --partition=compute
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00

# ... environment block from Section 2 ...
cd "${SLURM_SUBMIT_DIR}"

# (1) ground state — per-k independent, cheap; one rank is plenty
mpirun ./trout2 < Si_prim_epm_gs_9x9x9.inp > gs_9.log

# (2) real time, eph ON
mpirun ./trout2 < Si_frozen_phonon_rt_9x9x9.inp > on_9.log
cp Si_prim_sbe_nex.data on_nex.data; cp Si_prim_sbe_nex_k_real.data on_nkr.data

# (3) real time, eph OFF (reference) — sed-flip the two eph switches
sed "s/yn_sbe_eph               = 'y'/yn_sbe_eph               = 'n'/;\
     s/yn_sbe_eph_acoustic      = 'y'/yn_sbe_eph_acoustic      = 'n'/" \
    Si_frozen_phonon_rt_9x9x9.inp > off_9.inp
mpirun ./trout2 < off_9.inp > off_9.log
cp Si_prim_sbe_nex.data off_nex.data; cp Si_prim_sbe_nex_k_real.data off_nkr.data
```

The GS step **must finish before** the RT step (RT reads its `*_k/_eigen/_tm.data`).
Two SBE runs in the same dir overwrite each other's `Si_prim_sbe_*.data` — the
`cp` between steps captures each; for parallel field points give each its own
subdir.

### MPI rank count
`nproc` from `--nodes × --ntasks-per-node`. k-points are split across ranks
**evenly (±1)** by `split_num`, so `nproc` need **not** divide `nk`.
- Use **`nproc ≤ nk`** (more ranks than k-points idles ranks).
- For perfect balance pick a divisor: `nk = 9³ = 729 = 3⁶` → **27 / 81 / 243**.
- The inter-k ring does `nproc` communication hops/step, so 27–81 is the sweet
  spot at 9³; huge `nproc` trades compute for communication.

#### Recommended node counts for 9³ (14-core nodes)

Since `14·N` never divides `729 = 3⁶` (no factors of 2 or 7), cap `mpirun -np`
to a **divisor** for perfect ±0 k-balance (a few idle cores is fine):

| `mpirun -np` | nodes (14 cores) | k/rank | when |
|---|---|---|---|
| 27 | 2 | 27 | economical, longer walltime |
| **81** | **6** | **9** | **balanced default** |
| 243 | 18 | 3 | faster turnaround (Coulomb/ring-heavy) |

Above ~52 nodes (`729` ranks = 1 k/rank) is the ceiling — the `O(nproc)` ring
hops dominate well before that. Running all `14·N` ranks (e.g. 84 on 6 nodes)
also works, with a minor ±1 imbalance; `-np 81` is just cleaner.

#### Coulomb HF cost — can it be off?

`yn_sbe_coulomb='y'` (Σ^HF dynamic-gap renormalization) is an `O(nk²)`-per-step
all-gather (frozen over dt) — the single most expensive term at 9³. It **can be
turned off** (`yn_sbe_coulomb='n'`):
- **Physically safe at LOW carrier density** (weak fields / small `nex`, e.g. a
  ~100 kV/cm THz pulse): the density-driven gap shift is negligible → drop it
  for a large speedup.
- **Keep it ON for high density** (strong fields, near-degenerate plasma), where
  the bandgap renormalization matters — or use the cheap stand-in
  `yn_sbe_bgr_threshold='y'` (II threshold only).
- Turning Coulomb off does **not** disable the ring (e-ph/II/Auger ride
  `yn_sbe_superres`). Those gathers are themselves `O(nk²)`/step, so the node
  guidance above still applies unless you also cut channels.

---

## 3b. Checkpoint & restart (long / walltime-limited runs)

A production RT run can outlast the queue's walltime. The B4 checkpoint lets a
run resume the physics state (the density matrix) from where it stopped instead
of starting over.

### The two inputs (`&sbe`)

| input | default | meaning |
|---|---|---|
| `sbe_checkpoint_step` | `0` (off) | write a checkpoint **every N time steps**. Each rank streams its own state to `SYSNAME_sbe_ckpt_r<rank>.bin`. |
| `yn_sbe_checkpoint_restart` | `'n'` | `'y'` = at startup, **resume** from the checkpoint set instead of starting at t=0. |

### What is in a checkpoint

Per MPI rank, one binary stream file `SYSNAME_sbe_ckpt_rNNNNN.bin`
(`NNNNN` = 5-digit rank), **overwritten in place** each checkpoint (only the
latest survives — it is not a history). It holds:
- the step index `it` and the accumulated `energy`,
- the per-channel ledgers (`led_dn`, `led_de`),
- this rank's local density-matrix block `rho(:, :, ik_min:ik_max)`,
- the Houston branch `X_branch(:, ik_min:ik_max)` (if the run tracks it).

The **field is NOT saved** — it is recomputed deterministically from the input
each step, so nothing else is needed to continue exactly.

### How to use it — the two-phase resubmit

1. **First run:** set `sbe_checkpoint_step = 500` (say) and keep
   `yn_sbe_checkpoint_restart = 'n'`. The run writes `*_sbe_ckpt_r*.bin` every
   500 steps.
2. **If it is killed** (walltime, node failure): resubmit the **same input**
   with `yn_sbe_checkpoint_restart = 'y'` (leave `sbe_checkpoint_step` on so it
   keeps checkpointing). On startup rank 0 prints:
   ```
   # B4: resumed from checkpoint, continuing at step <it0>
   ```
   and the run continues from the last checkpointed step to `nt`.

In the x12 sbatch this is a one-line `sed` on the RT input, e.g.:
```bash
# resume variant of the RT input
sed "s/yn_sbe_checkpoint_restart = 'n'/yn_sbe_checkpoint_restart = 'y'/" \
    Si_frozen_phonon_rt_9x9x9.inp > Si_..._rt_resume.inp
mpirun ./trout2 < Si_..._rt_resume.inp > on_9_resume.log
```
(If `sbe_checkpoint_step` is not already in your input, add both lines to `&sbe`;
they are absent by default.)

### Hard constraints (all must match between the original run and the resume)

- **Same number of MPI ranks.** Each `.bin` holds only that rank's `ik_min:ik_max`
  k-block, and the partition depends on `nproc`. A different rank count reads the
  wrong slab → wrong physics (or a read error). This is the #1 rule.
- **Same input** — grid (`num_kgrid`), `nstate`, `dt`, the field, and the channel
  set. The field is recomputed, so any field/`dt`/grid change desynchronises the
  resume.
- **Same working directory** (shared scratch): the `.bin` files must be found
  where they were written. `yn_sbe_checkpoint_restart='y'` with no checkpoint
  file present is a hard `error stop` ("checkpoint file missing").
- The `.bin` is a raw per-rank memory dump — **not portable** across a different
  machine/compiler/endianness. Resume on the same build.

### ⚠️ Output-file caveat (important for analysis)

On restart the time-series outputs (`SYSNAME_sbe_rt.data`, `_sbe_nex.data`,
`_sbe_channels.data`, …) are **reopened for writing (truncated)** and then
written only for the `t=0` block and steps `it0 … nt`. **The rows for steps
`1 … it0-1` are lost from those files** even though the physics state resumed
correctly. If you need the *complete* `nex(t)` / current time series across a
restart, **rename or copy the `*.data` files before resubmitting**, then stitch
the pre- and post-resume files together (drop the duplicated `t=0` row). The
final-state quantities and any `_k` snapshot at/after the resume are unaffected.

### Verifying a resume
- The banner line `# B4: resumed from checkpoint, continuing at step <it0>` with
  the expected `it0`.
- `electrons = 8.000` from the first resumed step (the density matrix carried the
  correct trace across the restart).

---

## 3c. Time step (dt) — by material

The binding constraint on `dt` is **not the field** (THz pulses are slow) but the
**interband coherence at the direct gap** (and the active-band energy spread).
The velocity-gauge CF4 propagator diagonalises the instantaneous Hamiltonian, so
it is *exact* for a static H; the per-step error comes from the field rotating
the Houston basis over `dt`. Resolve the fastest interband beat:

> **dt ≲ (2πħ / E_gap^direct) / 10**   (ħ = 0.6582 eV·fs)

| material | direct gap E₀ (Γ) | interband period 2πħ/E₀ | recommended `dt` (moderate field) |
|---|---|---|---|
| GaAs | 1.42 eV | 2.9 fs | ~0.15–0.20 fs (6–8 a.u.) |
| CdS | 2.55 eV | 1.6 fs | ~0.12 fs (5 a.u.) |
| **Si** | ~3.4 eV (E₀′, Γ₂₅′→Γ₁₅) | 1.2 fs | **~0.10 fs (4 a.u.)** |
| graphene | 0 (gapless Dirac) | bandwidth-limited | ~0.05 fs (2 a.u.) or finer |

**Validated numerically for Si** (weak 100 kV/cm THz, clean VG, 4³, all 24 bands,
to the field peak): the **current** converges already at `dt = 0.25 fs` (~1 %),
but the **absorbed energy / nex** is under-resolved at 0.25 fs (it comes out with
the *wrong sign*) and converges only at `dt ≈ 0.1 fs`. So:
- **Carrier / absorption studies (`nex`):** use `dt ≈ 0.1 fs` for Si (scale by
  the gap for other materials, per the table).
- **Linear current / conductivity only:** `dt = 0.25 fs` is acceptable.

Caveats:
- **Stronger fields drive faster dynamics** → reduce `dt` proportionally (halve
  near tunneling / ≥ few MV/cm).
- **graphene is gapless** → always the smallest `dt` (no gap protection; the
  highest active band + the field-driven cone crossing set it).
- **Dissipation rarely binds `dt`** (`ν_sat⁻¹ ~ 8 fs ≫ dt`); the ring CPTP
  limiter warns and auto-scales if a step's scattering flux is too large.
- The values above are **starting points** — confirm with a **2× `dt` halving**
  on your field until `nex`/current stop changing (as done for Si above).

---

## 3d. Performance tuning (measured levers)

Where the wall-clock time goes, and which knobs actually move it. Measured on
**Si, 4³, weak 100 kV/cm THz** (`100.txt`), 1200 steps, `dt = 0.25 fs`,
**Coulomb off**:

| channels | nstate | s/step | rel. speed | `nex` [cm⁻³] |
|---|---|---|---|---|
| full (eph+ac+II+Auger+holes+eeh) | 24 | 0.216 | 1.0× | 3.48×10²¹ |
| full | 16 | 0.200 | 1.08× | 3.47×10²¹ |
| full | 12 | 0.196 | 1.10× | 3.46×10²¹ |
| **eph+ac only** | 24 | 0.044 | **4.9×** | 1.83×10²¹ |
| **eph+ac only** | 16 | 0.029 | **7.4×** | 1.83×10²¹ |

### Cost model (what scales how)
- **Unitary** (velocity-gauge CF4, `cf4_unitary_step`) runs on the **full
  `nstate` basis**: `O(nstate³ · nk)` per step. (Frozen-core shrinks only the
  *dissipators*, not the unitary — see §4.)
- **Ring dissipators** (e-ph, II, Auger, hole-II) draw the partner from the
  whole BZ: `O(nk² · n_active²)` per step, one pass **per channel**.
- So at a fixed grid the ring share **grows as `nk²`** while the unitary grows as
  `nk` — the ring dominates more and more with grid size.

### The two levers

**1. Channel set — the dominant cost, but NOT free here.** Dropping
`II + Auger + ii_holes + eeh` (keeping `eph + acoustic`) is **~80 % of the
per-step time at 4³** (0.216 → 0.044 s/step, **4.9×**) — and **larger at 9³**,
where the ring is a bigger share (ring ∝ nk²). **But at this field these
channels are not inert:** `nex` drops **3.48 → 1.83×10²¹ (~1.9×)** — they carry
about half the generation. So this is a **physics decision**, not a free
speedup:
- For a **100 kV/cm THz** field, impact ionization (threshold field ~MV/cm) is
  physically questionable; the ~2× it adds here is plausibly the long
  multi-cycle pulse + the `ν_sat` calibration inflating it (the "absolute yield
  is an upper estimate" caveat, wiki/00 / x12). **Scrutinise whether II/Auger
  belong at your field.** If the study is **phonon-assisted generation only**,
  `eph + ac` gives the 5–7× speedup honestly.
- If II/Auger *are* part of the physics, keep them and pay the cost (or cut the
  grid / dt instead).

**2. `nstate` — physically ~free, modest speed gain.** `nex` is already
**basis-converged at `nstate = 12`** (3.46 vs 3.48×10²¹ at 24) for this weak
field, so lowering it (with a matching GS) costs no physics — it only removes VG
headroom a *stronger* field would need. But the **speedup is modest with full
channels** (24 → 12 = only ~10 %, ring-dominated) and larger only once the ring
channels are off (eph-only 24 → 16 = ~35 %). At 9³ the ring dominates further,
so `nstate` helps full-channel runs even less. Verify `nex(nstate)` and keep it
as low as basis-sufficiency allows, but don't expect big gains while II/Auger
are on.

### Ranking for a 9³ production run (Coulomb already off)
1. **Channel set** — decide II/Auger on physical grounds (biggest lever; ring ∝ nk²).
2. **Grid `nk`** — ring ∝ nk²; 7³ vs 9³ ≈ 4.5× on the ring (X-valley coarser).
3. **`dt`** — the largest your observable allows (§3c: ~0.1 fs for `nex`, 0.25 fs current-only).
4. **`nstate`** — lower to basis-sufficiency (free physics, modest speed).
5. **`out_projection_k_step`** — k-resolved output is I/O ∝ nk; write it less often.

---

## 4. Reading a healthy run banner

A correct RT run prints (rank 0 stdout → `on_9.log`):

```
   k-points =    729,   bands =   16,   active =    7
  Fermi energy (eV)      =      11.8901 eV
  Gamma reference: k-point    365 at (  0.0000  0.0000  0.0000) [reduced]
  n_active_bands         =    7 /   16
 ...
   <step>  <t[fs]>  <Jx Jy Jz>   <electrons>   <energy>
```
and ends with `end SALMON`.

**The five things to verify:**
1. `Gamma reference: k-point N at (0,0,0)` — the frozen-core window is anchored
   at Γ. On an **odd** grid (9³) the coords must be exactly `(0,0,0)`; on an even
   grid (4³) it is the nearest-to-Γ point (small non-zero coords, expected).
2. `n_active_bands = N / M` with `N ≥ 2`. For Si at Γ with the `−3/+6 eV` window
   it is **7/16** (VBM triplet Γ₂₅′ + CB manifold). A `0/16` or a `< 2` count now
   aborts with a diagnostic — if you ever see `0/16` *proceed*, the binary is
   pre-fix (stale).
3. **`electrons = 8.000` at every step** (Si, nelec=8). The frozen-window ring is
   trace-exact; drift means a real bug.
4. `end SALMON` at the tail (no `SIGSEGV`, no `N processes killed`).
5. The **VG basis-edge** line on stderr (Section 6) — physics convergence, not a
   crash, but it gates whether the *absolute* yield is trustworthy.

---

## 5. Failure catalog (what we hit, and the fix)

| Symptom in the log | Root cause | Fix |
|---|---|---|
| `SIGSEGV` in `compute_coulomb_selfenergy_ring` at step 1, one rank | frozen-core count/mask broadcast disagreed on a non-synchronized start → `active_idx` out of bounds | PR #93: derive count from the mask, cross-rank consistency guard |
| `n_active_bands = 24/24` but per-band flags show only 8 active | same broadcast disagreement | PR #93 |
| intermittent field-file read errors / segfault reading `file_input1` | every rank opened the same text file on shared FS; a lagging/partial file gave a short read | PR #93: rank 0 reads once, `comm_bcast(Ac_ext_t)`; iostat guards |
| `n_active_bands = 0/16`, correct band energies printed, then `5 processes killed` | **stale binary** (impossible from the merged source) | **clean rebuild** (`rm -rf build/`); PR #94 also removes the mask broadcast and anchors at Γ so it cannot recur silently |
| window keyed off the wrong k (k=1 is the grid corner, not Γ, on 9³) | window was built from `gs%eigen(:,1)` | PR #94: find Γ = nearest-to-origin k, broadcast its energies, build the window from that (identical on all ranks) |

The throughline: on a distributed start these were **compiler-agnostic** MPI
issues, all fixed with standard collectives. If a *new* multi-node-only symptom
appears, first rule out the stale binary, then check whether the run banner's Γ
reference and `n_active_bands` agree with a serial run of the same input.

---

## 6. Physics convergence — the VG basis edge (`P_top`)

Separate from correctness: the RT monitor prints, on stderr,
```
WARNING: VG basis edge reached -- P_top = 1.16E-02 > 1.0E-03 (top band 8).
Increase N_b (nstate) and re-check convergence.
```
when population reaches the **top of the active window**. It means the band
budget is too small for this field, so the **absolute** `nex` is an *upper
estimate* (carriers pile at the window edge instead of dispersing higher). See
`wiki/06_vg_basis_nb_convergence.md`.

- The **on/off contrast** and the **distributions** (which k, which band) stay
  robust even with the warning — both runs share the same basis.
- For a converged **absolute** yield (needed for the Keldysh-bracket verdict):
  raise `nstate` (e.g. 24) and/or widen `frozen_free_threshold_ev`, run 2–3
  values, and confirm `nex` plateaus and the `P_top` warning clears.

---

## 7. Run log (this project)

| Date | Grid | Ranks | Result |
|---|---|---|---|
| 2026-07-21 | 9³ (729 k) | ~84 (6×14) | ✅ GS→ON→OFF complete; Γ ref = k365 (0,0,0); `n_active_bands = 7/16`; electrons = 8.000; **nex(eph ON)/nex(OFF) = ×231** (phonon-assisted indirect Γ→X). ⚠️ `P_top ≈ 1.2e-2` throughout → absolute nex (1.6e22 cm⁻³) is an upper estimate; needs an N_b-convergence sweep before the bracket verdict. |

Append each production run here (grid, ranks, active count, electron
conservation, the physics number, and any warning) so the cluster history lives
with the code.
