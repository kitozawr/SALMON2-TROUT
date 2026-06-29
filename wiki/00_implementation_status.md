# Implementation Status — live progress tracker

**Read this first on every resume.** It records what is done, what is next, key decisions, and the test inventory. Update it in the same commit as any code change.

Roadmap: (1) Si + nonlocal-super-compute (Parts A–G, all done, merged via PR #44 into `develop-2.0.0`); (2) new materials — **GaAs, Si/Si_cb, wurtzite CdS, monolayer graphene** — each a validated Python EPM reference + cited CPTP dissipation channels. The maintainer drives one bounded increment per session ("continue"). Test grid: **4×4×4 (or 5×5×5), scalar (no spinor)**. Current test count: **14/14** (`python3 tests/run_all.py`).

---

## ⭐ CURRENT WORK — PRIMITIVE-CELL (non-orthogonal, NO folding) pivot — branch `claude/sbe-nonorthogonal`

**Read this block FIRST on resume — it is the live frontier.** Off merged `develop-2.0.0` (PR #46 merged the plotter/real-carrier/intra-current work). All runs are in the session scratchpad `…/scratchpad/{gaas_prim,gaas_prim_odd,gaas_prim_so,gaas_prim_so7,si_prim,cds_prim}`; the SALMON binary is `build/salmon`.

**WHY the pivot:** the "L over-populates Γ" pattern was proven a **band-folding artifact**. Folded cubic gives L/Γ≈760 (anti-Zener); the **primitive FCC** gives Γ≫L≫X (correct Kane ordering, L/Γ≈0.01, matches Zener P). Removing folding swings L/Γ by 10⁵. So we now run the SBE on the **primitive cell directly** (non-orthogonal, no cosets, no unfold map, no sublattice projection — the primitive cell IS the irreducible problem).

**DONE & verified (this branch):**
- `epm_gaas_primitive.py` — FCC 2-atom non-orthogonal (plane waves = all-same-parity G's = BCC reciprocal). Reuses `epm_gaas_reference` H/momentum. Scalar: Γ 1.386 / L 2.677 / X 3.943 eV. **Spinor** via `INCLUDE_SPIN_ORBIT=True` (reuses ref SO machinery, mu calibrated on the primitive basis): Δ₀=0.341 eV, gap 1.273 eV, Kramers doublets, spin-split zero along ⟨100⟩/⟨111⟩ (Dresselhaus). Scalar path byte-identical when SO off. CLI: `gap`/`gs`/`bandpath`.
- `epm_si_primitive.py` — reconfigures `epm_gaas_primitive` for Si (V^A=0, Kunikiyo, a=10.26). **Cutoff RAISED to 27 Ry** (Si Δ-valley camel-back needs it; 11.1 mislands CBM at X) → indirect 1.059 eV @ 0.850·X. Cutoff only sizes GS, not the SBE.
- `epm_cds_primitive.py` — **4-atom HEXAGONAL wurtzite** primitive (genuinely new non-orthogonal geometry, a1/a2 at 120°). Reuses `epm_wurtzite_cds` geometry/H/BC1967 form factors + the proven non-orthogonal GS writer. Gap 2.547 eV vs BC1967 2.58. nelec=16, nstate=16.
- `plot_sbe_results.py` — **(a)** Cartesian-BZ heatmap (`_cartesian_bz_grid` un-shears the triclinic grid into a Wigner-Seitz Cartesian volume; reads `# b1/# b2/# b3` from the GS `k.data` header via `_read_bmatrix`/`_bmatrix_for`; only the non-orthogonal datasets get it, cubic legacy untouched). **(b)** `plot_primitive_spectral` = per-frame **movie** into `spectral_frames/` (thin skeleton, fixed colour scale, for ffmpeg). **(c)** spinor auto-detect + spin-splitting plot already worked.
- `src/ssbe/gs_info_ssbe.f90 read_k_data` — now skips **any** number of `#`/blank header lines (was hard-coded 5) so the EPM can write the reciprocal vectors into `k.data`. **Backward-compatible; GaAs reference re-verified intact** (8³ rerun: Γ dominates ~128×). This is the ONLY SBE/Fortran change on the branch.
- SBE runs verified end-to-end on the primitive cells: GaAs scalar (8³, **9³ odd → explicit Γ**), GaAs spinor (4³, **7³ odd**), Si (8³ + **dissipators+super-mode**: electrons=8.0000 exact CPTP, energy relaxes 12.5% post-pulse), CdS (7×7×5).

**KEY GOTCHAS (do not relitigate):**
- **`yn_sbe_spinor='y'` is MANDATORY in `&sbe` for a spinor GS dataset.** Default 'n' makes the solver read the 16-band spinor GS as scalar (`nb_vb=nelec/2=4`, `ib_lcb=5` = a *valence* band, `occ_max=2`) → garbage population (the bogus fc=1.59). With the flag: `nb_vb=nelec=8`, `ib_lcb=9`, `occ_max=1` → physical (fc≤1).
- **ODD k-grids sample Γ explicitly** (even grids straddle it at e.g. 0.1875); use odd (7³,9³) for a clean Γ spot. The odd grid resolves the sharp near-Γ excitation pocket the even grid misses (so totals rise but converge).
- **CdS with a c-axis (E∥c) field:** the top valence Γ9→CB transition is **dipole-forbidden at Γ** (only Γ7 couples to E∥c) → excitation is suppressed at Γ and peaks off-axis. This is real wurtzite selection-rule physics, NOT an artifact. (Use E⊥c to populate Γ.)
- **MPI load balance:** `split_num` (util_ssbe) gives a contiguous partition, max imbalance = 1 k-point/rank (`P/nk`). Non-orthogonality is balance-neutral (per-k cost is geometry-independent); odd cubes just factor awkwardly vs power-of-2 ranks. Super-mode ring sizes buffers to `maxn` (comm balanced); `n_active_bands` is global so dissipator cost is data-independent (no carrier-pileup imbalance). Pick nproc | nk.

**STANDING TODOs from the maintainer (2026-06-29, verbatim intent — sync target of this block):**
0. **Keep THIS wiki block + tasks in sync** so a fresh session resumes after a context/limit cutover. (readme+wiki = long-term memory.)
1. ✅ **DONE.** All four material EPM Hamiltonians are now vectorized. Finding: GaAs/Si and graphene were O(npw²) Python double-loops (only CdS used einsum). Vectorized `build_hamiltonian_sc` (GaAs/Si, pairwise dG-broadcast + per-shell form-factor fill) and graphene `build_hamiltonian` — both verified **byte-identical** to the loop (GaAs/Si 0.0e0; graphene 3.5e-15) with all EPM tests passing (Dirac cone + 4-coset/2-coset folding intact).
2. Add **Fortran EPM for CdS and graphene** — *without folding* (primitive cell) it should work.
3. `--spectral` must colour the **FOUR** levels (VB-1,VB,CB1,CB2) on the bandmap, not just one — needs the SBE to emit the 4-level primitive populations (currently only LCB).
4. ⚠️ **INVESTIGATED — real limitation found, gate needed.** `apply_eph_relaxation` (bloch_solver ~1784) is **k-local**: for each carrier band `ia` at a k-point it finds the partner **band** `ib` at the **SAME k** whose energy matches `evals(ia) ∓ ħω` (lines 1824-1833) — it uses band energies, NOT hardcoded valley coordinates, so there is nothing coordinate-wise to "fix". BUT the physics only closes in the **folded** picture: there the L/X/Δ valleys fold onto higher *bands at the same supercell k*, so the intra-k band-to-band transfer *is* intervalley scattering. In the **primitive** cell the valleys live at **different k-points**, so the same-k energy-matched partner is NOT the other valley — true intervalley e-ph is unrepresented (the search silently relaxes band-to-band at fixed k). **Consequence:** CdS **Fröhlich polar-optical** (intra-valley, q≈0) is fine on the primitive cell; **Si/GaAs intervalley** (g/f, Γ–L–X) is NOT — it needs an **inter-k** energy-matched final-state search (over k AND bands), analogous to the nonlocal-II ring (the "nonlocal e-ph" the C5 note already lists as pending). **Action (maintainer principle — gate the unphysical):** on a primitive (non-orthogonal / `n_coset=1`) dataset, intervalley e-ph modes should be gated off (or e-ph restricted to the polar-optical mode) until the inter-k search exists. The earlier Si-primitive dissipator run DID conserve trace and relax energy, but via intra-k band-to-band transitions, not true intervalley transfer.

   **e-ph IS material-aware from `&epm epm_material` — verified per material** (the phonon table comes from `get_material_params(epm_material)`, NOT hard-wired): Si → **6 intervalley g/f modes** {10,19,63,19,51,57 meV}; CdS → **1 Fröhlich polar mode** {38 meV}; GaAs → polar LO (mode 1) + 5 intervalley. **Per-material primitive-cell e-ph validity:**
   | Material (`&epm`) | e-ph modes | Primitive-cell validity |
   |---|---|---|
   | **CdS** | 1 polar-optical (Fröhlich, q≈0) | ✅ **valid** (intra-valley; relaxes within band/valley) |
   | **GaAs** | polar LO + 5 intervalley | ⚠️ polar mode valid; the 5 intervalley modes need inter-k |
   | **Si / Si_cb** | 6 intervalley g/f | ⚠️ **all intervalley** → intra-k search is the *folded* picture only; needs inter-k on the primitive cell |

   So **CdS dissipators run correctly on the primitive cell** (its only e-ph mode is intra-valley); Si/GaAs intervalley e-ph on the primitive cell awaits the inter-k final-state search. Maintainer instruction (2026-06-29): "e-ph should use the material from `&epm`" — confirmed it does.

   **🟢 MAINTAINER DECISION (2026-06-29, corrected): inter-k e-ph goes THROUGH THE RING. "If the ring (`yn_sbe_superres`) is on, inter-k goes through it."** So the activation is gated on the **systolic ring** (Part D), NOT on `yn_sbe_eph` alone:
   - `yn_sbe_eph='y'` **without** the ring → stays **k-local (intra-k)** — the current behaviour, valid for folded cells and for polar-optical (CdS). Byte-unchanged.
   - `yn_sbe_eph='y'` **with** `yn_sbe_superres='y'` (ring on) → the e-ph inter-k transfer **rides the same systolic ring** that already circulates the k-blocks for the Coulomb HF: as each q-block transits, run the energy-matched e-ph final-state search of the local carriers against the transiting block and accumulate the inter-k population transfer (amp_damp), exactly alongside `form_sigma_hf_ring`. This is the pending "nonlocal e-ph (C5/**D**)" — and **D = the ring**, so the design is: add an e-ph accumulation pass inside the ring transit (`compute_coulomb_selfenergy_ring` hop loop — the code comment already says "extra nonlocal accumulators (II, e-ph, e-e) can be added here without new communication"). 🚧 NEXT.

   **CPTP algorithm (the correctness key — implement exactly this):** gather `eval(nba,nk)` once per step (cheap; the ring already circulates the population blocks `transit`). In the hop loop, for each local carrier `(ik,a)` and each transiting `(iq,b)`, form the symmetric energy-matched pair rate `w = nu(eps)·w_mode·gaussian(|E(ik,a)−E(iq,b)|∓ħω)·Pauli`, and accumulate the **net** diagonal change `Δf(ik,a) += [in from (iq,b)] − [out to (iq,b)]`. Because rank-A's "(ik→iq) out" uses the **same** pairwise `w` as rank-B's "(iq→ik) in" (detailed balance, emission↔absorption with N_B/(N_B+1)), the transfers match across ranks ⇒ **global trace conserved without any bidirectional deposit** — each rank only ever writes its OWN local `Δf`. Apply `Δf` locally after the ring (clamp to [0,occ_max]) plus the source-coherence damping at rate `Γ_out(ik,a)` (local). Intra-k (same-k bands) is included automatically as the `iq==ik` term, so the ring path SUBSUMES and replaces the intra-k call when the ring is on. Needs a CPTP test (trace-conserving, PSD, γ=0 identity) mirroring `test_eph_cptp.f90` but for the inter-k pair transfer.
5. Keep working examples + docs updated (long-term memory).
6. ✅ **VERIFIED (already implemented).** `plot_conductivity` = σ(ω)=J(ω)·conj(E)/(|E|²+floor) (`_sigma_ratio`), Hann-windowed, **default 0–4 THz**, Re+Im. `plot_conductivity_stft` = Re σ(ω,t) 2-D map, **hop defaults to 1 sample → N−1 of N overlap** (the requested smoothness), 0–4 THz; effective hop only rises for the render cap (max_cols) while keeping overlap maximal. Note: needs a ps-scale run for true THz resolution (short test runs warn "trace too short for 0–4 THz").
7. Find the un-merged branch with **deeppseudodot** copied in; set up DFT(primitive, in-salmon, Si example)→deep-EPM coefficient fitting; drop a ready Si example in `samples/` with a DFT-compatibility layer (user runs the long calc; EPM-compat not yet).
8. Refresh the **SBE console output header/banner** (outdated).

**PRE-LIST material-pipeline tasks (predate the 8-item list, still open — order: map → physics → next material):**
- **A. GaAs (reference, 100% done):** scalar + spinor GS/SBE/maps/dissipators+super-mode all ✅.
- **B. Si (material 1) ✅:** GS (indirect 1.059 eV @ 0.85·X), SBE, maps, dissipators+super-mode (CPTP, 12.5% relax) all done. *Caveat: its e-ph is all-intervalley → awaits inter-k (TODO-4 decision).*
- **C. CdS (material 2) ✅:** GS (2.547 eV), SBE, maps, **dissipators+super-mode all done** — ran on the hexagonal primitive cell with polar e-ph + Coulomb + Auger (eeh forbidden; II skipped — uncited prefactor): **electrons=16.00000000 exact (CPTP)**, real LCB carriers excited (peak 0.69) and acted on by the channels. CdS e-ph is polar/intra-valley → valid on the primitive cell as-is.
- **D. graphene (material 3) ✅ (clean):** `epm_graphene_primitive.py` — 2-atom hexagonal primitive (non-orthogonal, 2D-in-vacuum, no folding), reuses `epm_graphene` + the non-orthogonal GS writer. GS gapless Dirac at K, **SBE ran with electrons=2.00000 exact (CPTP)**, in-plane current driven, and the Cartesian-BZ heatmap shows carriers localized at the **six K/K′ Dirac corners** (textbook). Maps (bandpath/snapshot/Cartesian) render. **Dissipators still pending** — the graphene e-ph/Auger registry entry `epm_material='graphene'` (E2g 196 meV, A1′ 160 meV, gapless carrier-multiplication, no-Kuhn-Zurek) is NOT yet added (folded-era item 3).
- **E. Examples + wiki HUGE update ⬜:** `samples/` recipes + wiki for the whole primitive-cell pipeline (all materials, spinor, Cartesian map, spectral movie, `yn_sbe_spinor='y'`, odd grids). Ongoing; wiki/00 kept live.

---

---

## Legend
✅ done & tested  · 🚧 in progress · ⬜ not started · 🔭 future / out of current scope

---

## Status by part

| Part | Item | State | Branch / PR | Notes |
|---|---|---|---|---|
| pre | CF4/Yoshida/Strang propagator, Kuhn-Zurek dephasing, frozen core | ✅ | merged | do NOT touch the integrator |
| pre | k-local impact ionization (Stobbe quartic), GaAs | ✅ | merged | Houston basis, HF-factorized amplitude damping |
| pre | Coulomb TDHF exchange Σ^HF (Golde-Kira-Meier-Koch) | ✅ | merged (#41/#43) | non-k-local, all-gather/step; δρ subtraction; Hermitian |
| pre | EPM GaAs + 4-fold FCC folding + unfold pipeline | ✅ | merged | `epm_gaas_reference.py`, `SYSNAME_unfold.data` |
| **A** | Silicon EPM (`epm_material='Si'`, V^A=0, Kunikiyo) | ✅ | #44 | gap 1.059 eV conv (Kunikiyo calc 1.068), CBM 0.86·2π/a |
| **B** | II fit-form switch (`sbe_ii_form`, `sbe_ii_exponent`) | ✅ | #44 | d**a, prefactor au_ev**a; GaAs a=4 unchanged |
| **E** | HF sublattice-block projection (`yn_sbe_hf_sublattice_proj`) | ✅ | #44 | proj_ij=Σ_s w_s(i)w_s(j) off-diag, diag kept; Hermitian; default 'y' |
| **C** | Nonlocal super-compute (`yn_sbe_superres`) | ✅ | #44 | C0-C8 all done; ring (D) + nl-II + nl-eph + e-e all on the super-mode path |
| C0 | Flags scaffolding (10 params, default OFF) | ✅ | #44 | namelist read+bcast+log; default run byte-for-byte unchanged |
| C-prim | Pure rate/search primitives module `sbe_superres_ssbe.f90` | ✅ | #44 | ν(ε), N_B, bins, Fröhlich asinh, II rate, BGR, amp-damp map, golden-rule prefactor + unit conv, thermal split, Si/GaAs tables |
| C2 | Houston-basis adiabatic populations reuse | ✅ | #44 | e-ph runs in the existing houston_dissipate ZHEEV basis (t2 = U†ρU) |
| C3 | Energy-bin final-state SEARCH (energy_partner_weights) | ✅ | #44 | windowed broadened-delta weights over candidates; deterministic, no MC; unit-tested |
| C4 | Nonlocal impact ionization (momentum exchange) | ✅ | #44 | valence partner + Pauli factors from BZ-averaged occupation (gathered/step); CPTP; end-to-end stable. Full momentum-resolved final states = refinement |
| C5 | e-ph population-relaxing Lindblad (k-local, full phonon table) | ✅ | #44 | Si 6 intervalley / GaAs LO+5; golden-rule weights D²/ħω (norm.), detailed-balance emis/abs, ν(ε) cap, Pauli-clamped, CPTP; trace=32 conserved. Nonlocal version pending (C4/D) |
| C6 | CPTP gate (amplitude-damping map test) | ✅ | #44 | trace, qubit positivity det≥0, transfer formulas, Hermiticity, γ=0 identity |
| C7 | BGR-gated II threshold | ✅ | #44 | running n(t) shifts E_th=E_th0−|K n^⅓| above gate; end-to-end stable |
| C8 | Dissipator sub-cycling | ✅ | #44 | m_sub from eph_numax·τ; II+e-ph split into m CPTP sub-steps; trace conserved at ν_sat=1e18. m=1 (unchanged) when e-ph off |
| **D** | Ring/pipeline MPI (systolic ring, one fused pass) | ✅ | #44 | Σ^HF via ring in super-mode; O(Nk/P) mem; ring==all-gather (6e-21), MPI 1==2 ranks (per-k 0.0). nl-II/eph/e-e accumulators slot in |
| **G** | Screening primitives (TF/Debye, Lindhard/RPA, LOPC) | ✅ | #44 | pure functions + GaAs/Si dielectric constants; unit-tested |
| **F** | Carrier-carrier (e-e/e-h) CPTP thermalization (`yn_sbe_eeh`) | ✅ | #44 | intra-k FD relaxation: conserves number AND energy exactly, EID coherence damping, CPTP; unit-tested + end-to-end. Inter-k momentum-resolved on ring = refinement |
| **MS** | Maxwell-SBE multiscale adaptation (`theory='multiscale_experiment'`) | ✅ | #44 | all A–G channels usable per-macropoint via the shared `init_sbe_bloch_solver`/`dt_evolve_bloch_cf4`; `sbe%icomm` set unconditionally (was Coulomb-only → MPI_COMM_NULL for BGR/nl-II-only); banner-print gated to one macropoint (`verbose`) |
| **VG** | Band budget: VG basis sufficiency & N_b convergence | ✅ | #44 | separate axis from PW cutoff, NOT cured by Houston basis; primitives + test (interlacing/shift formula) + runtime P_top monitor (warns on error channel, continues) |
| doc | Wiki pages 01–06 committed as long-term memory | ✅ | #44 | maintained per increment |

---

## Decisions log (gotchas that bit us / must not be re-litigated)
- **Si gap:** a 3-parameter local EPM gives ~1.06 eV converged = Kunikiyo's own calc (1.068), NOT the 1.12 eV experimental value. The spec's "30 meV from 1.12" is against experiment; treat **1.068 eV (Kunikiyo calc)** as the real target. CBM at 0.86·2π/a along ⟨100⟩.
- **Si diamond:** V^A ≡ 0 for all shells → the existing zincblende `VS·cos + i·VA·sin` machinery is reused verbatim with V^A=0 (purely real diamond). τ=(a/8)(1,1,1) for both.
- **Coulomb Σ^HF uses δρ = ρ − ρ₀** (ρ₀ = diag(occup)) so Σ(t=0)=0 and the equilibrium EPM gap is unchanged. Verified Hermitian; 1-rank vs 2-rank populations bitwise-identical.
- **nex density** = (excited e per cell)/V_cell; the /N_k is already in `calc_trace` via /Σ(kweight) (BZ average). Do NOT divide by N_k again.
- **Max ionization ceiling** (32 e per cubic cell) = 1.77e23 cm⁻³ (full valence inversion).
- **GitHub CI is serial-only by request** (compiler/syntax gate). MPI verified locally.
- **Commits are signed-off** (`git commit -s`), trailers Co-Authored-By + Claude-Session.
- **Part G = screening** (TF/Debye, static Lindhard/RPA = default, dynamic LOPC = GaAs-only); **Part F = carrier-carrier (e-e/e-h)** collision channel that USES G.
- **No HF double-counting:** carrier-carrier (F) is the correlation (2nd-Born/GW) self-energy — dissipative only. The static screened-exchange energy shift stays SOLELY in HF (Σ^HF). Do not add it twice.
- **Carrier-carrier invariants:** conserves Σf_k (number) AND ΣE_k f_k (energy) within the carrier subsystem (it thermalizes to a hot Fermi-Dirac, does NOT relax energy to the lattice). Use both as machine-precision validation invariants.
- **e-e/e-h rate scale:** 1e13–1e14 s⁻¹ at n=1e17–1e19 cm⁻³ (thermalization ~10–200 fs). Static screening under-estimates; dynamic (LOPC) needed for sub-100-fs.
- **Maxwell-SBE multiscale:** the multiscale driver runs ONE independent SBE cell per macropoint, driven by that point's macroscopic Maxwell A(t). It calls the SAME `init_sbe_bloch_solver`/`dt_evolve_bloch_cf4` as the single-cell solver, so every A–G channel flows through the `&sbe` namelist automatically — no per-channel wiring in `multiscale_ssbe`. The nonlocal collectives (Coulomb all-gather/ring, BGR + nl-II global reductions) reduce over the **per-macropoint** group `icomm_macro`, so momentum exchange is correctly confined to each cell's own k-grid. `sbe%icomm` MUST be set for every run (not only Coulomb) — it is now set unconditionally in init.
- **Provenance rule (STRICT — agreed with maintainer):** a channel may be enabled for a material ONLY if its constants are backed by a cited source for THAT material. **No source ⇒ invalid ⇒ forbidden**, and the SBE init `error stop`s — constants are NEVER transferred from another material. The registry carries per-channel gates `ii_ok/eph_ok/eeh_ok/coulomb_ok`. (This corrected an earlier mistake where CdS had estimated/copied dissipation constants — those were removed.)
- **CdS (wurtzite P6₃mc):** EPM band structure VALIDATED; dissipation channels CITED & enabled (see effect-support matrix in wiki/02). Cell = orthorhombic `al(1:3)=(a,a√3,c)` (not a single cubic constant). Form factors are the REAL cited **wurtzite** values from **Bergstresser & Cohen, Phys. Rev. 164, 1069 (1967), Table II** (LOCAL potential — spherical atomic potentials, no nonlocal term ⇒ `rvnl_tm=0`; no cited CdS nonlocal parameter exists, so none is fabricated). **Validated:** direct gap at Γ converges to **2.55 eV vs BC1967's 2.58 eV** (|Δ|≈0.03 eV); structure factors match Table II (002/101/102). Two bugs that gave 13 eV were fixed: (1) potential normalized by **total atoms/cell** (1/n, the BC1967 "volume per atom" normalization — not per species); (2) use the **wurtzite** Table II form factors, not the zinc-blende anchors. **Folding to the SBE cell is EXACT:** the orthorhombic supercell (8 atoms, the al-vector cell) block-diagonalizes over the 2 cosets (off-coset |H|≈8e-17, analogue of the cubic 4-fold FCC folding) and reproduces the same 2.54 eV gap; coset 0 = Γ_hex (direct gap), coset 1 = zone-edge (6.2 eV). **CdS dissipation — three CITED & enabled, carrier-carrier FORBIDDEN:** enabled = Fröhlich polar-optical e-ph (primary; ħω_LO=38 meV [Raman], ν_sat=α·ω_LO=2.89e13 s⁻¹ from α=0.5 [cyclotron]), Coulomb (**ε₀=8.9** [Berlincourt 1963]), impact ionization (E_th=1.5·E_g=3.6 eV; **prefactor is a fit parameter** → the run aborts unless `sbe_ii_prefactor` is set). **Carrier-carrier e-e is FORBIDDEN (`eeh_ok=.false.`)** — the strict-provenance resolution of the earlier inconsistency: there is **no cited CdS e-e rate**, so the FD-thermalization channel would borrow the generic 1e14 s⁻¹ scale cited only for GaAs/Si (Goodnick-Lugli; Fischetti-Laux), which the provenance rule forbids. The CdS literature only fixes a *timescale* (sub-100fs @ n≥1e18 [Shah 1986; Elsaesser 1991]) and an *Auger coefficient* (C=2.0e-30 cm⁶/s [Haury 1998]) — both belong to the **density-gated Auger Lindblad channel (wiki Sec 13), which is NOT yet implemented** (the constants `CDS_AUGER_C`/`CDS_EE_ACT_N` are declared but unused). A user with their own rate may opt in via `sbe_eeh_nu_sat` (same escape hatch as the II prefactor). Piezoelectric/deformation-acoustic cited but not yet SBE channels. `epm_wurtzite_cds.py`, `tests/test_wurtzite_cds_epm.py`.
- **Material registry (single source of per-material constants):** `get_material_params(name)` in `sbe_superres_ssbe.f90` returns one `s_material_params` struct (dielectric ε0/ε∞, II fit form/exponent/prefactor/threshold, e-ph phonon table + ν_sat, lattice, diamond flag) assembled from the cited module constants. Every channel auto-selects through it; the namelist defaults are sentinels (`sbe_ii_form='auto'`, `sbe_ii_exponent/prefactor/threshold ≤ 0`, `sbe_coulomb_epsilon ≤ 0`) that resolve to the material value, and any explicit namelist value overrides. GaAs resolves to the exact legacy numbers (byte-for-byte). Unknown material + a material-dependent channel ⇒ `error stop` with the supported list. **Adding a material = one `case` block + its name in `MAT_SUPPORTED`** — no edits scattered across channels. (EPM form factors are the other per-material table, already case-based in `epm_cohen_bergstresser.f90`.)
- **Fortran MPI EPM == Python reference (cubic), VERIFIED by build+run:** the Python EPM is the source of truth; the in-SALMON `src/epm` solver was using the **FCC primitive cell** (4 valence bands, Cartesian k.data, an a.u.² cutoff) and so did NOT reproduce the Python's **simple-cubic 8-atom supercell + FCC-in-cubic parity band-folding** (16 valence bands, reduced k.data, `(2π/a)²`-unit cutoff). Three concrete mismatches at the same `epm_pw_cutoff_ry`: 181 vs 171 PW, gap 1.43 vs 1.38 eV at Γ, FCC-grid vs cubic-grid k-points. **Fixed** by switching the Fortran to the simple-cubic supercell (`cb_lattice_vectors_sc`), the integer-shell cutoff, the explicit parity selection rule in `build_hamiltonian`, and reduced-coordinate `k.data`. **Verified**: built SALMON (gfortran+OpenMPI, `-fallow-argument-mismatch`) and diffed `theory='epm'` output vs `epm_gaas_reference.py` for **GaAs and Si (scalar)** — k-points exact, band energies to **5e-11 Ha**, occupations identical, valence/optical momentum to **~1e-10**; the only residual is the basis-arbitrary degenerate-subspace coupling of the single top (truncation-boundary) conduction band, which differs even between LAPACK/scipy and is physically irrelevant. Contract locked by `tests/test_epm_cubic_folding_contract.py`. **Spinor and non-cubic (CdS/graphene) Fortran EPM still TODO** — use the Python refs there.
- **Auger recombination (Sec 13) wired for CdS — standard CPTP Lindblad, doesn't break the code:** new `apply_auger_recombination` (bloch_solver) — gap-edge mean-field closure: a CB electron (ic1) recombines with a VB hole (iv1, destroying an e–h pair) and the released E_g promotes a second ic1 electron to the energy-matched hot state ic_hot. Two `amp_damp_channel` GKLS maps at the per-carrier rate γ = C·n² (R = C·n³), with occ_max-normalized [0,1]-clamped Pauli factors; trace-preserving ⇒ total carrier NUMBER conserved. Density-gated (inert below n_gate). Registry gate `auger_ok` + `auger_c_cm6s`/`auger_n_gate_cm3`; namelist `yn_sbe_auger` / `sbe_auger_c_cm6s` / `sbe_auger_n_gate_cm3`. Enabled for **CdS** (C=2.0e-30 [Haury 1998], n_gate=1e18 [Shah 1986]); **GaAs/Si/graphene forbidden** (no cited C → error stop). **Key physics (maintainer):** Auger acts on the REAL (Houston/adiabatic) populations, not the virtual driving polarization, so with the tiny cited C it is a RARE event — at normal fields the dynamics are byte-identical with the channel on vs off (verified CdS 4³: trace=32.000, Jz identical, no NaN). Its job is to be present and exactly CPTP, not to dominate. A gotcha fixed in passing: `sbe%nv_act` (gap-edge valence count) was computed only inside the impact-ionization init block; hoisted to run unconditionally so an Auger-only run has it set. `tests/test_auger_cptp.f90` (number-conserving, recombining, PSD/CPTP, Hermitian, γ=0 identity, Pauli∈[0,occ]).
- **BUG FIXED — e-ph / carrier-carrier silently skipped unless decoh or impact also on (any mode, scalar or spinor):** the Strang dissipative half-steps called `houston_dissipate` only under `if (flag_decoh .or. flag_impact)` (bloch_solver lines 910/939), but `houston_dissipate` is exactly where e-ph (`apply_eph_relaxation`) and carrier-carrier (`apply_carrier_carrier`) are applied (its INTERNAL gate `flag_impact .or. flag_eph .or. flag_eeh` is correct). So a run with **only** `yn_sbe_eph='y'` or **only** `yn_sbe_eeh='y'` never entered the block → the channel was set up (banner printed) but **never applied** (eph-only/eeh-only output was byte-identical to clean). The earlier "validated" smokes always paired e-ph with `yn_sbe_superres`/decoh/impact, masking it; the maintainer's **scalar-mode dissipation check** surfaced it. **Fix:** both gates now read `flag_decoh .or. flag_impact .or. flag_eph .or. flag_eeh`. Clean/decoh/impact/coulomb configs are byte-unchanged; only the buggy eph-only/eeh-only cases change (now active). Verified scalar GaAs 4³: eph-only nex 6.744e20 vs clean 6.757e20; eeh-only 2.398e19 vs clean 2.340e19 (both differ now, electron number conserved). **Occupation-max Pauli check (also requested):** all population-changing channels normalize their Pauli factors by `sbe%occ_max` (`merge(1,2, spinor)` — 2 scalar / 1 spinor): e-ph `1-ρ/f` (1712/1718), impact ionization `ρ/f` & `1-ρ/f` (1549-1556, 1623-1630), carrier-carrier `ρ/occ` (sbe_superres 701) with target `α·occ·f_FD`; off-diagonal coherences are damped occupation-independently (correct — coherence damping is not Pauli-blocked). So the 1-vs-2 maximum is handled correctly in both modes.
- **CdS + graphene Python EPM now emit folded scalar GS files (EPM→SBE closed, VERIFIED end-to-end):** the Python refs (`epm_wurtzite_cds.py`, `epm_graphene.py`) gained a `main_gs()` that builds the MP grid, diagonalizes the FOLDED supercell Hamiltonian, and writes `SYSNAME_k/_eigen/_tm.data` via a shared writer `epm_io.py` (byte-compatible with gs_info_ssbe; reduced-k, Hartree, rvnl_tm=0 local). **CdS:** orthorhombic 8-atom cell (al=a,a√3,c), nelec=32/nstate=32, 2-fold folded; ran in the SBE binary (clean 20-step run, trace=32 conserved, field-driven current). **graphene:** NEW rectangular 4-atom cell + 2-fold fold (PART G2), al=(a,√3a,vacuum), nelec=4/nstate=8; ran in the SBE (trace=4 conserved, in-plane current). Two graphene gotchas fixed: (1) **structure factor must be normalized by the number of primitive cells** (`struct_norm`) — the supercell sum over 4 atoms was 2× too strong vs the validated 2-atom primitive, so the folded bands didn't match; with `struct_norm=2` the rect folded-K reproduces the primitive Dirac energy 1.034 eV exactly and the cone stays gapless. (2) **It is a minimal π-model** (Ramanujam 3 form factors → Dirac cone = lowest band pair, 1 π e⁻/atom), so nelec=4 (not 16) puts the Fermi level at the Dirac point. **Fortran `src/epm` is cubic-only — NOT used for CdS/graphene.** Unfold map (`_unfold.data`) still pending (needs the SBE 4→N coset reader). `tests/test_epm_folded_gs.py`.
- **Unfold map generalized 4 cosets → N (2-coset CdS/graphene), backward-compatible:** the unfolded-population pipeline was hardcoded to the **4 FCC cosets** (`gs%unfold_w(1:4)`, `do s=1,4`, `pop_lev(1:4,1:4)`). Generalized by adding `gs%n_coset` (read from a 4th header field in `_unfold.data`; legacy GaAs files with a 3-field header fall back to 4) and changing the three `1,4` loops (`read_unfold_data` offsets, the `bloch_solver` weight loop, the `datafile`/`realtime` block writer) to `1,n_coset`. The arrays stay dimensioned (1:4) and slots > n_coset are zero, so the **GaAs 4-coset path is byte-unchanged** (verified: legacy GaAs map → trace=32, 4 cosets/k). Python side: `epm_io.compute_unfold_map` (per-k coset spectral weights + dominant coset + primitive rank from the coset-block spectrum) + `write_unfold_file`; CdS uses `orth_coset` (2-fold), graphene `rect_coset` (2-fold). Graphene end-to-end: the unfolded CB population correctly localizes at the Dirac points (zero near Γ), confirming the map places carriers in the right primitive-BZ sectors. Also added clean primitive **band-path** emission (`_bandpath.data`) for CdS (A-Γ-M-K-Γ, direct gap 2.55 eV at Γ) and graphene (Γ-M-K-Γ, gapless Dirac at K — fixed K=(2/3,1/3) reduced, not the non-corner (1/3,1/3)). `tests/test_epm_folded_gs.py`.
- **VG band budget ≠ PW cutoff:** the band count N_b carried into the dynamics is a correctness axis SEPARATE from the EPM plane-wave cutoff, and the Houston/adiabatic basis does NOT fix an insufficient N_b — it diagonalizes the already-truncated H_VG^(N_b), inheriting the truncation error (Hylleraas-Undheim-MacDonald). Monitor `P_top = max_k ρ̃_{N_b,N_b}` (criterion a): the real-time solver warns on the **error channel** and **continues** when P_top > 1e-3 (it does not stop — the user re-runs with more bands + an N_b convergence study). Re-verify N_b for every new material / driver wavelength; the converged value does not transfer. See wiki/06.
- **a.u. conversion audit (verified):** rates s⁻¹→a.u.⁻¹ via `×(au_fs·1e-15)`; energies eV→Ha via `/au_ev`; `mev_to_ha=E·1e-3/au_ev`; deformation potentials `d_evcm_to_au=D·BOHR_CM/au_ev`, `d_evang_to_au=D·BOHR_ANG/au_ev`; mass density `rho_gcm3_to_au=ρ·BOHR_CM³/ME_G`; carrier density a.u.⁻³→cm⁻³ via `1e24/BOHR_ANG³`. The e-ph `eph_wrel=D²/ħω` is normalized to a dimensionless per-mode weight, so GaAs (eV/Å) vs Si (1e8 eV/cm) deformation-potential units only need within-material consistency (the absolute scale is `ν_sat`).

## CPTP invariants (must hold for every new dissipator)
- Each exp(τD) a genuine GKLS map; clamp all (1−ρ) Pauli factors to [0,1].
- Yoshida wraps ONLY the unitary part; dissipator always τ>0.
- Hermitian-only dissipator must conserve populations exactly (Γ_aa=0).

---

## Test inventory (`tests/`)
| Test | Covers | How |
|---|---|---|
| `test_si_epm_gap.py` | Part A | primitive sublattice-block Si gap ≈1.06 eV, CBM≈0.86·2π/a |
| `test_ii_form_switch.py` | Part B | rate γ=P(ε−E_th)^a scaling for a=2 vs 4 |
| `test_hf_sublattice_proj.py` | Part E | proj zeroes inter-sublattice, keeps diag, Hermitian; real GaAs weights |
| `test_superres_rates.f90` | Part C primitives | ν(ε) limits, N_B, energy-bin area/peak, Fröhlich asinh, II 2^a scaling, BGR −19/−41 meV, data tables (standalone gfortran) |
| `test_eph_cptp.f90` | Part C5/C6 | amplitude-damping map: trace, qubit positivity (det≥0), transfer formulas, coherence damping, Hermiticity, γ=0 identity |
| `test_superres_rates.f90` (extended) | Part C5 | + unit conversions (meV/eV·cm/eV·Å/g·cm⁻³→a.u.), golden_rule_prefactor, eph_thermal_split (fe+fa=1, fe/fa=(N+1)/N) |
| `test_screening.f90` | Part G | eps_TF limits, Lindhard F(0)=1/F(1)=½/monotone, eps_Lindhard→TF at small q, plasmon ω_p², LOPC Vieta sum/product + anticrossing |
| `test_carrier_carrier.f90` | Part F | carrier_carrier_relax conserves number+energy, CPTP positivity, EID damping, converges to FD, inversion no-op; fit_fermi_dirac self-consistency |
| `test_vg_basis_nb.f90` | VG basis sufficiency | eta admixture + threshold, conv-error metric, P_top gate; Hylleraas-Undheim-MacDonald interlacing/upper-bound + 2nd-order truncation-shift formula (self-contained Jacobi solver, no LAPACK) |
| `test_material_registry.f90` | material registry + provenance gates | GaAs/Si all gates `.true.`; CdS gates (e-ph/Coulomb/II `.true.`, **e-e `.false.` forbidden**, II prefactor sentinel); unknown→not-found |
| `test_wurtzite_cds_epm.py` | CdS EPM geometry + cited FFs | orthorhombic cell from al, 8 atoms (4 Cd+4 S), Hermitian H, broken inversion, exact BC1967 Table II anchors; band solve reported NOT-yet-validated (not asserted) |
| `test_epm_cubic_folding_contract.py` | Fortran↔Python EPM convention (cubic) | simple-cubic basis count (171 PW @ cutoff 11.1), exact 4-coset (FCC-in-cubic) folding block-diagonality, folded spectrum = union of coset spectra, reduced-k MP grid (first pt −3/8) — the contract the Fortran `src/epm` solver must mirror |
| `test_auger_cptp.f90` | Auger recombination map (Sec 13) | replicates apply_auger_recombination's two-amp_damp map: number (trace) conserved, excited population decreases (recombination), PSD/CPTP, Hermitian, γ=0 identity, populations in [0,occ], recombination rate→0 with no hole (standalone gfortran) |
| `test_epm_folded_gs.py` | CdS+graphene folded GS + bandpath + 2-coset unfold (Python) | epm_io writer round-trip + non-orthogonal-cell rejection; graphene rect 4-atom 2-fold exact + gapless folded Dirac + struct_norm reproduces primitive Dirac + bandpath-K is the gapless Dirac point; CdS bandpath Γ gap ~2.5 eV; main_gs band count/occupations/reduced-k; graphene 2-coset unfold map (n_coset=2, weights sum to 1, exact one-coset folding) |
| _(add per increment)_ | | |

End-to-end smoke (manual): scalar GaAs 4³, `yn_sbe_superres='y'` + `yn_sbe_eph='y'`
→ runs stable, correct diagnostics (ħω=36 meV, N_B=0.33 @300 K), trace = 32.000
conserved over all steps (CPTP at the dynamics level).

End-to-end smoke (manual, not in run_all): scalar GaAs 4³, `yn_sbe_coulomb='y'`
+ `yn_sbe_hf_sublattice_proj='y'`, weak field → runs stable, finite, projection
diagnostic confirmed. Full "zero spurious Γ→X transfer" weak-field validation
is a longer run, deferred to a dedicated validation pass.

Run all: `python3 tests/run_all.py` (each test prints PASS/FAIL and exits nonzero on failure).

---

## Materials — Python EPM references (THE source of truth) ✅
The **Python EPM references are primary** and each is validated against a cited
benchmark (run `python3 tests/run_all.py`, currently 11/11). The fast in-SALMON
**MPI EPM (`src/epm/`) is SECONDARY** — calibrated against the Python refs and
**deprecated for the non-cubic materials (CdS, graphene) until debugged**;
generate those GS with the Python references.

| Material | Module | Band-validated | Folding | SBE GS files (eigen/tm/k + unfold) |
|---|---|---|---|---|
| GaAs | `epm_gaas_reference.py` | CB 1966 | ✅ 4-fold FCC (exact) | ✅ emitted (full EPM→SBE pipeline) |
| Si (Kunikiyo, default) | `epm_si_reference.py` | gap 1.059 eV, CBM 0.850 | ✅ 4-fold FCC (reuses GaAs) | ✅ emitted (full pipeline) |
| Si_cb (Cohen-Bergstresser) | `epm_si_reference.py --variant Si_cb` | gap 0.818 eV, CBM 0.850 | ✅ 4-fold FCC | ✅ emitted |
| CdS | `epm_wurtzite_cds.py` | gap 2.55 vs 2.58 eV (BC1967) | ✅ 2-fold orth←hex, **verified exact** (off-coset \|H\|≈8e-17) | ✅ **Python emits eigen/tm/k + bandpath + 2-coset unfold**; **ran end-to-end in the SBE** (trace=32 conserved) |
| graphene | `epm_graphene.py` | zero gap at K, v_F 9.6e5, Γ −7.78, M −2.70 | ✅ **rectangular 4-atom + 2-fold, verified exact** (off-coset=0; reproduces primitive Dirac) | ✅ **Python emits eigen/tm/k + bandpath + 2-coset unfold**; **ran end-to-end** (trace=4; unfolded carriers localize at the Dirac points) |

**HONEST status:** **all four materials now emit SBE GS files and run in the
SBE.** GaAs/Si via the Fortran `theory='epm'` path (verified == Python) and the
Python ref; **CdS and graphene via the Python references** on their folded
supercells (CdS orthorhombic 2-fold; graphene rectangular 4-atom 2-fold — both
verified exact and both ran end-to-end in the SBE binary with trace conserved).
**The Fortran `src/epm` solver is cubic-only and is NOT ready for CdS/graphene**
— use the Python refs for those. Still pending for the 2-fold materials: the
band-**unfolding** map (`_unfold.data`), which needs the SBE reader generalized
from the hardcoded 4 cosets to N (gs_info_ssbe read_unfold_data + the
bloch_solver unfold loops); the spinor Fortran path; and the Auger Lindblad
channel. **graphene is a minimal π-model** (Ramanujam 3 form factors): Dirac
cone = lowest band pair, 1 π e⁻/atom → 4-atom cell uses nelec=4, nstate=8,
occ 2/π-band, Fermi at the Dirac point (not a full 4-electron valence model).

**Si vs Si_cb:** identical machinery (diamond, V^A≡0, a=10.26 Bohr, 4-fold fold);
ONLY the V^S triplet differs (Kunikiyo vs CB). See README "Supported materials".

## Next action on resume — STANDING TODOs (context restarted)
Branch `develop-2.0.0` (merged from PR #44). Order of priority per the maintainer:

1. **Python EPM refs for all 4 materials FIRST, MPI EPM second.** Done: all four
   Python refs validated (table above). **MPI Fortran EPM (cubic GaAs/Si/Si_cb)
   now VERIFIED byte-equivalent to the Python ref** (built SALMON + diffed; see the
   decisions-log entry "Fortran MPI EPM == Python reference"). **Still TODO:** the
   Fortran EPM for the **spinor** GaAs path and the **non-cubic** materials (CdS
   orthorhombic 2-fold, graphene) — use the Python refs there until debugged.

2. **Sublattice-projection-with-Coulomb is the ONLY correct HF mode — verify on
   ALL materials.** `apply_hf_sublattice_projection` zeroes the spurious
   inter-coset Σ^HF created by folding. The cosets DIFFER per material (GaAs/Si
   = 4-fold FCC; CdS/graphene = 2-fold, different coset vectors). TODO: confirm
   `yn_sbe_hf_sublattice_proj='y'` + `yn_sbe_coulomb='y'` uses the correct
   per-material coset/unfold weights (gs%unfold_w) and is block-diagonal to
   machine precision for each. Currently wired for the FCC 4-fold path; the
   CdS/graphene 2-fold unfold maps must feed the same projector.

3. **GRAPHENE dissipators (G4/G5) + no-Kuhn-Zurek (G6)** — registry entry
   `epm_material='graphene'` NOT yet added. Needed (constants in wiki, PROMPT_graphene_full):
   - e-ph: E2g 196 meV (g²=0.0405 eV²), A1' 160 meV (g²=0.0994 ×2 GW), acoustic
     D_ac=16 eV (tunable); total ~1e10 s⁻¹ @300K, 196 meV threshold step.
   - e-e/Auger: α=2.2/ε_eff, static|dynamical RPA switch, gapless CM up to ~2.
   - **assert: graphene + Kuhn-Zurek flag → error** (coherence is many-body only).
   - is_diamond=.true. (V_A=0, centrosymmetric); odd-only HHG (linear), 6m±1 (circular).

4. **AUGER CPTP Lindblad primitive** (wiki Section 13) — ✅ **DONE for CdS**
   (`apply_auger_recombination`, `yn_sbe_auger`; see the decisions-log entry
   "Auger recombination wired for CdS"). Number-conserving, density-gated, exactly
   CPTP (two `amp_damp_channel` maps), rare by construction (acts on real Houston
   carriers). `tests/test_auger_cptp.f90`. **Still TODO:** graphene Auger/CM (G5)
   needs the graphene registry entry (item 3) + a cited gapless-CM coefficient;
   GaAs/Si Auger is gated forbidden (no cited C). The provenance resolution below
   stands (the Haury/Shah constants feed THIS channel, NOT carrier-carrier):
   **RESOLVED (provenance, conservative):**
   the cited Haury/Shah constants describe THIS Auger channel, not the existing
   FD-thermalization (`yn_sbe_eeh`) channel — which has no cited CdS *rate* and
   would borrow the GaAs/Si 1e14. So CdS `eeh_ok` is now **`.false.` (forbidden)**,
   matching README/wiki/02 and the test (`test_material_registry.f90` asserts
   forbidden). The constants `CDS_AUGER_C`/`CDS_EE_ACT_N` are declared (cited) but
   unused until this Auger channel is written; only then does a CdS density-gated
   carrier-multiplication channel switch on (it will NOT re-enable the borrowed-rate
   FD channel). A user with their own CdS e-e rate can still opt in via
   `sbe_eeh_nu_sat`.

5. **WIKI COMPLETENESS sweep (audited 2024-restart).** Fixed: wiki/00 line 58
   CdS "all dissipation FORBIDDEN" contradiction; wiki/00 roadmap line; the
   `[md]` shorthand refs in wiki/02 → real source labels; wiki/02 matrix + §12 +
   §13 are current. **✅ DONE (this doc pass):**
   - **wiki/04** (config examples) now has **CdS and graphene recipes** (Python-EPM
     workflow, al-vector cells), the Auger params (`yn_sbe_auger`,
     `sbe_auger_c_cm6s`, `sbe_auger_n_gate_cm3`), `epm_material` extended to
     CdS/graphene with the Python-reference note, and the unfold/bandpath plotting.
   - **wiki/05** (folding/unfolding) now has §6 (**CdS 2-fold orth←hex** +
     **graphene 2-fold rect←hex**) and §7 (the **N-coset unfold** generalization,
     4 cubic / 2 wurtzite-rectangular).
   - **wiki/01** (physics models) now has §12 (**Auger recombination** CPTP
     Lindblad), the e-ph/eeh gate-fix note (§6), the carrier-carrier ✅ status, and
     the **occ_max Pauli normalization** note (§8). STILL missing: the graphene
     dissipation channels (e-ph E2g/A1', gapless-CM Auger, no-Kuhn-Zurek G6) — they
     await the graphene registry entry (item 3).

6. **GRAPHENE FOLDING (G2) — ✅ DONE (rectangular 4-atom cell).** `epm_graphene.py`
   now has the orthorhombic 4-atom rectangular cell (a × √3a, zigzag x / armchair
   y), the 2-fold hex→rect fold (`rect_folding_check`: off-coset |H|=0, exact),
   and the cell reproduces the primitive Dirac cone (gapless) once the structure
   factor is normalized per primitive cell (`struct_norm`). GS files emit on this
   folded cell. **Remaining:** the `graphene_unfold.data` MAP (see item 7) — only
   the folding + GS emission are done, not the inverse unfold map.

7. **WIRE CdS + graphene EPM → SBE GS files — ✅ DONE (eigen/tm/k + bandpath +
   2-coset unfold map).** `epm_wurtzite_cds.py` and `epm_graphene.py` emit
   `SYSNAME_k/_eigen/_tm.data`, `_bandpath.data` (clean primitive bands), and
   `_unfold.data` (2-coset) via the shared `epm_io.py`; all run end-to-end in the
   SBE binary (trace conserved). **The 4→N coset generalization is DONE:**
   `gs_info_ssbe::read_unfold_data` reads `n_coset` from the map header (4 = cubic
   FCC, 2 = wurtzite/rectangular; legacy 3-field GaAs headers fall back to 4),
   `gs%n_coset` is bcast, and the `bloch_solver` population loop + the
   `datafile`/`realtime` unfold-block writer run over `1..n_coset`. Verified:
   GaAs 4-coset (legacy file) trace=32 + 4 cosets/k; graphene 2-coset trace=4 +
   2 cosets/k with carriers correctly localized at the Dirac points; the plotter
   renders both the unfolded and folded k–t maps.

8. Optional refinements (unchanged): inter-k momentum-resolved F/C4 on the ring;
   dynamic LOPC; absolute golden-rule e-ph prefactor; Chefonov Si bleaching run.

Also pending (validation, not code): a longer Si super-mode run to check the
e-ph cooling actually reduces the Drude conductivity (bleaching) -- needs a Si
dataset (epm_material='Si', scalar) and the carrier-density/current diagnostics.
