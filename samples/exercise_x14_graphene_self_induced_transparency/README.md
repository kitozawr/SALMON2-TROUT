# exercise x14 — graphene **self-induced transparency** at 1–100 kV/cm: field before / field after

**What this exercise delivers.** The readiness check the maintainer asked for
before the graphene self-induced-transparency (SIT) study: which effects are
ready, how the **transmission coefficient is obtained from the field before and
after the sheet**, the **level check** of the Dirac cone the SBE runs on, and
the **population-saturation check** (Auger recombination and impact ionization /
carrier multiplication balance each other above a density — F. Rana, PRB 76,
155431 (2007)). Two things were found and fixed on the way: the production
graphene EPM basis had a **spurious 0.21 eV gap at the Dirac point**, and
graphene lacked its own virtual-population filter — it now has the **2D colmem
analog** ([`wiki/10` §8.11](../../wiki/10_open_quantum_systems_literature.md)).

Everything below was run in this repository (24² smoke grid, serial build);
the shell-resolving production runs are a cluster job whose recipe and budget
are given in §5.

---

## 1. Physics in one paragraph

A free-standing graphene sheet absorbs the universal fraction **πα = 2.29 %**
of normally incident light whenever the interband transition at the photon
energy is Pauli-allowed (ħω ≫ 2E_F, k_BT). The resonant shell in k-space has
radius `k_res = ħω/2ħv_F` around K/K′. A strong pulse fills that shell
(Pauli blocking / Rabi saturation) and the sheet **bleaches** — its
transmission rises: self-induced transparency. What limits the bleaching is
how fast the shell refills/empties: optical-phonon emission (E2g 196 meV,
A1′ 160 meV, ~20 fs) and the **Coulomb channel**, where on the gapless cone
impact ionization (carrier multiplication, CM) and Auger recombination are
exact time-reverses. Their balance `R = G` holds when the electron and hole
quasi-Fermi levels coincide; for a symmetric e–h population that is μ = 0, i.e.
the **intrinsic thermal density**

    n₀ = n_i(T) = (π/6) (k_B T / ħ v_F)²  = 8.1×10¹⁰ cm⁻² at 300 K  (∝ T²)

— below it the plasma net-multiplies, above it net-recombines: **the pair
population saturates at n_i(T)** (unit-tested, `tests/test_rana_saturation.f90`).
In the field range 1–100 kV/cm the intensities are 1.3×10³–1.3×10⁷ W/cm² and
the vector-potential excursion `A₀ = E₀/ω` is tiny (6.6×10⁻⁴ a.u. at 100 kV/cm,
0.8 eV): this is the **perturbative-interband, weak-dressing** regime — the
Rabi area of the 8-cycle pulse reaches 0.25 rad at 100 kV/cm, so saturation
sets in at the top of the range. (The THz regime of the same fields, where
doped CVD graphene bleaches by Drude-carrier heating, needs two increments the
code does not have yet — §7.)

## 2. Readiness matrix (graphene, in-plane field, 1–100 kV/cm)

| effect | status | note |
|---|---|---|
| EPM ground state (π/π* Dirac pair) | ✅ **fixed** | 43-PW basis (`epm_pw_cutoff_ry = 29.4`); the old 7-PW basis is gapped at K by 0.21 eV (§4) |
| coherent velocity-gauge SBE (CF4 + Houston) | ✅ | 2 bands, the VG "basis-edge" monitor fires for ANY excitation (top band = the conduction band) — expected, ignore |
| e-ph optical E2g/A1′ (inter-k ring) | ✅ | Piscanec/Lazzeri couplings, ring mandatory |
| e-ph acoustic (D = 16 eV, TF-screened, grid-resolved q) | ✅ | Hwang–Das Sarma |
| 2D Rana Auger / carrier multiplication (ring) | ✅ | R−G on the gathered n, p; **evaluated at the e-ph bath T** (no carrier T_e yet — §7) |
| 2D-sheet Σ^HF (`yn_sbe_coulomb`) | ✅ (off by default here) | renormalizes v_F → shifts k_res; enable deliberately |
| Option A dressed reference (`yn_sbe_dressed_ref`) | ✅ | material-agnostic; on the exactly gapless K point the frozen state maps to equilibrium exactly |
| **virtual-population filter for graphene — the 2D colmem analog** | ✅ **NEW** | `yn_sbe_colmem` + `yn_sbe_colmem_pop` now run on graphene: phonon lines for the e-ph sectors, the **2D Dirac-plasmon line** ω_pl(n,p; Q_TF) for the Rana source (§3.3, `wiki/10` §8.11) |
| impact ionization (`yn_sbe_impact_ionization`) | 🚫 physics | gap-threshold law meaningless on a gapless cone — CM *is* the Rana channel |
| carrier-carrier FD fit (`yn_sbe_eeh`) | 🚫 physics | no cited graphene rate |
| Kuhn–Zurek dephasing (`sbe_decoh_*`) | 🚫 physics | many-body coherence loss (error stop) |
| **transmission from field before/after** | ✅ post-processing | sheet boundary condition on the `Jm` column (§3); the single-cell driver itself has no self-consistent field — the `S_rr` diagnostic tells when that matters |
| doped / finite-T initial occupation | — limit | not available (`gs%occup` = integer filling) — §7 |
| carrier temperature in the Rana rates | — limit | uses the lattice T — §7 |

## 3. Field before → field after → transmission

### 3.1 The sheet boundary condition (`transmission.py`)
The single-cell driver writes `E_ext` (the incident field it was driven with)
and the matter current `Jm` (electron current per cell volume, vacuum included)
to `*_sbe_rt.data`; its `E_tot` column equals `E_ext`. For a 2D sheet at normal
incidence, E is continuous across the sheet and H jumps by the sheet current, so
(Hartree a.u., `Z₀ = 4π/c`)

    J_s = −Jm · L_z                       physical charge current per unit width
    E_t = (2 E_inc − Z₀ J_s)/(1 + n_sub)   field AFTER   (n_sub = 1 free-standing)
    E_r = E_t − E_inc                      reflected

(`J_phys = −Jm` because the driver's energy ledger is `dW = −E·Jm V dt`.) Power,
**fluence-integrated** (Parseval-exact — the primary numbers):

    T = n_sub ∫E_t² / ∫E_inc²,   R = ∫E_r² / ∫E_inc²,   A = 1 − T − R
    (c/4π)(E_inc² − E_t² − E_r²) = E_t J_s   pointwise  ⇒  A = Joule absorption in the transmitted field

Linear universal sheet (σ = e²/4ħ = ¼ a.u.): `T = 1/(1+π/2c)² = 0.97746`,
`A = πα/(1+πα/2)² = 0.02241`, `R = 1.3×10⁻⁴` — reproduced by
`tests/test_sheet_transmission.py`, together with the Fresnel substrate limit
`T = 4n/(1+n)²` and the pointwise energy identity for an arbitrary current.

### 3.2 The self-consistency (radiation-reaction) diagnostic
The SBE was driven by `E_inc`, not by the local field `E_t`. Its own energy
ledger gives `A_E = ∫E_inc J_s / F_inc`, and **exactly** `A = A_E − S_rr` with

    S_rr = (Z₀/2) ∫J_s² / F_inc .

On a mesh that resolves the shell `S_rr ≈ (Z₀σ)²/2 ≈ 3×10⁻⁴` (negligible against
2.3 %). On an *unresolved* mesh a few near-resonant k-points carry a huge
**reactive** current and `S_rr` becomes comparable to `A` — the script prints it
as the reliability flag. When it is not ≪ A, use an exact route: (a) a
radiation-reaction term in the single-cell driver, `dA_ind/dt = (2π/c)·J_s`
(one line in `realtime_ssbe`, not yet added — then `E_tot` *is* the transmitted
field), or (b) `theory='maxwell_sbe'` (1D FDTD, incident/reflected/transmitted
written to `*_sbe_wave.data`) with **`hx_m = L_z = 37.79 a.u.`** so the FDTD
cell carries the same sheet current as the SBE cell, `dt ≤ hx_m/c = 0.28 a.u.`
(CFL) and enough `nxvac_m` to hold the pulse — a recipe, not shipped/tested here.

A single-FFT-bin "T at the carrier" is deliberately **not** reported: pulse
reshaping moves spectral weight between bins and such a ratio is not bounded
by 1 (it came out 1.08 on the smoke mesh). The script gives the band-integrated
`T_band` over the incident FWHM band and `Re σ/σ_univ` (= 1 for the universal
sheet) as secondary diagnostics.

### 3.3 The graphene 2D colmem analog (what "ready" required)
The maintainer's directive: graphene needs its own virtual-population filter
("2d colmem аналог"). The 2026-07-20 graphene exclusion of `yn_sbe_colmem/_pop`
was a guard, not physics; it is lifted, and the gapless cone gets the one
ingredient the 3D machinery lacked. Graphene's population channel is the Rana
**Coulomb** rate model on the *global* densities n, p — so the memory line of
that sector is not a phonon energy but the plasma response, the build-up time
of screening. Implemented (`wiki/10` §8.11; no new inputs, no new parameters):

| sector | memory line |
|---|---|
| e-ph coherence / e-ph source | graphene phonon table (E2g, A1′, acoustic) — the standard `colmem_lines` |
| **Rana (Coulomb) source n, p** | **2D Dirac plasmon** `ω_pl² = 2(W_c+W_v)Q_TF/ε_r`, `W(μ) = 2k_BT ln[2cosh(μ/2k_BT)]` (Falkovsky–Varlamov Drude weight), at the collision's own screening momentum `Q_TF` [R07 Eq. 13] — Hwang & Das Sarma, PRB 75, 205418 (2007) |

At 300 K, ε_r = 10: ω_pl = 31 meV (intrinsic), 133 meV at 10¹² cm⁻² — the phonon
scale, all ≪ 2ħω. The filter passes a constant density exactly (calibrated R07
rates untouched) and transmits the 2ω breathing at |R(2ω)| = 0.16
(`tests/test_colmem_2d.f90`). Together with `yn_sbe_dressed_ref` (which removes
the Dirac-point rotation background at the source) this is the graphene version
of the three-sector fix of `wiki/10` §8.10 — the `mem` variant of this exercise.

## 4. Level check — the Dirac cone the SBE actually runs on

`python3 tests/test_graphene_dirac_levels.py` (numpy only):

| basis (`epm_pw_cutoff_ry`) | PW | gap at K | v_F | Γ bottom / M dip |
|---|---|---|---|---|
| 2.94 (the old x11 input) | 7 | **0.2125 eV — spurious** | — | −8.5 / −2.8 eV |
| **29.4 (x11 + x14 now)** | 43 | 6.5×10⁻⁶ eV (Python), **0.0000 eV (Fortran bandpath)** | 0.960×10⁶ m/s | −7.78 / −2.70 eV (thesis windows) |

The 7-vector set {0, first shell} is not closed under the little group C₃ᵥ of K
(a rotation about K maps a first-shell G onto b₂−b₁, a second-shell vector), so
the truncation breaks the symmetry protection of the Dirac degeneracy. 43 PW is
shell-complete to n = 12 and restores it; the SBE still uses `nstate = 2`, so the
larger GS basis costs nothing. Also: the SALMON Monkhorst–Pack mesh is
half-shifted, `(2i−N−1)/2N` — **K = (2/3, 1/3) is on the mesh only for ODD
multiples of 3** (9, 15, …, 147, 153); 12, 24 and 150 straddle K by half a step.
The EPM v_F (0.96×10⁶ m/s) is 4 % below the 10⁸ cm/s the Rana constants assume
(n_i ∝ 1/v_F²: 8 %).

## 5. k-mesh resolution and the production budget

The resonant shell must be sampled: mesh points per shell radius
`= (ħω/2v_F)/(|b|/N)`:

| N | 12 | 24 | 48 | 96 | **150** | 300 |
|---|---|---|---|---|---|---|
| 0.4 eV | 0.13 | 0.26 | 0.51 | 1.0 | 1.6 | 3.2 |
| **0.8 eV** | 0.26 | 0.53 | 1.0 | 2.1 | **3.2** | 6.4 |
| 1.5 eV | 0.48 | 0.97 | 1.9 | 3.9 | 6.0 | 12 |

Below ~1 the mesh sees a few discrete near-resonant k-points instead of a shell
(the 24² smoke mesh has exactly 4 points at 0.777 eV, 40 meV under the 0.816 eV
carrier — hence Re σ/σ_univ = 4 there). **Production at 0.8 eV: N ≥ 150** (odd
multiple of 3 if the exactly gapless K point itself is wanted: 147 or 153).

Cost: the graphene ring is e-ph + Rana only — **O(nk²)**, exp-bound (no O(nk³)
II/Auger kernel). Measured: ≈ 40 ms/step at 24² (576 k) on 4 threads ⇒ ≈ 1.3
s/step at 150² on 48 threads ⇒ **≈ 25 min per 100-fs run** (1048 steps at
dt = 4 a.u.); the coherent runs are seconds. The 15-run scan ≈ 6 h; the
essential set (100 kV/cm × {coh, diss, mem} + dark) ≈ 1.5 h.

## 6. Run

```bash
cd samples/exercise_x14_graphene_self_induced_transparency
python3 make_inputs.py --nk 24 --outdir smoke_nk24      # pipeline smoke (seconds per run)
python3 make_inputs.py --nk 150 --outdir prod_nk150     # production (cluster; see wiki/11)
cd smoke_nk24 && cp ../run_scan.sh . && OMP_NUM_THREADS=4 SALMON=../../../build/salmon bash run_scan.sh
python3 ../transmission.py runs/E*/graphene_sit_sbe_rt.data --plot        # -> transmission_scan.png
python3 ../saturation_check.py runs/E100kVcm_{coh,diss,mem}/graphene_sit runs/dark_mem/graphene_sit --plot
```
`make_inputs.py` is the single source of truth (GS + RT share nk/cell/sysname;
`--hw-ev`, `--fields`, `--cycles`, `--sigma-ev`, `--temp-k`, `--eps-r`); the
variants are `coh` (no dissipation), `diss` (ring: e-ph + acoustic + Rana,
Markovian — as x11) and `mem` (diss + the 2D colmem analog + dressed reference).
`rt_dark_*` are the zero-field controls.

## 7. What the 24² smoke scan shows (pipeline validation — the mesh is UNRESOLVED, magnitudes are not physics)

17 runs, **electrons = 2.000 at every step in every run**, all banners present
(`# graphene: 2D collisional-memory analog`, `# 2D colmem analog: R − G evaluated
on n, p memory-filtered with the 2D Dirac-plasmon line …`, field audit
`peak |E| = 0.100 MV/cm`).

| E₀ [kV/cm] | variant | T | R | A (sheet BC) | A_E (ledger) | S_rr | Re σ/σ_univ |
|---|---|---|---|---|---|---|---|
| 1 … 100 | coh | 0.95600 | 0.0132 | 0.03078 | 0.05722 | 0.0264 | 4.11 |
| 1 … 100 | diss | 0.95446 | 0.0133 | 0.03227 | 0.05881 | 0.0265 | 4.18 |
| 1 … 100 | mem | 0.95592 | 0.0132 | 0.03085 | 0.05731 | 0.0265 | 4.12 |

- **Linear regime:** T, R, A are field-independent from 1 to 100 kV/cm to 10⁻⁴
  (the four detuned discrete points do not saturate like a shell; the coherent
  100 kV/cm run bleaches by only 7×10⁻⁴ of A). Saturation needs the resolved
  shell (§5).
- **Bookkeeping closes:** 100 kV/cm makes 9.5×10¹⁰ pairs/cm² (`saturation_check`),
  which is exactly the ledger's 5.7 % of the 2.05×10⁻⁷ J/cm² fluence at 0.816 eV
  (πα would give 3.6×10¹⁰ — the 4-point mesh over-absorbs ×2.5); `S_rr ≈ A`
  flags the unresolved mesh as intended.
- **The mechanism of `wiki/10` in miniature:** `diss` adds +5 % absorption —
  Markovian e-ph dephasing converts the reversible resonant polarization into
  real population (B25's "dephasing ionization"); `mem` returns to the coherent
  value (+0.2 %): the 2D colmem analog removes the fabricated part. At these
  fields the *dressing* itself is negligible (virtual fraction at peak 0.000 —
  A₀ ≪ Δk), so the population filters change little; their job is the strong-field
  / THz end.
- **Saturation direction is right:** 100 kV/cm gives 9.5×10¹⁰ cm⁻² > n_i(300 K)
  = 8.08×10¹⁰ → the Rana ledger ends in **net recombination** (−1.5×10⁸ pairs/cm²
  over the tail); 10 kV/cm gives 1.0×10⁹ < n_i → **net carrier multiplication**
  (+7.6×10⁷, +7 %). The dark runs stay exactly at zero: the vacuum is a fixed
  point of CVCC generation (a hot third carrier is needed), so there is no
  spurious dark drift.

Figures from the smoke scan: `transmission_scan.png`, `saturation_check.png`.

## 8. What to look for at production (150², 0.8 eV)

- `Re σ/σ_univ → 1` and `S_rr ≪ A` in the `coh` 1 kV/cm run — the linear
  universal-sheet check (T → 0.9775). If not, the shell is still unresolved.
- **T(E₀) rising** toward 100 kV/cm: the coherent bleaching fraction scales with
  the Rabi area² (≈ 0.06 at 100 kV/cm ⇒ ΔA/A ≈ −3 %); `diss` recovers part of
  it (phonon emission empties the shell in ~20 fs); `mem` is the physically
  correct dissipative answer.
- The post-pulse density: above `n_i(T)` it decays on the R07 lifetime (~ps at
  10¹² cm⁻²; only the first ~60 fs of the tail are in the run — the slope, not
  the plateau, is what you see). Remember the plateau the code relaxes to is
  `n_i(T_bath)`; a hot plasma (T_e ~ 1000–3000 K ⇒ n_i ~ 10¹²–10¹³) would sit
  higher (§9).

## 9. Limits recorded (not blockers for the near-IR study)

1. **No doping / finite-T initial state.** `gs%occup` is integer filling. The
   THz-regime SIT of doped CVD graphene (E_F ≈ 0.1–0.4 eV, Drude heating,
   Hwang 2013 / Mics 2015 / Hafez 2018) needs (a) an FD(E_F, T) initial state
   (touches the dressed-reference formula, which assumes full/empty bands) and
   (b) a k-mesh refined around K by ~10² beyond what a uniform MP grid affords.
2. **T_e in the Rana channel.** The balance density is evaluated at the bath T;
   a carrier temperature from the two moments (n, energy) of the gathered Dirac
   populations would move the plateau to `n_i(T_e)`. Bounded increment.
3. **Self-consistent sheet field.** The radiation-reaction term in the single-cell
   driver (§3.2) — one line plus a flag; until then `S_rr` is the flag.
4. The VG basis-edge monitor warns on the 2-band model whenever the conduction
   band is populated; the warning is expected here.

## 10. Tests added

`tests/test_graphene_dirac_levels.py` (43-PW cone: gapless, v_F, e–h symmetry,
linearity, thesis windows; the 7-PW trap; the resolution advisory),
`tests/test_sheet_transmission.py` (Fresnel, universal sheet, energy identity,
spectral/fluence), `tests/test_rana_saturation.f90` (n₀ = n_i(T), T² law,
two-sided monotone CPTP saturation), `tests/test_colmem_2d.f90` (plasmon-line
limits, filter fixed point and |R(2ω)|, Rana source overrides). Suite:
`python3 tests/run_all.py`.
