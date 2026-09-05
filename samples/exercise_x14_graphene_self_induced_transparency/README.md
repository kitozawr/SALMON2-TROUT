# exercise x14 — graphene under the DAST THz transient (1–100 kV/cm, normal incidence): field before / field after, Zener pair creation, two-temperature Coulomb sector

**What this exercise delivers.** The cluster-ready study the maintainer asked for
before launching the graphene self-induced-transparency (SIT) runs: the **dense
production grid**, the **self-consistent sheet field** (so the solver's total field
*is* the transmitted field), the **DAST single-cycle THz drive** rescaled to
1–100 kV/cm, the **two-temperature model** of the Coulomb (Auger / carrier-
multiplication) sector with cooling through the phonon channel, the **level
check** that exposed a spurious 0.21 eV Dirac-point gap in the old basis, and
the **population-saturation check** (Auger and impact ionization balance at
n_i(T), Rana 2007). The methods write-up is
[`wiki/12_graphene_sheet_solver.md`](../../wiki/12_graphene_sheet_solver.md);
the collisional-memory design is `wiki/10` §8.11.

Geometry: THz beam at **normal incidence** on a free-standing monolayer, E in the
sheet plane (x); everything in `A_eV_fs` units.

---

## 1. Physics and the honest expectation (read before running)

**The drive.** `DAST_singlecycle_100kV.txt` (exercise x12): a Gaussian A(t), so
E = −dA/dt is one clean cycle at **3.36 THz** (ħω = 13.9 meV, period 298 fs),
285 fs long. `make_inputs.py` rescales it to the exact peak field (1, 3, 10, 30,
100 kV/cm), removes the initial offset and cos²-windows the ends (the SALMON
file reader returns 0 outside the file — an offset would be a δ-spike in E).
Peak k-excursion `A₀ = E₀/ω`: **0.062 a.u. at 100 kV/cm** = 6 spacings of the
147² mesh (0.0106 a.u.); 6×10⁻⁴ a.u. at 1 kV/cm.

**What the field does to the cone.** At THz the photon is far below any
interband transition the mesh can resolve; the physics is the field *sweeping*
k + A(t) through the Dirac point: massless **Landau–Zener / Schwinger pair
creation** with P(k⊥) = exp(−π v_F k⊥²/E) [Allor–Cohen–McGady 2008; Dóra–Moessner
2010], rate per area Γ = (g/4π²) E^{3/2}/v_F^{1/2} (g = 4) — see `wiki/12` §7:

| E₀ [kV/cm] | A₀ [a.u.] | A₀ / Δk (147²) | LZ tube k⊥ = √(E/πv_F) [a.u.] | pairs per passage (analytic) |
|---|---|---|---|---|
| 1 | 6.2×10⁻⁴ | 0.06 | 3.7×10⁻⁴ | ~1.5×10⁹ cm⁻² |
| 10 | 6.2×10⁻³ | 0.6 | 1.2×10⁻³ | ~5×10¹⁰ |
| 30 | 1.9×10⁻² | 1.8 | 2.0×10⁻³ | ~2.5×10¹¹ |
| 100 | 6.2×10⁻² | 5.9 | 3.7×10⁻³ | ~1.5×10¹² |

(two passages per cycle, Stückelberg interference neglected). The mesh
*resolves* this process when A₀ spans ≥ 2 spacings — **E₀ ≳ 30 kV/cm at 147²**;
below that the created density is the analytic one and the mesh only samples
the K point itself (147 = odd×3 ⇒ K exactly on the half-shifted MP mesh).

**Expected transmission range.** The sheet starts *empty* (T = 0 filling, no
thermal or doping carriers — §9), so at low field it is transparent at THz:
T → 1 (the interband σ_univ response sits at k = ħω/2ħv_F ≈ 6×10⁻⁴ a.u., far
inside one mesh cell; the universal 2.3 % is *not* available to the mesh at
3.4 THz). As the field grows, the created pairs absorb: coherently only their
creation energy (~2v_F k⊥ ≈ 0.1 eV per pair → a few per cent of the pulse energy
at 100 kV/cm), with phonon scattering also the intraband (Drude) energy they
acquire while accelerated to v_F A₀ ≈ 0.7 eV. **For the intrinsic sheet the
model predicts THz-induced *absorption* that grows with the field** — ΔT from
≈ 0 at 1–10 kV/cm to −(several %) coherent and up to −(tens of %) with
dissipation at 100 kV/cm — as observed in undoped graphene [Tani et al. 2012].
The self-induced *transparency* (bleaching) reported for doped CVD graphene
[Hwang 2013; Paul 2013; Mics 2015] is the Drude-weight reduction of
*pre-existing* carriers by heating and requires a doped / finite-temperature
initial occupation, which this solver does not have yet (§9). The numbers from
the local runs are in §7.

For comparison, the near-IR interband regime (0.8 eV, `--field acos2`) at the
same fields is perturbative: pulse area 0.25 rad at 100 kV/cm ⇒ ΔA/A ≈ −θ²/12
≈ −0.5 %, ΔT ~ 10⁻⁴.

## 2. Readiness matrix (graphene, in-plane field, 1–100 kV/cm)

| effect | status | note |
|---|---|---|
| EPM ground state (π/π* Dirac pair) | ✅ **fixed** | 43-PW basis (`epm_pw_cutoff_ry = 29.4`); the old 7-PW basis is gapped at K by 0.21 eV (§4) |
| coherent velocity-gauge SBE (CF4 + Houston) | ✅ | 2 bands; the VG "basis-edge" monitor fires for ANY excitation (top band = the conduction band) — expected |
| **2D-sheet self-consistent field** (`yn_sbe_sheet_field`) | ✅ **NEW** | radiation reaction in the single-cell driver: `E_tot` = transmitted field, ledger in the local field, checkpointed (§3) |
| e-ph optical E2g/A1′ + acoustic (inter-k ring) | ✅ | Piscanec/Lazzeri; Hwang–Das Sarma acoustic, TF-screened |
| 2D Rana Auger / carrier multiplication (ring) | ✅ | R−G on the gathered n, p |
| **two-temperature Coulomb sector** (`yn_sbe_rana_te`) | ✅ **NEW** | T_e and quasi-Fermi levels from the cone moments; R−G, Q_TF, plasmon line at T_e; lattice = e-ph bath (cooling via phonons); `*_sbe_te.data` (§5) |
| **2D collisional-memory analog** (`yn_sbe_colmem`, `_pop`) | ✅ | phonon lines for e-ph, 2D Dirac-plasmon line for the Rana source (`wiki/10` §8.11) |
| Option A dressed reference (`yn_sbe_dressed_ref`) | ✅ | removes the Dirac-point rotation background |
| 2D-sheet Σ^HF (`yn_sbe_coulomb`) | ✅ (off by default) | renormalizes v_F; enable deliberately |
| impact ionization / eeh / Kuhn–Zurek | 🚫 physics | gapless cone (error stop) |
| doped / finite-T initial occupation; hot phonons; substrate in the sheet field | — limits | §9 |

## 3. Field before → field after → transmission

The driver writes `E_ext` (incident) and, with `yn_sbe_sheet_field = 'y'`,
propagates in the **local** field of the sheet (Hartree a.u., Z₀ = 4π/c):

    E_t = E_inc − (Z₀/2) J_s,   J_s = −Jm·L_z,   dA_ind/dt = −(2π/c) L_z Jm,
    A_tot = A_ext + A_ind,      E_tot = E_ext + (2π/c) L_z Jm

so the `E_tot`/`Ac_tot` columns of `*_sbe_rt.data` **are** the transmitted field
and the energy ledger is the work of the local field (`wiki/12` §5). Then

    T = ∫E_t²/∫E_inc²,  R = ∫E_r²/∫E_inc² (E_r = E_t − E_inc),  A = 1 − T − R,

fluence-integrated (Parseval-exact; a single-FFT-bin "T at the carrier" is not
bounded by 1 for a reshaped pulse and is not used). `transmission.py` detects
the self-consistent mode (`SC`), prints the deviation of `E_tot` from the
boundary-condition reconstruction on `Jm` as a consistency number (`dEt`, the
explicit-Euler lag), and — for runs without the flag (`pert`) — the
radiation-reaction term `S_rr = (Z₀/2)∫J_s²/F = A_E − A` that tells whether the
perturbative estimate is trustworthy. Linear universal sheet: T = 0.97746,
A = 0.02241 (`tests/test_sheet_transmission.py`).

## 3a. Basis: the velocity-gauge f-sum rule and the pure-gauge restoration

In the velocity gauge every electron of the filled π band carries the diamagnetic
current A·N_e/V; a *complete* basis cancels it exactly (a uniform A is a pure
gauge), an `nstate`-band basis only to the fraction S = ⟨Σ_m 2|p_nm|²/Δε⟩ (0.70 for
2 bands, 0.90 for 4, 0.964 for 8, 0.970 for 16 — the rest sits > 10 eV up). The
remainder η N_e A/V is a reactive current ∝ E/ω: negligible in the near-IR,
decisive at 3 THz, where the bare 2-band sheet reflects 85 % (a plasma mirror) and
even 16 bands leave R ≈ 8 %. The solver prints S and η at start-up. With
`yn_sbe_vg_sumrule = 'y'` it subtracts, at every step, the **adiabatic ground-state
current of the same truncated H_k(A(t))** (one ZHEEV per k): identically zero in a
complete basis, exact at every A and for any population, no fitted quantity
(wiki/12 §6a, `tests/test_vg_sumrule.f90`). A linear static form (−η N_e A) was
tried first and withdrawn — it over-corrects at 100 kV/cm and the anti-inductive
sheet runs away under the self-consistent field (T > 1). `sumrule_check.py` reports
S, η and the residual A-projection of any run (must be ≈ 0). With the
restoration the transmission is independent of `nstate` (2/3/4/8 agree at 24²,
§7.1–7.2), so production runs on the cheapest basis, `nstate = 2`.

## 3b. Doped / finite-temperature initial state (a real sample)

`sbe_ef_ev` (Fermi level from the Dirac point, eV) and `sbe_temp_init_k` set the
initial occupation to occ_max·f_FD(ε; μ, T) instead of integer filling. That is
what turns the sheet from an intrinsic semimetal into a *metal*: a partially
filled band, a Drude weight D = 2k_BT ln[2cosh(μ/2k_BT)] → E_F, an intraband
conductivity, and the 60–70 % THz transmission real CVD graphene has. Three
things go with it:

- **The pure-gauge reference stays undoped.** The f-sum-rule restoration of §3a
  subtracts the adiabatic ground-state current of the truncated basis. Evaluated
  on the *doped* occupation it would also subtract the physical intraband
  current, which is the whole point of a doped run. The solver therefore keeps
  `gs%occup_ref` (the undoped, T→0 filling) for that subtraction alone. The
  doped carriers' own truncation error is smaller than their Drude weight by
  η/⟨∂²ε/∂k²⟩ ≈ 0.28/13.5 ≈ 2 %.
- **The mesh has to resolve the Fermi surface.** k_F = E_F/ħv_F must exceed a few
  mesh spacings, else the density is set by a handful of points. At E_F = 0.2 eV
  (k_F = 0.0167 a.u.) that means N ≳ 147 (9 points per valley inside the Fermi
  disc, 3 % density error); N = 24 gives *zero* and 1.8×10¹⁰ cm⁻² instead of
  3.4×10¹². The solver counts the partially occupied k-points at start-up and
  warns below 20; `plot_occupation.py` draws the initial occupation of the doped
  and the undoped sheet on the actual mesh and prints the sheet density, k_F and
  that count — run it before committing to a mesh. Cheap alternative for a smoke test: raise E_F (E_F = 0.6 eV is
  resolved at N = 48).
- **The dressed reference and the basis-edge monitor follow the doped baseline.**
  `dressed_ref_delta` now takes the initial occupation vector (adiabatic rotation
  of ρ₀ minus ρ₀, still trace-neutral), and the VG basis-edge warning measures
  the *excess* over the initial occupation — a doped metal has its top band full
  at every k inside the Fermi surface, which is not a basis failure.

## 4. Level check — the Dirac cone the SBE runs on

`python3 tests/test_graphene_dirac_levels.py`:

| basis (`epm_pw_cutoff_ry`) | PW | gap at K | v_F | Γ bottom / M dip |
|---|---|---|---|---|
| 2.94 (old x11 input) | 7 | **0.2125 eV — spurious** (Python = Fortran bandpath) | — | −8.5 / −2.8 eV |
| **29.4 (x11 + x14 now)** | 43 | 6.5×10⁻⁶ eV (Python), 0.0000 (Fortran) | 0.960×10⁶ m/s | −7.78 / −2.70 eV |

The 7-vector set is not closed under the little group C₃ᵥ of K, so the
truncation breaks the symmetry protection of the Dirac degeneracy. Also: the
half-shifted Monkhorst–Pack mesh (2i−N−1)/2N contains K = (2/3, 1/3) **only for
odd multiples of 3** (…, 147, 153); 12, 24, 150 straddle K by half a step —
hence nk = 147 here.

## 5. Population saturation and the two-temperature model

`tests/test_rana_saturation.f90`: the Coulomb balance R = G holds at the intrinsic
density n_i(T) = (π/6)(k_BT/ħv_F)² = **8.08×10¹⁰ cm⁻² at 300 K** (∝ T²); the
CPTP channel relaxes to it monotonically from above (Auger) and below (carrier
multiplication). Evaluated at the lattice temperature this under-estimates the
plateau of a hot plasma by (T_e/T_L)². With `yn_sbe_rana_te = 'y'` the common
carrier temperature and the quasi-Fermi levels are **read from the first two
moments (n, p, energy) of the gathered cone populations** at every ring step
(`dirac_fit_te`; `tests/test_dirac_te_fit.f90` recovers T to 10⁻⁴ from
explicit-mesh moments); R−G, Q_TF and the plasmon line run at T_e, the phonon
Bose factors at the lattice T — the two-temperature model with cooling through
the e-ph channel, without a separate T_e rate equation. `*_sbe_te.data`
records t, T_e, μ_c, μ_h, n, p, n_i(T_e), T_bath; `saturation_check.py` plots it.

## 6. Run

**The production scans, shipped ready to run on the cluster.** Three input sets sit
in the exercise; they are the server runs this exercise exists for.

| set | mesh | occupation | variants | fields [kV/cm] | what it gives | cost |
|---|---|---|---|---|---|---|
| `prod_nk300_doped/` | 300² | E_F = 0.2 eV, 300 K | coh | 1…1000 (7) | **the converged T(E₀) curve — the headline result** | O(N_k): minutes/field per node, few core-hours total |
| `prod_nk300_intrinsic/` | 300² | undoped | coh | 1…1000 (7) | the control the doping is measured against | same |
| `prod_nk147/` | 147² | E_F = 0.2 eV, 300 K | diss, mem | 1…300 (5) + dark | τ, the mean free path, T_e, the absolute absorption | ring is O(N_k²): ≈7.5 h/run on 48 threads |

```bash
cd samples/exercise_x14_graphene_self_induced_transparency
# 1) the converged transmission curve (cheap, this is the main result)
for d in prod_nk300_doped prod_nk300_intrinsic; do
  ( cd $d && cp ../run_scan.sh . && OMP_NUM_THREADS=48 SALMON=../../../build/salmon bash run_scan.sh )
done
python3 field_scan_plot.py --doped 'prod_nk300_doped/runs/*/graphene_sit_sbe_rt.data'                            --intrinsic 'prod_nk300_intrinsic/runs/*/graphene_sit_sbe_rt.data'                            --t-meas 0.60 0.70 --n-sub 1.65 --out T_of_field_nk300.png
python3 plot_occupation.py prod_nk300_doped/graphene_sit --ef-ev 0.2   # pre-flight: 140 partial k-points
python3 drift_saturation.py prod_nk300_doped/graphene_sit --ef-ev 0.2  # the saturation curve vs continuum

# 2) the dissipative half (tau and the absolute absorption)
cd prod_nk147 && cp ../run_scan.sh . && OMP_NUM_THREADS=48 SALMON=../../../build/salmon bash run_scan.sh
python3 ../drude_check.py 'runs/*/graphene_sit_sbe_rt.data' --t-meas 0.60 0.70 --n-sub 1.65
python3 ../saturation_check.py runs/E100kVcm_mem/graphene_sit runs/dark_mem/graphene_sit --plot
```
Everything can also be driven end to end by `run_field_scan.sh` with `NK`, `EF`,
`FIELDS` and `VARIANTS`. A 24² smoke set is in `smoke_nk24/`.

**Why 300² for the transmission and 147² for the dissipators.** The Fermi surface
must be resolved (§3b: k_F ≳ 3 mesh spacings means N ≳ 280 at E_F = 0.2 eV), and at
300² the drift-saturation curve is converged onto the continuum (§7.12). The unitary
propagation is O(N_k), so that mesh is cheap. The graphene ring is **O(N_k²)**, so
the same mesh costs (90000/21609)² ≈ 17× the 147² dissipative run; the dissipative
production therefore runs at 147², where the Fermi surface is marginal but the
scattering physics is not mesh-critical. A converged dissipative 300² scan needs MPI
over k across nodes (`wiki/11`).
Variants: `coh` (no dissipation), `diss` (ring e-ph + acoustic + Rana at the lattice
T, Markovian — as x11), `mem` (diss + 2D colmem analog + dressed reference + T_e).
All carry the sheet field (`--no-sheet` to switch it off) and the velocity-gauge
**pure-gauge restoration** `yn_sbe_vg_sumrule = 'y'` (`--no-sumrule` to switch it
off; §3a). `--field acos2 --hw-ev 0.8 --cycles 8` gives the near-IR pulse of the
first x14 version. Other switches: `--nstate 2` (default; ≤ 4 at dt = 0.1 fs, ≥ 8 needs dt ≤ 0.05 fs — §7.2),
`--pol x|y` (in-plane polarisation; the beam is at normal incidence, k ∥ z, so E is
always in the sheet plane), `--n-layers 2` (two electronically decoupled sheets in
the same local field — incoherent/large-angle-twisted bilayer; `sbe_sheet_nlayers`),
`--snap-fs 50` (k-resolved level-population snapshots for `plot_levels.py`),
`--ef-ev 0.2 --temp-init-k 300` (a doped sheet with Drude carriers, §3b).

**The headline figure in one command** (doped + intrinsic scans on the same mesh,
then the T(E₀) / σ(E₀) plot of §7.9):

```bash
NK=147 EF=0.2 OMP_NUM_THREADS=48 SALMON=../../build/salmon bash run_field_scan.sh
# -> scan_nk147_ef0.2/doped_vs_intrinsic.png
```
`NK` must resolve the doping: k_F = E_F/ħv_F needs ≈ 3 mesh spacings, so NK ≳ 280
at E_F = 0.2 eV, 140 at 0.4, 93 at 0.6 (wiki/12 §4a.0). NK = 147 at 0.2 eV is the
cheapest setting that resolves the Fermi surface at all and is what the shipped
figure uses; the Drude weight is then ≈ 30 % low, a field-independent scale error
that leaves the shape of T(E₀) intact. Coherent runs only, so the whole pair of
scans is minutes per field.

```bash
python3 ../plot_levels.py runs/E100kVcm_mem/graphene_sit --times 100,150,200,300   # band populations vs t + k-maps around K
python3 ../sumrule_check.py graphene_sit --runs 'runs/*/graphene_sit_sbe_rt.data'  # basis f-sum rule + residual reactive current
python3 ../drude_check.py runs/*/graphene_sit_sbe_rt.data --t-meas 0.60 0.70 --n-sub 1.65  # doped sheet: D, tau, mean free path
```

**Cost (measured 2026-09-04, 4 threads).** Coherent runs: seconds to minutes at
147². Dissipative runs: the graphene ring (e-ph inter-k + Rana) is O(N_k²) and
exponential-bound — ≈ 60 ms/step at 24² with the THz drive (all cone states are
sources) ⇒ ≈ 84 s/step at 147² on 4 threads ⇒ **≈ 7 s/step on 48 threads ⇒ ≈ 7.5 h
per 3844-step run**. Levers: `--dt-fs 0.2` (halves it; the CF4 unitary is exact
per step, check `nex` against 0.1 fs), MPI ranks over k (the ring gather is
MPI-parallel, `wiki/11`). Essential set = 100 kV/cm × {coh, diss, mem} + `dark_mem`.

**Do not economise on `--tail-fs`.** The 100 fs ring-down after the 284 fs transient
looks like 25 % of the cost for nothing, but the sheet current is still ringing when
the drive ends, and every reported number is a fluence integral. Measured at 72²,
E_F = 0.6 eV, coherent, with and without it:

| E₀ [kV/cm] | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|
| T, `--tail-fs 0` | 0.6846 | 0.6807 | 0.6439 | 0.6454 | 0.7495 |
| T, `--tail-fs 100` | 0.7093 | 0.7053 | 0.6679 | 0.6713 | 0.7739 |
| A, `--tail-fs 0` | 0.0589 | 0.0593 | 0.0592 | 0.0554 | 0.0775 |
| A, `--tail-fs 100` | 0.0095 | 0.0100 | 0.0112 | 0.0035 | 0.0287 |

T comes out 0.024 low at every field (the shape survives, the absolute value does
not) and the absorption of a *coherent* run — which has no dissipation at all and
must be near zero away from pair creation — is inflated six-fold, because the
transmitted fluence that is still in flight at the last step is simply not counted.

**Dissipative runs are the exception, and for the reason you would guess:** the ring
damps the current, so by the end of the drive there is little left in flight. Same
mesh and doping, `diss`:

| | T, no tail | T, +100 fs | A, no tail | A, +100 fs |
|---|---|---|---|---|
| 10 kV/cm | 0.77773 | 0.77864 | 0.20363 | 0.20182 |
| 100 kV/cm | 0.82761 | 0.82785 | 0.15297 | 0.15248 |

— a shift of 0.0009 and 0.0002 against the coherent run's 0.0246. So the tail is
cheap insurance in general and indispensable for `coh`; if the ring is on and the
budget is tight, it is the one place the 25 % can honestly be saved.

`transmission.py` also refuses to analyse a record shorter than the run's own `nt`
(`--allow-partial` overrides), so a job that is still going cannot quietly contribute
a number.

## 7. Local validation (2026-09-04, 4 threads, container)

All runs: DAST single-cycle 3.36 THz proxy, x-polarised, normal incidence,
sheet field ON, `coh` unless stated; T, R, A fluence-integrated from
`transmission.py`; "resid." = residual A-projection of the sheet current in
units of N_e A/A_2D (`sumrule_check.py`, must be ≈ 0).

**7.1 The basis artifact and its removal (24², 100 kV/cm).**

| nstate | restoration | T | R | A | A_tot/A_ext | after the pulse | resid. |
|---|---|---|---|---|---|---|---|
| 2 | off | 0.147 | 0.851 | 0.002 | 0.20 | A_ind relaxes | 0.2935 |
| 4 | off | 0.646 | 0.321 | 0.033 | 0.53 | relaxes | — |
| 8 | off | 0.937 | 0.013 | 0.049 | 0.90 | relaxes | — |
| 16 | off | 0.945 | 0.005 | 0.050 | 0.92 | relaxes | — |
| 2 | linear static −ηN_eA (withdrawn) | 1.213 | 0.499 | −0.712 | 1.71 | **runaway** (A_ind, J grow together) | −0.023 |
| 8 | linear static (withdrawn), `diss` | — | — | — | 1.40 | **runaway** | — |
| 2 | pure gauge (Eq. 10) | 0.999996 | 1.3e-6 | 0.000002 | 1.000 | stationary | 0.0000 |
| 3 | pure gauge | 0.999996 | 1.4e-6 | 0.000003 | 1.000 | stationary | 0.0000 |
| 4 | pure gauge | 0.999996 | 1.4e-6 | 0.000003 | 1.000 | stationary | 0.0000 |
| 8 | pure gauge | 0.972 | 3.4e-4 | 0.027 | 0.986 | slow DC drift | 0.0001 |

The 24² mesh has no interband channel at 14 meV (smallest π–π* gap 0.78 eV, K
half a spacing off the mesh, no k-point in the Landau–Zener tube), so the
physical answer at this resolution is T ≃ 1 — which nstate = 2, 3 and 4 now
give identically. The 2.7 % of nstate = 8 is population leaked into bands
4–8 (22–58 eV above π; 5.6×10⁻⁶ per cell carrying 2.5×10⁻⁴ eV per cell, the
whole ledger) — a 14 meV field cannot do that; see 7.2 for the time-step test.
Linear checks of the mechanism (24², nstate = 2, 1 kV/cm, no sheet): the SBE
current equals the exact adiabatic response of the truncated H(A) to 1×10⁻⁴
(dt 0.1/0.05/0.02 fs), in-phase coefficient η = 0.2995 = Eq. (9); at 0.1 and
0.4 eV drive: 0.2985 / 0.2808 vs the dispersive prediction 0.2985 / 0.2806.

**7.2 Time step at nstate = 8 (24², 100 kV/cm).** With dt = 0.05 fs the nstate = 8
run gives T = 0.999985, A = 1.4×10⁻⁵, ledger 1.1×10⁻⁷ eV per cell — the same
answer as nstate = 2/3/4, and 2400× less deposited energy than at dt = 0.1 fs.
The 2.7 % at dt = 0.1 fs was the S4/CF4 step applied to bands 34–90 eV above
the cone (13.7 rad of phase per step): population leaks into them. Rule: at
dt = 0.1 fs use nstate ≤ 4; nstate ≥ 8 needs dt ≤ 0.05 fs. **Production recipe
therefore: nstate = 2 + pure-gauge restoration** (`make_inputs.py` default) —
the restoration makes 2, 3 and 4 bands agree to 10⁻⁶ in T, the THz physics
lives on the cone, and the ring cost is lowest.

**7.3 147² (K on the mesh), nstate = 2, pure gauge, coherent — the field scan.**
Static η = 0.2777 (x = y: the 147² mesh restores the hexagonal isotropy the 24²
mesh breaks at the 10 % level) — removed by the restoration; what remains in the
A-projection is the *physical* Drude response of the pairs the field creates.

| E₀ [kV/cm] | T | R | A | ledger A_E | pairs after the pulse [cm⁻²] | peak diabatic n_c (per cell) | A_tot/A_ext | Drude resid. |
|---|---|---|---|---|---|---|---|---|
| 1 (K arbitrary, see below) | 0.816 | 1.82 | −1.63 | −2.08 | 1.7×10³ | 8×10⁻¹¹ | 0.002 | 60.6 |
| 10 (K arbitrary) | 0.9569 | 0.0348 | 0.0083 | 0.0083 | 2.1×10¹⁰ (4.2×10⁻⁵ per cell) | 2.2×10⁻⁴ | 0.818 | 0.0100 |
| 30 | 0.9936 | 0.0039 | 0.0025 | 0.0025 | 3.5×10¹⁰ | 2.4×10⁻³ | 0.938 | 0.0029 |
| 100 | 0.9824 | 0.0043 | 0.0134 | 0.0134 | 6.7×10¹¹ (3.5×10⁻⁴ per cell) | 2.0×10⁻² (99 % virtual, returns) | 0.941 | 0.0048 |
| 1 (K averaged, rerun) | 0.999995 | 1.2×10⁻⁶ | 4×10⁻⁶ | 4×10⁻⁶ | 0 net | 9.3×10⁻⁵ (= baseline) | 1.000 | 0.0000 |
| 10 (K averaged, rerun) | 0.999993 | 2.4×10⁻⁶ | 4×10⁻⁶ | 4×10⁻⁶ | 0 net (returns to the 0.9×10⁻⁴ baseline) | 3.8×10⁻⁴ | 1.000 | −0.0001 |
| 30 (K averaged, rerun) | 0.9995 | 2.4×10⁻⁴ | 2.6×10⁻⁴ | 2.6×10⁻⁴ | 2×10⁹ net (1×10⁻⁶ per cell) | 2.6×10⁻³ | 0.995 | 0.0005 |
| 100 (K averaged, rerun) | 0.9867 | 0.0030 | 0.0103 | 0.0103 | 4.0×10¹¹ net (3.0−0.9×10⁻⁴ per cell; 0.9×10⁻⁴ = the two half-filled K points) | 2.0×10⁻² | 0.958 | 0.0041 |
| 10, 150² (K off mesh) | 0.9925 | 0.0004 | 0.0071 | 0.0059 | 1×10¹⁰ (5×10⁻⁶ per cell) | 3.8×10⁻⁴ | 1.001 | −0.0003 |
| 100, 150² (K off mesh) | 0.9563 | 0.0067 | 0.0370 | 0.0365 | 1.1×10¹² (5.6×10⁻⁴ per cell) | 2.0×10⁻² | 0.956 | 0.0035 |
| 100, 300² (K off mesh) | 0.9688 | 0.0046 | 0.0266 | 0.0264 | 7.7×10¹¹ (4.0×10⁻⁴ per cell) | 2.0×10⁻² | 0.950 | 0.0041 |

**The Dirac point on the mesh (found in this scan).** With K exactly on the
mesh the two π/π* levels at K are degenerate and the integer ground-state
filling (band 1 = 2, band 2 = 0) picks LAPACK's arbitrary basis inside the
degenerate pair: a broken-symmetry state with a velocity expectation of order
v_F at that k-point — a spurious current ~ 2v_F/N_k per valley that does not
depend on the field. At 100 kV/cm it is buried under the physical response;
at 10 kV/cm it is the whole "Drude residual" (0.010) and the 3.5 % reflection;
at 1 kV/cm, together with the sheet field, it acts as a relay (the pure-gauge
reference jumps by ±2v_F w_K at A = 0 and the sheet pins A_tot to zero: the
nonsense row above, R > 1). Fix (`gs_info_ssbe.f90`, all materials): a
degenerate, partially filled level group is given its **group-average
occupation** — the T → 0 density matrix, the identity on the block, current-free
and rotation-invariant; the pure-gauge reference becomes continuous. The rows
marked "K averaged" are the rerun with this fix; the 150² rows have K half a
spacing off the mesh (no degeneracy at all) and bracket the result from the
other side. Below ~30 kV/cm the Landau–Zener tube is thinner than one mesh
cell (0.01 / 0.05 / 0.3 / 1.6 / 9.6 cells at 1 / 3 / 10 / 30 / 100 kV/cm on
147²), so the pair creation there is Eq. (8) of wiki/12 analytically —
1.5×10⁹ / 4.7×10¹⁰ / 2.5×10¹¹ / 1.5×10¹² cm⁻² per passage at 1 / 10 / 30 / 100
kV/cm — rather than mesh-resolved; the 30 and 100 kV/cm rows are resolved.

At 100 kV/cm the sheet current is 96 % anti-correlated with A_tot *after* the
restoration: that is the kinetic inductance of the ~10¹² cm⁻² pairs created by
the Landau–Zener sweep through K (Eq. 8 of wiki/12 predicts 1.5×10¹² per
passage; two passages with Stückelberg interference give the 6.7×10¹¹ left),
their creation energy is the 1.3 % absorbed (ledger = fluence to 3 digits), and
the induced field relaxes after the pulse (A_ind −6.1 → −4.0×10⁻³ a.u. over the
tail, J decaying) — a stable, passive sheet. **Expected range of the
transmission change for the intrinsic sheet, coherent limit: T ≈ 1.000 up to
30 kV/cm (induced absorption < 0.1 %), T = 0.97 ± 0.01 at 100 kV/cm (induced
absorption 1–4 % across the 147²/150²/300² meshes and the two polarisations,
2.7 % on 300²; reflection ≤ 0.7 %)** — i.e. the intrinsic sheet *darkens* by a
few per cent at the top of the DAST range and does not bleach. With phonon
dissipation the carriers are also heated while accelerated, and the
`diss`/`mem` production runs will raise A (README §8). The mesh spread is the
transverse sampling of the Landau–Zener strip (0.7 cells at 147², 1.4 at 300²);
convergence to ≲ 0.5 % needs N ≳ 600 or a K-refined mesh.

**7.4 Polarisation and two decoupled layers (147², 100 kV/cm).** The `input`
field file carries its own three Cartesian components — `epdir_re1` is *not*
applied to it (a first "y" run silently reproduced the x result to six digits;
`make_inputs.py --pol y` now writes the waveform into the A_y column). True y
run (147², K averaged, 100 kV/cm): T = 0.9735, R = 0.0056, A = 0.0209 against
x: T = 0.9867, R = 0.0030, A = 0.0103. This factor of 2 is **not** the crystal's
anisotropy (see below) but the mesh's: the Landau–Zener tube is a strip
2k_⊥ = 7.5×10⁻³ a.u. wide (0.7 mesh cells at 147²) and 2A₀ = 0.12 long (12
cells), and how many mesh points fall inside it depends on its orientation on
the hexagonal mesh — the same sampling spread as 147² vs 150² for one
polarisation (1.0 % vs 3.7 %). The transverse tube width sets the convergence
requirement: ≳ 3 cells across the tube means N ≳ 600 at 100 kV/cm (or a
K-refined mesh); the 300² run below is the first step of that series.

Two electronically decoupled sheets in the same local field
(`sbe_sheet_nlayers = 2`, 147², K averaged, 100 kV/cm): T = 0.9585, R = 0.0112,
A = 0.0304 against the single layer's T = 0.9867, R = 0.0030, A = 0.0103. The
naive estimate T₁² = 0.9736 is 1.5 % (absolute) too high: in the coherent
non-linear regime the two layers do not simply multiply — they share one local
field (A_tot/A_ext = 0.92 instead of 0.96) and the re-radiated field reshapes
the waveform each layer sees. For the *linear* regime of real CVD samples the
sum-of-conductances formula applies (d ≪ λ): with z = Z₀σ complex,

    T₁ = |2/(2+z)|²,   T₂ = |2/(2+2z)|²,   T₁²/T₂ = 16|1+z|²/|2+z|⁴
                                                  = 1 + ½[(Im z)² − (Re z)²] + O(z³)

**The sign of the T·T error depends on what kind of sheet it is.** For a purely
dissipative (real σ) sheet T₁² is LOW: 0.03 % at T₁ = 0.98, 0.6 % at 0.89,
9.6 % at the sample's z = 0.565 (T₁²/T₂ = 0.904), 25 % at T₁ = 0.40. For a
purely reactive (inductive, imaginary σ) sheet — which is what the coherent
doped runs here are — T₁² is HIGH: the same |z| = 0.565 gives T₁²/T₂ = 1.133,
+13 %. A doped THz sheet is a mixture, so |z| alone does not even fix the sign;
only |z| ≤ 0.15 (T₁ ≥ 0.87) keeps the error under 1 % either way. The
substrate's own Fresnel factor must also be divided out *before* squaring. For Bernal (AB) bilayer or small-angle moiré the electronic structure
itself changes and neither estimate applies.

**The same comparison on a DOPED sheet (72², +100 fs tail, coherent), which is what
a real CVD bilayer is** — figure `wiki/figures/graphene_layers_1_vs_2.png`, rebuilt by

```bash
python3 layers_plot.py --set 0.6 'ef0.6_L1/runs/*/*_rt.data' 'ef0.6_L2/runs/*/*_rt.data' \
                       --set 0.4 'ef0.4_L1/runs/*/*_rt.data' 'ef0.4_L2/runs/*/*_rt.data' \
                       --t-meas 0.60 0.70 --out layers.png
```


| E₀ [kV/cm] | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|
| E_F = 0.6 eV: T₁ | 0.7093 | 0.7053 | 0.6679 | 0.6713 | 0.7739 |
| T₂ | 0.4362 | 0.4341 | 0.4128 | 0.3778 | 0.4715 |
| T₁²/T₂ | 1.153 | 1.146 | 1.081 | 1.193 | 1.270 |
| E_F = 0.4 eV: T₁ | 0.8487 | 0.8422 | 0.7939 | 0.8358 | 0.8772 |
| T₂ | 0.6493 | 0.6434 | 0.5847 | 0.5894 | 0.7226 |
| T₁²/T₂ | 1.109 | 1.103 | 1.078 | 1.185 | 1.065 |

T₁² is HIGH at every field and both dopings — the reactive branch, because a coherent
doped sheet at 14 meV is an inductor (band-averaged z = 0.017 − 1.195i at 0.6 eV).
Inverting a measured bilayer transmission as √T₂ therefore understates the monolayer
conductance by ~5 %. Two further checks: `--predict-layers 2`, which uses the run's
own σ(ω), reproduces T₂ to 1.5 % in the linear regime (against 14 % for the
single-frequency formula) and degrades to 27 % at u = 5.6, always predicting the
stack too bright; and bin by bin across the driven band the two-layer run returns
exactly twice the one-layer σ(ω) (median ratio 1.997 and 2.007), which is the check
that `sbe_sheet_nlayers` itself is right. The band-averaged Re σ ratio is 1.83 rather
than 2 because the two-layer transmitted field is filtered more strongly where |σ| is
largest — a property of that diagnostic, not of the solver.

**Answer to "does the direction matter?"** For the linear and the χ⁽³⁾ response
the hexagonal sheet is isotropic in-plane (any 2nd- and 4th-rank tensor of C₆ᵥ
is); anisotropy (zigzag vs armchair) enters at χ⁽⁵⁾ through trigonal warping,
i.e. at relative order (A₀/|K−M|)² ≈ (0.06/0.78)² ≈ 0.6 % of the already-small
non-linear part — below 10⁻⁴ in T at 100 kV/cm. On a mesh the x/y difference is
a resolution artifact: the 24² mesh breaks C₆ (S_x = 0.700 vs S_y = 0.728), the
147² mesh does not (0.7223 both).

**7.5 Near-IR πα check (147², acos2 0.8 eV, 8 cycles, nstate = 2, pure gauge, sheet).**
1 kV/cm: T = 0.98128, R = 6.8×10⁻⁴, A = 0.01804 (ledger 0.01804); 100 kV/cm:
T = 0.98136, A = 0.01796. The sheet absorbs 0.80 of the universal πα value
(A = 0.0224; Re σ/σ_univ = 0.90 over the pulse band) — the 3.2 mesh points per
resonance-shell radius and the finite-basis matrix elements at k_res, not a
solver error (the reactive residual is −0.013, i.e. the interband *capacitive*
response above resonance, the correct sign). The field dependence is the
predicted coherent bleaching: ΔA/A = −0.4 % (θ²/12 ≈ −0.5 %), ΔT = +8×10⁻⁵ —
near-IR self-induced transparency at 1–100 kV/cm is a 10⁻⁴ effect.

**7.6 Ring with T_e (24², nstate = 2, 100 kV/cm, `diss` / `mem`).** The 24² mesh
has no real pair channel at this field (coherent: A = 2×10⁻⁶), so whatever the
dissipative runs absorb is what the ring does to the *virtual* dressing (peak
diabatic n_c = 2.3×10⁻² per cell, all of it returning in the coherent run).
`diss` (Markovian ring, lattice T): T = 0.9988, A = 0.0012; ring-visible density
7.6×10⁹ cm⁻² generated out of the dressing, Rana ledger +1.9×10⁹ (net
multiplication, n < n_i) — the "dephasing ionization" of `wiki/10` §8.7 in its
graphene form. `mem` (2D colmem analog + dressed reference + T_e): T = 0.999996,
A = 3×10⁻⁶ — identical to the coherent run; ring-visible density 5×10⁴ cm⁻²
(10⁵ times less), Rana ledger +3.5×10⁶: the memory filters remove the
fabricated generation completely, as they do for Si/GaAs/CdS (`wiki/10` §8.11).
The two-temperature fit is meaningless when there are no carriers (it returned
T_e ≈ 4×10⁴ K for 10⁸ cm⁻²); the solver now holds T_e at the lattice value while
n + p < 10⁻³ n_i(T_lattice), so `*_sbe_te.data` reads 300 K there. Real carrier
heating (T_e in the 10³ K range, cooling on the phonon timescale) needs the
147² production runs, where pairs are actually created.

**7.7 Above 100 kV/cm — the intrinsic sheet keeps darkening (147² K-averaged and
150², nstate = 2, pure gauge, self-consistent sheet, coherent).**

| E₀ [kV/cm] | T (147²) | R | A | Re σ/σ_univ | T (150²) | A (150²) |
|---|---|---|---|---|---|---|
| 1 | 0.999995 | 1.2e-6 | 4e-6 | 0.000 | — | — |
| 10 | 0.999993 | 2.4e-6 | 4e-6 | 0.000 | 0.9925 | 0.0070 |
| 30 | 0.9995 | 2.4e-4 | 2.6e-4 | 0.03 | — | — |
| 100 | 0.9867 | 0.0030 | 0.0103 | 0.78 | 0.9563 | 0.0370 |
| 200 | 0.9505 | 0.0099 | 0.0396 | 2.83 | 0.9515 | 0.0389 |
| 300 | 0.9336 | 0.0178 | 0.0486 | 4.02 | 0.9404 | 0.0459 |
| 500 | 0.9369 | 0.0200 | 0.0430 | 3.94 | 0.9230 | 0.0535 |
| 1000 | 0.9163 | 0.0335 | 0.0503 | 5.52 | 0.9103 | 0.0546 |

Three things settle here. **(i)** The transmission of the *intrinsic* sheet falls
monotonically over three decades of field and never rises: there is no
self-induced transparency of undoped graphene in this model. **(ii)** The
absorbed fraction saturates near 5 % above ~300 kV/cm, which is the analytic
Landau-Zener value (wiki/12 §7: the pairs created per passage scale as E^{3/2}
and their creation energy as √E, so the absorbed *energy* scales as the fluence
and the *fraction* is field-independent). **(iii)** The two meshes, which
differ by a factor 3.7 in A at 100 kV/cm, agree to 1–2 % above 200 kV/cm: once
the Landau-Zener strip is several mesh cells wide the sampling ambiguity of
§7.3 is gone. The growing reflection (3.4 % at 1000 kV/cm) and Re σ = 5.5 σ_univ
are the Drude response of the created plasma, not an artifact — the residual
A-projection stays below 0.005 throughout.

**7.8 A doped sheet — what the measured 60 % → 70 % needs.** A sample on PET
transmitting 60 % (substrate included) is not the intrinsic sheet of §7.7. With
n_PET = 1.65 and two incoherent faces the bare substrate passes 88.3 %, so the
graphene itself accounts for 0.60/0.883 = 0.68 and 0.70/0.883 = 0.79, i.e. a
sheet conductance Z₀σ of 0.565 → 0.327, **σ = 24.7 → 14.3 σ_univ (−42 %)**.
That is a Drude conductance: σ_dc = Dτ/π with the Dirac-cone Drude weight
D = 2k_BT ln[2cosh(μ/2k_BT)] → E_F. In transport units it is R_s = 666 Ω/sq at low
field and 1153 Ω/sq at high field — an ordinary as-transferred CVD monolayer.
Splitting it across dopings, each with the τ that reproduces the same σ:

| E_F [eV] | n_2D [cm⁻²] | τ [fs] | mobility [cm²/V s] | mean free path [nm] | onset A₀ = k_F [kV/cm] |
|---|---|---|---|---|---|
| 0.1 | 8.0×10¹¹ | 127 | 11700 | 122 | 13 |
| **0.2** | **3.2×10¹²** | **64** | **2900** | **61** | **27** |
| 0.3 | 7.2×10¹² | 43 | 1300 | 41 | 40 |
| 0.4 | 1.3×10¹³ | 32 | 730 | 31 | 54 |
| 0.6 | 2.9×10¹³ | 21 | 330 | 20 | 81 |

E_F = 0.2–0.4 eV is where an ordinary sample sits; 0.1 eV would need a mobility
CVD-on-polymer does not have, and 0.6 eV a mobility that is low even for
chemically doped film. Independent corroboration: Hassanpour Amiri *et al.*,
*Doping free transfer of graphene using aqueous ammonia flow*, RSC Adv. **10**,
1127 (2020), attribute the unintentional doping of transferred CVD graphene to
ionic etch residue of **typical density 4×10¹² cm⁻²** — on the Dirac cone that is
E_F = 0.224 eV and a saturation onset of 30 kV/cm, within 12 % of what the
transmission alone gives. Their ammonia-washed (doping-free) films should follow
the *intrinsic* curve instead and darken with field: a cheap test of the mechanism. The field at which the transmission starts to rise picks
between them (last column). The E_F = 0.6 eV used in the 48² runs below is a
mesh-affordable proxy, not the sample.

`sbe_ef_ev` / `sbe_temp_init_k` now put that state into the solver (§3b). Which
of the two factors, D or τ, can give −42 %?

- **D cannot.** At *fixed* carrier density, heating lowers μ but adds thermal
  carriers; D passes a shallow minimum and comes back up. Over 300–3000 K the
  deepest excursion is 0.88 at n = 10¹³ cm⁻² and 0.90 at 3×10¹² (table in
  `drude_check.py`, pinned by `tests/test_doped_drude.py`). Carrier heating on
  its own is worth ~10 %, not 42 %.
- **τ can, and so can the drift itself.** Two field scales control it. The
  Fermi sea is displaced by the vector potential, A₀ = 6.2×10⁻⁴ a.u. per kV/cm
  for this transient, against k_F = E_F/ħv_F = 0.0167 a.u.: they are equal at
  **27 kV/cm**. Above that the displacement exceeds the Fermi radius every
  half-cycle, the drift velocity saturates at v_F, and the differential
  conductivity falls like k_F/A₀. The collisional excursion eEτ/ħ crosses k_F at
  the same place for τ = 63 fs. Simultaneously the carriers are pushed past the
  optical-phonon threshold (E₂g 196 meV, A₁′ 160 meV; v_F A₀ = 0.74 eV at
  100 kV/cm), and every emission randomises momentum — τ drops.

So the mechanism the measurement is showing is the one you suspected: the
conductivity falls because the carriers stop responding linearly (mobility /
mean free path), not because the Drude weight disappears. The onset the model
puts at ~27 kV/cm for E_F = 0.2 eV is inside the measured range, and the
predicted saturation is *sub-linear in the drift*, which is why the transmission
creeps up by ten points rather than jumping.

**7.10 At the sample's own doping (147², E_F = 0.2 eV) — the signature is the ratio
to the intrinsic control.** A 3×10¹² cm⁻² sheet without dissipators is only a weak
inductor, so the raw transmission is quiet (T = 0.968 at 1 kV/cm, not the measured
0.68 — a coherent run has no momentum relaxation and cannot produce absorption,
only reactive screening). Divide by the intrinsic control on the same mesh and the
doping carriers stand alone:

| E₀ [kV/cm] | 1 | 10 | 30 | 100 | 300 | 1000 |
|---|---|---|---|---|---|---|
| T doped (E_F = 0.2 eV) | 0.9682 | 0.9666 | 0.9595 | 0.9629 | 0.9442 | 0.9090 |
| T intrinsic | 0.999995 | 0.999993 | 0.999498 | 0.986694 | 0.933550 | 0.916273 |
| extinction added by the doping | 3.18 % | 3.34 % | **4.00 %** | 2.41 % | −1.14 % | 0.79 % |
| Re σ/σ_univ (doped) | 2.73 | 2.86 | 3.45 | 2.81 | 3.73 | 6.11 |

The doping-induced extinction **peaks at 30 kV/cm**, against the 27 kV/cm that
A₀ = k_F predicts from k_F alone (k_F = E_F/ħv_F = √(πn) is the radius of the
occupied disc, and in the velocity gauge A is literally its displacement in
reciprocal space — wiki/12 §4a.5.1), then collapses by a factor of five. At 300 kV/cm
it is *negative*: the doped sheet transmits better than the undoped one, because its
Drude response has saturated while its occupied states Pauli-block part of the
Landau–Zener pair creation. Figure: `doped_vs_intrinsic.png` (three panels: T(E₀),
the extinction ratio, σ(E₀)); the 48²/E_F = 0.6 eV proxy of §7.9 is the same physics
with a sheet conductance large enough to move the raw transmission.

**7.9 The doped sheet, calculated: transmission rises above the saturation
onset.** 48² mesh, E_F = 0.6 eV (n = 3.0×10¹³ cm⁻², 8 mesh points per valley
inside the Fermi disc), T_init = 300 K, nstate = 2, pure gauge, self-consistent
sheet, **coherent** (no dissipators at all — so everything below is mechanism 3
alone). D_fit and τ_fit come from fitting the run's own current to
dJ_s/dt = (D/π)E_tot − J_s/τ.

| E₀ [kV/cm] | T | R | A | D_fit [eV] | Re σ/σ_univ |
|---|---|---|---|---|---|
| 1 | 0.7286 | 0.2497 | 0.022 | 0.533 | 23.6 |
| 3 | 0.7284 | 0.2498 | 0.022 | 0.533 | 23.6 |
| 10 | 0.7257 | 0.2519 | 0.022 | 0.537 | 23.8 |
| 30 | 0.6990 | 0.2737 | 0.027 | 0.580 | 25.8 |
| 100 | 0.6622 | 0.3304 | 0.007 | 0.664 | 30.3 |
| 200 | 0.7202 | 0.2736 | 0.006 | 0.546 | 25.4 |
| 300 | 0.7716 | 0.2066 | 0.022 | 0.424 | 20.5 |
| 500 | 0.8122 | 0.1455 | 0.042 | 0.325 | 15.8 |
| 1000 | 0.8263 | 0.1124 | 0.061 | 0.290 | 13.4 |

The **intrinsic control on the same mesh**, same pulse, same everything except the
initial occupation:

| E₀ [kV/cm] | 1 | 10 | 30 | 60 | 100 | 200 | 300 | 500 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| T, doped (E_F = 0.6 eV) | 0.7286 | 0.7257 | 0.6990 | 0.6414 | 0.6622 | 0.7202 | 0.7716 | 0.8122 | 0.8263 |
| T, intrinsic | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.9970 | 0.9772 | 0.9519 | 0.9230 | 0.9026 |

(figure: `doped_vs_intrinsic.png`). The two curves run in opposite directions and
cross near 700 kV/cm.

Read from the bottom up, this is the measurement. **T is flat at 0.726 through
the linear regime (1–10 kV/cm), dips to 0.66 near 100 kV/cm, then rises
monotonically to 0.83.** The turn is where predicted: A₀ = k_F at 81 kV/cm for
E_F = 0.6 eV. Above it the fitted Drude weight collapses, 0.66 → 0.29 eV, and
with it the sheet conductivity, **30.3 → 13.4 σ_univ (−56 %)** — against the
−42 % (24.7 → 14.3 σ_univ) the PET measurement implies. Same sign, same
mechanism, same order. **The dip before the turn is not physics.** §7.12 shows it is the
discretization bump of a Fermi disc holding only 12 mesh points: the continuum
drift-saturation curve is monotone, and the T minimum at 60 kV/cm falls exactly on
u = A₀/k_F = 0.386, where the 48² mesh curve peaks at 1.171. Converged, the doped
sheet is flat and then brightens, with no darkening in between.

Note what is **not** in this table: no phonons, no Auger, no heating — the
dissipators are off. Mechanism 3 (current saturation on the Dirac cone) accounts
for the whole effect on its own. Mechanisms 1 and 2 add to it: heating supplies
≈10 % more, and optical-phonon emission above 196 meV shortens τ, which lowers σ
further in the same direction. So the measured rise does not need all three —
but none of them makes the sheet *darker*, and the model has no mechanism that
would.

*Drude-weight accuracy — measured convergence.* The fitted D is below the analytic
D = E_F because the mesh under-resolves ∂²ε/∂k² ∝ 1/k near the Dirac point. At
1 kV/cm (linear regime):

| mesh | E_F [eV] | partially occupied k-points | n_2D vs analytic | D_fit/E_F |
|---|---|---|---|---|
| 147² | 0.2 | 36 | 3.27 vs 3.36×10¹² (−2.6 %) | 0.659 |
| 300² | 0.2 | 116 | 3.345 vs 3.36×10¹² (−0.5 %) | **0.930** |
| 147² | 0.4 | 72 | 1.284×10¹³ | 0.894 |
| 48² | 0.6 | 8 | 3.04 vs 2.89×10¹³ (+5 %) | 0.888 |

The density needs one mesh shell inside the Fermi circle, the Drude weight three
or four. The deficit is a scale error common to all fields, so the *shape* of the
T(E₀) curve above is unaffected; absolute conductivities want k_F ≳ 3 mesh
spacings (N ≳ 280 at E_F = 0.2 eV).

**7.11 With the phonon ring on: the scattering time and the mean free path.**
Same 48² mesh and doping, `diss` (e-ph ring at a 300 K lattice), against the
coherent runs:

| E₀ [kV/cm] | τ [fs] | mean free path [nm] | T_e | Re σ/σ_univ | T | R | A |
|---|---|---|---|---|---|---|---|
| 10 | 141 | 136 | 361 K | 18.1 | 0.6442 | 0.061 | 0.294 |
| 100 | 60 | 58 | 2050 K | 14.5 | 0.7213 | 0.048 | 0.231 |
| 100, coherent | 7930 | 7613 | — | 30.3 | 0.6622 | 0.330 | 0.007 |

Two things are robust here. **(i)** Scattering turns the sheet from an inductive
mirror (R = 0.33, A = 0.007) into a Drude absorber (R = 0.048, A = 0.231) — a real
sample, not a mirror. **(ii)** T rises 0.644 → 0.721 and the high-field conductivity
is 14.5 σ_univ, against a measured 0.68 → 0.79 and 14.3 σ_univ.

**Caveat: neither the absolute τ nor its field dependence is mesh-converged.**
Repeating the scan on 72² (2.25× the k-points):

| mesh | E₀ [kV/cm] | τ [fs] | mean free path [nm] | Re σ/σ_univ | T |
|---|---|---|---|---|---|
| 48² | 10 | 141 | 136 | 18.1 | 0.644 |
| 48² | 100 | 60 | 58 | 14.5 | 0.721 |
| 72² | 10 | 33 | 31 | 10.6 | 0.779 |
| 72² | 100 | 72 | 69 | 8.5 | 0.828 |

(72² rows: the tail-inclusive repeat. A second, independent route agrees — inverting
σ(ω) bin by bin, τ = −Im σ/(ω Re σ), with no time-domain fit and no driven window,
gives 21.5 fs at 10 kV/cm and 58.0 fs at 100. Same rise, from data the fit never
touches, so the reversal against 48² is not an artifact of the fitting window.)

The magnitude moves by 4× at fixed field, and the *trend* reverses: at 48² τ falls
2.3× from 10 to 100 kV/cm, at 72² it rises by 2.2× over the same interval. An earlier
version of this file read the 48² fall as the carriers crossing the optical-phonon
thresholds (v_F A₀ = 0.74 eV at 100 kV/cm against E₂g 196 meV, A₁′ 160 meV); the
finer mesh does not reproduce it, so **that reading is withdrawn** — it was mesh
noise, not a mobility drop that the calculation can claim. The e-ph ledger rate
doubles between the two meshes for the same reason: the inter-k golden-rule sum is
unresolved, the σ_E = 0.1 eV energy window holding too few final states, so the rate
samples the mesh. What is robust is the sign
of every effect at fixed mesh when the ring is switched on: scattering shortens τ
against the collisionless run, converts reflection into absorption, lowers σ and
raises T. Quoting an absolute τ, or comparing the absolute transmission with the
measurement quantitatively, needs the ring on a Fermi-surface-resolving mesh: the
O(N_k²) production run of §6.

The carrier temperature comes from the total electronic energy (361 K and 2050 K),
*not* from the per-channel `dE_eph` column, which for a doped initial state is a
gross exchange counter: it grows at 7×10⁻⁶ eV/cell/fs from t = 0 regardless of
field, while the electronic energy of the doped sea moves by only 2.5×10⁻⁵ eV/cell
over 384 fs at 10 kV/cm and the carrier number is conserved to 0.14 %. The doped
Fermi sea is stationary under the dissipators; the ledger column is not a net loss.

*Still to run (cluster).* The same with `--ef-ev 0.2 --temp-init-k 300 --variants
diss,mem` on a Fermi-surface-resolving mesh (147² or denser, §3b), which replaces
the E_F = 0.6 eV proxy by the sample's own doping. Cost as in §6.

**7.12 Drift saturation, and separating it from the mesh.** The rise of T is the
kinematics of the Dirac cone: the band velocity has fixed modulus v_F, so a Fermi
disc displaced by A saturates at J = n e v_F, and the differential conductivity is

    sigma_eff/sigma_lin = G(u)/u,  u = A_0/k_F,   G(u)/u -> k_F/A_0  for u >> 1

with G(u) the exact displaced-disc integral (wiki/12 Eqs. 4a.12-4a.13). This is the
CHORD response -- the current the same peak field produces, which is what a field
scan measures. The differential response a weak probe on a strong pump would see is
G'(u), and it falls faster, as 1/(4u^3) (wiki/12 Eq. 4a.16); the two must not be
swapped. At 300 K the continuum values are

| u = A₀/k_F | 0.05 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| σ_eff/σ_lin (continuum) | 1.000 | 0.995 | 0.979 | 0.952 | 0.912 | 0.851 | 0.630 | 0.487 | 0.332 | 0.200 |
| 300², E_F = 0.2 eV (140 partial pts) | 0.930 | 1.062 | 1.015 | 0.978 | — | 0.878 | — | 0.543 | — | 0.249 |
| 147², E_F = 0.2 eV (36) | 6.2 | 1.225 | 1.000 | 0.984 | — | 0.829 | — | 0.532 | — | 0.243 |
| 48², E_F = 0.6 eV (12) | 0.901 | 1.015 | 1.171 | 1.062 | — | 1.003 | — | 0.635 | — | 0.333 |

The continuum curve is **monotone** — no darkening at any field. A Fermi disc holding
one shell of mesh points reproduces the linear limit to ~10 % but develops a spurious
15–20 % bump near u ≈ 0.2–0.5, and at very small u on the coarsest discs the ratio
breaks down entirely (the 147² entry at u = 0.05 is nine points, two of them the
degenerate K pair). That bump is the dip of §7.9: the 48² transmission minimum sits
at 60 kV/cm = u 0.386, exactly the peak of the 48² mesh curve. It shrinks as the disc
fills (300²: 1.06 at worst).

Measured on the runs instead of on the adiabatic sum, the artifact is the spurious
peak Re σ develops before saturation, and it is governed by k_F/Δk alone — not by the
doping (same transient, +100 fs tail, coherent, E_F = 0.6 eV unless noted):

| run | k_F/Δk | σ plateau | σ peak | spurious bump | T plateau | T minimum | dip in T |
|---|---|---|---|---|---|---|---|
| 48², E_F = 0.6 eV | 1.54 | 23.6 | 31.1 | **+32 %** | 0.7275 | 0.6414 | −11.8 % |
| 72², E_F = 0.4 eV | 1.54 | 12.9 | 17.4 | **+35 %** | 0.8487 | 0.7939 | −6.5 % |
| 72², E_F = 0.6 eV | 2.32 | 25.9 | 30.8 | **+19 %** | 0.7081 | 0.6571 | −7.2 % |
| 111², E_F = 0.6 eV | 3.57 | 25.6 | 27.8 | **+8.6 %** | 0.7149 | 0.6900 | −3.5 % |
| 147², E_F = 0.6 eV | 4.73 | 27.0 | 29.0 | **+7.6 %** | 0.6997 | 0.6763 | −3.3 % |

Two different meshes at two different dopings but the same k_F/Δk give the same bump;
filling the disc from 1.5 to 3.6 spacings takes it from a third to under a tenth —
and then it stops, 111² and 147² differing by 78 % in k-points and by nothing in the
dip. The four-mesh figure is `wiki/figures/graphene_T_of_field_mesh.png`; rebuild it
with

```bash
V=<scan root>
python3 field_scan_plot.py \
  --doped "$V/n147/runs/*/graphene_sit_sbe_rt.data" --doped-label '147² ($k_F/\Delta k$ = 4.7)' \
  --series "111² (3.6):$V/n111/runs/*/graphene_sit_sbe_rt.data" \
           "72² (2.3):$V/n72/runs/*/graphene_sit_sbe_rt.data" \
           "48² (1.5):$V/n48/runs/*/graphene_sit_sbe_rt.data" \
  --continuum --t-meas 0.60 0.70 --n-sub 1.65 --out T_of_field_mesh.png
```

The four curves agree above A₀ = k_F to better than 0.02 and disagree below it, which
is what makes the disagreement a discretization error rather than a field scale. How
much of the bump reaches T depends on the sheet impedance, which is why the same
artifact makes a 11.8 % dent at 0.6 eV and 6.5 % at 0.4 eV — the dip depth is not the
artifact, the conductivity bump is. Reproduce the comparison with
`field_scan_plot.py --continuum`, which draws the parameter-free continuum curve
beside the data and prints the residual.

**Most of the darkening before the rise is numerical, but not all of it.** The mesh
series settles the first half: three quarters of the dip disappears between
k_F/Δk = 1.5 and 3.6, and the disappearance is governed by that ratio alone, not by
the doping. It does not settle the second: the dip stops shrinking after that —
−3.5 % at 111² and −3.3 % at 147², with 78 % more k-points between them — so a
residual few-per-cent darkening near A₀ ≈ k_F survives mesh refinement. It is not the
basis either: the same 147² scan at nstate = 4 dips by 3.41 % against 3.39 % at
nstate = 2, unchanged to two digits, though the basis moves T by 0.035 and Re σ by
11 % (§7.13). Both numerical suspects are therefore excluded and the dip is left
unexplained; the remaining candidate is the departure of the real EPM band from an
ideal cone over the excursion A₀ ≈ k_F, which the continuum derivation assumes away. The honest shape is therefore: flat while A₀ ≲ 0.5 k_F, a residual few-per-cent
dip around A₀ ≈ k_F, the drift-saturation rise as 1/E₀, and Landau–Zener darkening at
several hundred kV/cm. An under-resolved Fermi surface multiplies that residual dip by
three or four, which is why Eq. (4a.3) still has to be satisfied.
`drift_saturation.py` produces the table and `graphene_drift_saturation.png`.

**7.13 A doped sheet is not nstate-converged the way an intrinsic one is.** The
pure-gauge restoration (§7.1) makes the *undoped* filling basis-independent to 10⁻⁶ in
T, because it subtracts the adiabatic ground-state current of exactly the same
truncated Hamiltonian. Carriers added on top of that reference get no such
subtraction. Repeating the 147², E_F = 0.6 eV scan at nstate = 4 (dt = 0.1 fs, energy
ledger 5×10⁻⁹ eV/cell, nex baseline identical to nstate = 2 — no leakage):

| E₀ [kV/cm] | u | T (nb = 2) | T (nb = 4) | Re σ/σ_univ (nb = 2) | (nb = 4) |
|---|---|---|---|---|---|
| 1 | 0.01 | 0.70004 | 0.73477 | 26.92 | 23.77 |
| 60 | 0.74 | 0.67629 | 0.70975 | 29.00 | 25.95 |
| 100 | 1.24 | 0.68124 | 0.71716 | 28.69 | 25.46 |

An 11 % change in the sheet response from the basis alone, at 1 kV/cm, where the
vector potential (6×10⁻⁴ a.u.) cannot mix anything — so this is not field-induced.
The natural reading is the second-order repulsion of the π* level by the bands above
it, ΔE ∝ A², which is precisely a shift of the Drude weight: D_spec falls from
0.608 eV to 0.548 eV against the analytic E_F = 0.600.

Two consequences. **(i)** The D → E_F agreement of §7.12's table is a consistency
check, not an independent confirmation: near K the two-band EPM *is* the Dirac model,
so the adiabatic band derivative the restoration leaves behind for the doped carriers
is the exact Dirac velocity and D = E_F follows almost by construction. **(ii)**
Absolute doped conductances from this solver carry an ~10 % basis systematic; ratios,
shapes and field dependences at fixed nstate do not. Notably the dip of §7.12 does
**not**: 3.41 % at nb = 4 against 3.39 % at nb = 2.

## 8. What to look for at production (147², DAST)

- `coh` 1 kV/cm: T ≈ 1 (empty sheet, THz far below the resolvable interband
  window); the pair density after the pulse follows the analytic Γ ∝ E^{3/2}
  table of §1 only from ~30 kV/cm up (mesh-resolved LZ tube).
- `coh` 100 kV/cm: pairs ~10¹² cm⁻² after the pulse; absorption a few per cent
  (creation energy only); `R` small.
- `diss`/`mem` 100 kV/cm: the created carriers are accelerated to ~0.5 eV and
  cooled by optical-phonon emission — the ledger's absorption grows to tens of
  per cent; T_e in `*_sbe_te.data` rises to thousands of K during the cycle and
  decays toward 300 K on the phonon timescale (tens of fs); the Rana channel
  *multiplies* carriers while n < n_i(T_e) and recombines once T_e has dropped.
  `mem` differs from `diss` by the removed dephasing-ionization share (`wiki/10`
  §8.7 logic) and by the T_e-consistent balance.
- Reflection stays ≪ absorption (R ~ (Z₀σ/2)² — a few 10⁻³ even for σ ~ 5 σ_univ).

## 9. Limits recorded (not blockers for this study)

1. **Initial state at T = 0, undoped** (`gs%occup` = integer filling): no thermal
   / doping Drude background, hence no *bleaching* channel; the FD(E_F, T)
   initial occupation is the next increment (it also generalizes the dressed
   reference to fractional filling).
2. Below ~30 kV/cm at 147² the LZ tube is thinner than the mesh spacing — the
   pair creation is then analytic rather than mesh-resolved.
3. Hot phonons are not included (fixed-temperature bath).
4. The Coulomb sector and the T_e fit assume quasi-thermal branch distributions.
5. The sheet field is free-standing; a substrate index enters the boundary
   condition trivially (`wiki/12` Eq. 4) but is not yet a driver option.

## 10. Tests added by this exercise

`test_graphene_dirac_levels.py`, `test_sheet_transmission.py`,
`test_rana_saturation.f90`, `test_colmem_2d.f90`, `test_dirac_te_fit.f90`,
`test_vg_sumrule.f90` (velocity-gauge f-sum rule and the pure-gauge restoration,
§3a), `test_doped_drude.py` (Fermi-Dirac occupation on a k-mesh and its
resolution requirement, Dirac-cone Drude weight and its weak temperature
dependence at fixed density, measured transmission → sheet conductance,
current-saturation field scale, §3b/§7.9) — `python3 tests/run_all.py`, 31/31.

Analysis scripts: `transmission.py` (T/R/A, sheet BC, energy ledger, Re σ),
`saturation_check.py` (populations, Rana ledger, T_e), `sumrule_check.py`
(basis f-sum rule, residual reactive current), `plot_levels.py` (level
populations vs t and k-maps), `drude_check.py` (doped sheet: D, τ, mean free
path, experiment conversion), `field_scan_plot.py` (the T(E₀)/σ(E₀) figure),
`run_field_scan.sh` (the whole scan-and-plot pipeline), `plot_occupation.py`
(the initial level occupation of a doped vs an undoped sheet on the actual mesh:
sheet density, k_F, partially occupied k-points — the pre-flight check for a doped
run), `drift_saturation.py` (the universal drift-saturation curve σ_eff/σ_lin
against A₀/k_F, continuum against any mesh — the diagnostic that separates the
physical brightening from the discretization bump, wiki/12 §4a.5).
