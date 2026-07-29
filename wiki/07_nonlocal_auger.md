# Nonlocal Auger recombination — implementation spec

**Status: MERGED into `develop-2.0.0` (2026-07-02→03, PRs #52/#57/#58/#62).**
The 3D detailed-balance kernel reuse (§0/§5) is done and merged:
`auger_interk_dpop` (the exact time-reverse of `ii_interk_dpop` — same
quadruples, |V(q)|², broadened δ and cited II magnitude; reversed occupations;
hot-gain capped), wired as *one gather, two kernels* inside
`apply_ii_interk_ring`, ring-gated, requiring the II channel (shared constants
— no separate C), with the k-local C·n³ Auger off when the ring is on. Unit
test incl. the **Fermi-Dirac detailed-balance fixed point** (net II+Auger dpop
< 1e-6). The **BGR↔Σ^HF guard** (§0.2b) is in (mutual exclusion, error stop).
The |V(q)|² weight now carries the **CDRB model ε(q) + Cartesian umklapp
G-sum** (§3, PR #57); the **graphene 2D Rana** branch (§6) is implemented &
validated against the cited lifetimes (PR #58); the effective-C(n)
order-of-magnitude check vs S14/L90 is done (§7, PR #62). All four journal
sources (R07/K15/S14/L90) are **source-verified** against the PDFs.
**The one remaining piece:** ~~wire the graphene 2D Rana rate into a live CPTP
SBE channel~~ — **DONE 2026-07-04 (see the UPDATE below)**; graphene
`auger_ok` is now `.true.` (with the `auger_2d_rana` discriminator). GaAs
dynamic λ(n(t)) free-carrier screening — **also DONE 2026-07-04**. The
Bloch-overlap factors I(G) (currently →1) are **resolved as redundant at the
current calibration level** (see the wiki/00 decisions entry): the kernels'
absolute rate scale is pinned to the cited Stobbe/Keldysh magnitude, so the
sub-unity I(G) factors renormalise a quantity that is *calibrated away*; the
residual shape effect on the q-distribution is second order per [L90] (the G-sum and
ε(q) are the big effects — both in), and carrying the plane-wave coefficients
into the SBE dataset would be a large data-plumbing change for it.

**UPDATE 2026-07-04 (wiki/00 TODO sweep): the graphene 2D Rana CPTP channel is
LIVE** (`rana_auger_dpop` in sbe_superres + `apply_rana_auger_ring` in
bloch_solver, ring-gated like graphene e-ph). Net pair relaxation R − G =
(CCCV+CVVV) − (their Eq.-17 reverses) on the INSTANTANEOUS quasi-Fermi levels
(`dirac_mu_2d` inversion of the gathered CB-electron / VB-hole sheet
densities); applied as a uniform-fractional CB→VB (or VB→CB when G > R —
thresholdless carrier multiplication) population transfer with smooth
`cap·(1−exp(−ΔN/cap))` saturation against BOTH the source population and the
destination phase space; trace exactly conserved by construction; coherences
amplitude-damped √(f_new/f_old). ε_r default = 10 (the R07 Fig. 4 benchmark
substrate), override via `sbe_coulomb_epsilon`; T = the e-ph bath
(`sbe_eph_temperature_k`); cell area from V·|b3|/2π (guards: b3 must be the
out-of-plane axis and the b-matrix must be in the GS header). Unit test
`test_rana_auger_cptp` (trace, bounds, the equilibrium μ=0 fixed point, the
implied τ_r = 1.48 ps R07 benchmark, τ→∞ saturation, empty no-op).
**Calc-validated** (graphene 12²×1, serial build): electrons = 2.000 every
step; below the 300 K thermal density the channel gently GENERATES toward
equilibrium (detailed balance — expected, not a bug), above it the post-pulse
nex decays monotonically (−4.6% over the 3 fs tail at n₂d ≈ 6×10¹³ cm⁻²,
consistent with the steep R07 τ_r(n) drop). The GaAs λ²(n(t)) also went live
in the ring kernels: min(Debye, degenerate-TF) crossover on the gathered
excited density (each formula overestimates outside its regime), registry-
gated `dyn_lambda_ok` (GaAs only; Si stays λ=0 per Burt [L90]); one-time
diagnostic print at first activation; verified live on GaAs 4³ impact+ring
(λ² = 1.02e-4 a.u. at n_exc = 6.7e17 cm⁻³ — hand-checked against the Debye
formula; electrons = 8.000 conserved).

**UPDATE 2026-07-03 (branch `claude/auger-eps-q-umklapp`): the K15 ε(q) +
umklapp piece of §3 is IN** — both ring kernels now weight the quadruples with
`interk_vq` = Σ over the 27 neighbouring reciprocal images G of
1/[ε_CDRB(|q+G|)·(|q+G|²+λ²+q²_reg)], with q in the **Cartesian metric** of the
actual (non-orthogonal) cell (bmat rows) and ε(q) the CDRB model dielectric
[K15 Eq. (8); CDRB PRB 47, 9892 (1993)] built from the registry ε∞ and the
valence-gas q_TF/ω_p (n = nelec/V). This replaces the old reduced-coordinate
bare `1/(q²_red + κ²=0.05)` (κ² was an uncited fixed regulariser; the new
q²_reg = (half the smallest grid spacing)² refines with the k-grid). λ²=0
(Si: Burt's dynamical argument [L90]; density-dependent Debye/TF λ(n(t)) for
GaAs = refinement). Both [L90] "must-not-drop" pieces are thus in at the
**overlap-free level (I(G)→1)**: the Bloch-overlap factors I₁₃(G)I₂₄(G′−G)
need the plane-wave coefficients the SBE does not carry — still a refinement.
Same weight in both kernels ⇒ detailed balance preserved exactly (the FD
fixed-point test passes with the linear-regime τ=1e-7; the CDRB weight raised
Γ so the old τ=1e-6 left a visible 1.7e-6 nonlinear tail — not a balance
violation).

**UPDATE 2026-07-03 (branch `claude/graphene-rana-2d`): the graphene 2D Rana
branch (§6) is IMPLEMENTED AND VALIDATED against the cited [R07] lifetime
benchmarks** (the maintainer supplied the journal texts — legal academic
access). `dirac_mu_2d` (quasi-Fermi of the g=4 Dirac gas, [R07 Eq. 7]),
`rana_qtf` ([R07 Eq. 13]; in Hartree a.u. Q_TF = 4kT·log[...]/(ε_r v²)) and
`rana_rcccv` ([R07 Eq. 14] verbatim — **the collinear-collapse
√((k1+Q)(Q−k2)k1k2) is a NUMERATOR factor** (the §6 transcription below had it
inverted; fixed against the paper); M_d/M_e per Eqs. (10)-(12), overlaps = 1
on the collinear line, in a.u. M_d = 2π/(ε_r(Q+Q_TF)); `reverse=.true.` = the
CVCC generation partner, Eq. (17)). `test_rana_2d` asserts the CITED numbers:
τ_r(10¹² cm⁻², 300 K, ε_r=10) = **1.48 ps** (Fig. 4 ~1–2 ps; minority
lifetime 1.1 ps), τ_r(10¹¹) = 7.4 ps (cited > 5 ps), τ_r(10¹², 77 K) =
1.34 ps (cited > 1 ps), plus equilibrium detailed balance (exact), the
CVVV mirror, ε-ordering (Fig. 5) and the [R07] T-trends. The **K15 model
ε(q)** used by the ring kernels (branch `claude/auger-eps-q-umklapp`) was
also verified VERBATIM against the supplied PRB 92, 035207 text (α = 1.563,
same functional form; the Debye-Hückel/Thomas-Fermi λ² formulas match the
Part-G screening module). **Still to wire:** the CPTP SBE channel that
applies R−G = (n−n₀)/τ_r to the graphene cone populations (graphene
`auger_ok` stays `.false.` until then); Remaining otherwise: the overlap factors I(G), per-material C
validation against [S14]/[L90] (§7). Originally recorded 2026-06-30 from the maintainer's
prompt. This is the next major dissipation-physics task: rebuild Auger
recombination from the current **k-local** approximation (initial states pinned
at the band extrema — the first-generation impact-ionization approximation) into
a fully **nonlocal** channel, reusing the machinery of the **already-nonlocal
impact ionization**.

> **Workflow (maintainer, 2026-06-30):** this is a **large SBE edit → do it on
> its own feature branch**, not on `develop-2.0.0`. *"Важно не сломать
> правильный код."* **Validation is by calculation only** — CPTP/number
> conservation AND occupation maps (the carrier "spot" sitting at the
> band-energy-minimum point), not by unit tests alone. The maintainer **approves
> being sent plots** (population maps, lifetime curves). Build `build/salmon`,
> run the small grids below, and send the maintainer the CPTP trace + population
> maps before merging.

## 0. The key idea that makes this tractable

Auger recombination and impact ionization are **time-reverses of each other**
(detailed balance). The squared matrix element |M|² is **identical**; only the
occupation factors differ (§5). Since nonlocal impact ionization already works,
nonlocal Auger is the **same kernel** with the Fermi factors swapped (and the
energy δ-function direction/sign flipped). **Do not duplicate** the k-sampling or
the umklapp G-sum — reuse them from the nonlocal impact-ionization path.

## 0.1 Computation reuse — what is genuinely shared (and what is NOT)

(Maintainer, 2026-06-30: *"эффективно использовать все расчёты — не только
симметрию с impact ionization"*.) Auger should ride the existing stack and add
little but the Fermi-factor reweighting — BUT the sharing must respect the
self-energy decomposition (§0.2), or it double-counts. **What is genuinely shared
(verified against the current code, 2026-06-30):**

1. **The metric-aware *bare* Coulomb building block** `coulomb_kernel(sbe,gs,ik,iq)`
   (bloch_solver):
   $$V_{\rm bare}(\mathbf q)=\frac{\texttt{strength}\cdot 4\pi}{\varepsilon_{\rm bg} \Omega_{\rm cell} N_k (q^2+\kappa_0^2)},$$
   metric-aware **minimum-image** $|\mathbf k-\mathbf q|$ via `gs%b_matrix` (correct
   on the non-orthogonal primitive cells), $q{=}k$ self-term excluded, prefactor
   `sbe%coul_pref`, Yukawa regulariser $\kappa_0^2=$`sbe%coul_screen2`. Today it is
   used **only by Σ^HF** (as the bare exchange). The min-image/metric/self-term
   machinery is what Auger/II reuse — **not** a ready-made screened $\tilde W$.
2. **The ring gather** of the Houston spectrum + populations (`comm_summation` over
   `sbe%icomm`; `sbe%glob_occ`, the gathered `eval,f`) — the SAME one nonlocal-II /
   inter-k e-ph already do once per step. $n$, $E_F$, $k_F$ (hence $\lambda^2$) and
   the occupation factors $f_{n\mathbf k}$ are read from it. No second gather.
3. **The momentum-conserving k-index map** `mp_grid_triple` / `mp_partner_triple`
   (the O(1) `klut` for $\mathbf k_1+\mathbf k_2=\mathbf k_1'+\mathbf k_2'\
   (\mathrm{mod}\ \mathbf G)$ from nonlocal-II) — the umklapp $\mathbf G$-sum is the
   loop over the lattice shifts folding $\mathbf k_4$ into the first BZ.
4. **The GKLS/Lindblad primitive** `amp_damp_channel` — Auger is two such maps
   (recombination $c_1 \to v_1$ + hot promotion $c_1 \to c_{\rm hot}$), exactly
   as the current k-local `apply_auger_recombination`; only the source/target
   $(\mathbf k,\text{band})$ pairs and the **rate** change.

**What is NOT yet built and must be added by THIS task** (do not claim it exists):
- The **screened** interaction $\tilde W(\mathbf q)=V_{\rm bare}(\mathbf q)/\varepsilon(\mathbf q,\omega)$.
  The Part-G dielectric primitives (`eps_lindhard_static`, `eps_thomas_fermi`,
  `tf_kappa2_degenerate`, `debye_kappa2`, `plasmon_freq2`, `lopc_branches`) exist
  and are unit-tested **but are NOT consumed by any channel yet** — wiring them
  into the collision kernel is part of this work.
- A *microscopic* screened e-e: the present `carrier_carrier_relax` is a **flat-rate
  FD-thermalisation** ($\alpha=1-e^{-\nu\tau}$ toward a fitted Fermi–Dirac) — it does
  **not** use $\tilde W$, so there is no screened-e-e kernel to inherit. The screened
  $\tilde W$ Auger needs is **new shared infrastructure** (which a future microscopic
  e-e and the refined II would then also use).

## 0.2 Consistency rules — no double-counting (answers the maintainer's questions)

**(a) Exchange vs correlation — Σ^HF must NOT be screened like the collision channels.**
The self-energy splits (GW-like) into the *coherent* Hartree–Fock exchange
$\Sigma^{\rm HF}= i G v$ (uses the **bare** $v$ = `coulomb_kernel`; a Hermitian
mean-field that renormalises the bands/Rabi, README §8) and the *dissipative*
correlation $\Sigma^{\rm corr}= iG(W-v)$ (uses the **screened** $W=v/\varepsilon$;
this is the e-e/II/Auger scattering). They are **different objects**: the
screened-exchange energy shift lives **solely in Σ^HF**, the scattering uses
$W-v$. So you may reuse the bare `coulomb_kernel` building block, but the **÷ε
screening is applied ONLY to the collision channels (II/Auger/e-e), never to
Σ^HF**. Reusing a single "global screened $\tilde W$" for both would double-count
the screening and smear the exchange shift into the scattering — the standing
**"No HF double-counting"** rule (decisions log: carrier-carrier is the 2nd-Born/GW
correlation, dissipative only; the static screened-exchange shift stays in Σ^HF).

**(b) BGR vs HF — they renormalise the SAME gap; do not apply both.**
Σ^HF already renormalises the band energies dynamically — its diagonal gives
$\tilde\varepsilon^\lambda_k=\varepsilon^\lambda_k-\sum_q V_{k-q} \delta f^\lambda_q$,
so when `yn_sbe_coulomb='y'` the impact-ionization rate is evaluated in the
Σ^HF-renormalized Houston basis and **already sees the carrier-shrunk gap**. The
**BGR** model (Part C7) ADDS a phenomenological $E_{\rm th}(n)=E_{\rm th0}-|K n^{1/3}|$
[Vashishta–Kalia] on top of the threshold. Vashishta–Kalia's $n^{1/3}$ law contains
**both** exchange and correlation; Σ^HF captures **only** exchange. So:
- **Σ^HF on ⇒ keep BGR off** (its exchange part is already in the renormalized
  bands — running both double-counts the exchange gap shrinkage);
- **Σ^HF off ⇒ BGR is the cheap stand-in** for the whole density-dependent gap
  shift (the intended use of `yn_sbe_bgr_threshold`).
The code currently lets both toggle independently → **this is a consistency hazard;
the Auger task should add a guard / documented mutual-exclusion** (and, for Auger,
the gap that enters $E_g$ in the hot-promotion energy must come from ONE source —
the Σ^HF-renormalized bands if Coulomb is on, else the BGR-shifted gap, never both).

**(c) Dynamic screening is material-specific — and the Auger frequency is high.**
The collision screening $\varepsilon(\mathbf q,\omega)$ must dispatch on the material:
**static Lindhard/Thomas–Fermi** (`eps_lindhard_static`, default, all materials) vs
**dynamic LO-phonon–plasmon coupled (LOPC)** which is **GaAs-only** (`lopc_branches`,
`plasmon_freq2`; Part G model c). Crucially, Auger transition energies are **~1 eV
≫ the plasma/phonon frequencies**, so the *static* free-carrier screening is
inappropriate at the Auger energy scale [Burt]:
- **Si:** take $\lambda^2=0$ (no free-carrier screening; §3.6) — screen through the
  background $\varepsilon_{\rm bg}(q)$ (K15 model $\varepsilon(q)$) only;
- **GaAs:** use the dynamic/plasmon-pole form where the LOPC machinery provides it;
- **graphene:** the 2D static↔dynamical RPA switch [R07/Tomadin] (§6).
So "reuse the screening" means **reuse the Part-G *primitives*, dispatched per
material to the right (static vs dynamic) model at the Auger frequency** — not a
single static kernel applied globally. Where no cited dynamic model exists, fall
back to the background $\varepsilon_{\rm bg}(q)$ + $\lambda^2=0$ rather than a wrong
static free-carrier screen.

Current code state to build on:
- Impact ionization: **nonlocal** (`apply_ii_interk_ring` + `ii_interk_dpop` +
  MP map `mp_grid_triple`/`mp_partner_triple`), ring-gated.
- Electron–phonon: matrix element already in the e-ph coupling; inter-k e-ph
  through the ring (`apply_eph_interk_ring`).
- Coulomb / screened exchange: Σ^HF machinery (`compute_coulomb_selfenergy_ring`).
- **Auger: still k-local** (`apply_auger_recombination`, initial states pinned at
  Γ / valley minima) — this is what to replace.

**Gate:** like the inter-k e-ph and nonlocal II, the nonlocal Auger is enabled
**only when the ring is on** (`yn_sbe_superres='y'`); ring-off stays the k-local
behaviour (byte-unchanged).

## 1. Sources (verified — use only these)

| tag | reference | what we take |
|---|---|---|
| **[K15]** | E. Kioupakis, D. Steiauf, P. Rinke, K. T. Delaney, C. G. Van de Walle, *Phys. Rev. B* **92**, 035207 (2015) | full direct + phonon-assisted (indirect) Auger formalism; model ε(q); free-carrier screening λ |
| **[L90]** | D. B. Laks, G. F. Neumark, S. T. Pantelides, *Phys. Rev. B* **42**, 5176 (1990) | exact **nonlocal** pure-Auger in Si; criticality of the umklapp G-sum and q-dependent ε(q); thresholds |
| **[S14]** | D. Steiauf, E. Kioupakis, C. G. Van de Walle, *ACS Photonics* **1**, 643 (2014) | GaAs coefficients (direct/phonon, eeh/hhe) |
| **[R07]** | F. Rana, *arXiv:0705.1204v2* [cond-mat] (2007) | graphene: 2D Auger/impact-ionization on the Dirac spectrum |

### ⚠️ Rejected source — DO NOT USE
**"Haury et al., PRB 57, 11513 (1998), C = 2.0×10⁻³⁰ cm⁶/s" for CdS is a
hallucination.** The real Haury et al. paper is **PRL 79, 511 (1997)** on the
ferromagnetic transition in CdMnTe quantum wells — wrong journal, volume, page,
year, and topic. The coefficient is unconfirmed anywhere. **It has already been
removed from the code** (CdS Auger is gated off). Never re-introduce this number.

## 2. What "nonlocal" means for Auger

- **k-local (current):** initial states 1,2,3 fixed at the band extrema (Γ for
  direct-gap GaAs/CdS; the valley minima for Si). [K15] Eqs. (19)–(20). Valid for
  wide-gap materials where carriers cluster at the band edges.
- **Nonlocal (target):** sum over **all occupied** initial states 1,2,3 across the
  occupied part of the BZ. The final state is fixed by crystal-momentum
  conservation
  $$\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_2 - \mathbf{k}_3 + \mathbf{G}$$
  with **G** an umklapp reciprocal-lattice vector that folds k₄ into the first BZ.
  **Sum over all G**, and evaluate ε(q) at the **actually transferred** momentum
  q = k₁ − k₃ + G.

**Two nonlocal pieces that must NOT be dropped [L90]:**
1. **The umklapp G-sum.** Dropping it underestimates the rate by ~**an order of
   magnitude**.
2. **q-dependent ε(q)** instead of static ε₀. Replacing ε(q)→ε₀ underestimates
   by ~**an order of magnitude** (static ε₀ over-screens large q, but large-q
   umklapp transitions dominate, where true ε(q)→1). Conversely ε=1 over-estimates
   by an order. → both the **full G-sum** and the **full ε(q)** are mandatory.

## 3. Direct (pure-Coulomb) Auger — 3D formalism

[K15] Eqs. (1)–(8); [L90] Eqs. (1)–(8). Composite indices **1**≡(n₁,k₁).

- **Rate (Fermi golden rule):**
  $$R = 2 \frac{2\pi}{\hbar}\sum_{1234} P |M_{1234}|^2 \delta(\epsilon_1+\epsilon_2-\epsilon_3-\epsilon_4)$$
  (the leading 2 = spin; momentum conservation is enforced inside M).
- **Occupation factor:** $P = f_1 f_2 (1-f_3)(1-f_4)$. (eeh ⇒ ∝ n²p; hhe ⇒ ∝ np².
  At n=p: R = Cₙn³ (eeh) and R = Cₚn³ (hhe).)
- **Antisymmetrized |M|²:**
  $$|M_{1234}|^2=|M^d-M^x|^2+|M^d|^2+|M^x|^2,$$
  with $M^d=\langle\psi_1\psi_2|W|\psi_3\psi_4\rangle$, $M^x=\langle\psi_1\psi_2|W|\psi_4\psi_3\rangle$.
- **Screened-Coulomb matrix element (with the umklapp sum):**
  $$\langle\psi_1\psi_2|W|\psi_3\psi_4\rangle=\frac1V\sum_{\mathbf G}\delta_{\mathbf k_1+\mathbf k_2, \mathbf k_3+\mathbf k_4+\mathbf G'} \tilde W(\mathbf k_1-\mathbf k_3+\mathbf G) I_{1,3}(\mathbf G) I_{2,4}(\mathbf G'-\mathbf G)$$
  $$\tilde W(\mathbf q)=\frac{1}{\varepsilon(\mathbf q)}\frac{4\pi e^2}{q^2+\lambda^2},\qquad I_{\alpha\beta}(\mathbf G)=\sum_{\mathbf G_1}c_\alpha^*(\mathbf G_1)c_\beta(\mathbf G_1-\mathbf G)$$
  (the overlap integrals I use the periodic parts u of the Bloch functions — the
  EPM plane-wave coefficients we already have).
- **Model ε(q) [K15] (Cappellini–Del Sole–Reining–Bechstedt):**
  $$\varepsilon(\mathbf q)=1+\Big[(\varepsilon_\infty-1)^{-1}+\alpha(q/q_{\rm TF})^2+\tfrac{\hbar^2q^4}{4m^2\omega_p^2}\Big]^{-1},\quad \alpha=1.563$$
  (limits: ε(0)=ε_∞, ε(∞)→1).
- **Free-carrier screening λ [K15]:** electron part — Debye–Hückel
  $\lambda_e^2=4\pi n e^2/(\varepsilon_\infty k_BT)$ (non-degenerate,
  $E_F-E_{\rm CBM}<\tfrac32k_BT$) or Thomas–Fermi
  $\lambda_e^2=6\pi n e^2/[\varepsilon_\infty(E_F-E_{\rm CBM})]$ (degenerate);
  hole part analogous; total $\lambda^2=\lambda_e^2+\lambda_h^2$.
  **For Si [L90]: take λ=0** (Burt's dynamical argument: Auger transition
  frequencies ~1 eV ≫ plasma ~0.1 eV → static free-carrier screening invalid);
  keep screening only through ε(q).
- **Coefficient:** $C(n)\equiv R(n)/(n^3V)$.

## 4. Phonon-assisted (indirect) Auger — 3D formalism

[K15] Eqs. (9)–(21). **Reuse the existing e-ph matrix elements g**; just remove
the Γ-pinning of the initial states.

- **Rate (2nd order):**
  $$R=2\frac{2\pi}{\hbar}\sum_{1234,\nu\mathbf q}\tilde P |\tilde M_{1234;\nu\mathbf q}|^2 \delta(\epsilon_1+\epsilon_2-\epsilon_3-\epsilon_4\mp\hbar\omega_{\nu\mathbf q})$$
  (upper/lower sign = phonon emission/absorption).
- **Factor:** $\tilde P=f_1f_2(1-f_3)(1-f_4)(n_{\nu\mathbf q}+\tfrac12\pm\tfrac12)$,
  $n_{\nu\mathbf q}$ = Bose. (The emission channel survives at T=0 — needed for
  low-T validation.)
- **8-diagram element:**
  $$|\tilde M|^2=|\tilde M_1{+}\tilde M_2{+}\tilde M_3{+}\tilde M_4-\tilde M_5{-}\tilde M_6{-}\tilde M_7{-}\tilde M_8|^2+|\tilde M_1{+}..{+}\tilde M_4|^2+|\tilde M_5{+}..{+}\tilde M_8|^2$$
  Terms 1–4 = direct-Coulomb M^d with a g insertion on each of the 4 lines (energy
  denominators $\epsilon_m-\epsilon_i\pm\hbar\omega+i\eta$); terms 5–8 = the same
  with exchange M^x. **All eight on the same wavefunction set** (phases matter).
- **e-ph element (already implemented):**
  $g_{n\mathbf k,m\mathbf k+\mathbf q;\nu}=(\hbar/2M_0\omega_{\nu\mathbf q})^{1/2}\langle\psi_{n\mathbf k}|(\partial_{\nu\mathbf q}V)^*|\psi_{m\mathbf k+\mathbf q}\rangle$.
- **local→nonlocal change:** in the Γ-approximation [K15] Eq. (19) the prefactor
  is N³/8; the nonlocal version restores the full k₁,k₂,k₃ sum (with Fermi
  factors) and fixes k₄ via the umklapp relation. Γ-approx error ~
  (k_F/|q|_dominant)² — negligible for nitrides, but it **is** the nonlocality for
  narrow-gap GaAs/Si.

## 5. The Auger ⟷ impact-ionization bridge (reuse the working kernel)

[R07] §2: generation processes are the time-reverses of recombination; same |M|²,
only the occupation factors differ.

| process | occupation factor |
|---|---|
| **Auger (recombination)** | $f_1 f_2 (1-f_3)(1-f_4)$ |
| **Impact ionization (generation)** | $(1-f_1)(1-f_2) f_3 f_4$ |

→ Compute the nonlocal Auger rate with the **same kernel** as the working
nonlocal impact ionization, swapping in the recombination Fermi factors (and
flipping the energy-δ sign/direction). Detailed balance holds at equilibrium
(total generation = total recombination).

## 6. Graphene — a SEPARATE 2D branch (not the 3D machinery)

[R07]. Gapless Dirac spectrum; the 3D umklapp/overlap machinery does not apply —
implement as its own branch.
- Spectrum $E_s(\mathbf k)=s\hbar v|\mathbf k|$, $v=10^8$ cm/s, s=±1; densities
  carry a factor **4** (spin × 2 valleys K,K′).
- **Phase-space restriction (key):** CCCV energy conservation forces k₁,k₂,Q
  **collinear**; overlaps → 1; the process is nearly forbidden (lifetimes > 1 ps).
- Overlaps $|\langle u_{s'k'}|u_{sk}\rangle|^2=\tfrac12[1+ss'\cos\theta]\xrightarrow{\rm collinear}1$.
- Recombination rate $R_{\rm CCCV}$, matrix elements $M_d,M_e$ with Thomas–Fermi
  vector [R07 Eq. (13)] $Q_{\rm TF}=\frac{e^2k_BT}{\pi\varepsilon_\infty\hbar^2v^2}\log[(e^{E_{f+1}/k_BT}+1)(e^{-E_{f-1}/k_BT}+1)]$
  (in Hartree a.u. with $\varepsilon_\infty=\varepsilon_r/4\pi$: $Q_{\rm TF}=4k_BT\log[\cdots]/(\varepsilon_r v^2)$).
- δ-collapsed to a 3D integral [R07 Eq. (14), the paper's main result — **the √ is
  a NUMERATOR factor**; an earlier draft of this line had it inverted, corrected
  2026-07-03 against the journal text, `rana_rcccv`]:
  $$R_{\rm CCCV}=\frac{1}{\hbar^2v} \int_0^\infty \frac{dk_1}{2\pi} \int_0^\infty \frac{dk_2}{2\pi} \int_{k_2}^\infty \frac{dQ}{2\pi} |M(k_1,k_2,Q)|^2 \sqrt{(k_1{+}Q)(Q{-}k_2)k_1k_2} [1{-}f_{-1}(Q{-}k_2)][1{-}f_{+1}(k_1{+}Q)]f_{+1}(k_1)f_{+1}(k_2)$$
  $R=R_{\rm CCCV}(n,p)+R_{\rm CVVV}(n,p)$, $R_{\rm CVVV}(n,p)=R_{\rm CCCV}(p,n)$,
  $1/\tau_r=R/\min(n,p)$. Units: **cm⁻²·s⁻¹** (2D!). **No single C [cm⁶/s]** —
  validate by lifetime curves, not by C. ✅ **IMPLEMENTED & validated**
  (`rana_rcccv`/`dirac_mu_2d`/`rana_qtf`, `test_rana_2d`): τ_r(10¹², 300 K,
  ε_r=10)=1.48 ps (Fig. 4), τ_r(10¹¹)=7.4 ps (>5), τ_r(10¹², 77 K)=1.34 ps
  (>1), G=R at equilibrium. **Still to wire:** the CPTP SBE channel applying
  R−G=(n−n₀)/τ_r to the cone (graphene `auger_ok` stays `.false.` until then).

## 7. Per-material parameters & validation targets (T=300 K unless noted)

> **SOURCE-VERIFIED 2026-07-03** against the journal PDFs supplied by the
> maintainer: the [S14] GaAs C-table below matches the paper verbatim
> (incl. ref-18 direct eeh 5×10⁻³⁴ and the Lush/Strauss experimental
> cross-checks) and the [L90] Si numbers match (Dziewior–Schmid
> Cₙ=2.8×10⁻³¹, Cₚ=0.99×10⁻³¹; and L90's own wording: the neglect of the
> umklapp G-sum and the static-ε₀ replacement are "the two worst
> approximations", each costing ~an order of magnitude — the direct
> motivation of the merged CDRB-ε(q)+umklapp kernel weight, PR #57).
> The §6 graphene targets were re-derived from [R07] directly (PR #58).
>
> ✅ **C(n) ORDER-OF-MAGNITUDE CHECK DONE (2026-07-03, `tests/validate_auger_c.f90`).**
> The effective ring-Auger coefficient C_eff = R/(n²p) was extracted from the
> `auger_interk_dpop` kernel on the real Si primitive **4³** EPM spectrum with
> FD electron/hole populations at n=p, T=300 K (quasi-Fermi levels by
> bisection), in the linear regime (τ=1e-2 a.u.t):
>
> | n=p [cm⁻³] | C_eff (σ=0.01 Ha) | C_eff (σ=0.02 Ha) |
> |---|---|---|
> | 10¹⁸ | ~0 (grid-unresolved: <1e-4 e⁻/cell) | ~0 |
> | 10¹⁹ | 2.6×10⁻³⁰ | 1.2×10⁻²⁹ |
> | 10²⁰ | 2.8×10⁻³⁰ | 1.2×10⁻²⁹ |
>
> **Findings:** (1) C_eff is **n-independent** across 10¹⁹→10²⁰ (2.6 vs 2.8e-30),
> confirming the **R ∝ n³** pair-density scaling the extraction assumes — the
> single most important structural check. (2) The magnitude lands **within
> ~1 order of magnitude** of the source-verified Dziewior–Schmid Si Cₙ=2.8×10⁻³¹
> (ratio ~9–10× at the tighter σ=0.01 Ha) — exactly the order-of-magnitude
> agreement [L90]/[S14] describe for this class of calculation. (3) The ×4
> σ-sensitivity (0.01 vs 0.02 Ha) is the energy-conservation broadening on the
> coarse 4³ mesh; a finer grid + smaller σ converges downward (σ=0.01 already
> halves it). **Honest caveat:** the magnitude is tied to the cited Keldysh II
> fit (P=2×10¹² s⁻¹ eV⁻², a=2) by detailed balance, NOT a first-principles
> pure-Auger C — this validates that the II↔Auger detailed-balance construction
> yields a **physically reasonable** recombination coefficient, and that its C is
> the right order of magnitude, not a first-principles reproduction of the exp C.
> Run: `gfortran -O2 src/ssbe/sbe_superres_ssbe.f90 tests/validate_auger_c.f90 -o v`
> then `./v` in a dir with `Si_prim_eigen.data`/`Si_prim_k.data` (4³, nstate=20).

Keep the gap a tunable parameter (rigid scissor of the conduction band) for
alloys / DFT-gap correction.

- **GaAs — ADD [S14].** Direct-gap, Eg=1.43 eV, n=p=10¹⁸ cm⁻³. No spin-orbit;
  exclude the split-off band (350 meV below VBM) as initial and (by central-zone
  energy conservation) final hole state.

  | process | direct C [cm⁶/s] | phonon C [cm⁶/s] |
  |---|---|---|
  | eeh | < 10⁻³³ (negligible) | 1.1×10⁻³¹ |
  | hhe | 2.2×10⁻³¹ | 3.1×10⁻³¹ |

  Checks: hhe ≈ 5× eeh at Eg=1.43; direct eeh negligible; direct hhe ~ phonon hhe.
  Experiment cross-check: (7±4)×10⁻³⁰ (Strauss); ≤1.6×10⁻²⁹ (Lush) — order-of-mag
  agreement. Grids: 80³ (direct), 56³ (phonon) — but here use the small grids in §8.
- **Si — ADD [L90].** Indirect, CBM at k≈0.85·(2π/a) along Δ, 6 valleys. **λ=0**
  (§3); ε(q) mandatory; umklapp G-sum mandatory. Fermi–Dirac for majority,
  Boltzmann for minority. For hhe include split-off (full 27 hh/lh/so transitions).

  | process | C (exp, Dziewior–Schmid) [cm⁶/s] | pure-Auger expectation |
  |---|---|---|
  | eeh (n-Si) | Cₙ=2.8×10⁻³¹ | pure-AR reproduces → agreement |
  | hhe (p-Si) | Cₚ=0.99×10⁻³¹ | pure-AR ~10× low → remainder is the phonon channel |

  Geometry: same-valley (×6) and orthogonal-valley (×12) contribute at the same
  order; opposite-valley (×3) ~2 orders smaller — droppable. Thresholds: eeh≈8 meV,
  hhe≈76 meV.
- **graphene — UPDATE (2D branch §6) [R07].** No single C. ε_∞=10ε₀ (Al₂O₃) or 4ε₀
  (SiO₂); v=10⁸ cm/s. Lifetime targets: τ_r>1 ps at n,p<10¹² cm⁻² (all T); τ_r>5 ps
  at n,p<10¹¹ cm⁻²; τ_r≈1.1 ps at n=10¹² cm⁻², T=300 K, ε_∞=10ε₀. Qualitative:
  generation/recombination curves cross at the equilibrium densities; at n=p,
  R_CCCV=R_CVVV; smaller ε_∞ ⇒ larger rates.
- **CdS — DISABLE.** The only coefficient ever cited (Haury PRB 57, 11513) is a
  hallucination (§1) — no verified value. Physically CdS is wide-gap (Eg≈2.42 eV)
  ⇒ direct Auger exponentially suppressed; the dominant channel would be the
  indirect (phonon/alloy/defect) one. **Already disabled in the code.** If revived
  later: use the nitride scheme [K15] (direct negligible, indirect dominant) and
  **only a verified** coefficient — never the Haury number.

## 8. Numerics (general)

- **δ-functions** → finite-width Gaussians; converge the width and the η in the §4
  denominators separately. [K15] guides: δ≈0.1 eV (direct), 0.3 eV (phonon); result
  stable for η in 10–500 meV (coefficients vary ≤~20–40 %, intermediate states are
  off-resonance).
- **k-grids: SMALL only** — odd **7×7×7 for GaAs**, even **8×8×8 for Si** (heavy).
- **k-point cutoff [K15]:** keep k within $E_{\rm cutoff}=E_F+M k_BT$ of the band
  edge; pick integer M so the post-cutoff carrier density differs from the full one
  by <1 %.
- **Statistics:** Fermi–Dirac at high density (degeneracy). C becomes
  density-dependent (phase-space filling) — the power-law exponent drops below 3.
  Not a bug: report C(n)=R(n)/(n³V), not a constant.

## 9. Acceptance checklist

- [ ] **GaAs:** direct+phonon, eeh/hhe reproduce the §7 table (hhe/eeh≈5; direct
      eeh<10⁻³³).
- [ ] **Si:** pure-AR eeh≈Cₙ; pure-AR hhe≈×10 below Cₚ; thresholds 8/76 meV; G-sum
      and ε(q) on; λ=0.
- [ ] **graphene:** 2D branch; lifetimes meet §7; R_CCCV=R_CVVV at n=p.
- [ ] **CdS:** disabled; the Haury number is nowhere in the code. *(done — this
      session)*
- [ ] Auger reuses the nonlocal impact-ionization kernel with swapped Fermi factors
      (§5); k-sampling and G-sum are **not** duplicated.
- [ ] **Ring-gated:** nonlocal Auger active only with `yn_sbe_superres='y'`;
      ring-off behaviour byte-unchanged.
- [ ] **Validation by calculation:** CPTP/number conservation + population maps
      (carrier spot at the band-minimum) sent to the maintainer; large edit on its
      own feature branch.
