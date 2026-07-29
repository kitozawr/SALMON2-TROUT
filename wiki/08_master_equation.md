# The master equation — complete mathematical specification

**Purpose of this page.** One self-contained mathematical statement of *everything*
the solver propagates, so that no textbook excavation is needed to see how each
physical effect enters. Every term below is written out in full, tied to the exact
namelist flag that switches it on, and to the routine that implements it. Long
formulas are the point. Companion pages: per-channel modelling assumptions in
[`01_physics_models.md`](01_physics_models.md), every constant with its citation in
[`02_constants.md`](02_constants.md), the integrator in
[`03_numerical_methods.md`](03_numerical_methods.md), the nonlocal Auger/II theory in
[`07_nonlocal_auger.md`](07_nonlocal_auger.md).

Units: Hartree atomic units throughout ($\hbar = m_e = |e| = 1$); user-facing
namelist knobs are in eV/fs/K/cm⁻³ and converted internally.

---

## 0. The object and the full equation

The dynamical object is the **one-particle reduced density matrix in the fixed
(field-free) Bloch band basis**, one matrix per crystal-momentum grid point:

$$
\rho_{nm}(\mathbf k, t) = \langle \hat c^\dagger_{m\mathbf k} \hat c_{n\mathbf k}\rangle,
\qquad n,m = 1\ldots N_b,\quad \mathbf k \in \text{MP grid},
$$

with the equilibrium initial condition $\rho_{nm}(\mathbf k,0)=\delta_{nm} f^0_n$,
$f^0_n \in \{0, \texttt{occ}\}$, where $\texttt{occ} = 2$ (scalar) or $1$
(spinor, `yn_sbe_spinor='y'`). Diagonal elements are band populations,
off-diagonals are interband coherences (microscopic polarizations).

The **full master equation** is a Gorini–Kossakowski–Sudarshan–Lindblad (GKLS)
quantum master equation, k-point by k-point, with the nonlocal channels coupling
the k-points through gathered quantities:

$$
\boxed{ 
\partial_t \rho(\mathbf k)  = 
\underbrace{- i\big[ H_{VG}(\mathbf k,t)  +  \Sigma^{\rm HF}[\rho](\mathbf k,t) , \rho(\mathbf k)\big]}_{\text{coherent (unitary) part, §1–§2}} +
 \sum_{c} \mathcal D_c[\rho](\mathbf k)
 }
$$

where the sum runs over the enabled CPTP dissipation channels $c$:

| $\mathcal D_c$ | physics | flag(s) | section | routine |
|---|---|---|---|---|
| $\mathcal D_{\rm KZ}$ | Kuhn–Zurek/Caldeira–Leggett wave-packet dephasing | `sbe_decoh_temperature_k`, `sbe_decoh_tau_m_fs` (>0) | §4.1 | Hadamard kernel in `houston_dissipate` |
| $\mathcal D_{\rm e\text{-}ph}$ | electron–phonon population relaxation | `yn_sbe_eph` (+ ring: `yn_sbe_superres`) | §4.2 | `apply_eph_relaxation` / `eph_interk_dpop` |
| $\mathcal D_{\rm II}$ | impact ionization (carrier multiplication) | `yn_sbe_impact_ionization` (+ ring) | §4.3 | `apply_impact_ionization` / `ii_interk_dpop` |
| $\mathcal D_{\rm Aug}$ | Auger recombination (3 representations) | `yn_sbe_auger` (+ ring) | §4.4 | `auger_interk_dpop` / `apply_auger_recombination` / `rana_auger_dpop` |
| $\mathcal D_{\rm ee}$ | carrier–carrier thermalization | `yn_sbe_eeh` | §4.5 | `carrier_carrier_relax` |

Two more physical effects are **not** dissipators and enter elsewhere:
$\Sigma^{\rm HF}$ (Coulomb exchange, `yn_sbe_coulomb`) is a *Hermitian* mean field
inside the commutator (§2); the band-gap renormalization
(`yn_sbe_bgr_threshold`) is a *parameter modifier* of $\mathcal D_{\rm II}$
(§4.3.4) and is **mutually exclusive** with $\Sigma^{\rm HF}$ (§2.3).

Trace convention: $N_e = \sum_{\mathbf k} w_{\mathbf k} \mathrm{Tr} \rho(\mathbf k)$
with MP weights $w_{\mathbf k}=1/N_k$ is conserved **exactly** by every term above
(the unitary part trivially; each $\mathcal D_c$ by construction — the per-channel
invariants are tabulated in §6).

---

## 1. The coherent part: velocity gauge, $H_{VG}$, and the Houston basis

### 1.1 From minimal coupling to $H_{VG}$

Minimal coupling with a spatially uniform vector potential $\mathbf A(t)$ (dipole
approximation, length scale of the field ≫ cell):

$$
\hat H(t) = \tfrac12\big(\hat{\mathbf p} + \mathbf A(t)\big)^2 + \hat V
= \hat H_0 + \mathbf A(t)\cdot\hat{\mathbf p} + \tfrac12 A^2(t).
$$

Projected onto the field-free Bloch eigenbasis $\{|n\mathbf k\rangle\}$ of
$\hat H_0$ ($\hat H_0|n\mathbf k\rangle = \varepsilon_n(\mathbf k)|n\mathbf k\rangle$),
and with the crystal momentum **fixed on the grid** (this is the defining property
of the velocity gauge — no moving grid), the matrix propagated is

$$
\big(H_{VG}\big)_{nm}(\mathbf k,t) = \varepsilon_n(\mathbf k) \delta_{nm} +
 \mathbf A(t)\cdot \boldsymbol\pi_{nm}(\mathbf k),
$$

where the scalar $\tfrac12A^2(t)$ is a global (k- and band-independent) phase and
is dropped from the propagation — but **restored exactly** wherever an absolute
kinetic energy is compared to a threshold (impact ionization §4.3, e-ph saturation
§4.2), because by the Houston identity
$\varepsilon^{\rm kin} = E_a(\mathbf k,t) + \tfrac12A^2 - E_{\rm CBM}$ is the
field-free kinetic energy of the accelerated carrier,
$E_n(\mathbf k + \mathbf A) \approx E_a$.

### 1.2 The velocity operator $\boldsymbol\pi$

$$
\boldsymbol\pi_{nm}(\mathbf k) = \mathbf p_{nm}(\mathbf k) + \big(\hat v_{\rm NL}\big)_{nm}(\mathbf k),
\qquad \hat v_{\rm NL} = - i [\hat{\mathbf r}, \hat V_{\rm NL}].
$$

For a **local** EPM pseudopotential $\hat v_{\rm NL}\equiv 0$ (the `rvnl_tm` block
of `_tm.data` is zero). Two genuinely nonlocal cases carry it:

* **spin-orbit (spinor GaAs)**: $\hat v_{SO} = \nabla_{\mathbf k}\hat H_{SO}(\mathbf k)$,
  computed analytically and written into block 2; the spinor problem doubles the
  basis, $H^{\rm spinor}_0 = H^{\rm loc}\otimes\mathbb 1_2 + H_{SO}$, and the Bloch
  equation stays a **single** $2N_b\times2N_b$ equation (spin channels do not
  factorize). Requires `yn_sbe_spinor='y'` + `yn_vnl_correction='y'`.
* **real DFT input** (`theory='dft'` + `yn_out_tm='y'`): the Kleinman–Bylander
  nonlocal projector velocity, genuinely nonzero.

Everything downstream — $H_{VG}$, the Houston diagonalizations, the branch
velocities of §4.1, the current of §5 — uses the same $\boldsymbol\pi$
consistently.

### 1.3 The Houston (adiabatic) basis — where all dissipators act

At each step and each $\mathbf k$ the instantaneous eigenproblem

$$\big(H_{VG}(\mathbf k,t) + \Sigma^{\rm HF}(\mathbf k,t)\big) W(\mathbf k,t) = W(\mathbf k,t) \mathrm{diag}\{E_a(\mathbf k,t)\}$$

defines the Houston (adiabatic) basis $W$ and quasi-energies $E_a$. The
**physical** (real-carrier) populations are

$$f_a(\mathbf k,t) = \big(W^\dagger \rho  W\big)_{aa},$$

while the fixed-basis diagonals $\rho_{nn}$ contain *virtual* (field-dressed)
admixtures during the pulse. **Every dissipator is applied in the Houston
basis**: rotate $\tilde\rho = W^\dagger\rho W$, act on $\tilde\rho$, rotate back.
This is the gauge-consistent choice (Wismer–Yakovlev PRB 97, 144302): the
collision partners and Pauli factors are evaluated on *real* carriers, and when
$\Sigma^{\rm HF}$ is on, the Houston basis is automatically the
*Coulomb-renormalized* adiabatic basis.

The interband transition (Zener tunnelling / multiphoton injection) is contained
**entirely in the coherent part**: sweeping $\mathbf k + \mathbf A(t)$ through an
avoided crossing transfers population between Houston branches — no injection
term is (or may be) added by hand.

---

## 2. Coulomb exchange $\Sigma^{\rm HF}$ (`yn_sbe_coulomb`)

### 2.1 The mean field

Time-dependent Hartree–Fock exchange (Golde–Kira–Meier–Koch, Phys. Status Solidi
B 248, 863 (2011), Eqs. 4–5), evaluated on the **deviation from equilibrium**:

$$
\Sigma^{\rm HF}_{nm}(\mathbf k)
= -\sum_{\mathbf q\neq \mathbf k} V(\mathbf k-\mathbf q) \delta\rho_{nm}(\mathbf q),
\qquad \delta\rho = \rho-\rho_0,\quad \rho_0 = \mathrm{diag}(f^0),
$$

$$
V(\mathbf p) = \frac{\texttt{sbe\_coulomb\_strength}\cdot 4\pi}
{\varepsilon \Omega_{\rm cell} N_k \big(|\mathbf p|^2+\kappa^2\big)},
$$

with $\varepsilon$ the material dielectric (registry; override
`sbe_coulomb_epsilon`) and $\kappa$ = `sbe_coulomb_screen_au`. The $\delta\rho$
subtraction makes $\Sigma^{\rm HF}(t{=}0)=0$ exactly: the EPM/DFT bands already
contain the equilibrium exchange, so only the **carrier-induced** renormalization
is added.

### 2.2 What the single commutator contains

Substituting into $-i[H_{VG}+\Sigma^{\rm HF},\rho]$ and separating diagonal /
off-diagonal parts reproduces **both** textbook renormalizations at once —
no separate terms are needed:

* diagonal → renormalized single-particle energies (dynamic band-gap shrinkage):
  $\tilde\varepsilon^\lambda_{\mathbf k} = \varepsilon^\lambda_{\mathbf k} - \sum_{\mathbf q} V_{\mathbf k-\mathbf q}  f^\lambda_{\mathbf q};$
* off-diagonal → renormalized Rabi frequency (excitonic enhancement of the drive):
  $\Omega_{\mathbf k} = \mathbf d_{\mathbf k} \cdot \mathbf E + \sum_{\mathbf q} V_{\mathbf k-\mathbf q}  p_{\mathbf q},$
  with the Pauli factor $(1-f^e-f^h)$ emerging automatically from the commutator
  structure.

$\Sigma^{\rm HF}$ is Hermitian ⇒ it is a *unitary* generator: trace and
positivity are untouched (CPTP-safe). It is frozen at $\rho(t)$ over each step
(mean-field predictor) and folded into **both** the CF4 Hamiltonians and the
Houston diagonalization of every dissipator. On folded supercells the spurious
inter-coset exchange is removed by the sublattice projection
(`yn_sbe_hf_sublattice_proj`, see [`05_folding_unfolding.md`](05_folding_unfolding.md)).

### 2.3 BGR: the stand-in, never the companion

`yn_sbe_bgr_threshold='y'` lowers the impact-ionization threshold with density,
$E_{\rm th}(t) = E_{\rm th,0} - |K n^{1/3}(t)|$ (Vashishta–Kalia form), as a
*cheap stand-in* for exactly the diagonal shift of §2.2 when Coulomb is off.
Enabling both **double-counts the same physics** — the init aborts
(`error stop`) by design. Maximally-loaded runs use $\Sigma^{\rm HF}$ and keep
BGR off (see `samples/exercise_x11_full_dissipation_showcase/`).

---

## 3. The dissipative superstructure: GKLS form, CPTP maps, splitting

### 3.1 GKLS generators and the two primitives

Every channel is a Gorini–Kossakowski–Sudarshan–Lindblad generator

$$\mathcal D[\rho] = \sum_j \Big( L_j \rho L_j^\dagger - \tfrac12\{L_j^\dagger L_j, \rho\} \Big),$$

whose finite-step map $e^{\tau\mathcal D}$ is **completely positive and
trace-preserving for every $\tau\ge0$** — the defining correctness property of
this fork (no positivity clipping anywhere). Concretely only two primitive maps
are ever used, both in the Houston basis:

**(a) Amplitude damping** $a \to b$ with step probability
$\gamma = 1-e^{-\Gamma\tau}$: the Kraus map
$K_0=\mathrm{diag}(\ldots,\sqrt{1-\gamma},\ldots)$, $K_1=\sqrt{\gamma} |b\rangle\langle a|$
gives exactly

$$
\tilde\rho_{aa}\to(1-\gamma)\tilde\rho_{aa},\qquad
\tilde\rho_{bb}\to\tilde\rho_{bb}+\gamma\tilde\rho_{aa},\qquad
\tilde\rho_{ab'}\to\sqrt{1-\gamma} \tilde\rho_{ab'}\ (b'\neq a).
$$

This is why every population-moving channel damps the coherences of a level that
lost population by $\sqrt{f^{\rm new}_a/f^{\rm old}_a}$ — that factor **is**
$\sqrt{1-\gamma}$, i.e. the *exact* amplitude-damping Kraus coherence factor,
not a phenomenological add-on. Populations always move with the saturating form
$f (1-e^{-\Gamma\tau})\le f$ (never negative) into Pauli-blocked destinations
clamped to $[0,\texttt{occ}]$.

**(b) Pure dephasing** by a Hadamard (entrywise) product with a positive
semidefinite kernel matrix $M$: $\tilde\rho \to M \circ \tilde\rho$. By the Schur
product theorem this is CPTP whenever $M\succeq0$ with unit diagonal (§4.1).

### 3.2 Operator splitting (why the composition is still CPTP)

The full step over $h$ (see [`03_numerical_methods.md`](03_numerical_methods.md)
for the integrator itself) is the Strang sandwich

$$
\rho(t+h) = D \left(\tfrac h2\right)\circ
\Big[S_2(p_1h)\circ S_2(p_2h)\circ S_2(p_1h)\Big]\circ
D \left(\tfrac h2\right)[\rho(t)],
$$

with $S_2$ the CF4/Magnus unitary sub-step on Gauss–Legendre nodes (exact matrix
exponentials by eigendecomposition — unitary to machine precision) and $D$ the
composition of all enabled dissipators. The Suzuki–Yoshida triple-jump
coefficients wrap **only** the unitary part: the middle sub-step $p_2h<0$ is a
harmless backward *unitary* rotation, but a negative-time dissipator would flip
the Gaussian kernel into $e^{+\lambda\Delta X^2|\tau|}$ — not positive
semidefinite (violates the Schoenberg/Bochner criterion) — and break complete
positivity. A composition of CPTP maps is CPTP, so the whole step is.

The nonlocal (ring) channels — inter-k e-ph, nonlocal II/Auger, 2D Rana — are
applied **once per step after the unitary block** (they are MPI-collective:
one gather of $\{E_a,f_a\}$ over all $\mathbf k$, then a pure kernel, then the
local application), which is the same Strang-consistent first-order placement as
$D(h/2)$ with the identical CPTP guarantees per map.

---

## 4. The channels, term by term

### 4.1 Kuhn–Zurek / Caldeira–Leggett dephasing $\mathcal D_{\rm KZ}$ *(optional; off in collision-loaded runs)*

Wave-packet decoherence by phonon-bath position monitoring. Each Houston branch
$a$ carries a wave-packet position $X_a(t)$ advanced by its group velocity,

$$
V_a(\mathbf k,t) = \big(W^\dagger\boldsymbol\pi W\big)_{aa} + \mathbf A(t)\cdot\hat{\mathbf e},\qquad
X_a \mathrel{+}= \tfrac12\big(V_a(t)+V_a(t+h)\big)h,
$$

and the coherences decay with the **squared branch separation**:

$$
\tilde\rho_{ab}  \longleftarrow  \exp \big[-\lambda (X_a-X_b)^2 \tau\big] \tilde\rho_{ab},
\qquad \lambda = \frac{k_B T}{\tau_m}
$$

(`sbe_decoh_temperature_k`, `sbe_decoh_tau_m_fs`). The kernel matrix
$M_{ab}=e^{-\lambda(X_a-X_b)^2\tau}$ is a Gaussian RBF kernel — positive
semidefinite by Schoenberg's theorem — so the Hadamard map is exactly CPTP for
any $\tau\ge0$ and any parameters (§3.1b). This is the Caldeira–Leggett
high-temperature position-coupling limit: $\Gamma_{ab} \propto k_BT \Delta X^2/\tau_m$
= Zurek's decoherence rate for spatially separated packets.

**When not to use it:** with collision channels on, they already amplitude-damp
the coherences (§3.1a) — adding KZ double-counts decoherence. For graphene it is
forbidden outright (`error stop`): gapless Dirac coherence loss is many-body.

### 4.2 Electron–phonon $\mathcal D_{\rm e\text{-}ph}$ (`yn_sbe_eph`)

**Rate scale.** A saturating collision rate versus carrier kinetic energy
$\varepsilon$ (measured from the nearest band edge, $A^2/2$ restored):

$$\nu(\varepsilon) = \nu_{\rm sat}\Big[1 - e^{-(\varepsilon/\varepsilon_0)^{n}}\Big]$$

(smooth by design — a hard min() derivative kink destabilizes the stiff solver;
$\nu_{\rm sat}$, $\varepsilon_0$, $n$ = `sbe_eph_nu_sat` (registry default),
`sbe_eph_eps0_ev`, `sbe_eph_n`; cited per material in `02_constants.md`).

**Phonon table.** Each material carries its cited mode list
$\{\hbar\omega_p, w_p\}$: GaAs = polar-LO + 5 intervalley deformation modes; Si =
6 intervalley g/f modes; CdS = single Fröhlich LO (38 meV); graphene = E2g(Γ,
196 meV) + A1′(K, 160 meV, ×2 GW). Relative weights $w_p \propto D_p^2/\hbar\omega_p$
(or $\langle g^2\rangle/\hbar\omega$), normalized. Thermal split per mode with the
Bose factor $N_B(\hbar\omega_p/k_BT_{\rm ph})$:

$$
f^{\rm em}_p = \frac{N_B+1}{2N_B+1},\qquad f^{\rm ab}_p = \frac{N_B}{2N_B+1},
\qquad \frac{f^{\rm em}}{f^{\rm ab}} = \frac{N_B+1}{N_B}\ \text{(detailed balance)}.
$$

**Inter-k (ring) form** — the physical one on primitive cells (intervalley final
states live at *different* $\mathbf k$). For each source Houston state $(a,\mathbf k)$
the partial rate to every destination $(b,\mathbf q)$ over the **whole BZ**:

$$
\Gamma_{(a\mathbf k)\to(b\mathbf q)} =
\nu(\varepsilon_{a\mathbf k})\sum_p w_p
\Big[ f^{\rm em}_p \delta_\sigma \big(E_a(\mathbf k)-E_b(\mathbf q)-\hbar\omega_p\big) +
    f^{\rm ab}_p \delta_\sigma \big(E_b(\mathbf q)-E_a(\mathbf k)-\hbar\omega_p\big)\Big]
\Big[1-\frac{f_b(\mathbf q)}{\texttt{occ}}\Big]_+ ,
$$

with $\delta_\sigma$ a normalized Gaussian of width $\sigma$
(`sbe_search_sigma_e_ev`, grid-matched default 0.2 eV). The total out-rate
$\Gamma_{\rm out}=\sum_{b\mathbf q}\Gamma_{(a\mathbf k)\to(b\mathbf q)}$ moves

$$
\Delta f_{a\mathbf k} = - f_{a\mathbf k}\big(1-e^{-\Gamma_{\rm out}\tau}\big),\qquad
\Delta f_{b\mathbf q} = + \big|\Delta f_{a\mathbf k}\big| \frac{\Gamma_{(a\mathbf k)\to(b\mathbf q)}}{\Gamma_{\rm out}},
$$

so $\sum\Delta f = 0$ **identically** (trace exact), and the source's coherences
damp by $e^{-\frac12(\Gamma_{ {\rm out},a}+\Gamma_{ {\rm out},b})\tau}$
(§3.1a). Implemented in `eph_interk_dpop` (pure, unit-tested), gathered/applied
once per step through the ring. The k-local variant (`.not.` ring) restricts the
sum to $\mathbf q=\mathbf k$ — correct for folded supercells (valleys fold onto
same-k bands) and for intra-valley polar-optical (CdS); it is gated off when the
ring is on (no double count).

### 4.3 Impact ionization $\mathcal D_{\rm II}$ (`yn_sbe_impact_ionization`)

#### 4.3.1 The cited magnitude (both representations share it)

$$
\gamma(\varepsilon^{\rm kin}) = P \big(\varepsilon^{\rm kin} - E_{\rm th}\big)^{a} 
\Theta_{\rm ramp}\big(\varepsilon^{\rm kin}-E_{\rm th}\big),\qquad
\varepsilon^{\rm kin} = E_h(\mathbf k,t)+\tfrac12A^2 - E_{\rm CBM},
$$

Stobbe–Redmer–Schattke quartic for GaAs ($a=4$, $P=2\times10^{12} \mathrm{s^{-1}eV^{-4}}$,
$E_{\rm th}=2.1$ eV), Keldysh quadratic for Si ($a=2$, $E_{\rm th}=1.1$ eV) and CdS
($E_{\rm th}=3.6$ eV cited; $P$ = **fit parameter**, must be given explicitly).
$\Theta_{\rm ramp}$ = the step smoothed over `sbe_ii_ramp_ev` (the fit's energy
resolution).

#### 4.3.2 Momentum-conserving nonlocal (ring) form — the physical one

The true two-particle event: hot conduction electron $(\mathbf k_1, h)$ +
valence electron $(\mathbf k_2, v)$ → two conduction electrons
$(\mathbf k_1', c)$, $(\mathbf k_2', c)$ + a hole, with **exact crystal-momentum
conservation** on the MP grid,

$$
\mathbf k_2' = \mathbf k_1 + \mathbf k_2 - \mathbf k_1' \ (\mathrm{mod}\ \mathbf G)
\quad\Longleftrightarrow\quad
m_2' = \big(m_1 + m_2 - m_1'\big)\bmod n \ \text{per axis (integer index map)}.
$$

Partial rate of the quadruple:

$$
\Gamma_{\mathbf k_1 h}^{(\mathbf k_1' \mathbf k_2)} =
\gamma(\varepsilon^{\rm kin}_{h\mathbf k_1}) 
\big|V(\mathbf k_1-\mathbf k_1')\big|^2 
\delta_\sigma \big(E_h(\mathbf k_1)+E_v(\mathbf k_2)-E_c(\mathbf k_1')-E_c(\mathbf k_2')\big) 
\frac{f_{v\mathbf k_2}}{\texttt{occ}}
\Big[1-\frac{f_{c\mathbf k_1'}}{\texttt{occ}}\Big]_+
\Big[1-\frac{f_{c\mathbf k_2'}}{\texttt{occ}}\Big]_+ .
$$

The screened Coulomb weight carries the **umklapp G-sum in the Cartesian metric
of the actual (possibly non-orthogonal) cell** and the **CDRB model dielectric**
[K15 Eq. (8); CDRB PRB 47, 9892 (1993)]:

$$
\big|V(\mathbf q)\big|^2  \to  \sum_{\mathbf G \in \{-1,0,1\}^3}
\frac{1}{\varepsilon_{\rm CDRB}\big(|\mathbf q+\mathbf G|^2\big) 
\big(|\mathbf q+\mathbf G|^2 + \lambda^2 + q^2_{\rm reg}\big)},
$$

$$
\varepsilon_{\rm CDRB}(q^2) = 1 + \left[\frac{1}{\varepsilon_\infty - 1} +
\alpha \frac{q^2}{q_{\rm TF}^2} + \frac{q^4}{4 \omega_p^2}\right]^{-1},
\qquad \alpha = 1.563,
$$

built from the **valence** gas ($n_v = n_{\rm elec}/\Omega_{\rm cell}$,
$k_F=(3\pi^2 n_v)^{1/3}$, $q_{\rm TF}^2 = 4k_F/\pi$, $\omega_p^2=4\pi n_v$);
$q^2_{\rm reg}$ = (half the smallest k-grid spacing)² regularizes the discrete
$q\to0$ term and refines with the grid. Dropping either the G-sum or
$\varepsilon(q)$ underestimates the rate by ~10× [L90]. The Bloch overlap factors
are kept at $I(\mathbf G)\to1$ — *deliberately*: the absolute magnitude is pinned
to the cited $\gamma(\varepsilon)$ anyway, so sub-unity overlaps renormalize a
calibrated-away quantity (full argument in `00`/`07`).

Each realized event applies the exactly-trace-conserving stencil

$$
\Delta f_{h\mathbf k_1} = -w,\quad \Delta f_{c\mathbf k_1'} = +w,\quad
\Delta f_{v\mathbf k_2} = -w,\quad \Delta f_{c\mathbf k_2'} = +w,
$$

with the total out capped at $f_{h\mathbf k_1}(1-e^{-\Gamma_{\rm tot}\tau})$ and
distributed $\propto$ partial rates. **One electron promoted per event, one hole
left behind: carrier multiplication with the pair landing in the *correct*
(momentum-resolved) valleys** — on Si this resolves the indirect gap explicitly.
Cost O($N_k^3$), exactly additive over the source $\mathbf k_1$ ⇒ MPI-distributed
over each rank's k-range and `comm_summation`-ed (O($N_k^3/P$)).

#### 4.3.3 Dynamic free-carrier screen $\lambda^2(n(t))$ (GaAs only)

$$
\lambda^2(t) = \min\Big[\underbrace{\tfrac{4\pi  n_{\rm exc}(t)}{\varepsilon_0 k_BT}}_{\text{Debye}}, 
\underbrace{\tfrac{4}{\varepsilon_0}\Big(\tfrac{3 n_{\rm exc}(t)}{\pi}\Big)^{1/3}}_{\text{degenerate TF}}\Big],
$$

evaluated on the gathered excited density each step — each formula
*over*estimates $\kappa^2$ outside its own regime, so min() selects the valid
branch at every $n$. Registry-gated `dyn_lambda_ok`: **GaAs only**. For Si,
$\lambda=0$ is the *correct* physics (Burt's dynamical argument [L90]: the ~1 eV
Auger/II transition frequency far exceeds the carrier plasma frequency — the
static free-carrier screen cannot follow).

#### 4.3.4 k-local fallback and BGR

Ring off ⇒ the k-local Lindblad closure (Rosati–Iotti–Dolcini–Rossi
factorization into two frozen-rate amplitude-damping channels at the same
$\mathbf k$; partner populations as clamped scalar factors) — the
first-generation approximation, kept for folded cells. `yn_sbe_bgr_threshold`
moves its $E_{\rm th}$ with density (§2.3) — only when $\Sigma^{\rm HF}$ is off.

### 4.4 Auger recombination $\mathcal D_{\rm Aug}$ (`yn_sbe_auger`) — three representations

#### 4.4.1 Nonlocal ring Auger = the exact time-reverse of §4.3.2

Same quadruples, same $\gamma |V(\mathbf q)|^2 \delta_\sigma$ weight, **reversed
occupation factors**:

$$
P^{\rm rev} = \Big[1-\frac{f_{v\mathbf k_2}}{\texttt{occ}}\Big]_+ 
\frac{f_{c\mathbf k_1'}}{\texttt{occ}} \frac{f_{c\mathbf k_2'}}{\texttt{occ}},
$$

and negated stencil ($+$hot, $-c_1'$, $+$valence, $-c_2'$): two conduction
electrons meet a hole; one recombines, the released energy promotes the other to
the hot state. **No separate Auger coefficient exists or is needed** — the rate
scale *is* the cited II magnitude, because Auger and II share $|M|^2$ and differ
only in occupations (microreversibility). **Detailed balance is exact**: for
Fermi–Dirac occupations on an energy-conserving quadruple,

$$
f_1 f_2 (1-f_3)(1-f_4) = (1-f_1)(1-f_2) f_3 f_4
\quad\Longrightarrow\quad \mathcal D_{\rm II} + \mathcal D_{\rm Aug}\big|_{\rm FD} = 0,
$$

the equilibrium-fixed-point unit test. Hot-state gain capped at
$(\texttt{occ}-f)(1-e^{-\Gamma\tau})$ (no overfill).

#### 4.4.2 Graphene: the 2D Rana channel (gapless — its own mathematics)

No gap ⇒ no threshold law; the collinear-collapsed 2D Coulomb integrals of
[R07 = Rana, PRB 76, 155431 (2007)] on the Dirac spectrum
$E_\pm(\mathbf k)=\pm v|\mathbf k|$ ($v = 1\times10^8$ cm/s, $g=4$):

carrier density inversion (per branch):
$ n(\mu) = \frac{g}{2\pi}\int_0^\infty k f \big(\tfrac{vk-\mu}{k_BT}\big) dk
 \Rightarrow  \mu_c(n),\ \mu_v(p)$ (bisection);

screening [R07 Eq. 13]:
$ Q_{\rm TF} = \frac{4 k_BT}{\varepsilon_r v^2}\Big[\ln \big(e^{\mu_c/k_BT}+1\big)+\ln \big(e^{-\mu_v/k_BT}+1\big)\Big];$

CCCV recombination rate per area [R07 Eq. 14, the √ in the numerator]:

$$
R_{\rm CCCV} = \frac{1}{v}\int_0^\infty \frac{dk_1}{2\pi}\int_0^\infty \frac{dk_2}{2\pi}\int_{k_2}^\infty \frac{dQ}{2\pi} 
|M|^2 \sqrt{(k_1{+}Q)(Q{-}k_2) k_1 k_2} 
f_c(k_1)f_c(k_2)\big[1{-}f_c(k_1{+}Q)\big]\big[1{-}f_v(Q{-}k_2)\big],
$$

$$
|M|^2 = M_d^2 + M_e^2 + (M_d-M_e)^2,\qquad
M_d = \frac{2\pi}{\varepsilon_r (Q+Q_{\rm TF})},\quad
M_e = \frac{2\pi}{\varepsilon_r (|Q+k_1-k_2|+Q_{\rm TF})};
$$

CVVV = the hole mirror ($n\leftrightarrow p$); the generation partners
$G$ [Eq. 17] have the occupations reversed. The live channel applies the **net**

$$\frac{dn}{dt} = -(R_{\rm CCCV}+R_{\rm CVVV}) + (G_{\rm CVCC}+G_{\rm VCCC})  \equiv  -\frac{n-n_0(T)}{\tau_r(n)}$$

as a uniform-fractional CB→VB population transfer (VB→CB when $G>R$ —
thresholdless **carrier multiplication**), quasi-Fermi levels re-inverted from
the gathered populations every step, saturation-capped against both source and
destination, trace exact. At equilibrium $R=G$ identically. Validated:
$\tau_r(10^{12} {\rm cm^{-2}}, 300 {\rm K}, \varepsilon_r{=}10) = 1.48$ ps
(R07 Fig. 4).

#### 4.4.3 k-local $Cn^3$ legacy

Ring off: $\gamma = C n^2$ per carrier ($R=Cn^3$), density-gated, two
amplitude-damping maps (recombination + hot promotion). **No material ships a
verified default $C$** — requires explicit `sbe_auger_c_cm6s` (provenance gate).

### 4.5 Carrier–carrier thermalization $\mathcal D_{\rm ee}$ (`yn_sbe_eeh`)

The e-e/e-h collision integral drives the carrier subsystem toward a hot
Fermi–Dirac **without changing its particle number or energy** (it cannot relax
energy to the lattice — that is e-ph's job). Implemented as the exactly-conserving
relaxation: per $\mathbf k$, fit $(\beta,\mu)$ of

$$
f^{\rm FD}_a = \frac{1}{e^{\beta(E_a-\mu)}+1}\quad\text{s.t.}\quad
\sum_a f^{\rm FD}_a = \sum_a \frac{f_a}{\texttt{occ}},\qquad
\sum_a E_a f^{\rm FD}_a = \sum_a E_a \frac{f_a}{\texttt{occ}},
$$

then mix with $\alpha = 1-e^{-\nu_{cc}\tau}$:

$$
\tilde\rho_{aa} \to (1-\alpha) \tilde\rho_{aa} + \alpha \texttt{occ} f^{\rm FD}_a,
\qquad \tilde\rho_{ab} \to (1-\alpha) \tilde\rho_{ab}\ (a\neq b),
$$

a convex combination of CPTP maps ⇒ CPTP; the off-diagonal factor is the
excitation-induced dephasing (EID). Rate scale $\nu_{cc}$ =
`sbe_eeh_nu_sat` (cited $10^{13}$–$10^{14} \rm s^{-1}$ at
$n=10^{17}$–$10^{19} \rm cm^{-3}$, Goodnick–Lugli / Fischetti–Laux); screening
from Part G (§4.6). Provenance-gated (forbidden for CdS/graphene — no cited rate).

### 4.6 The screening library (Part G — used by §2, §4.3, §4.4, §4.5)

$$
\varepsilon_{\rm TF}(q) = 1+\frac{\kappa^2}{q^2};\qquad
\kappa^2_{\rm Debye} = \frac{4\pi n}{\varepsilon_{\rm bg} k_BT};\qquad
\kappa^2_{\rm TF,deg} = \frac{4}{\varepsilon_{\rm bg}}\Big(\frac{3n}{\pi}\Big)^{1/3};
$$

$$
\varepsilon_{\rm Lind}(q) = 1 + \frac{\kappa^2}{q^2} F \Big(\frac{q}{2k_F}\Big),\quad
F(x) = \frac12 + \frac{1-x^2}{4x}\ln\Big|\frac{1+x}{1-x}\Big|
\ \text{(static Lindhard, default e-e screen)};
$$

$\varepsilon_{\rm CDRB}(q)$ as in §4.3.2; $\omega_p^2 = 4\pi n/(\varepsilon_\infty m^*)$;
LOPC phonon-plasmon branches (GaAs-dynamic, material-specific). The static screen
is *wrong* at the ~1 eV Auger frequency scale — hence Si's $\lambda=0$ (§4.3.3)
and the CDRB *model* dielectric (valence response) rather than a free-carrier
static screen inside $|V(q)|^2$.

---

## 5. Observables

**Current (gauge-invariant):**
$ \mathbf J(t) = \frac{1}{\Omega_{\rm cell}}\sum_{\mathbf k} w_{\mathbf k} 
\mathrm{Tr}\big[(\boldsymbol\pi(\mathbf k) + \mathbf A(t)) \rho(\mathbf k,t)\big]$ —
the $\mathbf A$ diamagnetic term is internal, inter/intra-band compensation is
exact in the velocity gauge (no perturbative splitting). HHG spectrum
$\propto|\omega J(\omega)|^2$; the intra/inter split is meaningful only in the
Houston basis.

**Excited carriers:** $n_{\rm ex}(t) = \frac{1}{\Omega_{\rm cell}}\sum_{\mathbf k}
w_{\mathbf k}\sum_{a>N_v} f_a(\mathbf k,t)$ (Houston populations; per *primitive*
cell volume — the al_vec caveat in the wiki/00 decisions log).

**k-resolved maps:** Houston lowest-CB population per k (`_sbe_nex_k`), the
4-level real-carrier file (`_sbe_nex_k_lev_real`), coset-unfolded populations on
folded cells; rendered by `plot_sbe_results.py` (see
[`09_plotting_and_analysis.md`](09_plotting_and_analysis.md)).

---

## 6. Invariants and the flag map (the contract every channel obeys)

| term | trace $\sum\mathrm{Tr}\rho$ | positivity | energy | number of carriers |
|---|---|---|---|---|
| unitary $[H,\rho]$ | exact | exact | driven by the field | conserved |
| $\Sigma^{\rm HF}$ | exact (Hermitian) | exact | mean-field shift | conserved |
| KZ dephasing | exact | exact (PSD kernel) | untouched (diagonals untouched) | conserved |
| e-ph | exact by stencil | §3.1a maps | **relaxes to the lattice** (ħω per hop, detailed balance) | conserved |
| impact ionization | exact by stencil | §3.1a maps | conserved within $\delta_\sigma$ | **+1 pair per event** |
| ring Auger / Rana | exact by stencil | §3.1a maps | conserved within $\delta_\sigma$ / rate-model | **−1 pair per event** (or +, CM side) |
| carrier–carrier | exact | convex mix | **conserved** (FD fit constraint) | conserved |

Master flag map: clean run = *no* `&sbe` flags (pure unitary, the baseline);
`yn_sbe_superres='y'` is the backbone that upgrades e-ph/II/Auger to their
momentum-resolved nonlocal forms and enables the graphene channels. The
per-material availability matrix (provenance gates) lives in
[`02_constants.md`](02_constants.md); the maximally-loaded valid configuration of
every material is `samples/exercise_x11_full_dissipation_showcase/`.
