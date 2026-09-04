# Graphene sheet solver: gapless-cone dissipation with collisional memory, a two-temperature Coulomb sector, and a self-consistent sheet field

**Status: implemented and unit-tested (2026-09-04); calc-validation in exercise x14. Maintained together with `wiki/10` §8.11.**
*Written in the form of a methods paper so that it can be lifted into a manuscript; every equation is the one the code evaluates, every constant is cited, every claim carries its test.*

---

## Abstract

We describe the graphene branch of the SALMON2-TROUT semiconductor-Bloch-equation (SBE) solver, built to study field-induced transparency / absorption of a monolayer graphene sheet under 1–100 kV/cm single-cycle THz transients (the maintainer's DAST source) and near-infrared pulses. Three elements were added to the general velocity-gauge, completely-positive (GKLS) SBE machinery: (i) a **collisional-memory (non-Markovian) treatment of the ring dissipators adapted to the gapless Dirac cone** — the electron–phonon sectors keep their phonon-line kernels while the Coulomb (Auger / carrier-multiplication) sector, which on the cone is a *global* rate model, receives the 2D Dirac-plasmon line of the instantaneous electron–hole plasma as its memory kernel; (ii) a **two-temperature description of the Coulomb sector**, in which a carrier temperature $T_e$ and the quasi-Fermi levels are read from the first two moments of the gathered Dirac-cone populations while the lattice stays at the phonon-bath temperature and cools the carriers through the phonon channel; (iii) a **self-consistent sheet field** (radiation reaction) in the single-cell driver, so that the total field written by the solver *is* the transmitted field and the transmission coefficient follows from the field before and after the sheet. A fourth element turned out to be indispensable at THz: a **parameter-free pure-gauge restoration of the velocity-gauge current** — a truncated basis cancels only part of the diamagnetic current of the filled π band (70 % for two bands), the remainder, proportional to $A=E/\omega$, turns the sheet into a plasma mirror, and the fitted linear correction first tried over-corrects at high field and makes the self-consistent sheet unstable; subtracting the adiabatic ground-state current of the same truncated Hamiltonian removes the artifact exactly at every field with no adjustable quantity (§6a). We give the equations, the numerical realization, the cost scaling ($O(N_k^2)$ for the graphene ring), the k-mesh rules (resonance-shell resolution; Dirac point on the half-shifted Monkhorst–Pack mesh only for odd multiples of 3; Zener excursion $A_0$ versus mesh spacing), the unit tests that pin each piece, and the calculation-level validations. A level check performed on the way exposed and removed a spurious 0.21 eV gap at the Dirac point of the previously used 7-plane-wave empirical pseudopotential basis.

---

## 1. Scope and notation

Monolayer graphene is represented by the π/π* pair of the Ramanujam local empirical pseudopotential (EPM) in the 2-atom hexagonal primitive cell embedded in a 20 Å vacuum slab; the SBE runs on `nstate = 2` bands. Hartree atomic units are used throughout ($e=\hbar=m_e=1$, $c = 137.036$); the sheet lies in the $xy$ plane, the driving field is in-plane. $N_k$ denotes the number of k-points of the $N\times N\times1$ Monkhorst–Pack (MP) mesh; $A_{2D}=(\sqrt3/2)a^2$ is the primitive-cell area, $L_z$ the slab height along the vacuum axis, so the cell volume is $V=A_{2D}L_z$.

The ground-state basis is the shell-complete set $|\mathbf G|^2\le 29.4$ a.u. (43 plane waves). The earlier 7-plane-wave set ($|\mathbf G|^2\le 2.94$) is **not closed under the little group $C_{3v}$ of K** (the rotation about K maps a first-shell vector onto the second-shell vector $\mathbf b_2-\mathbf b_1$ absent from the set); the symmetry protection of the Dirac degeneracy is then lost and a spurious gap of 0.2125 eV opens at K (identical in the Python reference and in the Fortran band path). With 43 plane waves the gap is $<10^{-5}$ eV, $v_F = 0.960\times10^6$ m/s, and the Γ-bottom / M-dip values fall inside the thesis acceptance windows (`tests/test_graphene_dirac_levels.py`).

## 2. Master equation and channels

The density matrix $\rho_{\mathbf k}(t)$ obeys the velocity-gauge SBE in GKLS form,
$$
\dot\rho_{\mathbf k} = -i\,[H_{\mathbf k}(t),\rho_{\mathbf k}] + \sum_c \mathcal D_c[\rho],\qquad
H_{\mathbf k}(t) = \varepsilon_{\mathbf k} + \mathbf A(t)\cdot\boldsymbol\pi_{\mathbf k} + \tfrac12 A^2,
$$
propagated with the fourth-order commutator-free Magnus (CF4) scheme in Strang splitting with the dissipators $\mathcal D_c$, which act in the instantaneous Houston (field-dressed) basis of $H_{\mathbf k}(t)$ (`wiki/03`, `wiki/08`). On graphene the "ring" (inter-k) channels are:

* **electron–phonon**, the two Kohn-anomaly optical modes $E_{2g}$ (Γ, 196 meV) and $A_1'$ (K, 160 meV) with $\langle g^2_\Gamma\rangle=0.0405$, $\langle g^2_K\rangle=0.0994$ eV² and the GW factor 2 on the K mode [1,2], plus the quasi-elastic acoustic deformation-potential mode ($D=16$ eV, $v_{ph}=2\times10^6$ cm/s [3], grid-resolved $q$, Thomas–Fermi screened);
* **Coulomb (Auger recombination / carrier multiplication)**, the Rana rate model [4]: CCCV/CVVV recombination $R$ and their CVCC/VCCC time-reverses $G$ evaluated on the instantaneous quasi-Fermi levels of the gathered sheet densities $n$, $p$, applied as a uniform-fractional, trace-exact, bounded population transfer $R-G$. Impact ionization with a gap threshold, the carrier–carrier Fermi–Dirac fit and Kuhn–Zurek dephasing are not defined on a gapless cone and are refused by the code.

The Coulomb balance $R=G$ holds iff the electron and hole quasi-Fermi levels coincide; for a symmetric pair population this is $\mu=0$, i.e. the intrinsic density of the cone at the temperature $T$ the rates are evaluated at,
$$
n_0(T)=n_i(T)=\frac{\pi}{6}\left(\frac{k_BT}{\hbar v_F}\right)^2 \;=\; 8.08\times10^{10}\ \mathrm{cm^{-2}}\ \text{at 300 K}. \tag{1}
$$
Below $n_i$ the plasma net-multiplies, above it net-recombines: the pair population saturates at $n_i(T)$ (`tests/test_rana_saturation.f90`: root of $R-G$ at $n_i$ to $10^{-4}$, two-sided monotone CPTP relaxation, $T^2$ law).

## 3. Collisional memory on the gapless cone (the "2D colmem analog")

### 3.1 Motivation
The Markovian ring dissipators scatter the reversible field-induced admixture of conduction character ("dressing") as if it were real population (`wiki/10` §1). For gapped materials this was cured in three sectors (§8.7–8.10 of `wiki/10`): a memory kernel in the coherence damping, a memory filter on the populations that feed the collision kernels, and a dressed reference for the carrier measure. Graphene had been excluded by a guard. Two features of the cone required a dedicated version: its population channel is the **Coulomb rate model on the global densities $n,p$**, so the virtual share inflates the very quantities the R07 rates are evaluated on; and the memory line of a Coulomb collision is not a phonon energy but the **plasma response** — screening builds up on the inverse plasma frequency [5].

### 3.2 Kernels
All memory filters have the form used in `wiki/10` §8.6–8.9: a set of Lorentzian lines $\{c_j,\mu_j\}$, $\mu_j = 1/\tau_c \pm i\omega_j$ with the thermal split $(N_j+1):N_j$, the common width $1/\tau_c=\sigma_E$ (the ring's energy-matching width), and the discrete Markov anchor $R(0)=1$, so that a constant quantity is a machine-exact fixed point (calibrated rates untouched) while a modulation at frequency $\omega$ is transmitted with $|R(\omega)|=|\sum_j c_j/(\mu_j+i\omega)|$.

| sector | line set |
|---|---|
| e-ph coherence damping (`yn_sbe_colmem`) | graphene phonon table: $E_{2g}$, $A_1'$, acoustic |
| e-ph population source (`yn_sbe_colmem_pop`) | same |
| **Coulomb source densities $n,p$** | **2D Dirac plasmon** $\omega_{pl}(n,p)$ at $q=Q_{TF}$ |

The long-wavelength plasmon of the Dirac gas [6], $\omega_{pl}^2(q)=2e^2E_Fq/(\kappa\hbar^2)$ (Gaussian units), is generalized to the two-component electron–hole plasma by adding the Drude weights, each branch's $E_F$ becoming the finite-temperature intraband Drude weight [7]
$$
W(\mu)=2k_BT\,\ln\!\big[2\cosh(\mu/2k_BT)\big]\;\to\;|\mu|\ (\text{degenerate}),\quad 2k_BT\ln2\ (\text{intrinsic}),
$$
and it is evaluated at the collision's own screening momentum, the Thomas–Fermi vector of R07 Eq. (13),
$$
\omega_{pl}^2=\frac{2\,[W(\mu_c)+W(\mu_h)]\,Q_{TF}}{\varepsilon_r},\qquad
Q_{TF}=\frac{4k_BT}{\varepsilon_r v_F^2}\ln\!\big[(e^{\mu_c/k_BT}+1)(e^{\mu_h/k_BT}+1)\big]. \tag{2}
$$
At 300 K and $\varepsilon_r=10$: $\hbar\omega_{pl}=31$ meV (intrinsic), 35 meV at $n=p=10^{11}$, 133 meV at $10^{12}$ cm⁻² — the phonon scale, all $\ll 2\hbar\omega_{\rm laser}$ for near-IR drives (`tests/test_colmem_2d.f90`: degenerate limit $\omega_{pl}^2=2E_FQ_{TF}/\varepsilon_r$ to $10^{-6}$; a $2\omega$ breathing at 0.8 eV transmitted at $0.17=|R(2\omega)|$). No new input variables and no new free parameters enter: $T$, $\varepsilon_r$ and the $\mu$'s are the Rana channel's own.

### 3.3 Composition
Per ring step: dressed-reference subtraction (basis level) → ring gate → the raw Coulomb source moments $(n,p,\varepsilon)$ are captured → the phonon-line filter replaces the gathered populations for the e-ph kernel → the Coulomb channel filters $n,p$ with the plasmon line and evaluates $R-G$ on them (`rana_auger_dpop(n2d_in,p2d_in)`), while its transfer stencils and the CPTP limiter still use the instantaneous populations. Trace is exact and positivity unchanged. Kuhn–Zurek dephasing remains forbidden on graphene (many-body coherence loss).

## 4. Two-temperature Coulomb sector

### 4.1 Model
The R07 rates are quasi-equilibrium expressions in $(\mu_c,\mu_h,T)$. Evaluating them at the lattice temperature under-estimates the balance density of a hot plasma by $(T_e/T_L)^2$ (Eq. 1). We therefore read a **common carrier temperature** from the gathered populations. With the Dirac-point energy $E_D$ (midpoint of the instantaneous band edges, exact by electron–hole symmetry), the moments per unit area are
$$
n=\frac1{N_kA_{2D}}\sum_{\mathbf k} f_{c\mathbf k},\quad
p=\frac1{N_kA_{2D}}\sum_{\mathbf k}(2-f_{v\mathbf k}),\quad
\varepsilon=\frac1{N_kA_{2D}}\Big[\sum_{\mathbf k} f_{c\mathbf k}(\epsilon_{c\mathbf k}-E_D)+\sum_{\mathbf k}(2-f_{v\mathbf k})(E_D-\epsilon_{v\mathbf k})\Big],
$$
and $(k_BT_e,\mu_c,\mu_h)$ solve
$$
n=n_D(\mu_c,T_e),\quad p=n_D(\mu_h,T_e),\quad \varepsilon=\varepsilon_D(\mu_c,T_e)+\varepsilon_D(\mu_h,T_e), \tag{3}
$$
$$
n_D(\mu,T)=\frac{g}{2\pi}\!\int_0^\infty\!k\,f\!\Big(\frac{v_Fk-\mu}{k_BT}\Big)dk,\qquad
\varepsilon_D(\mu,T)=\frac{g}{2\pi}\!\int_0^\infty\!k\,(v_Fk)\,f\!\Big(\frac{v_Fk-\mu}{k_BT}\Big)dk,\quad g=4,
$$
with the closed forms $n_D(0,T)=\frac{\pi}{6}(k_BT/v_F)^2$ and $\varepsilon_D(0,T)=\frac{2}{\pi}\cdot\frac{3\zeta(3)}{2}\,(k_BT)^3/v_F^2$. At fixed $(n,p)$ the energy is monotone in $T$, so Eq. (3) is solved by a bisection in $\ln T$ with the density inversions nested inside (`dirac_fit_te`). $T_e$ is clamped from below at the lattice temperature (carriers are not colder than the bath) and the fit falls back to the bath when there are no carriers.

### 4.2 Use
$k_BT_e$ replaces the bath temperature in every carrier-side quantity of the Coulomb sector: the R07 integrals, $Q_{TF}$, the plasmon line and its thermal split. The phonon Bose factors keep the lattice temperature. **No separate rate equation for $T_e$ is integrated**: heating by the field and cooling by optical/acoustic phonon emission are already contained in the SBE kinetics; the two-temperature model is realized by *reading* $T_e$ from the distribution at every ring step. The time series $(T_e,\mu_c,\mu_h,n,p,n_i(T_e))$ is written to `*_sbe_te.data`. Interpretation: for a hot, non-thermal photo-excited shell the fit returns an effective temperature at which the R07 channel *generates* pairs toward $n_i(T_e)$ — the carrier-multiplication regime of graphene [8,9] — and switches to recombination as phonon cooling lowers $T_e$.

Unit test `tests/test_dirac_te_fit.f90`: the closed form of $\varepsilon_D(0,T)$; explicit 2D-mesh moments of thermal populations reproduce $n_D,\varepsilon_D$ to $0.5\%$ and the fit recovers $T$ to $10^{-4}$ and $\mu$ to $10^{-4}k_BT$; the degenerate case $E_F=0.3$ eV at 300 K; fallbacks; monotonicity.

## 5. Sheet electrodynamics

### 5.1 Boundary condition
At normal incidence the tangential $E$ is continuous across a current sheet and $H$ jumps by the sheet current $J_s$; with $Z_0=4\pi/c$ and a substrate of index $n_s$ behind the sheet,
$$
E_t=\frac{2E_{\rm inc}-Z_0J_s}{1+n_s},\qquad E_r=E_t-E_{\rm inc},\qquad
J_s=-J_m L_z, \tag{4}
$$
where $J_m$ is the electron current per cell volume written by the solver (its energy ledger is $dW=-\mathbf E\!\cdot\!\mathbf J_m\,V\,dt$, so the charge current is $-J_m$). Fluence-integrated,
$$
T=n_s\frac{\int E_t^2}{\int E_{\rm inc}^2},\quad R=\frac{\int E_r^2}{\int E_{\rm inc}^2},\quad A=1-T-R,\qquad
\frac{c}{4\pi}\big(E_{\rm inc}^2-E_t^2-E_r^2\big)=E_tJ_s\ \ (n_s=1) \tag{5}
$$
pointwise: $A$ is the Joule absorption of the sheet in the *local* field. For the universal conductance $\sigma=e^2/4\hbar$: $T=(1+\pi/2c)^{-2}=0.97746$, $A=\pi\alpha/(1+\pi\alpha/2)^2=0.02241$ (`tests/test_sheet_transmission.py`).

### 5.2 Radiation reaction in the velocity gauge
If the solver is driven by $E_{\rm inc}$ alone, its absorption $A_E=\int E_{\rm inc}J_s/F$ exceeds the sheet's by exactly
$$
S_{rr}=\frac{Z_0}{2}\frac{\int J_s^2\,dt}{F_{\rm inc}},\qquad A=A_E-S_{rr}, \tag{6}
$$
which is $O((Z_0\sigma)^2/2)\approx3\times10^{-4}$ on a mesh that resolves the response but becomes comparable to $A$ when a few discrete near-resonant k-points carry a large reactive current, and is not small at all for a THz-driven plasma ($Z_0\sigma/2\sim0.1$). The single-cell driver therefore propagates in the **local** field (`yn_sbe_sheet_field`): with $E=-\dot A$,
$$
\frac{dA_{\rm ind}}{dt}=-\frac{2\pi}{c}L_z\,J_m(t),\qquad A_{\rm tot}=A_{\rm ext}+A_{\rm ind},\qquad E_{\rm tot}=E_{\rm ext}+\frac{2\pi}{c}L_zJ_m, \tag{7}
$$
integrated explicitly with the current of the previous step (lag error $O(\Delta t\,\omega\,Z_0\sigma/2)$). `Ac_tot/E_tot` in `*_sbe_rt.data` are then the transmitted field, `E_ext` the incident one, the energy ledger uses $E_{\rm tot}$, and $A_{\rm ind}$ is part of the checkpoint state. The post-processor `transmission.py` detects this mode, uses $E_{\rm tot}$ directly and reports the deviation from the boundary-condition reconstruction (Eq. 4) as a consistency number.

## 6. Numerical realization

1. Ground state: in-SALMON EPM, 43 plane waves, $N\times N\times1$ half-shifted MP mesh; **K is on the mesh only for odd multiples of 3** ($(2i-N-1)/2N$ contains $2/3$ iff $N$ is odd) — 147 or 153 for production, 12/24/150 straddle K by half a spacing.
2. Per time step: CF4 unitary half-steps around the Strang-split dissipators; one Houston pass + gather per ring step; channels in order e-ph (with the phonon-line memory), Coulomb (with the plasmon-line source filter and $T_e$); sheet field update by Eq. (7); outputs.
3. Cost: the graphene ring is e-ph + Rana only — **$O(N_k^2)$**, exponential-bound (no $O(N_k^3)$ impact-ionization kernel). Measured $\approx40$ ms/step at $24^2$ on 4 threads $\Rightarrow\approx1.3$ s/step at $147^2$ on 48 threads. The 285-fs single-cycle THz transient + 100 fs tail at $\Delta t=0.1$ fs is 3844 steps: $\approx1.4$ h per dissipative run, minutes per coherent run.
4. Mesh rules. *Near-IR*: the resonance shell $k_{\rm res}=\hbar\omega/2\hbar v_F$ needs $\gtrsim3$ mesh points per radius ($N\ge150$ at 0.8 eV). *THz*: the relevant scale is the Zener excursion $A_0=E_0/\omega$ ($0.062$ a.u. at 100 kV/cm, 3.36 THz) against the spacing $|\mathbf b|/N$ ($0.0106$ at $N=147$) and the Landau–Zener tube $k_\perp^{LZ}=\sqrt{E/\pi v_F}=3.7\times10^{-3}$ a.u.; $A_0\ge2\,|\mathbf b|/N$ holds for $E_0\gtrsim30$ kV/cm at $N=147$; below that the pair creation is analytic (Eq. 8) rather than mesh-resolved.
5. Time step: the CF4 exponential is exact per step; $\Delta t$ is set by the field interpolation and the dissipator splitting — 0.1 fs for the THz drive, checked against 0.05 fs (x14).

## 6a. The velocity-gauge f-sum rule, the diamagnetic sheet current and its parameter-free removal (2026-09-04)

**Observation.** The first self-consistent THz run (24², $n_b=2$, 100 kV/cm) gave $T=0.15$, $R=0.85$: the empty sheet behaved as a plasma mirror. The sheet current was perfectly anti-correlated with the vector potential, $\mathrm{corr}(J_s,A_{\rm tot})=-1.000$ and uncorrelated with $E$, i.e. purely reactive, $J_s\simeq-\eta\,(N_e/A_{2D})\,A_{\rm tot}$ with $\eta=0.30$.

**Origin.** In the velocity gauge the current $J=\mathrm{Tr}[(\mathbf p+\mathbf A)\rho]/V$ carries the diamagnetic term $\mathbf A N_e/V$ for *every* electron of the filled π band. A uniform static $\mathbf A$ is a pure gauge, so this term must be cancelled exactly by the paramagnetic (interband) response — which happens only if the basis is complete: per band $n$ and direction $a$,
$$
S_n^a(\mathbf k)=\sum_{m\ne n}\frac{2|p^a_{nm}(\mathbf k)|^2}{\varepsilon_m-\varepsilon_n}=1-\frac{\partial^2\varepsilon_n}{\partial k_a^2},\qquad \big\langle S_n^a\big\rangle_{\rm full\ band}=1, \tag{9}
$$
and with $n_b$ bands the occupied-state average $\langle S^{(n_b)}\rangle<1$. The uncancelled fraction $\eta_a=1-\langle S^{a,(n_b)}\rangle$ multiplies $A=E/\omega$: harmless in the near-IR, decisive at 3 THz where $A_0$ is 60 times larger for the same field. This is the `wiki/06` basis-sufficiency issue in its most violent form — a 2D sheet whose whole valence band responds as if it were free.

**It is a property of the truncated basis, not of the mesh, the time step or adiabaticity.** This was checked directly (24², $n_b=2$, 1 kV/cm, no sheet field): the SBE current reproduces the *exact adiabatic response of the same truncated Hamiltonian* $H_{\mathbf k}(A)=\varepsilon_{\mathbf k}+A\,p^x_{\mathbf k}$ (ground state of the 2×2 problem at the instantaneous $A(t)$, Hellmann–Feynman current) to $1\times10^{-4}$ at every instant, both at 1 and at 100 kV/cm; the in-phase coefficient is $\eta_x=0.2995$ for $\Delta t=0.1,\,0.05,\,0.02$ fs alike, equal to Eq. (9) evaluated on the ground-state data; and driving at 0.1 and 0.4 eV instead of 14 meV moves it to 0.2985 and 0.2808, exactly the dispersive value $1-\big\langle\sum_m 2|p_{nm}|^2\Delta\varepsilon/(\Delta\varepsilon^2-\omega^2)\big\rangle$ (0.2985, 0.2806). The solver integrates its truncated model correctly; the model's ground state carries a current a complete basis would not.

**Convergence with $n_b$** (24², DAST 100 kV/cm, coherent, self-consistent sheet, no correction; static $\eta$ from the ground-state data, `sumrule_check.py`):

| $n_b$ | $\eta_x$ static | $A_{\rm tot}/A_{\rm ext}$ | $T$ | $R$ | $A$ |
|---|---|---|---|---|---|
| 2 | 0.300 | 0.20 | 0.147 | 0.851 | 0.002 |
| 4 | 0.097 | 0.53 | 0.646 | 0.321 | 0.033 |
| 8 | 0.036 | 0.90 | 0.937 | 0.013 | 0.049 |
| 16 | 0.030 | 0.92 | — | — | — |

The residual mirror at $n_b=16$ ($R\approx(Z_0\sigma/2)^2\approx(9.6\,\eta)^2\approx8\%$ at 3.36 THz) is still larger than the physical signal, so basis enlargement alone does not converge fast enough in the THz regime; the missing weight sits in bands $\gtrsim10$ eV up.

**The linear static correction is not admissible.** A first remedy, $J_a\to J_a-\eta_a(N_e/V)A_a(t)$ with the static $\eta$, removes the artifact at small $A$ but *over*-corrects at 100 kV/cm: the adiabatic ground-state current of the truncated model is **non-linear** in $A$ once $A$ is comparable to the k-distance of the mesh points nearest to K ($\partial J_{gs}/\partial A$ falls from $0.2995\,N_e/V$ to $0.268$ at $A=0.03$ and $0.276$ at $A_0=0.062$ a.u. on the 24² mesh — the level repulsion of the near-K pairs grows as $A(p_{cc}-p_{vv})$ closes their gap). An over-corrected sheet has a *negative* kinetic inductance, $J_{\rm phys}=+\kappa A$; with the self-consistent field this mode is unstable — after the pulse $A_{\rm ind}$ and $J$ grow together and the ledger shows gain ($T=1.21$, $R=0.50$, $A=-0.71$ at $n_b=2$; the same runaway appeared at $n_b=8$ in the dissipative runs). A fitted coefficient also sits badly in a first-principles code. It was removed.

**Remedy adopted: pure-gauge restoration, no parameter.** `yn_sbe_vg_sumrule='y'` now subtracts, inside `calc_current_bloch`, the adiabatic ground-state current of the *same* truncated Hamiltonian at the instantaneous vector potential,
$$
J_a(t)\ \to\ J_a(t)-J^{\rm gs}_a\big(\mathbf A(t)\big),\qquad
J^{\rm gs}_a(\mathbf A)=\frac{1}{V}\sum_{\mathbf k}w_{\mathbf k}\sum_n f^{(0)}_{n}\,\big\langle\phi_{n\mathbf k}(\mathbf A)\big|\,p_a+A_a\,\big|\phi_{n\mathbf k}(\mathbf A)\big\rangle
=\frac{\partial}{\partial A_a}\frac1V\sum_{\mathbf k}w_{\mathbf k}\sum_n f^{(0)}_n\Big[E_{n\mathbf k}(\mathbf A)+\tfrac12A^2\Big], \tag{10}
$$
where $\phi_{n\mathbf k}(\mathbf A)$, $E_{n\mathbf k}(\mathbf A)$ are the eigenvectors and eigenvalues of $H_{\mathbf k}(\mathbf A)=\varepsilon_{\mathbf k}+\mathbf A\cdot\mathbf p_{\mathbf k}$ (the propagator's own coupling, including any nonlocal or coset correction) and $f^{(0)}_n$ the ground-state occupations in energy order (adiabatic continuation). This is the statement that a uniform $\mathbf A$ is a pure gauge, imposed within the truncated space: in a complete basis $E_{n\mathbf k}(\mathbf A)=\varepsilon_n(\mathbf k+\mathbf A)$ and the BZ sum of Eq. (10) vanishes identically for every $\mathbf A$, so the correction is zero there; in the truncated basis its linear limit is the static $\eta_aN_eA_a/V$ and beyond it the exact non-linear ground-state current. Because the k-trace of $\rho_{\mathbf k}$ is conserved, the artifact resides in the ground-state sum alone — the subtraction is exact for **any** population (real pairs, their Drude current and the polarization currents are untouched: $J-J^{\rm gs}=\sum_{\mathbf k}\sum_n\big(f_n-f^{(0)}_n\big)\langle p+A\rangle_n+$ coherences). Cost: one $n_b\times n_b$ diagonalization per k per step, beside the six of the S4 propagator. The energy ledger, the sheet self-field (Eq. 7) and the outputs all use the restored current. Unit test `tests/test_vg_sumrule.f90` (Hellmann–Feynman identity at $A=10^{-5},0.05,0.3$, linear limit $=\eta$, non-linear departure); `sumrule_check.py` recomputes Eq. (10) from the ground-state files and reports the residual A-projection of any run.

**With the restoration** (24², DAST, coherent, self-consistent sheet): the post-pulse field is stationary ($A_{\rm ind}\to$ const, $J\to0$), $A_{\rm tot}/A_{\rm ext}=1.000$, and the residual A-projection of the sheet current is $<10^{-4}$ at all $n_b$:

| $n_b$ | $E_0$ [kV/cm] | $T$ | $R$ | $A$ | $\eta_{\rm dyn}$ residual |
|---|---|---|---|---|---|
| 2 | 1 | 1.000000 | 6×10⁻⁸ | 0.0000 | 0.0000 |
| 2 | 100 | 0.999996 | 1.3×10⁻⁶ | 0.0000 | 0.0000 |
| 4 | 100 | 0.999996 | 1.4×10⁻⁶ | 0.0000 | 0.0000 |
| 8 | 100 | 0.972 | 3.4×10⁻⁴ | 0.027 | 0.0001 |

The 24² mesh (smallest π–π* gap 0.78 eV, K half a spacing off the mesh) has no interband channel at 14 meV and no k-point inside the Landau–Zener tube, so the physical answer at this resolution is $T\simeq1$ — the earlier $A=0.049$ at $n_b=8$ was part of the artifact, not absorption. The absorption physics (Eq. 8, §7) lives on the 147² mesh with K on it (§8 and x14 README §7). Production recipe: $n_b=8$ **with** the restoration; the dissipators keep acting on the π/π* window only (`frozen_core_threshold_ev=-15`, `frozen_free_threshold_ev=14` around the Γ-anchored $E_F$).

**Remark for the experiment.** Real CVD samples *are* strongly absorbing at THz — but through the Drude conductance of doping-induced carriers ($\sigma_{dc}\sim20$–$50\,\sigma_{\rm univ}$, sheet transmission 0.5–0.7 already in the linear regime), not through the artifact above; reproducing that requires the FD$(E_F,T)$ initial state (§9). A measured transmission below 90 % "without subtracting the substrate" also contains the substrate's own Fresnel loss ($\approx11\%$ per face for $n\approx2$); the sheet contribution is the ratio to the bare-substrate reference, Eq. (4) with $n_s$. Two electronically decoupled layers (a large-angle twisted or incoherently stacked bilayer) at $d\ll\lambda$ sit in the same local field and add their sheet currents: `sbe_sheet_nlayers = 2` (Eq. 7 with $2L_zJ_m$); the naive estimate $T_2\approx T_1^2$ holds to second order in $Z_0\sigma$ ($T_1^2/T(2\sigma)=1-(Z_0\sigma)^2/2+\dots$, i.e. $\lesssim1\%$ low for $T_1\ge0.9$) but ignores that both layers see the field reduced by *both* currents, which matters for the field-dependent (non-linear) part.

## 7. Expected physics (analytic anchors for the validation)

*Near-IR, 1–100 kV/cm.* Perturbative interband regime: $A_0\ll|\mathbf b|/N$; pulse area $\theta\approx A_0v_F\tau_p/2=0.25$ rad at 100 kV/cm (8 cycles), so the coherent bleaching of the resonant shell is $\Delta A/A\approx-\theta^2/12\approx-0.5\%$ and, equivalently, the Pauli factor of the shell within the pulse bandwidth, $n/N_{\rm shell}\approx3.6\times10^{10}/6\times10^{12}$, gives $\Delta A/A\approx-1\%$: the transmission changes by $\sim10^{-4}$ absolute — the sheet stays at the universal $T=0.9775$.

*THz (DAST, 3.36 THz), 1–100 kV/cm, intrinsic sheet.* The field sweeps $\mathbf k+\mathbf A(t)$ through the Dirac point; pair creation is the massless Landau–Zener/Schwinger process with $P(k_\perp)=\exp(-\pi v_Fk_\perp^2/E)$ [10,11], i.e. a rate per area
$$
\Gamma=\frac{g}{4\pi^2}\,\frac{E^{3/2}}{v_F^{1/2}}\quad(g=4), \tag{8}
$$
which for the single-cycle transient gives $n\approx(1/\pi^2)A_0\sqrt{E_{\rm peak}/v_F}\approx1.5\times10^{12}$ cm⁻² per passage at 100 kV/cm (two passages per cycle, Stückelberg interference neglected), scaling as $E^{3/2}/\omega$: $\sim5\times10^{10}$ at 10 kV/cm. The created carriers absorb: coherently only their creation energy ($\sim2v_Fk_\perp^{LZ}\approx0.1$ eV per pair, a few per cent of the pulse energy at 100 kV/cm), with phonon scattering also the intraband (Drude) energy they acquire while accelerated to $v_FA_0\approx0.7$ eV. **For the intrinsic sheet the model therefore predicts THz-induced *absorption* growing with the field** (as observed in undoped graphene [12]); the self-induced *transparency* of doped CVD graphene [13–15] is the Drude-weight reduction of pre-existing carriers by heating and requires a doped / finite-temperature initial occupation, which the present solver does not have (§9).

## 8. Validation

| item | test / run | result |
|---|---|---|
| Dirac levels | `test_graphene_dirac_levels` | 43 PW: gap $6.5\times10^{-6}$ eV, $v_F=0.960\times10^6$ m/s, e–h asymmetry 0.8 %, linear to 0.5 eV; 7 PW: 0.2125 eV spurious gap (Fortran bandpath: 0.2125 / 0.0000 eV) |
| Rana balance | `test_rana_saturation` | $n_0=n_i(T)$ to $10^{-4}$, $T^2$ law, two-sided monotone CPTP saturation |
| plasmon line + filter | `test_colmem_2d` | limits, fixed point $10^{-12}$, $|R(2\omega)|$ transmission, Rana overrides |
| two-temperature fit | `test_dirac_te_fit` | mesh moments 0.5 %, $T_e$ to $10^{-4}$, degenerate limit, fallbacks |
| sheet BC | `test_sheet_transmission` | universal sheet $T=0.97746$, $A=0.02241$, energy identity, Fresnel |
| pipeline, 24² | x14 smoke (17 runs) | electrons = 2.000 every step; 9.5×10¹⁰ pairs/cm² at 100 kV/cm ⇔ ledger 5.7 %; Rana sign flips at $n_i$; Markovian ring +5 % "dephasing ionization", removed by the memory analog |
| VG f-sum rule / pure gauge | `test_vg_sumrule`; §6a tables | SBE = exact adiabatic response of the truncated $H(A)$ to $10^{-4}$ ($\Delta t$, $\omega$ checks); plasma mirror at $n_b=2$; linear static correction over-corrects and runs away with the sheet field; Eq. (10) restoration: residual $<10^{-4}$, sheet stationary, $A_{\rm tot}/A_{\rm ext}=1.000$ |
| calc-level (THz/near-IR, 147²) | x14 README §7 | filled from the local runs of 2026-09-04 (see the README table) |

## 9. Limitations and outlook

1. **Initial state at $T=0$, undoped.** `gs%occup` is integer filling; the thermal/doped Drude background of real samples and the bleaching physics of doped graphene need an FD$(E_F,T)$ initial occupation (which also generalizes the dressed-reference formula from full/empty bands to fractional ones).
2. **Mesh at THz.** Below $\sim30$ kV/cm the Landau–Zener tube is thinner than the $147^2$ spacing; the pair creation is then given by Eq. (8) analytically rather than by the mesh. A K-refined (non-uniform) mesh would need the momentum-map machinery of the acoustic mode to be generalized.
3. **Hot phonons.** The lattice is a fixed-temperature bath; the optical-phonon bottleneck of graphene (hot $A_1'$/$E_{2g}$ populations) is not included.
4. **Quasi-equilibrium form of the Coulomb sector.** The R07 rates and the two-temperature fit assume thermal branch distributions; the early, strongly non-thermal shell is mapped onto an effective $T_e$.
5. The sheet field is free-standing (vacuum both sides); a substrate index enters Eq. (4) trivially and can be added to the driver.

## References

1. S. Piscanec, M. Lazzeri, F. Mauri, A. C. Ferrari, J. Robertson, Phys. Rev. Lett. **93**, 185503 (2004).
2. M. Lazzeri, C. Attaccalite, L. Wirtz, F. Mauri, Phys. Rev. B **78**, 081406(R) (2008).
3. E. H. Hwang, S. Das Sarma, Phys. Rev. B **77**, 195412 (2008).
4. F. Rana, Phys. Rev. B **76**, 155431 (2007).
5. H. Haug, A.-P. Jauho, *Quantum Kinetics in Transport and Optics of Semiconductors* (Springer), build-up of screening.
6. E. H. Hwang, S. Das Sarma, Phys. Rev. B **75**, 205418 (2007).
7. L. A. Falkovsky, A. A. Varlamov, Eur. Phys. J. B **56**, 281 (2007).
8. T. Winzer, A. Knorr, E. Malic, Nano Lett. **10**, 4839 (2010).
9. D. Brida *et al.*, Nat. Commun. **4**, 1987 (2013).
10. D. Allor, T. D. Cohen, D. A. McGady, Phys. Rev. D **78**, 096009 (2008).
11. B. Dóra, R. Moessner, Phys. Rev. B **81**, 165431 (2010).
12. S. Tani, F. Blanchard, K. Tanaka, Phys. Rev. Lett. **109**, 166603 (2012).
13. H. Y. Hwang *et al.*, Phys. Rev. B **87**, 115413 (2013).
14. M. J. Paul *et al.*, New J. Phys. **15**, 085019 (2013).
15. Z. Mics *et al.*, Nat. Commun. **6**, 7655 (2015).
16. N. Boroumand *et al.*, Rep. Prog. Phys. **88**, 070501 (2025) — the non-Markovian framing (`wiki/10` §6).
17. R. Ramanujam, M.S. thesis, Arizona State University (2015) — the π-EPM form factors.
