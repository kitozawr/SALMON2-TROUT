# SALMON2-SBE — SALMON fork with CF4/Yoshida SBE propagation, CPTP decoherence & EPM ground states

This repository is a fork of the original [SALMON project](http://salmon-tddft.jp/), an open-source software package for *ab-initio* quantum-mechanical calculations of light-matter interactions. 

This fork extends SALMON's Semiconductor Bloch Equations (SBE) module with a **commutator-free Magnus 4 (CF4) / Suzuki-Yoshida exponential propagator**, a **strictly CPTP Kuhn-Zurek/Caldeira-Leggett decoherence model**, an optional **k-local impact-ionization Lindblad channel**, an optional **Coulomb (time-dependent Hartree–Fock) renormalization** for the extreme-THz regime (Golde–Kira–Meier–Koch), **frozen-core optimizations**, and a self-contained **local Empirical Pseudopotential Method (EPM)** ground-state solver (Cohen-Bergstresser, GaAs) that closes the EPM → SBE pipeline without external scripts.

## Contents

- [Key Fork Features](#key-fork-features)
- [The CF4 + Suzuki-Yoshida + CPTP Operator Splitting](#the-cf4--suzuki-yoshida--cptp-operator-splitting)
- [Configuration Parameters](#configuration-parameters)
  - [Real-time output frequency (`&analysis`)](#real-time-output-frequency-analysis)
  - [Plotting the real-time output (`plot_sbe_results.py`)](#plotting-the-real-time-output-plot_sbe_resultspy)
  - [EPM ground-state solver (`&epm`)](#epm-ground-state-solver-epm)
- [Examples](#examples)
  - [Minimal SBE Input Example](#minimal-sbe-input-example)
  - [Minimal EPM → SBE Pipeline Example](#minimal-epm--sbe-pipeline-example)
  - [Spinor (spin-orbit) EPM → SBE Pipeline Example](#spinor-spin-orbit-epm--sbe-pipeline-example)
  - [Band-structure calculation (`theory='dft_band'`)](#band-structure-calculation-theorydft_band)
- [Building & Continuous Integration](#building--continuous-integration)
- [References & Theoretical Background](#references--theoretical-background)
- [License](#license)


## Key Fork Features

### 1. Strictly CPTP Kuhn-Zurek/Caldeira-Leggett Decoherence
Phenomenological dephasing schemes (whether $-\rho_{nm}/T_2$ in the field-free basis, or double-commutator dissipators built from the instantaneous dressed Hamiltonian $H_{eff}=H_0+\mathbf{A}(t)\cdot\mathbf{p}$) generally fail to be **completely positive and trace preserving (CPTP)** for arbitrary parameters/timesteps, and can introduce artifacts such as spurious adiabatic-following relaxation in strong fields.

This fork instead implements a **Kuhn-Zurek/Caldeira-Leggett wave-packet dephasing model** that is *exactly* CPTP by construction:
* At every step the instantaneous (Houston/adiabatic) eigenbasis $U(t)$ of $H_{VG}(t)=H_0+\mathbf{A}(t)\cdot\mathbf{p}$ is computed, together with the branch (wave-packet) positions $X_a(t)$, propagated via their group velocities $V_a = (U^\dagger\boldsymbol{\pi}U)_{aa}+\mathbf{A}(t)$.
* The density matrix is rotated into this basis, $\tilde\rho = U^\dagger\rho U$, and dephased through an **exact Hadamard/Gaussian (RBF) kernel** $\tilde\rho_{ab}\leftarrow e^{-\lambda (X_a-X_b)^2\,\tau}\,\tilde\rho_{ab}$.
* By the Schoenberg/Bochner positive-definiteness of the Gaussian kernel and the Schur product theorem, this Hadamard map is CPTP for **any** $\tau\ge 0$ — no positivity violations, no ad-hoc clipping.
* The decoherence rate is set physically via $\lambda = k_B T/\tau_m$ (`sbe_decoh_temperature_k`, `sbe_decoh_tau_m_fs`).

### 2. CF4 + Suzuki-Yoshida Exponential Propagator
Replaces the previous ETDRK4/Taylor-4 propagators with a **commutator-free Magnus 4th-order (CF4)** exponential integrator evaluated on two-point Gauss-Legendre quadrature nodes, composed into a 4th-order scheme via the **Suzuki-Yoshida triple-jump** ($p_1=1.35120719196$, $p_2=-1.70241438392$):
* Each CF4 sub-step is realized as **two exact unitary rotations** $\rho \to e^{-i\Omega_2}e^{-i\Omega_1}\rho\, e^{+i\Omega_1}e^{+i\Omega_2}$, with $\Omega_{1,2}$ Hermitian combinations of the Hamiltonian sampled at the Gauss-Legendre nodes, exponentiated *exactly* via eigendecomposition (no Padé/Krylov truncation error — unitary to machine precision).
* The full step combines the coherent and dissipative parts via **Strang splitting**, $D(h/2)\circ\big[S_2(p_1h)\circ S_2(p_2h)\circ S_2(p_1h)\big]\circ D(h/2)$, with the Suzuki-Yoshida composition wrapping **only** the unitary part. This is essential for CPTP: a negative sub-step ($p_2 h<0$) is a harmless backward-time *unitary* rotation, but applying it to the dissipator would invert the sign of the Hadamard kernel exponent and break positive semi-definiteness.
* Branch positions $X_a$ are advanced using the **midpoint (average of endpoint) velocities**, matching the overall 4th-order accuracy of CF4 (a forward-Euler update would degrade the scheme to 1st order).

### 3. Frozen Core / Active Subspace Optimization
For systems with many deep bands (e.g., 80 bands where 60 lie below -20 eV), evaluating the nonlinear commutator $[V, \rho]$ is computationally wasteful. This fork introduces an **Active Subspace** projection:
* Deep core and high-energy free bands are "frozen" and evolve purely under the exact linear operator $L$.
* Nonlinear light-matter interactions are computed exclusively in the active subspace (e.g., $20 \times 20$ instead of $80 \times 80$ ZGEMM calls).
* Yields an additional **~30× speedup** for the nonlinear evaluation step without sacrificing physical accuracy.

### 4. Exact Current Operator
Computes the total current as $J = \text{Tr}[(\mathbf{p} + \mathbf{A}) \rho]$ (in atomic units) without relying on perturbative expansions, ensuring proper inter/intra-band compensation in the Velocity Gauge. Hermiticity stabilization (`ρ = ρ†`) is enforced at each step to maintain real-valued currents and FFT stability.

### 5. Local Empirical Pseudopotential Method (EPM) ground states
A self-contained local-EPM ground-state solver (`theory='epm'`, `src/epm`) that computes the Cohen-Bergstresser band structure and momentum matrix elements for zincblende GaAs directly in SALMON, and writes `SYSNAME_k.data`/`SYSNAME_eigen.data`/`SYSNAME_tm.data` in exactly the format read by `gs_info_ssbe` — closing the EPM → SBE pipeline end-to-end without external scripts (`rvnl_tm` is written as identically zero, since a local pseudopotential has no nonlocal velocity correction).

### 6. Spinor (spin-orbit split) EPM input + `yn_sbe_spinor`
The Python EPM reference (`epm_gaas_reference.py`, hardcoded flag `INCLUDE_SPIN_ORBIT = True`) promotes the scalar $N_{PW}$ problem to the **spinor $2N_{PW}$ problem**: the plane-wave basis is doubled to $|G,s\rangle$, and

$$\hat H_0^{\rm spinor}(k) = \hat H^{\rm loc}(k)\otimes\mathbb 1_2 + \hat H_{SO}(k),$$

with the projected Weisz/Bloom-Bergstresser spin-orbit operator (Chelikowsky-Cohen form) whose single strength constant $\mu$ is auto-calibrated at $\Gamma$ to the GaAs split-off gap $\Delta_0 = 0.341$ eV ($\Gamma_8$–$\Gamma_7$). Because $\hat H_{SO}$ is **nonlocal**, the velocity acquires the mandatory correction $\hat v_{SO} = -i[\hat r,\hat H_{SO}] = \nabla_k\hat H_{SO}$, which the script computes **analytically** (verified against finite differences) and writes into block 2 (`rvnl_tm`) of the `_tm.data` file — so a spinor SBE run must set `yn_vnl_correction='y'` to use the full $\hat\pi^{\rm spinor} = \hat p\otimes\mathbb 1_2 + \hat v_{SO}$ consistently in $H_{VG}$, the decoherence branch velocities and the current.

On the SBE side the new `&sbe` flag **`yn_sbe_spinor`** switches the solver to such spin-orbit split input files: occupations become **1 per spinor band** over the first `nelec` bands (instead of 2 per band over `nelec/2`), and every `nelec/2`-derived index (Fermi level for the frozen core, lowest conduction band for `_sbe_nex_k`, valence trace for `_sbe_nex`, automatic minimum-gap search) consistently uses `nelec` valence bands. The spinor Bloch equation stays a **single** $2N_b\times 2N_b$ equation — spin-orbit couples the spin channels, so it does not factorize.

### 7. k-local impact ionization (optional Lindblad channel, `yn_sbe_impact_ionization`)
A fully **optional** carrier-multiplication channel for GaAs (off by default — it can be slow, and for a first estimate the purely coherent + dephasing code is sufficient). Hot conduction electrons with kinetic energy above the threshold ionize a valence electron, creating a "cold" e–h pair at the band edges, with the rate taken from the verified **Stobbe–Redmer–Schattke fit** (PRB 49, 4494) to a full Fermi-golden-rule calculation on the *same* Cohen-Bergstresser band structure as our EPM:

$$\gamma_{\rm St}(\varepsilon^{\rm kin}) = P\,(\varepsilon^{\rm kin}-E_{\rm th})^4\,\Theta(\varepsilon^{\rm kin}-E_{\rm th}),\qquad P = 2\times10^{12}\ {\rm s^{-1}eV^{-4}},\quad E_{\rm th}=2.1\ {\rm eV},$$

with $\varepsilon^{\rm kin}$ measured from the field-free CBM (the $\tfrac12 A^2$ scalar dropped in $H_{VG}$ is restored exactly here; by the Houston identity $\varepsilon^{\rm kin}=E_h(k+A)-E_{\rm CBM}$ — the scale on which the fit is defined). The quartic two-particle event $\hat A_h=\sqrt{\gamma_{\rm St}}\,c^\dagger_{h'}c^\dagger_{c_1}c_{v_1}c_h$ is closed **k-locally** (no momentum transfer — k-points stay independent, VG parallelism intact) and Hartree-Fock-factorized (two-particle closure, Rosati–Iotti–Dolcini–Rossi PRB 90, 125140) into two **frozen-rate amplitude-damping channels** in the same Houston basis as the Kuhn-Zurek dephasing (no extra ZHEEV): primary relaxation $h\to h'$ (the conduction branch closest to $\varepsilon_h-E_g$) and cold-pair creation $v_1\to c_1$, with partner populations and Pauli blockers entering as scalar factors clamped to $[0,1]$ — every map is **exactly CPTP** for any step. In the Boltzmann limit the diagonals reproduce the canonical impact-ionization collision integral $W_h=\gamma_{\rm St}\tilde\rho_{hh}\tilde\rho_{v_1v_1}(1-\tilde\rho_{c_1c_1})(1-\tilde\rho_{h'h'})$ (one e–h pair per event, $\dot n_c = \gamma_{\rm St} n_{\rm hot}$ in the dilute limit), and ionization additionally **destroys the coherences** of the participating branches — a decoherence channel of its own. A threshold gate keeps the cost at $O(N_C)$ comparisons per k-point while no populated branch exceeds $E_{\rm th}$ ("rare impact events").

**Declared limitations of the fit** (Stobbe): direction-averaged (their matrix elements are nearly isotropic — energy is the dominant variable); electron-initiated channel only (hole-initiated omitted); no phonon-assisted ionization, collisional broadening, or field-induced threshold softening (Quade–Schöll–Rossi: at MV/cm there is strictly no fixed threshold — near-threshold rates are underestimated, a known limitation); fit energy resolution $\delta E = 0.2$ eV (the $\Theta$ step is smoothed by a linear ramp of this width). Electron–electron scattering ($O(N_k^2)$, expensive) and Auger recombination ($\gamma_{\rm Auger}\sim10^6$ s$^{-1}$, negligible on sub-ps scales) are deliberately excluded.

### 8. Coulomb (time-dependent Hartree–Fock) renormalization (optional, `yn_sbe_coulomb`)
A fully **optional** Coulomb mean field for the extremely nonlinear THz regime, following **Golde–Kira–Meier–Koch** (*Phys. Status Solidi B* **248**, 863 (2011), Eqs. 4–5). At peak fields of a few MV/cm, carriers are driven across a large fraction of the Brillouin zone and the Coulomb interaction renormalizes both the band energies and the field coupling. In the multiband density-matrix form this is the time-dependent Hartree–Fock **exchange (Fock) self-energy**

$$\Sigma^{\rm HF}_{nm}(k) = -\sum_{q\neq k} V(k-q)\,\delta\rho_{nm}(q),\qquad V(p)=\frac{\texttt{strength}\cdot 4\pi}{\varepsilon\,\Omega_{\rm cell}\,N_k\,(|p|^2+\kappa^2)},$$

added to $H_{VG}$. The single commutator $-i[\,H_{VG}+\Sigma^{\rm HF},\rho\,]$ reproduces **both** of the paper's renormalizations at once: the diagonal part gives the renormalized single-particle energies $\tilde\varepsilon^\lambda_k=\varepsilon^\lambda_k-\sum_{q}V_{k-q}f^\lambda_q$, and the off-diagonal part gives the renormalized Rabi frequency $\Omega_k=\mathbf d_k\!\cdot\!\mathbf E_{\rm THz}+\sum_{q}V_{k-q}p_q$, with the $(1-f^e_k-f^h_k)$ Pauli-blocking factor emerging automatically from the commutator structure. Key design choices:

* **Basis.** $\Sigma^{\rm HF}$ is built and stored in the **velocity-gauge stationary-Bloch basis** in which $\rho_{nm}(k)$ is propagated. The convolution is gauge-covariant under the uniform Peierls shift $k\to k-\mathbf A(t)$ (the $\mathbf A$ cancels in $k-q$), so it is evaluated directly on the grid-$k$ density matrix with no transformation. Because $\Sigma^{\rm HF}$ is **added to $H_{VG}$**, the **Houston basis** (eigenbasis of $H_{VG}+\Sigma^{\rm HF}$) that the dissipative channels diagonalize automatically becomes the Coulomb-renormalized adiabatic basis — consistent with the paper's $\mathbf E$-renormalized (Houston) picture.
* **Equilibrium subtraction.** The EPM bands carry no explicit exchange, so the convolution uses the **deviation** $\delta\rho=\rho-\rho_0$ from the ground state ($\rho_0=\mathrm{diag}(\texttt{occup})$). All pieces of $\delta\rho$ vanish at $t=0$ (no excited electrons/holes, no polarization), so $\Sigma^{\rm HF}(t{=}0)=0$: the equilibrium gap stays exactly the EPM gap and **only** the carrier-induced (dynamical) renormalization is added.
* **CPTP-safe.** $\Sigma^{\rm HF}$ is Hermitian ($\rho$ Hermitian, $V$ real), so it enters as a coherent (unitary) generator and preserves trace and positivity; it is **frozen at $\rho(t)$** over each $h$ (mean-field predictor) and re-evaluated once per step.
* **Cost.** Unlike the k-local dephasing/ionization channels, the exchange sum is **non-k-local** — it couples all k-points, costing $O(N_k^2\,n_{\rm act}^2)$ compute (a screened-Coulomb convolution, parallelized across MPI ranks by the k-partition) and **one MPI all-gather of the active-band $\rho$ per step** ($O(n_{\rm act}^2 N_k)$ words). Both scale as $n_{\rm act}^2$, so always pair it with the **frozen-core active subspace** (§4: set `frozen_core_threshold_ev`/`frozen_free_threshold_ev`) to keep $n_{\rm act}$ small — this shrinks the per-step traffic to a few MB and the convolution to the gap-edge bands. The compute/communication ratio is $\sim N_k$ (the transfer is relatively cheaper on denser grids, where the convolution dominates), and the all-gather is verified MPI-correct (1-rank vs 2-rank populations are bitwise-identical). It is off by default; the screening $\kappa$ (`sbe_coulomb_screen_au`) regularizes the $q\to0$ tail and the $q=0$ self-term is excluded.

---

## The CF4 + Suzuki-Yoshida + CPTP Operator Splitting

The full (optionally Coulomb-renormalized) master equation advanced by the propagator is

$$\partial_t\rho(k,t) = -i\big[\,H_{VG}(k,t) + \Sigma^{\rm HF}[\rho](k,t)\,,\,\rho(k,t)\big] + \mathcal{D}_{\rm KZ}[\rho] + \mathcal{D}_{\rm II}[\rho],$$

with $H_{VG}=H_0(k)+\mathbf A(t)\cdot\boldsymbol\pi$ the velocity-gauge band Hamiltonian, $\Sigma^{\rm HF}$ the optional Coulomb exchange mean field (§8), $\mathcal{D}_{\rm KZ}$ the strictly-CPTP Kuhn–Zurek dephasing (§1) and $\mathcal{D}_{\rm II}$ the optional impact-ionization Lindblad channel (§7). It is advanced over a step $h$ as

$$\rho(t+h) = D(h/2)\circ\Big[S_2(p_1 h)\circ S_2(p_2 h)\circ S_2(p_1 h)\Big]\circ D(h/2)\,[\rho(t)]$$

The Coulomb self-energy $\Sigma^{\rm HF}$ is **frozen at $\rho(t)$** over the step (a non-k-local mean field; re-evaluated once per $h$) and folded into the Hamiltonians $H_1,H_2,H_{VG}$ of **both** the unitary $S_2$ and the dissipative $D$ blocks, so the Houston basis the dissipators diagonalize is the Coulomb-renormalized one.

**Unitary part $S_2(\tau)$ — CF4 on Gauss-Legendre nodes** (nodes $c_{1,2}=\tfrac12\mp\tfrac{\sqrt3}{6}$, weights $\alpha_{1,2}=\tfrac14\pm\tfrac{\sqrt3}{6}$):
* $H_1=H_{VG}(t+c_1\tau)$, $H_2=H_{VG}(t+c_2\tau)$
* $\Omega_1=\tau(\alpha_1 H_1+\alpha_2 H_2)$, $\Omega_2=\tau(\alpha_2 H_1+\alpha_1 H_2)$
* $\rho \to e^{-i\Omega_2}e^{-i\Omega_1}\,\rho\,e^{+i\Omega_1}e^{+i\Omega_2}$, each exponential built *exactly* from an eigendecomposition of the Hermitian generator (unitary to machine precision).

**Dissipative part $D(\tau)$ — Strang/Hadamard Kuhn-Zurek dephasing** (always applied with $\tau=+h/2 > 0$):
* Diagonalize $H_{VG}(t)\to U(t),\ \{E_a\}$ (Houston/adiabatic basis); $\tilde\rho = U^\dagger\rho U$
* $\tilde\rho_{ab}\leftarrow \exp[-\lambda(X_a-X_b)^2\tau]\,\tilde\rho_{ab}$ (PSD Hadamard/Gaussian kernel ⇒ exactly CPTP for $\tau\ge0$)
* Rotate back $\rho = U\tilde\rho U^\dagger$; update $X_a \mathrel{+}= \tfrac12(V_a(t)+V_a(t+h))\,h$, with $V_a=(U^\dagger\boldsymbol\pi U)_{aa}+\mathbf{A}(t)$.

**Why Yoshida wraps only the unitary part:** the middle Yoshida sub-step has $p_2 h<0$. For $S_2$ this is merely a unitary rotation run backwards in time — exact and always valid. For $D$, however, a negative $\tau$ would turn the Gaussian kernel $e^{-\lambda\Delta X^2\tau}$ into $e^{+\lambda\Delta X^2|\tau|}$, which is not positive semi-definite (violates the Schoenberg/Bochner criterion and the Schur product theorem) and would break CPTP. Hence $D$ is applied only with $\tau=+h/2$, via Strang splitting around the (always-safe) Yoshida-composed unitary block.

---

## Configuration Parameters

The `&sbe` namelist now accepts the following parameters:

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `sbe_decoh_temperature_k` | K | `-1.0d0` | Bath temperature $T$ for the Kuhn-Zurek/Caldeira-Leggett dephasing rate $\lambda=k_B T/\tau_m$. Both this and `sbe_decoh_tau_m_fs` must be `> 0` to enable decoherence. |
| `sbe_decoh_tau_m_fs` | fs | `-1.0d0` | Wave-packet momentum-relaxation time $\tau_m$ entering $\lambda=k_B T/\tau_m$. |
| `frozen_core_threshold_ev` | eV | `0.0d0` | Freeze bands below $E_F + \text{threshold}$. (Use negative values, e.g., `-15.0`). |
| `frozen_free_threshold_ev` | eV | `0.0d0` | Freeze bands above $E_F + \text{threshold}$. (Use positive values, e.g., `+20.0`). |
| `yn_sbe_spinor` | — | `'n'` | `'y'`: ground-state input files come from a **spinor (spin-orbit split)** system — occupation 1 per spinor band, `nelec` valence bands instead of `nelec/2`. Combine with `yn_vnl_correction='y'` when the dataset carries the $\hat v_{SO}=\nabla_k\hat H_{SO}$ correction in `rvnl_tm`. |
| `yn_sbe_impact_ionization` | — | `'n'` | `'y'`: enable the **k-local impact-ionization** Lindblad channel (Stobbe rate fit, GaAs). Fully optional; threshold-gated, so it costs ~nothing while no populated branch exceeds $E_{\rm th}$. |
| `sbe_ii_prefactor` | s⁻¹eV⁻⁴ | `2.0d12` | Stobbe fit prefactor $P$ in $\gamma_{\rm St}=P(\varepsilon^{\rm kin}-E_{\rm th})^4$. |
| `sbe_ii_threshold_ev` | eV | `2.1d0` | Ionization threshold $E_{\rm th}$ above the field-free CBM. |
| `sbe_ii_ramp_ev` | eV | `0.2d0` | Linear $\Theta$-smoothing width (the fit's energy resolution); `<= 0` gives a hard step. |
| `yn_sbe_coulomb` | — | `'n'` | `'y'`: enable the **Coulomb (time-dependent Hartree–Fock / exchange) renormalization** (§8, Golde–Kira–Meier–Koch). Non-k-local mean field, $O(N_k^2)$ per step — off by default, best on modest grids. |
| `sbe_coulomb_epsilon` | — | `12.9d0` | Background dielectric constant $\varepsilon$ screening the exchange kernel (GaAs default). |
| `sbe_coulomb_strength` | — | `1.0d0` | Overall scaling of the exchange kernel (set `0` to disable while leaving the flag on; `>1` to enhance). |
| `sbe_coulomb_screen_au` | Bohr⁻¹ | `0.0d0` | Yukawa screening $\kappa$ regularizing $V(q)\propto1/(q^2+\kappa^2)$; `0` = bare Coulomb with the $q=0$ self-term excluded. |

*Note: Internal conversions to atomic units (Hartree) are handled automatically (`kB_au`, `au_fs`).*

### Real-time output frequency (`&analysis`)

Real-time SBE propagation writes three diagnostic files (`SYSNAME_sbe_rt_energy.data`, `SYSNAME_sbe_nex.data`, `SYSNAME_sbe_nex_k.data`), each on its own cadence selectable in the `&analysis` namelist. The k-resolved file in particular can grow to gigabytes for dense k-grids/long runs, so its default stride is ten times coarser than the band-projection output:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `out_rt_energy_step` | `10` | Stride (in time steps) for `SYSNAME_sbe_rt_energy.data` (total energy / trace) and the stdout progress log. |
| `out_projection_step` | `100` | Stride for `SYSNAME_sbe_nex.data` (number of excited electrons/holes, summed over k). |
| `out_projection_k_step` | `1000` | Stride for `SYSNAME_sbe_nex_k.data` (Houston-basis population of the lowest conduction band, resolved per k-point). Defaults to 10× `out_projection_step` to avoid producing terabyte-scale output on dense k-grids; increase the stride (larger value) further for very large `nk`/`nt`. |

`SYSNAME_sbe_nex_k.data` reports, for every saved time `t`, one block of `nk` lines `ik, kx, ky, kz, population_lcb`, where `population_lcb = (W^\dagger \rho W)_{aa}` is the diagonal element of the density matrix rotated into the instantaneous Houston (adiabatic) eigenbasis $W$ of $H_{VG}(t)$ for the lowest conduction band $a = N_{elec}/2+1$ — i.e. the same gauge-independent basis used internally by the CPTP dephasing step. With `yn_sbe_spinor='y'` the lowest conduction band is $a = N_{elec}+1$ (the lower spin sub-band of the first conduction level).

### Plotting the real-time output (`plot_sbe_results.py`)

The repository root contains a self-contained `plot_sbe_results.py` script (matplotlib + numpy, not part of the Fortran build — copy it into the calculation directory and run it there). It scans the directory for `SYSNAME_sbe_rt_energy.data`, `SYSNAME_sbe_nex.data` and `SYSNAME_sbe_nex_k.data`, and produces (with no interactive windows, `Agg` backend):
* line plots of total energy and excited-electron/hole counts vs time;
* for `SYSNAME_sbe_nex_k.data`, one PNG per saved time step (the time value is encoded in the file name), each showing the Houston-basis lowest-conduction-band population as three 2D heatmap slices of the k-grid ($k_x$-$k_y$, $k_x$-$k_z$, $k_y$-$k_z$);
* a band-structure plot from `SYSNAME_k.data` + `SYSNAME_eigen.data` along a high-symmetry path (`--band-path`, default `L Γ X W K`), energies shifted to VBM = 0;
* a band-structure plot from `band.dat` (a `theory='dft_band'` run) vs path distance, energies shifted to the `--band-vbm` band index (default `nb//2`).

```sh
cp plot_sbe_results.py /path/to/calculation/
cd /path/to/calculation/
python3 plot_sbe_results.py            # writes PNGs into ./sbe_plots/
```

**Spinor (spin-orbit split) datasets:** the plotter detects a spinor `_eigen.data` automatically from the occupation column (1 electron per band instead of 2) and then **sums the spins of each level**: adjacent (Kramers partner) spin sub-bands are merged into one level — occupations summed ($1+1=2$ per valence level), level energy = mean of the pair. The band plot draws one solid curve per level on top of the faint spin-resolved sub-bands, so tiny Dresselhaus splittings don't render as doubled lines while the real spin-orbit splittings ($\Gamma_8$/$\Gamma_7$, $\Delta_0 = 0.341$ eV for GaAs) stay visible between levels. Control with `--spin-sum {auto,on,off}` (default `auto`):

```sh
python3 plot_sbe_results.py                 # auto-detects spinor input, sums spins per level
python3 plot_sbe_results.py --spin-sum off  # raw 2*Nb spin-resolved bands
```

**Folded vs unfolded bands.** The cubic 8-atom cell is a supercell of 4 primitive FCC cells, so the MP-grid band plot shows the primitive bands **folded 4-fold**: every cubic k-point carries the states of 4 primitive BZ points, and the conduction manifold appears as 4 overlaid copies of CB1/CB2/CB3 — dense crossings that are an artifact of the supercell representation, not of the physics (cf. band unfolding, Quan-Rybin-Scheffler-Carbogno, *PRB* 113, 085112 (2026)). The folding is **exact** here (the parity selection rule makes the Hamiltonian block-diagonal over the 4 FCC reciprocal sublattices to machine precision — asserted at runtime), so the clean primitive picture can be recovered exactly:

```sh
python3 epm_gaas_reference.py bandpath   # fast: writes SYSNAME_bandpath.data (no MP dataset)
python3 plot_sbe_results.py --only-bands # -> bandpath_*.png + bandpath_spin_splitting_*.png
```

The `bandpath` mode diagonalizes the FCC-sublattice blocks of the cubic Hamiltonian along the primitive path `L-Γ-X-W-K-Γ` (configurable constants at the top of the script) and the plotter renders the **unfolded** primitive-cell bands (CB1/CB2/CB3 individually resolved) plus, for spinor data, the **Dresselhaus spin splitting** $\Delta_j(k)$ of the levels around the gap in meV (zero along the [100]/[111] axes by symmetry, ~10–140 meV peaks near W/K for GaAs — directly comparable with published spin-splitting panels). High-symmetry labels for the folded MP plot are specified in FCC-primitive reduced coordinates and converted/wrapped into the cubic BZ internally.

**Unfolded k-resolved band populations.** The supercell branch index `nelec+1` mixes *different physical primitive bands* from k to k, so the SBE can instead report the population of the **physical bands closest to the gap** (VB-1, VB, CB1, CB2; spins summed) of every folded primitive BZ point. This is a three-stage pipeline — **EPM → SALMON → plotter** — that does **not** require regenerating the ground state.

**1. EPM: build the unfold map (once, cheap).** Re-diagonalize the MP-grid cubic Hamiltonians and assign every cubic band to the 4 FCC sublattices. This writes the *spectral weights* $w_s=\lvert\langle\psi\lvert P_s\rvert\psi\rangle\rvert^2$ ($\sum_s w_s=1$) of each band on each sublattice, plus the energy-ranked primitive-band index:

```sh
python3 epm_gaas_reference.py unfoldmap   # -> SYSNAME_unfold.data (next to the GS files)
```

This must be run on the **same k-grid as the ground state** (`nk`/`nb` in the file header must match the GS dataset; the SBE stops with a clear message otherwise). The GS dataset itself (`eigen.data`, `tm.data`, `k.data`) is **not** touched — the eigenvalues are unchanged, so no full `main()` rerun is needed. For the A(k,E) spectral skeleton, also generate the clean primitive dispersion once:

```sh
python3 epm_gaas_reference.py bandpath    # -> SYSNAME_bandpath.data (L-Γ-X-W-K-Γ)
```

**2. SALMON: run the SBE dynamics.** With `SYSNAME_unfold.data` present, the SBE automatically writes `SYSNAME_sbe_nex_k_unfold.data` alongside `_sbe_nex_k.data`. Per saved time it records the crystal-gauge population of the four physical levels (VB-1, VB, CB1, CB2; spins summed) at each primitive point $k_{\rm prim} = k_{\rm sc} + G_0(s)$, $s=1..4$ — i.e. 4·nk unfolded points covering the primitive FCC BZ. The population of each cubic band is **distributed over the sublattices by the spectral weights** $w_s$ (not a hard argmax), so at a symmetry degeneracy it splits *equally* among the equivalent primitive points and the result is symmetric by construction. (The optional impact-ionization channel is likewise sublattice-resolved when the map is present — the secondary pair is created in the primitive sector of the primary.)

**3. Plotter: visualise.** `plot_sbe_results.py` picks up `_sbe_nex_k_unfold.data` automatically and produces, in `./sbe_plots/`:

* `nex_k_unfold_*` — time–k maps and snapshots over the **primitive FCC BZ** (points wrapped into the first BZ). These legitimately carry satellite peaks at the zone boundary: "CB1" is the lowest conduction band of *each* valley, so sublattices 2/3/4 place their X-valley population at the X points;
* `nex_k_fold_*` — the **folded cubic-zone** view, summing the four sublattices back onto the regular cubic grid ($k_{\rm sc}=k_{\rm prim}-G_0(s)$): a single clean, hole-free zone showing the per-cubic-k total of the lowest conduction band;
* `nex_k_unfold_spectral_*` (with `--spectral`, needs `SYSNAME_bandpath.data`) — A(k,E) band-structure views, **one pair per output time** in `spectral_frames/` (path view + $k_x$ projection), with the band coloured by population and broadened by carrier kinetic energy. The colour scale is fixed across frames, so the frame sequence is a ready-to-assemble movie of the band dynamics (`ffmpeg -pattern_type glob -i 'spectral_frames/*path*.png' bands.mp4`).

```sh
python3 plot_sbe_results.py --spectral    # unfolded + folded k-maps + per-frame A(k,E)
```

In short: rerun `unfoldmap` (and `bandpath`) only, **rebuild the Fortran**, rerun the dynamics, then plot — the ground state is reused as-is.

### EPM ground-state solver (`&epm`)

The `&epm` namelist configures the local-EPM ground-state solver (`theory='epm'`):

| Parameter | Units | Default | Description |
| :--- | :--- | :--- | :--- |
| `epm_material` | — | `'GaAs'` | Material whose tabulated Cohen-Bergstresser local form factors are used (currently `'GaAs'`). |
| `epm_lattice_constant_au` | Bohr | `10.68d0` | Zincblende lattice constant $a$. |
| `epm_pw_cutoff_ry` | Ry | `11.1d0` | Plane-wave cutoff $|\mathbf{k}+\mathbf{G}|^2$ for the basis set. |

---

## Examples

### Minimal SBE Input Example

```fortran
&calculation
  theory = 'sbe'
/

&sbe
  ! ... standard SALMON SBE system parameters ...

  ! ---------------------------------------------------------
  ! 1. Kuhn-Zurek/Caldeira-Leggett Decoherence (strictly CPTP)
  ! ---------------------------------------------------------
  ! lambda = kB*T / tau_m;  enabled only when both are > 0
  sbe_decoh_temperature_k = 300.0d0
  sbe_decoh_tau_m_fs      = 10.0d0

  ! ---------------------------------------------------------
  ! 2. Frozen Core / Active Subspace Optimization
  ! ---------------------------------------------------------
  ! Only evolve bands within ±15 eV of the Fermi level non-linearly.
  ! Deep core bands will only undergo exact linear phase oscillation.
  frozen_core_threshold_ev = -15.0d0
  frozen_free_threshold_ev =  15.0d0
/
```

**Reverting to default behavior:**
* Set `sbe_decoh_temperature_k` and/or `sbe_decoh_tau_m_fs` to a non-positive value to recover the original purely-coherent (no dephasing, $D\equiv 0$, trivially CPTP) behavior.
* Set both `frozen_core_threshold_ev` and `frozen_free_threshold_ev` to `0.0d0` to force all bands into the active nonlinear subspace.

### Minimal EPM → SBE Pipeline Example

#### Standalone Python reference (`epm_gaas_reference.py`)

For quick debugging without building/running SALMON, the repository root also contains `epm_gaas_reference.py` -- a monolithic, single-machine NumPy/SciPy reimplementation of the GaAs Cohen-Bergstresser local-EPM solver (no MPI/OpenMP). It builds the same lattice/plane-wave basis/Hamiltonian/momentum matrices as `src/epm`, and writes byte-compatible `SYSNAME_k.data`/`_eigen.data`/`_tm.data` files that `gs_info_ssbe` can read directly -- so its output can be diffed against the Fortran `theory='epm'` run, or fed straight into an SBE real-time calculation. All parameters (lattice constant, plane-wave cutoff, k-grid, number of bands/electrons, sysname) are hardcoded constants at the top of the script -- including the spinor switch `INCLUDE_SPIN_ORBIT` (see the spinor pipeline example below) -- edit them there and run:

```sh
python3 epm_gaas_reference.py
```

This is a debugging aid only -- `theory='epm'` in SALMON remains the primary, MPI/OpenMP-parallel ground-state path.

```fortran
! Step 1: ground state via local EPM (writes GaAs_k/_eigen/_tm.data)
&calculation
  theory = 'epm'
/
&epm
  epm_material            = 'GaAs'
  epm_lattice_constant_au = 10.68d0
  epm_pw_cutoff_ry        = 11.1d0
/
```
```fortran
! Step 2: real-time SBE propagation reading the files generated above
&calculation
  theory = 'sbe'
/
&system
  ! sysname, lattice vectors, num_kgrid, nstate, nelec must match the EPM run
/
```

### Spinor (spin-orbit) EPM → SBE Pipeline Example

Step 1 — generate the spin-orbit split ground state with the Python reference (the spinor switch is a hardcoded constant at the top of the script):

```sh
# epm_gaas_reference.py:  INCLUDE_SPIN_ORBIT = True   (default)
python3 epm_gaas_reference.py
# writes GaAs_cubic_so_k.data / _eigen.data / _tm.data:
#   64 spin-orbit split bands (occupation 1 per band),
#   mu auto-calibrated at Gamma to Delta0 = 0.341 eV (Gamma8-Gamma7),
#   v_SO = grad_k H_SO written analytically into block 2 (rvnl_tm)
```

Step 2 — real-time SBE propagation on the spinor dataset (note `nstate` doubled, `yn_sbe_spinor` and `yn_vnl_correction` both `'y'`):

```fortran
&calculation
  theory = 'sbe'
/
&control
  sysname = 'GaAs_cubic_so'
/
&units
  unit_system = 'au'
/
&system
  yn_periodic = 'y'
  al(1:3) = 10.68d0, 10.68d0, 10.68d0   ! must match the EPM run
  nelec  = 32
  nstate = 64                            ! 2*Nb spinor bands
/
&kgrid
  num_kgrid(1:3) = 4, 4, 4               ! must match the EPM run
/
&tgrid
  dt = 0.05d0
  nt = 20000
/
&emfield
  ae_shape1 = "Acos2"
  epdir_re1(1:3) = 0.0d0, 0.0d0, 1.0d0
  I_wcm2_1 = 1.0d+11
  tw1 = 500.0d0
  omega1 = 0.056d0
/
&sbe
  yn_sbe_spinor     = 'y'   ! spinor input: occupation 1/band, nelec valence bands
  yn_vnl_correction = 'y'   ! use pi = p + v_SO from rvnl_tm everywhere
/
```

Step 3 — plot (spin pairs are summed into levels automatically):

```sh
python3 plot_sbe_results.py
```

Setting `INCLUDE_SPIN_ORBIT = False` in the script restores the scalar pipeline (`GaAs_cubic`, 32 bands, occupation 2 per band) byte-for-byte; the SBE input then keeps `yn_sbe_spinor = 'n'` (default).

### Band-structure calculation (`theory='dft_band'`)

`theory='dft_band'` diagonalizes the **converged** Kohn-Sham Hamiltonian at k-points along a high-symmetry path and writes the eigenvalues to `band.dat`. It is a post-processing step: run a normal `theory='dft'` ground state first, then restart from it. A ready-to-run pair lives in `samples/exercise_04_bulkSi_gs/` (`Si_gs.inp` + `Si_band.inp`).

```sh
cd samples/exercise_04_bulkSi_gs

# 1. Ground state (writes the restart directory data_for_restart/)
salmon < Si_gs.inp

# 2. dft_band restarts from ./restart — point it at the GS output
ln -s data_for_restart restart

# 3. Band structure along L-G-X-M-G  ->  band.dat
salmon < Si_band.inp
```

The path is given explicitly in the `&band` namelist (reduced reciprocal coordinates):

```fortran
&calculation
  theory = 'dft_band'
/
&control
  sysname    = 'Si'
  yn_restart = 'y'      ! restart from the ground-state density in ./restart
/
&band
  lattice         = 'non'              ! use the explicit kpt/ndiv_segment path below
  nref_band       = 20                 ! converge eigenvalues up to this band index
  tol_esp_diff    = 1.0d-5             ! per-band convergence tolerance on |dE| (a.u.)
  num_of_segments = 4                  ! L-G-X-M-G : 4 segments, 5 end points
  ndiv_segment(1:4) = 16, 16, 16, 16   ! k-points per segment
  kpt(1:3,1) = 0.5d0, 0.5d0, 0.5d0     ! L
  kpt(1:3,2) = 0.0d0, 0.0d0, 0.0d0     ! G
  kpt(1:3,3) = 0.5d0, 0.0d0, 0.0d0     ! X
  kpt(1:3,4) = 0.5d0, 0.5d0, 0.0d0     ! M
  kpt(1:3,5) = 0.0d0, 0.0d0, 0.0d0     ! G
  kpt_label(1) = 'L'
  kpt_label(2) = 'G'
  kpt_label(3) = 'X'
  kpt_label(4) = 'M'
  kpt_label(5) = 'G'
/
```

`band.dat` starts with a small header (`Number_of_Bands`, `Number_of_kpt_in_each_block`, `Number_of_blocks`), then one `ik  k_red(1:3)  k_cart(1:3)` line per k-point, followed by `ik  ib  energy(spin...)` eigenvalue lines (energies in Hartree). For the sample above the silicon valence-band top sits at $\Gamma$ with the conduction-band minimum near $X$ (indirect gap), as expected for an LDA silicon band structure.

`plot_sbe_results.py` plots `band.dat` directly (it is picked up automatically alongside the other band-structure files): energies are converted to eV and shifted to a valence-band-maximum reference, with vertical guides drawn at the detected path nodes (direction changes). Since `band.dat` carries no occupations, the VBM band index defaults to `nb//2` (half filling); override it with `--band-vbm IDX`.

```sh
cp plot_sbe_results.py /path/to/band_calculation/
cd /path/to/band_calculation/
python3 plot_sbe_results.py --only-bands --energy-range -13 7   # -> sbe_plots/band_dat_band.png
```

| `&band` parameter | Default | Description |
| :--- | :--- | :--- |
| `lattice` | `''` | `'non'`: take the path from `kpt`/`ndiv_segment` below. `'sc'`/`'fcc'`/`'bcc'`/`'hex'`: use a built-in default path for that Bravais lattice. |
| `nref_band` | `0` | Eigenvalues are converged (and convergence is checked) up to this band index. |
| `tol_esp_diff` | `1.0d-5` | Per-band convergence tolerance on the eigenvalue change between iterations (a.u.). |
| `num_of_segments` | `0` | Number of straight segments in the path (a path of `N` segments has `N+1` end points). |
| `ndiv_segment(:)` | `0` | Number of k-points sampled along each segment. |
| `kpt(1:3,:)` | `0` | Segment end points in **reduced reciprocal** coordinates (one more than `num_of_segments`). |
| `kpt_label(:)` | `''` | Optional labels for the end points (`'G'`, `'X'`, ...). |

---

## Building & Continuous Integration

The simplest serial (no-MPI) build, which the GitHub Actions workflow (`.github/workflows/build.yml`) also runs on every push and pull request to catch compilation errors:

```sh
cmake -B build -S . \
  -D CMAKE_BUILD_TYPE=Release \
  -D USE_MPI=OFF \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch -fallow-invalid-boz"
cmake --build build -j "$(nproc)"
# -> ./build/salmon
```

The `-fallow-argument-mismatch -fallow-invalid-boz` flags are needed for the bundled serial-communication stub under modern gfortran (≥ 10). For production runs configure with `-D USE_MPI=ON` and an MPI Fortran compiler, or use `configure.py` as in the upstream SALMON documentation.

---

## References & Theoretical Background

1. **Commutator-Free Magnus Integrators:** Blanes, S., & Moan, P. C. "Practical symplectic partitioned Runge-Kutta and Runge-Kutta-Nyström methods." *J. Comput. Appl. Math.* 142, 313-330 (2002); Alvermann, A., & Fehske, H. "High-order commutator-free exponential time-propagation of driven quantum systems." *J. Comput. Phys.* 230, 5930-5956 (2011).
2. **Suzuki-Yoshida Composition:** Yoshida, H. "Construction of higher order symplectic integrators." *Phys. Lett. A* 150, 262-268 (1990).
3. **CPTP / Lindblad & RBF-kernel positivity:** Schoenberg, I. J. "Metric spaces and completely monotone functions." *Ann. Math.* 39, 811-841 (1938) (Bochner/Schoenberg PSD criterion for Gaussian/RBF kernels); Schur product theorem (Hadamard maps of PSD matrices are PSD).
4. **Caldeira-Leggett / Kuhn-Zurek Decoherence:** Caldeira, A. O., & Leggett, A. J. "Path integral approach to quantum Brownian motion." *Physica A* 121, 587-616 (1983); Zurek, W. H. "Decoherence, einselection, and the quantum origins of the classical." *Rev. Mod. Phys.* 75, 715 (2003).
5. **Cohen-Bergstresser Local Pseudopotentials:** Cohen, M. L., & Bergstresser, T. K. "Band Structures and Pseudopotential Form Factors for Fourteen Semiconductors of the Diamond and Zinc-blende Structures." *Phys. Rev.* 141, 789 (1966).
6. **Spin-Orbit in EPM:** Weisz, G. "Band Structure and Fermi Surface of White Tin." *Phys. Rev.* 149, 504 (1966); Bloom, S., & Bergstresser, T. K. "Band structure of α-Sn, InSb and CdTe including spin-orbit effects." *Solid State Commun.* 6, 465 (1968); Chelikowsky, J. R., & Cohen, M. L. "Nonlocal pseudopotential calculations for the electronic structure of eleven diamond and zinc-blende semiconductors." *Phys. Rev. B* 14, 556 (1976).
7. **Velocity-Gauge SBE / Houston Basis:** Wismer, M. S., & Yakovlev, V. S. "Gauge-independent decoherence models for solids in external fields." *Phys. Rev. B* 97, 144302 (2018).
8. **Original SALMON SBE:** Sato, S. A. et al. "Multiscale computational approach for light-matter interactions." *Phys. Rev. B* 92, 115145 (2015).

## License

SALMON is available under the Apache License version 2.0.

    Copyright 2017-2026 SALMON developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.