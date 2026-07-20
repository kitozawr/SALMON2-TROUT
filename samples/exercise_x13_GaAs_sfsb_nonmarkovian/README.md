# exercise x13 — SFSB non-Markovian memory-kernel ionization (GaAs Γ–L)

**What this teaches:** replacing the relaxation-time approximation (a constant
dephasing time T₂) with a **non-Markovian heat bath** — the strong-field
spin-boson (SFSB) model of **Boroumand et al., Rep. Prog. Phys. 88, 070501
(2025)** [B25], transcribed and implemented in
[`wiki/10`](../../wiki/10_open_quantum_systems_literature.md) §6–7.

The entire environment enters through ONE object, the bath correlation
function C(t₁−t₂) [B25 Eq. (5)]; the conduction population is the double-time
**memory integral** [B25 Eq. (3)]

$$n_c(K,t) = \tfrac12\,\mathrm{Re}\!\int^t\! dt_1\!\int^{t_1}\! dt_2\,
\Omega^*(K_{t_1},t_1)\,\Omega(K_{t_2},t_2)\,
e^{\,iS(t_1,t_2) + C(t_1-t_2)}$$

The kernel $e^{C(\tau)}$ does not factorize across time steps — that
non-factorizable factor IS the memory. The Markovian relaxation-time
approximation is the special case C = −τ/T₂.

## Run

```bash
../../build/salmon < GaAs_line_epm_gs.inp     # 1D Γ–L line ground state (~seconds)
../../build/salmon < GaAs_sfsb_rt.inp         # SFSB run (~1 min, OpenMP over k)
```

Outputs: `GaAs_line_sfsb_nex.data` (nc(t), same normalization as
`_sbe_nex.data`) and `GaAs_line_sfsb_nck.data` (k-resolved blocks).

## The physics story (validated numbers, 192-pt line, γ≈2.15, 3.5-photon)

| bath | nex(t_end) [cm⁻³] | η | B25 claim reproduced |
|---|---|---|---|
| none (coherent) | 3.0e18 | 1 | baseline; k-grid converged 96/192/384 to 0.2 %/0.005 % |
| RTA T₂ = 6 fs | 9.1e19 | **×30** | **dephasing ionization** — the RTA converts returning *virtual* polarization into *real* carriers [Fig 1] |
| ohmic j₀=1, ω_c=2.1ω₀, 300 K | 1.2e17 | **×25 suppressed** | dephasing-**suppressed** ionization at low T / strong coupling [Fig 3] |
| … same, `yn_sbe_bath_imc='n'` | 4.7e20 | ×158 | killing Im C flips suppression → enhancement: **the bath phase (dynamic gap addition) is the suppressor** [Fig 3(c)] |
| … same bath, 20 000 K | 3.0e20 | ×101 | enhancement only at extreme temperatures [Fig 3] |
| ohmic/debye j₀=0.1, ω_c=0.1ω₀, T ≤ 3000 K | ≈3.2e18 | ≈1.08 | physical baths barely change ionization at realistic T [Fig 2(a)] |
| stronger field (E₀=10⁹ V/m, γ≈0.43) | — | RTA ×30→×1.4 | bath influence shrinks in the tunneling regime [Fig 5(b)] |

Full figure: [`sfsb_validation.png`](sfsb_validation.png) (ships with this exercise).

## Validity & gotchas (honest scope — wiki/10 §7)

- **2nd order in the drive**: small ionized fraction (multiphoton, γ > 1);
  no depletion feedback. Not a replacement for the density-matrix SBE — a
  physically-correct **reference bracket** for how much ionization a *real*
  environment can add or remove.
- **1D line**: `num_kgrid = (N,1,1)` and **E ∥ b1** (the code error-stops
  otherwise). B25 use the same 1D reduction and verify it against 3D.
- **Two-band reduction**: gap = band-edge E_c−E_v,top (continuous by
  eigenvalue sorting); coupling = quadrature (bright-state) sum over the top
  `sbe_sfsb_nv` valence bands — individual |d_vc| of degenerate members are
  gauge-random per k, only the quadrature sum is invariant.
- **Selection-rule walls**: if the channel is strictly forbidden over part of
  the line (e.g. CdS Γ–M with E⊥c: allowed pocket |q₁| ≲ 0.1 then exactly 0),
  no smooth two-band line exists — the code prints a WARNING and the result
  is not k-converged. GaAs Γ–L is clean; that is why this exercise uses it.
