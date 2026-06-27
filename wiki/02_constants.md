# Constants & Coefficients

> Every default numerical constant with its primary source. If you change a default, update the citation here AND in the code comment. ✅ = value live in code.

## 0. Material registry — the single per-material source

All per-material constants that the SBE **dissipation channels** consume are assembled in one place: `get_material_params(name)` in [`../src/ssbe/sbe_superres_ssbe.f90`](../src/ssbe/sbe_superres_ssbe.f90), returning an `s_material_params` struct (dielectric ε₀/ε∞, impact-ionization fit form/exponent/prefactor/threshold, electron-phonon phonon table + ν_sat, lattice constant, diamond flag). The channels auto-select through it; the `&sbe` namelist defaults are **sentinels** (`sbe_ii_form='auto'`, `sbe_ii_exponent/prefactor/threshold ≤ 0`, `sbe_coulomb_epsilon ≤ 0`) that resolve to the material value, and any explicit value overrides. GaAs resolves to the exact legacy numbers.

**To add a material** (e.g. Ge, GaN, …):
1. Add the cited raw constants (dielectric, phonon table, II fit) to the parameter block at the top of `sbe_superres_ssbe.f90`, *and a row in this page*.
2. Add one `case ('<name>')` block to `get_material_params` filling the struct from those constants.
3. Add the name to `MAT_SUPPORTED`.
4. (If the material needs EPM bands) add its form factors to `epm_cohen_bergstresser.f90` (also case-based) and to `epm_gaas_reference.py`.
5. Add the expected struct values to `tests/test_material_registry.f90`.

No channel code changes — every dissipation channel reads the struct. A material-dependent channel requested for a name not in the registry stops with `error stop` and the supported list.

### Provenance rule (STRICT) — no source ⇒ forbidden ✅

**A channel may be enabled for a material ONLY if its constants are backed by a cited source for THAT material. No source ⇒ the constant is invalid ⇒ the channel is forbidden, and the SBE init aborts (`error stop`). Constants are NEVER transferred from another material.** The registry carries per-channel provenance gates `ii_ok / eph_ok / eeh_ok / coulomb_ok` (default `.false.`); a material sets a gate `.true.` only when it supplies the cited constants.

### Effect-support matrix (which modes are allowed per material) ✅

| Effect | GaAs | Si | CdS | graphene | CdS source |
|---|---|---|---|---|---|
| EPM band structure | ✅ CB1966 | ✅ Kunikiyo 1994 | ✅ **validated** (2.55 vs 2.58 eV) | ✅ **validated** (Dirac cone, v_F=9.6e5) | BC1967; exact 2-fold folding to the al-vector SBE cell |
| Kuhn-Zurek decoherence | ✅ | ✅ | ✅ | ⛔ **disabled by policy** (G6) | material-independent (λ=k_B T/τ_m, user-supplied) |
| Coulomb HF | ✅ ε=12.9 | ✅ ε=11.7 | ✅ ε=8.9 | 🚧 (2D V(q)=2πe²/εq) | static dielectric (isotropic avg) [Berlincourt 1963] |
| Electron-phonon | ✅ Fischetti 1e14 | ✅ Meng 1.3e14 | ✅ Fröhlich LO 38 meV, ν_sat=2.9e13 | 🚧 E2g/A1'/acoustic (G4) | ħω_LO [Raman], α≈0.5 [cyclotron]; ν_sat = α·ω_LO [md] |
| Impact ionization | ✅ Stobbe quartic | ✅ Keldysh quadratic | ✅ E_th=3.6 eV; **prefactor = user fit param** | 🚧 Auger (gapless) | E_th=1.5·E_g [md (3/2)E_g]; prefactor scarce → must set `sbe_ii_prefactor` |
| Carrier-carrier (e-e/e-h) | ✅ | ✅ | ✅ | 🚧 (G5) | **sub-100fs thermalization at n > 1e18, Shah 1986; Elsaesser 1989** |
| Auger recombination | 🚧 | 🚧 | ✅ | 🚧 (gapless CM) | **Haury 1998 (C = 2.0e-30 cm⁶/s)** |

CdS dissipation: Fröhlich polar-optical e-ph is the **primary** room-T channel. The impact-ionization **prefactor** is a fit parameter (no cited CdS value) — enabling II without an explicit `sbe_ii_prefactor` aborts. **Carrier-carrier (e-e) and Auger recombination are now cited and enabled** — they switch on **dynamically when the carrier density exceeds the activation thresholds** (e-e at n ≳ 1×10¹⁸ cm⁻³ → sub-100 fs thermalization [Shah 1986; Elsaesser 1989]; Auger via C·n³ with C = 2.0×10⁻³⁰ cm⁶/s [Haury 1998]). The **piezoelectric acoustic** (e₃₃/e₃₁/e₁₅ [Berlincourt 1963]) and **deformation-potential acoustic** channels are cited but not yet SBE Lindblad channels (future). ⛔ = enabling it for that material aborts; nothing is borrowed between materials. 🚧 = the constants are cited and documented but the SBE channel/material wiring is in progress.

### What the code selects per material (the complete `s_material_params`) ✅

These are **exactly** the fields `get_material_params(name)` fills — i.e. every constant a channel uses when you set `epm_material` and leave the `&sbe` knobs at their `'auto'`/sentinel defaults. Deeper provenance for each number is in the per-effect sections below.

**`epm_material = 'GaAs'`** (zincblende, polar)
| Struct field | Value | Used by | Source |
|---|---|---|---|
| `a_lattice_au` | 10.68 Bohr (5.65 Å) | reference | std GaAs |
| `is_diamond` | `.false.` (V^A≠0) | EPM structure factor | — |
| `eps0` / `eps_inf` | 12.9 / 10.89 | Coulomb HF, screening | std GaAs; NSM/Ioffe |
| `ii_form` | `stobbe_quartic` | impact ionization | Stobbe-Redmer-Schattke, PRB 49, 4494 (1994) |
| `ii_exponent` | 4 | impact ionization | same |
| `ii_prefactor` | 2×10¹² s⁻¹eV⁻⁴ | impact ionization | same (Eq. 11) |
| `ii_threshold_ev` | 2.1 | impact ionization | same |
| `eph_nu_sat_si` | 1.0×10¹⁴ s⁻¹ | e-ph rate cap | Fischetti, IEEE TED 38, 634 (1991) |
| `eph_polar` | `.true.` (Fröhlich LO) | e-ph | — |
| `eph_nph` | 6 | e-ph | — |
| `eph_hw_mev(1:6)` | 36 (LO), 27.8, 29.9, 29.0, 29.3, 29.9 | e-ph phonon energies | LO: std GaAs/Adachi; IV: Fischetti-Laux, PRB 38, 9721 (1988) |
| `eph_wraw(1)` | Σ of the 5 IV weights (polar-LO dominant ≈50% after norm) | e-ph weight | project convention |
| `eph_wraw(2:6)` | D²/ħω with D[eV/Å]=10,10,10,5,7 (Γ→L,Γ→X,L→L,L→X,X→X) | e-ph weight | Fischetti-Laux 1988 |

**`epm_material = 'Si'` or `'Si_cb'`** (diamond, non-polar)
| Struct field | Value | Used by | Source |
|---|---|---|---|
| `a_lattice_au` | 10.26 Bohr (5.431 Å) | reference | std Si; Kunikiyo JAP 75, 297 (1994) |
| `is_diamond` | `.true.` (V^A≡0) | EPM structure factor | diamond symmetry |
| `eps0` / `eps_inf` | 11.7 / 11.7 (non-polar ⇒ equal) | Coulomb HF, screening | std Si |
| `ii_form` | `keldysh_quadratic` | impact ionization | Keldysh, JETP 21, 1135 (1965); Cartier, APL 62, 3339 (1993) |
| `ii_exponent` | 2 (set 4.6 for full-band) | impact ionization | Cartier 1993; full-band Kamakura, JAP 75, 3500 (1994) |
| `ii_prefactor` | 2×10¹² s⁻¹eV⁻² | impact ionization | project (same scale as GaAs; tune to 1e13–1e14 @ 3–4 eV) |
| `ii_threshold_ev` | 1.1 | impact ionization | near the 1.12 eV gap (NOT 3/2·E_g); Cartier 1993 |
| `eph_nu_sat_si` | 1.3×10¹⁴ s⁻¹ | e-ph rate cap | Meng, PRB 91, 075201 (2015); Fischetti-Laux 1988 |
| `eph_polar` | `.false.` (no Fröhlich LO) | e-ph | — |
| `eph_nph` | 6 (intervalley g/f) | e-ph | — |
| `eph_hw_mev(1:6)` | 10, 19, 63, 19, 51, 57 (g-TA,g-LA,g-LO,f-TA,f-LA,f-TO) | e-ph phonon energies | Jacoboni-Reggiani, RMP 55, 645 (1983); Pop set |
| `eph_wraw(1:6)` | D²/ħω with D[10⁸eV/cm]=0.3,1.5,6.0,0.5,3.5,1.5 | e-ph weight | Pop set (Jacoboni-Lugli); Canali, PRB 15, 3994 (1977) |

Notes: `eph_wraw` is the **un-normalized** D²/ħω weight; the channel normalizes Σ=1 at runtime, so the GaAs (eV/Å) vs Si (10⁸ eV/cm) deformation-potential units only need within-material consistency. The absolute e-ph magnitude is set by `eph_nu_sat_si`. `ii_prefactor`/`ii_exponent` carry matching units (s⁻¹eV⁻ᵃ).

**`epm_material = 'CdS'`** (wurtzite P6₃mc, polar, non-centrosymmetric)
| Struct field | Value | Used by | Source |
|---|---|---|---|
| `cell_au` (al box) | (7.816, 13.538, 12.685) Bohr = (a, a√3, c) | EPM geometry / `&system al(1:3)` | a=4.136 Å, c/a=1.623; BC1967 Table I |
| `is_diamond` | `.false.` (V^A≠0, broken inversion) | EPM structure factor | wurtzite, u=3/8 |
| `eps0 / eps_inf` | 8.9 / 5.3 | Coulomb HF / screening | static/∞ dielectric (isotropic avg) [Berlincourt 1963] |
| `eph_polar / eph_nph / eph_hw_mev(1)` | `.true.` / 1 / 38 meV | e-ph | Fröhlich LO [Raman; md] |
| `eph_nu_sat_si` | 2.89e13 s⁻¹ | e-ph rate scale | α·ω_LO from α=0.5, ħω_LO=38 meV [md] |
| `ii_form / ii_exponent / ii_threshold_ev` | keldysh_quadratic / 2 / 3.6 eV | impact ionization | Keldysh soft; E_th=1.5·E_g [md] |
| `ii_prefactor` | **sentinel (−1)** → user must set `sbe_ii_prefactor` | impact ionization | no cited CdS prefactor (fit parameter) [md] |
| `auger_coeff` | 2.0e-30 cm⁶/s | Auger recombination | Haury et al., PRB 57, 11513 (1998) |
| `ee_activation_n` | 1.0e18 cm⁻³ | e-e scattering threshold | Shah et al., IEEE JQE 22, 1728 (1986) |
| `coulomb_ok / eph_ok / ii_ok / eeh_ok` | `.true.` | provenance gates | all cited |

Every CdS constant carries its source. The impact-ionization **prefactor** is the one genuinely uncited quantity (a fit parameter), so it stays a sentinel and the run aborts unless the user supplies `sbe_ii_prefactor`. **Carrier-carrier (e-e) and Auger are now cited and enabled**, switching on dynamically above their density thresholds (`ee_activation_n`; Auger via C·n³). Piezoelectric (e₃₃/e₃₁/e₁₅ [Berlincourt 1963]) and deformation-acoustic are cited but not yet SBE channels.

## 1. EPM form factors (local pseudopotential)

### GaAs (zincblende, Rydberg, |G|² in (2π/a)² units) ✅
| Quantity | Value (Ry) | Source |
|---|---|---|
| V^S(3) | −0.23 | Cohen & Bergstresser, Phys. Rev. 141, 789 (1966) |
| V^S(8) | +0.01 | CB 1966 |
| V^S(11) | +0.06 | CB 1966 |
| V^A(3) | +0.07 | CB 1966 |
| V^A(4) | +0.05 | CB 1966 |
| V^A(11) | +0.01 | CB 1966 |
| a | 5.65 Å = 10.68 Bohr | std GaAs; CB 1966 |

### Silicon (diamond, Rydberg) — DEFAULT Kunikiyo ✅
| Quantity | Kunikiyo (default) | Cohen-Bergstresser (alt `Si_cb`) | Source |
|---|---|---|---|
| V^S(3) | −0.2258 | −0.21 | Kunikiyo et al., JAP 75, 297 (1994) Table I; CB 1966 |
| V^S(8) | +0.05698 | +0.04 | same |
| V^S(11) | +0.070709 | +0.08 | same |
| V^A(all) | 0 (exact) | 0 (exact) | diamond: two identical atoms → antisymmetric structure factor vanishes |
| a | 5.431 Å = 10.26 Bohr | same | std Si; Kunikiyo 1994 |

**Validation (this fork):** converged indirect gap 1.059 eV (Kunikiyo's own calc 1.068; exp 1.12), CBM at 0.86·2π/a along ⟨100⟩. The 3-parameter local EPM intentionally does not reach the experimental gap. **Structure factor:** τ=(a/8)(1,1,1) for both; diamond purely real cos(G·τ); zincblende cos + i·sin.

### CdS (wurtzite P6₃mc, **LOCAL** EPM) ✅ validated against BC1967
**Source: T. K. Bergstresser & M. L. Cohen, *Phys. Rev.* 164, 1069 (1967).** The paper assumes spherically-symmetric atomic potentials (NO angular/nonlocal term) → the pseudopotential is **purely local** (so `rvnl_tm = 0`, like local GaAs/Si; there is **no** cited CdS nonlocal parameter, so none is added). BC1967 give the **wurtzite** form factors directly in Table II (their tuned values); the cited values (Ry, keyed by |G|² in the reduced units (√2π/a_W)² = (2π/a_ZB)², a_ZB=√2·a_W):

| shell (\|G\|²) | V^S | V^A | Source |
|---|---|---|---|
| 002 (3.04) | −0.26 | +0.23 | BC1967 Table II (wurtzite CdS) |
| 101 (3.43) | −0.24 | +0.18 | same |
| 102 (5.70) | −0.20 | +0.08 | same |
| 103 (9.50) | +0.04 | +0.05 | same |
| 200/201 (10.7/11.4) | +0.04 | +0.05 | same |
| 202 (13.3) | +0.02 | +0.03 | same |
| a, c/a, u | 4.136 Å, 1.623, 3/8 | | BC1967 Table I |
| gap | **2.58 eV**, direct at Γ | | BC1967 Table I (validation target) |

**Validation ✅:** implemented in [`../epm_wurtzite_cds.py`](../epm_wurtzite_cds.py) (hexagonal primitive cell for the band check; orthorhombic √3×1×1 cell from `al(1:3)` for the SBE). The **structure factors match Table II** (002: |S^S|,|S^A|=0.71,0.71; 101: 0.33,0.80; 102: 0.35,0.35) and the **converged direct gap at Γ = 2.55 eV reproduces the paper's 2.58 eV** (|Δ|≈0.03 eV, inside the paper's ~0.27 eV form-factor accuracy). Two things had to be right (both were initially wrong, giving 13 eV): the potential is normalized by **total atoms per cell** (1/n, the BC1967 "volume per atom" normalization — *not* per species), and the **wurtzite** Table II form factors are used (not the zinc-blende anchors).

**Folding to the SBE cell (EXACT) ✅:** the SBE uses the orthorhombic 8-atom `al(1:3)` cell, a √3×1×1 supercell of the hexagonal primitive cell (folding factor 2). Because the supercell potential is primitive-periodic, the supercell Hamiltonian at Γ is **block-diagonal over the 2 cosets** of the primitive reciprocal lattice (off-coset |H| ≈ 8×10⁻¹⁷ — the analogue of the cubic 4-fold FCC folding's parity selection rule). Coset 0 = Γ_hex carries the 2.54 eV direct gap; coset 1 = the zone-edge partner (6.2 eV, so the CBM stays at Γ). The orthorhombic gap equals the primitive gap → the folding is consistent and exactly unfoldable. Unit-tested in [`../tests/test_wurtzite_cds_epm.py`](../tests/test_wurtzite_cds_epm.py).

## 2. Spin-orbit (GaAs spinor mode) ✅
| Quantity | Value | Source |
|---|---|---|
| Δ₀ (Γ8-Γ7 split-off) | 0.341 eV (calibration target) | Vurgaftman-Meyer-Ram-Mohan, JAP 89, 5815 (2001) |
| E_g(GaAs) direct, 0 K | 1.519 eV | Blakemore, JAP 53, R123 (1982) |
| SO operator | projected Weisz/Bloom-Bergstresser | Weisz PR 149, 504 (1966); Bloom-Bergstresser SSC 6, 465 (1968); Chelikowsky-Cohen PRB 14, 556 (1976) |

## 3. Impact-ionization fits: γ = P(ε−E_th)^a Θ(ε−E_th) ✅ (form switch = Part B)
| Material | a | E_th (eV) | P | Threshold | Source |
|---|---|---|---|---|---|
| GaAs (default) | 4 | 2.1 | 2e12 s⁻¹eV⁻⁴ | hard | Stobbe-Redmer-Schattke, PRB 49, 4494 (1994) Eq. 11 |
| Si (default) | 2 | 1.1 | tuned ~1e13–1e14 at 3–4 eV | soft | Keldysh, JETP 21, 1135 (1965); Cartier et al., APL 62, 3339 (1993) |
| Si (full-band opt.) | 4.6 | 1.15 | Kamakura fit | soft | Kamakura et al., JAP 75, 3500 (1994) |

Si e⁻ threshold 1.14 eV, hole 1.37 eV; near the 1.12 eV gap, NOT 3/2·Eg. Ramp σ_E smooths Θ.

## 4. Silicon electron-phonon deformation potentials ✅ (Part C5)
Six intervalley phonons (Pop "new" set default; D in 1e8 eV/cm, E in meV). [Jacoboni & Reggiani, RMP 55, 645 (1983); Pop set in Jacoboni-Lugli; anchor Canali et al., PRB 15, 3994 (1977)]
| Phonon | type | E (meV) | D (Pop, default) | D (J-R alt) |
|---|---|---|---|---|
| g-TA | g | 10 | 0.3 | 0.5 |
| g-LA | g | 19 | 1.5 | 0.8 |
| g-LO | g | 63 | 6.0 | 11 |
| f-TA | f | 19 | 0.5 | 0.3 |
| f-LA | f | 51 | 3.5 | 2.0 |
| f-TO | f | 57 | 1.5 | 2.0 |

g-type couple same-⟨100⟩-axis valleys; f-type orthogonal-axis. **Project decision: empirical values, NOT ab-initio** (ab-initio Cheng et al., PRB 104, 195201 (2021) is substantially lower — reference only).

| Other Si quantity | Value | Source |
|---|---|---|
| Acoustic def. potential Ξ_d | ~9 eV | Jacoboni-Reggiani 1983 |
| D_LA / D_TA | 6.39 / 3.01 eV | J-R 1983 |
| Mass density ρ | 2.33 g/cm³ | std Si |
| v_LA / v_TA | 9.01e5 / 5.23e5 cm/s | std Si; J-R 1983 |
| Non-parabolicity α | 0.5 eV⁻¹ | J-R 1983 |
| Δ-valley position | 0.85·2π/a along ⟨100⟩ | Kunikiyo 1994 |
| m_l / m_t | 0.916 / 0.19 mₑ | std Si |

## 5. GaAs electron-phonon (polar: Fröhlich + intervalley) ✅ (Part C5)
### Fröhlich polar-optical
| Quantity | Value | Source |
|---|---|---|
| ħω_LO | 36 meV | std GaAs; Adachi |
| Fröhlich α | 0.068 | Adachi |
| ε₀ / ε_inf | 12.9 / 10.89 | std GaAs |
| m*_Γ | 0.067 mₑ | std GaAs |
| Rate form | (1/√E)·asinh(√(E/ħω₀−1)) | Fawcett-Boardman-Swain, JPCS 31, 1963 (1970) |

### GaAs intervalley deformation potentials (D in eV/Å, E in meV)
| Process | D | E | Source |
|---|---|---|---|
| Γ→L | 10 | 27.8 | Fischetti-Laux, PRB 38, 9721 (1988); IEEE TED 38, 634 (1991) |
| Γ→X | 10 | 29.9 | same |
| L→L | 10 | 29 | same |
| L→X | 5 | 29.3 | same |
| X→X | 7 | 29.9 | same |
| Acoustic Ξ_dΓ | 7.01 eV | — | Fischetti-Laux |

### GaAs valley structure
| Quantity | Value | Source |
|---|---|---|
| Γ-L separation | 0.29 eV | Adachi; Jiang & Wu 0.296 |
| Γ-X separation | 0.48 eV | Fischetti-Laux |
| Equivalent L / X valleys | 4 / 3 | std |
| L-valley DOS mass (per valley, in Lindblad) | 0.22 mₑ | Adachi; project decision (per-valley, not 0.55 averaged) |

## 6. Collision-rate saturation ν(ε) = ν_sat[1 − exp(−(ε/ε₀)^n)] ✅
| Quantity | Value | Source |
|---|---|---|
| ν_sat (Si) | 1.3e14 s⁻¹ | project (Fischetti priority); Meng et al., PRB 91, 075201 (2015); Fischetti-Laux 1988 |
| ν_sat (GaAs) | 1e14 s⁻¹ | Fischetti, IEEE TED 38, 634 (1991) (ab-initio Bernardi PNAS 112, 5291 ~1e14) |
| n | 2 (start) | parametrization |
| ε₀ | ~0.8 eV | parametrization |

**Never a hard min(ν, ν_sat)** — derivative discontinuity destabilizes the stiff solver.

## 7. Bandgap renormalization (II threshold) ✅ (Part C7)
| Quantity | Value | Source |
|---|---|---|
| EHP cube-root law | ΔE_gap[eV] = −1.9e-8 (n[cm⁻³])^(1/3) | Vashishta & Kalia, PRB 25, 6492 (1982) |
| at 1e18 / 1e19 | −19 / −41 meV | same |
| Gate density | 5e18 cm⁻³ | project (below = within fit uncertainty) |
| K ambiguity | [1.9, 3.8]e-8 eV·cm | factor-2 convention |
| Dopant BGN (distinct) | 22.5(n/1e18)^(1/2) meV | Lanyon & Tuft, IEEE TED ED-26, 1014 (1979) |

## 8. Coulomb / Hartree-Fock ✅
| Quantity | Value/form | Source |
|---|---|---|
| Σ^HF_nm(k) | −Σ_{q≠k} V(k−q) δρ_nm(q) | Golde-Kira-Meier-Koch, PSS B 248, 863 (2011) |
| V(p) | strength·4π / (ε Ω_cell N_k (|p|²+κ²)) | same |
| ε (GaAs) | 12.9 | std GaAs |
| screening κ | `sbe_coulomb_screen_au` (0 = bare, q=0 excluded) | regularization |
| sublattice projection | block-diagonal over 4 FCC cosets | Popescu & Zunger, PRB 85, 085201 (2012); Ku-Berlijn-Lee, PRL 104, 216401 (2010) |

## 9. Integrator constants ✅
| Quantity | Value | Source |
|---|---|---|
| CF4 nodes c1,c2 | ½ ∓ √3/6 | Blanes & Moan, JCAM 142, 313 (2002); Alvermann & Fehske, JCP 230, 5930 (2011) |
| CF4 weights α1,α2 | ¼ ± √3/6 | same |
| Yoshida p1 | 1/(2−2^(1/3)) ≈ 1.35120719196 | Yoshida, PLA 150, 262 (1990) |
| Yoshida p2 | −2^(1/3)/(2−2^(1/3)) ≈ −1.70241438392 | Yoshida 1990 |

## 10. Energy bins / search (super-mode) ✅ (Part C3)
| Quantity | Value/rule | Source |
|---|---|---|
| δ(ΔE) replacement | normalized Gaussian or unit-area rectangle | Stobbe 1994 (0.2 eV rect); Kunikiyo 1994 (5 meV bins) |
| σ_E | ~ mean inter-level spacing of the grid | grid-matched |
| pair search | deterministic / energy-windowed; NO Monte-Carlo | project (MC breaks CPTP) |
| sub-cycling trigger | ν_max(h/2) ≳ 0.2 → m sub-steps, (τ/2m)ν_max ≲ 0.1 | Lindblad 1976 |

## 11. Ring/pipeline MPI (super-mode) ✅ (Part D)
| Quantity | Value | Source |
|---|---|---|
| memory per rank | O(N_k/P + one transit block) | Plimpton, JCP 117, 1 (1995) |
| communication | O(P) sends of size O(N_k/P) | same |
| force-decomposition upgrade | N/√P | Plimpton & Hendrickson, JCC 17, 326 (1996) |

---

## 12. Dielectric screening (Part G ✅) and carrier-carrier (Part F)
| Quantity | Value | Source |
|---|---|---|
| GaAs ε₀ / ε_∞ | 12.9 / 10.89 | std GaAs; NSM/Ioffe |
| GaAs ħω_LO / ħω_TO | 36 / 33.6 meV | std GaAs (300 K NSM 36.1/33.2) |
| GaAs m*_Γ / m_hh / m_lh | 0.067 / 0.51 / 0.082 mₑ | std GaAs |
| GaAs ω_p (n=1e18 / 1e19) | ≈43.5 / 137 meV (uses ε_∞) | Mahan; Mooradian-McWhorter PR 177, 1231 |
| GaAs LOPC anticrossing (ω_p=ω_LO) | n ≈ 7×10¹⁷ cm⁻³ | Mooradian-McWhorter PR 177, 1231 (1969) |
| Si ε (static, non-polar) | 11.7 | std Si |
| Si m_l / m_t / m_hh / m_lh | 0.98 / 0.19 / 0.49 / 0.16 mₑ | Si effective-mass tables |
| Si LOPC | N/A (non-polar, no Fröhlich LO-plasmon) | Fischetti PRB 44, 5527 (1991) |
| **CdS ε₀ / ε_∞ (isotropic avg)** | **8.9 / 5.3** | Berlincourt et al., PR 129, 1009 (1963) |
| **CdS m*_e / m*_h** | **0.21 / 0.70 mₑ** | Chiang & Himpsel (per arXiv:1011.2793) |
| **CdS Auger coeff (n-type)** | **2.0e-30 cm⁶/s** | Haury et al., PRB 57, 11513 (1998) |
| **CdS e-e thermalization** | **< 100 fs at n ≥ 1e18 cm⁻³** | Shah et al., IEEE JQE 22, 1728 (1986); Elsaesser et al. PRL 66, 1757 (1991) |
| **CdS Thomas-Fermi screening** | **dynamic q_TF(n, T_e), ε_∞ = 5.3** (fast e-e process) | Hwang-Das Sarma / Lindhard |
| κ_TF² (degenerate, a.u.) | 4(3n/π)^⅓/ε = 4k_F/(πε) | Ashcroft-Mermin; arXiv:2312.13059 |
| κ_D² (Debye, a.u.) | 4πn/(εk_BT) | Ashcroft-Mermin |
| Lindhard F(x) | ½+(1−x²)/(4x)ln\|(1+x)/(1−x)\| | Lindhard, Mat.-Fys. Medd. 28, 8 (1954) |
| carrier-carrier rate scale | 1e13–1e14 s⁻¹ @ 1e17–1e19 cm⁻³ | Goodnick-Lugli PRB 37, 2578; Fischetti-Laux PRB 38, 9721 |
| EID slope (qualitative) | γ=γ₀+a·n, linear; e-carrier ≈8× e-exciton | Honold et al. PRB 40, 6442 (1989) |
| screening default | static Lindhard/RPA (b); LOPC (c) GaAs-only n≳5e17 | recommendation |

**No HF double-counting:** carrier-carrier (F) = correlation (2nd-Born/GW) self-energy, dissipative only; the static screened-exchange shift stays in Σ^HF. Carrier-carrier conserves Σf_k and ΣE_k f_k (validation invariants).

## 13. Lindblad jump operators (Auger / impact-ionization-recombination) 🚧

Auger recombination and its inverse (impact ionization / carrier multiplication) are added as a **density-gated, number-conserving CPTP channel** built from GKLS jump operators on the active level set — the same Taj-Rossi → Kossakowski machinery used for the other scattering channels. Reusable across materials (written for graphene's gapless CM, reused for CdS).

**Process & operators.** Auger recombination annihilates an electron–hole pair and promotes a third carrier; impact ionization is the inverse. On a level set with adiabatic populations `ρ̃`, the jump operators are the products `L_recomb ∝ a†_h' a_e1 a_e2` (two carriers → one hot carrier + recombination) and `L_ii ∝ a†_e1 a†_e2 a_h'` (one hot carrier → two), cast into amplitude-transfer maps that **conserve total carrier number** (Auger is number-conserving within the e+h subsystem) and apply Pauli blocking (final-state `(1−f)` factors clamped to [0,1]). Each finite step `exp(τD)` is a genuine GKLS map (τ>0).

**Rate / gating.** The Auger rate per unit volume is `R_Auger = C·n³` (electron-density form) with the cited **C = 2.0×10⁻³⁰ cm⁶/s for CdS** [Haury 1998]; the channel is **density-gated** — it is inert below `ee_activation_n` (CdS 1×10¹⁸ cm⁻³ [Shah 1986]) and turns on as the photoexcited density crosses it. For **graphene** the process is **nearly thresholdless** (gaplessness → carrier multiplication, CM factor up to ~2 [Plötzing 2014]); the collinear-scattering phase space is regularized by the lifetime broadening generated self-consistently by the e-ph + e-e rates (not an ad-hoc width). Screening choice (`static`|`dynamical` RPA) is a modeling switch — static maximizes CM/AR, dynamical can suppress it [Tomadin et al., PRB 88, 035430 (2013)].

**CPTP invariants (tested where implemented):** trace (carrier number) conserved; populations stay in [0, occ]; Hermiticity preserved; γ=0 → identity. Sources: Taj-Rossi PRA 78, 052113 (2008); Rosati et al. PRB 90, 125140 (2014); Auger coeff Haury et al. PRB 57, 11513 (1998).
