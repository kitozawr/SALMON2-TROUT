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

### Registry values (live) ✅
| Material | ε₀ | ε∞ | II form / a / E_th[eV] | e-ph ν_sat[s⁻¹] | phonons | a[Bohr] |
|---|---|---|---|---|---|---|
| GaAs | 12.9 | 10.89 | stobbe_quartic / 4 / 2.1 | 1.0e14 | Fröhlich-LO 36 meV + 5 IV | 10.68 |
| Si / Si_cb | 11.7 | 11.7 | keldysh_quadratic / 2 / 1.1 | 1.3e14 | 6 intervalley g/f | 10.26 |

Sources for these numbers are in the per-material sections below (dielectric, intervalley deformation potentials, impact-ionization fits).

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

## 4. Silicon electron-phonon deformation potentials 🚧 (Part C5)
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

## 5. GaAs electron-phonon (polar: Fröhlich + intervalley) 🚧 (Part C5)
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

## 6. Collision-rate saturation ν(ε) = ν_sat[1 − exp(−(ε/ε₀)^n)] 🚧
| Quantity | Value | Source |
|---|---|---|
| ν_sat (Si) | 1.3e14 s⁻¹ | project (Fischetti priority); Meng et al., PRB 91, 075201 (2015); Fischetti-Laux 1988 |
| ν_sat (GaAs) | 1e14 s⁻¹ | Fischetti, IEEE TED 38, 634 (1991) (ab-initio Bernardi PNAS 112, 5291 ~1e14) |
| n | 2 (start) | parametrization |
| ε₀ | ~0.8 eV | parametrization |

**Never a hard min(ν, ν_sat)** — derivative discontinuity destabilizes the stiff solver.

## 7. Bandgap renormalization (II threshold) 🚧 (Part C7)
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

## 10. Energy bins / search (super-mode) 🚧 (Part C3)
| Quantity | Value/rule | Source |
|---|---|---|
| δ(ΔE) replacement | normalized Gaussian or unit-area rectangle | Stobbe 1994 (0.2 eV rect); Kunikiyo 1994 (5 meV bins) |
| σ_E | ~ mean inter-level spacing of the grid | grid-matched |
| pair search | deterministic / energy-windowed; NO Monte-Carlo | project (MC breaks CPTP) |
| sub-cycling trigger | ν_max(h/2) ≳ 0.2 → m sub-steps, (τ/2m)ν_max ≲ 0.1 | Lindblad 1976 |

## 11. Ring/pipeline MPI (super-mode) 🚧 (Part D)
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
| κ_TF² (degenerate, a.u.) | 4(3n/π)^⅓/ε = 4k_F/(πε) | Ashcroft-Mermin; arXiv:2312.13059 |
| κ_D² (Debye, a.u.) | 4πn/(εk_BT) | Ashcroft-Mermin |
| Lindhard F(x) | ½+(1−x²)/(4x)ln\|(1+x)/(1−x)\| | Lindhard, Mat.-Fys. Medd. 28, 8 (1954) |
| carrier-carrier rate scale | 1e13–1e14 s⁻¹ @ 1e17–1e19 cm⁻³ | Goodnick-Lugli PRB 37, 2578; Fischetti-Laux PRB 38, 9721 |
| EID slope (qualitative) | γ=γ₀+a·n, linear; e-carrier ≈8× e-exciton | Honold et al. PRB 40, 6442 (1989) |
| screening default | static Lindhard/RPA (b); LOPC (c) GaAs-only n≳5e17 | recommendation |

**No HF double-counting:** carrier-carrier (F) = correlation (2nd-Born/GW) self-energy, dissipative only; the static screened-exchange shift stays in Σ^HF. Carrier-carrier conserves Σf_k and ΣE_k f_k (validation invariants).
