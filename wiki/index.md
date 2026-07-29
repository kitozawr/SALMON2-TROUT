# Wiki

Full documentation for the SALMON2-TROUT fork. Start at the
[README](../) for the pitch, the master equation and the quick start;
these pages are the detail.

| page | what's in it |
|---|---|
| [00 — Implementation status](00_implementation_status.md) | **Read first on resume.** What is done, what is next, and the decisions log (gotchas that must not be re-litigated). |
| [01 — Physics models](01_physics_models.md) | The physical content of each term: bands, coupling, the channel inventory. |
| [02 — Constants](02_constants.md) | Every material constant with its literature source. Strict provenance — no fabricated numbers. |
| [03 — Numerical methods](03_numerical_methods.md) | Propagators, the `dt` criterion, the frozen-core window, MPI/k-point layout. |
| [04 — Configuration examples](04_configuration_examples.md) | Worked `&sbe` inputs + the canonical parameter reference (all flags, defaults, guards). |
| [05 — Folding & unfolding](05_folding_unfolding.md) | Cubic supercell ↔ primitive cell, the unfold pipeline and its spectral weights. |
| [06 — VG basis & N_b convergence](06_vg_basis_nb_convergence.md) | Velocity-gauge basis sufficiency, the `dt` artefact, and the **band budget vs field strength**. |
| [07 — Nonlocal Auger](07_nonlocal_auger.md) | The momentum-conserving II/Auger pair: kernel, screening, umklapp, and the C(n) validation. |
| [08 — Master equation](08_master_equation.md) | The complete mathematical specification — every effect as an explicit term. |
| [09 — Plotting & analysis](09_plotting_and_analysis.md) | `plot_sbe_results.py` (maps, BZ views, spectral movies, `--levels`) and the injection probes. |
| [10 — Open quantum systems literature](10_open_quantum_systems_literature.md) | Conspects of the cited open-system papers + the non-Markovian memory work. |
| [11 — Supercomputer (Lomonosov-2)](11_supercomputer_lomonosov2.md) | Building and running at scale: cost model, measured levers, OpenMP × MPI split. |

---

**Conventions used throughout.** Defaults are OFF — every dissipation channel is
opt-in. Every constant is cited. Guards abort rather than silently double-count
(e.g. BGR vs Σ<sup>HF</sup>). Where a number is an upper estimate rather than a
converged result, the page says so.
