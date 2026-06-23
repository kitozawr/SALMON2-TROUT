# Vendored snapshot: DeePseudopot

This directory is a **vendored, read-only snapshot** of the upstream project
**DeePseudopot** by Kaixuan (Tommy) Lin et al. (Rabani group, UC Berkeley / LBNL).

- Upstream:      https://github.com/TommyLinkl/DeePseudopot
- Snapshot commit: `803e6a400e621a7a584862f70e8823afbbe60b8c`
- Snapshot date:  2026-06-11

It is included here purely as **reference material / provenance** for the
`theory='dft'` → EPM form-factor extraction workflow added to this fork
(see `tools/dft_to_epm/`): DeePseudopot is the machine-learning generalisation
of exactly the semi-empirical pseudopotential (SEPM/EPM) local form factors that
`src/epm/` consumes, so its Hamiltonian builders, k-path conventions and
`*_pot.par` form-factor files document the same physical objects our extractor
produces.

## How to cite the upstream work

> K. Lin, M. J. Coley-O'Rourke & E. Rabani,
> "Deep-learning atomistic semi-empirical pseudopotential model for nanomaterials,"
> *npj Comput. Mater.* **11**, 381 (2025).

## License note

The upstream repository did **not** ship an explicit license file at the
snapshot commit above. This vendored copy is therefore included on a
reference-only basis; redistribution/reuse terms are governed by the upstream
authors. If you intend to build on the DeePseudopot code itself (as opposed to
the SALMON-side extractor in `tools/dft_to_epm/`), please contact the upstream
authors for licensing. Nothing in this snapshot is compiled into or linked
against the SALMON binary.
