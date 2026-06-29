"""
Silicon PRIMITIVE-CELL EPM ground state -- 2-atom FCC (diamond), NON-orthogonal,
NO folding. The Si companion to epm_si_reference.py (cubic 8-atom supercell),
built on exactly the same primitive-FCC machinery as epm_gaas_primitive.py.

Silicon is the diamond structure: two IDENTICAL atoms per primitive cell, so the
antisymmetric structure factor vanishes (V^A == 0 for every shell) -- only V^S
enters. Geometrically Si is the same FCC Bravais lattice as zincblende GaAs, so
the primitive-cell recipe is identical: the plane-wave basis is the cubic
(2pi/a)(h,k,l) set restricted to one parity class (h,k,l all even or all odd =
the BCC reciprocal lattice). No cosets, no folding, no unfold map.

We REUSE epm_gaas_primitive verbatim and just re-point its module globals at Si
(MATERIAL='Si' -> Kunikiyo form factors, a = 10.26 Bohr). This mirrors how
epm_si_reference reuses epm_gaas_reference, so there is one validated code path.

Unlike GaAs (direct gap at Gamma), Si is INDIRECT: the conduction minimum sits
near X (the 6 Delta valleys at ~0.85*(2pi/a) along <100>). The primitive-cell
SBE therefore populates the X/Delta region, not Gamma -- the unfolded test of
whether the folded supercell distorts that ordering too.

Usage:
    python3 epm_si_primitive.py            # report gaps + write GS + bandpath
    python3 epm_si_primitive.py gap        # just the Gamma/X/L vertical gaps
    python3 epm_si_primitive.py gs         # GS files + bandpath
    python3 epm_si_primitive.py bandpath   # bandpath only
"""
import sys

import epm_gaas_primitive as prim


def configure_for_si(variant='Si'):
    """Re-point the primitive-FCC machinery at Silicon (diamond, V^A=0)."""
    assert variant in ('Si', 'Si_cb'), "variant must be 'Si' (Kunikiyo) or 'Si_cb'"
    prim.MATERIAL = variant
    prim.SYSNAME = 'Si_prim' if variant == 'Si' else 'Si_cb_prim'
    prim.A_LATTICE_AU = 10.26        # Si lattice constant [Bohr] = 5.431 Ang
    # Si's conduction minimum is the Delta valley (~0.85*X), whose camel-back
    # needs a larger plane-wave basis than GaAs's Gamma minimum: at 27 Ry the
    # indirect gap is 1.059 eV @ 0.850*(2pi/a) (converged Kunikiyo value); at the
    # GaAs cutoff (11.1) the dip is unresolved and the CBM mislands at X. The
    # cutoff only sizes the one-off GS diagonalization -- the SBE uses 8 bands.
    prim.PW_CUTOFF_RY = 27.0
    prim.NUM_KGRID = (8, 8, 8)
    prim.NELEC = 8                   # 2 atoms x 4 valence e- -> 8 -> 4 filled bands
    prim.NSTATE = 8
    return prim


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    configure_for_si('Si')
    if mode == 'gap':
        prim.report_gap()
    elif mode == 'bandpath':
        prim.main_bandpath()
    elif mode == 'gs':
        prim.main_gs()
        prim.main_bandpath()
    else:
        prim.report_gap()
        print()
        prim.main_gs()
        prim.main_bandpath()
