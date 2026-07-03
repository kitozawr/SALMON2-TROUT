# Exercise x10 — Si primitive-cell SBE with **ALL dissipation channels**

The open-system template: same Si FCC 2-atom primitive cell, field files and
HHG tooling as **exercise_x8** (which is the dissipation-free basis/HHG study —
read it first), but the SBE step enables **every dissipation channel**, all
riding the nonlocal ring. Intended as the starting point for long open-system
experiments (e.g. the free-carrier **bleaching** run: pump → hot carriers →
e-ph cooling → Drude/intra-band current response).

| channel | flag | what it is |
|---|---|---|
| ring / super-mode | `yn_sbe_superres='y'` | makes e-ph and II momentum-resolved **inter-k** (Si e-ph is all-intervalley → needs the ring) |
| electron-phonon | `yn_sbe_eph='y'` | Si 6 intervalley g/f modes (cited registry values) |
| impact ionization | `yn_sbe_impact_ionization='y'` | nonlocal, crystal-momentum conserving on the ring |
| Auger | `yn_sbe_auger='y'` | exact time-reverse of the II (detailed balance — no separate C; **requires** the II channel) |
| carrier-carrier | `yn_sbe_eeh='y'` | e-e/e-h Fermi-Dirac thermalization (cited Si rate) |
| decoherence | `sbe_decoh_temperature_k=300`, `sbe_decoh_tau_m_fs=30` | Kuhn-Zurek gauge-covariant dephasing (phenomenological; either ≤ 0 disables) |
| BGR threshold | `yn_sbe_bgr_threshold='y'` | density-dependent II threshold. **Mutually exclusive with `yn_sbe_coulomb`** (gap double-count guard): Σ^HF off → BGR is the stand-in |

`yn_sbe_coulomb` (Σ^HF mean-field renormalization) is **off** by default here —
it is not a dissipation channel; to use it set it `'y'` **and**
`yn_sbe_bgr_threshold='n'`, or the run aborts by design.

## Cost / grid

The ring II + Auger pair is **O(nk³)/step** — fine at 4³–8³, prohibitive
around 20³. Default here is **8³** (x8 uses 16³, but it is dissipation-free).
Drop both inputs to `4, 4, 4` for a fast smoke test; raise `nt` (and MPI ranks
dividing nk) for production runs.

## Run

```bash
cd samples/exercise_x10_Si_full_dissipation

# (1) ground state: in-SALMON Fortran EPM (also emits Si_prim_bandpath.data)
../../build/salmon < Si_prim_epm_gs.inp > gs.log

# (2) real-time SBE, all channels on (edit file_input1 to pick the amplitude)
../../build/salmon < Si_prim_sbe_rt.inp | tee rt.log

# (3) plots + HHG (identical tooling to x8)
python3 ../../plot_sbe_results.py -i . -o plots --spectral --valleys
python3 hhg_spectrum.py
```

## What to check

* **CPTP**: the `electrons` column stays at `8.000` for every step — every
  channel is a trace-preserving GKLS map; any drift is a bug, stop and report.
* The startup banner lists each enabled channel with its cited constants
  (6 phonon modes for Si, II magnitude, eeh rate, decoherence λ).
* `Si_sbe_nex.data`: with dissipation on, nex should be *lower* at late times
  than an x8-style clean run at the same field (Auger recombination), and the
  energy (`Si_sbe_rt_energy.data`) relaxes after the pulse (e-ph cooling).
* For bleaching-type analysis use the intra-band (Drude) current
  (`Si_sbe_intra_current.data`, `yn_out_intraband_current='y'`) and the
  conductivity plots from `plot_sbe_results.py` (needs ps-scale runs for THz
  resolution — see the plotter's note).

## Field files

Same amplitude ladder as x8 (A/eV/fs; col 1 = time [fs], never rescale it):
`100.txt` A_max=2.76 (VG overflow demo) / `10.txt` 0.28 (default) / `1.txt`
0.028 / `0.1.txt` 0.0028. Basis-sufficiency methodology (P_top monitor,
nex ∈ [0, nelec]) — see exercise_x8's README; with `nstate=20` the 10.txt
drive is clean.
