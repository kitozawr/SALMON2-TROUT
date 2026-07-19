# Exercise x12 — Si frozen-core **phonon-assisted indirect** (Γ→X) generation

Si is an **indirect-gap** semiconductor: the valence-band maximum (Γ) and the
conduction-band minimum (the Δ/X valley) sit at different crystal momenta, so
band-edge absorption **cannot** be a vertical optical transition — it needs a
**phonon** to carry the momentum. This example drives Si with a near-tunneling
IR field in a **frozen-core 2V+2C Houston window** and shows the e-ph ring
opening exactly that indirect channel: turn the phonon off and the real-carrier
density collapses; turn it on and carriers are generated in the **lowest**
conduction band at the **X-valley** and cooled toward the band minimum.

It is the open-system counterpart to the frozen-core discussion in
[`wiki/03` §12](../../wiki/03_numerical_methods.md) and a direct validation that
the dissipators are **exactly electron-number conserving on a frozen window**.

## The active window (2 valence + 2 conduction)

`frozen_core_threshold_ev`/`frozen_free_threshold_ev` select the gap-edge
Houston window; the deep valence and high conduction bands are **frozen from
dissipation** but still evolve under the **full-basis VG unitary** (so a strong
field can push virtual population up through them and back — basis sufficiency).
For the 4³ Si primitive cell at Γ the window is

| active band | character | E−E_F |
|---|---|---|
| 3, 4 | VBM (Γ₂₅′ doublet) — **2 valence** | −1.70 eV |
| 5 | CBM (Γ₂′ singlet) | +1.70 eV |
| 6, 7 | next CB (Γ₁₅ doublet) | +3.51 eV |

so the conduction side is the CBM singlet **+** the degenerate Γ₁₅ pair (the two
lowest conduction *levels*). Bands 1–2 and 8–16 are frozen. `nstate=16` gives
the unitary its headroom; the ring/dissipators only ever see the 5 active bands.

## Run

```bash
cd samples/exercise_x12_Si_frozen_phonon_indirect

# (1) ground state (in-SALMON EPM)
../../build/salmon < Si_prim_epm_gs.inp > gs.log

# (2) phonon ON  (indirect channel open)
../../build/salmon < Si_frozen_phonon_rt.inp > on.log
cp Si_prim_sbe_nex.data on_nex.data ; cp Si_prim_sbe_nex_k_real.data on_nkr.data

# (3) phonon OFF (indirect channel closed) — flip yn_sbe_eph -> 'n'
sed "s/yn_sbe_eph               = 'y'/yn_sbe_eph               = 'n'/;\
     s/yn_sbe_eph_acoustic      = 'y'/yn_sbe_eph_acoustic      = 'n'/" \
    Si_frozen_phonon_rt.inp > off.inp
../../build/salmon < off.inp > off.log
cp Si_prim_sbe_nex.data off_nex.data ; cp Si_prim_sbe_nex_k_real.data off_nkr.data

# (4) plot: nex(t), energy distribution, cooling curve
python3 phonon_analysis.py         # -> phonon_assisted_demo.png
```

## What you should see (`phonon_assisted_demo.png`)

| | eph OFF (field only) | eph ON (phonon) |
|---|---|---|
| final `nex` | 5.9×10¹⁸ cm⁻³ | 1.1×10²¹ cm⁻³ (**≈180×**) |
| ⟨E−E_CBM⟩ (lowest CB) | 1.0 eV (hot) | 0.42 eV (cooled to the valley) |
| electrons (every step) | **8.000** | **8.000** |

The phonon opens the indirect Γ→X channel (≈180× more real carriers — most of
Si's band-edge absorption is phonon-assisted), lands them in the **lowest**
conduction band (band-5 / band-6 population ratio ≈ 33×, i.e. **not** thrown to
the top band), and cools them from 1.23 → 0.42 eV above the CBM by phonon
emission.

> Reference numbers regenerated **2026-07-19** after the CP extension of the
> frozen-window dissipators (`wiki/00`): the active↔frozen coherences now carry
> the same loss-Kraus factors as the active block, so a real scattering event
> decoheres the reversible frozen-band dressing (collision-assisted generation)
> instead of leaving ρ non-positive. The eph-OFF column is bit-identical to the
> pre-fix build; the eph-ON carrier yield is larger.

## Adding impact ionization (+ phonon-assisted II)

Impact ionization rides the **same ring** and is verified to stay exactly
electron-number conserving on this frozen window (`sbe_ii_phassist` enables the
phonon-assisted branch). In `Si_frozen_phonon_rt.inp` set both:

```
  yn_sbe_impact_ionization = 'y'
  yn_sbe_ii_holes          = 'y'     ! must be paired with impact_ionization
  sbe_ii_phassist          = 1.0d0
```

Electrons stay **8.000** (checked at I=10¹¹ and 10¹² W/cm², with Auger on, and
with the FULL channel set eph+ac+II+holes+phassist+Auger+eeh+Coulomb at
dt=10 a.u. — 8.000 at every step, no negative nelec/nhole).

## Notes / caveats

- **Basis convergence.** This field (I=10¹¹ W/cm², Keldysh γ≈1) is
  basis-converged: `nex` changes < 4 % from `nstate` 12 → 16. The built-in
  VG basis-edge monitor (`P_top` warning on stderr) stays quiet. **Deeper
  tunneling** (e.g. I=10¹² at this ω) sweeps virtual population across many
  bands — the monitor will warn and you must raise `nstate` (frozen core keeps
  that cheap, since only the active window is dissipated).
- **Grid.** The ring is O(nk³)/step; 4³ is the fast default. The indirect
  X-valley is only coarsely sampled at 4³ — raise `num_kgrid` (and MPI ranks
  dividing nk) to resolve the valley shell.
- **Trace exactness.** The frozen-window ring is trace-exact (`wiki/00`, the
  `ring_apply_dpop` fix): electron number is conserved to machine precision at
  every `dt`, so any residual `nex` is real carrier generation, not a leak.
