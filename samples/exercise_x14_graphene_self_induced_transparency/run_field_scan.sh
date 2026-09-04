#!/usr/bin/env bash
# x14 headline figure, end to end: a DOPED and an INTRINSIC field scan on the same
# mesh and the same pulse, then transmission / conductivity vs field.
#
#   NK=147 EF=0.2 FIELDS=1,10,30,100,300,1000 OMP_NUM_THREADS=48 \
#   SALMON=../../build/salmon bash run_field_scan.sh
#
# Environment (all optional):
#   NK        k-mesh side, default 147. The doping must be RESOLVED: k_F = E_F/hbar v_F
#             needs ~3 mesh spacings, i.e. NK >~ 4.68/k_F(a.u.) -- 280 at E_F = 0.2 eV,
#             140 at 0.4, 93 at 0.6 (wiki/12 sec. 4a.0). NK = 147 with E_F = 0.2 eV is
#             the cheapest setting that still resolves the Fermi surface at all
#             (9 mesh points per valley inside the Fermi disc; the Drude weight is
#             then ~30 % low, a field-independent scale error).
#   EF        Fermi level from the Dirac point [eV], default 0.2 (n = 3.4e12 cm^-2,
#             the CVD-on-PET range). Pick it from the measurement: wiki/12 sec. 4a.0.
#   TINIT     temperature of the initial occupation [K], default 300.
#   FIELDS    comma-separated peak fields [kV/cm], default 1,10,30,100,300,1000.
#   VARIANTS  coh (default; the T(E0) shape, cheap) or diss / mem (adds the phonon ring,
#             which gives tau and the absolute absorption -- O(N_k^2), see the cost note
#             in README sec. 6 before choosing NK).
#   T_MEAS    measured transmissions to draw (substrate included), default "0.60 0.70".
#   N_SUB     substrate index for T_MEAS, default 1.65 (PET).
#   OUTDIR    working directory, default scan_nk<NK>_ef<EF>.
#   SALMON, OMP_NUM_THREADS as usual.
#
# Cost: with VARIANTS=coh this is the cheap half of x14 -- the unitary propagation is
# O(N_k), so NK = 300 (90000 k-points) is minutes per field on a node and the whole
# converged T(E0) curve, which is the main result, costs a few core-hours. With
# VARIANTS=diss/mem the graphene ring is O(N_k^2) and the cost jumps by (N_k/N_k')^2:
# see README sec. 6 -- 147^2 is a 7.5 h run on 48 threads, 300^2 needs MPI over k.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NK=${NK:-147}; EF=${EF:-0.2}; TINIT=${TINIT:-300}
FIELDS=${FIELDS:-1,10,30,100,300,1000}; VARIANTS=${VARIANTS:-coh}
T_MEAS=${T_MEAS:-"0.60 0.70"}; N_SUB=${N_SUB:-1.65}
OUTDIR=${OUTDIR:-scan_nk${NK}_ef${EF}}
SALMON=${SALMON:-$HERE/../../build/salmon}

mkdir -p "$OUTDIR"; cd "$OUTDIR"
echo "# x14 field scan: nk = $NK, E_F = $EF eV, T_init = $TINIT K, fields = $FIELDS"
python3 "$HERE/make_inputs.py" --nk "$NK" --nstate 2 --fields "$FIELDS" --variants "$VARIANTS" \
        --ef-ev "$EF" --temp-init-k "$TINIT" --outdir doped
python3 "$HERE/make_inputs.py" --nk "$NK" --nstate 2 --fields "$FIELDS" --variants "$VARIANTS" \
        --outdir intrinsic

for set in doped intrinsic; do
  ( cd "$set" && cp -f "$HERE/run_scan.sh" . && SALMON="$SALMON" bash run_scan.sh )
done

python3 "$HERE/field_scan_plot.py" \
        --doped "doped/runs/*/graphene_sit_sbe_rt.data" \
        --intrinsic "intrinsic/runs/*/graphene_sit_sbe_rt.data" \
        --t-meas $T_MEAS --n-sub "$N_SUB" \
        --title "Same ${NK}^2 mesh, same pulse: only the initial occupation differs" \
        --out doped_vs_intrinsic.png
python3 "$HERE/drude_check.py" "doped/runs/*/graphene_sit_sbe_rt.data" \
        --t-meas $T_MEAS --n-sub "$N_SUB"
python3 "$HERE/plot_occupation.py" doped/graphene_sit --ef-ev "$EF" --temp-init-k "$TINIT" \
        --out occupation.png
echo "# figures: $OUTDIR/doped_vs_intrinsic.png  $OUTDIR/occupation.png"
