#!/usr/bin/env bash
#
# End-to-end DFT -> EPM form-factor extraction for Silicon (conventional cubic
# cell), using the existing low-precision DFT sample. Produces a form-factor
# table that theory='epm' can read via epm_material='file'.
#
# Usage:  SALMON=/path/to/build/salmon  bash run_dft_to_epm.sh
#
set -euo pipefail

SALMON="${SALMON:-salmon}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
GS_SAMPLE="$REPO/samples/exercise_04_bulkSi_gs"
TOOL="$REPO/tools/dft_to_epm/dft_to_epm.py"

WORK="$HERE/work"
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
cp "$GS_SAMPLE/Si_gs.inp" "$GS_SAMPLE/Si_band.inp" "$GS_SAMPLE/Si_rps.dat" .

echo "==> [1/3] DFT ground state (theory='dft')"
"$SALMON" < Si_gs.inp > dft_gs.out 2>&1
ln -sfn data_for_restart restart

echo "==> [2/3] Band structure (theory='dft_band') -> band.dat"
"$SALMON" < Si_band.inp > dft_band.out 2>&1

echo "==> [3/3] Fit EPM local form factors from band.dat"
# a = 5.43 Angstrom = 10.2626 Bohr; diamond => V^A = 0 (no --shells-a);
# nval = 32 electrons / 2 = 16 valence bands.
python3 "$TOOL" \
    --dft band.dat --format band_dat --cell cubic \
    --a-lattice-au 10.2626 --cutoff-ry 11.1 \
    --material-name Si_fromDFT --shells-s 3,8,11 \
    --nval 16 --nbands-fit 18 --weight-valence 3.0 \
    --out-prefix Si_fromDFT

echo
echo "==> done. Form factors written to:"
echo "    $WORK/Si_fromDFT_epm_formfactors.data"
echo
echo "    To build the EPM ground state from them, copy that file next to"
echo "    Si_epm_fromdft.inp (in the sample dir) and run:"
echo "        \$SALMON < Si_epm_fromdft.inp"
