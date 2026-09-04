#!/usr/bin/env bash
# x14: ground state once, then every rt_*.inp in its own runs/<name>/ directory
# (the GS files are symlinked in). Usage:  ./run_scan.sh [pattern]   e.g. ./run_scan.sh 'rt_E10kVcm_*.inp'
set -euo pipefail
SALMON=${SALMON:-../../build/salmon}
SYS=graphene_sit
PATTERN=${1:-rt_*.inp}
if [ ! -f ${SYS}_k.data ]; then
    echo "# ground state ..."
    $SALMON < ${SYS}_epm_gs.inp > ${SYS}_gs.log
fi
mkdir -p runs
for inp in $PATTERN; do
    name=${inp%.inp}; name=${name#rt_}
    d=runs/$name
    mkdir -p "$d"
    for f in ${SYS}_k.data ${SYS}_eigen.data ${SYS}_tm.data; do ln -sf ../../$f "$d/$f"; done
    ( cd "$d" && $SALMON < ../../$inp > run.log 2>&1 ) && echo "done  $name  ($(grep -c . "$d/${SYS}_sbe_rt.data") rows)" \
        || { echo "FAILED $name -- see $d/run.log"; }
done
echo "# analyse:  python3 transmission.py runs/*/${SYS}_sbe_rt.data --plot"
echo "#           python3 saturation_check.py runs/E100kVcm_mem/${SYS}"
