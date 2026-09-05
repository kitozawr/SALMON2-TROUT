#!/usr/bin/env bash
# x14: ground state once, then every rt_*.inp (or the given pattern) in its own runs/<name>/
# directory; the GS files and the field files are symlinked in.
#   ./run_scan.sh                      # everything
#   ./run_scan.sh 'rt_E100kVcm_*.inp'  # one field
set -euo pipefail
SALMON=${SALMON:-../../../build/salmon}
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
    for f in DAST_*.txt; do [ -f "$f" ] && ln -sf ../../$f "$d/$f"; done
    if ( cd "$d" && $SALMON < ../../$inp > run.log 2>&1 ); then
        msg="done  $name  ($(grep -c . "$d/${SYS}_sbe_rt.data") rows)"
    else
        msg="FAILED $name -- see $d/run.log"
    fi
    echo "$msg"
    # Also straight to a file, unbuffered by construction. A driver that pipes this
    # script through grep/tee block-buffers stdout, so a watcher tailing the driver's
    # log sees nothing until the whole scan ends; runs/scan_progress.log is immediate.
    echo "$(date +%FT%T)  $msg" >> runs/scan_progress.log
done
echo "# analyse:  python3 ../transmission.py runs/*/${SYS}_sbe_rt.data --plot"
echo "#           python3 ../saturation_check.py runs/E100kVcm_mem/${SYS}"
