#!/bin/bash
#==============================================================================
# x12 Keldysh-bracket benchmark -- 9x9x9 Si, SLURM / sbatch (Lomonosov-2).
#
# Submits ONE job that runs the whole workflow in order:
#   (1) EPM ground state   (theory='epm', cheap)        -> Si_prim_{k,eigen,tm}.data
#   (2) SBE real time, eph ON   (indirect channel open)  -> on_*.data
#   (3) SBE real time, eph OFF  (indirect channel closed) -> off_*.data
# Then run the Python analysis on the login node (see the README).
#
# SALMON reads its namelist from STDIN. Only MPI rank 0 reads stdin (it copies
# the input to .namelist.tmp, which every rank then reads), so `srun ./salmon <
# input.inp` is the correct launch -- there is NO input-file command-line flag.
# Run from a SHARED-filesystem directory (scratch): rank 0 writes .namelist.tmp
# and the GS *.data there and every node reads them back.
#
# k-points at 9^3 = 729. They are split across ranks evenly (+/-1); nproc need
# NOT divide 729. Use nproc <= 729; for perfect load balance pick a divisor of
# 729 = 3^6 (27, 81, or 243). This template uses 81 (9 k-points/rank).
#==============================================================================

#SBATCH --job-name=x12_si_9
#SBATCH --output=x12_9x9x9_%j.out
#SBATCH --error=x12_9x9x9_%j.err
### EDIT for your machine (these are Lomonosov-2 "compute" defaults) ###
#SBATCH --partition=compute
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=14          # Lomonosov-2 nodes = 14 cores
#SBATCH --time=02:00:00
### #SBATCH --account=<your_project>   # uncomment + set if your site requires it

set -euo pipefail

#------------------------------------------------------------------------------
# EDIT: the SALMON binary and the module environment you build with.
#------------------------------------------------------------------------------
SALMON=./salmon                       # path to the built binary (gfortran+OpenMPI)
# module purge
# module load openmpi gcc              # <- your Lomonosov-2 modules (OpenMPI + gfortran)

# Total MPI ranks = nodes * ntasks-per-node. With the SBATCH above: 6*14 = 84.
# We ask srun for exactly 81 (a clean divisor of 729); the 3 spare slots idle.
NRANKS=81

# Run from the directory that holds the .inp files (a shared scratch path).
cd "${SLURM_SUBMIT_DIR}"

echo "=== (1) ground state (9x9x9 EPM) ==="
srun --ntasks=1 "${SALMON}" < Si_prim_epm_gs_9x9x9.inp > gs_9.log
#   ^ GS is per-k independent and tiny; one rank is plenty and avoids any
#     multi-rank cost for a few-second diagonalization. (You MAY use more.)

echo "=== (2) SBE real time -- eph ON ==="
srun --ntasks=${NRANKS} "${SALMON}" < Si_frozen_phonon_rt_9x9x9.inp > on_9.log
cp Si_prim_sbe_nex.data        on_nex.data
cp Si_prim_sbe_nex_k_real.data on_nkr.data

echo "=== (3) SBE real time -- eph OFF (reference) ==="
sed "s/yn_sbe_eph               = 'y'/yn_sbe_eph               = 'n'/;\
     s/yn_sbe_eph_acoustic      = 'y'/yn_sbe_eph_acoustic      = 'n'/" \
    Si_frozen_phonon_rt_9x9x9.inp > off_9.inp
srun --ntasks=${NRANKS} "${SALMON}" < off_9.inp > off_9.log
cp Si_prim_sbe_nex.data        off_nex.data
cp Si_prim_sbe_nex_k_real.data off_nkr.data

echo "=== done. Now on the login node: ==="
echo "    python3 phonon_analysis.py    # -> phonon_assisted_demo.png"
echo "    python3 rate_benchmark.py     # -> rate_benchmark.png (Keldysh bracket)"
echo "    python3 calibration_scan.py   # if you swept I_wcm2_1 (see README)"
