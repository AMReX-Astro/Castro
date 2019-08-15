#!/bin/bash
#SBATCH --job-name=acoustic_pulse_general
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --qos=regular
#SBATCH --time=01:00:00
#SBATCH --constraint=haswell

echo "starting..."

export OMP_NUM_THREADS=1

export CASTRO_EXEC=./Castro2d.intel.haswell.MPI.ex

# Cori Haswell has 2 Intel "Haswell" processors each with 16 cores.
# number of nodes = # of MPI * number of threads / 32 (cores / node)

RUNPARAMS="castro.sdc_order=4 castro.time_integration_method=2 castro.limit_fourth_order=1"

srun ${CASTRO_EXEC} inputs.64 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_sdc4_plt &> 64.out
srun ${CASTRO_EXEC} inputs.128 ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_sdc4_plt &> 128.out
srun ${CASTRO_EXEC} inputs.256 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_sdc4_plt &> 256.out
srun ${CASTRO_EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_sdc4_plt &> 512.out
