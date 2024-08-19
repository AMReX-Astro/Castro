#!/bin/bash
#SBATCH --job-name=react_conv
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --qos=regular
#SBATCH --time=04:00:00
#SBATCH --constraint=haswell

echo "starting..."

export OMP_NUM_THREADS=1


slack_job_start.py "starting reacting_convergence problem @ NERSC" @michael

export CASTRO_EXEC=./Castro2d.intel.haswell.MPI.ex

# Cori Haswell has 2 Intel "Haswell" processors each with 16 cores.
# number of nodes = # of MPI * number of threads / 32 (cores / node)


# SDC-4
RUNPARAMS="
castro.sdc_order=4
castro.time_integration_method=2
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
castro.sdc_solver=1"

srun ${CASTRO_EXEC} inputs.64 amr.plot_file=react_converge_64_sdc4_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.128 amr.plot_file=react_converge_128_sdc4_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.256 amr.plot_file=react_converge_256_sdc4_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.512 amr.plot_file=react_converge_512_sdc4_plt ${RUNPARAMS}

# SDC
RUNPARAMS="
castro.sdc_order=2
castro.time_integration_method=2
castro.ppm_type=0
castro.sdc_solve_for_rhoe=1
castro.sdc_solver_tol=1.e-10
castro.sdc_solver=1"

srun ${CASTRO_EXEC} inputs.64 amr.plot_file=react_converge_64_sdc_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.128 amr.plot_file=react_converge_128_sdc_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.256 amr.plot_file=react_converge_256_sdc_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.512 amr.plot_file=react_converge_512_sdc_plt ${RUNPARAMS}

# Strang
RUNPARAMS="castro.time_integration_method=0"

srun ${CASTRO_EXEC} inputs.64 amr.plot_file=react_converge_64_strang_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.128 amr.plot_file=react_converge_128_strang_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.256 amr.plot_file=react_converge_256_strang_plt ${RUNPARAMS}
srun ${CASTRO_EXEC} inputs.512 amr.plot_file=react_converge_512_strang_plt ${RUNPARAMS}
