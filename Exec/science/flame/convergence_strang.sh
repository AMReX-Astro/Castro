#!/bin/bash

# echo the commands
set -x

DIM=1
EXEC=./Castro${DIM}d.gnu.MPI.TRUESDC.ex

RUNPARAMS="
castro.time_integration_method=0
castro.ppm_type=1"

mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=1024 amr.plot_file=flame_1024_strang_plt &> strang_1024.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=2048 amr.plot_file=flame_2048_strang_plt &> strang_2048.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=4096 amr.plot_file=flame_4096_strang_plt &> strang_4096.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=8192 amr.plot_file=flame_8192_strang_plt &> strang_8192.out
mpiexec -n 16 ${EXEC} inputs.1d.sdc ${RUNPARAMS} amr.n_cell=16384 amr.plot_file=flame_16384_strang_plt &> strang_16384.out



