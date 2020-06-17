#!/bin/bash

# echo the commands
set -x

DIM=3
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS="castro.ppm_type=0 castro.time_integration_method=0"

mpiexec -n 8 ${EXEC} inputs.3d.64 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_plm_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs.3d.128 ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_plm_plt &> 128.out
mpiexec -n 8 ${EXEC} inputs.3d.256 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_plm_plt &> 256.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_plm_plt00081 mediFile=acoustic_pulse_128_plm_plt00161 fineFile=acoustic_pulse_256_plm_plt00321 > convergence.${DIM}d.lo.plm.out

#mpiexec -n 4 ${EXEC} inputs.2d.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_plm_plt &> 512.out

#RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_plm_plt00161 mediFile=acoustic_pulse_256_plm_plt00321 fineFile=acoustic_pulse_512_plm_plt00641 > convergence.${DIM}d.hi.plm.out

