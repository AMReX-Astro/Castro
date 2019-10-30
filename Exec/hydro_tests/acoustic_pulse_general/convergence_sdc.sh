#!/bin/bash

# echo the commands
set -x

DIM=2
EXEC=./Castro${DIM}d.gnu.MPI.ex

RUNPARAMS=" castro.sdc_order=2 castro.time_integration_method=2"

mpiexec -n 8 ${EXEC} inputs.64 ${RUNPARAMS} amr.plot_file=acoustic_pulse_64_sdc_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs.128 ${RUNPARAMS} amr.plot_file=acoustic_pulse_128_sdc_plt &> 128.out
mpiexec -n 8 ${EXEC} inputs.256 ${RUNPARAMS} amr.plot_file=acoustic_pulse_256_sdc_plt &> 256.out
mpiexec -n 8 ${EXEC} inputs.512 ${RUNPARAMS} amr.plot_file=acoustic_pulse_512_sdc_plt &> 512.out

RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_64_sdc_plt00101 mediFile=acoustic_pulse_128_sdc_plt00201 fineFile=acoustic_pulse_256_sdc_plt00401 > convergence.${DIM}d.lo.sdc.out
RichardsonConvergenceTest${DIM}d.gnu.ex coarFile=acoustic_pulse_128_sdc_plt00201 mediFile=acoustic_pulse_256_sdc_plt00401 fineFile=acoustic_pulse_512_sdc_plt00801 > convergence.${DIM}d.hi.sdc.out

