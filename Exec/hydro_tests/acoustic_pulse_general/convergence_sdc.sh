#!/bin/bash

mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.64 castro.sdc_order=2 castro.time_integration_method=2
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.128 castro.sdc_order=2 castro.time_integration_method=2
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.256 castro.sdc_order=2 castro.time_integration_method=2
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.2d.512 castro.sdc_order=2 castro.time_integration_method=2

RichardsonConvergenceTest2d.gnu.ex coarFile=acoustic_pulse_64_plt00081 mediFile=acoustic_pulse_128_plt00161 fineFile=acoustic_pulse_256_plt00321 mediError=med.out coarError=coar.out > convergence_sdc.1.out
RichardsonConvergenceTest2d.gnu.ex coarFile=acoustic_pulse_128_plt00161 mediFile=acoustic_pulse_256_plt00321 fineFile=acoustic_pulse_512_plt00641 mediError=med.out coarError=coar.out > convergence_sdc.2.out

