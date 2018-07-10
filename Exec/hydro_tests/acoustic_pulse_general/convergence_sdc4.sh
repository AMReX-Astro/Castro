#!/bin/bash

mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.64 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.128 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.256 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1

RichardsonConvergenceTest2d.gnu.ex coarFile=acoustic_pulse_64_plt00101 mediFile=acoustic_pulse_128_plt00201 fineFile=acoustic_pulse_256_plt00401 mediError=med.out coarError=coar.out > convergence_sdc4.1.out

mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.512 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1

RichardsonConvergenceTest2d.gnu.ex coarFile=acoustic_pulse_128_plt00201 mediFile=acoustic_pulse_256_plt00401 fineFile=acoustic_pulse_512_plt00801 mediError=med.out coarError=coar.out > convergence_sdc4.2.out

