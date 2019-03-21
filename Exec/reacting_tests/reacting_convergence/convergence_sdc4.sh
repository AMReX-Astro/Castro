#!/bin/bash

mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs.64 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.128 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.256 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1

RichardsonConvergenceTest2d.gnu.ex coarFile=react_converge_64_plt00301 mediFile=react_converge_128_plt00601 fineFile=react_converge_256_plt01201 mediError=med.out coarError=coar.out > convergence.1.sdc4.out

mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.512 castro.sdc_order=4 castro.time_integration_method=2 castro.fourth_order=1 castro.limit_fourth_order=1 castro.use_reconstructed_gamma1=1

RichardsonConvergenceTest2d.gnu.ex coarFile=react_converge_128_plt00601 mediFile=react_converge_256_plt01201 fineFile=react_converge_512_plt02401 mediError=med.out coarError=coar.out > convergence.2.sdc2.out



