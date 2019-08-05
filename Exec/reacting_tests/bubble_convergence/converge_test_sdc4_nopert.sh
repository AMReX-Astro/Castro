#!/bin/bash

RUNPARAMS="
castro.sdc_order=4
castro.time_integration_method=2
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1"

mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs_2d.64 ${RUNPARAMS} castro.do_react=0 amr.probin_file=probin.64.nopert
mpiexec -n 8 ./Castro2d.gnu.MPI.ex inputs_2d.128 ${RUNPARAMS} castro.do_react=0 amr.probin_file=probin.128.nopert
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs_2d.256 ${RUNPARAMS} castro.do_react=0 amr.probin_file=probin.256.nopert
