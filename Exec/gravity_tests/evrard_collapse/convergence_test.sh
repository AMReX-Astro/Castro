#!/bin/sh

EXEC=./Castro3d.gnu.MPI.ex
mpiexec -n 8 ${EXEC} inputs amr.plot_file=evrard_64_plt &> 64.out
mpiexec -n 8 ${EXEC} inputs amr.plot_file=evrard_128_plt amr.n_cell=128 128 128 &> 128.out
mpiexec -n 8 ${EXEC} inputs amr.plot_file=evrard_256_plt amr.n_cell=256 256 256 &> 256.out

