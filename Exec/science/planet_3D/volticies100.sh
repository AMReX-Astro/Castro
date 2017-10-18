#!/bin/bash
#PBS -l nodes=1:ppn=12,walltime=00:02:00
#PBS -N voltice_100
#PBS -q debug

module load  mpiexec/0.84_432
module load gcc/6.1.0
module load anaconda/2


cd /gpfs/home/taryu/project7_HJ/Castro/Exec/science/planet/
make realclean
make clean
make -j12
#mpiexec -n 12 ./Castro2d.gnu.MPI.ex inputs_2d

echo End Job
