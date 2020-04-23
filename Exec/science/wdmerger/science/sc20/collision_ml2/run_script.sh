#!/bin/bash
#BSUB -P AST106
#BSUB -W 120
#BSUB -nnodes 48
#BSUB -alloc_flags smt1
#BSUB -J Castro
#BSUB -o Castro.%J
#BSUB -e Castro.%J

cd $LS_SUBCWD

inputs=inputs

n_mpi=288
n_omp=1
n_gpu=1
n_cores=1
n_rs_per_node=6

Castro_ex=Castro3d.pgi.TPROF.MPI.CUDA.ex

jsrun -n $n_mpi -r $n_rs_per_node -c $n_cores -a 1 -g $n_gpu -X 1 -brs $Castro_ex $inputs
