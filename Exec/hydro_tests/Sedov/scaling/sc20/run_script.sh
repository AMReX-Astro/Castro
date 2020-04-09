#!/bin/bash
#BSUB -P AST106
#BSUB -W 15
#BSUB -nnodes 1
#BSUB -alloc_flags smt1
#BSUB -J Castro
#BSUB -o Castro.%J
#BSUB -e Castro.%J

cd $LS_SUBCWD

inputs=inputs

n_mpi=1
n_omp=1
n_gpu=1
n_cores=1
n_rs_per_node=1

export OMP_NUM_THREADS=$n_omp
export OMP_PLACES="threads"
export OMP_SCHEDULE="dynamic"

Castro_ex=Castro3d.ex

jsrun -n $n_mpi -r $n_rs_per_node -c $n_cores -a 1 -g $n_gpu -X 1 -brs $Castro_ex $inputs
