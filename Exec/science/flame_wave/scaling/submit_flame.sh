#!/bin/bash
#BSUB -P ast106
#BSUB -W 50
#BSUB -nnodes 64
#BSUB -alloc_flags smt1
#BSUB -J flame_gpu
#BSUB -o flame_gpu.%J
#BSUB -e flame_gpu.%J

cd $LS_SUBCWD

inputs_file=scaling/inputs.scaling.3d

n_mpi=384 # 16 nodes * 6 gpu per node
n_omp=1
n_gpu=1
n_cores=1
n_rs_per_node=6

export OMP_NUM_THREADS=$n_omp

Castro_ex=./Castro3d.pgi.MPI.CUDA.ex

jsrun -n $n_mpi -r $n_rs_per_node -c $n_cores -a 1 -g $n_gpu $Castro_ex $inputs_file
