#!/bin/bash
#BSUB -P ast106
#BSUB -W 25
#BSUB -nnodes 4
#BSUB -alloc_flags smt1
#BSUB -J Sedov_cpu
#BSUB -o Sedov_cpu.%J
#BSUB -e Sedov_cpu.%J

cd $LS_SUBCWD

inputs_file=inputs.3d.sph_1level  

n_mpi=24 # 4 nodes * 6 mpi per node
n_omp=7
n_gpu=0
n_cores=7
n_rs_per_node=6

export OMP_NUM_THREADS=$n_omp

Castro_ex=./Castro3d.pgi.MPI.OMP.ex

jsrun -n $n_mpi -r $n_rs_per_node -c $n_cores -a 1 -g $n_gpu $Castro_ex $inputs_file
