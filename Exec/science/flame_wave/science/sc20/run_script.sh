#!/bin/bash
#BSUB -P AST106
#BSUB -W 360
#BSUB -nnodes 768
#BSUB -alloc_flags smt1
#BSUB -J Castro
#BSUB -o Castro.%J
#BSUB -e Castro.%J

cd $LS_SUBCWD

jsrun -n 4608 -r 6 -c 1 -a 1 -g 1 -X 1 -brs ./Castro3d.pgi.TPROF.MPI.CUDA.ex inputs
