#!/bin/bash
#SBATCH -N 32
#SBATCH -C knl,quad,cache
#SBATCH -p regular
#SBATCH -t 00:45:00
#SBATCH -J test3
#SBATCH -A m1938


export OMP_PROC_BIND=close 
export OMP_PLACES=threads

export OMP_NUM_THREADS=16 

cd $SLURM_SUBMIT_DIR



srun -n 512 -c 16 --cpu_bind=cores ./Castro3d.gnu.mic-knl.MPI.OMP.ex inputs_test_wdmerger_3D



