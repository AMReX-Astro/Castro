#!/bin/bash
#SBATCH -A AST106
#SBATCH -J testing
#SBATCH -o %x-%j.out
#SBATCH -t 00:30:00
#SBATCH -p batch
# here N is the number of compute nodes
#SBATCH -N 64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest

EXEC=Castro3d.hip.x86-trento.MPI.HIP.ex
INPUTS=inputs.He.25cm.static.1000Hz

module load cpe
module load PrgEnv-gnu
module load cray-mpich
module load craype-accel-amd-gfx90a
module load rocm/6.3.1

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
export NMPI_PER_NODE=8
export TOTAL_NMPI=$(( ${SLURM_JOB_NUM_NODES} * ${NMPI_PER_NODE} ))

srun -n${TOTAL_NMPI} -N${SLURM_JOB_NUM_NODES} --ntasks-per-node=8 --gpus-per-task=1 ./$EXEC $INPUTS max_step=25
