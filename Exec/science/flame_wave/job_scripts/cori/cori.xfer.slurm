#!/bin/ksh
#SBATCH -M esedison
#SBATCH -J xrb-hpss-xfer
#SBATCH -A m1938
#SBATCH -p xfer
#SBATCH -t 12:00:00

cd $SLURM_SUBMIT_DIR

# do our archiving
pidfile=process.pid

./process.xrb

PID=$!
trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL

rm -f process.pid
