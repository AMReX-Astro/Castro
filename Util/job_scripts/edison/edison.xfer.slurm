#!/bin/ksh
#SBATCH --qos=xfer
#SBATCH -J xrb-hpss-xfer
#SBATCH -t 12:00:00
#SBATCH --licenses=SCRATCH

# note: you need to do
#
#  module laod esslurm
#  sbatch edison.xfer.slurm
#
# to launch a transfer job

cd $SLURM_SUBMIT_DIR

# do our archiving
pidfile=process.pid

./process.xrb

PID=$!
trap 'kill -s TERM $PID' EXIT TERM HUP XCPU KILL

rm -f process.pid
