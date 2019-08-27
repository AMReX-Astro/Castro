#!/bin/bash
#SBATCH --job-name=f3anw-strang
#SBATCH --account=m3018
#SBATCH --nodes=20
#SBATCH --tasks-per-node=32
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --constraint=haswell

# Cori Haswell has 2 Intel "Haswell" processors each with 16 cores.
# number of nodes = # of MPI * number of threads / 32 (cores / node)

export OMP_NUM_THREADS=1

export CASTRO_EXEC=./Castro2d.intel.haswell.MPI.ex


function find_chk_file {
    # find_chk_file takes a single argument -- the wildcard pattern
    # for checkpoint files to look through
    chk=$1

    # find the latest 2 restart files.  This way if the latest didn't
    # complete we fall back to the previous one.
    temp_files=$(find . -maxdepth 1 -name "${chk}" -print | sort | tail -2)
    restartFile=""
    for f in ${temp_files}
    do
        # the Header is the last thing written -- check if it's there, otherwise,
        # fall back to the second-to-last check file written
        if [ ! -f ${f}/Header ]; then
            restartFile=""
        else
            restartFile="${f}"
        fi
    done

}

# look for 7-digit chk files
find_chk_file "*chk???????"

if [ "${restartFile}" = "" ]; then
    # look for 6-digit chk files
    find_chk_file "*chk??????"
fi

if [ "${restartFile}" = "" ]; then
    # look for 5-digit chk files
    find_chk_file "*chk?????"
fi

# restartString will be empty if no chk files are found -- i.e. new run
if [ "${restartFile}" = "" ]; then
    restartString=""
else
    restartString="amr.restart=${restartFile}"
fi

slack_job_start.py "starting flame_wave problem @ NERSC; ${restartFile}" @michael

srun ${CASTRO_EXEC} inputs.boost_10_10 ${restartString}

