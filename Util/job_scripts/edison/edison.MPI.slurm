#!/bin/bash
#SBATCH -J flame_wave
#SBATCH -A m1938
#SBATCH -N 6
#SBATCH -p regular
#SBATCH -t 04:00:00

# -N = number of nodes = # of MPI * number of threads / 24 (cores / node)

export OMP_NUM_THREADS=1

slack_job_start.py "starting flame_wave problem @ NERSC" xrb


#cd $SLURM_SUBMIT_DIR


function find_chk_file {
    # find_chk_file takes a single argument -- the wildcard pattern
    # for checkpoint files to look through
    chk=$1

    # find the latest 2 restart files.  This way if the latest didn't
    # complete we fall back to the previous one.
    temp_files=$(find . -type d -name "${chk}" -print | sort | tail -2)
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


# here, we'll rely on restartFile being global

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


# Edison has 2 Intel "IvyBridge" processors each with 12 cores.  Each processor
# connects to 1/2 of the node's memory, and is called a NUMA node.  So there
# are 2 NUMA nodes per Edison node.

# -n is the number of MPI tasks
# -c sets the number of OMP threads per MPI -- this should be the same as OMP_NUM_THREADS
srun -n 144 -c ${OMP_NUM_THREADS} ./Castro2d.intel.sandybridge.MPI.ex inputs_2d ${restartString}



