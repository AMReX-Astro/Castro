#!/bin/bash
#SBATCH --job-name=fstrang
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --qos=regular
#SBATCH --time=04:00:00
#SBATCH --constraint=haswell

# Cori Haswell has 2 Intel "Haswell" processors each with 16 cores.
# number of nodes = # of MPI * number of threads / 32 (cores / node)

export OMP_NUM_THREADS=1


slack_job_start.py "starting flame strang problem @ NERSC" @michael

export CASTRO_EXEC=./Castro1d.intel.haswell.MPI.ex


function find_chk_file {
    # find_chk_file takes a single argument -- the wildcard pattern
    # for checkpoint files to look through
    chk=$1

    # find the latest 2 restart files.  This way if the latest didn't
    # complete we fall back to the previous one.
    #temp_files=$(find . -type d -name "${chk}" -print | sort | tail -2)
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

GLOBAL_PARAMS="
castro.dtnuc_e=0.25
amr.check_int=100"

NEEDED_FILES="
${CASTRO_EXEC}
inputs.1d.sdc
probin.sdc
helm_table.dat"

NZONES="
256
512
1024"


# Strang
RUNPARAMS="
castro.time_integration_method=0
castro.ppm_type=1"


for nz in ${NZONES}
do
    rdir=flame_${nz}_strang
    if [ ! -d ${rdir} ]; then
        mkdir ${rdir}
    fi
    if [ ! -f ${rdir}/run_done ]; then
        cd ${rdir}
        for nf in ${NEEDED_FILES}
        do
            if [ ! -f ${nf} ]; then
                cp ../${nf} .
            fi
        done

        slack_job_start.py "checking ${rdir} for restart" @michael

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

        slack_job_start.py "starting ${rdir} ${restartFile}" @michael
        srun ${CASTRO_EXEC} inputs.1d.sdc amr.plot_file=${rdir}_plt ${RUNPARAMS} amr.n_cell=${nz} ${GLOBAL_PARAMS} ${restartString} &> run.out
        touch run_done
        cd ..
    fi
done


