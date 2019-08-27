#!/bin/bash
#SBATCH --job-name=fsdc2
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --qos=regular
#SBATCH --time=06:00:00
#SBATCH --constraint=haswell

# Cori Haswell has 2 Intel "Haswell" processors each with 16 cores.
# number of nodes = # of MPI * number of threads / 32 (cores / node)

export OMP_NUM_THREADS=1


slack_job_start.py "starting flame problem @ NERSC" @michael

export CASTRO_EXEC=./Castro1d.intel.haswell.MPI.ex


function find_chk_file {
    # find_chk_file takes a single argument -- the wildcard pattern
    # for checkpoint files to look through
    chk=$1

    # find the latest 2 restart files.  This way if the latest didn't
    # complete we fall back to the previous one.
    temp_files=$(find . -type d -name "${chk}" -print | sort | tail -2)
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
1024
2048
4096
8192"


# SDC-2
RUNPARAMS="
castro.sdc_order=2
castro.time_integration_method=2
castro.sdc_solve_for_rhoe=1
castro.sdc_solver_tol_dens=1.e-8
castro.sdc_solver_tol_spec=1.e-8
castro.sdc_solver_tol_ener=1.e-6
castro.sdc_solver_atol=1.e-10
castro.sdc_solver=1
castro.sdc_use_analytic_jac=1
castro.sdc_solver_relax_factor=1"


for nz in ${NZONES}
do
    rdir=flame_${nz}_sdc2
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

        find_chk_file "*chk?????"
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

