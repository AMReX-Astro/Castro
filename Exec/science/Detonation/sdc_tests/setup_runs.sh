CASTRO_EXEC=./Castro1d.gnu.ex

NEEDED_FILES="
${CASTRO_EXEC}
inputs.1d.sdc
probin.sdc
helm_table.dat"

NZONES="
256
512
1024"

CFL="
0.02
0.1
0.5"

GLOBAL_RUNPARAMS="
stop_time=3.e-6"

#=============================================================================
# Lobatto SDC-4
#=============================================================================

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=4
castro.sdc_quadrature=1
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
castro.sdc_solve_for_rhoe=1
castro.sdc_solver_tol_dens=1.e-8
castro.sdc_solver_tol_spec=1.e-8
castro.sdc_solver_tol_ener=1.e-8
castro.sdc_solver=2"

for c in ${CFL}
do

    for nz in ${NZONES}
    do
        rdir=det_z${nz}_c${c}_lobatto_sdc4
        if [ ! -d ${rdir} ]; then
            mkdir ${rdir}
        fi

        cd ${rdir}
        for nf in ${NEEDED_FILES}
        do
            if [ ! -f ${nf} ]; then
                cp ../${nf} .
            fi
        done

        nohup ${CASTRO_EXEC} inputs.1d.sdc amr.plot_file=${rdir}_plt ${GLOBAL_RUNPARAMS} ${RUNPARAMS} castro.cfl=${c} amr.n_cell=${nz} >& out &
        cd ..
    done
done


#=============================================================================
# Radau SDC-4
#=============================================================================

RUNPARAMS="
castro.time_integration_method=2
castro.sdc_order=4
castro.sdc_quadrature=2
castro.limit_fourth_order=1
castro.use_reconstructed_gamma1=1
castro.sdc_solve_for_rhoe=1
castro.sdc_solver_tol_dens=1.e-8
castro.sdc_solver_tol_spec=1.e-8
castro.sdc_solver_tol_ener=1.e-8
castro.sdc_solver=2"

for c in ${CFL}
do

    for nz in ${NZONES}
    do
        rdir=det_z${nz}_c${c}_radau_sdc4
        if [ ! -d ${rdir} ]; then
            mkdir ${rdir}
        fi

        cd ${rdir}
        for nf in ${NEEDED_FILES}
        do
            if [ ! -f ${nf} ]; then
                cp ../${nf} .
            fi
        done

        nohup ${CASTRO_EXEC} inputs.1d.sdc amr.plot_file=${rdir}_plt ${GLOBAL_RUNPARAMS} ${RUNPARAMS} castro.cfl=${c} amr.n_cell=${nz} >& out &
        cd ..
    done
done



#=============================================================================
# Strang CTU
#=============================================================================

RUNPARAMS="
castro.time_integration_method=0
castro.ppm_type=1
castro.ppm_reference_eigenvectors=1"

for c in ${CFL}
do

    for nz in ${NZONES}
    do
        rdir=det_z${nz}_c${c}_strang_ctu
        if [ ! -d ${rdir} ]; then
            mkdir ${rdir}
        fi

        cd ${rdir}
        for nf in ${NEEDED_FILES}
        do
            if [ ! -f ${nf} ]; then
                cp ../${nf} .
            fi
        done

        nohup ${CASTRO_EXEC} inputs.1d.sdc amr.plot_file=${rdir}_plt ${GLOBAL_RUNPARAMS} ${RUNPARAMS} castro.cfl=${c} amr.n_cell=${nz} >& out &
        cd ..
    done
done


