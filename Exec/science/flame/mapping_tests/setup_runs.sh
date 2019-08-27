CASTRO_EXEC=./Castro1d.gnu.MPI.ex

NEEDED_FILES="
${CASTRO_EXEC}
inputs.1d.sdc.map
probin.sdc.map
flame_4096_sdc4_plt323427.slice
helm_table.dat"

NZONES="
1024
2048
4096"

# SDC-4
RUNPARAMS=""

for nz in ${NZONES}
do
    rdir=flame_${nz}_sdc4
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

    nohup mpiexec -n 4 ${CASTRO_EXEC} inputs.1d.sdc.map amr.plot_file=${rdir}_plt ${RUNPARAMS} amr.n_cell=${nz} >& out &
    cd ..
done

