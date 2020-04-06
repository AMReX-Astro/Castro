
dir=scaling_results

mkdir -p $dir

ncell_list="256 384 512 1024 2048"
ngpu_list="6 12 24 48 72 96 144 192 384 768 1536 3072"
grid_size_list="32 48 64 96 128"

inputs=inputs
probin=probin

Castro_ex=Castro3d.pgi.MPI.CUDA.ex

if [ ! -e $dir/$Castro_ex ]; then
    cp $Castro_ex $dir/
fi

for ngpu in $ngpu_list
do

    if [ $ngpu -gt 6 ]; then
        nnodes=$(echo "$ngpu / 6" | bc)
    else
        nnodes=1
    fi

    for ncell in $ncell_list
    do

        ncell_min=0
        ncell_max=1048576

        if   [ $ngpu -eq 6 ]; then
            ncell_min=256
            ncell_max=256
        elif [ $ngpu -le 48 ]; then
            ncell_min=256
            ncell_max=512
        elif [ $ngpu -le 384 ]; then
            ncell_min=512
            ncell_max=1024
        elif [ $ngpu -le 3072 ]; then
            ncell_min=1024
            ncell_max=2048
        fi

        if [ $ncell -gt $ncell_max ] || [ $ncell -lt $ncell_min ]; then
            continue
        fi

        for grid_size in $grid_size_list
        do

            grid_size_min=0
            grid_size_max=1048576

            if [ $grid_size -gt $grid_size_max ] || [ $grid_size -lt $grid_size_min ]; then
                continue
            fi

            suffix=ncell.$ncell.ngpu.$ngpu.grid.$grid_size

            new_inputs=inputs.$suffix
            new_probin=probin.$suffix
            new_run_script=run_script_$suffix.sh

            if [ -e $dir/$new_inputs ]; then
                continue
            fi

            if [ $ngpu -lt 6 ]; then
                n_rs_per_node=$ngpu
            else
                n_rs_per_node=6
            fi

            cp $inputs $dir/$new_inputs
            cp $probin $dir/$new_probin
            sed -i "s/amr.n_cell.*/amr.n_cell = $ncell $ncell $ncell/g" $dir/$new_inputs
            sed -i "s/amr.max_grid_size.*/amr.max_grid_size = $grid_size/g" $dir/$new_inputs
            sed -i "s/amr.probin_file.*/amr.probin_file = $new_probin/g" $dir/$new_inputs

            cp run_script.sh $dir/$new_run_script
            sed -i "s/#BSUB -o.*/#BSUB -o Castro.$suffix.out/g" $dir/$new_run_script
            sed -i "s/#BSUB -e.*/#BSUB -e Castro.$suffix.out/g" $dir/$new_run_script
            sed -i "s/#BSUB -nnodes.*/#BSUB -nnodes $nnodes/g" $dir/$new_run_script
            sed -i "0,/Castro_ex.*/s//Castro_ex=$Castro_ex/" $dir/$new_run_script
            sed -i "0,/inputs.*/s//inputs=$new_inputs/" $dir/$new_run_script
            sed -i "0,/n_mpi.*/s//n_mpi=$ngpu/" $dir/$new_run_script
            sed -i "0,/n_gpu.*/s//n_gpu=1/" $dir/$new_run_script
            sed -i "0,/n_rs_per_node.*/s//n_rs_per_node=$n_rs_per_node/" $dir/$new_run_script

            cd $dir
            echo "Submitting job with suffix "$suffix
            bsub $new_run_script
            cd -

        done

    done
done
