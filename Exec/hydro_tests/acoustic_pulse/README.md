# acoustic_pulse

This is the acoustic pulse problem from McCorquodale & Colella 2011.
We use this to measure convergence.

Note: for DIM > 1, you can still run this in a "1-d" mode by
setting the probin parameter init_as_1d to 1, 2, or 3, to just
do the pulse in that coordinate direction.


# Convergence testing

Here we detail the procedure used for convergence testing at NERSC and
processing the results for the pretty table included in the Castro SDC
paper.

  in `Castro/Exec/hydro_tests/acoustic_pulse`
    ```
    make COMP=intel -j 4
    ```

  in the output directory:
    ```
    mkdir acoustic_pulse
    cd acoustic_pulse
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/inputs.2d.* .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/probin .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/job_scripts/cori.acoustic_pulse_convergence.slurm .
    sbatch cori.acoustic_pulse_convergence.slurm
    ```

  to do the analysis, in the output directory:
    ```
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse/job_scripts/check_convergence.sh .
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse/job_scripts/create_pretty_tables.py .

  you may need to edit `check_convergence.sh` to use bash and have a
  different executable name.  Then
    ```
    ./check_convergence.sh
    ```
  This outputs `convergence.2d.lo.sdc4.out` and
 `convergence.2d.hi.sdc4.out`.

  To process these into a set of LaTeX table rows:
    ```
    python3 ./create_pretty_tables.py
    ```
