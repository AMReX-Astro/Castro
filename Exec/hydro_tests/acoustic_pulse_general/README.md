# acoustic_pulse

This is an alternate version of the acoustic pulse problem from
McCorquodale & Colella 2011.  In this version, we perturb the pressure
and use a constant entropy via the EOS to find the density.

# Convergence testing

Here we detail the procedure used for convergence testing at NERSC and
processing the results for the pretty table included in the Castro SDC
paper.

  in `Castro/Exec/hydro_tests/acoustic_pulse_general`
    ```
    make COMP=intel -j 4
    ```

  in the output directory:
    ```
    mkdir acoustic_pulse_general
    cd acoustic_pulse_general
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/inputs.2d.* .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/probin .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/helm_table.dat .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/job_scripts/cori.acoustic_pulse_convergence.slurm .
    sbatch cori.acoustic_pulse_convergence.slurm
    ```

  to do the analysis, in the output directory:
    ```
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse_general/job_scripts/check_convergence.sh .
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse_general/job_scripts/create_pretty_tables.py .

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


