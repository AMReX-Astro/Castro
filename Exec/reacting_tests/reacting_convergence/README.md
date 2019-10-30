# reacting_convergence

This is a reaction test used for measuring convergence in reacting hydro,
built from the acoustic_pulse_general test.

# convergence testing

* in `Castro/Exec/reacting_tests/reacting_convergence`:
    ```
    make COMP=intel -j 4
    ```

* in output directory:
    ```
    mkdir reacting_convergence
    cd reacting_convergence
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/inputs.* .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/probin .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/helm_table.dat .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/job_scripts/cori.reacting_convergence.slurm .
    sbatch cori.reacting_convergence.slurm
    ```

* to do the analysis, in the output directory:
    ```
    cp ~/development/Castro//Exec/reacting_tests/reacting_convergence/analysis/check_convergence.sh .
    cp ~/development/Castro//Exec/reacting_tests/reacting_convergence/analysis/create_pretty_tables.py .
    ```

  edit it to use bash and give the proper executable name.  Then:
    ```
    ./check_convergence.sh
    ```
  This outputs 3 sets of convergence data files.  To process these
  into a set of LaTeX table rows::
    ```
    python3 create_pretty_tables.py
    ```
