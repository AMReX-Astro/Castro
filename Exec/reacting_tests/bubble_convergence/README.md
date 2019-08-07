# bubble_convergence

This problem was developed to test the convergence of the 4th order
SDC solver with gravity and reactions.

# convergence testing

* in `Castro/Exec/reacting_tests/bubble_convergence`
    ```
    make COMP=intel -j 4
    ```

* in output directory:
    ```
    mkdir bubble_convergence
    cd bubble_convergence
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/helm_table.dat .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/inputs_2d.* .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/probin.* .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/job_scripts/cori.bubble_convergence.slurm .
    sbatch cori.bubble_convergence.slurm
    ```

* to measure convergence:

  in the `bubble_convergence` output directory:
    ```
    cp ~/development/Castro//Exec/reacting_tests/bubble_convergence/job_scripts/check_convergence.sh .
    cp ~/development/Castro//Exec/reacting_tests/bubble_convergence/job_scripts/create_pretty_tables.py .
    ```

  edit it to use bash and give the proper executable name.  Then:
    ```
    ./check_convergence.sh
    ```
  This outputs 3 sets of convergence data files.  To process these
  into a set of LaTeX table rows:
    ```
    python3 create_pretty_tables.py
    ```

# HSE convergence

This setup can also be used just to test HSE convergence.  

* Build as:
  ```
  make USE_REACT=FALSE -j 20
  ```
  Then run the tests as:
  ```
  ./converge_test_sdc4_nopert.sh
  ```

  We can get the maximum velocity as:
  ```
  fextrema.Linux.gfortran.exe bubble_64_plt00667 | grep -i magvel
  fextrema.Linux.gfortran.exe bubble_128_plt01334 | grep -i magvel
  fextrema.Linux.gfortran.exe bubble_256_plt02667 | grep -i magvel

