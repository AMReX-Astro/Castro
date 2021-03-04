# bubble_convergence

This problem was developed to test the convergence of the 4th order
SDC solver with gravity, reactions, and reflecting boundary
conditions.

This was featured in the Castro SDC paper:
https://ui.adsabs.harvard.edu/abs/2019ApJ...886..105Z


## Measuring convergence

The script `converge_test.sh` will test the convergence for the bubble
rising problem by running tests at 5 different resolutions (the 32**2
resolution was not shown in the paper).

* Build as:
  ```
  make USE_TRUE_SDC=TRUE -j 20
  ```

* Run the tests via: 
  ```
  ./converge_test.sh
   ```

  This relies on having built the RichardsonConvergenceTest too in
  `amrex/Tools/C_util/Convergence/` and putting the executable in your
  path.

  3 output files will be produced: `sdc_converge.lo.out`,
  `sdc_converge.mid.out`, and `sdc_converge.hi.out`, giving the order
  of convergence for the 32-64-128, 64-128-256, and 128-256-512 runs
  respectively.

* Create a summary table using `create_pretty_tables.py`.  E.g., for
  the mid and hi resolution cases, you would do:
  ```
  python create_pretty_tables.py --simple sdc_converge.mid.out sdc_converge.hi.out
  ```

  This gives:
  ```
  $\rho$                           3.590942e+15   3.264       3.739197e+14   3.713       2.852085e+13
  $\rho u$                          1.11992e+24   3.794       8.071755e+22   3.930       5.296302e+21
  $\rho v$                          1.31454e+24   3.544       1.127087e+23   3.839       7.878383e+21
  $\rho E$                         3.701644e+32   2.947       4.801223e+31   3.646       3.834012e+30
  $\rho e$                         3.701199e+32   2.947       4.800721e+31   3.646       3.833886e+30
  $T$                              1.438083e+18   3.508       1.264438e+17   3.829       8.898941e+15
  $\rho X(\isotm{He}{4})$          3.589696e+15   3.266       3.732104e+14   3.711       2.849884e+13
  $\rho X(\isotm{C}{12})$          1.519871e+13   2.544       2.605906e+12   3.797       1.874147e+11
  $\rho X(\isotm{O}{16})$              35896250   3.262            3741870   3.714           285088.2
  $\rho X(\isotm{Fe}{56})$             35908410   3.264            3739049   3.713           285208.5
  ```

  which is essentially the same values shown in table 11 of the Castro SDC paper.

## HSE convergence

This setup can also be used just to test HSE convergence.  

* Build as:
  ```
  make USE_TRUE_SDC=TRUE USE_REACT=FALSE -j 20
  ```
  Then run the tests as:
  ```
  ./converge_test_sdc4_nopert.sh
  ```

  We can get the maximum velocity as:
  ```
  fextrema.gnu.ex bubble_64_plt00667 | grep -i magvel
  fextrema.gnu.ex bubble_128_plt01334 | grep -i magvel
  fextrema.gnu.ex bubble_256_plt02667 | grep -i magvel
  ```

  This gives (for max |U|):
  ```
   64: 0.018842065993
  128: 0.0011912049874
  256: 8.8777896776e-05
  ```
  demonstrating nearly 4th order convergence of the HSE state.


## Original Cori convergence testing

The results in the Castro SDC paper were run on Cori, using
the file in `job_scripts/`.

* in `Castro/Exec/reacting_tests/bubble_convergence`
    ```
    make USE_TRUE_SDC=TRUE COMP=intel -j 4
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
