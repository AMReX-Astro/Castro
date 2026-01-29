# hse_convergence

This is meant to be a simple 1-d test for assessing the convergence of
hydro + gravity in maintaining HSE.

A simple hydrostatic, isentropic atmosphere is setup (using the
Helmholtz stellar equation of state).  For 4th order, the initial
model is integrated to be 4th-order accurate at cell-centers.
The base conditions are ρ = 1.e7 g/cm**3 and T = 1.e8 K.

By default, HSE boundary conditions are used at the top and bottom.

Convergence can be measure either with the RichardsonConvergenceTest
(from AMReX) tool or by looking at the max |U| in the plotfiles.

The following describes the results for various solvers.

## CTU + PLM

Build this simply as `make`, and run the suite via:

```
./convergence_plm.sh
```

This explores

  * the default HSE BCs and `use_pslope` -- results are in
    `plm.converge.out`

    The maximum velocity for different zoning is:

    ```
     64   1016837.3104
    128    512627.30938
    256    257385.3507
    512    128964.94588
    ```

    This looks first-order.

  * the HSE BCs with reflection and `use_pslope` -- results are in
    `plm-hsereflect.converge.out`

    The maximum velocity for different zoning is:

    ```
     64      4814.029903
    128      1214.1915455
    256       304.83736032
    512        76.376611631
    ```

    This looks second-order.

  * reflecting BCs (instead of HSE BCs) without `use_pslope` --
    results are in `plm-reflect-nopslope.converge.out`

    The maximum velocity for different zoning is:

    ```
     64   1324129.692
    128    663443.58475
    256    332063.12615
    512    166116.20661
    ```

    This looks first-order.

  * reflect BCs with `use_pslope` -- results are in
    `plm-reflect-pslope.converge.out`

    The maximum velocity for different zoning is:

    ```
     64       639.41980735
    128       154.6195799
    256        37.100875875
    512         8.898682755
    ```

    This looks second-order.

These tests show that the best results (by far) come from
`use_pslope=1` and reflecting BCs

## CTU + PPM

Again, build simply as `make` and run the suite via:

```
./convergence_ppm.sh
```

This explores a similar set of conditions as the PLM:

  * the default HSE BCs (no pslope) -- results are in
    `ppm.converge.out`

    The maximum velocity for different zoning is:

    ```
     64    681956.19954
    128    341814.31565
    256    171107.26187
    512     85605.18186
    ```

    This looks first-order.

  * the default HSE BCs but with velocity reflection
    and temperature interpolation -- results are in
    `ppm-hsereflect.converge.out`

    The maximum velocity for different zoning is:

    ```
     64      2997.6612833
    128       757.50207433
    256       190.3748635
    512        47.721329387
    ```

    This looks second-order.

  * reflecting BCs -- results are in `ppm-reflect.converge.out`

    The maximum velocity for different zoning is:

    ```
     64   1333801.695
    128    669021.94809
    256    334511.29859
    512    167523.07073
    ```

    This looks first-order.

  * reflecting BCs + `use_pslope` -- results are in
    `ppm-reflect-pslope.converge.out`

    The maximum velocity for different zoning is:

    ```
     64       566.56312924
    128       145.4205635
    256        36.291205477
    512         8.9841497662
    ```

    This looks second-order.

  * reflecting BCs + `ppm_well_balanced` -- results are in
    `ppm-reflect-wellbalanced.converge.out`

    The maximum velocity for different zoning is:

    ```
     64         0.00048977043343
    128         0.001728828207
    256         0.00031968121456
    512         2.5713543173e-05
    ```

    This looks better than second-order.

These tests show that the best results (by far) come from reflecting
BCs with `ppm_well_balanced=1`.  There is no equivalent for this in
the PLM solver.

## SDC

This uses the true SDC solver, looking at both 2nd- and 4th-order
solvers.

Build as:

```
make USE_TRUE_SDC=TRUE
```

run the suite as:

```
./convergence_sdc.sh
```

The following variations are explored:

  * SDC-2 + PLM with `use_pslope` and reflecting BCs -- results are in
    `sdc2-reflect.converge.out`

    The maximum velocity for different zoning is:

    ```
     64         0.00048049855069
    128         0.0017153949315
    256         0.0003098853472
    512         2.4050544423e-05
    ```

    This looks better than second-order, and almost the same as
    CTU + PPM with well-balancing.

  * SDC-2 + PPM and reflecting BCs -- results are in
    `sdc2-ppm-reflect.converge.out`

    The maximum velocity for different zoning is:

    ```
     64   1970628.8788
    128    987751.15483
    256    494405.42935
    512    247328.03663
    ```

    This looks first-order.

  * SDC-2 + PLM with `use_pslope` and HSE BCs -- results are in
    `sdc2.converge.out`

    The maximum velocity for different zoning is:

    ```
     64     39132.815939
    128      9770.7837714
    256      2441.4880416
    512       610.27386207
    ```

    This looks second-order.

  * SDC-2 + PPM with HSE BCs -- results are in `sdc2-ppm.converge.out`

    The maximum velocity for different zoning is:

    ```
     64      8004.7751989
    128      9770.7837714
    256      2441.4880416
    512       610.27386207
    ```

    This looks second-order.

  * SDC-4 with reflecting BCs -- results are in
    `sdc4-reflect.converge.out`

    The maximum velocity for different zoning is:

    ```
     64         0.018841266999
    128         0.0011894966175
    256         8.7006924189e-05
    512         1.2416964741e-05
    ```

    This looks like it gets down to roundoff, but converges almost 4th
    order until then.


These tests show that the PLM + reflect (which uses the well-balanced
`use_pslope`) and the SDC-4 + reflect give the lowest errors and
expected (or better) convergence.
