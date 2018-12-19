***************************
Offloading a routine to GPU
***************************

In order to offload a routine to the GPU, we insert a ``#pragma gpu`` statement on the line before the call. This tells the preprocessor to generate a CUDA device version of the function, with all the additional CUDA GPU management. In addition to this, there are a few additional modifications required to ensure the function operates correctly on the GPU. We outline these below.

C++
---

- Make sure that the box ``lo`` and ``hi`` are the first two arguments and use `AMREX_INT_ANYD` as the macro wrapping ``bx.loVect()`` and ``bx.hiVect()``
- likewise, inplace of ``ZFILL()``, use ``AMREX_REAL_ANYD``
- Scalars must be passed by value to Fortran functions (not by reference)
- Make sure that you change the variable definitions in the Fortran code to include value so that the Fortran knows the variables are being passed by value rather than by reference. If you don’t do this, the code is liable to segfault
- Use ``BL_TO_FORTRAN_ANYD()`` to wrap MultiFab arguments (and in the header file wrap with ``BL_FORT_FAB_ARG_3D``)

Fortran
-------

- the routine must only operate on a single zone in ``lo:hi``.
- mark the routine with ``!$gpu``
- If you get the error ``nvlink error: Undefined reference to 'foo_device' in 'tmp_build_dir/o/3d.pgi.CUDA.EXE/bar.o'``, then you’ve probably forgotten to do this
even for temporary arrays, you cannot write to an ``i+1`` zone (e.g. when doing limiting) -- this will cause a race condition.  If necessary do extra computation to avoid the temporaries
- If a module defines its own variables, these variables need to be ``allocatable`` (even if they are scalars) and marked as ``attributes(managed)`` (see ``meth_params.F90`` for examples)
    * Additional routines may be needed to allocate these variables before they’re used and deallocate them when they’re no longer needed
- Temporary variables must be defined outside of function calls
    * E.g. if a function calls ``contain foo(x(a:b)/y)``, you need to define a new variable ``z = x(a:b)/y`` then pass this into the function as ``foo(z)``
    * If you don’t do this, may see the error ``Array reshaping is not supported for device subprogram calls``
- If importing a function from another module, make sure to put the import within the function/subroutine, and put ``! function`` at the end of the line, e.g.
```
use my_module, only: my_func ! function
```
- Individual functions should be imported individually (so not ``use my_module, only: func1, func2 ! function``) and there must be a space either side of the `!`
- Make sure the fortran file is ``.F90`` rather than ``.f90`` (and remember to update the ``Make.xx`` file to reflect this).
    * If you don’t do this you will see the error ``Label field of continuation line is not blank``


To be documented
----------------

when do we need to mark stuff as attributes(managed)?


To check if we launched a kernel
--------------------------------

Run
```
nvprof ./Castro.xxx
```


How to debug
------------

- Run under ``cuda-memcheck``
- Run under ``cuda-gdb``
- Turn off GPU offloading for some part of the code with
```
Device::endDeviceLaunchRegion();
... ;
Device::beginDeviceLaunchRegion();
```
