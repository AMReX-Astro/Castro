***************************
Offloading a routine to GPU
***************************

In order to offload a routine to the GPU, we insert a ``#pragma gpu``
statement on the line before the call. This tells the preprocessor to
generate a CUDA device version of the function, with all the
additional CUDA GPU management. In addition to this, there are a few
other modifications required to ensure the function operates correctly
on the GPU. We outline these below.

C++
---

- Make sure that the box ``lo`` and ``hi`` are the first two arguments
  and use ``AMREX_INT_ANYD`` as the macro wrapping ``bx.loVect()`` and
  ``bx.hiVect()``

- Likewise, inplace of ``ZFILL()``, use ``AMREX_REAL_ANYD``

- Scalars must be passed by value to Fortran functions (not by
  reference)

- Make sure that you change the variable definitions in the Fortran
  code to include value so that the Fortran knows the variables are
  being passed by value rather than by reference. If you don’t do
  this, the code is liable to segfault

- FArrayBox functions like `setVal()` and `saxpy()` are not on the
  GPU, so you should explicitly write out these operations in C++.
  The MultiFab counterparts are on the GPU.

- Use ``BL_TO_FORTRAN_ANYD()`` to wrap MultiFab arguments (and in the header file wrap with ``BL_FORT_FAB_ARG_3D``)

To illustrate these modifications, consider the function ``divu``. To call the function on the CPU, we would write

.. code-block:: c++

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& obx = mfi.growntilebox(1);

        divu(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             ZFILL(dx),
             BL_TO_FORTRAN_ANYD(div[mfi]));
     }

Implementing the changes described above to offload this to GPU, this becomes

.. code-block:: c++

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& obx = mfi.growntilebox(1);

    #pragma gpu
        divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             AMREX_REAL_ANYD(dx),
             BL_TO_FORTRAN_ANYD(div[mfi]));
     }

Fortran
-------

- the routine must only operate on a single zone in ``lo:hi``.

  - even for temporary arrays, you cannot write to an ``i+1`` zone
    (e.g. when doing limiting) -- this will cause a race condition.
    If necessary do extra computation to avoid the temporaries

- mark the routine with ``!$gpu``

  - If you get the error ``nvlink error: Undefined reference to
    'foo_device' in 'tmp_build_dir/o/3d.pgi.CUDA.EXE/bar.o'``, then
    you’ve probably forgotten to do this.

- If a module defines its own variables, these variables need to be
  ``allocatable`` (even if they are scalars) and marked as
  ``attributes(managed)``. Additional routines may be needed to
  allocate these variables before they’re used and deallocate them
  when they’re no longer needed.

  - Examples of this can be seen for the sponge parameters. These are
    marked as ``allocatable`` and ``attributes(managed)`` in
    ``sponge_nd.F90``, and therefore must be allocated (by
    ``ca_allocate_sponge_params``) and deallocated (by
    ``ca_deallocate_sponge_params``) at the beginning and end of the

- Temporary variables must be defined outside of function calls. E.g. if a
  function call contains ``foo(x(a:b)/y)``, you need to define a new variable
  ``z = x(a:b)/y`` then pass this into the function as ``foo(z)``.

  - If you don’t do this, you may see the error ``Array reshaping is
    not supported for device subprogram calls``

- If importing a function from another module, make sure to put the
  import within the function/subroutine, and put ``! function`` at the
  end of the line, e.g.

.. code-block:: fortran

   use my_module, only: my_func ! function

- Individual functions should be imported individually (so not ``use
  my_module, only: func1, func2 ! function``) and there must be a
  space either side of the ``!``

- Make sure the fortran file is ``.F90`` rather than ``.f90`` (and
  remember to update the ``Make.xx`` file to reflect this). If you
  don’t do this you will see the error ``Label field of continuation
  line is not blank``

  - This is required as we use the convention that ``.F90`` files are
    processed by the preprocessor, ``.f90`` files are not. The
    preprocessor will therefore only generate the required device
    function if the file has the correct extension.

We can see some of the above modifications by looking at the
subroutine ``derangmomx`` in ``Derive_nd.F90``:

.. code-block:: fortran

   subroutine derangmomx(L,L_lo,L_hi,ncomp_L, &
                         u,u_lo,u_hi,ncomp_u, &
                         lo,hi,domlo,domhi, &
                         dx,xlo) bind(C, name="derangmomx")

      use amrex_constants_module, only: HALF
      use math_module, only: cross_product ! function
      use amrex_fort_module, only : rt => amrex_real
      use prob_params_module, only: center

      implicit none

      integer, intent(in), value :: ncomp_L, ncomp_u
      integer, intent(in) :: L_lo(3), L_hi(3)
      integer, intent(in) :: u_lo(3), u_hi(3)
      integer, intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
      real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
      real(rt), intent(in) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
      real(rt), intent(in) :: dx(3), xlo(3)

      integer          :: i, j, k
      real(rt)         :: loc(3), mom(3), ang_mom(3), rho

      !$gpu

      do k = lo(3), hi(3)
         loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)

      ...

   end subroutine derangmomx

- Here, we can see that the ``cross-product`` function from the
  ``math_module`` is marked as ``! function``, which tells the
  preprocessor to generate a device version of this function.

- The scalars ``ncomp_L`` and ``ncomp_u`` are both passed in by value.

- The ``!$gpu`` directive has been inserted after the definition of
  all the variables passed into the routine and all the local
  variables, but before the main body of the function.

- The routine only operates on values in a single zone of ``lo:hi``.


.. To be documented
.. ----------------
..
.. when do we need to mark stuff as attributes(managed)?


To check if we launched a kernel
--------------------------------

Run

.. code-block:: sh

   nvprof ./Castro.xxx



How to debug
------------

- Run under ``cuda-memcheck``

- Run under ``cuda-gdb``

- Turn off GPU offloading for some part of the code with

.. code-block:: c++

    Gpu::setLaunchRegion(0);
    ... ;
    Gpu::setLaunchRegion(1);

- Run with ``CUDA_LAUNCH_BLOCKING=1``.  This means that only one
  kernel will run at a time.  This can help identify if there are race
  conditions.

