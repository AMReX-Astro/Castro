*********************
GPU Programming Model
*********************

CPUs and GPUs have separate memory, which means that working on both
the host and device may involve managing the transfer of data between
the memory on the host and that on the GPU.

In Castro, the core design when running on GPUs is that all of the compute
should be done on the GPU.

When we compile with ``USE_CUDA=TRUE`` or ``USE_HIP=TRUE``, AMReX will allocate
a pool of memory on the GPUs and all of the ``StateData`` will be stored there.
As long as we then do all of the computation on the GPUs, then we don't need
to manage any of the data movement manually.

.. note::

   We can tell AMReX to allocate the data using managed-memory by
   setting:

   ::

      amrex.the_arena_is_managed = 1

   This is generally not needed.

The programming model used throughout Castro is C++-lambda-capturing
by value.  We access the ``FArrayBox`` stored in the ``StateData``
``MultiFab`` by creating an ``Array4`` object.  The ``Array4`` does
not directly store a copy of the data, but instead has a pointer to
the data in the ``FArrayBox``.  When we capture the ``Array4`` by
value in the GPU kernel, the GPU gets access to the pointer to the
underlying data.


Most AMReX functions will work on the data directly on the GPU (like
``.setVal()``).

In rare instances where we might need to operate on the data on the
host, we can force a copy to the host, do the work, and then copy
back.  For an example, see the reduction done in  ``Gravity.cpp``.

.. note::

   For a thorough discussion of how the AMReX GPU offloading works
   see :cite:`amrex-ecp`.


Runtime parameters
------------------

The main exception for all data being on the GPUs all the time are the
runtime parameters.  At the moment, these are allocated as managed
memory and stored in global memory.  This is simply to make it easier
to read them in and initialize them on the CPU at runtime.


