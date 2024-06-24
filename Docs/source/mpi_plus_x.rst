.. _ch:mpiplusx:

******************************
Running Options: CPUs and GPUs
******************************

Castro uses MPI for coarse parallelization, distributing boxes across
compute nodes.  For fine-grained parallelism, OpenMP is used for
CPU-based computing and CUDA is used for GPUs.

Running on CPUs
===============

The preferred was of running on CPUs is to use MPI+OpenMP, compiling as::

  USE_MPI=TRUE
  USE_OMP=TRUE

Castro uses tiling to divide boxes into smaller tiles and distributes
these tiles to the OpenMP threads.  This is all managed at the MFIter
level -- no OpenMP directives need to be present in the compute
kernels themselves.  See `MFIter with Tiling
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#sec-basics-mfiter-tiling>`_
for more information.

The optimal number of OpenMP threads depends on the computer
architecture, and some experimentation is needed.  Tiling works best
with larger boxes, so increasing ``amr.max_grid_size`` can benefit
performance.


Running on GPUs
===============

Castro's compute kernels can run on GPUs and this is the preferred way
to run on supercomputers with GPUs.  The exact same compute kernels
are used on GPUs as on CPUs.

.. note::

   Almost all of Castro runs on GPUs, with the main exception being
   the true SDC solver (``USE_TRUE_SDC = TRUE``).

When using GPUs, almost all of the computing is done on the GPUs.  In
the ``MFIter`` loops over boxes, the loops put a single zone on each
GPU thread, to take advantage of the massive parallelism.  The
Microphysics routines (EOS, nuclear reaction networks, etc.) also take
advantage of GPUs, so entire simulations can be run on the GPU.

Best performance is obtained with bigger boxes, so setting
``amr.max_grid_size = 128`` and ``amr.blocking_factor = 32`` can give
good performance.



Castro / AMReX have an option to use managed memory for the GPU --
this means that the data will automatically be migrated from host to
device (and vice versa) as needed, whenever a page fault is
encountered.  This can be enabled via:
``amrex.the_arena_is_managed=1``.

By default, Castro will abort if it runs out of GPU memory.  You can
disable this via ``amrex.abort_on_out_of_gpu_memory=0`` -- together
with running with managed memory, this can allow the memory to be
swapped off of the GPU to make more room available.  This is not
recommended -- oversubscribing the GPU memory will severely impact
performance.

.. index:: castro.hydro_memory_footprint_ratio

The CTU hydrodynamics scheme creates a lot of temporary ``FAB`` 's
when doing the update.  This can lead to the code oversubscribing the
GPU memory during the hydro advance.  To alleviate this, Castro can
break a box into tiles and work on one tile at a time (this is the
approach we use with OpenMP).  In the hydro solver, this is controlled
by ``hydro_tile_size``.  By setting
``castro.hydro_memory_footprint_ratio`` to a number > 0, Castro will
dynamically estimate a good tile size to use for the hydro during the
first timestep and then use this subsequently.  This larger the
number, the more local memory we will allow the hydro solver to use
(generally this means a larger ``hydro_tile_size``).  This can allow
you to run a problem on a smaller number of GPUs if the hydro
temporary memory was the cause of oversubscription.  Current
recommendations are to try ``castro.hydro_memory_footprint_ratio``
between ``2.0`` and ``4.0``.


NVIDIA GPUs
-----------

With NVIDIA GPUs, we use MPI+CUDA, compiled with GCC and the NVIDIA compilers.
To enable this, compile with::

  USE_MPI = TRUE
  USE_OMP = FALSE
  USE_CUDA = TRUE

.. note::

   For recent GPUs, like the NVIDIA RTX 4090, you may need to change
   the default CUDA architecture.  This can be done by adding:

   .. code::

      CUDA_ARCH=89

   to the ``make`` line or ``GNUmakefile``.

.. note::

   CUDA 11.2 and later can do link time optimization.  This can
   increase performance by 10-30% (depending on the application), but
   may greatly increase the compilation time.  This is disabled by
   default.  To enable link time optimization, add:

   .. code::

      CUDA_LTO=TRUE

   to the ``make`` line of ``GNUmakefile``.

AMD GPUs
--------

For AMD GPUs, we use MPI+HIP, compiled with the ROCm compilers.
To enable this, compile with::

  USE_MPI = TRUE
  USE_OMP = FALSE
  USE_HIP = TRUE


Working at Supercomputing Centers
=================================

Our best practices for running any of the AMReX Astrophysics codes
at different supercomputing centers is produced in our workflow
documentation: https://amrex-astro.github.io/workflow/
