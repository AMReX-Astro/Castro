.. _ch:mpiplusx:

******************************
Running Options: CPUs and GPUs
******************************

Castro uses MPI for coarse parallelization, distributing boxs across
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
to run on supercomputers with GPUs.  At the moment, offloading is
handled using CUDA and managed memory.  The exact same compute kernels
are used on GPUs as on CPUs.

.. note::

   Almost all of Castro runs on GPUs, with the main exception being
   the true SDC solver (``USE_TRUE_SDC = TRUE``).

To enable GPU computing, compile with::

  USE_MPI = TRUE
  USE_OMP = FALSE
  USE_CUDA = TRUE

Currently, we support the PGI compilers.

When using GPUs, almost all of the computing is done on the GPUs.  In
the MFIter loops over boxes, the loops put a single zone on each GPU
thread, to take advantage of the massive parallelism.  The Microphysics
in StarKiller also takes advantage of GPUs, so entire simulations can
be run on the GPU.

Best performance is obtained with bigger boxes, so setting
``amr.max_grid_size = 128`` and ``amr.blocking_factor = 32`` can give
good performance.


Working at Supercomputing Centers
=================================

Our best practices for running any of the AMReX Astrophysics codes
at different supercomputing centers is produced in our workflow
documentation: https://amrex-astro.github.io/workflow/

