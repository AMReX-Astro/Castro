****************
Performance Tips
****************

Getting good performance is a mix of single-processor performance and
parallel efficiency.

Parallel Efficiency
===================

Load Balancing
--------------

Good parallel performance occurs when there are the same number of boxes
on each MPI task.  If too many MPI tasks are used, then there will not be
enough work to go around.

Use a large ``amr.max_grid_size`` and aim for one box per MPI task one
each level.

Usually it's a good idea to do a small scaling study before running a
large simulation to find the optimal number of nodes to run your
problem.

Gravity / Multigrid
-------------------

.. index:: amr.blocking_factor

Multigrid performance is best when the grid can be coarsened down to
the smallest possible size.  This means that the number of zones in
each dimension of a box should be a power-of-two.  This can be controlled
by ``amr.blocking_factor``.  See also the `AMReX docs on linear solvers <https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#>`_.

.. index:: amr.subcycling_mode

It can also be faster to run without subcycling in some instances,
which is controlled by ``amr.subcycling_mode``, as described
in :ref:`sec:subcycling`.


.. index:: gravity.max_multipole_order

For isolated systems, we use a multipole expansion for constructing
the Dirichlet boundary conditions for the Poisson solve.  The
construction of the multipole expansion can be expensive.  If the
system is roughly spherical and the boundaries are far from the mass,
then you probably can get away with fewer multipole moments.  For
reference, for the white dwarf mergers in :cite:t:`katz:2016`,
$l_\mathrm{max} = 6$ was used.  The maximum multipole moment is
set by ``gravity.max_multipole_order``.


Global Diagnostics
------------------

.. index:: castro.sum_interval

:ref:`sec:global_diag` are output at regular intervals.  These can
require reduction operations across all processors, which can be
expensive for big calculations.  Often we don't need the diagnostics
every step, so setting ``castro.sum_interval`` to a larger value (like
``10``) could improve performance.

Metrics
=======

To understand where the simulation is spending the most time, it is suggested
to run using the `AMReX Tiny Profiler <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools.html#tiny-profiling>`_.

Simply build as:

.. prompt:: bash

   make TINY_PROFILE=TRUE

and then run your simulation (any number of processors, CPU or GPU).  Typically
running for 10--100 steps is enough to get a sense of where the time is spent.
Upon completion, a detailed report will be output to stdout showing which
functions / kernels are using the most time.

This can help you determine where to focus your optimization efforts.

Alternately, you can use the GNU Profiler on CPUs to get line-by-line
performance details.  This can be enabled by building with

.. prompt:: bash

   make USE_GPROF=TRUE

running, and then using the ``gprof`` command line tool.


GPUs
====

For GPU performance, see the discussion in :ref:`sec:running_on_gpus`.

The main parameter to explore when running on GPUs is
``castro.hydro_memory_footprint_ratio``, which can use tiling in the
hydrodynamics solver to prevent oversubscription of GPU memory.



Reactions
=========

Reactions are often the most time-consuming part of a simulation.  The
following are some things to try to improve the performance:

* Try using the Runge-Kutta-Chebyshev integrator (see the Microphysics
  `ODE integrators docs <https://amrex-astro.github.io/Microphysics/docs/ode_integrators.html>`_.

  This is an explicit integrator that can work with moderately-stiff
  networks.  Experience shows that it can work well with flames
  (sometimes being twice as fast as the VODE integrator), but probably
  not very efficiently with detonations.

* Use the analytic Jacobian (selected via ``integrator.jacobian=1``).

  The analytic Jacobian is faster to evaluate than the
  difference-approximation.  Note that if the integration fails in a zone,
  by default Castro will record this failure and reject the step
  triggering the Castro :ref:`ch:retry`.  This can be expensive, so
  something to try instead is to catch the failure in the burner and
  retry the burn with the alternate Jacobian (e.g., numerical differencing),
  since sometimes that helps the integrator get through.  This can
  be enabled via:

  ::

    integrator.use_burn_retry = 1
    integrator.retry_swap_jacobian = 1

  sometimes it can be advantageous to force a burn failure if the first
  pass it taking too many integration steps.  This can be done by setting
  ``integrator.ode_max_steps`` to a small value (like ``5000`` instead of
  the default ``150000``).

* Try a single-precision Jacobian.

  On GPUs, the Jacobian can consume a lot of memory, and the linear algebra
  solve is often the most expensive part of the ODE integration.  Since the
  role of the Jacobian is simply to point the nonlinear solver in the right
  direction for computing the correction, we can sometimes benefit from using
  a single-precision Jacobian.  This needs to be set at compile time
  by building as:

  .. prompt:: bash

     make USE_SINGLE_PRECISION_JACOBIAN=TRUE

* Disable burning where it is not needed.

  The parameters ``castro.react_rho_min``, ``castro.react_rho_max``,
  ``castro.react_T_min``, and ``castro.react_T_max`` can be used to
  control the density and temperature where we burn.  Disabling
  reactions in low density regions (where you don't expect appreciable
  energy generation) can help with performance.

* Experiment with SDC

  For explosive flows, the simplified-SDC solver can be more
  efficient, since it eliminates some stiffness from the system.  See
  :ref:`sec:flowchart` for information on how to enable the SDC
  solver.

* Use self-consistent NSE

  For explosive flows, when the temperature reaches 5 GK or more, the
  nuclei can enter nuclear statistical equilibrium.  The cancellation
  of the forward and reverse rates can be challenging for the reaction
  integrator to handle.  Instead, we can detect that we are entering
  NSE and use the NSE solution instead of integration in these cases.
  This is described in the `Microphysics self-consistent NSE docs
  <https://amrex-astro.github.io/Microphysics/docs/nse_net.html>`_.
