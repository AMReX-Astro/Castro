.. _ch:sdc:

*****************************
Spectral Deferred Corrections
*****************************

The Castro SDC solver couples the hydrodynamics tightly together,
iteratively improving the convergence of the solution.  This is the
basis of the 4th order accurate Castro solver.  The algorithm is described
in :cite:`castro-sdc`.

.. note::

   Here we are referring to the full SDC time integration scheme
   (``castro.time_integration_method = 2``), not the simplified-SDC solver.


The options that describe the quadrature and iterations are:

* ``castro.sdc_order`` : the desired spatial and temporal order.  2 and 4 are supported.

* ``castro.sdc_quadrature`` : the quadrature scheme used for the
  time-integration.  This determines the number and location of the
  temporal nodes.  Supported values are 0 for Gauss-Lobatto and 1 for
  Radau IIA.

* ``castro.sdc_extra`` : the number of extra iterations to take.  By
  default the number of iterations used is equal to the value of
  ``sdc_order``.


The options that affect the nonlinear solve are:

* ``sdc_solver`` : the method we use to do the nonlinear solution of
  the reaction system.  Values are:

  * 1 : pure Newton iteration (we subdivide the time interval if
    needed to get the Newton method to converge).

  * 2 : use VODE to solve the nonlinear system by expressing it as an ODE system.

  * 3 : use VODE for the first iteration and then Newton for the
    subsequent iterations.

  The tolerances for both the Newton and VODE solver are controlled by
  the usual Microphysics parameters: ``integrator.rtol_spec``,
  ``integrator.atol_spec``, ``integrator.rtol_enuc``,
  ``integrator.atol_enuc``.

  In all cases, the type of Jacobian (analytic or numerical) is determined by
  ``integrator.jacobian``.






