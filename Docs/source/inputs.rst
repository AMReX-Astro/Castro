***********
Input Files
***********

The Castro executable uses an inputs file at runtime to set and
alter the behavior of the algorithm and initial conditions.

Runtime parameters take the form ``namespace.parameter = value`` ,
where *namespace* identifies the major code component.  Some
parameters take multiple values separated by spaces.

Typically
named ``inputs``, it

  * Sets the AMReX parameters for gridding, refinement, etc., through the
    ``geometry`` and ``amr`` namespaces.

  * Enables different physics behaviors through the ``castro`` namespace.

  * Sets the problem-specific runtime parameters through the ``problem`` namespace.

  * Sets any Microphysics runtime parameters, through the various namespaces
    defined in Microphysics (like ``eos``, ``integrator``, ``network``, ...).

.. warning:: Because the inputs file is handled by the C++ portion
   of the code, any quantities you specify in scientific notation,
   must take the form ``1.e5`` and not ``1.d5``—the ``d``
   specifier is not recognized.


.. note::

   Additionally, note that in Castro, all quantities are in CGS units.


Common inputs Options
=====================


Problem Geometry
----------------

The ``geometry`` namespace is used by AMReX to define the
computational domain. The main parameters here are:

  * ``geometry.prob_lo``: physical location of low corner of the
    domain (type: ``Real``; must be set)

    Note: a number is needed for each dimension in the problem.

  * ``geometry.prob_hi``: physical location of high corner of the
    domain (type: ``Real``; must be set)

    Note: a number is needed for each dimension in the problem.

  * ``geometry.coord_sys``: coordinate system, 0 = Cartesian,
    1 = :math:`r`-:math:`z` (2-d only), 2 = spherical (1-d only) (must be set)

  * ``geometry.is_periodic``: is the domain periodic in this direction?
    0 if false, 1 if true (default: ``0 0 0``)

    Note: an integer is needed for each dimension in the problem.

  * ``castro.center``: physical location of problem center on the
    domain (type: ``Real``; default: ``0.0 0.0 0.0``). The problem
    center is used for gravity, rotation, and some other quantities.
    This is not necessarily the geometric center of the domain—often
    you should choose it to coincide with the center of mass of your
    system.

    Note: a number is needed for each dimension in the problem.

As an example, the following::

    geometry.prob_lo = 0 0 0
    geometry.prob_hi = 1.e8 2.e8 2.e8
    geometry.coord_sys = 0
    geometry.is_periodic = 0 1 0
    castro.center = 5.e7 1.e8 1.e8

This defines the domain to run from :math:`(0,0,0)` at the lower left to
:math:`(10^8,\, 2\times 10^8,\, 2\times 10^8)` at the upper right in physical
space, specifies a Cartesian geometry, and makes the domain periodic
in the :math:`y`-direction only. The problem center is set to be halfway in
between the lower left and upper right corners.

Domain Boundary Conditions
--------------------------

Boundary conditions are specified using integer keys that are interpreted
by AMReX. The runtime parameters that we use are:

  * ``castro.lo_bc``: boundary type of each low face (must be set)

  * ``castro.hi_bc``: boundary type of each high face (must be set)

The valid boundary types are:

.. table:: boundary condition types
   :align: center

   +------------------------+-----------------+
   | 0: Interior / Periodic | 3: Symmetry     |
   +------------------------+-----------------+
   | 1: Inflow              | 4: Slip Wall    |
   +------------------------+-----------------+
   | 2: Outflow             | 5: No Slip Wall |
   +------------------------+-----------------+

.. note:: ``castro.lo_bc`` and ``castro.hi_bc`` must be consistent
   with ``geometry.is_periodic``—if the domain is periodic in a
   particular direction then the low and high bc’s must be set to 0
   for that direction.

As an example, the following::

    castro.lo_bc = 1 4 0
    castro.hi_bc = 2 4 0

    geometry.is_periodic = 0 0 1

This defines a problem with inflow (1) in the low-\ :math:`x` direction,
outflow (2) in the high-\ :math:`x` direction, slip wall (4) on
the low and high :math:`y`-faces, and periodic in the :math:`z`-direction.
See § :ref:`soft:phys_bcs`.

Resolution
----------

The grid resolution is specified by defining the resolution at the
coarsest level (level 0) and the number of refinement levels and
factor of refinement between levels. The relevant parameters are:

  * ``amr.n_cell``: number of cells in each direction at the coarsest
    level (integer :math:`> 0`; must be set)

  * ``amr.max_level``: number of levels of refinement above the
    coarsest level (integer :math:`\geq 0`; must be set)

  * ``amr.ref_ratio``: ratio of coarse to fine grid spacing
    between subsequent levels (2 or 4; must be set)

  * ``amr.regrid_int``: how often (in terms of number of steps) to
    regrid (integer; must be set)

  * ``amr.regrid_on_restart``: should we regrid immediately after
    restarting? (0 or 1; default: 0)

.. note:: if ``amr.max_level = 0`` then you do not need to set
   ``amr.ref_ratio`` or ``amr.regrid_int``.

Some examples::

    amr.n_cell = 32 64 64

would define the domain to have 32 cells in the :math:`x`-direction, 64 cells
in the :math:`y`-direction, and 64 cells in the :math:`z`-direction *at the
coarsest level*. (If this line appears in a 2D inputs file then the
final number will be ignored.)

::

    amr.max_level = 2

would allow a maximum of 2 refined levels in addition to the coarse
level. Note that these additional levels will only be created only if
the tagging criteria are such that cells are flagged as needing
refinement. The number of refined levels in a calculation must be
:math:`\leq` ``amr.max_level``, but can change in time and need not
always be equal to ``amr.max_level``.

::

    amr.ref_ratio = 2 4

would set factor of 2 refinement between levels 0 and 1, and factor of 4
refinement between levels 1 and 2. Note that you must have at least
``amr.max_level`` values of ``amr.ref_ratio`` (Additional values
may appear in that line and they will be ignored).

::

    amr.regrid_int = 2 2

tells the code to regrid every 2 steps. Thus in this example, new
level 1 grids will be created every 2 level-0 time steps, and new
level 2 grids will be created every 2 level-1 time steps. If
``amr.regrid_int`` :math:`<` 0 for any level, then regridding starting at that
level will be disabled. If ``amr.regrid_int = -1`` only, then we
never regrid for any level. Note that this is not compatible with
``amr.regrid_on_restart = 1``.


Other parameters
----------------

There are a large number of solver-specific runtime parameters. We describe these
together with the discussion of the physics solvers in later chapters.
