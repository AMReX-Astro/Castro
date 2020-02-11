**********
Regridding
**********

The details of the regridding strategy are described in
§ :ref:`sec:tagging`; here we cover how the input parameters can
control the gridding.

As described later, the user defines Fortran subroutines which tag
individual cells at a given level if they need refinement. This list
of tagged cells is sent to a grid generation routine, which uses the
Berger-Rigoutsos algorithm :cite:`br-refine` to create rectangular
grids that contain the tagged cells.

The relevant runtime parameters are:

  * ``amr.regrid_file``: name of file from which to read the grids
    (text; default: no file)

    If set to a filename, e.g. ``fixed_girds``, then list of grids at
    each fine level are read in from this file during the gridding
    procedure. These grids must not violate the ``amr.max_grid_size``
    criterion. The rest of the gridding procedure described below will
    not occur if ``amr.regrid_file`` is set.

  * ``amr.n_error_buf``: radius of additional tagging
    around already tagged cells (integer :math:`\geq 0`; default: 1)

  * ``amr.max_grid_size``: maximum size of a grid in any
    direction (integer :math:`> 0`; default: 128 (2-d), 32 (3-d))

    Note: ``amr.max_grid_size`` must be even, and a multiple of
    ``amr.blocking_factor`` at every level.

  * ``amr.blocking_factor``: grid size must be a multiple of this
    (integer :math:`> 0`; default: 2)
    ``amr.blocking_factor`` at every level must be a power of 2
    and the domain size must be a multiple of ``amr.blocking_factor``
    at level 0.

    .. note:: This can be very important for elliptic problems with
       multigrid. A higher blocking factor allows the multigrid
       algorithm to coarsen more at the lowest level, reducing the
       amount of work required by the bottom solver.

  * ``amr.grid_eff``: grid efficiency (Real :math:`>0` and :math:`<1`;
    default: 0.7)

    When creating a refined grid, do we make boxes that only include
    the coarse cells that were explicitly tagged for refinement? or do
    we allow ourselves to encompass nearby, untagged cells in order to
    make larger and more regular boxes? This is the grid efficiency.

    When ``blocking_factor = 1``, *grid efficiency* is exactly the
    fraction of refined cells in the fine ``BoxArray`` which
    correspond to coarse cells which were tagged. For other blocking
    factors, we actually apply ``grid_eff`` at the level which has been
    coarsened by ``blocking_factor``, so it is no longer strictly this
    fraction, but the idea is still the same.

  * ``amr.refine_grid_layout``: refine grids more if # of
    processors :math:`>` # of grids (0 if false, 1 if true; default: 1)

Note also that ``amr.n_error_buf``, ``amr.max_grid_size`` and
``amr.blocking_factor`` can be read in as a single value which is
assigned to every level, or as multiple values, one for each level.

As an example, consider::

    amr.grid_eff = 0.9
    amr.max_grid_size = 64
    amr.blocking_factor} = 32

The grid efficiency, ``amr.grid_eff``, means that during the grid
creation process, at least 90% of the cells in each grid at the level
at which the grid creation occurs must be tagged cells. A higher
grid efficiency means fewer cells at higher levels, but may result
in the production of lots of small grids, which have inefficient cache
and OpenMP performance and higher communication costs.

The ``amr.max_grid_size`` parameter means that the final grids will be
no longer than 64 cells on a side at every level.  Alternately, we
could specify a value for each level of refinement as
``amr.max_grid_size = 64 32 16`` in which case our final grids will be
no longer than 64 cells on a side at level 0, 32 cells on a side at
level 1, and 16 cells on a side at level 2. The
``amr.blocking_factor`` means that all of the final grids will be
multiples of 32 at all levels.  Again, this can be specified on a
level-by-level basis, like ``amr.blocking_factor = 32 16 8``, in which
case the dimensions of all the final grids will be multiples of 32 at
level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

Getting good performance
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters can have a large impact on the performance
of Castro, so taking the time to experiment with is worth the effort.
Having grids that are large enough to coarsen multiple levels in a
V-cycle is essential for good multigrid performance in simulations
that use self-gravity.



How grids are created
~~~~~~~~~~~~~~~~~~~~~

The gridding algorithm proceeds in this order:

#. Grids are created using the Berger-Rigoutsos clustering algorithm
   modified to ensure that all new fine grids are divisible by
   ``amr.blocking_factor``.

#. Next, the grid list is chopped up if any grids are larger than
   ``max_grid_size``.  Note that because ``amr.max_grid_size`` is a
   multiple of ``amr.blocking_factor`` the ``amr.blocking_factor``
   criterion is still satisfied.

#. Next, if ``amr.refine_grid_layout = 1`` and there are more
   processors than grids, and if ``amr.max_grid_size`` / 2 is a
   multiple of ``amr.blocking_factor``, then the grids will be
   redefined, at each level independently, so that the maximum length
   of a grid at level :math:`\ell`, in any dimension, is
   ``amr.max_grid_size`` [:math:`\ell`] / 2.

#. Finally, if ``amr.refine_grid_layout = 1``, and there are still
   more processors than grids, and if ``amr.max_grid_size`` / 4 is a
   multiple of ``amr.blocking_factor``, then the grids will be
   redefined, at each level independently, so that the maximum length
   of a grid at level :math:`\ell`, in any dimension, is
   ``amr.max_grid_size`` [:math:`\ell`] / 4.
