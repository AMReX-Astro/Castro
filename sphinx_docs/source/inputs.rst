***********
Input Files
***********

The Castro executable uses two inputs files at runtime to set and
alter the behavior of the algorithm and initial conditions.

The main inputs file, typically named ``inputs`` is used to set AMReX
parameters and the control flow in the C++ portions of Castro. Each
parameter here has a namespace (like ``amr.optionname`` or
``castro.optionname``).  Parameters set here are read using the AMReX
``ParmParse`` class infrastructure.

.. warning:: Because the inputs file is handled by the C++ portion of
   the code, any quantities you specify in scientific notation, must
   take the form ``1.e5`` and not ``1.d5``—the ``d`` specifier is not
   recognized.

Additionally, note that in Castro, all quantities are in CGS units.

The second inputs file, typically named ``probin`` is used by the
Fortran code that initializes the problem setup.  It is read at
problem initialization (via a Fortran namelist) and the
problem-specific quantities are stored in a Fortran module
``probdata_module`` defined in the problem’s ``probdata.f90`` file.

Only the ``inputs`` file is specified on the commandline. The
associated ``probin`` file is specified in the inputs file
using the ``amr.probin_file`` parameter, e.g.::

    amr.probin_file = my_special_probin

for example, has the Fortran code read a file called ``my_special_probin``.

Working with probin Files
=========================

There are three several Fortran namelists that can be defined in the
``probin`` file:

  * ``&fortin`` is the main namelist read by the problem’s
    ``probinit`` subroutine in the ``Prob_nd.F90`` file.

  * ``&extern`` is used to set different microphysics options

  * ``&tagging`` is used to get the parameters (defined in
    ``tagging_module``) that affect how we tag for refinement.

  * ``&sponge`` is used to set the parameters for the sponge source
    term that is used to damp velocities.


Common inputs Options
=====================


Problem Geometry
----------------

The ``geometry`` namespace is used by BoxLib to define the
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
    system. See § :ref:`soft:prob_params` for more details.

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
by BoxLib. The runtime parameters that we use are:

  * ``castro.lo_bc``: boundary type of each low face (must be set)

  * ``castro.hi_bc``: boundary type of each high face (must be set)

The valid boundary types are:

+-------------------------+------------------+
| 0 – Interior / Periodic | 3 – Symmetry     |
+-------------------------+------------------+
| 1 – Inflow              | 4 – Slip Wall    |
+-------------------------+------------------+
| 2 – Outflow             | 5 – No Slip Wall |
+-------------------------+------------------+

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

Regridding
----------

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

Simulation Time
---------------

There are two paramters that can define when a simulation ends:

  * ``max_step``: maximum number of level 0 time steps (integer
    :math:`\geq 0`; default: -1)

  * ``stop_time``: final simulation time (Real :math:`\geq 0`; default:
    -1.0)

To control the number of time steps, you can limit by the maximum
number of level 0 time steps (``max_step``) or by the final
simulation time (``stop_time``), or both. The code will stop at
whichever criterion comes first.

Note that if the code reaches ``stop_time`` then the final time
step will be shortened so as to end exactly at ``stop_time``, not
past it.

As an example::

    max_step  = 1000
    stop_time  = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level 0 steps taken equals 1000, whichever comes first.

Time Step
---------

If ``castro.do_hydro = 1``, then typically
the code chooses a time step based on the CFL number:

.. math::
   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\min\left\{\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right\}\right]
   :label: eq:cfl

If SDC integration is used instead, then we have

.. math::

   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\left(\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right)^{-1}\right]^{-1}

(If we are simulating in 1D or 2D, the extraneous parts related to :math:`v` and/or :math:`w` are removed.)

The following parameters affect the timestep choice:

  * ``castro.cfl``: CFL number (Real :math:`> 0` and :math:`\leq 1`;
    default: 0.8)

  * ``castro.init_shrink``: factor by which to shrink the initial
    time step (Real :math:`> 0` and :math:`\leq 1`; default: 1.0)

  * ``castro.change_max``: factor by which the time step can
    grow in subsequent steps (Real :math:`\geq 1`; default: 1.1)

  * ``castro.fixed_dt``: level 0 time step regardless of cfl
    or other settings (Real :math:`> 0`; unused if not set)

  * ``castro.initial_dt``: initial level 0 time
    step regardless of other settings (Real :math:`> 0`; unused if not set)

  * ``castro.dt_cutoff``: time step below which calculation
    will abort (Real :math:`> 0`; default: 0.0)

  * ``castro.hard_cfl_limit``: whether or not to abort the
    simulation if the hydrodynamics update creates velocities that
    violate the CFL criterion (Integer; default: 1)

As an example, consider::

    castro.cfl = 0.9
    castro.init_shrink = 0.01
    castro.change_max = 1.1
    castro.dt_cutoff = 1.e-20

This defines the :math:`\mathtt{cfl}` parameter in :eq:`eq:cfl` to be
0.9, but sets (via ``init_shrink``) the first timestep we take to
be 1% of what it would be otherwise. This allows us to ramp up to
the hydrodynamic timestep at the start of a simulation. The
``change_max`` parameter restricts the timestep from increasing by
more than 10% over a coarse timestep. Note that the time step can
shrink by any factor; this only controls the extent to which it can
grow. The ``dt_cutoff`` parameter will force the code to abort if
the timestep ever drops below :math:`10^{-20}`. This is a safety
feature—if the code hits such a small value, then something likely
went wrong in the simulation, and by aborting, you won’t burn through
your entire allocation before noticing that there is an issue.

If we know what we are doing, then we can force a particular timestep::

    castro.fixed_dt = 1.e-4

This sets the level 0 time step to be 1.e-4 for the entire simulation,
ignoring the other timestep controls. Note that if
``castro.init_shrink`` :math:`\neq 1` then the first time step will in fact
be ``castro.init_shrink`` :math:`\cdot` ``castro.fixed_dt``.

::

    castro.initial_dt = 1.e-4

sets the *initial* level 0 time step to be :math:`10^{-4}` regardless of
``castro.cfl`` or ``castro.fixed_dt``. The time step can
grow in subsequent steps by a factor of castro.change_max each step.

*diffusion*: If diffusion is enabled, the timestep will also be limited by:

.. math::

   \Delta t = \frac{1}{2}\min_{i,j,k}\left[\min\left\{\frac{\Delta x^2}{D_{i,j,k}},
                                                      \frac{\Delta y^2}{D_{i,j,k}},
                                                      \frac{\Delta z^2}{D_{i,j,k}}\right\}\right]

where :math:`D \equiv k / (\rho c_V)` if we are diffusing temperature,
and :math:`D \equiv k / (\rho c_P)` if we are diffusing enthalpy. No
input parameter is necessary to enable this constraint. See Chapter
:ref:`ch:diffusion` for more details.

*reactions*: If reactions are enabled, the timestep will also
be limited by two constraints:

.. math:: \Delta t = \mathtt{dtnuc\_e}\, \min_{i,j,k} \left\{\frac{e_{i,j,k}}{\dot{e}_{i,j,k}}\right\}

.. math:: \Delta t = \mathtt{dtnuc\_X}\, \min_{i,j,k} \left\{\min_n\frac{X^n_{i,j,k}}{\dot{X}^n_{i,j,k}}\right\}

where :math:`e` is the internal energy, and :math:`X^n` is the mass fraction of
the :math:`n`\ th species. The safety factors correspond to the runtime parameters
``castro.dtnuc_e`` and ``castro.dtnuc_X``. These limiters
say that the timestep must be small enough so that no zone can change
its internal energy by more than the fraction in one
step, and so that no zone can change the abundance of any isotope by
more than the fraction in one step. The time derivatives
:math:`\dot{e}` and :math:`\dot{X}^n` are estimated by calling the right-hand-side
of the nuclear network given the state at the time the timestep limiter
is being calculated. (We use a small number floor to prevent division by zero.)
To prevent the timestep from being dominated by trace species, there is
an additional option ``castro.dtnuc_X_threshold`` which is the
mass fraction threshold below which a species will not be considered in
the timestep constraint. and are set to
a large number by default, effectively disabling them. Typical choices
for these values in the literature are :math:`\sim 0.1`.

Subcycling
----------

Castro supports a number of different modes for subcycling in time,
set via amr.subcycling_mode.

  * ``amr.subcycling_mode`` = ``Auto`` (default): the code will run with
    equal refinement in space and time. In other words, if level
    :math:`n+1` is a factor of 2 refinement above level :math:`n`,
    then :math:`n+1` will take 2 steps of half the duration for every
    level :math:`n` step.

  * If ``amr.subcycling_mode`` = ``None``: the code will not refine in
    time. All levels will advance together with a timestep dictated by
    the level with the strictest :math:`dt`. Note that this is
    identical to the deprecated command ``amr.nosub = 1``.

  * If ``amr.subcycling_mode`` = ``Manual``: the code will subcycle
    according to the values supplied by ``amr.subcycling_iterations``.

In the case of ``amr.subcycling_mode`` = Manual, we subcycle in
manual mode with largest allowable timestep. The number of iterations
at each level is then specified as::

    amr.subcycling_iterations = 1 2 1 2

Here, we take 1 level-0 timestep at a time (required). Take 2 level-1
timesteps for each level-0 step, 1 timestep at level-2 for each
level-1 step, and take 2 timesteps at level-3 for each level-2 step.

Alternately, we could do::

    amr.subcycling_iterations = 2

which will subcycle twice at every level (except level 0).

Restart Capability
------------------

Castro has a standard sort of checkpointing and restarting capability.
In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

  * ``amr.check_file``: prefix for restart files (text;
    default: chk)

  * ``amr.check_int``: how often (by level 0 time steps) to
    write restart files (integer :math:`> 0`; default: -1)

  * ``amr.check_per``: how often (by simulation time) to
    write restart files (Real :math:`> 0`; default: -1.0)

    Note that ``amr.check_per`` will write a checkpoint at the first
    timestep whose ending time is past an integer multiple of this
    interval.  In particular, the timestep is not modified to match
    this interval, so you won’t get a checkpoint at exactly the time
    you requested.

  * ``amr.restart``: name of the file (directory) from which to
    restart (Text; not used if not set)

  * ``amr.checkpoint_files_output``: should we write
    checkpoint files? (0 or 1; default: 1)

    If you are doing a scaling study then set
    ``amr.checkpoint_files_output`` = 0 so you can test scaling of the
    algorithm without I/O.

  * ``amr.check_nfiles``: how parallel is the writing of
    the checkpoint files? (Integer :math:`\geq 1`; default: 64)

    See the chapter :ref:`ch:io` for more details on parallel I/O and the
    ``amr.check_nfiles`` parameter.

  * ``amr.checkpoint_on_restart``: should we write a
    checkpoint immediately after restarting? (0 or 1; default: 0)

  * ``castro.grown_factor``: factor by which domain has been
    grown (Integer :math:`\geq 1`; default: 1)

.. note:: You can specify both ``amr.check_int`` or ``amr.check_per``,
   if you so desire; the code will print a warning in case you did
   this unintentionally. It will work as you would expect – you will
   get checkpoints at integer multiples of ``amr.check_int`` timesteps
   and at integer multiples of ``amr.check_per`` simulation time
   intervals.

   ``amr.plotfile_on_restart`` and ``amr.checkpoint_on_restart``
   require amr.regrid_on_restart to be in effect.

As an example::

    amr.check_file = chk_run
    amr.check_int = 10

means that restart files (really directories) starting with the prefix
“chk_run” will be generated every 10 level-0 time steps. The
directory names will be ``chk_run00000``, ``chk_run00010``,
``chk_run00020``, etc.

If instead you specify::

    amr.check_file = chk_run
    amr.check_per = 0.5

then restart files (really directories) starting with the prefix
“chk_run” will be generated every 0.1 units of
simulation time. The directory names will be ``chk_run00000``,
``chk_run00043``, ``chk_run00061``, etc, where :math:`t = 0.1` after
43 level-0 steps, :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from ``chk_run00061``, for example, then set::

    amr.restart = chk_run00061

.. _sec:PlotFiles:

Controlling Plotfile Generation
-------------------------------

The main output from Castro is in the form of plotfiles (which are
really directories). The following options in the inputs file control
the generation of plotfiles:

  * ``amr.plot_file``: prefix for plotfiles (text; default:
    “plt”)

  * ``amr.plot_int``: how often (by level-0 time steps) to
    write plot files (Integer :math:`> 0`; default: -1)

  * ``amr.plot_per``: how often (by simulation time) to write
    plot files (Real :math:`> 0`; default: -1.0)

   .. note:: ``amr.plot_per`` will write a plotfile at the first
      timestep whose ending time is past an integer multiple of this
      interval.  In particular, the timestep is not modified to match
      this interval, so you won’t get a checkpoint at exactly the time
      you requested.

  * ``amr.plot_vars``: name of state variables to include in plotfiles
    (valid options: ALL, NONE or a list; default: ALL)

  * ``amr.derive_plot_vars``: name of derived variables to include in
    plotfiles (valid options: ALL, NONE or a list; default: NONE

  * ``amr.plot_files_output``: should we write plot files?
    (0 or 1; default: 1)

    If you are doing a scaling study then set
    ``amr.plot_files_output`` = 0 so you can test scaling of the
    algorithm without I/O.

  * ``amr.plotfile_on_restart``: should we write a plotfile
    immediately after restarting? (0 or 1; default: 0)

  * ``amr.plot_nfiles``: how parallel is the writing of the
    plotfiles? (Integer :math:`\geq 1`; default: 64)

    See the Software Section for more details on parallel I/O and the
    ``amr.plot_nfiles`` parameter.

  * ``castro.plot_X``: include all the species mass
    fractions in the plotfile (0 or 1; default: 0)

All the options for ``amr.derive_plot_vars`` are kept in
``derive_lst`` in ``Castro_setup.cpp``. Feel free to look at
it and see what’s there.

.. note:: You can specify both ``amr.plot_int`` or ``amr.plot_per``,
   if you so desire; the code will print a warning in case you did
   this unintentionally. It will work as you would expect – you will
   get plotfiles at integer multiples of amr.plot_int timesteps and at
   integer multiples of amr.plot_per simulation time intervals.

As an example::

    amr.plot_file = plt_run
    amr.plot_int = 10

means that plot files (really directories) starting with the prefix
“plt_run” will be generated every 10 level-0 time steps. The
directory names will be ``plt_run00000``, ``plt_run00010``,
``plt_run00020``, etc.

If instead you specify::

    amr.plot_file = plt_run
    amr.plot_per = 0.5

then restart files (really directories) starting with the prefix
“plt_run” will be generated every 0.1 units of simulation time. The
directory names will be ``plt_run00000``, ``plt_run00043``,
``plt_run00061``, etc, where :math:`t = 0.1` after 43 level-0 steps, :math:`t =
0.2` after 61 level-0 steps, etc.

Screen Output
-------------

There are several options that set how much output is written to the
screen as Castro runs:

  * ``amr.v``: verbosity of ``Amr.cpp`` (0 or 1; default: 0)

  * ``castro.v``: verbosity of ``Castro.cpp`` (0 or 1; default: 0)

  * ``gravity.v``: verbosity of ``Gravity.cpp`` (0 or 1; default: 0)

  * ``diffusion.v``: verbosity of ``Diffusion.cpp`` (0 or 1;
    default: 0)

  * ``mg.v``: verbosity of multigrid solver (for gravity) (allow
    values: 0, 1, 2, 3, 4; default: 0)

  * ``amr.grid_log``: name of the file to which the grids are
    written (text; not used if not set)

  * ``amr.run_log``: name of the file to which certain output is
    written (text; not used if not set)

  * ``amr.run_log_terse``: name of the file to which certain
    (terser) output is written (text; not used if not set)

  * ``amr.sum_interval``: if :math:`> 0`, how often (in level-0 time
    steps) to compute and print integral quantities (Integer; default: -1)

    The integral quantities include total mass, momentum and energy in
    the domain every ``castro.sum_interval`` level-0 steps.  The print
    statements have the form::

           TIME= 1.91717746 MASS= 1.792410279e+34

   for example. If this line is commented out then
   it will not compute and print these quanitities.

  * ``castro.do_special_tagging``: allows the user to set a special
    flag based on user-specified criteria (0 or 1; default: 1)

    ``castro.do_special_tagging`` = 1 can be used, for example, to
    calculate the bounce time in a core collapse simulation; the
    bounce time is defined as the first time at which the maximum
    density in the domain exceeds a user-specified value. This time
    can then be printed into a special file as a useful diagnostic.

As an example::

    amr.grid_log = grdlog
    amr.run_log = runlog

Every time the code regrids it prints a list of grids at all relevant
levels. Here the code will write these grids lists into the file
``grdlog``. Additionally, every time step the code prints certain
statements to the screen (if ``amr.v`` = 1), such as::

    STEP = 1 TIME = 1.91717746 DT = 1.91717746
    PLOTFILE: file = plt00001

The ``run_log`` option will output these statements into
*runlog* as well.

Terser output can be obtained via::

    amr.run_log_terse = runlogterse

This file, ``runlogterse`` differs from ``runlog``, in that it
only contains lines of the form::

    10  0.2  0.005

in which “10” is the number of steps taken, “0.2” is the
simulation time, and “0.005” is the level-0 time step. This file
can be plotted very easily to monitor the time step.

Other parameters
----------------

There are a large number of solver-specific runtime parameters. We describe these
together with the discussion of the physics solvers in later chapters.
