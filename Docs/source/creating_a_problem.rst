***************************
Setting Up Your Own Problem
***************************

Castro problems are organized loosely into groups describing their
intent (e.g., science, hydro tests, ...).  These groups are
sub-directories under the ``Castro/Exec/`` directory.  Each problem is
then placed in a sub-directory of the appropriate group (for example,
``Castro/Exec/hydro_tests/Sedov`` holds the Sedov test problem).

To create a new problem, you will create a new directory under one
of the groups and place in it the following files:

  * ``GNUmakefile`` : the makefile for this problem.  This will tell
    Castro what options to use and what network and EOS to build.

  * ``problem_initialize.H`` and
    ``problem_initialize_state_data.H`` : this holds the problem
    initialization routines.  MHD and radiation problems require
    an additional file.

  * ``_prob_params`` (optional) : a list of runtime parameters that
    you problem will read.  These parameters are controlled by the
    set in the inputs file as ``problem.param`` where ``param`` is
    one of the parameters listed in ``_prob_params``.

  * ``Make.package`` : this is a makefile fragment that is included
    during the build process.  It tells the build system about any
    problem-specific files that need to be compiled.

  * ``inputs`` : this is the main inputs file that controls Castro and
    AMReX's behavior.

The best way to get started writing a new problem is to copy an
existing problem setup into a new directory.

.. index:: _prob_params

Runtime Parameters
------------------

The problem-specific runtime parameters are defined in ``_prob_params``.
This has the form::

   name       datatype    default      namelist?     size

Here:

* `name` is the name of the variable to put into ``probdata_module``

* `datatype` is one of ``real``, ``integer``, ``character``, or
  ``logical``.

* `default` is the default value of the runtime parameter.  It may be
  overridden at runtime by reading from the namelist.

* `namelist` indicates if the variable should be able to be set as
  a runtime parameter.  If it is empty or marked with
  ``n``, then the variable will still be creates in the ``problem`` namespace,
  but it will not be able to be set via the commandline or inputs file.
  A common usage of this is to define global variables that might be set
  at problem initialization that are used elsewhere in Castro.

* `size` is for arrays, and gives their size.  It can be any integer
  or variable that is known to the ``probdata_module``.  If you need a
  variable from a different module to specify the size (for example,
  ``nspec`` from ``network``), then you express the size as a tuple:
  ``(nspec, network)`` and the appropriate use statement will be
  added.

The variables will all be initialized for the GPU as well.


Problem Initialization
----------------------

Here we describe the main problem initialization routines.

.. index:: initialize_problem

* ``initialize_problem()``

  This C++ routine is primarily responsible for doing any one-time
  initialization for the problem (like reading in an
  initial model through the ``model_parser.H`` functionality—see the
  ``toy_convect`` problem setup for an example.

  This routine can also postprocess any of the parameters defined
  in the ``_prob_params`` file, which are defined in ``problem`` namespace.

  .. note:: many problems set the value of the ``problem::center[]`` array
     from the ``problem`` namespace.  This is used to note the
     center of the problem (which does not necessarily need to be
     the center of the domain, e.g., for axisymmetric problems).
     ``center`` is used in source terms (including rotation and
     gravity) and in computing some of the derived variables (like
     angular momentum).

  If you need coordinate information, it can be obtained
  by constructing a ``Geometry`` object using ``DefaultGeometry()``
  and accessing its ``ProbLo()`` and ``ProbHi()`` methods.


* ``problem_initialize_state_data()``:

  This routine will initialize the state data in a given zone.
  The arguments passed are:

  - ``i``, ``j``, ``k``: the index of the zone to fill the data in

  - ``state``: an array containing the simulation state data

  - ``GeomData``: a ``GeometryData`` object that can be used for obtaining
    ``dx``, ``problo``, ``probhi``, etc.

  Filling data is done by simply writing to ``state(i,j,k,URHO)``, etc.

.. _create:bcs:

Boundary Conditions
-------------------

.. index:: boundary conditions

Standard boundary conditions, including outflow (zero-gradient), periodic,
and symmetry (reflect) are handled by AMReX directly.  Castro has a special
hydrostatic boundary condition that can be used for the lower boundary.  It
is accessed by setting the ``castro.lo_bc`` flag to 1 in the vertical coordinate
direction, e.g., for 2-d as::

   castro.lo_bc       =  0   1

The flag value 1 is traditionally named "inflow" by AMReX, but generally means that
the boundary implementation is left to the user.  To tell Castro to use the
hydrostatic boundary condition here, we set::

   castro.yl_ext_bc_type = 1
   castro.hse_interp_temp = 1
   castro.hse_reflect_vels = 1

The first parameter tells Castro to use the HSE boundary condition for the lower
y direction.
In filling the ghost cells, hydrostatic equilibrium will be integrated
from the last interior zone into the boundary.  We need one more
equation for this integration, so we either interpolate the density or
temperature into the ghost cells, depending on the value of
``castro.hse_interp_temp``.  Finally, ``castro.hse_reflect_vels``
determines how we treat the velocity.  The default is to give is a
zero gradient, but in tests we've found that reflecting the velocity
while integrating the HSE profile can be better.  For modeling a
plane-parallel hydrostatic atmosphere, using the hydrostatic boundary
conditions instead of a simple symmetry boundary is essential when
using the standard CTU PPM solver.

.. index:: castro.fill_ambient_bc, castro.ambient_fill_dir, castro.ambient_outflow_vel

A different special boundary condition, based on outflow, is available at
the upper boundary.  This works together with the ``model_parser``
module to fill the ghost cells at the upper boundary with the initial
model data.  You set this as::

   castro.hi_bc = 2 2

   castro.fill_ambient_bc = 1
   castro.ambient_fill_dir = 1
   castro.ambient_outflow_vel = 1

where ``ambient_fill_dir`` is the 0-based direction to fill using an
ambient state defined by the problem setup.  In this example, we will
override the outflow (2) boundary condition in the y-direction.  That
problem setup needs to fill the ``ambient_state[:]`` array defined in
``ambient.H``.  An example of using this boundary is in the
``flame_wave`` problem.

The implementations of these boundary conditions is found in
``Castro/Source/problems/Castro_bc_fill_nd.cpp``.

Optional Files
--------------

The follow problem-specific files are optional. There are stubs for
each of these in the main source tree.

-  ``problem_checkpoint.H``, ``problem_restart.H`` :

   These provides two routines, respectively ``problem_checkpoint`` and
   ``problem_restart`` that can be used to add information to the
   checkpoint files and read it in upon restart. This is useful for
   some global problem-specific quantities. For instance, the
   ``wdmerger`` problem uses this to store center of mass position and
   velocity information in the checkpoint files that are used for
   runtime diagnostics.

   The name of the checkpoint directory is passed in as an argument.

-  ``problem_tagging.H``

   This implements problem-specific tagging for refinement, through a
   the function ``problem_tagging``. The full hydrodynamic state (State_Type)
   is passed in, and the problem can mark zones for refinement by setting the
   tag variable for a zone to set. An example is provided by the ``toy_convect``
   problem which refines a rectangular region (fuel layer) based on
   a density parameter and the H mass fraction.

   .. _problem_derives:

-  ``Problem_Derives.H``, ``Problem_Derive.H``, and ``Problem_Derive.cpp``

   Together, these provide a mechanism to create derived quantities
   that can be stored in the plotfile. ``Problem_Derives.H``
   provides the C++ code that defines these new plot variables. It
   does this by adding them to the ``derive_lst``—a list of
   derived variables that Castro knows about. When adding new
   variables, a descriptive name, a C++ routine that does the
   deriving, and component of ``StateData`` are specified.

   The other two files provide the header and implementation of the
   function that computes the derived variable.  A example is provided
   by the ``reacting_bubble`` problem, which derives several new
   quantities (perturbations against a background one-dimensional
   model, in this case).

-  ``Prob.cpp``, ``Problem.H``

   These files provide problem-specific routines for computing global
   diagnostic information through the sum_integrated_quantities
   functionality that is part of the ``Castro`` class.

   An example is provided by ``toy_flame``, where an estimate
   of the flame speed is computed by integrating the mass of fuel on
   the grid.


Model Parser
------------

.. index:: USE_MODEL_PARSER, MAX_NPTS_MODEL

Many problem setups begin with a 1-d initial model that is mapped onto
the grid.  The ``model_parser.H`` provides the functions that read in
the initial model and map it on the Castro grid.  To enable this, add::

  USE_MODEL_PARSER = TRUE

to the problem ``GNUmakefile``.  There are 2 other parameters that can
be set in the makefile to control the initial model storage:

  * ``MAX_NPTS_MODEL``: is the maximum number of data points in the
    1-d initial model.  This needs to be known at compile time so we
    can make the data managed for GPUs.

  * ``NUM_MODELS``: this is the number of different initial models we
    want to managed.  Typically we only want 1, but some problems,
    like ``flame_wave`` use 2, applied to different portions of the
    domain.

Two different model formats are allowed:

* *"Old" model format*:  This has two header lines at the top giving
  the number of points in the model and number of variables, and
  then the name of each variable (excluding the coordinate) is
  given separately on the next lines:

  ::

      # npts = 896
      # num of variables = 6
      # density
      # temperature
      # pressure
      # carbon-12
      # oxygen-16
      # magnesium-24
      195312.5000  5437711139.  8805500.952  0.4695704813E+28  0.3  0.7  0
      585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0

  The data begins with the coordinate and then the variables in the
  model, with one data point per line.

* *New model format*: This has a single header line that gives the column
  names, including the coordinate (first column), followed by the
  data.  For example:

  ::

      # radius  density temperature  pressure  carbon-12  oxygen-16  magnesium-24
      195312.5000  5437711139.  8805500.952  0.4695704813E+28  0.3  0.7  0
      585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0

.. note::

   For both formats, the variable names should be the same as the
   names used by Castro.  You can see the names in the `_variables
   <https://github.com/AMReX-Astro/Castro/blob/main/Source/driver/_variables>`_
   file.

When the model is read, the variables listed in the file are matched
to the ones that Castro knows about.  If the variable is recognized,
then it is stored in the model data, otherwise, it is ignored.

The data can then be mapped onto the grid using the ``interpolate()``
function, e.g., ::

    Real dens = interpolate(height, model::idens);

This fills ``dens`` with the density at the position ``height``.  In
addition to density, you can specify temperature (``model::itemp``),
pressure (``model::ipres``), species (indexed from ``model::ispec``),
or an auxiliary quantity (indexed from ``model::iaux``).
