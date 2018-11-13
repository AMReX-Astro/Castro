***************************
Setting Up Your Own Problem
***************************

To define a new problem, we create a new directory in one
of the subdirectories of ``Exec/``,
and place in it a ``Prob_2d.f90`` file (or 1d/3d,
depending on the dimensionality of the problem), a ``probdata.f90``
file, the ``inputs`` and ``probin`` files, and a
``Make.package`` file that tells the build system what problem-specific
routines exist. Finally, if you need custom boundary conditions, a
``bc_fill_2d.F90`` (or 1d/3d) file is needed. The
simplest way to get started is to copy these files from an existing
problem. Here we describe how to customize your problem.

The purpose of these files is:

-  ``probdata.f90``: this holds the ``probdata_module`` Fortran module
   that allocates storage for all the problem-specific runtime parameters that
   are used by the problem (including those that are read from the ``probin``
   file.

-  ``Prob_?d.f90``: this holds the main routines to
   initialize the problem and grid and perform problem-specific boundary
   conditions:

   -  ``probinit()``:

      This routine is primarily responsible for reading in the
      ``probin`` file (by defining the ``&fortin`` namelist and
      reading in an initial model (usually through the
      ``model_parser_module``—see the ``toy_convect`` problem
      setup for an example). The parameters that are initialized
      here are those stored in the ``probdata_module``.

   -  ``ca_initdata()``:

      This routine will initialize the state data for a single grid.
      The inputs to this routine are:

      -  ``level``: the level of refinement of the grid we are filling

      -  ``time``: the simulation time

      -  ``lo()``, ``hi()``: the integer indices of the box’s
         *valid data region* lower left and upper right corners. These
         integers refer to a global index space for the level and
         identify where in the computational domain the box lives.

      -  ``nscal``: the number of scalar quantities—this is not typically
         used in Castro.

      -  ``state_l1``, ``state_l2``, (``state_l3``): the
         integer indices of the lower left corner of the box in each
         coordinate direction. These are for the box as allocated in memory,
         so they include any ghost cells as well as the valid data regions.

      -  ``state_h1``, ``state_h2``, (``state_h3``): the
         integer indices of the upper right corner of the box in each
         coordinate direction. These are for the box as allocated in memory,
         so they include any ghost cells as well as the valid data regions.

      -  ``state()``: the main state array. This is dimensioned as::

             double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

         (in 2-d), where ``NVAR`` comes from the ``meth_params_module``.

         When accessing this array, we use the index keys provided by
         meth_params_module (e.g., ``URHO``) to refer to specific
         quantities

      -  ``delta()``: this is an array containing the zone width (:math:`\Delta x`)
         in each coordinate direction: :math:`\mathtt{delta(1)} = \Delta x`,
         :math:`\mathtt{delta(2)} = \Delta y`, :math:`\ldots`.

      -  ``xlo()``, ``xhi()``: these are the physical coordinates of the
         lower left and upper right corners of the *valid region*
         of the box. These can be used to compute the coordinates of the
         cell-centers of a zone as::

               do j = lo(2), hi(2)
                  y = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)
                  ...

         .. note:: this method works fine for the problem
            initialization stuff, but for routines that implement
            tiling, as discussed below, ``lo`` and ``xlo`` may not
            refer to the same corner, and instead coordinates should
            be computed using ``problo()`` from the
            ``prob_params_module``.

-  ``bc_fill_?d.F90``:

   These routines handle how Castro fills ghostcells
   *at physical boundaries* for specific data. Most problem
   setups won’t need to do anything special here, and inclusion
   of this file is optional–only use it if you need to set
   specific boundary conditions.

   These routines are registered in ``Castro_setup.cpp``, and
   called as needed. By default, they just
   pass the arguments through to ``filcc``, which handles all of
   the generic boundary conditions (like reflecting, extrapolation,
   etc.). The specific ‘fill’ routines can then supply the
   problem-specific boundary conditions, which are typically just
   Dirichlet boundary conditions (usually this means looking to see
   if the ``bc()`` flag at a boundary is ``EXT_DIR``). The
   problem-specific code implementing these specific conditions
   should *follow* the ``filcc`` call.

   -  ``ca_hypfill``:
      This handles the boundary filling for the hyperbolic system.

   -  ``ca_denfill``: At times, we need to fill just the density
      (always assumed to be the first element in the hyperbolic state)
      instead of the entire state. When the fill patch routine is called
      with ``first_comp`` = ``Density`` and ``num_comp`` = 1, then we
      use ``ca_denfill`` instead of ``ca_hypfill``.

      (Note: it seems that this may be used for more than just
      density, but it is only used for tagging and the plotfile)

   -  ``ca_grav?fill``: These routines fill will the ghostcells
      of the gravitational acceleration grids with the gravitational
      acceleration.

      Note: for constant gravity, these routines will never be called.
      For one of the Poisson-type gravities, you only need to do
      something special here if you are implementing an Interior
      boundary type (which you can test for by comparing
      ``bc(:,:,:)`` to ``EXT_DIR``.

      For the other standard physical boundary types, the ghost cell
      filling will be handled automatically by the default ``filcc``
      call in these routines.

      The gravitational acceleration in the ghost cells is used during
      the hydrodynamics portion of the code in predicting the
      interface states.

   -  ``ca_reactfill``: This handles boundary filling for
      any ``Reactions_Type`` ``MultiFAB``s, which are sometimes used to interface
      with the nuclear burning module. It stores the normal state data
      in addition to components for the energy release and species change.

   These routines take the following arguments:

   -  ``adv_l1``, ``adv_l2``, (``adv_l3``): the indicies of
      the lower left corner of the box holding the data we are working on.
      These indices refer to the entire box, including ghost cells.

   -  ``adv_h1``, ``adv_h2``, (``adv_h3``): the indicies of
      the upper right corner of the box holding the data we are working on.
      These indices refer to the entire box, including ghost cells.

   -  ``adv()``: the array of data whose ghost cells we are filling.
      Depending on the routine, this may have an additional index refering
      to the variable.

      This is dimensioned as::

            double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

   -  ``domlo()``, ``domhi()``: the integer indices of the lower
      left and upper right corners of the valid region of the *entire
      domain*. These are used to test against to see if we are filling
      physical boundary ghost cells.

      This changes according to refinement level: level-0 will
      range from 0 to ``castro.max_grid_size``,
      and level-n will range from 0 to
      :math:`\mathtt{castro.max\_grid\_size} \cdot \prod_n \mathtt{castro.ref\_ratio(n)}`.

   -  ``delta()``: is the zone width in each coordinate direction,
      as in ``initdata()`` above.

   -  ``xlo()``: this is the physical coordinate of the lower
      left corner of the box we are filling—including the ghost cells.

      .. note:: this is different than how ``xlo()`` was defined in
         ``initdata()`` above.

   -  ``time``: the simulation time

   -  ``bc()``: an array that holds the type of boundary conditions
      to enforce at the physical boundaries for ``adv``.

      Sometimes it appears of the form ``bc(:,:)`` and sometimes
      ``bc(:,:,:)``—the last index of the latter holds the variable
      index, i.e., density, pressure, species, etc.

      The first index is the coordinate direction and the second index
      is the domain face (1 is low, 2 is hi), so
      ``bc(1,1)`` is the lower :math:`x` boundary type, ``bc(1,2)`` is
      the upper :math:`x` boundary type, ``bc(2,1)`` is the lower
      :math:`y` boundary type, etc.

      To interpret the array values, we test against the quantities
      defined in ``bc_types.fi`` included in each subroutine,
      for example, ``EXT_DIR``, ``FOEXTRAP``, :math:`\ldots`. The
      meaning of these are explained below.

Optional Files
--------------

The follow problem-specific files are optional. There are stubs for
each of these in the main source tree.

-  ``Problem.f90`` :

   This provides two routines, ``problem_checkpoint`` and
   ``problem_restart`` that can be used to add information to the
   checkpoint files and read it in upon restart. This is useful for
   some global problem-specific quantities. For instance, the
   wdmerger [5]_ problem uses this
   to store center of mass position and velocity information in the
   checkpoint files that are used for runtime diagnostics.

   The name of the checkpoint directory is passed in as an argument.
   ``Problem_F.H`` provides the C++ interfaces for these routines.

-  ``problem_tagging_?d.F90``, ``problem_tagging_nd.F90``

   This implements problem-specific tagging for refinement, through a
   subroutine ``set_problem_tags``. The full hydrodynamic state
   (State_Type) is passed in, and the problem can mark zones for
   refinement by setting the tag variable for a zone to
   set. An example is provided by the ``toy_convect``
   problem which refines a rectangular region (fuel layer) based on
   a density parameter and the H mass fraction.

-  ``Problem_Derive_F.H``, ``Problem_Derives.H``, ``problem_derive_nd.f90``

   Together, these provide a mechanism to create derived quantities
   that can be stored in the plotfile. ``Problem_Derives.H``
   provides the C++ code that defines these new plot variables. It
   does this by adding them to the ``derive_lst``—a list of
   derived variables that Castro knows about. When adding new
   variables, a descriptive name, Fortran routine that does the
   deriving, and component of ``StateData`` are specified.

   The Fortran routine that does the deriving is put in the
   problem-specific ``problem_derive_nd.f90`` (and a prototype for
   C++ is put in ``Problem_Derives.H``). A example is provided by
   the ``reacting_bubble`` problem, which derives several new
   quantities (perturbations against a background one-dimensional
   model, in this case).

-  ``Prob.cpp``, ``Problem.H``, ``Problem_F.H``

   These files provide problem-specific routines for computing global
   diagnostic information through the sum_integrated_quantities
   functionality that is part of the ``Castro`` class.

   An example is provided by ``toy_flame``, where an estimate
   of the flame speed is computed by integrating the mass of fuel on
   the grid.

Dimension Agnostic Problem Initialization
-----------------------------------------

Most of the problem setups have separate implementations for 1-, 2-,
and 3D. A new method exists that allows you to write just a single
set of files for any dimensionality (this is called the *dimension
agnostic* format). To use this mode, set
``DIMENSION_AGNOSTIC`` = ``TRUE`` in your ``GNUmakefile``.
Then write you problem initialization in ``Prob_nd.F90``.
Analogous routines exist for tagging and boundary conditions. See the
``rotating_torus`` and ``Noh`` problem setups for an
example.

.. _software:io:


.. [5]
   available separately at
   https://github.com/BoxLib-Codes/wdmerger
