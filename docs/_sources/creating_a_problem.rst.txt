***************************
Setting Up Your Own Problem
***************************

Castro problems are organized loosely into groups describing their
intent (e.g., science, hydro tests, ...).  These groups are
sub-directories under the `Castro/Exec/` directory.  Each problem is
then placed in a sub-directory of the appropriate group (for example, 
``Castro/Exec/hydro_tests/Sedov`` holds the Sedov test problem).

To create a new problem, you will create a new directory under one
of the groups and place in it the following files:

  * ``GNUmakefile`` : the makefile for this problem.  This will tell
    Castro what options to use and what network and EOS to build.

  * ``Prob_nd.F90`` OR ``problem_initialize.H`` and
    ``problem_initialize_state_data.H`` : this holds the problem
    initialization routines, which may be implemented either in Fortran
    or in C++.

  * ``_prob_params`` (optional) : a list of runtime parameters that
    you problem will read.  These parameters are controlled by the
    ``probin`` file.

  * ``Make.package`` : this is a makefile fragment that is included
    during the build process.  It tells the build system about any
    problem-specific files that need to be compiled.

  * ``inputs`` : this is the main inputs file that controls Castro and
    AMReX's behavior.

  * ``probin`` : this is the problem-specific inputs file that
    contains Fortran namelists with problem parameters as well as any
    paramters for external physics (e.g., from the StarKiller
    Microphysics)

The best way to get started writing a new problem is to copy an
existing problem setup into a new directory.

.. index:: _prob_params

Runtime parameters
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

* `namelist` indicates if the variable should be in the namelist and
  controlled at runtime.  The namelist column should have a ``y`` if
  you want the parameter to be read from the ``&fortin`` namelist in
  the ``probin`` file at runtime.  If it is empty or marked with
  ``n``, then the variable will still be put into the
  ``probdata_module`` but it will not be initialized via the namelist.

* `size` is for arrays, and gives their size.  It can be any integer
  or variable that is known to the ``probdata_module``.  If you need a
  variable from a different module to specify the size (for example,
  ``nspec`` from ``network``), then you express the size as a tuple:
  ``(nspec, network)`` and the appropriate use statement will be
  added.

The variables will all be initialized for the GPU as well.


``Problem Initialization``
---------------

Here we describe the main problem initialization routines. There are
two implementations, in C++ (``problem_setup.H``) and Fortran (``Prob_nd.F90``),
and you can pick either but not both (C++ is recommended since eventually
we will switch the whole code to C++).

.. index:: probdata

* ``amrex_probinit()`` (Fortran) or ``initialize_problem()`` (C++):

  This routine is primarily responsible for doing any one-time
  initialization for the problem (like reading in an
  initial model through the ``model_parser_module``—see the
  ``toy_convect`` problem setup for an example).

  This routine can also postprocess any of the parameters defined
  in the ``_prob_params`` file, which are defined in ``probdata_module``.

  .. note:: many problems set the value of the ``center()`` array
     from ``prob_params_module`` here (in C++, the ``center[]`` variable
     from the ``problem`` namespace).  This is used to note the
     center of the problem (which does not necessarily need to be
     the center of the domain, e.g., for axisymmetric problems).
     ``center`` is used in source terms (including rotation and
     gravity) and in computing some of the derived variables (like
     angular momentum).

  In Fortran the arguments include the name and length of the probin file
  as well as the physical values of the domain's lower-left corner
  (``problo``) and upper-right corner (``probhi``). The C++ version
  accepts no arguments, but for example ``problo`` and ``probhi`` can
  be obtained by constructing a ``Geometry`` object using ``DefaultGeometry()``
  and accessing its ``ProbLo()`` and ``ProbHi()`` methods.


* ``ca_initdata()`` (Fortran) or ``initialize_problem_state_data()`` (C++):

  This routine will initialize the state data for a single grid.
  In Fortran the inputs to this routine are:

  -  ``level``: the level of refinement of the grid we are filling

  -  ``time``: the simulation time

  -  ``lo()``, ``hi()``: the integer indices of the box’s
     *valid data region* lower left and upper right corners. These
     integers refer to a global index space for the level and
     identify where in the computational domain the box lives.

  -  ``nscal``: the number of scalar quantities—this is not typically
     used in Castro.

  -  ``state()``: the main state array. This is dimensioned as::

       real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

     where ``NVAR`` comes from the ``meth_params_module``.  The
     spatial dimensions of the array come in through the arguments
     ``s_lo`` and ``s_hi``.

     When accessing this array, we use the index keys provided by
     meth_params_module (e.g., ``URHO``) to refer to specific
     quantities

  -  ``delta()``: this is an array containing the zone width (:math:`\Delta x`)
     in each coordinate direction: ``delta(1)`` = :math:`\Delta x`,
     ``delta(2)`` = :math:`\Delta y`, ...

  -  ``xlo()``, ``xhi()``: these are the physical coordinates of the
     lower left and upper right corners of the *valid region*
     of the box.  These should not be used, and will be removed in a future
     version of Castro.

Filling data is typically done in a loop like::

     do k = lo(3), hi(3)
        z = (dble(k)+HALF)*delta(3) + problo(3)

        do j=lo(2),hi(2)
           y = (dble(j)+HALF)*delta(2) + problo(2)

           do i=lo(1),hi(1)
              x = (dble(i)+HALF)*delta(1) + problo(1)

              state(i,j,k,URHO) = ...

           end do
        end do
     end do

Here, we compute the coordinates of the zone center, ``x``, ``y``, and ``z``
from the zone indices, ``i``, ``j``, and ``k``.

  In C++, the arguments passed are:

  - ``i``, ``j``, ``k``: the index of the zone to fill the data in

  - ``state``: an array containing the simulation state data

  - ``GeomData``: a ``GeometryData`` object that can be used for obtaining
    ``dx``, ``problo``, ``probhi``, etc.

  Filling data is done by simply writing to ``state(i,j,k,URHO)``, etc.

.. _create:bcs:

Boundary conditions
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

   castro.yl_ext_bc_type = "hse"
   castro.hse_interp_temp = 1
   castro.hse_reflect_vels = 1

The first parameter tells Castro to use the HSE boundary condition.
In filling the ghost cells, hydrostatic equilibrum will be integrated
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
problem setup needs to fill the ``ambient_state(:)`` array defined in
``ambient_module``.  An example of using this boundary is in the
``flame_wave`` problem.

The implementations of these boundary conditions is found in
``Castro/Source/problems/bc_ext_fill_nd.F90``.

If a problem requires different initial conditions, then they should
put a version of ``bc_ext_fill_nd.F90`` into the problem directory and
modify it as needed.  See the ``double_mach_reflection`` problem for
an example of this.

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

-  ``problem_tagging_nd.F90`` OR ``problem_tagging.H``

   This implements problem-specific tagging for refinement, through a
   subroutine ``set_problem_tags`` (Fortran) or function ``problem_tagging``
   (C++). The full hydrodynamic state (State_Type) is passed in, and the
   problem can mark zones for refinement by setting the tag variable for
   a zone to set. An example is provided by the ``toy_convect``
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

-  ``Prob.cpp``, ``Problem.H``

   These files provide problem-specific routines for computing global
   diagnostic information through the sum_integrated_quantities
   functionality that is part of the ``Castro`` class.

   An example is provided by ``toy_flame``, where an estimate
   of the flame speed is computed by integrating the mass of fuel on
   the grid.

