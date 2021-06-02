.. _ch:buildsystem:

*********************
Build System Overview
*********************


Make Parameters
---------------

These build parameters control the parallelism, included physics,
etc.  In general they can be set to ``TRUE`` or ``FALSE``.  The
Castro-specific ones are interpreted by ``Make.Castro``.

A lot of optional features are enabled at compile time.  This allows
Castro to reduce the memory footprint of the state arrays by not allocating
space for variables that are not used.

General Build Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: USE_ALL_CASTRO, USE_AMR_CORE, USE_HYPRE

These Parameters affect the build (parallelism, performance, etc.)
Most of these are parameters from AMReX.

  * ``USE_ALL_CASTRO``: compile all of the core Castro directories.
    This is the defailt (``TRUE``), and should not be changed for
    general simulations.  The purpose of this flag is for unit tests, which
    do not need all of the Castro directories compiled.  

  * ``USE_AMR_CORE``: compile all of the core AMReX directories, including
    ``Base/``, ``AmrCore/``, ``Amr/``, and ``Boundary/``.  This defaults
    to ``TRUE`` and should be left set for Castro simulations.  The purpose
    of this flag is for unit tests that don't need all of AMReX.

  * ``USE_MLMG``: use the AMReX multi-level multigrid solver for gravity
    and diffusion.  This should always be set to ``TRUE``.

  * ``USE_HYPRE``: compile in the Hypre library.  This will be automatically enabled
    for radiation.  You need to specify the path to the Hypre library via either
    ``HYPRE_DIR`` or ``HYPRE_OMP_DIR``.

Parallelization and GPUs
^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: USE_MPI, USE_OMP, USE_CUDA, USE_ACC

The following parameters control how work is divided across nodes, cores, and GPUs.

  * ``USE_CUDA``: compile with GPU support using CUDA. 

  * ``USE_ACC``: compile with OpenACC. Note: this is a work in
    progress and should not be used presently.


  * ``USE_MPI``: compile with the MPI library to allow for distributed parallelism.

  * ``USE_OMP``: compile with OpenMP to allow for shared memory parallelism.



General Physics Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: USE_SIMPLIFIED_SDC, USE_TRUE_SDC

The following parameters control how the coupling between hydro and reactions
is handled.

  * ``USE_SIMPLIFIED_SDC``: use the simplified spectral deferred corrections (SDC)
    solver for coupling hydro and reactions.  At the moment, this
    works with the CTU hydrodynamics solver.  This requires running with
    ``castro.time_integration_method = 3``.

  * ``USE_TRUE_SDC``: use the true SDC method to couple hydro and
    reactions.  This can do 2nd order or 4th order accuracy.  At the
    moment, this works on single level only.  This requires running
    with ``castro.time_integration_method = 2``.



Radiation Parameters
^^^^^^^^^^^^^^^^^^^^

  * ``USE_RAD``: use photon radiation diffusion.  Note, For
    multigroup radiation, you need to set the number of radiation
    groups.  This is controlled by the ``NGROUPS`` parameter.

    .. index:: USE_RAD, NGROUPS



Gravity Parameters
^^^^^^^^^^^^^^^^^^

  * ``USE_GRAV``: use gravity (this could be constant or self-gravity)

    .. index:: USE_GRAV

  * ``USE_SELF_GRAV``: use self-gravity.  At the moment, this is always set
    if ``USE_GRAV`` is enabled.

  * ``USE_GR``: use a post-Newtonian approximation for GR gravity for the monopole
    solver.

    .. index:: USE_GR

  * ``USE_POINTMASS``: include a pointmass source to the gravitational potential.

    .. index:: USE_POINTMASS

Microphysics Parameters
^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_DIFFUSION``: enable thermal diffusion.  The conductivity is
    set via ``CONDUCTIVITY_DIR``, which should be a directory in the
    Microphysics repo.

    .. index:: USE_DIFFUSION, CONDUCTIVITY_DIR

  * ``USE_REACT``: enable reactions.  When reactions are set, we need
    to specify a network and an integrator.  Typically these come from
    the Microphysics repo, but one common exception is the
    ``general_null`` network, which just defines a composition.  The
    parameters that come into play here are:

    * ``NETWORK_DIR``: the network to use.  This is expected to be a subdirectory
      in the Microphysics repo.

    * ``GENERAL_NET_INPUTS``: this is the text file that we read to define the
      composition if we are using the ``general_null`` network.

    * ``INTEGRATOR_DIR``: this is the ODE integrator to use to integrate the 
      reaction system.  This is expected to be a subdirectory in the Microphysics
      repo.

    .. index:: USE_REACT, general_null, GENERAL_NET_INPUTS, NETWORK_DIR, INTEGRATOR_DIR

  * ``USE_REACT_SPARSE_JACOBIAN``

  * ``USE_SPARSE_STOP_ON_OOB``

  * ``EOS_DIR``: the equation of state to use.  This will be a subdirectory under the
    Microphysics repo.

    .. index:: EOS_DIR


Hydrodynamics and Source Term Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_ROTATION``: include rotation sources

    .. index:: USE_ROTATION

  * ``USE_HYBRID_MOMENTUM``: have Castro evolve angular momentum in addition to linear
    momentum.

    .. index:: USE_HYBRID_MOMENTUM

  * ``USE_SHOCK_VAR``: include a variable in the State_Type StateData that marks the
    location of a shock.

    .. index:: USE_SHOCK_VAR


Simulation Flow Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_AUX_UPDATE``: some networks define auxillary quantities, which in general
    Castro will advect, but not otherwise change.  If we set ``USE_AUX_UPDATE=TRUE``
    then Castro will call a user-supplied routine ``advance_aux()`` that can
    change the auxillary quantities.

    .. index:: USE_AUX_UPDATE

  * ``USE_POST_SIM``: if this is defined, then Castro will call the user-defined 
    routine ``problem_post_simulation()`` after the full evolution of the problem
    has ended.

    .. index:: USE_POST_SIM

  * ``USE_MAESTRO_INIT``: this enables the code to allow Castro to restart from a 
    Maestro simulation.  This will need to be updated in the future to allow for 
    restarts from MAESTROeX.

    .. index:: USE_MAESTRO_INIT

  * ``USE_HDF5``: compile in support for HDF5.  This is needed for some tables used
    by Microphysics routines.

    .. index:: USE_HDF5

Tracer Particle Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_PARTICLES``: compile in support for tracer particles.





Build Process Procedure
-----------------------

.. note::

   At build time, there are a number of source files that are autogenerated based
   on the configuration of the problem.  Most of these files are output into
   ``tmp_build_dir/castro_sources/Nd.COMP.OPTIONS.EXE/``, where ``N`` is the 
   dimensionality, ``COMP`` is the compiler name, and ``OPTIONS`` can be any
   number of options (``MPI``, ``DEBUG``, ...).

This is the current build system process.

* ``set_variables.py`` is called

  .. index:: set_variables.py, _variables, state_indices_nd.F90, state_indices.H

  * This processes the Castro ``_variables`` file and writes
    ``state_indices_nd.F90`` and ``state_indices.H`` into the
    ``tmp_build_dir/castro_sources/`` directory.

    These are used to define the size of the various state arrays and
    the integer keys to index each state variable.

  * The hook for this is in ``Make.auto_source`` in the build rule for ``state_indices_nd.F90``

  * You can test this portion of the build system by doing ``make test_variables``

* (for ``general_null networks``), ``actual_network.F90`` is created

  .. index:: write_network.py

  * This is done by ``write_network.py``

  * The hook for this is in ``$(CASTRO_HOME)/Microphysics/networks/general_null/Make.package``

* Runtime parameter files for the microphysics routines are parsed by ``write_probin.py``

  .. index:: write_probin.py

  * This writes the routines that manage the Microphysics runtime
    parameters: ``extern_parameters.cpp``, ``extern_parameters.H``, and  ``extern.F90``.  This is output in
    ``tmp_build_dir/castro_sources/``.

  * The hook for this is in ``Make.Castro`` in the rule for ``extern.F90``

* Castro's runtime parameters are parsed by ``parse_castro_params.py``

  .. index:: parse_castro_params.py

  * This writes the C++ header files that manage and read the runtime parameters and also
    creates the Fortran module ``meth_params.F90``, which defines all
    of the runtime parameters available to Fortran, from the template
    ``meth_params.template`` in ``Source/driver``. The file is output in
    ``tmp_build_dir/castro_sources/``.  

  * The hook for this is in ``Make.Castro`` in the rule for ``meth_params.F90``

* Problem-specific runtime parameters are parsed by ``write_probdata.py``

  * If the problem directory defines a ``_prob_params`` then it is parsed
    and used to C++ header and source files ``prob_parameters.H`` and ``prob_parameters.cpp``.
    These handle reading the ``problem.*`` parameters from the inputs file.
    Even without a problem-specific ``_prob_params``, all of the 
    variables in ``Castro/Source/problems/_default_prob_params`` will be included.

  * The script ``Castro/Util/scripts/write_probdata.py`` is used

  * The hook for this is in ``Make.auto_source`` in the ``prob_params_auto.F90`` rule.

  * The ``prob_params_auto.F90`` file is output into ``tmp_build_dir/castro_sources/``.

* The Fortran dependencies file is created

  * This creates the ``f90.depends`` file in the ``tmp_build_dir``

  * The script ``amrex/Tools/F_scripts/dep.py`` is used

  * The hook for this is in ``amrex/Tools/GNUMake/Make.rules`` in the
    ``$(depEXETempDir)/f90.depends`` target

* The C/C++ dependencies file is created

  * This creates the individual ``.d`` files in ``tmp_build_dir``, one for each source file

  * A set of rules in ``Make.rules`` handles this. There is some
    description of what each line does in the comments of the make
    file

* Output to stdout the git version of the sources, via
  ``describe_sources.py``.  This doesn’t affect the build process

For all of this to work, we need the ``tmp_build_dir/s`` directory to
be first in the vpath, so our modified sources are found and used.


