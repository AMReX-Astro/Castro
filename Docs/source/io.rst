.. _ch:io:

**********
Outputting
**********

Restart Capability
------------------

.. index:: amr.check_file, amr.check_int, amr.check_per, amr.restart
.. index:: amr.checkpoint_files_output, amr.check_nfiles, amr.checkpoint_on_restart
.. index:: castro.output_at_completion, castro.grown_factor

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

    See the section :ref:`sec:parallel_io` for more details on parallel I/O and the
    ``amr.check_nfiles`` parameter.

  * ``amr.checkpoint_on_restart``: should we write a
    checkpoint immediately after restarting? (0 or 1; default: 0)

  * ``castro.output_at_completion``: should we write a final checkpoint/plotfile? (0 or 1)

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


Plotfile Outputting
-------------------

.. index:: amr.plot_files_output, amr.plotfile_on_restart, amr.write_plotfile_with_checkpoint

Castro has two levels of plotfiles, `regular` plotfiles and `small`
plotfiles.  The idea behind this distinction is that we can output a
small number of variables very frequently in the small plotfiles and
output a large number (or all variables) less frequently.  This helps
keep the data sizes down while allowing for fine-grained temporal
analysis of important quantities.


A few general controls determines whether we want to output plotfiles and when:

  * ``amr.plot_files_output`` : this is set to 1 to output plotfiles

  * ``amr.plotfile_on_restart`` : set this to 1 to dump out a plotfile
    immediately when we restart.

  * ``amr.write_plotfile_with_checkpoint`` : always output a plotfile
    when we dump a checkpoint file.

.. index:: amr.plot_file, amr.plot_per, amr.plot_int

The frequency of outputting and naming of regular plotfiles is
controlled by:

  * ``amr.plot_file`` : this is the base name for the plotfile,
    e.g. ``plt``.

  * ``amr.plot_per`` : this is the amount of simulation time between
    plotfile output

    .. note:: ``amr.plot_per`` will write a plotfile at the first
       timestep whose ending time is past an integer multiple of this
       interval.  In particular, the timestep is not modified to match
       this interval, so you won’t get a checkpoint at exactly the time
       you requested.

  * ``amr.plot_int`` this is the number of timesteps between plotfiles.
    Set this to -1 to rely on the simulation-time-based outputting.

.. index:: amr.small_plot_file, amr.small_plot_per, amr.small_plot_int

Similarly, the frequency of outputting and naming of small plotfiles
is controlled by:

  * ``amr.small_plot_file`` : this is the base name for the small plotfile,
    e.g. ``smallplt``.

  * ``amr.small_plot_per`` : this is the amount of simulation time between
    small plotfile output

  * ``amr.small_plot_int`` this is the number of timesteps between small plotfiles.
    Set this to -1 to rely on the simulation-time-based outputting.

Additional output options control how the I/O is done:

  * ``amr.plot_nfiles``: how parallel is the writing of the
    plotfiles? (Integer :math:`\geq 1`; default: 64)

    See the Software Section for more details on parallel I/O and the
    ``amr.plot_nfiles`` parameter.

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


Controlling What’s in the PlotFile
----------------------------------

.. index:: amr.plot_vars, amr.derive_plot_vars

There are a few options that can be set at runtime to control what
variables appear in the regular plotfile.

  * ``amr.plot_vars``: this controls which of the main
    state variables appear in the plotfile. The default is for all of
    them to be stored. But you can specify a subset by name, e.g.::

        amr.plot_vars = density

    to only store that subset.

  * ``amr.derive_plot_vars``: this controls which of the derived
    variables to be stored in the plotfile. Derived variables are
    created only when the plotfile is being created, using the
    infrastructure provided by AMReX to register variables and the
    associated C++ routine to do the deriving.

    By default, no derived variables are stored. You can store all
    derived variables that Castro knows about by doing::

       amr.derive_plot_vars = ALL

   or a subset by explicitly listing them, e.g.::

      amr.derive_plot_vars = entropy pressure

   To not output any derived variable,s this is set to ``NONE``.

.. index:: amr.small_plot_vars

For small plotfiles, the controls that lists the variables is:

  * ``amr.small_plot_vars`` : this is a list of which variables
    to include in the small plotfile.

  * ``amr.derive_small_plot_vars`` : this is a list of which derived
    variables to include in the small plotfile.


Plotfile Variables
------------------

Native variables
^^^^^^^^^^^^^^^^

These variables come directly from the ``StateData``, either the
``State_Type`` (for the hydrodynamic variables), ``Reactions_Type``
(for the nuclear energy generation quantities). ``PhiGrav_Type`` and
``Gravity_Type`` (for the gravity quantities), and ``Rad_Type`` (for
radiation quantities).


+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``density``                       | Mass density, :math:`\rho`                        | :math:`\gcc`                         |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``xmom``                          | x-momentum, :math:`(\rho u)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``ymom``                          | y-momentum, :math:`(\rho v)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``zmom``                          | z-momentum, :math:`(\rho w)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_E``                         | Total energy density                              | :math:`{\rm erg~cm^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_e``                         | Internal energy density                           | :math:`{\rm erg~cm^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Temp``                          | Temperature                                       | :math:`{\rm K}`                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_X``                         | Mass density of species X                         | :math:`\gcc`                         |
| (where X is any of the species    |                                                   |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``omegadot_X``                    | Creation rate of species X                        | :math:`{\rm s^{-1}}`                 |
| (where X is any of the species    | :math:`\omegadot_k = DX_k/Dt`                     |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_enuc``                      | Nuclear energy generation rate density            | :math:`{\rm erg~cm^{-3}~s^{-1}}`     |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiGrav``                       | Gravitational potential                           | :math:`{\rm erg~g^{-1}}`             |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``grav_x``, ``grav_y``,           | Gravitational acceleration                        | :math:`{\rm cm~s^{-2}}`              |
| ``grav_z``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rmom``                          | Radial momentum (defined for                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
|                                   | ``HYBRID_MOMENTUM``)                              |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``lmom``                          | Angular momentum (:math:`\theta`; defined for     | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
|                                   | ``HYBRID_MOMENTUM``)                              |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pmom``                          | z-momentum (defined for ``HYBRID_MOMENTUM``)      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Shock``                         | Shock flag (= 1 if a zone has a shock;            | --                                   |
|                                   | defined for ``SHOCK``)                            |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rad``, ``rad0``, ``rad1``,      | Radiation energy density                          |                                      |
| ...                               | (for multigroup radiation, each group has its     |                                      |
|                                   | own variable)                                     |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+



Derived variables
^^^^^^^^^^^^^^^^^

.. index:: castro.domain_is_plane_parallel

+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| variable name                     | description                                       | derive routine              | units                                   |
+===================================+===================================================+=============================+=========================================+
| ``abar``                          | Mean atomic mass                                  | ``derabar``                 | :math:`\amu`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``angular_momentum_x``,           | Angular momentum / volume in the x, y, or z dir   | ``derangmomx``,             | :math:`{\rm g~cm^{-1}~s^{-1}}`          |
| ``angular_momentum_y``,           | computed as :math:`[(\rho \ub) \times {\bf r}]_n` | ``derangmomy``,             |                                         |
| ``angular_momentum_z``            | where :math:`{\bf r}` is the distance from        | ``derangmomz``              |                                         |
|                                   | ``center`` and :math:`n` is either x, y, or z     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``diff_coeff``                    | Thermal diffusion coefficient,                    | ``derdiffcoeff``            | :math:`{\rm cm^2~s^{-1}}`               |
|                                   | :math:`\kth/(\rho c_v)`                           |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``diff_term``                     | :math:`\nabla\cdot(\kth\nabla T)`                 | ``derdiffterm``             | :math:`{\rm erg~cm^{-3}~s^{-1}}`        |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``divu``                          | :math:`\nabla \cdot \ub`                          | ``derdivu``                 | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_e``                        | Specific internal energy computed from the        | ``dereint2``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | conserved :math:`(\rho e)` state variable as      |                             |                                         |
|                                   | :math:`e = (\rho e)/\rho`                         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_E``                        | Specific internal energy computed from the        | ``dereint1``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | total energy and momentum conserved state as      |                             |                                         |
|                                   | :math:`e=[(\rho E)-\frac{1}{2}(\rho \ub^2)]/\rho` |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``entropy``                       | Specific entropy, :math:`s`, computed as          | ``derentropy``              | :math:`{\rm erg~g^{-1}~K^{-1}}`         |
|                                   | :math:`s = s(\rho, e, X_k)`, where `e` is         |                             |                                         |
|                                   | computed from :math:`(\rho e)`                    |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``enuc``                          | Nuclear energy generation rate / gram             | ``derenuc``                 | :math:`{\rm erg~g^{-1}~s^{-1}}`         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Ertot``                         | Total radiation energy density                    | ``derertot``                |                                         |
|                                   | (for multigroup radiation problems)               |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Frcomx``, ``Frcomy``,           | Comoving radiation flux                           | ``Radiation.cpp``           |                                         |
| ``Frcomz``                        |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Frlabx``, ``Frlaby``,           | Lab-frame radiation flux                          | ``Radiation.cpp``           |                                         |
| ``Frlabz``                        |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Gamma_1``                       | Adiabatic index,                                  | ``dergamma1``               | --                                      |
|                                   | :math:`d\log p/d\log \rho|_s`                     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``kineng``                        | Kinetic energy density,                           | ``derkineng``               | :math:`{\rm erg~cm^{-3}}`               |
|                                   | :math:`K = \frac{1}{2} |(\rho \ub)|^2`            |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``lambda``                        | Radiation flux limiter                            |                             | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``logden``                        | :math:`\log_{10} \rho`                            | ``derlogten``               | dimensionless, assuming :math:`\rho`    |
|                                   |                                                   |                             | is in CGS                               |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``MachNumber``                    | Fluid Mach number, :math:`|\ub|/c_s`              | ``dermachnumber``           | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``maggrav``                       | Gravitational acceleration magnitude              | ``dermaggrav``              | :math:`{\rm cm~s^{-2}}`                 |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magmom``                        | Momentum density magnitude,                       | ``dermagmom``               | :math:`{\rm g~cm^{-2}~s^{-1}}`          |
|                                   | :math:`|\rho \ub|`                                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvel``                        | Velocity magnitude, :math:`|\ub|`                 | ``dermagvel``               | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvort``                       | Vorticity magnitude, :math:`|\nabla\times\ub|`    | ``dermagvort``              | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``pressure``                      | Total pressure, including ions, electrons,        | ``derpres``                 | :math:`{\rm dyn~cm^{-2}}`               |
|                                   | and radiation (for non radhydro problems)         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``radvel``                        | Radial velocity (measured with respect to         | ``derradialvel``            | :math:`\cms`                            |
|                                   | ``center`` or vertical axis if                    |                             |                                         |
|                                   | ``domain_is_plane_parallel`` is set)              |                             |                                         |
|                                   | :math:`(xu + yv + zw)/r`                          |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``circvel``                       | Circumferential velocity (perpendicular to        | ``derradialvel``            | :math:`\cms`                            |
|                                   | ``radvel``.  If ``domain_is_plane_parallel`` is   |                             |                                         |
|                                   | set, then this is in the x-y plane                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``soundspeed``                    | Sound speed                                       | ``dersoundspeed``           | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``StateErr``                      |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``thermal_cond``                  | Thermal conductivity, :math:`\kth`                | ``dercond``                 | :math:`{\rm erg~cm^{-1}~s^{-1}~K^{-1}}` |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``t_sound_t_enuc``                |                                                   | ``derenuctimescale``        | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``uminusc``                       | (only for 1D) x-velocity :math:`-` sound          | ``deruminusc``              | :math:`\cms`                            |
|                                   | speed                                             |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``uplusc``                        | (only for 1D) x-velocity + sound speed            | ``deruplusc``               | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``X(q)``                          | Mass fraction of species q                        | ``derspec``                 | --                                      |
|                                   | :math:`X_k = (\rho X_k)/\rho`                     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``x_velocity``,                   | Fluid velocity,                                   | ``dervel``                  | :math:`\cms`                            |
| ``y_velocity``,                   | :math:`\ub = (\rho \ub)/\rho`                     |                             |                                         |
| ``z_velocity``                    |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+

problem-specific plotfile variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the section on :ref:`Problem_Derives.H <problem_derives>` for more details about defining your own plotfile variables.

+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``analytic``                      |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pi``                            |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pioverp0``                      |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``primarymask``                   |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``secondarymask``                 |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Terror``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Texact``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_angular_momentum_x``,  |                                                   |                                      |
| ``inertial_angular_momentum_y``,  |                                                   |                                      |
| ``inertial_angular_momentum_z``   |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_momentum_x``,          |                                                   |                                      |
| ``inertial_momentum_y``,          |                                                   |                                      |
| ``inertial_momentum_z``           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_radial_momentum_x``,   |                                                   |                                      |
| ``inertial_radial_momentum_y``,   |                                                   |                                      |
| ``inertial_radial_momentum_z``    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEff``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEffPM_P``                    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEffPM_S``                    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``tpert``                         |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+



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


.. _sec:global_diag:

Global Diagnostics
-------------------

.. index:: castro.sum_interval, global diagnostics, amr_diag.out, gravity_diag.out, grid_diag.out, species_diag.out

Castro can calculate integrals of quantities on the grid and other
global quantities and output them to both the screen and to a runtime
file at regular intervals.  By default, this capability is off.  To
enable it, one of the following runtime parameters can be set:

  * ``castro.sum_interval``: if :math:`> 0`, how often (in level-0 time
    steps) to compute and print integral quantities (Integer; default: -1)

    The integral quantities include total mass, momentum and energy in
    the domain every ``castro.sum_interval`` level-0 steps.  The print
    statements have the form::

           TIME= 1.91717746 MASS= 1.792410279e+34

    for example.

  * ``castro.sum_per``: how often in simulation time to output
    integral quantities (this is used as an alternate to
    ``castro.sum_interval``).

By default, 4 output files are created:

  * ``amr_diag.out`` : This includes timestep information, in the
    following columns:

    * timestep
    * time
    * dt
    * finest level
    * coarse timestep walltime

  * ``gravity_diag.out`` : For problems with Poisson gravity, this
    includes the gravitational wave amplitudes

  * ``grid_diag.out`` : This includes integrals of the state data:

    * time
    * mass
    * x-, y-, and z-momentum
    * x-, y-, and z-angular momentum
    * kinetic energy
    * internal energy
    * kinetic + internal energy
    * gravitational potential energy
    * total energy (including gravitational potential energy)

  * ``species_diag.out`` : This contains the mass of each of the nuclear species on the grid.

    .. note::

       The species masses are given in units of solar masses.

``Castro/Util/scripts/diag_parser.py`` contains Python code for parsing
these output files into Numpy arrays.  Usage instructions are included
in the file, along with an example script at
``Castro/Util/scripts/plot_species.py``.  This reads a
``species_diag.out`` file provided on the command line and makes a plot
of the total mass fractions over time.

Some problems have custom versions of the diagnostics with additional
information.  These are not currently supported by the Python parser.


.. _sec:parallel_io:

Parallel I/O
------------

Both checkpoint files and plotfiles are really directories containing
subdirectories: one subdirectory for each level of the AMR hierarchy.
The fundamental data structure we read/write to disk is a ``MultiFab``,
which is made up of multiple FAB’s, one FAB per grid. Multiple
``MultiFab`` s may be written to each directory in a checkpoint file.
``MultiFab`` s of course are shared across CPUs; a single ``MultiFab`` may be
shared across thousands of CPUs. Each CPU writes the part of the
``MultiFab`` that it owns to disk, but they don’t each write to their own
distinct file. Instead each MultiFab is written to a runtime
configurable number of files :math:`N` (:math:`N` can be set in the inputs file as the
parameter ``amr.checkpoint_nfiles`` and ``amr.plot_nfiles``; the
default is 64). That is to say, each ``MultiFab`` is written to disk
across at most :math:`N` files, plus a small amount of data that gets written
to a header file describing how the file is laid out in those :math:`N` files.

What happens is :math:`N` CPUs each opens a unique one of the :math:`N` files into
which the ``MultiFab`` is being written, seeks to the end, and writes
their data. The other CPUs are waiting at a barrier for those :math:`N`
writing CPUs to finish. This repeats for another :math:`N` CPUs until all the
data in the ``MultiFab`` is written to disk. All CPUs then pass some data
to CPU 0 which writes a header file describing how the ``MultiFab`` is
laid out on disk.

We also read ``MultiFabs`` from disk in a “chunky” manner, opening only :math:`N`
files for reading at a time. The number :math:`N`, when the ``MultiFab`` s were
written, does not have to match the number :math:`N` when the ``MultiFab`` s are
being read from disk. Nor does the number of CPUs running while
reading in the ``MultiFab`` need to match the number of CPUs running when
the ``MultiFab`` was written to disk.

Think of the number :math:`N` as the number of independent I/O pathways in
your underlying parallel filesystem. Of course a “real” parallel
filesystem should be able to handle any reasonable value of :math:`N`. The
value -1 forces :math:`N` to the number of CPUs on which you’re
running, which means that each CPU writes to a unique file, which can
create a very large number of files, which can lead to inode issues.
