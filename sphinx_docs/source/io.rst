.. _ch:io:

**********
Outputting
**********

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
filesytem should be able to handle any reasonable value of :math:`N`. The
value -1 forces :math:`N` to the number of CPUs on which you’re
running, which means that each CPU writes to a unique file, which can
create a very large number of files, which can lead to inode issues.
