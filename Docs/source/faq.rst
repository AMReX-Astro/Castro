**************************
Frequently Asked Questions
**************************

Compiling
=========

#. *Compiling fails giving me a cryptic message about a module not
   being found.*

   This usually indicates that the build system cannot find a source file.
   The source files are specified
   in the various ``Make.package`` files throughout the
   Castro directory hierarchy. make will look through the
   directories in the ``VPATH_LOCATIONS`` to find the files.

   There are 2 things you can do to check what’s happening. First, inspect
   the directories in ``VPATH_LOCATIONS``. This can be done via:

   .. prompt:: bash

      make print-VPATH_LOCATIONS

   Next, ask make to tell you where it is finding each of the source
   files. This is done through a script ``find_files_vpath.py``
   that is hooked into Castro’s build system. You can run this as:

   .. prompt:: bash

      make file_locations

   At the end of the report, it will list any files it cannot find in
   the vpath. Some of these are to be expected (like
   ``buildInfo.cpp``—these are written at compile-time). But any other
   missing files need to be investigated.

#. *I put a copy of one of the header files (e.g. ``problem_tagging.H``)
   in my problem setup but it does not seem to recognized / used by
   the build system.  Why doesn't my executable use my custom version
   of the header?*

   This is likely due to compiler caching / ccache.  You need to
   clear the cache and the build:

   .. prompt:: bash

      ccache -C
      make clean

   Then rebuild and it should be recognized.

#. *I’m still having trouble compiling. How can I find out what
   all of the make variables are set to?*

   Use:

   .. prompt:: bash

      make help

   This will tell you the value of all the compilers and their options.

.. _debugging_backtrace:

Debugging
=========

#. *Castro crashes with a floating point exception—how can
   I get more information?*

   The best thing to do is to recompile the code with ``TEST=TRUE``
   set in the ``GNUmakefile``. This will have AMReX catch the
   signals raised in C++ functions. Behind the
   scenes, this defines the ``AMREX_TESTING`` preprocessor flag, which
   will initialize memory allocated in fabs or multifabs to
   signaling NaNs (sNaN), and use the ``BLBackTrace::handler()``
   function to handle various signals raised in C++
   functions. This is a Linux/UNIX capability. This gives us a chance
   to print out backtrace information. The signals include seg fault,
   floating point exceptions (NaNs, divided by zero and overflow), and
   interruption by the user and system. What signals are handed to
   AMReX are controlled by AMReX(e.g., using interruption by the
   user, this was once used to find an MPI deadlock.) It also includes
   the ``AMREX_ASSERTION`` statements if ``USE_ASSERTION=TRUE`` or
   ``DEBUG=TRUE``.

   The AMReX parameters that affect the behavior are:

   -  ``amrex.fpe_trap_invalid``

   -  ``amrex.fpe_trap_zero``

   -  ``amrex.fpe_trap_overflow``

   For further capabilities, you can get
   more information than the backtrace of the call stack info by
   instrumenting the code.  Here is an
   example. You know the line ``Real rho = state(cell,0);`` is
   causing a segfault. You could add a print statement before that.
   But it might print out thousands (or even millions) of line before
   it hits the segfault. Instead, you could

   .. code:: c++

             std::ostringstream ss;
             ss << "state.box() = " << state.box() << " cell = " << cell;
             BL_BACKTRACE_PUSH(ss.str()); // PUSH takes std::string

             Real rho = state(cell,0);  // state is a Fab, and cell is an IntVect.

   The "print" prints to a stack of string, not stdout. When it hits
   the segfault, you will only see the last print out in the backtrace
   file (e.g. ``BackTrace.0``).

   You may need to include the header ``AMReX_BLBackTrace.H``.

#. *How can I monitor the state in a zone from the C side
   at various points in the evolution?*

   Given a MultiFab ``mf``, you can dump out the state as:

   .. code:: c++

           print_state(mf, IntVect(AMREX_D_DECL(10, 20, 30)));

   Here, the IntVect has the dimension that we were compiled with
   (and this is handled through the preprocessor ``AMREX_D_DECL``). In
   this case, we are inspecting zone (10, 20, 30), in the global index
   space. Note that since a multifab exists only on a single level, the
   integer indices here refer to the global index space on that level.

#. *What if I want to see all the data in a FArrayBox?*

   You can simply output a FAB to ``std::cout``. Imagine that you
   are in an MFIter loop, with a MultiFab ``mf``:

   .. code:: c++

           S = FArrayBox& mf[mfi];
           std::cout << S << std::endl;

   This will output the contents on the FAB, one zone per line.

Profiling
=========

#. *How can I get line-by-line profiling information?*

   With the GNU compilers, you can enabling profiling with gprof
   by compiling with

   ::

         USE_GPROF=TRUE

   in your ``GNUmakefile``.

   When you run, a file named ``gmon.out`` will be produced. This can
   be processed with gprof by running:

   .. prompt:: bash

      gprof exec-name

   where *exec-name* is the name of the executable. More detailed
   line-by-line information can be obtained by passing the -l
   argument to gprof.

#. *How can I use AMReX's profiling to see which functions dominate when
   running in parallel?*

   You can compile with:

   .. prompt:: bash

      TINY_PROFILE=TRUE

   then at the end of the simulation, a table of profiling information will
   be produced.  See `AMReX Profiling Tools <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools.html>`_ for more information.

Managing Runs
=============

#. *How can I force the running code to output, even it the plot or
   checkpoint interval parameters don’t require it?*

   Create a file called ``dump_and_continue``, e.g., as:

   .. prompt:: bash

      touch dump_and_continue

   This will force the code to output a checkpoint file that can be used
   to restart. Other options are ``plot_and_continue`` to output
   a plotfile, ``dump_and_stop`` to output a checkpoint file
   and halt the code, and ``stop_run`` to simply stop the code.


   .. note::

      The parameter ``amr.message_int`` controls how often the
      existence of these files is checked; by default it is 1, so the
      check will be done at the end of every timestep, but you can
      set it to some other integer to check only timesteps that are a
      multiple of that number.

#. *How can I output plotfiles in single precision?*

   The AMReX runtime parameter:

   ::

       fab.format = NATIVE_32

   controls this (put this in your inputs file). Note: checkpoint files are unaffected
   by this and will always be written out in the native precision (the ‘fab.format‘ parameter
   is overridden in the checkpoint code in AMReX).

#. *How can I check the compilation parameters of a Castro executable?*

   The build information (including git hashes, modules, EoS, network, etc.) can be displayed by running the executable as

   .. prompt:: bash

      ./Castro.exe --describe

.. _ch:faq:vis:

Runtime Errors
==============

.. index:: castro.limit_fluxes_on_small_dens, castro.state_interp_order,
           castro.abundance_failure_tolerance, castro.abundance_failure_rho_cutoff

#. *When running with retries, Castro requests too many substeps
   and crashes.*

   This can occur due to CFL violations or negative densities.  If
   there are density resets, try running with
   ``castro.limit_fluxes_on_small_dens = 1``.  This will use a flux
   limiter to prevent the density from going negative.

#. *There might be a problem when Castro tries to normalize mass fractions
   and encounters: ``Invalid mass fraction in Castro::normalize_species()``.*

   If the error happens at the beginning of the timestep, it is possible that
   something unexpected happened durng the interpolation from the coarse-level
   to the fine-level. Try to set ``castro.state_interp_order = 0`` in the
   input file. This allows piecewise constant refinement, but sacrifices
   some benefit of the refinement.

   If the error continues, try to increase the tolerance of determining
   specie abundance validity check by setting ``castro.abundance_failure_tolerance``
   to a higher value, or increasing the density floor below which this is
   ignored by changing ``castro.abundance_failure_rho_cutoff``.

Visualization
=============

#. *When I try to use Amrvis with the Nvidia driver, all I see is
   black—no data. How do I fix this?*

   You need to edit your xorg.conf file (usually found in /etc/X11/
   to enable the Dac8Bit option. The section will look like:

   ::

       Section "Device"
           Identifier     "Device0"
           Driver         "nvidia"
           VendorName     "NVIDIA Corporation"
           Option         "Dac8bit" "True"
       EndSection

   If you don’t already have an ``xorg.conf`` then you can create one
   by running ``nvidia-xconfig`` first.
