**************************
Frequently Asked Questions
**************************

Compiling
=========

#. *Compiling fails giving me a cryptic message about a module not
   being found, usually bl_types or bl_error_module, like:*

   ::

       mpif90 -fno-range-check -fno-second-underscore
        -Jo/3d.Linux.gcc.gfortran.MPI.EXE -I o/3d.Linux.gcc.gfortran.MPI.EXE
        -ffixed-line-length-0 -g -O3
        -I. -I/home/zingale/development/Microphysics/util
        -I../../Microphysics/EOS -I../../Microphysics/EOS/gamma_law
        -I../../Microphysics/networks
        -I../../Microphysics/networks/general_null
        -I. -I/home/zingale/development/BoxLib//Src/C_BaseLib
        -I/home/zingale/development/BoxLib//Src/C_AMRLib
        -I/home/zingale/development/BoxLib//Src/C_BoundaryLib -I../../Source
        -I../../Source/Src_3d -I../../Source/Src_nd -I../../constants
        -I/home/zingale/development/BoxLib//Src/F_BaseLib
        -I/home/zingale/development/BoxLib//Tools/C_scripts -c
        ../../Microphysics/EOS/eos.f90 -o
        o/3d.Linux.gcc.gfortran.MPI.EXE/eos.o
       ../../Microphysics/EOS/eos.f90:3:6:

          use bl_types
             1
       Fatal Error: Can't open module file ‘bl_types.mod’ for reading at (1):
        No such file or directory compilation terminated.

       /home/zingale/development/BoxLib//Tools/C_mk/Make.rules:122: recipe
        for target 'o/3d.Linux.gcc.gfortran.MPI.EXE/eos.o' failed

   This usually indicates that the build system cannot find a source file
   (note: the problem above is not bl_types, that just seems to be
   the way the error manifests itself). The source files are specified
   in the various Make.package files throughout the
   Castro directory hierarchy. make will look through the
   directories in the VPATH_LOCATIONS to find the files.

   There are 2 things you can do to check what’s happening. First, inspect
   the directories in VPATH_LOCATIONS. This can be done via:

   ::

       make print-VPATH_LOCATIONS

   Next, ask make to tell you where it is finding each of the source
   files. This is done through a script find_files_vpath.py
   that is hooked into Castro’s build system. You can run this as:

   ::

       make file_locations

   At the end of the report, it will list any files it cannot find in
   the vpath. Some of these are to be expected (like extern.f90
   and buildInfo.cpp—these are written at compiletime. But any
   other missing files need to be investigated.

#. *I’m still having trouble compiling. How can I find out what
   all of the make variables are set to?*

   Use:

   ::

       make help

   This will tell you the value of all the compilers and their options.

#. *How do I use a system’s BLAS library instead of compiling and
   linking the one that comes with the StarKiller microphysics?*

   To use a system’s BLAS library, set the Make variable
   USE_SYSTEM_BLAS to TRUE. This will then look at
   the Make variable BLAS_LIBRARY for the library to link
   (defaults to -lopenblas).

#. *How can I check to make sure the function signatures defined
   in C are consistent with their implementations in Fortran?*

   Use:

   ::

       make typecheck

   This will compile the code and report on any mismatched function signatures.

Debugging
=========

#. *Castro crashes with a floating point exception—how can
   I get more information?*

   The best thing to do is to recompile the code with TEST=TRUE
   set in the GNUmakefile. This will have AMReX catch the
   signals raised in both C and Fortran functions. Behind the
   scenes, this defines the AMREX_TESTING preprocessor flag, which
   will initialize memory allocated in fabs or multifabs to
   signaling NaNs (sNaN), and use the BLBackTrace::handler()
   function to handle various signals raised in both C and Fortran
   functions. This is a Linux/UNIX capability. This gives us a chance
   to print out backtrace information. The signals include seg fault,
   floating point exceptions (NaNs, divided by zero and overflow), and
   interruption by the user and system. What signals are handed to
   AMReX are controlled by AMReX(e.g., using interruption by the
   user, this was once used to find an MPI deadlock.) It also includes
   the BL_ASSERTION statements if USE_ASSERTION=TRUE or
   DEBUG=TRUE.

   The AMReX parameters that affect the behavior are:

   -  amrex.fpe_trap_invalid

   -  amrex.fpe_trap_zero

   -  amrex.fpe_trap_overflow

   For further capabilities, defining BACKTRACE=TRUE enables you
   to get more information than the backtrace of the call stack info by
   instrumenting the code. (This is in C code only). Here is an
   example. You know the line “Real rho = state(cell,0);” is
   causing a segfault. You could add a print statement before that.
   But it might print out thousands (or even millions) of line before
   it hits the segfault. With BACKTRACE, you could do

   ::

             #ifdef AMREX_BACKTRACING
                std::ostringstream ss;
                ss << ``state.box() = `` << state.box() << `` cell = `` << cell;
                BL_BACKTRACE_PUSH(ss.str()); // PUSH takes std::string
             #endif

             Real rho = state(cell,0);  // state is a Fab, and cell is an IntVect.

   The “print” prints to a stack of string, not stdout. When it
   hits the segfault, you will only see the last print out.

#. *How can I monitor the state in a zone from the C side
   at various points in the evolution?*

   Given a MultiFab mf, you can dump out the state as:

   ::

           print_state(mf, IntVect(D_DECL(10, 20, 30)));

   Here, the IntVect has the dimension that we were compiled with
   (and this is handled through the preprocessor D_DECL). In
   this case, we are inspecting zone (10, 20, 30), in the global index
   space. Note that since a multifab exists only on a single level, the
   integer indices here refer to the global index space on that level.

#. *What if I want to see all the data in a FArrayBox?*

   You can simply output a FAB to std::cout. Imagine that you
   are in an MFIter loop, with a MultiFab mf:

   ::

           S = FArrayBox& mf[mfi];
           std::cout << S << std::endl;

   This will output the contents on the FAB, one zone per line.

Profiling
=========

#. *How can I get line-by-line profiling information?*

   With the GNU compliers, you can enabling profiling with gprof
   by compiling with

   ::

         USE_GPROF=TRUE

   in your GNUmakefile.

   When you run, a file named gmon.out will be produced. This can
   be processed with gprof by running:

   ::

         gprof exec-name

   where *exec-name* is the name of the executable. More detailed
   line-by-line information can be obtained by passing the -l
   argument to gprof.

Managing Runs
=============

#. *How can I force the running code to output, even it the plot or
   checkpoint interval parameters don’t require it?*

   Create a file called dump_and_continue, e.g., as:

   ::

       touch dump_and_continue

   This will force the code to output a checkpoint file that can be used
   to restart. Other options are plot_and_continue to output
   a plotfile, dump_and_stop to output a checkpoint file
   and halt the code, and stop_run to simply stop the code.
   Note that the parameter amr.message_int controls how often
   the existence of these files is checked; by default it is 10, so the
   check will be done at the end of every timestep that is a multiple of 10.
   Set that to 1 in your inputs file if you’d like it to check every timestep.

#. *How can I output plotfiles in single precision?*

   The AMReX runtime parameter:

   ::

       fab.format = NATIVE_32

   controls this (put this in your inputs file). Note: checkpoint files are unaffected
   by this and will always be written out in the native precision (the ‘fab.format‘ parameter
   is overridden in the checkpoint code in AMReX).

#. *How can I check the compilation parameters of a Castro executable?*

   The build information (including git hashes, modules, EoS, network, etc.) can be displayed by running the executable as 

   ::

       ./Castro.exe --display

.. _ch:faq:vis:

Runtime Errors
==============

#. *When running with retries, Castro requests too many substeps
   and crashes.*

   This can occur due to CFL violations or negative densities.  If
   there are density resets, try running with
   ``castro.limit_fluxes_on_small_dens`` = 1.  This will use a flux
   limiter to prevent the density from going negative.

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

   If you don’t already have an xorg.conf then you can create one
   by running nvidia-xconfig first.
