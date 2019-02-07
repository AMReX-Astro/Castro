*********
Debugging
*********

There are several methods we typically use to debug issues in Castro.
Some descriptions are given below.

Using the compiler of other tools
=================================

Compiler checks
---------------

Recompile the code in debug-mode by setting::

   DEBUG := TRUE

in the ``GNUmakefile`` (or invoking ``make`` as ``make DEBUG=TRUE``).
This will create an executable with bounds checking and other compiler
options enabled.  Running this sometimes will show issues.


Check for NaNs
--------------

Compiling in debug mode also initializes uninitialized variables to
NaN.  For optimized code, the same can be done by setting::

   TEST := TRUE

in the ``GNUmakefile``.  To capture the NaNs, use the runtime parameter::

   amr.fpe_trap_invalid=1


Valgrind
--------

We frequently run Castro with valgrind to find illegal memory
accesses.  The valgrind documentation can give details on how to use
it.


Thread Sanitizer
----------------



Instrumenting the Code
======================



Physics issues
==============


