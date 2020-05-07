*********
Debugging
*********

There are several methods we typically use to debug issues in Castro.
Some descriptions are given below.

Using the compiler and other tools
==================================

Compiler checks
---------------

Recompile the code in debug-mode by setting::

   DEBUG := TRUE

in the ``GNUmakefile`` (or invoking ``make`` as ``make DEBUG=TRUE``).
This will create an executable with bounds checking and other compiler
options enabled.  Running this sometimes will show issues.


Check for invalid floating point exceptions
-------------------------------------------

Compiling in debug mode also initializes uninitialized variables to
NaN.  For optimized code, the same can be done by setting::

   TEST := TRUE

in the ``GNUmakefile``.  To capture the NaNs, use the runtime parameter::

   amrex.fpe_trap_invalid=1


Valgrind
--------

We frequently run Castro with valgrind to find illegal memory
accesses.  The valgrind documentation can give details on how to use
it.


Clang-tidy
----------

We run ```clang-tidy`` <https://clang.llvm.org/extra/clang-tidy/>`_ on all pull requests using a `GitHub action <https://github.com/AMReX-Astro/cpp-linter-action>`_. ``clang-tidy`` analyzes the source code, produces warnings for potential bugs and offers suggestions for performance improvements. 

It can also be run locally. This requires the ``clang-tidy`` and ``bear`` packages, and the python script
```run-clang-tidy.py` <https://github.com/AMReX-Astro/cpp-linter-action/blob/master/run-clang-tidy.py>`_. The analysis is performed by first compiling a problem using the ``bear`` package, then running the python script to analyze the source files. From within a problem directory, run

.. code-block:: bash

    bear make -j 20 USE_OMP=FALSE USE_MPI=FALSE DEBUG=TRUE 

    python3 run-clang-tidy.py -header-filter='Castro' -ignore-files='amrex|Microphysics' -j 20 > clang-tidy-report.txt

The compiler flags can be modified to suit the problem to be analyzed, but the ``DEBUG`` flag must be set to ``TRUE``. The ``header-filter`` option for the python script tells the script to only analyze header files containing the given regex pattern, and the ``ignore-files`` flag tells it to ignore any source files containing the given regex pattern. The ``-j`` option tells the script to run a given number of processes in parallel. The output is then redirected to a text file. 

Thread sanitizer
----------------



Instrumenting the Code
======================

Checking for NaNs
-----------------

In the C++ code, you can check whether a FAB contains NaNs using
the ``contains_nan()`` method:

.. code-block:: c++

   for (MFIter mfi(S, True); mfi.isValid(); ++mfi) {

     const Box& bx = mf.tilebox()

     // do some operations on S[mfi]

     if (S[mfi].contains_nan()) {
       amrex::Abort("S has NaNs")
     }
   }

There are other versions of ``contains_nan()`` that can take a Box
to operate over.

For Fortran code, the module ``nan_check`` has the function
``check_for_nan()`` that can be called in a Fortran routine to look
for NaNs:

.. code-block:: fortran

  subroutine check_for_nan(s, s_lo, s_hi, lo, hi, ncomp, comp)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: ncomp
    real(rt), intent(in) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)



Physics issues
==============


