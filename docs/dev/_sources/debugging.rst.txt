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

If the code crashes, it will produce one or more ``Backtrace.*``
files.  Looking over these files should pinpoint where the FPE
occurred allowing you to do further debugging.

You can output more information into the ``Backtrace.*`` files by
pushing it to the backtrace stack as described here:
:ref:`debugging_backtrace`.

Make sure your runtime options are valid
----------------------------------------

Castro can validate the runtime options specified in the inputs file
by running with ``castro.abort_on_invalid_params = 1``.


Valgrind
--------

We frequently run Castro with valgrind to find illegal memory
accesses.  The valgrind documentation can give details on how to use
it.


Clang-tidy
----------

.. index:: clang-tidy

We run `clang-tidy <https://clang.llvm.org/extra/clang-tidy/>`_ on all
pull requests using a `GitHub action
<https://github.com/AMReX-Astro/cpp-linter-action>`_. ``clang-tidy``
analyzes the source code, produces warnings for potential bugs and
offers suggestions for performance improvements.

It can also be run locally. Support for this is enabled via the AMReX build system,
and requires that you have  ``clang-tidy`` installed locally.  You build via:

.. prompt:: bash

   make USE_CLANG_TIDY=TRUE

and this will use the clang-tidy options set in
``Castro/.clang-tidy``.  If you do a parallel build, you should use
the ``-O`` flag to ensure that output is not mixed between files.

You can also use ask it to fix errors automatically via:

.. prompt:: bash

   make USE_CLANG_TIDY=TRUE CLANG_TIDY="clang-tidy --fix-errors"

and you can treat warnings as errors by adding ``CLANG_TIDY_WARN_ERROR=TRUE``.

.. note::

   Building a Castro problem with ``clang-tidy`` will suppress the
   checks in AMReX and Microphysics sources and headers.  This is set by
   the parameter ``CLANG_TIDY_IGNORE_SOURCES`` in ``Make.Castro``, and
   the ``HeaderFilterRegex`` whitelist in ``.clang-tidy``.


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



Physics issues
==============
