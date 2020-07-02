**************************
Development Best Practices
**************************

Coding Conventions
==================

Castro development should do its best to adhere to the following coding
style guidelines.

General
-------

* Indentation should be spaces, not tabs, with 4 spaces preferred.


C++
---

* Conditional should appear as:

  .. code:: c++

     if (condition)
     {
         ...
     }
     else if (condition)
     {
         ...
     }
     else
     {
         ...
     }


C++ to Fortran
--------------

* All C to Fortran interfaces should use the ISO C binding.  The
  Fortran subroutine should use

  .. code:: fortran

     subroutine subroutine_name(args) bind(C, name="subroutine_name)

* Data passed by reference from C++ should use ``*`` and not ``&``.

* Scalars should be passed by value, using the Fortran ``value`` attribute.

* Fortran routines that are called from C++ should being with ``ca_``.
  This is true even if these routines are also called by other
  Fortran.

* A separate tile ``lo`` and ``hi`` should be passed in to each
  subroutine, with the expectation that the Fortran subroutine will
  act exactly over the domain (lo, hi). If you want to include ghost
  zones in the calculation, include them in your box using
  ``mfi.growntilebox(ng)`` instead of ``mfi.validbox()``.


Fortran
-------

* Put all routines in a module—this ensures that argument lists are
  checked.

* Use the ``only`` clause in module ``use`` statements to explicitly
  make clear what is being accessed.

* In a module, there should be no "top-level" ``use`` statements (with
  the exception of getting access to the ``rt`` type).  Instead each
  function / subroutine in the module  should use what it needs directly.

* New Fortran files should have the .F90 file extension, not the .f90
  file extension, so that they can be preprocessed.


Documentation
-------------

C++ routines should use Doxygen style comments.

* C++ functions should be documented in the header file using the style:

  .. code:: c++

     ///
     /// Description of the function
     ///
     /// @param bar       Brief description of the variable
     ///
     void foo(int bar) { ...

* Member variables can either be documented using the above style of comment block or
  with a brief inline description:

  .. code:: c++

     int var; ///< Brief description after the variable

Fortran functions should use Sphinx style comments

* Fortran functions should be documented by placing a comment block
  immediately after their prototype (i.e. `without` a line in betwen ) using the style:

  .. code:: fortran

     subroutine foo(bar)
       ! Description of the function

       use some_module

       implicit none

       integer, intent(inout) :: bar   ! Brief description of bar
       ...

  Documentation for modules should be similarly formatted, with the comment block again
  coming `immediately` after the module definition.

Castro Releases
===============

This outlines the procedure for doing the monthly Castro release.

Castro uses submodules for dependencies, this means that, at a
minimum, we must update the AMReX and  Microphysics submodules monthly when we
issue new releases. The releases for AMReX and Microphysics must be done
first. Then navigate to each submodule directory, checkout the new
tag, and then from the top-level directory of Castro do a "git add" on
the ``external/`` directory to store the new tags. So, for example, at
the beginning of March 2020 we would first issue the ``20.03`` tag on
Microphysics, and wait for AMReX to release a ``20.03`` tag, then do::

   cd $CASTRO_HOME/external
   cd amrex
   git pull
   git checkout 20.03
   cd ..
   cd Microphysics
   git pull
   git checkout 20.03
   cd ..
   git add -u .
   git commit -m "Update AMReX and Microphysics to release 20.03"

Then we can proceed with issuing our own release.


Each month the ``development`` branch is merged into ``main`` and a
release is tagged.  This needs to be done in coordination with its
submodule dependencies.  The general procedure is:

  * Do 'git pull' in both main and development branches.  (Use `git
    checkout xxx` to change to branch xxx.

  * In main branch, do `git merge development`.  Fix any conflicts
    if there are any.  (There should not be any conflicts unless a
    commit is checked into main directly without going through
    development.)

  * In main branch, commit new release notes (``CHANGES.md``)
    summarizing changes since last major release.

  * Tag the new release: ``git tag -m "Castro YY.MM" YY.MM``

  * ``git push``

  * ``git push --tags``

  * ``git checkout development``

  * ``git merge main``

  * ``git push``


Interim updates
---------------

When breaking changes to Microphysics occur in its development branch
that Castro depends on, we must update the Microphysics submodule on
the Castro development branch in the same way, replacing the git
checkout statement with the latest commit hash on the Microphysics
development branch. (A git submodule always tracks a specific
commit/tag on the target repo -- it is not configured to automatically
track a particular branch.)  Since such breaking changes usually are
accompanied by a Castro change, it is best practice to ensure that
the PRs in both Microphysics and Castro have been approved, then
merge the Microphysics PR, then add the update to the Microphysics
submodule to the Castro PR, then merge. A similar process applies for AMReX.


Continuous Integration
======================

We use Travis CI to run integration tests on the code and to build and deploy the documentation. The current status of these tests on the development branch can be found here:

.. image:: https://travis-ci.com/AMReX-Astro/Castro.svg?branch=development
   :target: https://travis-ci.com/AMReX-Astro/Castro

Currently, travis runs the `clang static analyzer <https://clang-analyzer.llvm.org/>`_, which finds potential bugs in the code. It also runs a script to convert any tabs in the code into spaces. Both of these are run on pull requests to the Castro GitHub repo, and are run weekly on the development branch. 

The travis build settings can be found in the ``.travis.yml`` file.