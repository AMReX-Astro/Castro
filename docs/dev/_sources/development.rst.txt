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

We github actions to run integration tests on the code and to build and deploy the documentation.

Currently, we run the `clang static analyzer <https://clang-analyzer.llvm.org/>`_, which finds potential bugs in the code. It also runs a script to convert any tabs in the code into spaces. Both of these are run on pull requests to the Castro GitHub repo, and are run weekly on the development branch.

