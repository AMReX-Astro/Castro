***************
Castro Releases
***************

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


Each month the ``development`` branch is merged into ``master`` and a
release is tagged.  This needs to be done in coordination with its
submodule dependencies.  The general procedure is:

  * Do 'git pull' in both master and development branches.  (Use `git
    checkout xxx` to change to branch xxx.

  * In master branch, do `git merge development`.  Fix any conflicts
    if there are any.  (There should not be any conflicts unless a
    commit is checked into master directly without going through
    development.)

  * In master branch, commit new release notes (``CHANGES.md``)
    summarizing changes since last major release.

  * Tag the new release: ``git tag -m "Castro YY.MM" YY.MM``

  * ``git push``

  * ``git push --tags``

  * ``git checkout development``

  * ``git merge master``

  * ``git push``


Interim updates
---------------

When breaking changes to Microphysics occur in its development branch
that Castro depends on, we must update the Microphysics submodule on
the Castro development branch in the same way, replacing the git
checkout statement with the latest commit hash on the Microphysics
development branch. A git submodule always tracks a specific
commit/tag on the target repo -- it is not configured to automatically
track a particular branch. A similar process applies for AMReX.
