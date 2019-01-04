This is the current build system process, for a CUDA build.

* ``set_variables.py`` is called

  * this processes the Castro ``_variables`` file and writes
    ``set_conserved.H``, ``set_indices.F90``, and ``state_sizes.f90`` into the
    build directory.

    These are used to define the size of the various state arrays.

  * This is not GPU specific.

  * The hook for this is in ``Make.Castro`` in the build rule for ``set_indices.F90``

* (for ``general_null networks``), ``actual_network.F90`` is created

  * This is done by ``write_network.py``

  * The hook for this is in ``$(CASTRO_HOME)/Microphysics/networks/general_null/Make.package``

* Runtime parameter files are parsed by ``write_probin.py``

  * The hook for this is in ``Make.Castro`` in the rule for ``extern.F90``

* The Fortran dependencies file is created

  * This creates the ``f90.depends`` file in the ``tmp_build_dir``

  * The script ``amrex/Tools/F_scripts/dep.py`` is used

  * The hook for this is in ``amrex/Tools/GNUMake/Make.rules`` in the
    ``$(depEXETempDir)/f90.depends`` target

* The C/C++ dependencies file is created

  * This creates the individual ``.d`` files in ``tmp_build_dir``, one for each source file

  * A set of rules in ``Make.rules`` handles this. There is some
    description of what each line does in the comments of the make
    file

* Interpret the ``#pragma gpu``

  * The script ``write_cuda_headers.py`` (in ``amrex/Tools/F_scripts/``) is tasked with
    understanding our custom pragma.  Its flow is:

    * Loop over all C++ files, looking for routines that are marked
      with the pragma and return a dict keyed by the name of the
      function with values being a list of the arguments

    * Parse the headers

      * preprocess all of the ``.H`` files to the ``tmp_build_dir/s/``
        directory, giving them the prefix ``CPP-``.

      * now parse the preprocessed headers, grab the function
        signatures there, modify them with the CUDA launch, and insert
        them into a copy of the original, unpreprocessed
        header.  These new copies are also put in ``tmp_build_dir/s/``.

      * loop through the C++ files that had the pragma, and add the
        needed launch macro.  These new ``.cpp`` files are put in the same
        ``tmp_build_dir/s/`` directory.

* Output to stdout the git version of the sources, via
  ``describe_sources.py``.  This doesn’t affect the build process

* Create device and host versions of each needed Fortran file. This
  is done as each ``.F90`` file is compiled with a rule in ``Make.rules`` that
  invokes ``gpu_fortran.py`` and then directs the compilation to build
  that version.

  * We look for a ``!$gpu`` comment in routines, and use that as an
    indication to mark it up with a host and device version of the
    routine

  * The modified ``.F90`` files are placed in ``tmp_build_dir/s/``

For all of this to work, we need the ``tmp_build_dir/s`` directory to
be first in the vpath, so our modified sources are found and used.


