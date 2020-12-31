******************
Runtime Parameters
******************

Introduction to Runtime Parameters
==================================

Castro has 2 sets of runtime parameters—those controlled by
C++ and those controlled by Fortran. The C++ parameters are set
in the inputs file and managed by the AMReX ``ParmParse``
class. For Castro-specific parameters, we list the runtime
parameters in a file ``_cpp_parameters`` and generate the
C++ code and headers at compile time.

The behavior of the network, EOS, and other microphysics routines are
controlled by a different set of runtime parameters. These parameters are defined
in plain-text files ``_parameters`` located in the different
directories that hold the microphysics code. At compile time, a
script in the AMReX bulid system, ``findparams.py``, locates all
of the ``_parameters`` files that are needed for the given choice
of network, integrator, and EOS, and assembles all of the runtime
parameters into a module named ``extern_probin_module`` (using the
``write_probin.py`` script). The parameters are set in your
probin file in the ``&extern`` namelist.

C++ parameter format
--------------------

The C parameters take the form of::

    # comment describing the parameter
    name   type   default   need in Fortran?   ifdef

Here,

  * `name` is the name of the parameter that will be looked for
    in the inputs file.

    The name can actually take the form of ``(a, b)``, where ``a`` is
    the name to be used in the inputs file where the parameter is set
    and ``b`` is the name used within the Castro C++ class.  It is not
    recommended to name new parameters with this functionality—this
    was implemented for backwards compatibility.


  * `type` is one of int, Real, or string

  * `default` is the default value of the parameter.

The next columns are optional, but you need to fill in all of the
information up to and including any of the optional columns you need
(e.g., if you are going to provide the fortran name, you also need to
provide "need in Fortran?" and "ifdef".

  * `need in Fortran?` is ``y`` if the runtime parameter should be
    made available in Fortran (through ``meth_params_module``).

  * `ifdef` provides the name of a preprocessor name that should
    wrap this parameter definition—it will only be compiled in if that
    name is defined to the preprocessor.

Finally, any comment (starting with ``#``) immediately before the
parameter definition will be used to generate the documentation
describing the parameters.

Microphysics/extern parameter format
------------------------------------

The microphysics/extern parameter definitions take the form of::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the `priority` is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.


Removed Runtime Parameters
==========================

The following runtime parameters have been removed for Castro.

*  ``castro.ppm_flatten_before_integrals`` : this parameter
   controlled whether we applied the flattening of the parabolic
   profiles before we integrated under their profiles or afterwards.
   The default was switched to flattening before the integration,
   which is more consistent with the original PPM methodology. This
   parameter was removed since the variation enabled by this parameter
   was not that great.

   (removed in commit: 9cab697268997714919de16db1ca7e77a95c4f98)

* ``castro.ppm_reference`` and ``castro.ppm_reference_edge_limit`` :
  these parameters controlled whether we use the integral under the
  parabola for the fastest wave moving toward the interface for the
  reference state and whether in the case that the wave is moving away
  from the interface we use the cell-center value or the limit of the
  parabola on the interface.  These were removed because there is
  little reason to not use the reference state.

   (removed in commit: 213f4ffc53463141084551c7be4b37a2720229aa)


Parameters by Namespace
=======================

.. toctree::

   runtime_parameters
