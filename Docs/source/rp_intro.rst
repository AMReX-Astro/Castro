******************
Runtime Parameters
******************

Introduction to Runtime Parameters
==================================

Castro runtime parameters are set
in the inputs file and managed by the AMReX ``ParmParse``
class. For Castro-specific parameters, we list the runtime
parameters in a file ``_cpp_parameters`` and generate the
C++ code and headers at compile time.

The behavior of the network, EOS, and other microphysics routines are
controlled by a different set of runtime parameters. These parameters are defined
in plain-text files ``_parameters`` located in the different
directories that hold the microphysics code. At compile time, a
a make function locates all
of the ``_parameters`` files that are needed for the given choice
of network, integrator, and EOS, and creates the ``extern_parameters.H`` and
``extern_parameters.cpp`` files that manage the parameters.  These
are set at runtime via the inputs file.

Castro-specific parameters
--------------------------

The Castro parameters that control the behavior of the code and
physics modules are listed in ``_cpp_parameters`` and take the form of::

    # comment describing the parameter
    name   type   default   ifdef

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

The next column is optional:

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


Parameters by Namespace
=======================

.. toctree::

   runtime_parameters
