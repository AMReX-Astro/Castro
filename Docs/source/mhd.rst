.. _ch:mhd:

***
MHD
***

Castro implements a constrained transport (CT) corner transport upwind
(CTU) MHD scheme based on the work of Miniati & Martin
:cite:`miniati_martin`.  MHD is enabled by compiling with ``USE_MHD =
TRUE``.  This replaces the pure hydrodynamics solver, but uses the
same driver as the CTU hydrodynamics solver.  This means that all of
the source terms supported by the hydrodynamics solver are also
supported by MHD.

.. note::

   Currently the MHD solver is single-level only.  AMR support is forthcoming.
