Retry Mechanism
===============

Castro's Strang CTU solver has a retry mechanism that can discard a
time step on a level and restart with a smaller timestep, subcycling
within the level to make up the full time step needed for that level.
It is enabled by setting::

   castro.use_retry      = 1
   castro.hard_cfl_limit = 0

The number of subcycles to try in the level is controlled via the
``castro.max_subcycles`` parameter.  It is not really suggested to go
beyond ``16``---any more is usually an indication of a bigger problem.

A retry can be triggered by a number of conditions:

  * Exceeding the CFL condition for a level

  * A negative density is encountered

  * Integration failure in the burner

    Note: this requires that the following be set in your ``&extern``
    namelist::

      retry_burn = F
      abort_on_failure = F

  



