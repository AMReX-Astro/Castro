.. _ch:retry:

***************
Retry Mechanism
***************

Castro's Strang CTU solver has a retry mechanism that can discard a
time step on a level and restart with a smaller timestep, subcycling
within the level to make up the full time step needed for that level.
It is enabled by setting::

   castro.use_retry = 1

.. note::

   The Castro retry mechanism is enabled by default for CTU + Strang
   and Simplified SDC integration.

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

    This instructs the integration routine in Microphysics to not
    abort when the integration fails, but instead to tell the calling
    Castro routine that the integration failed so Castro can handle
    the retry itself.

    .. note::

       The combination of ``use_retry = 0`` and ``abort_on_failure = F``
       is unsafe and not supported.

       For true SDC, we disable retry and reset ``abort_on_failure`` to
       always be true, since retry is not supported for that integration.


