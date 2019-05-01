.. _sponge_section:

******
Sponge
******

Castro uses a sponge source term to damp velocities in regions where
we are not interested in the dynamics.  Often these are buffer regions
between the domain of interest and the boundary conditions.

There are three sponges, each controlled by different runtime parameters:

  * **radial sponge** : The radial sponge operates beyond some radius
    (define in terms of the ``center(:)`` position for the problem).

    It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & r < r_\mathrm{lower} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (r - r_\mathrm{lower})}{\Delta r} \right ]  & r_\mathrm{lower} <= r < r_\mathrm{upper} \\
                     f_\mathrm{upper} & r >= r_\mathrm{upper} \end{array} \right .



   The purpose of the sponge is to damp velocities outside of a star, to
   prevent them from dominating the timestep constraint. The sponge parameters
   are set in your ``probin`` file, in the ``&sponge`` namelist. You can sponge either
   on radius from the center (using ``sponge_lower_radius`` and
   ``sponge_upper_radius``) or on density (using ``sponge_lower_density``
   and ``sponge_upper_density``). The timescale of the damping is
   set through ``sponge_timescale``.
