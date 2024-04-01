.. _sponge_section:

******
Sponge
******

Castro uses a sponge source term to damp velocities in regions where
we are not interested in the dynamics, to prevent them from
dominating the timestep constraint.  Often these are buffer regions
between the domain of interest and the boundary conditions.

The sponge parameters are set in the inputs file. The timescale of the
damping is set through ``castro.sponge_timescale``, while factors such as
the radius/density/pressure at which the sponge starts to begin being applied
are described below.

The sponge value, :math:`f_\mathrm{sponge}` is computed as described below
and then the sponge factor, :math:`f` is computed as:

.. math::

   f = - \left [ 1 - \frac{1}{1 + \alpha f_\mathrm{sponge}}\right ]

for an implicit update or

.. math::

   f = -\alpha f_\mathrm{sponge}

.. index:: sponge_implicit, sponge_timescale

for an explicit update.  This choice is controlled by
``sponge_implicit``.  Here, :math:`\alpha` is constructed from the
``sponge_timescale`` parameter, :math:`t_\mathrm{sponge}` as:

.. math::

   \alpha = \frac{\Delta t}{t_\mathrm{sponge}}

The sponge source is then added to the momentum and total energy equations.

The general sponge parameters are:

       ==========================     ========================
         variable                       runtime parameter
       ==========================     ========================
       :math:`f_\mathrm{lower}`       ``sponge_lower_factor``
       :math:`f_\mathrm{upper}`       ``sponge_upper_factor``
       ==========================     ========================

There are three sponges, each controlled by different runtime parameters:

  * **radial sponge** : The radial sponge operates beyond some radius
    (define in terms of the ``center(:)`` position for the problem).
    It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & r < r_\mathrm{lower} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (r - r_\mathrm{lower})}{\Delta r} \right ) \right ]  & r_\mathrm{lower} \le r < r_\mathrm{upper} \\
                     f_\mathrm{upper} & r \ge r_\mathrm{upper}
             \end{array} \right .


    The parameters controlling the various quantities here are:

       ==========================     ========================
         variable                       runtime parameter
       ==========================     ========================
       :math:`r_\mathrm{lower}`       ``sponge_lower_radius``
       :math:`r_\mathrm{upper}`       ``sponge_upper_radius``
       ==========================     ========================

    and :math:`\Delta r = r_\mathrm{upper} - r_\mathrm{lower}` .


  * **density sponge** : The density sponge turns on based on density
    thresholds.  At high densities, we usually care about the
    dynamics, so we will want the sponge off, but as the density
    lowers, we gradually turn on the sponge.  It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & \rho > \rho_\mathrm{upper} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (\rho - \rho_\mathrm{upper})}{\Delta \rho} \right ) \right ]  & \rho_\mathrm{upper} \ge \rho > \rho_\mathrm{lower} \\
                     f_\mathrm{upper} & \rho < \rho_\mathrm{lower}
             \end{array} \right .


    The parameters controlling the various quantities here are:

       ============================     ==========================
         variable                          runtime parameter
       ============================     ==========================
       :math:`\rho_\mathrm{lower}`       ``sponge_lower_density``
       :math:`\rho_\mathrm{upper}`       ``sponge_upper_density``
       ============================     ==========================

    and :math:`\Delta \rho = \rho_\mathrm{upper} - \rho_\mathrm{lower}` .


  * **pressure sponge** : The pressure sponge is just like the density sponge,
    except it is keyed on pressure.  It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & p > p_\mathrm{upper} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (p - p_\mathrm{upper})}{\Delta p} \right ) \right ]  & p_\mathrm{upper} \ge p \ge p_\mathrm{lower} \\
                     f_\mathrm{upper} & \rho < \rho_\mathrm{lower}
             \end{array} \right .


    The parameters controlling the various quantities here are:

       ============================     ==========================
         variable                          runtime parameter
       ============================     ==========================
       :math:`p_\mathrm{lower}`         ``sponge_lower_pressure``
       :math:`p_\mathrm{upper}`         ``sponge_upper_pressure``
       ============================     ==========================

    and :math:`\Delta p = p_\mathrm{upper} - p_\mathrm{lower}` .

