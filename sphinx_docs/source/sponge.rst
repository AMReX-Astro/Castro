.. _sponge_section:

******
Sponge
******

Castro uses a sponge source term to damp velocities in regions where
we are not interested in the dynamics.  Often these are buffer regions
between the domain of interest and the boundary conditions.



   The purpose of the sponge is to damp velocities outside of a star, to
   prevent them from dominating the timestep constraint. The sponge parameters
   are set in your ``probin`` file, in the ``&sponge`` namelist. You can sponge either
   on radius from the center (using ``sponge_lower_radius`` and
   ``sponge_upper_radius``) or on density (using ``sponge_lower_density``
   and ``sponge_upper_density``). The timescale of the damping is
   set through ``sponge_timescale``.
