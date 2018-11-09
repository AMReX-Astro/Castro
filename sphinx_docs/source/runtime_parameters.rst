.. _ch:parameters:

 castro  Namespace
=================

.. raw:: latex

   \small

.. table:: castro : AMR
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we do the          | 1                     |
   |                       | hyperbolic reflux at  |                       |
   |    \endfoot           | coarse-fine           |                       |
   |                       | interfaces?           |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | how to do limiting of | 0                     |
   |                       | the state data when   |                       |
   |                       | interpolating 0: only |                       |
   |                       | prevent new extrema   |                       |
   |                       | 1: preserve linear    |                       |
   |                       | combinations of state |                       |
   |                       | variables             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | highest order used in | 1                     |
   |                       | interpolation         |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Number of ghost zones | 0                     |
   |                       | for state data to     |                       |
   |                       | have. Note that if    |                       |
   |                       | you are using         |                       |
   |                       | radiation, choosing   |                       |
   |                       | this to be zero will  |                       |
   |                       | be overridden since   |                       |
   |                       | radiation needs at    |                       |
   |                       | least one ghost zone. |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether to re-compute | 1                     |
   |                       | new-time source terms |                       |
   |    \rowcolor{tableSha | after a reflux        |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | should we have state  | 0                     |
   |                       | data for custom       |                       |
   |                       | load-balancing        |                       |
   |                       | weighting?            |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : diagnostics, I/O
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | abort if we exceed    | 1                     |
   |                       | CFL = 1 over the      |                       |
   |    \endfoot           | cource of a timestep  |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | a string describing   | ""                    |
   |                       | the simulation that   |                       |
   |                       | will be copied into   |                       |
   |                       | the plotfile’s        |                       |
   |                       | job_info file         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | write a final         | 1                     |
   |                       | plotfile and          |                       |
   |    \rowcolor{tableSha | checkpoint upon       |                       |
   | de}                   | completion            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | display warnings in   | (0, 1)                |
   |                       | Fortran90 routines    |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | display information   | (0, 1)                |
   |                       | about updates to the  |                       |
   |    \rowcolor{tableSha | state (how much mass, |                       |
   | de}                   | momentum, energy      |                       |
   |                       | added)                |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Do we want to reset   | -1                    |
   |                       | the number of steps   |                       |
   |                       | in the checkpoint?    |                       |
   |                       | This ONLY takes       |                       |
   |                       | effect if             |                       |
   |                       | amr.regrid_on_restart |                       |
   |                       | = 1 and               |                       |
   |                       | amr.checkpoint_on_res |                       |
   |                       | tart                  |                       |
   |                       | = 1, (which require   |                       |
   |                       | that max_step and     |                       |
   |                       | stop_time be less     |                       |
   |                       | than the value in the |                       |
   |                       | checkpoint) and you   |                       |
   |                       | set it to value       |                       |
   |                       | greater than this     |                       |
   |                       | default value.        |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do we want to reset   | -1.e200               |
   |                       | the time in the       |                       |
   |    \rowcolor{tableSha | checkpoint? This ONLY |                       |
   | de}                   | takes effect if       |                       |
   |                       | amr.regrid_on_restart |                       |
   |                       | = 1 and               |                       |
   |                       | amr.checkpoint_on_res |                       |
   |                       | tart                  |                       |
   |                       | = 1, (which require   |                       |
   |                       | that max_step and     |                       |
   |                       | stop_time be less     |                       |
   |                       | than the value in the |                       |
   |                       | checkpoint) and you   |                       |
   |                       | set it to value       |                       |
   |                       | greater than this     |                       |
   |                       | default value.        |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | display center of     | 0                     |
   |                       | mass diagnostics      |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | how often (number of  | -1                    |
   |                       | coarse timesteps) to  |                       |
   |    \rowcolor{tableSha | compute integral sums |                       |
   | de}                   | (for runtime          |                       |
   |                       | diagnostics)          |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | how often (simulation | -1.0e0                |
   |                       | time) to compute      |                       |
   |                       | integral sums (for    |                       |
   |                       | runtime diagnostics)  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | calculate losses of   | 0                     |
   |                       | material through      |                       |
   |    \rowcolor{tableSha | physical grid         |                       |
   | de}                   | boundaries            |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : diffusion
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | scaling factor for    | 1.0                   |
   |                       | conductivity          |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | set a cutoff density  | -1.e200               |
   |                       | for diffusion – we    |                       |
   |                       | zero the term out     |                       |
   |                       | below this density    |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | enable enthalpy       | 0                     |
   |                       | diffusion             |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | enable species        | 0                     |
   |                       | diffusion             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | enable thermal        | 0                     |
   |                       | diffusion             |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | enable velocity       | 0                     |
   |                       | diffusion             |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : embiggening
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the factor by which   | 1                     |
   |                       | to extend the domain  |                       |
   |    \endfoot           | upon restart for      |                       |
   |                       | embiggening           |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | used with the         | -1                    |
   |                       | embiggening routines  |                       |
   |                       | to determine how to   |                       |
   |                       | extend the domain     |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : gravity and rotation
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits gravity       | -1                    |
   |                       | calculation to be     |                       |
   |    \endfoot           | turned on and off     |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits rotation      | -1                    |
   |                       | calculation to be     |                       |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | determines how the    | 4                     |
   |                       | gravitational source  |                       |
   |    \rowcolor{tableSha | term is added to the  |                       |
   | de}                   | momentum and energy   |                       |
   |                       | state variables.      |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | we can do a implicit  | 1                     |
   |                       | solution of the       |                       |
   |                       | rotation update to    |                       |
   |                       | allow for better      |                       |
   |                       | coupling of the       |                       |
   |                       | Coriolis terms        |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | to we recompute the   | 0                     |
   |                       | center used for the   |                       |
   |    \rowcolor{tableSha | multipole gravity     |                       |
   | de}                   | solve each step?      |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | mass of the point     | 0.0                   |
   |                       | mass                  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we have a central  | 0                     |
   |                       | point mass, we can    |                       |
   |    \rowcolor{tableSha | prevent mass from     |                       |
   | de}                   | building up in the    |                       |
   |                       | zones adjacent to it  |                       |
   |                       | by keeping their      |                       |
   |                       | density constant and  |                       |
   |                       | adding their mass to  |                       |
   |                       | the point mass object |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the coordinate axis   | 3                     |
   |                       | (:math:`x=1`,         |                       |
   |                       | :math:`y=2`,          |                       |
   |                       | :math:`z=3`) for the  |                       |
   |                       | rotation vector       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | determines how the    | 4                     |
   |                       | rotation source terms |                       |
   |    \rowcolor{tableSha | are added to the      |                       |
   | de}                   | momentum and energy   |                       |
   |                       | equations             |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits the           | 1                     |
   |                       | centrifugal terms in  |                       |
   |                       | the rotation to be    |                       |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits the Coriolis  | 1                     |
   |                       | terms in the rotation |                       |
   |    \rowcolor{tableSha | to be turned on and   |                       |
   | de}                   | off                   |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits the           | 1                     |
   |                       | d(omega)/dt terms in  |                       |
   |                       | the rotation to be    |                       |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the rotation periods  | 0.0                   |
   |                       | time evolution—this   |                       |
   |    \rowcolor{tableSha | allows the rotation   |                       |
   | de}                   | rate to change        |                       |
   |                       | durning the           |                       |
   |                       | simulation time       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the rotation period   | -1.e200               |
   |                       | for the corotating    |                       |
   |                       | frame                 |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Which reference frame | 1                     |
   |                       | to measure the state  |                       |
   |    \rowcolor{tableSha | variables with        |                       |
   | de}                   | respect to. The       |                       |
   |                       | standard in the       |                       |
   |                       | literature when using |                       |
   |                       | a rotating reference  |                       |
   |                       | frame is to measure   |                       |
   |                       | the state variables   |                       |
   |                       | with respect to an    |                       |
   |                       | observer fixed in     |                       |
   |                       | that rotating frame.  |                       |
   |                       | If this option is     |                       |
   |                       | disabled by setting   |                       |
   |                       | it to 0, the state    |                       |
   |                       | variables will be     |                       |
   |                       | measured with respect |                       |
   |                       | to an observer fixed  |                       |
   |                       | in the inertial frame |                       |
   |                       | (but the frame will   |                       |
   |                       | still rotate).        |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | include a central     | 1                     |
   |                       | point mass            |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : hydrodynamics
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if true, define an    | 0                     |
   |                       | additional source     |                       |
   |    \endfoot           | term                  |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Whether or not to     | 1                     |
   |                       | allow the internal    |                       |
   |                       | energy to be less     |                       |
   |                       | than the internal     |                       |
   |                       | energy corresponding  |                       |
   |                       | to small_temp         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | for the Colella &     | 2                     |
   |                       | Glaz Riemann solver,  |                       |
   |    \rowcolor{tableSha | what to do if we do   |                       |
   | de}                   | not converge to a     |                       |
   |                       | solution for the star |                       |
   |                       | state. 0 = do         |                       |
   |                       | nothing; print        |                       |
   |                       | iterations and exit 1 |                       |
   |                       | = revert to the       |                       |
   |                       | original guess for    |                       |
   |                       | p-star 2 = do a       |                       |
   |                       | bisection search for  |                       |
   |                       | another 2 \*          |                       |
   |                       | cg_maxiter            |                       |
   |                       | iterations.           |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | for the Colella &     | 12                    |
   |                       | Glaz Riemann solver,  |                       |
   |                       | the maximum number of |                       |
   |                       | iterations to take    |                       |
   |                       | when solving for the  |                       |
   |                       | star state            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | for the Colella &     | 1.0e-5                |
   |                       | Glaz Riemann solver,  |                       |
   |    \rowcolor{tableSha | the tolerance to      |                       |
   | de}                   | demand in finding the |                       |
   |                       | star state            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Which method to use   | 1                     |
   |                       | when resetting a      |                       |
   |                       | negative/small        |                       |
   |                       | density 1 = Reset to  |                       |
   |                       | characteristics of    |                       |
   |                       | adjacent zone with    |                       |
   |                       | largest density 2 =   |                       |
   |                       | Use average of all    |                       |
   |                       | adjacent zones for    |                       |
   |                       | all state variables 3 |                       |
   |                       | = Reset to the        |                       |
   |                       | original zone state   |                       |
   |                       | before the hydro      |                       |
   |                       | update                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the coefficient of    | 0.1                   |
   |                       | the artificial        |                       |
   |    \rowcolor{tableSha | viscosity             |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we do the CTU      | 1                     |
   |                       | unsplit method or a   |                       |
   |                       | method-of-lines       |                       |
   |                       | approach?             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits hydro to be   | -1                    |
   |                       | turned on and off for |                       |
   |    \rowcolor{tableSha | running pure rad      |                       |
   | de}                   | problems              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits sponge to be  | 0                     |
   |                       | turned on and off     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Threshold value of (E | 1.0e0                 |
   |                       | - K) / E such that    |                       |
   |    \rowcolor{tableSha | above eta1, the       |                       |
   | de}                   | hydrodynamic pressure |                       |
   |                       | is derived from E -   |                       |
   |                       | K; otherwise, we use  |                       |
   |                       | the internal energy   |                       |
   |                       | variable UEINT.       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Threshold value of (E | 1.0e-4                |
   |                       | - K) / E such that    |                       |
   |                       | above eta2, we update |                       |
   |                       | the internal energy   |                       |
   |                       | variable UEINT to     |                       |
   |                       | match E - K. Below    |                       |
   |                       | this, UEINT remains   |                       |
   |                       | unchanged.            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | set the flattening    | 0                     |
   |                       | parameter to zero to  |                       |
   |    \rowcolor{tableSha | force the             |                       |
   | de}                   | reconstructed         |                       |
   |                       | profiles to be flat,  |                       |
   |                       | resulting in a        |                       |
   |                       | first-order method    |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       | 0                     |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we do fourth-order | 0                     |
   |                       | accurate MOL hydro?   |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing HSE   | 0                     |
   |                       | boundary conditions,  |                       |
   |                       | should we get the     |                       |
   |                       | temperature via       |                       |
   |                       | interpolation (using  |                       |
   |                       | model_parser) or hold |                       |
   |                       | it constant?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing HSE   | 0                     |
   |                       | boundary conditions,  |                       |
   |    \rowcolor{tableSha | how do we treat the   |                       |
   | de}                   | velocity? reflect? or |                       |
   |                       | outflow?              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing HSE   | 0                     |
   |                       | boundary conditions,  |                       |
   |                       | do we zero the        |                       |
   |                       | velocity?             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether to use the    | 0                     |
   |                       | hybrid advection      |                       |
   |    \rowcolor{tableSha | scheme that updates   |                       |
   | de}                   | z-angular momentum,   |                       |
   |                       | cylindrical momentum, |                       |
   |                       | and azimuthal         |                       |
   |                       | momentum (3D only)    |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we drop from our   | 0                     |
   |                       | regular Riemann       |                       |
   |                       | solver to HLL when we |                       |
   |                       | are in shocks to      |                       |
   |                       | avoid the odd-even    |                       |
   |                       | decoupling            |                       |
   |                       | instability?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Should we limit the   | 0                     |
   |                       | density fluxes so     |                       |
   |    \rowcolor{tableSha | that we do not create |                       |
   | de}                   | small densities?      |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | integration order for | 2                     |
   |                       | MOL integration 1 =   |                       |
   |                       | first order, 2 =      |                       |
   |                       | second order TVD, 3 = |                       |
   |                       | 3rd order TVD, 4 =    |                       |
   |                       | 4th order RK          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | for piecewise linear, | 2                     |
   |                       | reconstruction order  |                       |
   |    \rowcolor{tableSha | to use                |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we construct       | 0                     |
   |                       | :math:`\gamma_e = p/( |                       |
   |                       | \rho e) + 1`          |                       |
   |                       | and bring it to the   |                       |
   |                       | interfaces for        |                       |
   |                       | additional            |                       |
   |                       | thermodynamic         |                       |
   |                       | information (this is  |                       |
   |                       | the Colella & Glaz    |                       |
   |                       | technique) or do we   |                       |
   |                       | use :math:`(\rho e)`  |                       |
   |                       | (the classic          |                       |
   |                       | Castro behavior).     |                       |
   |                       | Note this also uses   |                       |
   |                       | :math:`\tau = 1/\rho` |                       |
   |                       | instead of            |                       |
   |                       | :math:`\rho`.         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we use the         | 0                     |
   |                       | reference state in    |                       |
   |    \rowcolor{tableSha | evaluating the        |                       |
   | de}                   | eigenvectors?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | various methods of    | 0                     |
   |                       | giving temperature a  |                       |
   |                       | larger role in the    |                       |
   |                       | reconstruction—see    |                       |
   |                       | Zingale & Katz 2015   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | reconstruction type:  | 1                     |
   |                       | 0: piecewise linear;  |                       |
   |    \rowcolor{tableSha | 1: classic Colella &  |                       |
   | de}                   | Woodward ppm; 2:      |                       |
   |                       | extrema-preserving    |                       |
   |                       | ppm                   |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | which Riemann solver  | 0                     |
   |                       | do we use: 0:         |                       |
   |                       | Colella, Glaz, &      |                       |
   |                       | Ferguson (a two-shock |                       |
   |                       | solver); 1: Colella & |                       |
   |                       | Glaz (a two-shock     |                       |
   |                       | solver) 2: HLLC       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the small density     | -1.e200               |
   |                       | cutoff. Densities     |                       |
   |    \rowcolor{tableSha | below this value will |                       |
   | de}                   | be reset              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the small specific    | -1.e200               |
   |                       | internal energy       |                       |
   |                       | cutoff. Internal      |                       |
   |                       | energies below this   |                       |
   |                       | value will be reset   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the small pressure    | -1.e200               |
   |                       | cutoff. Pressures     |                       |
   |    \rowcolor{tableSha | below this value will |                       |
   | de}                   | be reset              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the small temperature | -1.e200               |
   |                       | cutoff. Temperatures  |                       |
   |                       | below this value will |                       |
   |                       | be reset              |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | extrapolate the       | 0                     |
   |                       | source terms (gravity |                       |
   |    \rowcolor{tableSha | and rotation) to      |                       |
   | de}                   | :math:`n+1/2`         |                       |
   |                       | timelevel for use in  |                       |
   |                       | the interface state   |                       |
   |                       | prediction            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are using the   | 1                     |
   |                       | sponge, whether to    |                       |
   |                       | use the implicit      |                       |
   |                       | solve for it          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if the transverse     | 1                     |
   |                       | interface state       |                       |
   |    \rowcolor{tableSha | correction, if the    |                       |
   | de}                   | new density is        |                       |
   |                       | negative, then        |                       |
   |                       | replace all of the    |                       |
   |                       | interface quantities  |                       |
   |                       | with their values     |                       |
   |                       | without the           |                       |
   |                       | transverse            |                       |
   |                       | correction.           |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if the interface      | 0                     |
   |                       | state for             |                       |
   |                       | :math:`(\rho e)` is   |                       |
   |                       | negative after we add |                       |
   |                       | the transverse terms, |                       |
   |                       | then replace the      |                       |
   |                       | interface value of    |                       |
   |                       | :math:`(\rho e)` with |                       |
   |                       | a value constructed   |                       |
   |                       | from the              |                       |
   |                       | :math:`(\rho e)`      |                       |
   |                       | evolution equation    |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | after we add the      | 0                     |
   |                       | transverse correction |                       |
   |    \rowcolor{tableSha | to the interface      |                       |
   | de}                   | states, replace the   |                       |
   |                       | predicted pressure    |                       |
   |                       | with an EOS call      |                       |
   |                       | (using :math:`e` and  |                       |
   |                       | :math:`\rho`).        |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | should we use the EOS | 0                     |
   |                       | in the Riemann solver |                       |
   |                       | to ensure             |                       |
   |                       | thermodynamic         |                       |
   |                       | consistency?          |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | flatten the           | 1                     |
   |                       | reconstructed         |                       |
   |    \rowcolor{tableSha | profiles around       |                       |
   | de}                   | shocks to prevent     |                       |
   |                       | them from becoming    |                       |
   |                       | too thin              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | for the piecewise     | 1                     |
   |                       | linear                |                       |
   |                       | reconstruction, do we |                       |
   |                       | subtract off          |                       |
   |                       | :math:`(\rho g)` from |                       |
   |                       | the pressure before   |                       |
   |                       | limiting?             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing an    | ""                    |
   |                       | external -x boundary  |                       |
   |    \rowcolor{tableSha | condition, who do we  |                       |
   | de}                   | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing an    | ""                    |
   |                       | external +x boundary  |                       |
   |                       | condition, who do we  |                       |
   |                       | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing an    | ""                    |
   |                       | external -y boundary  |                       |
   |    \rowcolor{tableSha | condition, who do we  |                       |
   | de}                   | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing an    | ""                    |
   |                       | external +y boundary  |                       |
   |                       | condition, who do we  |                       |
   |                       | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if we are doing an    | ""                    |
   |                       | external -z boundary  |                       |
   |    \rowcolor{tableSha | condition, who do we  |                       |
   | de}                   | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | if we are doing an    | ""                    |
   |                       | external +z boundary  |                       |
   |                       | condition, who do we  |                       |
   |                       | interpret it?         |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : parallelization
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        |                       | 1                     |
   |                       |                       |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | determines whether we | -1                    |
   |                       | use accelerators for  |                       |
   |                       | specific loops        |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : particles
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | permits tracer        | 0                     |
   |                       | particle calculation  |                       |
   |    \endfoot           | to be turned on and   |                       |
   |                       | off                   |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : reactions
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | disable burning       | 0                     |
   |                       | inside hydrodynamic   |                       |
   |    \endfoot           | shock regions         |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | permits reactions to  | -1                    |
   |                       | be turned on and off  |                       |
   |                       | – mostly for          |                       |
   |                       | efficiency’s sake     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Limit the timestep    | 1.e200                |
   |                       | based on how much the |                       |
   |    \rowcolor{tableSha | burning can change    |                       |
   | de}                   | the species mass      |                       |
   |                       | fractions of a zone.  |                       |
   |                       | The timestep is equal |                       |
   |                       | to dtnuc              |                       |
   |                       | :math:`\cdot\,(X / \d |                       |
   |                       | ot{X})`.              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | If we are using the   | 1.e-3                 |
   |                       | timestep limiter      |                       |
   |                       | based on changes in   |                       |
   |                       | :math:`X`, set a      |                       |
   |                       | threshold on the      |                       |
   |                       | species abundance     |                       |
   |                       | below which the       |                       |
   |                       | limiter is not        |                       |
   |                       | applied. This helps   |                       |
   |                       | prevent the timestep  |                       |
   |                       | from becoming very    |                       |
   |                       | small due to changes  |                       |
   |                       | in trace species.     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Limit the timestep    | 1.e200                |
   |                       | based on how much the |                       |
   |    \rowcolor{tableSha | burning can change    |                       |
   | de}                   | the internal energy   |                       |
   |                       | of a zone. The        |                       |
   |                       | timestep is equal to  |                       |
   |                       | dtnuc                 |                       |
   |                       | :math:`\cdot\,(e / \d |                       |
   |                       | ot{e})`.              |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | limit the zone size   | 1.e200                |
   |                       | based on how much the |                       |
   |                       | burning can change    |                       |
   |                       | the internal energy   |                       |
   |                       | of a zone. The zone   |                       |
   |                       | size on the finest    |                       |
   |                       | level must be smaller |                       |
   |                       | than dxnuc            |                       |
   |                       | :math:`\cdot\, c_s\cd |                       |
   |                       | ot (e / \dot{e})`,    |                       |
   |                       | where :math:`c_s` is  |                       |
   |                       | the sound speed. This |                       |
   |                       | ensures that the      |                       |
   |                       | sound-crossing time   |                       |
   |                       | is smaller than the   |                       |
   |                       | nuclear energy        |                       |
   |                       | injection timescale.  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Disable limiting      | 1.e200                |
   |                       | based on dxnuc above  |                       |
   |    \rowcolor{tableSha | this threshold. This  |                       |
   | de}                   | allows zones that     |                       |
   |                       | have already ignited  |                       |
   |                       | or are about to       |                       |
   |                       | ignite to be          |                       |
   |                       | de-refined.           |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Disable limiting      | -1                    |
   |                       | based on dxnuc above  |                       |
   |                       | this AMR level.       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | maximum temperature   | 1.e200                |
   |                       | for allowing          |                       |
   |    \rowcolor{tableSha | reactions to occur in |                       |
   | de}                   | a zone                |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | minimum temperature   | 0.0                   |
   |                       | for allowing          |                       |
   |                       | reactions to occur in |                       |
   |                       | a zone                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | maximum density for   | 1.e200                |
   |                       | allowing reactions to |                       |
   |    \rowcolor{tableSha | occur in a zone       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | minimum density for   | 0.0                   |
   |                       | allowing reactions to |                       |
   |                       | occur in a zone       |                       |
   +-----------------------+-----------------------+-----------------------+

.. raw:: latex

   \small

.. table:: castro : refinement
parameters

   +--------------------------+--+---+
   |                          |  |   |
   +--------------------------+--+---+
   | Table —continued         |  |   |
   +--------------------------+--+---+
   |                          |  |   |
   +--------------------------+--+---+
   |                          |  |   |
   +--------------------------+--+---+
   | .. raw:: latex           |  | 0 |
   |                          |  |   |
   |    \endfoot              |  |   |
   |                          |  |   |
   | .. raw:: latex           |  |   |
   |                          |  |   |
   |    \hline                |  |   |
   |                          |  |   |
   | .. raw:: latex           |  |   |
   |                          |  |   |
   |    \endlastfoot          |  |   |
   |                          |  |   |
   | .. raw:: latex           |  |   |
   |                          |  |   |
   |    \rowcolor{tableShade} |  |   |
   +--------------------------+--+---+
   |                          |  | 0 |
   +--------------------------+--+---+

.. raw:: latex

   \small

.. table:: castro : timestep control
parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the effective Courant | 0.8                   |
   |                       | number to use—we will |                       |
   |    \endfoot           | not allow the         |                       |
   |                       | hydrodynamic waves to |                       |
   | .. raw:: latex        | cross more than this  |                       |
   |                       | fraction of a zone    |                       |
   |    \hline             | over a single         |                       |
   |                       | timestep              |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the maximum factor by | 1.1                   |
   |                       | which the timestep    |                       |
   |                       | can increase from one |                       |
   |                       | step to the next.     |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | If we do request more | 1                     |
   |                       | than the maximum      |                       |
   |    \rowcolor{tableSha | number of subcycles,  |                       |
   | de}                   | should we fail, or    |                       |
   |                       | should we clamp to    |                       |
   |                       | that maximum number   |                       |
   |                       | and perform that      |                       |
   |                       | many?                 |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the smallest valid    | 0.0                   |
   |                       | timestep—if we go     |                       |
   |                       | below this, we abort  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | a fixed timestep to   | -1.0                  |
   |                       | use for all steps     |                       |
   |    \rowcolor{tableSha | (negative turns it    |                       |
   | de}                   | off)                  |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | a factor by which to  | 1.0                   |
   |                       | reduce the first      |                       |
   |                       | timestep from that    |                       |
   |                       | requested by the      |                       |
   |                       | timestep estimators   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the initial timestep  | -1.0                  |
   |                       | (negative uses the    |                       |
   |    \rowcolor{tableSha | step returned from    |                       |
   | de}                   | the timestep          |                       |
   |                       | constraints)          |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the largest valid     | 1.e200                |
   |                       | timestep—limit all    |                       |
   |                       | timesteps to be no    |                       |
   |                       | larger than this      |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do not permit more    | 10                    |
   |                       | subcycled timesteps   |                       |
   |    \rowcolor{tableSha | than this parameter.  |                       |
   | de}                   | Set to a negative     |                       |
   |                       | value to disable this |                       |
   |                       | criterion.            |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | enforce that the AMR  | 0                     |
   |                       | plot interval must be |                       |
   |                       | hit exactly           |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | If we’re doing        | 1.e-1                 |
   |                       | retries, set the      |                       |
   |    \rowcolor{tableSha | target threshold for  |                       |
   | de}                   | changes in density if |                       |
   |                       | a retry is triggered  |                       |
   |                       | by a negative         |                       |
   |                       | density. If this is   |                       |
   |                       | set to a negative     |                       |
   |                       | number then it will   |                       |
   |                       | disable retries using |                       |
   |                       | this criterion.       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | When performing a     | 0.5                   |
   |                       | retry, the factor to  |                       |
   |                       | multiply the current  |                       |
   |                       | timestep by when      |                       |
   |                       | trying again.         |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Tolerance to use when | 0.02                  |
   |                       | evaluating whether to |                       |
   |    \rowcolor{tableSha | do a retry. The       |                       |
   | de}                   | timestep suggested by |                       |
   |                       | the retry will be     |                       |
   |                       | multiplied by (1 +    |                       |
   |                       | this factor) before   |                       |
   |                       | comparing the actual  |                       |
   |                       | timestep to it. If    |                       |
   |                       | set to some number    |                       |
   |                       | slightly larger than  |                       |
   |                       | zero, then this       |                       |
   |                       | prevents retries that |                       |
   |                       | are caused by small   |                       |
   |                       | numerical             |                       |
   |                       | differences.          |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Number of iterations  | 2                     |
   |                       | for the SDC advance.  |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | enforce that the AMR  | 0                     |
   |                       | small plot interval   |                       |
   |    \rowcolor{tableSha | must be hit exactly   |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Check for a possible  | 0                     |
   |                       | post-timestep regrid  |                       |
   |                       | if certain stability  |                       |
   |                       | criteria were         |                       |
   |                       | violated.             |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Retry a timestep if   | 0                     |
   |                       | it violated the       |                       |
   |    \rowcolor{tableSha | timestep-limiting     |                       |
   | de}                   | criteria over the     |                       |
   |                       | course of an advance. |                       |
   |                       | The criteria will     |                       |
   |                       | suggest a new         |                       |
   |                       | timestep that         |                       |
   |                       | satisfies the         |                       |
   |                       | criteria, and we will |                       |
   |                       | do subcycled          |                       |
   |                       | timesteps on the same |                       |
   |                       | level until we reach  |                       |
   |                       | the original target   |                       |
   |                       | time.                 |                       |
   +-----------------------+-----------------------+-----------------------+

.. _ch:parameters:

 diffusion  Namespace
====================

.. raw:: latex

   \small

.. table:: diffusion parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Use MLMG as the       | 4                     |
   |                       | operator              |                       |
   |    \endfoot           |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the level of          | 0                     |
   |                       | verbosity for the     |                       |
   |                       | diffusion solve       |                       |
   |                       | (higher number means  |                       |
   |                       | more output)          |                       |
   +-----------------------+-----------------------+-----------------------+

.. _ch:parameters:

 gravity  Namespace
==================

.. raw:: latex

   \small

.. table:: gravity parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | if doing constant     | 0.0                   |
   |                       | gravity, what is the  |                       |
   |    \endfoot           | acceleration          |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Check if the user     | 0                     |
   |                       | wants to compute the  |                       |
   |                       | boundary conditions   |                       |
   |                       | using the brute force |                       |
   |                       | method. Default is    |                       |
   |                       | false, since this     |                       |
   |                       | method is slow.       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | should we apply a     | 1                     |
   |                       | lagged correction to  |                       |
   |    \rowcolor{tableSha | the potential that    |                       |
   | de}                   | gets us closer to the |                       |
   |                       | composite solution?   |                       |
   |                       | This makes the        |                       |
   |                       | resulting fine grid   |                       |
   |                       | calculation slightly  |                       |
   |                       | more accurate, at the |                       |
   |                       | cost of an additional |                       |
   |                       | Poisson solve per     |                       |
   |                       | timestep.             |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | ratio of dr for       | 1                     |
   |                       | monopole gravity      |                       |
   |                       | binning to grid       |                       |
   |                       | resolution            |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | For non-Poisson       | 0                     |
   |                       | gravity, do we want   |                       |
   |    \rowcolor{tableSha | to construct the      |                       |
   | de}                   | gravitational         |                       |
   |                       | acceleration by       |                       |
   |                       | taking the gradient   |                       |
   |                       | of the potential,     |                       |
   |                       | rather than           |                       |
   |                       | constructing it       |                       |
   |                       | directly?             |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | what type             | "fillme"              |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the maximum mulitpole | 0                     |
   |                       | order to use for      |                       |
   |    \rowcolor{tableSha | multipole BCs when    |                       |
   | de}                   | doing Poisson gravity |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | For all gravity       | MAX_LEV-1             |
   |                       | types, we can choose  |                       |
   |                       | a maximum level for   |                       |
   |                       | explicitly            |                       |
   |                       | calculating the       |                       |
   |                       | gravity and           |                       |
   |                       | associated potential. |                       |
   |                       | Above that level, we  |                       |
   |                       | interpolate from      |                       |
   |                       | coarser levels.       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | Do agglomeration?     | 1                     |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       | 1                     |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | how many FMG cycles?  | 0                     |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | Do N-Solve?           | 0                     |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | do we do a composite  | 0                     |
   |                       | solve?                |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | do we perform the     | 0                     |
   |                       | synchronization at    |                       |
   |                       | coarse-fine           |                       |
   |                       | interfaces?           |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the level of          | 0                     |
   |                       | verbosity for the     |                       |
   |    \rowcolor{tableSha | gravity solve (higher |                       |
   | de}                   | number means more     |                       |
   |                       | output on the status  |                       |
   |                       | of the solve /        |                       |
   |                       | multigrid             |                       |
   +-----------------------+-----------------------+-----------------------+

.. _ch:parameters:

 particles  Namespace
====================

.. raw:: latex

   \small

.. table:: particles parameters

   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | Table —continued      |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the name of an input  | ""                    |
   |                       | file containing the   |                       |
   |    \endfoot           | total particle number |                       |
   |                       | and the initial       |                       |
   | .. raw:: latex        | position of each      |                       |
   |                       | particle.             |                       |
   |    \hline             |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \endlastfoot       |                       |                       |
   |                       |                       |                       |
   | .. raw:: latex        |                       |                       |
   |                       |                       |                       |
   |    \rowcolor{tableSha |                       |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the name of timestamp | ""                    |
   |                       | files.                |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | the name of a file    | ""                    |
   |                       | with new particles at |                       |
   |    \rowcolor{tableSha | restart               |                       |
   | de}                   |                       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | to restart from a     | 0                     |
   |                       | checkpoint that was   |                       |
   |                       | written with          |                       |
   |                       | USE_PARTICLES=FALSE   |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether the local     | 1                     |
   |                       | densities at given    |                       |
   |    \rowcolor{tableSha | positions of          |                       |
   | de}                   | particles are stored  |                       |
   |                       | in output files       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the name of a         | ""                    |
   |                       | directory in which    |                       |
   |                       | timestamp files are   |                       |
   |                       | stored.               |                       |
   +-----------------------+-----------------------+-----------------------+
   | .. raw:: latex        | whether the local     | 0                     |
   |                       | temperatures at given |                       |
   |    \rowcolor{tableSha | positions of          |                       |
   | de}                   | particles are stored  |                       |
   |                       | in output files       |                       |
   +-----------------------+-----------------------+-----------------------+
   |                       | the level of          | 0                     |
   |                       | verbosity for the     |                       |
   |                       | tracer particle (0 or |                       |
   |                       | 1)                    |                       |
   +-----------------------+-----------------------+-----------------------+
