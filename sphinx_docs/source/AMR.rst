************************
Adaptive Mesh Refinement
************************

Our approach to adaptive refinement in Castro uses a nested hierarchy
of logically-rectangular grids with simultaneous refinement of the
grids in both space and time. The integration algorithm on the grid
hierarchy is a recursive procedure in which coarse grids are advanced
in time, fine grids are advanced multiple steps to reach the same time
as the coarse grids and the data at different levels are then
synchronized.

During the regridding step, increasingly finer grids
are recursively embedded in coarse grids until the solution is
sufficiently resolved. An error estimation procedure based on
user-specified criteria (described in Section `1 <#sec:tagging>`__)
evaluates where additional refinement is needed
and grid generation procedures dynamically create or
remove rectangular fine grid patches as resolution requirements change.

A good introduction to the style of AMR used here is in Lecture 1
of the Adaptive Mesh Refinement Short Course at
https://ccse.lbl.gov/people/jbb/shortcourse/lecture1.pdf.

.. _sec:tagging:

Tagging for Refinement
======================

Castro determines what zones should be tagged for refinement at the
next regridding step by using a set of built-in routines that test on
quantities such as the density and pressure and determining whether
the quantities themselves or their gradients pass a user-specified
threshold. This may then be extended if amr.n_error_buf :math:`> 0`
to a certain number of zones beyond these tagged zones. This section
describes the process by which zones are tagged, and describes how to
add customized tagging criteria.

The routines for tagging cells are located in the
Tagging_nd.f90 file in the Source/driver/ directory. (These are
dimension-agnostic routines that loop over all three dimensional
indices even for 1D or 2D problems.) The main routines are
ca_denerror, ca_temperror, ca_presserror,
ca_velerror, and ca_raderror. They refine based on
density, temperature, pressure, velocity, and radiation energy density
(if enabled), respectively. The same approach is used for all of
them. As an example, we consider the density tagging routine. There
are four parameters that control tagging. If the density in a zone is
greater than the user-specified parameter denerr, then that
zone will be tagged for refinement, but only if the current AMR level
is less than the user-specified parameter max_denerr_lev.
Similarly, if the absolute density gradient between a zone and any
adjacent zone is greater than the user-specified parameter
dengrad, that zone will be tagged for refinement, but only
if we are currently on a level below
max_dengrad_lev. Note that setting denerr alone
will not do anything; you’ll need to set max_dengrad_lev :math:`>=
1` for this to have any effect.

All four of these parameters are set in the &tagging namelist
in your probin file. If left unmodified, they
default to a value that means we will never tag. The complete set of
parameters that can be controlled this way is the following:

-  density:

   -  value: denerr, max_denerr_lev

   -  gradient: dengrad, max_dengrad_lev

-  temperature:

   -  value: temperr, max_temperr_lev

   -  gradient: tempgrad, max_tempgrad_lev

-  velocity (magnitude):

   -  value: velerr, max_velerr_lev

   -  gradient: velgrad, max_velgrad_lev

-  pressure:

   -  value: presserr, max_presserr_lev

   -  gradient: pressgrad, max_pressgrad_lev

-  radiation energy density:

   -  value: raderr, max_raderr_lev

   -  gradient: radgrad, max_radgrad_lev

Since there are multiple algorithms for determining
whether a zone is tagged or not, it is worthwhile to specify
in detail what is happening to a zone in the code during this step.
We show this in the following pseudocode section. A zone
is tagged if the variable itag = SET, and is not tagged
if itag = CLEAR (these are mapped to 1 and 0, respectively).

::

    itag = CLEAR

    for errfunc[k] from k = 1 ... N
        // Three possibilities for itag: SET or CLEAR or remaining unchanged
        call errfunc[k](itag)
    end for

In particular, notice that there is an order dependence of this operation; if errfunc[2]
CLEARs a zone and then errfunc[3] SETs that zone, the final operation will
be to tag that zone (and vice versa). In practice by default this does not matter, because the
built-in tagging routines never explicitly perform a ``CLEAR``. However,
it is possible to overwrite the Tagging_nd.f90 file if you want to change how
ca_denerror, ca_temperror, etc. operate. This is not recommended, and if you do so
be aware that CLEARing a zone this way may not have the desired effect.

We provide also the ability for the user to define their own tagging criteria.
This is done through the Fortran function set_problem_tags in the
problem_tagging_d.f90 files. This function is provided the entire
state (including density, temperature, velocity, etc.) and the array
of tagging status for every zone. As an example of how to use this, suppose we
have a 3D Cartesian simulation where we want to tag any zone that has a
density gradient greater than 10, but we don’t care about any regions
outside a radius :math:`r > 75` from the problem origin; we leave them always unrefined.
We also want to ensure that the region :math:`r \leq 10` is always refined.
In our probin file we would set denerr = 10 and max_denerr_lev = 1
in the &tagging namelist. We would also make a copy of
problem_tagging_3d.f90 to our work directory and set it up as follows:

::

    subroutine set_problem_tags(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                state,state_l1,state_l2,state_l3, &
                                state_h1,state_h2,state_h3,&
                                set,clear,&
                                lo,hi,&
                                dx,problo,time,level)

      use bl_constants_module, only: ZERO, HALF
      use prob_params_module, only: center
      use meth_params_module, only: URHO, UMX, UMY, UMZ, UEDEN, NVAR

      implicit none

      integer         ,intent(in   ) :: lo(3),hi(3)
      integer         ,intent(in   ) :: state_l1,state_l2,state_l3, &
                                        state_h1,state_h2,state_h3
      integer         ,intent(in   ) :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      double precision,intent(in   ) :: state(state_l1:state_h1, &
                                              state_l2:state_h2, &
                                              state_l3:state_h3,NVAR)
      integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision,intent(in   ) :: problo(3),dx(3),time
      integer         ,intent(in   ) :: level,set,clear

      double precision :: x, y, z, r

      do k = lo(3), hi(3)
         z = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
         do j = lo(2), hi(2)
            y = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
            do i = lo(1), hi(1)
               x = problo(1) + (dble(i) + HALF) * dx(1) - center(2)

               r = (x**2 + y**2 + z**2)**(HALF)

               if (r > 75.0) then
                 tag(i,j,k) = clear
               elseif (r <= 10.0) then
                 tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

    end subroutine set_problem_tags

.. _sec:amr_synchronization:

Synchronization Algorithm
=========================

Here we present the AMR algorithm for the compressible equations with
self-gravity. The gravity component of the algorithm is closely
related to (but not identical to) that in Miniati and Colella, JCP,
2007. The content here is largely based on the content in the original
Castro paper (:cite:`castro_I`). The most significant difference is the
addition of a different strategy for when to employ the synchronization;
but regardless of whether the original or new strategy is used, the fundamental
synchronization step is identical.

.. _sec:synchronization_methodology:

Synchronization Methodology
---------------------------

Over a coarse grid time step we collect flux register information for
the hyperbolic part of the synchronization:

.. math:: \delta\Fb = -\Delta t_c A^c F^c + \sum \Delta t_f A^f F^f

Analogously, at the end of a coarse grid time step we store the
mismatch in normal gradients of :math:`\phi` at the coarse-fine interface:

.. math::

   \delta F_\phi =  - A^c \frac{\partial \phi^c}{\partial n}
   + \sum A^f \frac{\partial \phi^f}{\partial n}

We want the composite :math:`\phi^{c-f}` to satisfy the multilevel
version of (`[eq:Self Gravity] <#eq:Self Gravity>`__) at the synchronization time, just
as we want the coarse and fine fluxes at that time to match. So the goal
is to synchronize :math:`\phi` across levels at that time and then zero out
this mismatch register.

At the end of a coarse grid time step we can define
:math:`{\overline{\Ub}}^{c-f}` and :math:`\overline{\phi}^{c-f}` as the composite
of the data from coarse and fine grids as a provisional solution at
time :math:`n+1`. (Assume :math:`\overline{\Ub}` has been averaged down so that
the data on coarse cells underlying fine cells is the average of the
fine cell data above it.)

The synchronization consists of two parts:

-  Step 1: Hyperbolic reflux

   In the hyperbolic reflux step, we update the conserved variables with
   the flux synchronization and adjust the gravitational terms to reflect
   the changes in :math:`\rho` and :math:`\ub`.

   .. math:: {\Ub}^{c, \star} = \overline{\Ub}^{c} + \frac{\delta\Fb}{V},

   where :math:`V` is the volume of the cell and the correction from
   :math:`\delta\Fb` is supported only on coarse cells adjacent to fine grids.

   Note: this can be enabled/disabled via castro.do_reflux. Generally,
   it should be enabled (1).

   Also note that for axisymmetric or 1D spherical coordinates, the
   reflux of the pressure gradient is different, since it cannot be
   expressed as a divergence in those geometries. We use a separate
   flux register in the hydro code to store the pressure term in these
   cases.

-  Step 2: Gravitational synchronization

   In this step we correct for the mismatch in normal derivative in
   :math:`\phi^{c-f}` at the coarse-fine interface, as well as accounting for
   the changes in source terms for :math:`(\rho \ub)` and :math:`(\rho E)` due to the
   change in :math:`\rho.`

   On the coarse grid only, we define

   .. math:: (\delta \rho)^{c} =  \rho^{c, \star} - {\overline{\rho}}^{c}  .

   We then form the composite residual, which is composed of two
   contributions. The first is the degree to which the current :math:`\overline{\phi}^{c-f}` does not satisfy the original equation on a
   composite grid (since we have solved for :math:`\overline{\phi}^{c-f}`
   separately on the coarse and fine levels). The second is the response
   of :math:`\phi` to the change in :math:`\rho.` We define

   .. math::

      R \equiv  4 \pi G \rho^{\star,c-f} - \Delta^{c-f} \; \overline{\phi}^{c-f}
      = - 4 \pi G (\delta \rho)^c - (\nabla \cdot \delta F_\phi ) |_c   .

   Then we solve

   .. math::

      \Delta^{c-f} \; \delta \phi^{c-f} = R
      \label{eq:gravsync}

   as a two level solve at the coarse and fine levels.
   We define the update to gravity,

   .. math:: \delta {\bf g}^{c-f} = \nabla (\delta \phi^{c-f})  .

   Finally, we need to

   -  add :math:`\delta \phi^{c-f}` directly to
      to :math:`\phi^{c}` and :math:`\phi^{f}` and interpolate :math:`\delta \phi^{c-f}` to any finer
      levels and add to the current :math:`\phi` at those levels.

   -  if level :math:`c` is not the coarsest level in the calculation, then we must transmit the
      effect of this change in :math:`\phi` to the coarser levels by updating the flux register between
      level :math:`c` and the next coarser level, :math:`cc.` In particular, we set

      .. math::

         \delta {F_\phi}^{cc-c} = \delta F_\phi^{cc-c}
         + \sum A^c \frac{\partial (\delta \phi)^{c-f}}{\partial n}  .

   The gravity synchronization algorithm can be disabled with
   gravity.no_sync = 1. This should be done with care. Generally,
   it is okay only if he refluxing happens in regions of low density that
   don’t affect the gravity substantially.

.. _sec:synchronization_sources:

Source Terms
------------

After a synchronization has been applied, the state on the coarse grid
has changed, due to the change in fluxes at the coarse-fine boundary as
well as the change in the gravitational field. This poses a problem
regarding the source terms, all of which generally rely either on the
state itself, or on the global variables affected by the synchronization
such as the gravitational field. The new-time sources constructed on the
coarse grid all depended on what the state was after the coarse-grid
hydrodynamic update, but the synchronization and associated flux
correction step retroactively changed that hydrodynamic update. So one
can imagine that in a perfect world, we would have calculated the
hydrodynamic update first, including the coarse-fine mismatch
corrections, and only then computed the source terms at the new time.
Indeed, an algorithm that did not subcycle, but marched every zone along
at the same timestep, could do so – and some codes, like FLASH,
actually do this, where no new-time source terms are computed on any
level until the hydrodynamic update has been fully completed and the
coarse-fine mismatches corrected. But in Castro we cannot do this; in
general we assume the ability to subcycle, so the architecture is set up
to always calculate the new-time source terms on a given level
immediately after the hydrodynamic update on that level. Hence on the
coarse level we calculate the new-time source terms before any fine grid
timesteps occur.

One way to fix this, as suggested by Miniati and Colella for the case
of gravity, is to explicitly compute what the difference in the source
term is as a result of any flux corrections across coarse-fine
boundaries. They work out the form of this update for the case of a
cell-centered gravitational force, which has contributions from both
the density advected across the coarse-fine interfaces
(i.e. :math:`\delta \rho \mathbf{g}`, where :math:`\delta \rho` is the density
change due to the coarse-fine synchronization on the coarse rid), as
well as the global change in the gravitational field due to the
collective mass motion (see Miniati and Colella for the explicit form
of the source term). This has a couple of severe limitations. First,
it means that when the form of the source term is changed, the form of
the corrector term is changed too. For example, it is less easy to
write down the form of this corrector term for the flux-based
gravitational energy source term that is now standard in Castro.
Second, gravity is a relatively easy case due to its linearity in the
density and the gravitational acceleration; other source terms
representing more complicated physics might not have an easily
expressible representation in terms of the reflux contribution. For
example, for a general nuclear reaction network (that does not have an
analytic solution), it is not possible to write down an analytic
expression for the nuclear reactions that occur because of
:math:`\delta \rho`.

Instead we choose a more general approach. On the coarse level, we save
the new-time source terms that were applied until the end of the fine
timesteps. We also save the fine level new-time source terms. Then, when
we do the AMR synchronization after a fine timestep, we first subtract
the previously applied new-time source terms to both the coarse and the
fine level, then do the flux correction and associated gravitational
sync solve, and then re-compute the new-time source terms on both the
coarse and the fine level [1]_. In this way, we get almost
the ideal behavior – if we aren’t subcycling, then we get essentially
the same state at the end of the fine timestep as we would in a code
that explicitly had no subcycling. The cost is re-computing the new-time
source terms that second time on each level. For most common source
terms such as gravity, this is not a serious problem – the cost of
re-computing :math:`\rho \mathbf{g}` (for example, once you already know
:math:`\mathbf{g}`) is negligible compared to the cost of actually computing
:math:`\mathbf{g}` itself (say, for self-gravity). If you believe that the
error in not recomputing the source terms is sufficiently low, or the
computational cost of computing them too high, you can disable this
behavior [2]_ using the
code parameter castro.update_sources_after_reflux.

Note that at present nuclear reactions are not enabled as part of this
scheme, and at present are not automatically updated after an AMR
synchronization. This will be amended in a future release of Castro.

.. _sec:synchronization_timing:

Synchronization Timing
----------------------

The goal of the synchronization step is for the coarse and fine grid to
match at the end of a coarse timesteps, after all subcycled fine grid
timesteps have been completed and the two levels have reached the same
simulation time. If subcycling is disabled, so that the coarse and fine
grid take the same timestep, then this is sufficient. However, in the
general subcycling case, the situation is more complicated. Consider the
discussion about source terms in `2.2 <#sec:synchronization_sources>`__. If
we have a coarse level and one fine level with a refinement ratio of
two, then for normal subcycling the fine grid takes two timesteps for
every one timestep taken by the coarse level. The strategy advocated by
the original Castro paper (and Miniati and Colella) is to only do the
AMR synchronization at the actual synchronization time between coarse
and fine levels, that is, at the end of the second fine timestep.
Consequently, we actually only update the source terms after that second
fine timestep. Thus note that on the fine grid, only the *new-time*
source terms in the *second* fine timestep are updated. But a
moment’s thought should reveal a limitation of this. The first fine grid
timestep was also responsible for modifying the fluxes on the coarse
grid, but the algorithm as presented above didn’t take full account of
this information. So, the gravitational field at the old time in
the second fine timestep is actually missing information that would have
been present if we had updated the coarse grid already. Is there a way
to use this information? For the assumptions we make in Castro, the
answer is actually yes. If we apply the effect of the synchronization
not at the synchronization time but at the end of every fine
timestep, then every fine timestep always has the most up-to-date
information possible about the state of the gravitational field. Now, of
course, in fine timesteps before the last one, we have not actually
reached the synchronization time. But we already know at the end of the
first fine timestep what the synchronization correction will be from
that fine timestep: it will be equal to 1/2 of the coarse contribution
to the flux register and the normal contribution to the flux register
for just that timestep. This is true because in Castro, we assume that
the fluxes provided by the hydrodynamic solver are piecewise-constant
over the timestep, which is all that is needed to be second-order
accurate in time if the fluxes are time centered [3]_. So it is fair to say
that halfway through the coarse timestep, half of the coarse flux has
been advected, and we can mathematically split the flux register into
two contributions that have equal weighting from the coarse flux. (In
general, of course, the coarse flux contribution at each fine timestep
is weighted by :math:`1/R` where :math:`R` is the refinement ratio between the
coarse and fine levels.) So, there is nothing preventing us from
updating the coarse solution at the synchronization time :math:`t^{n+1}_c`
after this first fine timestep; we already know at that point how the
coarse solution will change, so why not use that information? We can
then update the gravitational potential at :math:`t^{n+1/2}_c` that is used to
construct the boundary conditions for the gravitational potential solve
on the fine grid at the beginning of the second fine timestep.

In practice, this just means calling the synchronization routine
described in `2.1 <#sec:synchronization_methodology>`__, with the only
modification being that the flux register contribution from the coarse
grid is appropriately weighted by the fine grid timestep instead of
the coarse grid timestep, and we only include the current fine step:

.. math:: \delta\Fb = -\Delta t_f A^c F^c + \Delta t_f A^f F^f

The form of the :math:`\phi` flux register remains unchanged, because the
intent of the gravity sync solve is to simply instantaneously correct
the mismatch between the fine and coarse grid. The only difference,
then, between the old strategy and this new method is that we call the
synchronization at the end of every fine timestep instead of only the
last subcycled one, and we change the weighting appropriately. This
new method is more expensive as currently implemented because we have
to do :math:`R` gravitational sync solves, refluxes, and source term
recalculations instead of only one. However, it results in maximal
possible accuracy, especially in cases where significant amounts of
material are crossing refinement boundaries. The reflux strategy is
controlled by the parameter castro.reflux_strategy. At present
the old method is still the default.

Note that one does not need to be using self-gravity for this to be
beneficial. Even in pure hydrodynamics this can matter. If a regrid
occurs on the fine level, new zones on the boundaries of the current
fine level are filled by interpolation from the coarse level. In the
old method, that interpolation is not using the most up-to-date data
that accounts for the synchronization.

For multiple levels of refinement, the scheme extends naturally. In
the old method, we always call the synchronization at the
synchronization time between any two levels. So for example with two
jumps in refinement by a factor of two, there is a synchronization at
the end of the first two timesteps on level 2 (between level 1 and
level 2), a synchronization after the next two timesteps on level 2
(again between level 1 and level 2), and then a synchronization
between level 0 and level 1. In the new method, we always call the
synchronization at the end of every timestep *on the finest level
only*, and we simultaneously do the synchronization *on every
level*. The timestep :math:`\Delta t_f` in the flux register is just the
timestep on the finest level. (If this is unclear, give it a sanity
check: when the sum of all flux register totals is added up, the level
0 contribution will have a factor of :math:`\Delta t` equal to the coarse
grid timestep since the sum of the timesteps on the finest level over
the entire advance must equal the level 0 timestep. So, the final
contribution from the flux register is the same as if we had saved up
the flux mismatch until the end of the level 0 timestep.) The
synchronization no longer needs to be called at the end of any coarser
level’s timestep because it will already be up to date as a result of
the synchronizations applied at the end of the fine level timesteps.

.. [1]
   In the absence of a global field like
   the gravitational potential, this would only need to be done on the
   coarse level, as we always assume that the solution on the fine grid is
   correct and average it down to the coarse grid. In Castro we do it by
   default on the fine level too in anticipation of the fact that gravity
   is a common component of many of our production science
   simulations. This could be generalized so that if you aren’t using any
   global force fields, you don’t bother updating the fine level. If this
   is important to the science you want to do, please let the Castro developers know and we can look into it.

.. [2]
   in general it may be desirable for this to be a
   source-term specific setting, so that some source terms that are cheap
   or physically important are re-computed after a synchronization can be
   set to update, while others can be disabled. If this is important for
   your science application, please let the developers know, as this would
   be a straightforward extension of the current architecture.

.. [3]
   If this scheme
   is generalized to higher-order methods, in principle all one would need
   to do is integrate the fluxes until :math:`\Delta t / 2`, which is what we are
   doing here for the constant-in-time flux case.
