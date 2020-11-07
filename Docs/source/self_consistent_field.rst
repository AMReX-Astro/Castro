************************************
Self-Consistent Field Initialization
************************************

Introduction
============

The Hachisu self-consistent field (SCF) method is a way to generate
equilibrium rotating (or non-rotating) initial configurations. It can
generate single or multiple star systems. The SCF method was originally
developed in the 1960s, but it was a variant proposed by Hachisu in 1986
(:cite:`hachisu:1986a,hachisu:1986b`) that ended up becoming the
most popular technique. The SCF method was originally developed
for rapidly rotating single stars, but was soon extended to apply
to rotating binary systems, and it has been applied to construct
initial conditions by several groups studying binary white dwarf
or neutron star systems using several types of rotation laws
(:cite:`newtohline:1997,swc:2000,motl:2002,dsouza:2006,motl:2007,even:2009,kadam:2018,yoshida:2018`).
It has also been used for constructing toroidal configurations
(:cite:`kim:2016`).
The technique assumes a uniform temperature and composition (more generally,
a barotropic equation of state), self-gravitation represented by the
Poisson equation, and a well-defined rotation law (often rigid-body rotation).
The user is required to specify three quantities: the maximum density of
the equilibrium star(s), and two points on the stellar surface. For a
single star this is usually a point on the equator and a point on a pole.
For a detached binary system this corresponds to the inner and outer points
of the teardrop-shape configuration along the axis joining the binary.

At present we only support generating a single star using the Hachisu
SCF method, but we plan to extend this in a later release to rotating
toroidal configurations and binary star systems.

We note that while the Hachisu SCF method iteratively solves the
integral form of the Euler equations (essentially, it solves the Bernoulli
equation from classical fluid dynamics), other approaches have been used
in the literature. Usually these directly solve the coupled Poisson and
Bernoulli equations rather than relying on indirect iterative coupling
between them. See, for example, :cite:`eriguchi:1985,fujisawa:2015` and
:cite:`clement:1974,aksenov:1994`. These alternative methods are generally
more powerful and promise faster convergence, but are also more difficult
to implement. See also :cite:`jackson:2005` for yet another approach.



Usage and Code Parameters
=========================

To use the SCF initialization technique, you must compile with
``USE_GRAV=TRUE`` and ``USE_ROTATION=TRUE``, and enable the method
with ``castro.do_scf_initial_model`` = ``1``. You are responsible
for initializing the grid with an initial guess at the mass distribution.
This guess does not need to be accurate -- for example, just initializing
with a uniform spherical mass distribution should be sufficient.
(See ``Exec/scf_tests/single_star``.) However, it is important that
the grid is isothermal and has a uniform composition. This temperature
and composition will be retained for the final equilibrium configuration.

Several code parameters are available for controlling problem initialization
with SCF:

- ``castro.scf_maximum_density``: the target maximum density on the domain
- ``castro.scf_equatorial_radius``: the target equatorial radius of the star
- ``castro.scf_polar_radius``: the target polar radius of the star
- ``castro.scf_relax_tol``: tolerance required for SCF convergence

The first three options are required and must be set. One limitation of this
method is that (to our knowledge) there is no known way to specify more natural
parameters such as the total mass of the star.


Single Star Algorithm
=====================

The Bernoulli equation tells us

.. math::
   \Phi + \phi + H = \mathrm{constant}

where :math:`\Phi` is the gravitational potential, :math:`\phi` is the
rotational potential (:math:`\phi = -\omega^2 R^2` for rigid-body rotation,
where :math:`R` is the radius in the equatorial plane), and :math:`H` is the
enthalpy.

We denote the prescribed maximum density as :math:`\rho_0` and the equatorial
and polar radii as :math:`r_A` and :math:`r_B` respectively (the A and B indices
are chosen to be consistent with Hachisu). At these radii the enthalpy should
vanish, so we have:

.. math::
   C = \Phi_A + \phi_A

.. math::
   C = \Phi_B + \phi_B

Note that :math:`C = C_A = C_B` is constant everywhere on the body.

Then, given the value of a Bernoulli constant, we can invert the Bernoulli equation to
obtain the enthalpy:

.. math::
   H = C - \Phi - \phi

Given an enthalpy, we can call the equation of state given the enthalpy (and
composition and temperature) as an input to obtain the density. (Note that
this inversion generally requires a Newton-Raphson iteration in realistic
equations of state.)

.. math::
   \rho = \rho(T, H, X)

As one more housekeeping item, we'll notate the rotational potential as

.. math::
   \phi = \omega^2 \psi

where :math:`\psi = -R^2` for rigid body rotation.

The above ingredients are all that is needed to construct the algorithm for
obtaining an equilibrium rotating single star. Given an initial density distribution,
:math:`\rho^n`, gravitational field, :math:`\Phi^n`, and rotational field,
:math:`\phi^n`, we first calculate an updated guess for the rotation frequency
:math:`\omega`:

.. math::
   \omega^{n+1} = \sqrt{\frac{\Phi_B^n - \Phi_A^n}{\psi_A^n - \psi_B^n}}

which simply involves finding :math:`\Phi` and :math:`\psi` at these vanishing points.

With the updated rotation frequency, we can reconstruct the rotational potential
:math:`\phi`, and then update the enthalpy everywhere on the domain as:

.. math::
   H^{n+1} = C - \Phi - \Phi_R

However, we want to guarantee that the maximum density on the domain is fixed. Given
that this maximum density corresponds to a maximum enthalpy,

.. math::
   H_0 = H(\rho_0, T, X)

we can rescale all of the updated enthalpies such that the maximum is fixed:

.. math::
   H^{n+1} \rightarrow H^{n+1} \left( \frac{H_0}{H^{n+1}_{\mathrm{max}}} \right)

and then invert the EOS to obtain :math:`\rho^{n+1}`. Given the new density
distribution, we can then update the gravitational potential, :math:`\Phi^{n+1}`,
by solving the Poisson equation. This procedure is iterated until no zone
changes its density by more than a factor of ``castro.scf_relax_tol``.
