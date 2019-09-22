Treatment of gas pressure interface state in Castro + radiation

Background
----------

The interface states for radiation work in the primitive variable
system, :math:`q = (\rho, u, p, (\rho e)_g, E_r)^\intercal`, where :math:`p` is
the gas pressure only, and :math:`(\rho e)_g` is the gas energy density.

Written in the form:

.. math:: q_t + A(q) q_x = 0

The matrix :math:`A` takes the form:

.. math::

   A = \left (
   \begin{matrix}
   u & \rho & 0 & 0 & 0\\
   0 & u & {1}/{\rho} & 0 & \lambda_f/{\rho}\\
   0 & c_{g}^{2} \rho & u & 0 & 0\\
   0 & h_{g} \rho & 0 & u & 0\\
   0 & E_{r} \left(\lambda_f + 1\right) & 0 & 0 & u
   \end{matrix}\right )

here, :math:`c_g` is the gas sound speed and :math:`h_g = e_g + p/\rho` is the
gas specific enthalpy.
A has eigenvalues, :math:`\lambda^{(-)}= u - c`, :math:`\lambda^{(0)}= u`,
:math:`\lambda^{(+)}= u + c`, where the total sound speed, :math:`c` is related
to the gas sound speed, :math:`c_g`, as:

.. math:: c^2 = c_g^2 + (\lambda_f + 1)\frac{\lambda_f E_r}{\rho}

In constructing the interface states, start with a reference state,
:math:`q_\mathrm{ref}` and define the jumps carried by each wave as the
integral under the parabolic profile with respect to this reference
state:

.. math:: \Delta q^{(\nu)}\equiv q_\mathrm{ref} - \mathcal{I}^{(\nu)}(q)

and then the interface states are:

.. math:: q_\mathrm{int} = q_\mathrm{ref} - \sum_\nu (l^{(\nu)}\cdot \Delta q^{(\nu)}) r^{(\nu)}

Defining:

.. math::

   \begin{aligned}
   \beta^{(-)}= ( l^{(-)}\cdot \Delta q^{(-)}) &= \frac{1}{2 c^{2}} \left(
       \Delta E_r^{(-)}\lambda_f + \Delta p^{(-)}- \Delta u^{(-)}c \rho\right)\\
   &= \frac{1}{2 c^{2}} \left(
       \Delta p_\mathrm{tot}^{(-)}- \Delta u^{(-)}c \rho\right)\end{aligned}

where we recognized that :math:`p_\mathrm{tot} = p + \lambda_f E_r`.
Similarly, we have:

.. math::

   \begin{aligned}
   \beta^{(+)}= ( l^{(+)}\cdot \Delta q^{(+)}) &= \frac{1}{2 c^{2}} \left(
       \Delta E_r^{(+)}\lambda_f + \Delta p^{(+)}+ \Delta u^{(+)}c \rho\right)\\
              &= \frac{1}{2 c^{2}} \left(
       \Delta p_\mathrm{tot}^{(+)}+ \Delta u^{(+)}c \rho\right)\end{aligned}

and for the 0-wave, we have three degenerate eigenvalues and
eigenvectors. We label these with the subscripts :math:`\rho`, :math:`\rho e_g`,
and :math:`E_r`, indicating which of the corresponding components in our
primitive variable vector, :math:`q`, is non-zero in the right eigenvector,
:math:`r^{(\nu)}`. Then we have:

.. math::

   \begin{aligned}
   \beta^{(0)}_\rho &= 
       \Delta\rho^{(0)}- \frac{\Delta p^{(0)}_\mathrm{tot}}{c^2} \\
   %
   \beta^{(0)}_{{\rho e}_g} &= \Delta(\rho e)^{(0)}_g - \frac{\Delta p_\mathrm{tot}^{(0)}}{c^2} h_g \\
   %
   \beta^{(0)}_{E_r} &= \Delta E_r^{(0)}- \frac{\Delta p_\mathrm{tot}^{(0)}}{c^2} h_r\end{aligned}

where :math:`h_r = (\lambda_f + 1)E_r/\rho`. Note, these match the derivation
done in the Jupyter/SymPy notebook:

The gas pressure update in these terms is:

.. math:: p_\mathrm{int} = p_\mathrm{ref} - (\beta^{(+)}+ \beta^{(-)}) c_g^2 + \lambda_f \beta^{(0)}_{E_r}

This matches the expression in Castro. Notice that this expression is
unusual for pressure, since it jumps across the :math:`{(0)}`-wave, whereas
the total pressure would not.

Castro computes the edge state of the radiation energy density as:

.. math:: {E_r}_\mathrm{int} = {E_r}_\mathrm{ref} - (\beta^{(+)}+ \beta^{(-)}) h_r - \beta^{(0)}_{E_r}

and we see that the total pressure can be constructued as
:math:`{p_\mathrm{tot}}_\mathrm{int} = p_\mathrm{int} + \lambda_f
{E_r}_\mathrm{int}`, giving:

.. math::

   {p_\mathrm{tot}}_\mathrm{int} = {p_\mathrm{tot}}_\mathrm{ref} -
      (\beta^{(+)}+ \beta^{(-)}) c^2

This looks, as it should, analogous to the pure hydrodynamics case.

The Interface States
--------------------

What is the interface state? for this we need to choose a reference
state. The choice of reference state should not matter if we counted
the contributions of all waves, then:

.. math:: q_\mathrm{int} = q_\mathrm{ref} - \sum_\nu [ l^{(\nu)}\cdot (q_\mathrm{ref} - \mathcal{I}^{(\nu)}(q)) ] r^{(\nu)}

and :math:`q_\mathrm{ref}` cancels out. However, when we do characteristic
tracing, this is not guaranteed.

In Castro, the quantity :math:`\Delta p^{(0)}_\mathrm{tot}` is defined as:

.. math::

   \Delta p^{(0)}_\mathrm{tot} =
     p_\mathrm{tot,ref} - \mathcal{I}^{(0)}(p_\mathrm{tot})

and we adopt the common strategy of picking the reference state to be
the :math:`\mathcal{I}` corresponding to the fastest wave moving toward the
interface.

Consider a zone, :math:`i`, and tracing to the right edge of the zone to
form the interface state, :math:`i+1/2,L`, where :math:`L` indicates that it is
immediately to the left of the :math:`i+1/2` interface. The fastest wave
that can potentially move toward that interface is the :math:`{(+)}` wave,
so we pick:

.. math::

   q_\mathrm{ref} = \left (
      \begin{array}{c}
        \mathcal{I}^{(+)}(\rho) \\
        \mathcal{I}^{(+)}(u) \\
        \mathcal{I}^{(+)}(p_\mathrm{tot}) \\
        \mathcal{I}^{(+)}((\rho e)_g) \\
        \mathcal{I}^{(+)}(E_r)
      \end{array}
   \right )

Looking at the gas pressure interface state, we have:

.. math::

   p_\mathrm{int} = p_\mathrm{ref} - \frac{1}{2} \frac{c_g^2}{c^2} \left \{
      \left [ \Delta p_\mathrm{tot}^{(+)}+ \Delta u^{(+)}c\rho \right ]
    + \left [ \Delta p_\mathrm{tot}^{(-)}- \Delta u^{(-)}c\rho \right ] \right \}
    + \lambda_f \left [ \Delta E_r^{(0)}- \frac{h_r}{c^2} \Delta p^{(0)}_\mathrm{tot} \right ]

Substituting in our choice of reference state, we have:

.. math::

   \begin{aligned}
   p_\mathrm{int} = \mathcal{I}^{(+)}(p) &- \underbrace{\frac{1}{2} \frac{c_g^2}{c^2} \left \{
    \left [ \mathcal{I}^{(+)}(p_\mathrm{tot}) - \mathcal{I}^{(-)}(p_\mathrm{tot}) \right ] 
      - \rho c \left [ \mathcal{I}^{(+)}(u) - \mathcal{I}^{(-)}(u) \right ] \right \}}_{\mbox{\footnotesize carried by the ${(-)}$ wave}} \\
    &+ \underbrace{\lambda_f \left \{ \left [ \mathcal{I}^{(+)}(E_r) - \mathcal{I}^{(0)}(E_r) \right ]
      - \frac{h_r}{c^2} \left [ \mathcal{I}^{(+)}(p_\mathrm{tot}) - \mathcal{I}^{(0)}(p_\mathrm{tot}) \right ] \right \}}_{\mbox{\footnotesize carried by the ${(0)}$ wave}}\end{aligned}

We see that the expression for
:math:`p_\mathrm{int}` starts with the gas pressure. In the
event that no other waves are moving toward the interface, then we find:

.. math:: p_\mathrm{int} \rightarrow \mathcal{I}^{(+)}(p)

(since in our algorithms, we set the :math:`\beta`\ ’s to :math:`0` if those waves
are not moving toward the interface.

Alternative?
------------

A, perhaps more consistent, way to handle this is to predict :math:`p_\mathrm{tot}`
and :math:`E_r` to the interfaces. This is consistent with our choice of
reference state, and using :math:`p_\mathrm{tot}` is more analogous to the
pure hydrodynamics case. We then algebraically construct the
gas-pressure interface state as:

.. math:: p_\mathrm{int} = {p_\mathrm{tot}}_\mathrm{int} - \lambda_f {E_r}_\mathrm{int}

Note that this extends naturally to multigroup radiation as well, simply by
summing up the independent :math:`E_r`.

:math:`\gamma_e` system
-----------------------

The alternate approach that Castro takes to incorporate auxiliary
thermodynamic information into the system is to evolve an equation for
:math:`\gamma_e`. We use the primitive variable system, :math:`q = (\tau, u, p,
\gamma_e, E_r)^\intercal`, where :math:`\tau = 1/\rho`. The matrix :math:`A` now
takes the form:

.. math::

   A = \left (
      \begin{matrix}
      u & - \tau & 0 & 0 & 0\\
      0 & u & \tau & 0 & \lambda_f \tau\\
      0 & \frac{c_{g}^{2}}{\tau} & u & 0 & 0\\
      0 & - \alpha & 0 & u & 0\\
      0 & E_{r} \left(\lambda_f + 1\right) & 0 & 0 & u
   \end{matrix}\right)

The eigenvalues are unchanged and the eigenvectors are derived in the
same Jupyter notebook as with the previous system.
Here, :math:`\alpha = (\gamma_e - 1)(\gamma_e - \Gamma_1)`.
Now we have

.. math::

   \begin{aligned}
   \beta^{(+)}= ( l^{(+)}\cdot \Delta q^{(+)}) &= -\frac{1}{2C}
      \left ( \Delta u^{(+)}+ \frac{\Delta p_\mathrm{tot}^{(+)}}{C} \right ) \\
   \beta^{(0)}_\tau = ( l^{(0)}_\tau \cdot \Delta q^{(0)}) &= 
      \Delta \tau^{(0)}+ \frac{\Delta p_\mathrm{tot}^{(0)}}{C^2} \\
   \beta^{(0)}_{\gamma_e} = ( l^{(0)}_{\gamma_e} \cdot \Delta q^{(0)}) &= 
      \Delta \gamma_E^{(0)}+ \alpha \frac{\Delta p_\mathrm{tot}^{(0)}}{\tau C^2} \\
   \beta^{(0)}_{E_r} = ( l^{(0)}_{E_r} \cdot \Delta q^{(0)}) &= 
      \Delta E_r^{(0)}- \frac{h_r}{c^2} \Delta p_\mathrm{tot}^{(0)}\\
   \beta^{(-)}= ( l^{(-)}\cdot \Delta q^{(-)}) &= \frac{1}{2C} 
      \left ( \Delta u^{(-)}- \frac{\Delta p_\mathrm{tot}^{(-)}}{C} \right ) \end{aligned}

Here we use the Lagrangian sound speed, :math:`C = \rho c = c/\tau`.

The interface states are then:

.. math::

   \begin{aligned}
   \tau_\mathrm{int} &= \tau_\mathrm{ref} - \beta^{(+)}- \beta^{(-)}- \beta_\tau^{(0)}\\
   u_\mathrm{int} &= u_\mathrm{ref} + C (\beta^{(+)}- \beta^{(-)}) \\
   p_\mathrm{int} &= p_\mathrm{ref} + \frac{c_g^2}{\tau^2} ( \beta^{(+)}+ \beta^{(-)}) + \beta_{E_r} \lambda_f \\
   {\gamma_e}_\mathrm{int} &= {\gamma_e}_\mathrm{ref} - \beta_{\gamma_e}^{(0)}
      - \frac{\alpha}{\tau} (\beta^{(+)}+ \beta^{(-)}) \\
   {E_r}_\mathrm{int} &= {E_r}_\mathrm{ref} + \frac{h_r}{\tau^2} (\beta^{(+)}+ \beta^{(-)}) - \beta^{(0)}_{E_r}\end{aligned}

Again, we can also construct the total pressure on the interface:

.. math::

   \begin{aligned}
   {p_\mathrm{tot}}_\mathrm{int} &= p_\mathrm{int} + \lambda_f {E_r}_\mathrm{int}\\
       &= p_\mathrm{ref} + \lambda {E_r}_\mathrm{ref} + \frac{h_r \lambda_f + c_g^2}{\tau^2} (\beta^{(+)}+ \beta^{(-)}) \\
       &= {p_\mathrm{tot}}_\mathrm{ref} + C^2 (\beta^{(+)}+ \beta^{(-)})\end{aligned}
