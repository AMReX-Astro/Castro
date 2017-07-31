# Problem Setup

This is a pure diffusion problem (no hydro).  It uses the explicit
diffusion solver in Castro to diffuse a Gaussian thermal pulse.
Because diffusion in Castro is incorporated in the energy equation, we
are solving:

d(rho e)/dt = k L T

where L is the Laplacian and k is the thermal conductivity.  Here we
assume that k is constant, but Castro does not require that
assumption.  For this problem, we take rho = constant, and using a
gamma-law EOS, we have e = c_v T, so this can take the form of a
diffusion equation as:

dT/dt = k/(rho c_v) L T = D L T

where D is the constant diffusion coefficient.

Because we are doing the diffusion explicitly, there is a constraint
on the timestep,

dt < 0.5 dx**2/D


# Testing

The problem creates a derived variable, `analytic`, which is the
analytic solution at the time of the plotfile output.  This allows you
to compare the current solution to the analytic solution to assess the
error.

It also will use the Castro problem-specific post-simulation hooks (in
`Prob.cpp`) to output the L-inf norm of the error (numerical vs
analyic solution) at the end of the simulation.  This can be used
for convergence testing.
