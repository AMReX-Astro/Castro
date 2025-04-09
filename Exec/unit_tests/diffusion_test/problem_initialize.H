#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <prob_util.H>
#include <conductivity.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // set center variable

    // Note problem::center = 0 by default for 2d spherical geometry.
    if ((problem::center[0] == 0.0_rt) && (coord_type == 0)) {
        problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    }
#if AMREX_SPACEDIM >= 2
    if ((problem::center[1] == 0.0_rt) && (coord_type <= 1)) {
        problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    }
#endif
#if AMREX_SPACEDIM == 3
    if (problem::center[2] == 0.0_rt) {
        problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
    }
#endif

    // the conductivity is the physical quantity that appears in the
    // diffusion term of the energy equation.  It is set in the inputs
    // file.  For this test problem, we want to set the diffusion
    // coefficient, D = k/(rho c_v), so the free parameter we have to
    // play with is rho.  Note that for an ideal gas, c_v does not depend
    // on rho, so we can call it the EOS with any density.
    Real X[NumSpec] = {0.0_rt};
    X[0] = 1.0_rt;

    eos_t eos_state;
    eos_state.T = problem::T1;
    eos_state.rho = 1.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = X[n];
    }

    eos(eos_input_rt, eos_state);

#ifdef DIFFUSION
    // get the conductivity
    conductivity(eos_state);

    // diffusion coefficient is D = k/(rho c_v). we are doing an ideal
    // gas, so c_v is constant, so find the rho that combines with
    // the conductivity
    problem::rho0 = eos_state.conductivity / (problem::diff_coeff * eos_state.cv);
#else
    problem::rho0 = 1.0_rt;
#endif

}

#endif
