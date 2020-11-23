#ifndef problem_setup_H
#define problem_setup_H

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

    if (coord_type == 0) {
        problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    } else if (coord_type >= 1) {
        problem::center[0] = 0.0;
    }
#if AMREX_SPACEDIM >= 2
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
#endif
#if AMREX_SPACEDIM == 3
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

    // the conductivity is the physical quantity that appears in the
    // diffusion term of the energy equation.  It is set in the probin
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


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    int coord_type = geomdata.Coord(); 

    Real r[3] = {0.0_rt};

    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    // set the composition

    Real X[NumSpec] = {0.0_rt};
    X[0] = 1.0_rt;

    state(i,j,k,URHO) = problem::rho0;

    Real temp = analytic(r, 0.0_rt, coord_type);

    state(i,j,k,UTEMP) = temp;

    // compute the internal energy and temperature

    eos_t eos_state;
    eos_state.T = temp;
    eos_state.rho = state(i,j,k,URHO);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = X[n];
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UEDEN) = problem::rho0 * eos_state.e;
    state(i,j,k,UEINT) = problem::rho0 * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = problem::rho0 * X[n];
    }

}
#endif
