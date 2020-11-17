#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // Read in the initial model

    read_model_file(problem::model_name);

    // Save some of the data locally

#if AMREX_SPACEDIM == 1
    // 1-d assumes spherical, with center at origin
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.0_rt;
    problem::center[2] = 0.0_rt;
#elif AMREX_SPACEDIM == 2
    // 2-d assumes axisymmetric
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.0_rt;
#else
    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

#if AMREX_SPACEDIM == 1
  if (problo[0] != 0.0_rt) {
      amrex::Error("ERROR: xmin should be 0!");
  }

#elif AMREX_SPACEDIM == 2
  if (problo[0] != 0.0_rt) {
      amrex::Error("ERROR: xmin should be 0!");
  }
#endif

  // store the state at the very top of the model for the boundary
  // conditions

  problem::hse_rho_top = model::state(model::npts-1, model::idens);
  problem::hse_T_top = model::state(model::npts-1, model::itemp);
  for (int n = 0; n < NumSpec; n++) {
      problem::hse_X_top[n] = model::state(model::npts-1, model::ispec+n);
  }

  // set hse_eint_top and hse_p_top via the EOS

  eos_t eos_state;
  eos_state.rho = problem::hse_rho_top;
  eos_state.T = problem::hse_T_top;
  for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = problem::hse_X_top[n];
  }

  eos(eos_input_rt, eos_state);

  problem::hse_eint_top = eos_state.e;
  problem::hse_p_top = eos_state.p;

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif


#if AMREX_SPACEDIM == 1
    Real dist = x;
#elif AMREX_SPACEDIM == 2
    Real dist = std::sqrt(x*x + y*y);
#else
    Real dist = std::sqrt(x*x + y*y + z*z);
#endif

    state(i,j,k,URHO) = interpolate(dist, model::idens);
    state(i,j,k,UTEMP) = interpolate(dist, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(dist, model::ispec+n);
    }

    // Compute energy from rho,T and X

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rt, eos_state);

    // we'll add the density weighting shortly

    state(i,j,k,UEINT) = eos_state.e;

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT);
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }
}
#endif


