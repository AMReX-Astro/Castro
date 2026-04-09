#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();

#if AMREX_SPACEDIM >= 2
    const Real* probhi = dgeom.ProbHi();
#endif

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

  problem::hse_rho_top = model::profile(0).state(model::npts-1, model::idens);
  problem::hse_T_top = model::profile(0).state(model::npts-1, model::itemp);
  for (int n = 0; n < NumSpec; n++) {
      problem::hse_X_top[n] = model::profile(0).state(model::npts-1, model::ispec+n);
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

  problem::ifuel = network_spec_index("helium-4");

  // make sure we have the minimum nuclei needed to understand the
  // different initial models depending on the prob_type

  if (problem::prob_type == 1) {
      if (network_spec_index("carbon-12") < 0 ||
          network_spec_index("oxygen-16") < 0 ||
          network_spec_index("magnesium-24") < 0) {
          amrex::Error("invalid network -- required nuclei not preset");
      }

  } else if (problem::prob_type == 3) {
      if (network_spec_index("helium-4") < 0 ||
          network_spec_index("carbon-12") < 0 ||
          network_spec_index("oxygen-16") < 0 ||
          network_spec_index("iron-56") < 0) {
          amrex::Error("invalid network -- required nuclei not preset");
      }
  }

}

#endif
