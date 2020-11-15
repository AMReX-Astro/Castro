#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <initial_model.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // set center variable

    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
#if AMREX_SPACEDIM >= 2
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
#endif
#if AMREX_SPACEDIM == 3
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

    // set local variable defaults

    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
#if AMREX_SPACEDIM >= 2
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
#endif
#if AMREX_SPACEDIM == 3
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

    // we'll generate the initial model at the needed resolution

    int nx = (2.0_rt * problem::center[AMREX_SPACEDIM-1] + 1.e-8_rt) /
        problem::dx_model;

    model_t model_params;
    model_params.dens_base = problem::dens_base;
    model_params.p_base = problem::pres_base;
    for (int n = 0; n < NumSpec; ++n) {
        model_params.xn[n] = 0.0_rt;
    }
    model_params.xn[0] = 1.0_rt;

    generate_initial_model(nx, problo[AMREX_SPACEDIM-1], probhi[AMREX_SPACEDIM-1],
                           model_params);

    if (problem::single == 1) {
        problem::left_bubble_x_center = problo[0] + 0.5_rt * (probhi[0] - problo[0]);
    } else {
        problem::left_bubble_x_center = problo[0] + (probhi[0] - problo[0]) / 3.0_rt;
        problem::right_bubble_x_center = problo[0] + 2.0_rt * (probhi[0] - problo[0]) / 3.0_rt;
    }

}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 1
    Real height = x;
#elif AMREX_SPACEDIM == 2
    Real height = y;
#else
    Real height = z;
#endif

    Real T_height = interpolate(height, model::itemp);
    Real dens = interpolate(height, model::idens);
    Real pres = interpolate(height, model::ipres);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(height, model::ispec+n);
    }

    if (problem::single == 0) {

#if AMREX_SPACEDIM == 1
        amrex::Error("Error: 1-d not supported");

#elif AMREX_SPACEDIM == 2
        Real x1 = problem::left_bubble_x_center;
        Real y1 = problem::r_pert_center;
        Real z1 = 0.0_rt;

        Real x2 = problem::right_bubble_x_center;
        Real y2 = problem::r_pert_center;
        Real z2 = 0.0_rt;

#else
        Real x1 = problem::left_bubble_x_center;
        Real y1 = problem::center[1];
        Real z1 = problem::r_pert_center;

        Real x2 = problem::right_bubble_x_center;
        Real y2 = center[1];
        Real z2 = problem::r_pert_center;
#endif

        Real r1 = std::sqrt((x - x1) * (x - x1) +
                            (y - y1) * (y - y1) +
                            (z - z1) * (z - z1)) / problem::pert_width;

        Real r2 = std::sqrt((x - x2) * (x - x2) +
                            (y - y2) * (y - y2) +
                            (z - z2) * (z - z2)) / problem::pert_width;

        eos_t eos_state;

        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = state(i,j,k,UFS+n);
        }

        state(i,j,k,UTEMP) = T_height;
        state(i,j,k,URHO) = dens;

        // which bubble are we in? -- we want their rho perturbations to be the
        // same so they have the same buoyancy

        if (r1 < 2.0_rt) {
            state(i,j,k,URHO) = dens * (1.0_rt - (problem::pert_factor * (1.0_rt + std::tanh(2.0_rt - r1))));
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = 0.0_rt;
            }
            eos_state.xn[1] = 1.0_rt;
        }

        if (r2 < 2.0_rt) {
            state(i,j,k,URHO) = dens * (1.0_rt - (problem::pert_factor * (1.0_rt + std::tanh(2.0_rt - r2))));
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = 0.0_rt;
            }
            eos_state.xn[2] = 1.0_rt;
        }

        eos_state.p = pres;
        eos_state.rho = state(i,j,k,URHO);
        eos_state.T = T_height;

        eos(eos_input_rp, eos_state);

        state(i,j,k,UEINT) = eos_state.e;
        state(i,j,k,UTEMP) = eos_state.T;

        // make state conservative

        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = state(i,j,k,URHO) * eos_state.xn[n];
        }

        state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT);

        // assumes ke=0

        state(i,j,k,UEDEN) = state(i,j,k,UEINT);

        state(i,j,k,UMX) = 0.0;
        state(i,j,k,UMY) = 0.0;
        state(i,j,k,UMZ) = 0.0;

    } else {

#if AMREX_SPACEDIM == 1
        amrex::Error("Error: 1-d not supported");
#elif AMREX_SPACEDIM == 2
        Real x1 = problem::left_bubble_x_center;
        Real y1 = problem::r_pert_center;
        Real z1 = problem::center[2];
#else
        Real x1 = problem::left_bubble_x_center;
        Real y1 = problem::center[1];
        Real z1 = problem::r_pert_center;
#endif

        Real r1 = std::sqrt((x - x1) * (x - x1) +
                            (y - y1) * (y - y1) +
                            (z - z1) * (z - z1)) / problem::pert_width;

        eos_t eos_state;

        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = state(i,j,k,UFS+n);
        }

        state(i,j,k,UTEMP) = T_height;
        state(i,j,k,URHO) = dens;


        // which bubble are we in? -- we want their rho perturbations to be the
        // same so they have the same buoyancy

        if (r1 < 2.0_rt) {
            state(i,j,k,URHO) = dens * (1.0_rt - (problem::pert_factor * (1.0_rt + std::tanh(2.0_rt - r1))));
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = 0.0_rt;
            }
            eos_state.xn[1] = 1.0_rt;
        }

        eos_state.p = pres;
        eos_state.rho = state(i,j,k,URHO);
        eos_state.T = T_height;

        eos(eos_input_rp, eos_state);

        state(i,j,k,UEINT) = eos_state.e;
        state(i,j,k,UTEMP) = eos_state.T;


        // make state conservative

        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = state(i,j,k,URHO) * eos_state.xn[n];
        }

        state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT);

        // assumes ke=0

        state(i,j,k,UEDEN) = state(i,j,k,UEINT);

        state(i,j,k,UMX) = 0.0;
        state(i,j,k,UMY) = 0.0;
        state(i,j,k,UMZ) = 0.0;

    }
}
#endif


