#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <initial_model.H>
#include <global.H>

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

    // we use the fine grid dx for the model resolution
    auto fine_geom = global::the_amr_ptr->Geom(global::the_amr_ptr->maxLevel());

    auto dx = fine_geom.CellSizeArray();
    auto dx_model = dx[AMREX_SPACEDIM-1];

    int nx = (2.0_rt * problem::center[AMREX_SPACEDIM-1] + 1.e-8_rt) /
        dx_model;

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

#endif
