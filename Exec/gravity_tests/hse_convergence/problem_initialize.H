#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <initial_model.H>
#include <network.H>
#include <global.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
    }


    int ihe4 = network_spec_index("helium-4");
    if (ihe4 < 0) {
        amrex::Error("Error: helium-4 not present");
    }

    model_t model_params;
    model_params.T_base = problem::temp_base;
    model_params.dens_base = problem::dens_base;
    for (int n = 0; n < NumSpec; n++) {
        model_params.xn[n] = 100 * network_rp::small_x;
    }
    model_params.xn[ihe4] = 1.0_rt - (NumSpec - 1) * 100 * network_rp::small_x;

    // we add some buffer to the model so we can use it to fill ghost
    // cells in the boundary conditions

    int nbuf = 8;

    // we use the fine grid dx for the model resolution
    auto fine_geom = global::the_amr_ptr->Geom(global::the_amr_ptr->maxLevel());

    auto dx = fine_geom.CellSizeArray();
    auto dx_model = dx[AMREX_SPACEDIM-1];

    int nx_model = static_cast<int>((probhi[AMREX_SPACEDIM-1] -
                                     problo[AMREX_SPACEDIM-1] + 1.e-8_rt) / dx_model);

    // generate the initial model -- it will be stored in the model
    // parser global data

    generate_initial_model(nx_model,
                           problo[AMREX_SPACEDIM-1], probhi[AMREX_SPACEDIM-1],
                           model_params, nbuf);

}

#endif
