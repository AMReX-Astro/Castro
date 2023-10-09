#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <model_parser.H>
#include <initial_model.H>
#include <ambient.H>
#include <global.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    problem::refine_cutoff_height = problo[1] + problem::refine_cutoff_frac * (probhi[1] - problo[1]);


    // get the species indices

    bool species_defined = true;

    int ifuel1 = network_spec_index(problem::fuel1_name);
    if (ifuel1 < 0) {
        species_defined = false;
    }

    int ifuel2;
    if (problem::fuel2_name != "") {
        ifuel2 = network_spec_index(problem::fuel2_name);
        if (ifuel2 < 0) {
            species_defined = false;
        }
    }

    int ifuel3;
    if (problem::fuel3_name != "") {
        ifuel3 = network_spec_index(problem::fuel3_name);
        if (ifuel3 < 0) {
            species_defined = false;
        }
    }

    int iash1 = network_spec_index(problem::ash1_name);
    if (iash1 < 0) {
        species_defined = false;
    }

    int iash2;
    if (problem::ash2_name != "") {
        iash2 = network_spec_index(problem::ash2_name);
        if (iash2 < 0) {
            species_defined = false;
        }
    }

    int iash3;
    if (problem::ash3_name != "") {
        iash3 = network_spec_index(problem::ash3_name);
        if (iash3 < 0) {
            species_defined = false;
        }
    }

    if (! species_defined) {
        std::cout << ifuel1 << " " << ifuel2 << " " << ifuel3 << std::endl;
        std::cout << iash1 << " " << iash2 << " "<< iash3 << std::endl;
        amrex::Error("ERROR: species not defined");
    }

    model_t model_params;

    // set the composition of the underlying star

    for (int n = 0; n < NumSpec; n++) {
        model_params.xn_star[n] = problem::smallx;
    }
    model_params.xn_star[iash1] = problem::ash1_frac;
    if (problem::ash2_name != "") {
        model_params.xn_star[iash2] = problem::ash2_frac;
    }
    if (problem::ash3_name != "") {
        model_params.xn_star[iash3] = problem::ash3_frac;
    }

    // and the composition of the accreted layer

    for (int n = 0; n < NumSpec; n++) {
        model_params.xn_base[n] = problem::smallx;
    }
    model_params.xn_base[ifuel1] = problem::fuel1_frac;
    if (problem::fuel2_name != "") {
        model_params.xn_base[ifuel2] = problem::fuel2_frac;
    }
    if (problem::fuel3_name != "") {
        model_params.xn_base[ifuel3] = problem::fuel3_frac;
    }

    // check if they sum to 1

    Real sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += model_params.xn_star[n];
    }
    if (std::abs(sumX - 1.0_rt) > NumSpec * problem::smallx) {
        amrex::Error("ERROR: ash mass fractions don't sum to 1");
    }

    sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += model_params.xn_base[n];
    }
    if (std::abs(sumX - 1.0_rt) > NumSpec * problem::smallx) {
        amrex::Error("ERROR: fuel mass fractions don't sum to 1");
    }

    // we are going to generate an initial model from problo(2) to
    // probhi(2) with nx_model zones.  But to allow for a interpolated
    // lower boundary, we'll add 4 ghostcells to this, so we need to
    // compute dx

    // we use the fine grid dx for the model resolution
    auto fine_geom = global::the_amr_ptr->Geom(global::the_amr_ptr->maxLevel());

    auto dx = fine_geom.CellSizeArray();
    auto dx_model = dx[AMREX_SPACEDIM-1];

    int nx_model = static_cast<int>((probhi[AMREX_SPACEDIM-1] -
                                     problo[AMREX_SPACEDIM-1] + 1.e-8_rt) / dx_model);

    int ng = 4;

    // now generate the initial models

    model_params.dens_base = problem::dens_base;
    model_params.T_star = problem::T_star;
    model_params.T_hi = problem::T_hi;
    model_params.T_lo = problem::T_lo;

    model_params.H_star = problem::H_star;
    model_params.atm_delta = problem::atm_delta;

    model_params.low_density_cutoff = problem::low_density_cutoff;

    generate_initial_model(nx_model + ng,
                           problo[AMREX_SPACEDIM-1] - ng * dx_model,
                           probhi[AMREX_SPACEDIM-1],
                           model_params, 0);

    // now create a perturbed model -- we want the same base conditions
    // a hotter temperature

    model_params.T_hi = model_params.T_hi + problem::dtemp;

    generate_initial_model(nx_model + ng,
                           problo[AMREX_SPACEDIM-1] - ng * dx_model,
                           probhi[AMREX_SPACEDIM-1],
                           model_params, 1);

}
#endif
