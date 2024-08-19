#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <global.H>
#include <network.H>
#include <model_parser.H>
#include <initial_model.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // check to make sure that small_dens is less than low_density_cutoff
    // if not, funny things can happen above the atmosphere

    if (small_dens >= 0.99_rt * problem::low_density_cutoff) {
        amrex::Error("ERROR: small_dens should be set lower than low_density_cutoff");
    }

    // make sure hse_fixed_temp is the same as T_star, if it's specified

    if (hse_fixed_temp > 0.0_rt && hse_fixed_temp != problem::T_star) {
        amrex::Error("ERROR: hse_fixed_temp should be the same as T_star");
    }

    // get the species indices

    bool species_defined = true;

    int ifuel1 = network_spec_index(problem::fuel1_name);
    if (ifuel1 < 0) {
        species_defined = false;
    }

    int ifuel2;
    if (!problem::fuel2_name.empty()) {
        ifuel2 = network_spec_index(problem::fuel2_name);
        if (ifuel2 < 0) {
            species_defined = false;
        }
    }

    int ifuel3;
    if (!problem::fuel3_name.empty()) {
        ifuel3 = network_spec_index(problem::fuel3_name);
        if (ifuel3 < 0) {
            species_defined = false;
        }
    }

    int ifuel4;
    if (!problem::fuel4_name.empty()) {
        ifuel4 = network_spec_index(problem::fuel4_name);
        if (ifuel4 < 0) {
            species_defined = false;
        }
    }

    int iash1 = network_spec_index(problem::ash1_name);
    if (iash1 < 0) {
        species_defined = false;
    }

    int iash2;
    if (!problem::ash2_name.empty()) {
        iash2 = network_spec_index(problem::ash2_name);
        if (iash2 < 0) {
            species_defined = false;
        }
    }

    int iash3;
    if (!problem::ash3_name.empty()) {
        iash3 = network_spec_index(problem::ash3_name);
        if (iash3 < 0) {
            species_defined = false;
        }
    }

    if (! species_defined) {
        std::cout << ifuel1 << " " << ifuel2 << " " << ifuel3 << " " << ifuel4 << std::endl;
        std::cout << iash1 << " " << iash2 << " "<< iash3 << std::endl;
        amrex::Error("ERROR: species not defined");
    }

    model_t model_params;

    // set the composition of the underlying star


    for (Real &X : model_params.xn_star) {
        X = problem::smallx;
    }
    model_params.xn_star[iash1] = problem::ash1_frac;
    if (!problem::ash2_name.empty()) {
        model_params.xn_star[iash2] = problem::ash2_frac;
    }
    if (!problem::ash3_name.empty()) {
        model_params.xn_star[iash3] = problem::ash3_frac;
    }

    // and the composition of the accreted layer

    for (Real &X : model_params.xn_base) {
        X = problem::smallx;
    }
    model_params.xn_base[ifuel1] = problem::fuel1_frac;
    if (!problem::fuel2_name.empty()) {
        model_params.xn_base[ifuel2] = problem::fuel2_frac;
    }
    if (!problem::fuel3_name.empty()) {
        model_params.xn_base[ifuel3] = problem::fuel3_frac;
    }
    if (!problem::fuel4_name.empty()) {
        model_params.xn_base[ifuel4] = problem::fuel4_frac;
    }

    // check if they sum to 1

    Real sumX = 0.0_rt;
    for (Real X : model_params.xn_star) {
        sumX += X;
    }
    if (std::abs(sumX - 1.0_rt) > NumSpec * problem::smallx) {
        amrex::Error("ERROR: ash mass fractions don't sum to 1");
    }

    sumX = 0.0_rt;
    for (Real X : model_params.xn_base) {
        sumX += X;
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

    int nx_model = static_cast<int>((probhi[AMREX_SPACEDIM-1] - problo[AMREX_SPACEDIM-1]) /
                                    dx_model);

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

    // set center

    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        problem::center[d] = 0.5_rt * (problo[d] + probhi[d]);
    }

#if AMREX_SPACEDIM == 2
    // for axisymmetry, put the x-center on the x-axis
    // and the y-center at 0, so the height computation is okay
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.0_rt;
#endif

    // set the ambient state for the upper boundary condition

    ambient::ambient_state[URHO] = model::profile(0).state(model::npts-1, model::idens);
    ambient::ambient_state[UTEMP] = model::profile(0).state(model::npts-1, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        ambient::ambient_state[UFS+n] =
            ambient::ambient_state[URHO] * model::profile(0).state(model::npts-1, model::ispec+n);
    }

    ambient::ambient_state[UMX] = 0.0_rt;
    ambient::ambient_state[UMY] = 0.0_rt;
    ambient::ambient_state[UMZ] = 0.0_rt;

    // make the ambient state thermodynamically consistent

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / eos_state.rho;
    }

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = eos_state.rho * eos_state.e;
    ambient::ambient_state[UEDEN] = eos_state.rho * eos_state.e;
}

#endif
