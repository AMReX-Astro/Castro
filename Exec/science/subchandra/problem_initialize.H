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

    // set center variable

#if AMREX_SPACEDIM == 2
    // axisymmetric
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.0_rt;

#elif AMREX_SPACEDIM == 3
    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

    // we are going to generate an initial model that is half the width of the
    // domain (or full width for axisymmetry)

    // we use the fine grid dx for the model resolution
    auto fine_geom = global::the_amr_ptr->Geom(global::the_amr_ptr->maxLevel());

    // we will use a factor of 1/2 in dx to give a bit more resolution
    // when mapping
    auto dx = fine_geom.CellSizeArray();
    auto dx_model = std::min(problem::model_min_res, 0.5_rt * dx[AMREX_SPACEDIM-1]);

    int nx_model = static_cast<int>(problem::model_r_max / dx_model);

    model_t model_params;

    // first the composition

    // get the species indices

    int ihe4 = network_spec_index("helium-4");
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");
    int ine20 = network_spec_index("neon-20");

    int in14 = network_spec_index("nitrogen-14");

    if (ihe4 < 0 || ic12 < 0 || io16 < 0 || ine20 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    if (in14 < 0 && problem::X_N14 > 0.0_rt) {
        amrex::Error("ERROR: N14 not defined");
    }

    // core mass fractions should sum to 1
    if (std::abs(problem::core_XC12 +
                 problem::core_XO16 +
                 problem::core_XNe20 - 1.0_rt) > 1.e-10) {
        std::cout << problem::core_XC12 << " " << problem::core_XO16 << " " << problem::core_XNe20 << std::endl;
        amrex::Error("ERROR: core mass fractions should sum to 1");
    }

    for (auto& X : model_params.xn_core) {
        X = network_rp::small_x;
    }
    for (auto& X : model_params.xn_he) {
        X = network_rp::small_x;
    }

    int ncore = 0;
    if (problem::core_XC12 > 0.0) {
        ncore++;
    }

    if (problem::core_XO16 > 0.0) {
        ncore++;
    }

    if (problem::core_XNe20 > 0.0) {
        ncore++;
    }

    AMREX_ALWAYS_ASSERT(ncore > 0);

    // we need to initialize the core to be 1 - other_xn * small_x.
    // We'll distribute this remainder over the core nuclei

    int other_xn = NumSpec - ncore;

    model_params.xn_core[ic12] = problem::core_XC12 - other_xn * network_rp::small_x / ncore;
    model_params.xn_core[io16] = problem::core_XO16 - other_xn * network_rp::small_x / ncore;;
    model_params.xn_core[ine20] = problem::core_XNe20 - other_xn * network_rp::small_x / ncore;;

    if (problem::X_N14 > 0.0_rt) {
        model_params.xn_he[in14] = problem::X_N14;
    }

    if (problem::X_C12 > 0.0_rt) {
        model_params.xn_he[ic12] = problem::X_C12;
    }

    if (problem::X_O16 > 0.0_rt) {
        model_params.xn_he[io16] = problem::X_O16;
    }

    model_params.xn_he[ihe4] = 1.0_rt -
        model_params.xn_he[ic12] -
        model_params.xn_he[in14] -
        model_params.xn_he[io16] -
        (NumSpec-4) * network_rp::small_x;

    model_params.ihe4 = ihe4;
    model_params.ic12 = ic12;
    model_params.in14 = in14;
    model_params.io16 = io16;
    model_params.ine20 = ine20;

    model_params.X_N14 = problem::X_N14;
    model_params.X_C12 = problem::X_C12;
    model_params.X_O16 = problem::X_O16;

    // convert the envelope and WD mass into CGS

    model_params.M_WD = problem::M_WD * C::M_solar;
    model_params.M_He = problem::M_He * C::M_solar;

    model_params.delta = problem::delta;
    model_params.T_core = problem::temp_core;
    model_params.T_base = problem::temp_base;
    model_params.low_density_cutoff = problem::low_density_cutoff;
    model_params.T_fluff = problem::temp_fluff;

    model_params.isothermal_layer = problem::isothermal_layer;

    model_params.tol_hse = problem::tol_hse;
    model_params.tol_WD = problem::tol_WD;
    model_params.tol_He = problem::tol_He;

    // generate the model

    generate_initial_model(nx_model,
                           problem::model_r_max,
                           model_params);

    // find the distance (starting from the center) where the He layer begins

    problem::R_He_base = 0.0;
    problem::ihe4 = ihe4;

    for (int n = 0; n < model::npts; n++) {
        if (model::profile(0).state(n, model::ispec+problem::ihe4) > 0.5_rt) {
            problem::R_He_base = model::profile(0).r(n);
            break;
        }
    }

    amrex::Print() << Font::Bold << FGColor::Green
                   << "Base of He layer found at r = " << problem::R_He_base << " cm"
                   << ResetDisplay << std::endl;
}

#endif
