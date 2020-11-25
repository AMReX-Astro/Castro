#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
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

    int ifuel4;
    if (problem::fuel4_name != "") {
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
        std::cout << ifuel1 << " " << ifuel2 << " " << ifuel3 << " " << ifuel4 << std::endl;
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
    if (problem::fuel4_name != "") {
        model_params.xn_base[ifuel4] = problem::fuel4_frac;
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

    int nx_model = static_cast<int>((probhi[AMREX_SPACEDIM-1] - problo[AMREX_SPACEDIM-1]) /
                                    problem::dx_model);

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
                           problo[AMREX_SPACEDIM-1] - ng * problem::dx_model,
                           probhi[AMREX_SPACEDIM-1],
                           model_params, 0);

    // now create a perturbed model -- we want the same base conditions
    // a hotter temperature

    model_params.T_hi = model_params.T_hi + problem::dtemp;

    generate_initial_model(nx_model + ng,
                           problo[AMREX_SPACEDIM-1] - ng * problem::dx_model,
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
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    // lateral distance and height

    Real r;
    Real height;

    if (AMREX_SPACEDIM == 1) {
        r = 1.0_rt;
        height = x;

    } else if (AMREX_SPACEDIM == 2) {
        r = x;
        height = y;

    } else if (AMREX_SPACEDIM == 3) {
        r = std::sqrt(x * x + y * y);
        height = z;
    }

    // blending factor

    Real f;

    if (r < problem::x_half_max) {
        f = 1.0_rt;

    } else if (r > problem::x_half_max + problem::x_half_width) {
        f = 0.0_rt;

    } else {
        f = -(r - problem::x_half_max) / problem::x_half_width + 1.0_rt;
    }

    state(i,j,k,URHO) = f * interpolate(height, model::idens, 1) +
             (1.0_rt - f) * interpolate(height, model::idens, 0);

    state(i,j,k,UTEMP) = f * interpolate(height, model::itemp, 1) +
              (1.0_rt - f) * interpolate(height, model::itemp, 0);

    Real temppres = f * interpolate(height, model::ipres, 1) +
         (1.0_rt - f) * interpolate(height, model::ipres, 0);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = f * interpolate(height, model::ispec+n, 1) +
                  (1.0_rt - f) * interpolate(height, model::ispec+n, 0);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    eos_state.p = temppres;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rp, eos_state);

    state(i,j,k,UTEMP) = eos_state.T;
    state(i,j,k,UEINT) = eos_state.rho * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);

    // Initial velocities = 0

    state(i,j,k,UMX) = 0.e0_rt;
    state(i,j,k,UMY) = 0.e0_rt;
    state(i,j,k,UMZ) = 0.e0_rt;

    // convert to partial densities

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }
}
#endif
