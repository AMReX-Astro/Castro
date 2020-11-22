#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{

    // get the species indices

    problem::ihe4 = network_spec_index("helium-4");
    problem::ic12 = network_spec_index("carbon-12");
    problem::io16 = network_spec_index("oxygen-16");

    if (problem::ihe4 < 0 || problem::ic12 < 0 || problem::io16 < 0) {
        amrex::Error("ERROR: species indices not found");
    }

    // make sure that the carbon fraction falls between 0 and 1

    if (problem::cfrac > 1.0_rt || problem::cfrac < 0.0_rt) {
        amrex::Error("ERROR: cfrac must fall between 0 and 1");
    }

    // make sure that the oxygen fraction falls between 0 and 1

    if (problem::ofrac > 1.0_rt || problem::ofrac < 0.0_rt) {
        amrex::Error("ERROR: ofrac must fall between 0 and 1");
    }

    // make sure that the C/O fraction sums to no more than 1

    if (problem::cfrac + problem::ofrac > 1.0_rt) {
        amrex::Error("ERROR: cfrac + ofrac cannot exceed 1.");
    }

    // set the default mass fractions

    for (int n = 0; n < NumSpec; n++) {
        problem::xn[n] = problem::smallx;
    }
    problem::xn[problem::ic12] = amrex::max(problem::cfrac, problem::smallx);
    problem::xn[problem::io16] = amrex::max(problem::ofrac, problem::smallx);
    problem::xn[problem::ihe4] = 1.0_rt - problem::cfrac - problem::ofrac - (NumSpec - 2) * problem::smallx;

    // Set the ambient material

    problem::ambient_dens = problem::dens;
    for (int n = 0; n < NumSpec; n++) {
        problem::ambient_comp[n] = problem::xn[n];
    }

    eos_t eos_state;
    eos_state.rho = problem::ambient_dens;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n]  = problem::ambient_comp[n];
    }

    eos_state.T = problem::T_l;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_l = eos_state.e;

    eos_state.T = problem::T_r;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_r = eos_state.e;

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real width = problem::w_T * (probhi[0] - problo[0]);
    Real c_T = problo[0] + problem::center_T * (probhi[0] - problo[0]);

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    state(i,j,k,URHO) = problem::dens;

    Real sigma = 1.0_rt / (1.0_rt + std::exp(-(c_T - xcen)/ width));

    state(i,j,k,UTEMP) = problem::T_l + (problem::T_r - problem::T_l) * (1.0_rt - sigma);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * problem::xn[n];
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = problem::xn[n];
    }

#ifdef NSE_THERMO
    // set the aux quantities -- we need to do this if we are using the NSE network
    set_nse_aux_from_X(eos_state);

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * eos_state.aux[n];
    }
#endif

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = state(i,j,k,URHO) * (problem::vel - 2 * problem::vel * (1.0_rt - sigma));
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * state(i,j,k,URHO) * problem::vel * problem::vel;

}
#endif
