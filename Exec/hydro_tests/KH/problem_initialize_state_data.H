#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();


    Real vel1 = -0.5_rt;
    Real vel2 =  0.5_rt;

    int sine_n;
    Real w0;
    Real sigma;
    Real delta_y;

    if (problem::problem == 1) {
        sine_n = 4;
        w0 = 0.1_rt;
        sigma = 0.05_rt / std::sqrt(2.0_rt);

    } else if (problem::problem == 2) {
        sine_n = 2;
        w0 = 0.1_rt;
        delta_y = 0.05_rt;

    } else if (problem::problem == 3) {
        sine_n = 4;
        w0 = 0.01_rt;
        delta_y = 0.025_rt;

    } else if (problem::problem == 4) {
        sine_n = 2;
        w0 = 0.01_rt;
        delta_y = 0.025_rt;

    } else if (problem::problem == 5) {
        sine_n = 2;
        w0 = 0.01_rt;
        delta_y = 0.05_rt;
        sigma = 0.2_rt;
        vel1 = 1.0_rt;
        vel2 = 1.0_rt;
    }

    Real y1 = problem::center[1] - (probhi[1] - problo[1]) * 0.25e0_rt;
    Real y2 = problem::center[1] + (probhi[1] - problo[1]) * 0.25e0_rt;

    Real velz = 0.0;

#if AMREX_SPACEDIM == 3
    Real zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#else
    Real zz = 0.0;
#endif
#if AMREX_SPACEDIM == 2
    Real yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#else
    Real yy = 0.0;
#endif

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    // Assume the initial y-velocity represents the bulk flow
    // which will be perturbed in the following step

    Real vely = problem::bulk_velocity;
    Real dye = 0.0_rt;
    Real dens{};
    Real velx{};

    if (problem::problem == 1) {

        if (std::abs(yy - 0.5_rt * (y1 + y2)) < 0.5_rt * (y2 - y1)) {
            dens = problem::rho2;
            velx = 0.5_rt;
        } else {
            dens = problem::rho1;
            velx = -0.5_rt;
        }

        vely += w0 * std::sin(sine_n * M_PI * xx) *
            (std::exp(-(yy - y1) * (yy - y1) / (2 * sigma * sigma)) +
             std::exp(-(yy - y2) * (yy - y2) / (2 * sigma * sigma)));

    } else if (problem::problem == 2) {

        Real ramp = 1.0_rt / ((1.0_rt + std::exp(-2.0_rt * (yy - y1) / delta_y)) *
                              (1.0_rt + std::exp(2.0_rt * (yy - y2) / delta_y)));

        dens = problem::rho1 + ramp * (problem::rho2 - problem::rho1);
        velx = vel1 + ramp * (vel2 - vel1);

        vely += w0 * std::sin(sine_n * M_PI * xx);

    } else if (problem::problem == 3 || problem::problem == 4) {

        if (yy < y1 ) {
            dens = problem::rho1 - (problem::rho1 - problem::rho2) / 2 * std::exp((yy - y1) / delta_y);
            velx = vel1 - (vel1 - vel2) / 2 * std::exp((yy - y1) / delta_y);
        } else if (yy <= 0.5_rt * (y1 + y2) ) {
            dens = problem::rho2 + (problem::rho1 - problem::rho2) / 2 * std::exp((y1 - yy) / delta_y);
            velx = vel2 + (vel1 - vel2) / 2 * std::exp((y1 - yy) / delta_y);
        } else if (yy < y2) {
            dens = problem::rho2 + (problem::rho1 - problem::rho2) / 2 * std::exp((yy - y2) / delta_y);
            velx = vel2 + (vel1 - vel2) / 2 * std::exp((yy - y2) / delta_y);
        } else {
            dens = problem::rho1 - (problem::rho1 - problem::rho2) / 2 * std::exp((y2 - yy) / delta_y);
            velx = vel1 - (vel1 - vel2) / 2 * std::exp((y2 - yy) / delta_y);
        }

        vely += w0 * std::sin(sine_n * M_PI * xx);

    } else if (problem::problem == 5) {

        dens = problem::rho1 + (problem::rho2 - problem::rho1) * 0.5_rt * (std::tanh((yy - y1) / delta_y ) - std::tanh((yy - y2) / delta_y));
        velx = vel1 * (std::tanh((yy - y1) / delta_y) - std::tanh((yy - y2) / delta_y ) - 1.0_rt);
        vely += w0 * std::sin(sine_n * M_PI * xx) * (std::exp(-(yy - y1) * (yy - y1) / (sigma * sigma)) + std::exp(-(yy - y2) * (yy - y2) / (sigma * sigma)));
        dye  = 0.5_rt * (std::tanh((yy - y2) / delta_y) - std::tanh((yy - y1) / delta_y ) + 2.0_rt);

    } else {
#ifndef AMREX_USE_GPU
        amrex::Error("Error: This problem choice is undefined.");
#endif
    }

    state(i,j,k,URHO) = dens;
    state(i,j,k,UMX)  = dens * velx;
    state(i,j,k,UMY)  = dens * vely;
    state(i,j,k,UMZ)  = dens * velz;
    state(i,j,k,UFA)  = dye;

    // Establish the thermodynamic quantities
    eos_t eos_state;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 1.0_rt / static_cast<Real>(NumSpec);
        eos_state.xn[n] = 1.0_rt / static_cast<Real>(NumSpec);
    }

    eos_state.rho = state(i,j,k,URHO);
    eos_state.p = problem::pressure;

    eos(eos_input_rp, eos_state);

    state(i,j,k,UTEMP) = eos_state.T;

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) +
        0.5_rt * dens * (velx * velx + vely * vely + velz * velz);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }
}
#endif
