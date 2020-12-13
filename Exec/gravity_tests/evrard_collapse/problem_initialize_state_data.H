#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // Loop through the zones and set the zone state depending on whether we are
    // inside the sphere or if we are in an ambient zone.

    Real loc[3] = {0.0_rt};

    loc[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
    loc[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
    loc[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real radius = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        radius += (loc[d] - problem::center[d]) * (loc[d] - problem::center[d]);
    }
    radius = std::pow(radius, 0.5_rt);

    eos_t eos_state;
    if (radius <= problem::sphere_radius) {
        eos_state.rho = (problem::sphere_mass * C::M_solar) /
            (2.0_rt * M_PI * (problem::sphere_radius * problem::sphere_radius) * radius);
    } else {
        eos_state.rho = problem::ambient_density;
    }

    eos_state.e = 0.05_rt * C::Gconst * problem::sphere_mass * C::M_solar / radius;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = 1.0_rt / static_cast<Real>(NumSpec);
    }

    eos(eos_input_re, eos_state);

    state(i,j,k,URHO) = eos_state.rho;
    state(i,j,k,UTEMP) = eos_state.T;
    state(i,j,k,UEINT) = eos_state.e * eos_state.rho;
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = eos_state.xn[n] * eos_state.rho;
    }

    state(i,j,k,UMX) = state(i,j,k,URHO) * problem::smallu;
    state(i,j,k,UMY) = state(i,j,k,URHO) * problem::smallu;
    state(i,j,k,UMZ) = state(i,j,k,URHO) * problem::smallu;

    state(i,j,k,UEDEN) = state(i,j,k,UEINT) +
        (state(i,j,k,UMX) * state(i,j,k,UMX) +
         state(i,j,k,UMY) * state(i,j,k,UMY) +
         state(i,j,k,UMZ) * state(i,j,k,UMZ)) / (2.0_rt * state(i,j,k,URHO));
}
#endif
