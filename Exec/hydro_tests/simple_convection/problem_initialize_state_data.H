#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata,
                                    amrex::RandomEngine const& engine)
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

#if AMREX_SPACEDIM == 2
    Real height = y;
#else
    Real height = z;
#endif

    state(i,j,k,URHO) = interpolate(height, model::idens);
    amrex::Real pressure = interpolate(height, model::ipres);
    state(i,j,k,UTEMP) = interpolate(height, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(height, model::ispec+n);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.p = pressure;
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rp, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

    // Initial velocities = 0
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    if (state(i,j,k,URHO) > 2.0_rt * problem::low_density_cutoff) {

        amrex::Real cs = std::sqrt(eos_rp::eos_gamma * eos_state.p / state(i,j,k,URHO));

        state(i,j,k,UMX) = state(i,j,k,URHO) * (2.0 * amrex::Random(engine) - 1.0_rt);
        state(i,j,k,UMY) = state(i,j,k,URHO) * (2.0 * amrex::Random(engine) - 1.0_rt);
        state(i,j,k,UMZ) = 0.0_rt;

        state(i,j,k,UEDEN) += 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                        state(i,j,k,UMY) * state(i,j,k,UMY) +
                                        state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
    }

}

#endif
