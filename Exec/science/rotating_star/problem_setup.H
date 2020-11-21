#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>

#ifdef NSE_THERMO
#include <nse.H>
#endif

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();


    // Read initial model

    read_model_file(problem::model_name);


#if AMREX_SPACEDIM == 1
    problem::center[0] = 0.0_rt;

#elif AMREX_SPACEDIM == 2
    //  assume axisymmetric
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);

#else
    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

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
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    Real dist = std::sqrt(x * x + y * y + z * z);

    state(i,j,k,URHO) = interpolate(dist, model::idens);
    state(i,j,k,UTEMP) = interpolate(dist, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(dist, model::ispec+n);
    }

    Real sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += state(i,j,k,UFS+n);
    }
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,UFS+n) / sumX;
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

#ifdef NSE_THERMO
    // set the aux quantities -- we need to do this if we are using the NSE network
    // this will fill eos_state.aux[]
    set_nse_aux_from_X(eos_state);

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = eos_state.aux[n];
    }

    if (in_nse(eos_state)) {
        Real abar;
        Real dq;
        Real dyedt;
        Real xn[NumSpec];

        nse_interp(state(i,j,k,UTEMP), state(i,j,k,URHO),
                   state(i,j,k,UFX+AuxZero::iye), abar, dq, dyedt, xn);

        state(i,j,k,UFX+AuxZero::iabar) = abar;
        state(i,j,k,UFX+AuxZero::ibea) = dq;
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = xn[n];
        }

        // we just got new X's, so we need to re-renormalize them
        sumX = 0.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            sumX += state(i,j,k,UFS+n);
        }
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = state(i,j,k,UFS+n) / sumX;
        }
    }
#endif

    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }
#ifdef NSE_THERMO
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state(i,j,k,UFX+n);
    }
#endif

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * state(i,j,k,UFX+n);
    }

    // initial velocities = 0
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

}
#endif
