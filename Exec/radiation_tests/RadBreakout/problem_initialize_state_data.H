#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    constexpr Real Tindex = 0.5_rt;
    constexpr int nsub = 16;

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real dx_sub = dx[0] / static_cast<Real>(nsub);

    Real rhowind0 = problem::rhowind1 *
        std::pow(problem::rwind1/problem::rwind0, 2);
    Real rlast = model::profile(0).r(model::npts-1);
    Real rholast = model::profile(0).state(model::npts-1, model::idens);
    Real rbase = rlast * problem::rbasefac;

    Real Twind0 = problem::Twind1 * std::pow(problem::rwind1 / problem::rwind0, Tindex);
    Real Tlast = model::profile(0).state(model::npts-1, model::itemp);

    Real xcl = problo[0] + dx[0] * static_cast<Real>(i);

    Real vtot = 0.0_rt;
    Real rho = 0.0_rt;
    Real T = 0.0_rt;
    Real u = 0.0_rt;
    Real Ye = 0.0_rt;
    Real Abar = 0.0_rt;

    for (int ii = 0; ii < nsub; ii++) {
        Real xx = xcl + (static_cast<Real>(ii) + 0.5_rt) * dx_sub;
        Real vsub = xx * xx;
        vtot += vsub;

        if (xx >= model::profile(0).r(model::npts-1)) {
            Real rho_tmp;
            if (xx >= problem::rwind0 ) {
                rho_tmp = problem::rhowind1 * (problem::rwind1 / xx) * (problem::rwind1 / xx);
            } else {
                rho_tmp = rholast * std::pow(rhowind0 / rholast,
                                             ((std::log(xx - rbase) - std::log(rlast - rbase)) /
                                              (std::log(problem::rwind0 - rbase) - std::log(rlast - rbase))));
            }
            rho += rho_tmp * vsub;

            Real T_tmp;
            if (xx >= problem::rwind0) {
                T_tmp = problem::Twind1 * std::pow(problem::rwind1 / xx, Tindex);
            } else {
                T_tmp = Tlast * std::pow(Twind0 / Tlast,
                                         ((std::log(xx - rbase) - std::log(rlast - rbase)) /
                                          (std::log(problem::rwind0 - rbase) - std::log(rlast - rbase))));
            }
            T += vsub * T_tmp;

            Ye += model::profile(0).state(model::npts-1, model::iaux+iye) * vsub;
            Abar += model::profile(0).state(model::npts-1, model::iaux+iinvmu) * vsub;

            // u += 0.e0_rt;

        } else if (xx <= model::profile(0).r(0)) {
            rho += model::profile(0).state(0, model::idens) * vsub;
            T += model::profile(0).state(0, model::itemp) * vsub;
            // u += 0.e0_rt;
            Ye += model::profile(0).state(0, model::iaux+iye) * vsub;
            Abar += model::profile(0).state(0, model::iaux+iinvmu) * vsub;

        } else {
            rho += interpolate(xx, model::idens) * vsub;
            T += interpolate(xx, model::itemp) * vsub;
            u += interpolate(xx, model::ivelr) * vsub;
            Ye += interpolate(xx, model::iaux+iye) * vsub;
            Abar += interpolate(xx, model::iaux+iinvmu) * vsub;
        }
    }

    rho /= vtot;
    T /= vtot;
    u /= vtot;
    Ye /= vtot;
    Abar /= vtot;

    // what we read in from the model was really 1/mu_I, which is why we called it Abar here.
    // now we compute that total 1/mu = 1/mu_I + 1/mu_e = 1/Abar + Ye

    Real invmu = (1.0_rt + Abar * Ye) / Abar;

    state(i,j,k,URHO) = rho;
    state(i,j,k,UTEMP) = T;
    state(i,j,k,UMX) = rho * u;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // set the composition to be all in the first species

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);

    state(i,j,k,UFX+iye) = Ye * rho;
    state(i,j,k,UFX+iinvmu) = invmu * rho;

    // set the internal energy via the EOS

    Real rhoInv = 1.0_rt / state(i,j,k,URHO);

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n) * rhoInv;
    }
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state(i,j,k,UFX+n) * rhoInv;
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) +
        0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                  state(i,j,k,UMY) * state(i,j,k,UMY) +
                  state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
}
#endif
