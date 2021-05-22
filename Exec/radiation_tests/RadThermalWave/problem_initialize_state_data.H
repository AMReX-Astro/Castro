#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    const int nsub = 64;

    Real Vexp = 4.0_rt / 3.0_rt * M_PI * std::pow(problem::rexp, 3);

    Real rhoeexp = problem::Eexp / Vexp;

    Real rho = 1.0_rt;
    Real T = 1.0_rt;

    eos_t eos_state;
    eos_state.rho = rho;
    eos_state.T = T;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = 0.0_rt;
    }
    eos_state.xn[0]  = 1.0_rt;

    eos(eos_input_rt, eos_state);

    Real cv = eos_state.cv;

    rho = problem::rhocv / cv;

    Real eint = rhoeexp / rho;
    T = eint / cv;

    Real rhoe0 = rho * cv * problem::T0;

    Real dx_sub = dx[0] / static_cast<Real>(nsub);
    Real xcl = problo[0] + dx[0] * static_cast<Real>(i);
    Real xcr = xcl + dx[0];
    Real xcmin = amrex::min(std::abs(xcl), std::abs(xcr));
    Real xcmax = amrex::max(std::abs(xcl), std::abs(xcr));

#if AMREX_SPACEDIM >= 2
    Real dy_sub = dx[1] / static_cast<Real>(nsub);
    Real ycl = problo[1] + dx[1] * static_cast<Real>(j);
    Real ycr = ycl + dx[1];
    Real ycmin = amrex::min(std::abs(ycl), std::abs(ycr));
    Real ycmax = amrex::max(std::abs(ycl), std::abs(ycr));
#else
    Real dy_sub = 0.0_rt;
    Real ycl = 0.0_rt;
    Real ycmin = 0.0_rt;
    Real ycmax = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    Real dz_sub = dx[2] / static_cast<Real>(nsub);
    Real zcl = problo[2] + dx[2] * static_cast<Real>(k);
    Real zcr = zcl + dx[2];
    Real zcmin = amrex::min(std::abs(zcl), std::abs(zcr));
    Real zcmax = amrex::max(std::abs(zcl), std::abs(zcr));
#else
    Real dz_sub = 0.0_rt;
    Real zcl = 0.0_rt;
    Real zcmin = 0.0_rt;
    Real zcmax = 0.0_rt;
#endif


    Real rcmin = std::sqrt(xcmin * xcmin + ycmin * ycmin + zcmin * zcmin);
    Real rcmax = std::sqrt(xcmax * xcmax + ycmax * ycmax + zcmax * zcmax);

    state(i,j,k,URHO) = rho;
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    Real T_zone;
    Real rhoe_zone;

    if (rcmin >= problem::rexp) {
        T_zone = problem::T0;
        rhoe_zone = rhoe0;

    } else if (rcmax <= problem::rexp) {
        T_zone = T;
        rhoe_zone = rhoeexp;

    } else {
        Real vol_pert = 0.0_rt;
        Real vol_ambient = 0.0_rt;

        for (int kk = 0; kk <= dg2 * (nsub-1); kk++) {
            Real zz = zcl + (static_cast<Real>(kk) + 0.5_rt) * dz_sub;

            for (int jj = 0; jj <= dg1 * (nsub-1); jj++) {
                Real yy = ycl + (static_cast<Real>(jj) + 0.5_rt) * dy_sub;

                for (int ii = 0; ii < nsub; ii++) {
                    Real xx = xcl + (static_cast<Real>(ii) + 0.5_rt) * dx_sub;

                    Real rr2 = xx * xx + yy * yy + zz * zz;

                    if (coord_type == 1) {
                        // we are axisymmetric so weight by the radial coordinate
                        if (rr2 <= problem::rexp * problem::rexp) {
                            vol_pert += xx;
                        } else {
                            vol_ambient += xx;
                        }
                    } else {
                        //  Cartesian
                        if (rr2 <= problem::rexp * problem::rexp) {
                            vol_pert += 1.0_rt;
                        } else {
                            vol_ambient += 1.0_rt;
                        }
                    }

                }
            }
        }

        T_zone = (vol_pert * T + vol_ambient * problem::T0) / (vol_pert + vol_ambient);
        rhoe_zone = (vol_pert * rhoeexp + vol_ambient * rhoe0) / (vol_pert + vol_ambient);

    }

    state(i,j,k,UTEMP) = T_zone;
    state(i,j,k,UEDEN) = rhoe_zone;
    state(i,j,k,UEINT) = rhoe_zone;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);
}
#endif
