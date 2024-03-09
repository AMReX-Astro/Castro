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

    const int nsub = 16;

    Real dx_sub = dx[0] / static_cast<Real>(nsub);

#if AMREX_SPACEDIM >= 2
    Real dy_sub = dx[1] / static_cast<Real>(nsub);
#else
    Real dy_sub = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    Real dz_sub = dx[2] / static_cast<Real>(nsub);
#else
    Real dz_sub = 0.0_rt;
#endif

    Real r2init = (problem::r1 * probhi[0]) * (problem::r1 * probhi[0]);

    Real xcl = problo[0] + dx[0] * (static_cast<Real>(i));

#if AMREX_SPACEDIM >= 2
    Real ycl = problo[1] + dx[1] * (static_cast<Real>(j));
#else
    Real ycl = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    Real zcl = problo[2] + dx[2] * (static_cast<Real>(k));
#else
    Real zcl = 0.0_rt;
#endif

    Real vol_pert    = 0.0_rt;
    Real vol_ambient = 0.0_rt;

    for (int kk = 0; kk <= nsub - 1; ++kk) {
        Real zl = zcl + (static_cast<Real>(kk)         ) * dz_sub;
        Real zr = zcl + (static_cast<Real>(kk) + 1.0_rt) * dz_sub;
        Real zz = 0.5_rt * (zl + zr);

        for (int jj = 0; jj <= nsub - 1; ++jj) {
            Real yl = ycl + (static_cast<Real>(jj)         ) * dy_sub;
            Real yr = ycl + (static_cast<Real>(jj) + 1.0_rt) * dy_sub;
            Real yy = 0.5_rt * (yl + yr);

            for (int ii = 0; ii <= nsub - 1; ++ii) {
                Real xl = xcl + (static_cast<Real>(ii)         ) * dx_sub;
                Real xr = xcl + (static_cast<Real>(ii) + 1.0_rt) * dx_sub;
                Real xx = 0.5_rt * (xl + xr);

                Real rr2 = xx * xx + yy * yy + zz * zz;

                // Note: we assume spherical in 1D and Cartesian in 2D.

#if AMREX_SPACEDIM == 1
                // the volume of a subzone is (4/3) pi (xr^3 - xl^3).
                // we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)
                // The (4/3) pi dr factor is common, so we can neglect it.
                if (rr2 <= r2init) {
                    vol_pert += (xr * xr + xl * xr + xl * xl);
                }
                else {
                    vol_ambient += (xr * xr + xl * xr + xl * xl);
                }
#elif AMREX_SPACEDIM == 2
                if (rr2 <= r2init) {
                    vol_pert += xx; // common factor in volume ignored
                }
                else {
                    vol_ambient += xx;
                }
#else
                if (rr2 <= r2init) {
                    vol_pert = vol_pert + 1.0_rt; // common factor in volume ignored
                }
                else {
                    vol_ambient = vol_ambient + 1.0_rt;
                }
#endif

            }
        }
    }

    Real T_zone = std::pow((vol_pert * std::pow(problem::T1, 4) +
                            vol_ambient * std::pow(problem::T0, 4)) /
                           (vol_pert + vol_ambient), 0.25_rt);

    Real rho_zone = (vol_pert * problem::rho1 +
                     vol_ambient * problem::rho0) / (vol_pert + vol_ambient);

    state(i,j,k,URHO) = rho_zone;
    state(i,j,k,UTEMP) = T_zone;

    // set the composition to be all in the first species
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // set the internal energy via the EOS
    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T   = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = state(i,j,k,UFS+n) / state(i,j,k,URHO);
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = rho_zone * eos_state.e;
    state(i,j,k,UEDEN) = rho_zone * eos_state.e;
}
#endif
