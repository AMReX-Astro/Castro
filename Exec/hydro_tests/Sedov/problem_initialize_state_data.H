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

#if AMREX_SPACEDIM == 1
    int nsubx = problem::nsub;
    int nsuby = 1;
    int nsubz = 1;
#elif AMREX_SPACEDIM == 2
    int nsubx = problem::nsub;
    int nsuby = problem::nsub;
    int nsubz = 1;
#else
    int nsubx = problem::nsub;
    int nsuby = problem::nsub;
    int nsubz = problem::nsub;
#endif

    Real ds[3] = {0.0, 0.0, 0.0};
    Real xmin = 0.0;
    Real ymin = 0.0;
    Real zmin = 0.0;

    ds[0] = dx[0] / nsubx;
    xmin = problo[0] + dx[0] * static_cast<Real>(i);

#if AMREX_SPACEDIM >= 2
    ds[1] = dx[1] / nsuby;
    ymin = problo[1] + dx[1] * static_cast<Real>(j);
#endif

#if AMREX_SPACEDIM == 3
    ds[2] = dx[2] / nsubz;
    zmin = problo[2] + dx[2] * static_cast<Real>(k);
#endif

    Real vol_pert = 0.0_rt;
    Real vol_ambient = 0.0_rt;

    for (int kk = 0; kk <= nsubz-1; ++kk) {
        Real zz = zmin + ds[2] * (static_cast<Real>(kk) + 0.5_rt);

        for (int jj = 0; jj <= nsuby-1; ++jj) {
            Real yy = ymin + ds[1] * (static_cast<Real>(jj) + 0.5_rt);

            for (int ii = 0; ii <= nsubx-1; ++ii) {
                Real xx = xmin + ds[0] * (static_cast<Real>(ii) + 0.5_rt);

                // This is distance squared away from the blast wave center.
                Real dist = (problem::center[0] - xx) * (problem::center[0] - xx) +
                            (problem::center[1] - yy) * (problem::center[1] - yy) +
                            (problem::center[2] - zz) * (problem::center[2] - zz);
#if AMREX_SPACEDIM ==2
                if (coord_type == 2) {
                    Real dist = xx*xx - 2.0*xx*problem::center[0]*std::cos(yy - problem::center[1]) +
                        problem::center[0]*problem::center[0];
                }
#endif
                bool pert = false;

                if (dist <= problem::r_init * problem::r_init) {
                    pert = true;
                }

                if (coord_type == 1) {

                    // The volume of a cell is a annular cylindrical region.
                    // The main thing that matters is the distance from the
                    // symmetry axis.
                    //   V = pi*dy*(x_r**2 - x_l**2) = pi*dy*dx*0.5*xx
                    // (where x_r is the coordinate of the x right edge,
                    //        x_l is the coordinate of the x left edge,
                    //    and xx  is the coordinate of the x center of the cell)
                    //
                    // since dx and dy are constant, they cancel out

                    if (pert) {
                        vol_pert = vol_pert + xx;
                    }
                    else {
                        vol_ambient = vol_ambient + xx;
                    }

                }
                else if (coord_type == 2) {

                    Real xl = xx - 0.5_rt * ds[0];
                    Real xr = xx + 0.5_rt * ds[0];

                    // the volume of a subzone is (4/3) pi (xr^3 - xl^3).
                    // we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)
                    // The (4/3) pi dr factor is common, so we can neglect it.
#if AMREX_SPACEDIM == 1
                    Real vol = xr * xr + xl * xr + xl * xl;
#elif AMREX_SPACEDIM == 2
                    Real yl = yy - 0.5_rt * ds[1];
                    Real yr = yy + 0.5_rt * ds[1];
                    Real vol = (std::cos(yr) - std::cos(yl)) * (xr * xr + xl * xr + xl * xl);
#endif
                    if (pert) {
                        vol_pert = vol_pert + vol;
                    }
                    else {
                        vol_ambient = vol_ambient + vol;
                    }

                }
                else {

                    // Cartesian -- equal volume of dx * dy * dz is factored out.

                    if (pert) {
                        vol_pert = vol_pert + 1.0_rt;
                    }
                    else {
                        vol_ambient = vol_ambient + 1.0_rt;
                    }

                }

            }
        }
    }

    Real e_zone = (vol_pert * problem::e_exp + vol_ambient * problem::e_ambient) / (vol_pert + vol_ambient);

    Real eint = problem::dens_ambient * e_zone;

    state(i,j,k,URHO) = problem::dens_ambient;
    state(i,j,k,UMX) = 0.e0_rt;
    state(i,j,k,UMY) = 0.e0_rt;
    state(i,j,k,UMZ) = 0.e0_rt;

    state(i,j,k,UEDEN) = eint +
                         0.5e0_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) / state(i,j,k,URHO) +
                                     state(i,j,k,UMY) * state(i,j,k,UMY) / state(i,j,k,URHO) +
                                     state(i,j,k,UMZ) * state(i,j,k,UMZ) / state(i,j,k,URHO));

    state(i,j,k,UEINT) = eint;

    // The temperature initialization will be done later
    // in the call to clean_state. We want to avoid
    // EOS calls in ca_initdata. So that there is valid
    // data for an initial guess for the EOS, we initialize
    // to the ambient temperature everywhere.

    state(i,j,k,UTEMP) = problem::temp_ambient;

    state(i,j,k,UFS) = state(i,j,k,URHO);
}

#endif
