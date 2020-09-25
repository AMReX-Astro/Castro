#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        if (coord_type == 1 || coord_type == 2) {
            center[n] = 0.0_rt;
        }
        else {
            center[n] = 0.5_rt * (problo[n] + probhi[n]);
        }
    }

    Real xn_zone[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        xn_zone[n] = 0.0_rt;
    }
    xn_zone[0] = 1.0_rt;

    eos_t eos_state;

    // override the pressure with the temperature
    if (temp_ambient > 0.0_rt) {

        eos_state.rho = dens_ambient;
        eos_state.T = temp_ambient;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn_zone[n];
        }

        eos(eos_input_rt, eos_state);

     p_ambient = eos_state.p;

    }

    // Calculate ambient state data

    eos_state.rho = dens_ambient;
    eos_state.p   = p_ambient;
    eos_state.T   = 1.e9_rt; // Initial guess for iterations
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_zone[n];
    }
        
    eos(eos_input_rp, eos_state);

    e_ambient = eos_state.e;
    temp_ambient = eos_state.T;

    // set explosion pressure -- we will convert the point-explosion energy into
    // a corresponding pressure distributed throughout the perturbed volume

    Real vctr = 0.0;

    if (coord_type == 0) {

#if AMREX_SPACEDIM == 1

        amrex::Abort("Sedov problem unsupported in 1D Cartesian geometry.");

#elif AMREX_SPACEDIM == 2

        // Cylindrical problem in Cartesian coordinates

        vctr = M_PI * r_init * r_init;

#else

        // Spherical problem in Cartesian coordinates

        vctr = (4.0_rt / 3.0_rt) * M_PI * r_init * r_init * r_init;

#endif

    }
    else if (coord_type == 1) {

#if AMREX_SPACEDIM == 1

        vctr = M_PI * r_init * r_init;

#elif AMREX_SPACEDIM == 2

        vctr = (4.0_rt / 3.0_rt) * M_PI * r_init * r_init * r_init;

#else

        amrex::Abort("Sedov problem unsupported in 3D axisymmetric geometry.");

#endif

    }
    else if (coord_type == 2) {

        // Must have AMREX_SPACEDIM == 1 for this coord_type.

        vctr = (4.0_rt / 3.0_rt) * M_PI * r_init * r_init * r_init;

    }

    e_exp = exp_energy / vctr / dens_ambient;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    
#if AMREX_SPACEDIM == 1
    int nsubx = nsub;
    int nsuby = 1;
    int nsubz = 1;
#elif AMREX_SPACEDIM == 2
    int nsubx = nsub;
    int nsuby = nsub;
    int nsubz = 1;
#else
    int nsubx = nsub;
    int nsuby = nsub;
    int nsubz = nsub;
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

                Real dist = (center[0] - xx) * (center[0] - xx) +
                            (center[1] - yy) * (center[1] - yy) +
                            (center[2] - zz) * (center[2] - zz);

                bool pert = false;

                if (dist <= r_init * r_init) {
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

                    if (pert) {
                        vol_pert = vol_pert + (xr * xr + xl * xr + xl * xl);
                    }
                    else {
                        vol_ambient = vol_ambient + (xr * xr + xl * xr + xl * xl);
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

    Real e_zone = (vol_pert * e_exp + vol_ambient * e_ambient) / (vol_pert + vol_ambient);

    Real eint = dens_ambient * e_zone;

    state(i,j,k,URHO) = dens_ambient;
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

    state(i,j,k,UTEMP) = temp_ambient;

    state(i,j,k,UFS) = state(i,j,k,URHO);
}

#endif
