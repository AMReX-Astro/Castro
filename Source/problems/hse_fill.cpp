#include <AMReX_BLFort.H>
#include <Castro.H>
#include <runtime_parameters.H>
#include <ext_bc_types.H>

using namespace amrex;


// a hydrostatic boundary conditions -- this relies on the assumption
// that the gravitation acceleration is constant


void
hse_fill(const Box& bx, Array4<Real> const& adv,
              Geometry const& geom, const Vector<BCRec>& bcr,
              const Real time)
{


    auto domlo = geom.Domain().loVect3d();
    auto domhi = geom.Domain().hiVect3d();

    auto problo = geom.ProbLoArray();

    auto lo = bx.loVect();
    auto hi = bx.hiVect();

    auto adv_bx = Box(adv);
    auto adv_lo = adv_bx.loVect3d();
    auto adv_hi = adv_bx.hiVect3d();

    auto dx = geom.CellSizeArray();

    //
    // x boundaries
    //

    // XLO

    if (bcr[URHO].lo(0) == EXT_DIR && lo[0] < domlo[0]) {

        if (xl_ext_bc_type == EXT_HSE) {

            // we need to integrate in the i direction from adv_lo[0]
            // to domlo[0]-1, but we want that to be handled by a
            // single thread on the GPU

            Box gbx(IntVect(D_DECL(domlo[0]-1, lo[1], lo[2])),
                    IntVect(D_DECL(domlo[0]-1, hi[1], hi[2])));

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {

                Real dens_above = adv(domlo[0],j,k,URHO);
                Real temp_above = adv(domlo[0],j,k,UTEMP);
                Real X_zone[NumSpec];
                for (int n = 0; n < NumSpec; n++) {
                    X_zone[n] = adv(domlo[0],j,k,UFS+n) / dens_above;
                }
#if NAUX_NET > 0
                Real aux_zone[NumAux];
                for (int n = 0; n < NumAux; n++) {
                    aux_zone[n] = adv(domlo[0],j,k,UFX+n) / dens_above;
                }
#endif

                //  keep track of the density at the base of the domain

                Real dens_base = dens_above;

                // get pressure in this zone (the initial above zone)

                eos_rep_t eos_state;
                eos_state.rho = dens_above;
                eos_state.T = temp_above;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = X_zone[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    eos_state.aux[n] = aux_zone[n];
                }
#endif

                eos(eos_input_rt, eos_state);

                Real pres_above = eos_state.p;

                for (int ii = domlo[0]-1; ii >= adv_lo[0]; ii--) {

                    // we are integrating along a column at constant i.
                    // Make sure that our starting state is well-defined

                    // HSE integration to get density, pressure

                    // initial guesses

                    Real dens_zone = dens_above;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(ii+1,j,k,UTEMP) - adv(ii+2,j,k,UTEMP);
                    } else {
                        temp_zone = temp_above;
                    }

                    bool converged_hse = false;

                    Real p_want;
                    Real drho;

                    for (int iter = 0; iter < hse::MAX_ITER; iter++) {

                        // pressure needed from HSE

                        p_want = pres_above -
                            dx[0] * 0.5_rt * (dens_zone + dens_above) * gravity::const_grav;

                        // pressure from EOS

                        eos_state.rho = dens_zone;
                        eos_state.T = temp_zone;
                        // xn is already set above

                        eos(eos_input_rt, eos_state);

                        Real pres_zone = eos_state.p;
                        Real dpdr = eos_state.dpdr;

                        // Newton-Raphson - we want to zero A = p_want - p(rho)
                        Real A = p_want - pres_zone;
                        drho = A / (dpdr + 0.5_rt * dx[0] * gravity::const_grav);

                        dens_zone = amrex::max(0.9_rt*dens_zone,
                                               amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

                        // convergence?

                        if (std::abs(drho) < hse::TOL * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                    }

#ifndef AMREX_USE_CUDA
                    if (! converged_hse) {
                        std::cout << "ii, j, k, domlo[0]: " << ii << " " << j << " " << k << " " << domlo[0] << std::endl;
                        std::cout << "p_want:    " << p_want << std::endl;
                        std::cout << "dens_zone: " << dens_zone << std::endl;
                        std::cout << "temp_zone: " << temp_zone << std::endl;
                        std::cout << "drho:      " << drho << std::endl;
                        std::cout << std::endl;
                        std::cout << "column info: " << std::endl;
                        std::cout << "   dens: " << adv(ii,j,k,URHO) << std::endl;
                        std::cout << "   temp: " << adv(ii,j,k,UTEMP) << std::endl;
                        amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in -X BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(ii,j,k,UMX) = 0.0_rt;
                       adv(ii,j,k,UMY) = 0.0_rt;
                       adv(ii,j,k,UMZ) = 0.0_rt;

                   } else {

                      if (hse_reflect_vels == 1) {
                          // reflect normal, zero gradient for transverse
                          // note: we need to match the corresponding
                          // zone on the other side of the interface
                          int ioff = domlo[0]-ii-1;
                          adv(ii,j,k,UMX) = -dens_zone * (adv(domlo[0]+ioff,j,k,UMX) / adv(domlo[0]+ioff,j,k,URHO));

                          adv(ii,j,k,UMY) = -dens_zone * (adv(domlo[0],j,k,UMY) / dens_base);
                          adv(ii,j,k,UMZ) = -dens_zone * (adv(domlo[0],j,k,UMZ) / dens_base);
                      } else {
                          // zero gradient
                          adv(ii,j,k,UMX) = dens_zone * (adv(domlo[0],j,k,UMX) / dens_base);
                          adv(ii,j,k,UMY) = dens_zone * (adv(domlo[0],j,k,UMY) / dens_base);
                          adv(ii,j,k,UMZ) = dens_zone * (adv(domlo[0],j,k,UMZ) / dens_base);
                      }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   // store the final state

                   adv(ii,j,k,URHO) = dens_zone;
                   adv(ii,j,k,UEINT) = dens_zone * eint;
                   adv(ii,j,k,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(ii,j,k,UMX) * adv(ii,j,k,UMX) +
                                 adv(ii,j,k,UMY) * adv(ii,j,k,UMY) +
                                 adv(ii,j,k,UMZ) * adv(ii,j,k,UMZ)) / dens_zone;
                   adv(ii,j,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(ii,j,k,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(ii,j,k,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_above = dens_zone;
                   pres_above = pres_zone;

                }
            });

        }

    }


    // XHI

    if (bcr[URHO].hi(0) == EXT_DIR && hi[0] > domhi[0]) {

       if (xr_ext_bc_type == EXT_HSE) {

            Box gbx(IntVect(D_DECL(domhi[0]+1, lo[1], lo[2])),
                    IntVect(D_DECL(domhi[0]+1, hi[1], hi[2])));

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {

                Real dens_below = adv(domhi[0],j,k,URHO);
                Real temp_below = adv(domhi[0],j,k,UTEMP);
                Real X_zone[NumSpec];
                for (int n = 0; n < NumSpec; n++) {
                    X_zone[n] = adv(domhi[0],j,k,UFS+n) / dens_below;
                }
#if NAUX_NET > 0
                Real aux_zone[NumAux];
                for (int n = 0; n < NumAux; n++) {
                    aux_zone[n] = adv(domhi[0],j,k,UFX+n) / dens_below;
                }
#endif

                // keep track of the density at the top of the domain

                Real dens_base = dens_below;

                // get pressure in this zone (the initial below zone)

                eos_rep_t eos_state;
                eos_state.rho = dens_below;
                eos_state.T = temp_below;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = X_zone[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    eos_state.aux[n] = aux_zone[n];
                }
#endif

                eos(eos_input_rt, eos_state);

                Real pres_below = eos_state.p;

                for (int ii = domhi[0]+1; ii <= adv_hi[0]; ii++) {

                    // HSE integration to get density, pressure

                    // initial guesses
                    Real dens_zone = dens_below;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(ii-1,j,k,UTEMP) - adv(ii-2,j,k,UTEMP);
                    } else {
                        temp_zone = temp_below;
                    }

                    bool converged_hse = false;

                    Real p_want;
                    Real drho;

                    for (int iter = 0; iter < hse::MAX_ITER; iter++) {

                        // pressure needed from HSE
                        p_want = pres_below +
                            dx[0] * 0.5_rt * (dens_zone + dens_below) * gravity::const_grav;

                        // pressure from EOS

                        eos_state.rho = dens_zone;
                        eos_state.T = temp_zone;
                        // xn is already set above

                        eos(eos_input_rt, eos_state);

                        Real pres_zone = eos_state.p;
                        Real dpdr = eos_state.dpdr;

                        // Newton-Raphson - we want to zero A = p_want - p(rho)
                        Real A = p_want - pres_zone;
                        drho = A / (dpdr - 0.5_rt * dx[0] * gravity::const_grav);

                        dens_zone = amrex::max(0.9_rt*dens_zone,
                                               amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

                        // convergence?

                        if (std::abs(drho) < hse::TOL * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                    }

#ifndef AMREX_USE_CUDA
                   if (! converged_hse) {
                       std::cout << "ii, j, k, domhi[0]: " << ii << " " << j << " " << k << " " << domhi[0] << std::endl;
                       std::cout << "p_want:    " << p_want << std::endl;
                       std::cout << "dens_zone: " << dens_zone << std::endl;
                       std::cout << "temp_zone: " << temp_zone << std::endl;
                       std::cout << "drho:      " << drho << std::endl;
                       std::cout << std::endl;
                       std::cout << "column info: " << std::endl;
                       std::cout << "   dens: " << adv(ii,j,k,URHO) << std::endl;
                       std::cout << "   temp: " << adv(ii,j,k,UTEMP) << std::endl;
                       amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in +X BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(ii,j,k,UMX) = 0.0_rt;
                       adv(ii,j,k,UMY) = 0.0_rt;
                       adv(ii,j,k,UMZ) = 0.0_rt;

                   } else {

                       if (hse_reflect_vels == 1) {
                           // reflect normal, zero gradient for transverse
                           // note: we need to match the corresponding
                           // zone on the other side of the interface
                           int ioff = ii-domhi[0]-1;
                           adv(ii,j,k,UMX) = -dens_zone * (adv(domhi[0]-ioff,j,k,UMX) / adv(domhi[0]-ioff,j,k,URHO));

                           adv(ii,j,k,UMY) = -dens_zone * (adv(domhi[0],j,k,UMY) / dens_base);
                           adv(ii,j,k,UMZ) = -dens_zone * (adv(domhi[0],j,k,UMZ) / dens_base);
                       } else {
                           // zero gradient
                           adv(ii,j,k,UMX) = dens_zone * (adv(domhi[0],j,k,UMX) / dens_base);
                           adv(ii,j,k,UMY) = dens_zone * (adv(domhi[0],j,k,UMY) / dens_base);
                           adv(ii,j,k,UMZ) = dens_zone * (adv(domhi[0],j,k,UMZ) / dens_base);
                       }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   //  store the final state

                   adv(ii,j,k,URHO) = dens_zone;
                   adv(ii,j,k,UEINT) = dens_zone * eint;
                   adv(ii,j,k,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(ii,j,k,UMX) * adv(ii,j,k,UMX) +
                                 adv(ii,j,k,UMY) * adv(ii,j,k,UMY) +
                                 adv(ii,j,k,UMZ) * adv(ii,j,k,UMZ)) / dens_zone;
                   adv(ii,j,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(ii,j,k,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(ii,j,k,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_below = dens_zone;
                   pres_below = pres_zone;

                }
            });

       }

    }


#if AMREX_SPACEDIM >= 2
    //
    // y boundaries
    //

    // YLO

    if (bcr[URHO].lo(1) == EXT_DIR && lo[1] < domlo[1]) {

        if (yl_ext_bc_type == EXT_HSE) {

            Box gbx(IntVect(D_DECL(lo[0], domlo[1]-1, lo[2])),
                    IntVect(D_DECL(hi[0], domlo[1]-1, hi[2])));

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {

                Real dens_above = adv(i,domlo[1],k,URHO);
                Real temp_above = adv(i,domlo[1],k,UTEMP);
                Real X_zone[NumSpec];
                for (int n = 0; n < NumSpec; n++) {
                    X_zone[n] = adv(i,domlo[1],k,UFS+n) / dens_above;
                }
#if NAUX_NET > 0
                Real aux_zone[NumAux];
                for (int n = 0; n < NumAux; n++) {
                    aux_zone[n] = adv(i,domlo[1],k,UFX+n) / dens_above;
                }
#endif

                // keep track of the density at the base of the domain

                Real dens_base = dens_above;

                // get pressure in this zone (the initial above zone)

                eos_rep_t eos_state;
                eos_state.rho = dens_above;
                eos_state.T = temp_above;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = X_zone[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    eos_state.aux[n] = aux_zone[n];
                }
#endif

                eos(eos_input_rt, eos_state);

                Real pres_above = eos_state.p;

                for (int jj = domlo[1]-1; jj >= adv_lo[1]; jj--) {

                    // HSE integration to get density, pressure

                    // initial guesses

                    Real dens_zone = dens_above;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(i,jj+1,k,UTEMP) - adv(i,jj+2,k,UTEMP);
                    } else {
                        temp_zone = temp_above;
                    }

                    bool converged_hse = false;

                    Real p_want;
                    Real drho;

                    for (int iter = 0; iter < hse::MAX_ITER; iter++) {

                        // pressure needed from HSE

                        p_want = pres_above -
                            dx[1] * 0.5_rt * (dens_zone + dens_above) * gravity::const_grav;

                        // pressure from EOS

                        eos_state.rho = dens_zone;
                        eos_state.T = temp_zone;
                        // xn is already set above

                        eos(eos_input_rt, eos_state);

                        Real pres_zone = eos_state.p;
                        Real dpdr = eos_state.dpdr;

                        // Newton-Raphson - we want to zero A = p_want - p(rho)
                        Real A = p_want - pres_zone;
                        drho = A / (dpdr + 0.5_rt * dx[1] * gravity::const_grav);

                        dens_zone = amrex::max(0.9_rt*dens_zone,
                                               amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

                        // convergence?

                        if (std::abs(drho) < hse::TOL * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                    }

#ifndef AMREX_USE_CUDA
                   if (! converged_hse) {
                       std::cout << "i, jj, k, domlo[1]: " << i << " " << jj << " " << k << " " << domlo[1] << std::endl;
                       std::cout << "p_want:    " << p_want << std::endl;
                       std::cout << "dens_zone: " << dens_zone << std::endl;
                       std::cout << "temp_zone: " << temp_zone << std::endl;
                       std::cout << "drho:      " << drho << std::endl;
                       std::cout << std::endl;
                       std::cout << "column info: " << std::endl;
                       std::cout << "   dens: " << adv(i,jj,k,URHO) << std::endl;
                       std::cout << "   temp: " << adv(i,jj,k,UTEMP) << std::endl;
                       amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in -Y BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(i,jj,k,UMX) = 0.0_rt;
                       adv(i,jj,k,UMY) = 0.0_rt;
                       adv(i,jj,k,UMZ) = 0.0_rt;

                   } else {

                       if (hse_reflect_vels == 1) {
                           // reflect normal, zero gradient for transverse
                           // note: we need to match the corresponding
                           // zone on the other side of the interface
                           int joff = domlo[1]-jj-1;
                           adv(i,jj,k,UMY) = -dens_zone*(adv(i,domlo[1]+joff,k,UMY) / adv(i,domlo[1]+joff,k,URHO));

                           adv(i,jj,k,UMX) = -dens_zone*(adv(i,domlo[1],k,UMX) / dens_base);
                           adv(i,jj,k,UMZ) = -dens_zone*(adv(i,domlo[1],k,UMZ) / dens_base);
                       } else {
                           // zero gradient
                           adv(i,jj,k,UMX) = dens_zone * (adv(i,domlo[1],k,UMX) / dens_base);
                           adv(i,jj,k,UMY) = dens_zone * (adv(i,domlo[1],k,UMY) / dens_base);
                           adv(i,jj,k,UMZ) = dens_zone * (adv(i,domlo[1],k,UMZ) / dens_base);
                       }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   // store the final state

                   adv(i,jj,k,URHO) = dens_zone;
                   adv(i,jj,k,UEINT) = dens_zone * eint;
                   adv(i,jj,k,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(i,jj,k,UMX) * adv(i,jj,k,UMX) +
                                 adv(i,jj,k,UMY) * adv(i,jj,k,UMY) +
                                 adv(i,jj,k,UMZ) * adv(i,jj,k,UMZ)) / dens_zone;
                   adv(i,jj,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(i,jj,k,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(i,jj,k,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_above = dens_zone;
                   pres_above = pres_zone;

                }
            });

       }


    }


    // YHI

    if (bcr[URHO].hi(1) == EXT_DIR && hi[1] > domhi[1]) {

        if (yr_ext_bc_type == EXT_HSE) {

            Box gbx(IntVect(D_DECL(lo[0], domhi[1]+1, lo[2])),
                    IntVect(D_DECL(hi[0], domhi[1]+1, hi[2])));

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {

                Real dens_below = adv(i,domhi[1],k,URHO);
                Real temp_below = adv(i,domhi[1],k,UTEMP);
                Real X_zone[NumSpec];
                for (int n = 0; n < NumSpec; n++) {
                    X_zone[n] = adv(i,domhi[1],k,UFS+n) / dens_below;
                }
#if NAUX_NET > 0
                Real aux_zone[NumAux];
                for (int n = 0; n < NumAux; n++) {
                    aux_zone[n] = adv(i,domhi[1],k,UFX+n) / dens_below;
                }
#endif

                // keep track of the density at the base of the domain

                Real dens_base = dens_below;

                // get pressure in this zone (the initial below zone)

                eos_rep_t eos_state;
                eos_state.rho = dens_below;
                eos_state.T = temp_below;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = X_zone[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    eos_state.aux[n] = aux_zone[n];
                }
#endif

                eos(eos_input_rt, eos_state);

                Real pres_below = eos_state.p;

                for (int jj = domhi[1]+1; jj <= adv_hi[1]; jj++) {

                    // HSE integration to get density, pressure

                    // initial guesses
                    Real dens_zone = dens_below;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(i,jj-1,k,UTEMP) - adv(i,jj-2,k,UTEMP);
                    } else {
                        temp_zone = temp_below;
                    }

                    bool converged_hse = false;

                    Real p_want;
                    Real drho;

                    for (int iter = 0; iter < hse::MAX_ITER; iter++) {

                        // pressure needed from HSE
                        p_want = pres_below +
                            dx[1] * 0.5_rt * (dens_zone + dens_below) * gravity::const_grav;

                        // pressure from EOS

                        eos_state.rho = dens_zone;
                        eos_state.T = temp_zone;
                        // xn is already set above

                        eos(eos_input_rt, eos_state);

                        Real pres_zone = eos_state.p;
                        Real dpdr = eos_state.dpdr;

                        // Newton-Raphson - we want to zero A = p_want - p(rho)
                        Real A = p_want - pres_zone;
                        drho = A / (dpdr - 0.5_rt * dx[1] * gravity::const_grav);

                        dens_zone = amrex::max(0.9_rt*dens_zone,
                                               amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

                        // convergence?

                        if (std::abs(drho) < hse::TOL * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                    }

#ifndef AMREX_USE_CUDA
                   if (! converged_hse) {
                       std::cout << "i, jj, k, domhi[1]: " << i << " " << jj << " " << k << " " << domhi[1] << std::endl;
                       std::cout << "p_want:    " << p_want << std::endl;
                       std::cout << "dens_zone: " << dens_zone << std::endl;
                       std::cout << "temp_zone: " << temp_zone << std::endl;
                       std::cout << "drho:      " << drho << std::endl;
                       std::cout << std::endl;
                       std::cout << "column info: " << std::endl;
                       std::cout << "   dens: " << adv(i,jj,k,URHO) << std::endl;
                       std::cout << "   temp: " << adv(i,jj,k,UTEMP) << std::endl;
                       amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in +Y BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(i,jj,k,UMX) = 0.0_rt;
                       adv(i,jj,k,UMY) = 0.0_rt;
                       adv(i,jj,k,UMZ) = 0.0_rt;

                   } else {

                       if (hse_reflect_vels == 1) {
                           // reflect normal, zero gradient for transverse
                           // note: we need to match the corresponding
                           // zone on the other side of the interface
                           int joff = jj-domhi[1]-1;
                           adv(i,jj,k,UMY) = -dens_zone * (adv(i,domhi[1]-joff,k,UMY) / adv(i,domhi[1]-joff,k,URHO));

                           adv(i,jj,k,UMX) = -dens_zone * (adv(i,domhi[1],k,UMX) / dens_base);
                           adv(i,jj,k,UMZ) = -dens_zone * (adv(i,domhi[1],k,UMZ) / dens_base);
                       } else {
                           // zero gradient
                           adv(i,jj,k,UMX) = dens_zone * (adv(i,domhi[1],k,UMX) / dens_base);
                           adv(i,jj,k,UMY) = dens_zone * (adv(i,domhi[1],k,UMY) / dens_base);
                           adv(i,jj,k,UMZ) = dens_zone * (adv(i,domhi[1],k,UMZ) / dens_base);
                       }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   // store the final state

                   adv(i,jj,k,URHO) = dens_zone;
                   adv(i,jj,k,UEINT) = dens_zone * eint;
                   adv(i,jj,k,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(i,jj,k,UMX) * adv(i,jj,k,UMX) +
                                 adv(i,jj,k,UMY) * adv(i,jj,k,UMY) +
                                 adv(i,jj,k,UMZ) * adv(i,jj,k,UMZ)) / dens_zone;
                   adv(i,jj,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(i,jj,k,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(i,jj,k,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_below = dens_zone;
                   pres_below = pres_zone;

                }
            });
        }

    }
#endif

#if AMREX_SPACEDIM == 3
    //
    // z boundaries
    //

    // ZLO

    if (bcr[URHO].lo(2) == EXT_DIR && lo[2] < domlo[2]) {

        if (zl_ext_bc_type == EXT_HSE) {

            Box gbx(IntVect(D_DECL(lo[0], lo[1], domlo[2]-1)),
                    IntVect(D_DECL(hi[0], hi[1], domlo[2]-1)));

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {

                Real dens_above = adv(i,j,domlo[2],URHO);
                Real temp_above = adv(i,j,domlo[2],UTEMP);
                Real X_zone[NumSpec];
                for (int n = 0; n < NumSpec; n++) {
                    X_zone[n] = adv(i,j,domlo[2],UFS+n) / dens_above;
                }
#if NAUX_NET > 0
                Real aux_zone[NumAux];
                for (int n = 0; n < NumAux; n++) {
                    aux_zone[n] = adv(i,j,domlo[2],UFX+n) / dens_above;
                }
#endif

                // keep track of the density at the base of the domain

                Real dens_base = dens_above;

                // get pressure in this zone (the initial above zone)

                eos_rep_t eos_state;
                eos_state.rho = dens_above;
                eos_state.T = temp_above;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = X_zone[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    eos_state.aux[n] = aux_zone[n];
                }
#endif

                eos(eos_input_rt, eos_state);

                Real pres_above = eos_state.p;

                for (int kk = domlo[2]-1; kk >= adv_lo[2]; kk--) {

                    // HSE integration to get density, pressure

                    // initial guesses

                    Real dens_zone = dens_above;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(i,j,kk+1,UTEMP) - adv(i,j,kk+2,UTEMP);
                    } else {
                        temp_zone = temp_above;
                    }

                    bool converged_hse = false;

                    Real p_want;
                    Real drho;

                    for (int iter = 0; iter < hse::MAX_ITER; iter++) {

                        // pressure needed from HSE

                        p_want = pres_above -
                            dx[2] * 0.5_rt * (dens_zone + dens_above) * gravity::const_grav;

                        // pressure from EOS

                        eos_state.rho = dens_zone;
                        eos_state.T = temp_zone;
                        // xn is already set above

                        eos(eos_input_rt, eos_state);

                        Real pres_zone = eos_state.p;
                        Real dpdr = eos_state.dpdr;

                        // Newton-Raphson - we want to zero A = p_want - p(rho)
                        Real A = p_want - pres_zone;
                        drho = A / (dpdr + 0.5_rt * dx[2] * gravity::const_grav);

                        dens_zone = amrex::max(0.9_rt*dens_zone,
                                               amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

                        // convergence?

                        if (std::abs(drho) < hse::TOL * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                    }

#ifndef AMREX_USE_CUDA
                   if (! converged_hse) {
                       std::cout << "i, j, kk, domlo[2]: " << i << " " << j << " " << kk << " " << domlo[2] << std::endl;
                       std::cout << "p_want:    " << p_want << std::endl;
                       std::cout << "dens_zone: " << dens_zone << std::endl;
                       std::cout << "temp_zone: " << temp_zone << std::endl;
                       std::cout << "drho:      " << drho << std::endl;
                       std::cout << std::endl;
                       std::cout << "column info: " << std::endl;
                       std::cout << "   dens: " << adv(i,j,kk,URHO) << std::endl;
                       std::cout << "   temp: " << adv(i,j,kk,UTEMP) << std::endl;
                       amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in -Z BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(i,j,kk,UMX) = 0.0_rt;
                       adv(i,j,kk,UMY) = 0.0_rt;
                       adv(i,j,kk,UMZ) = 0.0_rt;

                   } else {

                       if (hse_reflect_vels == 1) {
                           // reflect normal, zero gradient for transverse
                           // note: we need to match the corresponding
                           // zone on the other side of the interface
                           int koff = domlo[2]-kk-1;
                           adv(i,j,kk,UMZ) = -dens_zone * (adv(i,j,domlo[2]+koff,UMZ) / adv(i,j,domlo[2]+koff,URHO));

                           adv(i,j,kk,UMX) = -dens_zone * (adv(i,j,domlo[2],UMX) / dens_base);
                           adv(i,j,kk,UMY) = -dens_zone * (adv(i,j,domlo[2],UMY) / dens_base);
                       } else {
                           // zero gradient
                           adv(i,j,kk,UMX) = dens_zone * (adv(i,j,domlo[2],UMX) / dens_base);
                           adv(i,j,kk,UMY) = dens_zone * (adv(i,j,domlo[2],UMY) / dens_base);
                           adv(i,j,kk,UMZ) = dens_zone * (adv(i,j,domlo[2],UMZ) / dens_base);
                       }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   // store the final state

                   adv(i,j,kk,URHO) = dens_zone;
                   adv(i,j,kk,UEINT) = dens_zone * eint;
                   adv(i,j,kk,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(i,j,kk,UMX) * adv(i,j,kk,UMX) +
                                 adv(i,j,kk,UMY) * adv(i,j,kk,UMY) +
                                 adv(i,j,kk,UMZ) * adv(i,j,kk,UMZ)) / dens_zone;
                   adv(i,j,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(i,j,kk,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(i,j,kk,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_above = dens_zone;
                   pres_above = pres_zone;

                }
            });
        }

    }

    // ZHI

    if (bcr[URHO].hi(2) == EXT_DIR && hi[2] > domhi[2]) {

        if (zr_ext_bc_type == EXT_HSE) {
#ifndef AMREX_USE_CUDA
            amrex::Error("ERROR: HSE boundaries not implemented for +Z");
#endif
        }
    }
#endif

}



