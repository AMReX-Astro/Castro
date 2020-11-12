#ifndef HSE_FILL_H
#define HSE_FILL_H

#include <runtime_parameters.H>
#include <ext_bc_types.H>

// a hydrostatic boundary conditions -- this relies on the assumption
// that the gravitation acceleration is constant

namespace hse
{
    constexpr int MAX_ITER = 250;
    constexpr Real TOL = 1.e-8_rt;
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
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
    auto adv_lo = adv_bx.loVect();
    auto adv_hi = adv_bx.hiVect();

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
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
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

                eos_t eos_state;
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

                for (int i = domlo[0]-1; i >= adv_lo[0]; i--) {

                    // we are integrating along a column at constant i.
                    // Make sure that our starting state is well-defined

                    // HSE integration to get density, pressure

                    // initial guesses

                    Real dens_zone = dens_above;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(i+1,j,k,UTEMP) - adv(i+2,j,k,UTEMP);
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
                        std::cout << "i, j, k, domlo[0]: " << i << " " << j << " " << k << " " << domlo[0] << std::endl;
                        std::cout << "p_want:    " << p_want << std::endl;
                        std::cout << "dens_zone: " << dens_zone << std::endl;
                        std::cout << "temp_zone: " << temp_zone << std::endl;
                        std::cout << "drho:      " << drho << std::endl;
                        std::cout << std::endl;
                        std::cout << "column info: " << std::endl;
                        std::cout << "   dens: " << adv(i,j,k,URHO) << std::endl;
                        std::cout << "   temp: " << adv(i,j,k,UTEMP) << std::endl;
                        amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in -X BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(i,j,k,UMX) = 0.0_rt;
                       adv(i,j,k,UMY) = 0.0_rt;
                       adv(i,j,k,UMZ) = 0.0_rt;

                   } else {

                      if (hse_reflect_vels == 1) {
                          // reflect normal, zero gradient for transverse
                          // note: we need to match the corresponding
                          // zone on the other side of the interface
                          int ioff = domlo[0]-i-1;
                          adv(i,j,k,UMX) = -dens_zone * (adv(domlo[0]+ioff,j,k,UMX) / adv(domlo[0]+ioff,j,k,URHO));

                          adv(i,j,k,UMY) = -dens_zone * (adv(domlo[0],j,k,UMY) / dens_base);
                          adv(i,j,k,UMZ) = -dens_zone * (adv(domlo[0],j,k,UMZ) / dens_base);
                      } else {
                          // zero gradient
                          adv(i,j,k,UMX) = dens_zone * (adv(domlo[0],j,k,UMX) / dens_base);
                          adv(i,j,k,UMY) = dens_zone * (adv(domlo[0],j,k,UMY) / dens_base);
                          adv(i,j,k,UMZ) = dens_zone * (adv(domlo[0],j,k,UMZ) / dens_base);
                      }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   // store the final state

                   adv(i,j,k,URHO) = dens_zone;
                   adv(i,j,k,UEINT) = dens_zone * eint;
                   adv(i,j,k,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(i,j,k,UMX) * adv(i,j,k,UMX) +
                                 adv(i,j,k,UMY) * adv(i,j,k,UMY) +
                                 adv(i,j,k,UMZ) * adv(i,j,k,UMZ)) / dens_zone;
                   adv(i,j,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(i,j,k,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(i,j,k,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_above = dens_zone;
                   pres_above = pres_zone;

                }
            });

        } else {

#ifndef AMREX_USE_CUDA
            amrex::Error("invalid BC option");
#endif
        }

    }


    // XHI

    if (bcr[URHO].hi(0) == EXT_DIR && hi[0] > domhi[0]) {

       if (xr_ext_bc_type == EXT_HSE) {

            Box gbx(IntVect(D_DECL(domhi[0]+1, lo[1], lo[2])),
                    IntVect(D_DECL(domhi[0]+1, hi[1], hi[2])));

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
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

                eos_t eos_state;
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

                for (int i = domhi[0]+1; i <= adv_hi[0]; i++) {

                    // HSE integration to get density, pressure

                    // initial guesses
                    Real dens_zone = dens_below;

                    // temperature and species held constant in BCs

                    Real temp_zone;
                    if (hse_interp_temp == 1) {
                        temp_zone = 2*adv(i-1,j,k,UTEMP) - adv(i-2,j,k,UTEMP);
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
                       std::cout << "i, j, k, domhi[0]: " << i << " " << j << " " << k << " " << domhi[0] << std::endl;
                       std::cout << "p_want:    " << p_want << std::endl;
                       std::cout << "dens_zone: " << dens_zone << std::endl;
                       std::cout << "temp_zone: " << temp_zone << std::endl;
                       std::cout << "drho:      " << drho << std::endl;
                       std::cout << std::endl;
                       std::cout << "column info: " << std::endl;
                       std::cout << "   dens: " << adv(i,j,k,URHO) << std::endl;
                       std::cout << "   temp: " << adv(i,j,k,UTEMP) << std::endl;
                       amrex::Error("ERROR in bc_ext_fill_nd: failure to converge in +X BC");
                   }
#endif

                   // velocity

                   if (hse_zero_vels == 1) {

                       // zero normal momentum causes pi waves to pass through

                       adv(i,j,k,UMX) = 0.0_rt;
                       adv(i,j,k,UMY) = 0.0_rt;
                       adv(i,j,k,UMZ) = 0.0_rt;

                   } else {

                       if (hse_reflect_vels == 1) {
                           // reflect normal, zero gradient for transverse
                           // note: we need to match the corresponding
                           // zone on the other side of the interface
                           int ioff = i-domhi[0]-1;
                           adv(i,j,k,UMX) = -dens_zone * (adv(domhi[0]-ioff,j,k,UMX) / adv(domhi[0]-ioff,j,k,URHO));

                           adv(i,j,k,UMY) = -dens_zone * (adv(domhi[0],j,k,UMY) / dens_base);
                           adv(i,j,k,UMZ) = -dens_zone * (adv(domhi[0],j,k,UMZ) / dens_base);
                       } else {
                           // zero gradient
                           adv(i,j,k,UMX) = dens_zone * (adv(domhi[0],j,k,UMX) / dens_base);
                           adv(i,j,k,UMY) = dens_zone * (adv(domhi[0],j,k,UMY) / dens_base);
                           adv(i,j,k,UMZ) = dens_zone * (adv(domhi[0],j,k,UMZ) / dens_base);
                       }
                   }

                   eos_state.rho = dens_zone;
                   eos_state.T = temp_zone;

                   eos(eos_input_rt, eos_state);

                   Real pres_zone = eos_state.p;
                   Real eint = eos_state.e;

                   //  store the final state

                   adv(i,j,k,URHO) = dens_zone;
                   adv(i,j,k,UEINT) = dens_zone * eint;
                   adv(i,j,k,UEDEN) = dens_zone * eint +
                       0.5_rt * (adv(i,j,k,UMX) * adv(i,j,k,UMX) +
                                 adv(i,j,k,UMY) * adv(i,j,k,UMY) +
                                 adv(i,j,k,UMZ) * adv(i,j,k,UMZ)) / dens_zone;
                   adv(i,j,k,UTEMP) = temp_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       adv(i,j,k,UFS+n) = dens_zone * X_zone[n];
                   }
#if NAUX_NET > 0
                   for (int n = 0; n < NumAux; n++) {
                       adv(i,j,k,UFX+n) = dens_zone * aux_zone[n];
                   }
#endif

                   // for the next zone

                   dens_below = dens_zone;
                   pres_below = pres_zone;

                }
            });

       } else {

#ifndef AMREX_USE_CUDA
           amrex::Error("invalid BC option");
#endif
       }

    }


#if AMREX_SPACEDIM >= 2
    !-------------------------------------------------------------------------
    ! y boundaries
    !-------------------------------------------------------------------------

    ! YLO
    if (bc(2,1,1) == EXT_DIR .and. lo(2) < domlo(2)) then

       if (yl_ext == EXT_HSE) then

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                ! we are integrating along a column at constant i.
                ! Make sure that our starting state is well-defined
                dens_above = adv(i,domlo(2),k,URHO)
                temp_above = adv(i,domlo(2),k,UTEMP)
                X_zone(:) = adv(i,domlo(2),k,UFS:UFS-1+nspec)/dens_above
                aux_zone(:) = adv(i,domlo(2),k,UFX:UFX-1+naux)/dens_above

                ! keep track of the density at the base of the domain
                dens_base = dens_above

                ! get pressure in this zone (the initial above zone)
                eos_state%rho = dens_above
                eos_state%T = temp_above
                eos_state%xn(:) = X_zone(:)
                eos_state%aux(:) = aux_zone(:)

                call eos(eos_input_rt, eos_state)

                eint = eos_state%e
                pres_above = eos_state%p

                ! integrate downward
                jmin = adv_lo(2)
                jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
                ! For CUDA, this should only be one thread doing the work:
                ! we'll arbitrary choose the zone with index domlo(2) - 1.
                if (hi(2) /= jmax) then
                   jmax = jmin - 1
                end if
#endif
                do j = jmax, jmin, -1
                   y = problo(2) + delta(2)*(dble(j) + HALF)

                   ! HSE integration to get density, pressure

                   ! initial guesses
                   dens_zone = dens_above

                   ! temperature and species held constant in BCs
                   if (hse_interp_temp == 1) then
                      temp_zone = 2*adv(i,j+1,k,UTEMP) - adv(i,j+2,k,UTEMP)
                   else
                      temp_zone = temp_above
                   endif

                   converged_hse = .FALSE.


                   do iter = 1, MAX_ITER

                      ! pressure needed from HSE
                      p_want = pres_above - &
                           delta(2)*HALF*(dens_zone + dens_above)*const_grav

                      ! pressure from EOS
                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone(:)
                      eos_state%aux(:) = aux_zone(:)

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      dpdr = eos_state%dpdr
                      eint = eos_state%e

                      ! Newton-Raphson - we want to zero A = p_want - p(rho)
                      A = p_want - pres_zone
                      drho = A/(dpdr + HALF*delta(2)*const_grav)

                      dens_zone = max(0.9_rt*dens_zone, &
                           min(dens_zone + drho, 1.1_rt*dens_zone))

                      ! convergence?
                      if (abs(drho) < TOL*dens_zone) then
                         converged_hse = .TRUE.
                         exit
                      endif

                   enddo

#ifndef AMREX_USE_CUDA
                   if (.not. converged_hse) then
                      print *, "i, j, k,domlo(2): ", i, j, k, domlo(2)
                      print *, "p_want:    ", p_want
                      print *, "dens_zone: ", dens_zone
                      print *, "temp_zone: ", temp_zone
                      print *, "drho:      ", drho
                      print *, " "
                      print *, "column info: "
                      print *, "   dens: ", adv(i,j:domlo(2),k,URHO)
                      print *, "   temp: ", adv(i,j:domlo(2),k,UTEMP)
                      call castro_error("ERROR in bc_ext_fill_nd: failure to converge in -Y BC")
                   endif
#endif

                   ! velocity
                   if (hse_zero_vels == 1) then

                      ! zero normal momentum causes pi waves to pass through
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMY) = ZERO
                      adv(i,j,k,UMZ) = ZERO

                   else

                      if (hse_reflect_vels == 1) then
                         ! reflect normal, zero gradient for transverse
                         ! note: we need to match the corresponding
                         ! zone on the other side of the interface
                         joff = domlo(2)-j-1
                         adv(i,j,k,UMY) = -dens_zone*(adv(i,domlo(2)+joff,k,UMY)/adv(i,domlo(2)+joff,k,URHO))

                         adv(i,j,k,UMX) = -dens_zone*(adv(i,domlo(2),k,UMX)/dens_base)
                         adv(i,j,k,UMZ) = -dens_zone*(adv(i,domlo(2),k,UMZ)/dens_base)
                      else
                         ! zero gradient
                         adv(i,j,k,UMX) = dens_zone*(adv(i,domlo(2),k,UMX)/dens_base)
                         adv(i,j,k,UMY) = dens_zone*(adv(i,domlo(2),k,UMY)/dens_base)
                         adv(i,j,k,UMZ) = dens_zone*(adv(i,domlo(2),k,UMZ)/dens_base)
                      endif
                   endif
                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone(:)
                   eos_state%aux(:) = aux_zone(:)

                   call eos(eos_input_rt, eos_state)

                   pres_zone = eos_state%p
                   eint = eos_state%e

                   ! store the final state
                   adv(i,j,k,URHO) = dens_zone
                   adv(i,j,k,UEINT) = dens_zone*eint
                   adv(i,j,k,UEDEN) = dens_zone*eint + &
                        HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
                   adv(i,j,k,UTEMP) = temp_zone
                   adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                   ! for the next zone
                   dens_above = dens_zone
                   pres_above = pres_zone

                end do
             end do
          end do
#ifndef AMREX_USE_CUDA
       else
          call castro_error("invalid BC option")
#endif
       endif  ! yl_ext check


    endif


    ! YHI
    if (bc(2,2,1) == EXT_DIR .and. hi(2) > domhi(2)) then

       if (yr_ext == EXT_HSE) then

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                ! we are integrating along a column at constant i.
                ! Make sure that our starting state is well-defined
                dens_below = adv(i,domhi(2),k,URHO)
                temp_below = adv(i,domhi(2),k,UTEMP)
                X_zone(:) = adv(i,domhi(2),k,UFS:UFS-1+nspec)/dens_below

                ! keep track of the density at the base of the domain
                dens_base = dens_below

                ! get pressure in this zone (the initial below zone)
                eos_state%rho = dens_below
                eos_state%T = temp_below
                eos_state%xn(:) = X_zone(:)
                eos_state%aux(:) = aux_zone(:)

                call eos(eos_input_rt, eos_state)

                eint = eos_state%e
                pres_below = eos_state%p

                ! integrate upward
                jmin = domhi(2)+1
                jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
                ! For CUDA, this should only be one thread doing the work:
                ! we'll arbitrary choose the zone with index domlo(1) - 1.
                if (hi(2) /= jmax) then
                   jmax = jmin - 1
                end if
#endif
                do j = jmin, jmax
                   y = problo(2) + delta(2)*(dble(j) + HALF)

                   ! HSE integration to get density, pressure

                   ! initial guesses
                   dens_zone = dens_below

                   ! temperature and species held constant in BCs
                   if (hse_interp_temp == 1) then
                      temp_zone = 2*adv(i,j-1,k,UTEMP) - adv(i,j-2,k,UTEMP)
                   else
                      temp_zone = temp_below
                   endif

                   converged_hse = .FALSE.

                   do iter = 1, MAX_ITER

                      ! pressure needed from HSE
                      p_want = pres_below + &
                           delta(2)*HALF*(dens_zone + dens_below)*const_grav

                      ! pressure from EOS
                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone(:)
                      eos_state%aux(:) = aux_zone(:)

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      dpdr = eos_state%dpdr
                      eint = eos_state%e

                      ! Newton-Raphson - we want to zero A = p_want - p(rho)
                      A = p_want - pres_zone
                      drho = A/(dpdr - HALF*delta(2)*const_grav)

                      dens_zone = max(0.9_rt*dens_zone, &
                           min(dens_zone + drho, 1.1_rt*dens_zone))

                      ! convergence?
                      if (abs(drho) < TOL*dens_zone) then
                         converged_hse = .TRUE.
                         exit
                      endif

                   enddo

#ifndef AMREX_USE_CUDA
                   if (.not. converged_hse) then
                      print *, "i, j, k, domhi(2): ", i, j, k, domhi(2)
                      print *, "p_want:    ", p_want
                      print *, "dens_zone: ", dens_zone
                      print *, "temp_zone: ", temp_zone
                      print *, "drho:      ", drho
                      print *, " "
                      print *, "column info: "
                      print *, "   dens: ", adv(i,j:domhi(2),k,URHO)
                      print *, "   temp: ", adv(i,j:domhi(2),k,UTEMP)
                      call castro_error("ERROR in bc_ext_fill_nd: failure to converge in +Y BC")
                   endif
#endif

                   ! velocity
                   if (hse_zero_vels == 1) then

                      ! zero normal momentum causes pi waves to pass through
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMY) = ZERO
                      adv(i,j,k,UMZ) = ZERO

                   else

                      if (hse_reflect_vels == 1) then
                         ! reflect normal, zero gradient for transverse
                         ! note: we need to match the corresponding
                         ! zone on the other side of the interface
                         joff = j-domhi(2)-1
                         adv(i,j,k,UMY) = -dens_zone*(adv(i,domhi(2)-joff,k,UMY)/adv(i,domhi(2)-joff,k,URHO))

                         adv(i,j,k,UMX) = -dens_zone*(adv(i,domhi(2),k,UMX)/dens_base)
                         adv(i,j,k,UMZ) = -dens_zone*(adv(i,domhi(2),k,UMZ)/dens_base)
                      else
                         ! zero gradient
                         adv(i,j,k,UMX) = dens_zone*(adv(i,domhi(2),k,UMX)/dens_base)
                         adv(i,j,k,UMY) = dens_zone*(adv(i,domhi(2),k,UMY)/dens_base)
                         adv(i,j,k,UMZ) = dens_zone*(adv(i,domhi(2),k,UMZ)/dens_base)
                      endif
                   endif
                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone(:)
                   eos_state%aux(:) = aux_zone(:)

                   call eos(eos_input_rt, eos_state)

                   pres_zone = eos_state%p
                   eint = eos_state%e

                   ! store the final state
                   adv(i,j,k,URHO) = dens_zone
                   adv(i,j,k,UEINT) = dens_zone*eint
                   adv(i,j,k,UEDEN) = dens_zone*eint + &
                        HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
                   adv(i,j,k,UTEMP) = temp_zone
                   adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                   ! for the next zone
                   dens_below = dens_zone
                   pres_below = pres_zone

                end do
             end do
          end do
#ifndef AMREX_USE_CUDA
       else
          call castro_error("invalid BC option")
#endif
       end if  ! yr_ext check

    endif
#endif

#if AMREX_SPACEDIM == 3

    !-------------------------------------------------------------------------
    ! z boundaries
    !-------------------------------------------------------------------------

    ! ZLO
    if (bc(3,1,1) == EXT_DIR .and. lo(3) < domlo(3)) then

       if (zl_ext == EXT_HSE) then

          do j= lo(2), hi(2)
             do i = lo(1), hi(1)
                ! we are integrating along a column at constant i.
                ! Make sure that our starting state is well-defined
                dens_above = adv(i,j,domlo(3),URHO)
                temp_above = adv(i,j,domlo(3),UTEMP)
                X_zone(:) = adv(i,j,domlo(3),UFS:UFS-1+nspec)/dens_above
                aux_zone(:) = adv(i,j,domlo(3),UFX:UFX-1+naux)/dens_above

                ! keep track of the density at the base of the domain
                dens_base = dens_above

                ! get pressure in this zone (the initial above zone)
                eos_state%rho = dens_above
                eos_state%T = temp_above
                eos_state%xn(:) = X_zone(:)
                eos_state%aux(:) = aux_zone(:)

                call eos(eos_input_rt, eos_state)

                eint = eos_state%e
                pres_above = eos_state%p

                ! integrate downward
                kmin = adv_lo(3)
                kmax = domlo(3) - 1
#ifdef AMREX_USE_CUDA
                if (hi(3) /= kmax) then
                   kmax = kmin - 1
                end if
#endif
                do k = kmax, kmin, -1
                   z = problo(3) + delta(3)*(dble(k) + HALF)

                   ! HSE integration to get density, pressure

                   ! initial guesses
                   dens_zone = dens_above

                   ! temperature and species held constant in BCs
                   if (hse_interp_temp == 1) then
                      temp_zone = 2*adv(i,j,k+1,UTEMP) - adv(i,j,k+2,UTEMP)
                   else
                      temp_zone = temp_above
                   endif

                   converged_hse = .FALSE.


                   do iter = 1, MAX_ITER

                      ! pressure needed from HSE
                      p_want = pres_above - &
                           delta(3)*HALF*(dens_zone + dens_above)*const_grav

                      ! pressure from EOS
                      eos_state%rho = dens_zone
                      eos_state%T = temp_zone
                      eos_state%xn(:) = X_zone(:)
                      eos_state%aux(:) = aux_zone(:)

                      call eos(eos_input_rt, eos_state)

                      pres_zone = eos_state%p
                      dpdr = eos_state%dpdr
                      eint = eos_state%e

                      ! Newton-Raphson - we want to zero A = p_want - p(rho)
                      A = p_want - pres_zone
                      drho = A/(dpdr + HALF*delta(3)*const_grav)

                      dens_zone = max(0.9_rt*dens_zone, &
                           min(dens_zone + drho, 1.1_rt*dens_zone))

                      ! convergence?
                      if (abs(drho) < TOL*dens_zone) then
                         converged_hse = .TRUE.
                         exit
                      endif

                   enddo

#ifndef AMREX_USE_CUDA
                   if (.not. converged_hse) then
                      print *, "i, j, k,domlo(3): ", i, j, k, domlo(3)
                      print *, "p_want:    ", p_want
                      print *, "dens_zone: ", dens_zone
                      print *, "temp_zone: ", temp_zone
                      print *, "drho:      ", drho
                      print *, " "
                      print *, "column info: "
                      print *, "   dens: ", adv(i,j,k:domlo(3),URHO)
                      print *, "   temp: ", adv(i,j,k:domlo(3),UTEMP)
                      call castro_error("ERROR in bc_ext_fill_1d: failure to converge in -Z BC")
                   endif
#endif

                   ! velocity
                   if (hse_zero_vels == 1) then

                      ! zero normal momentum causes pi waves to pass through
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMY) = ZERO
                      adv(i,j,k,UMZ) = ZERO

                   else

                      if (hse_reflect_vels == 1) then
                         ! reflect normal, zero gradient for transverse
                         ! note: we need to match the corresponding
                         ! zone on the other side of the interface
                         koff = domlo(3)-k-1
                         adv(i,j,k,UMZ) = -dens_zone*(adv(i,j,domlo(3)+koff,UMZ)/adv(i,j,domlo(3)+koff,URHO))

                         adv(i,j,k,UMX) = -dens_zone*(adv(i,j,domlo(3),UMX)/dens_base)
                         adv(i,j,k,UMY) = -dens_zone*(adv(i,j,domlo(3),UMY)/dens_base)
                      else
                         ! zero gradient
                         adv(i,j,k,UMX) = dens_zone*(adv(i,j,domlo(3),UMX)/dens_base)
                         adv(i,j,k,UMY) = dens_zone*(adv(i,j,domlo(3),UMY)/dens_base)
                         adv(i,j,k,UMZ) = dens_zone*(adv(i,j,domlo(3),UMZ)/dens_base)
                      endif
                   endif
                   eos_state%rho = dens_zone
                   eos_state%T = temp_zone
                   eos_state%xn(:) = X_zone
                   eos_state%aux(:) = aux_zone

                   call eos(eos_input_rt, eos_state)

                   pres_zone = eos_state%p
                   eint = eos_state%e

                   ! store the final state
                   adv(i,j,k,URHO) = dens_zone
                   adv(i,j,k,UEINT) = dens_zone*eint
                   adv(i,j,k,UEDEN) = dens_zone*eint + &
                        HALF*sum(adv(i,j,k,UMX:UMZ)**2)/dens_zone
                   adv(i,j,k,UTEMP) = temp_zone
                   adv(i,j,k,UFS:UFS-1+nspec) = dens_zone*X_zone(:)

                   ! for the next zone
                   dens_above = dens_zone
                   pres_above = pres_zone

                end do
             end do
          end do
#ifndef AMREX_USE_CUDA
       else
          call castro_error("invalid BC option")
#endif
       endif  ! zl_ext check


    endif

    ! ZHI
    if (bc(3,2,1) == EXT_DIR .and. hi(3) > domhi(3)) then

       if (zr_ext == EXT_HSE) then
#ifndef AMREX_USE_CUDA
          call castro_error("ERROR: HSE boundaries not implemented for +Z")
#endif
#ifndef AMREX_USE_CUDA
       else
          call castro_error("invalid BC option")
#endif
       end if  ! zr_ext check

    end if
#endif

}

#if 0

#ifndef CXX_MODEL_PARSER
  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only: problo
    use model_parser_module, only: idens_model, interpolate_sub
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real(rt) :: x, y, z

    !$gpu

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    ! XLO
    if ( bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       imin = adv_lo(1)
       imax = domlo(1)-1
#ifdef AMREX_USE_CUDA
       ! For CUDA, this should only be one thread doing the work:
       ! we'll arbitrary choose the zone with index domlo(2) - 1.
       if (hi(1) /= imax) then
          imax = imin - 1
       end if
#endif
       do i = imin, imax
          x = problo(1) + delta(1)*(dble(i) + HALF)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                call interpolate_sub(adv(i,j,k), x, idens_model)
             end do
          end do
       end do
    end if

    ! XHI
    if ( bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       imin = domhi(1)+1
       imax = adv_hi(1)
#ifdef AMREX_USE_CUDA
       if (lo(1) /= imin) then
          imin = imax + 1
       end if
#endif
       do i = imin, imax
          x = problo(1) + delta(1)*(dble(i)+ HALF)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                call interpolate_sub(adv(i,j,k), x, idens_model)
             end do
          end do
       end do
    endif

#if AMREX_SPACEDIM >= 2
    ! YLO
    if ( bc(2,1) == EXT_DIR .and. lo(2) < domlo(2)) then
       jmin = adv_lo(2)
       jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
       ! For CUDA, this should only be one thread doing the work:
       ! we'll arbitrary choose the zone with index domlo(2) - 1.
       if (hi(2) /= jmax) then
          jmax = jmin - 1
       end if
#endif
       do j = jmin, jmax
          y = problo(2) + delta(2)*(dble(j) + HALF)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), y, idens_model)
             end do
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2) == EXT_DIR .and. hi(2) > domhi(2)) then
       jmin = domhi(2)+1
       jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
       if (lo(2) /= jmin) then
          jmin = jmax + 1
       end if
#endif
       do j = jmin, jmax
          y = problo(2) + delta(2)*(dble(j)+ HALF)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), y, idens_model)
             end do
          end do
       end do
    end if
#endif

#if AMREX_SPACEDIM == 3
    ! ZLO
    if ( bc(3,1) == EXT_DIR .and. lo(3) < domlo(3)) then
       kmin = adv_lo(3)
       kmax = domlo(3)-1
#ifdef AMREX_USE_CUDA
       if (hi(3) /= kmax) then
          kmax = kmin - 1
       end if
#endif
       do k = kmin, kmax
          z = problo(3) + delta(3)*(dble(k) + HALF)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), z, idens_model)
             end do
          end do
       end do
    end if

    ! ZHI
    if ( bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then
       kmin = domhi(3)+1
       kmax = adv_hi(3)
#ifdef AMREX_USE_CUDA
       if (lo(3) /= kmin) then
          kmin = kmax + 1
       end if
#endif
       do k = kmin, kmax
          z = problo(3) + delta(3)*(dble(k)+ HALF)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), z, idens_model)
             end do
          end do
       end do
    end if
#endif


#endif

#endif

#endif
