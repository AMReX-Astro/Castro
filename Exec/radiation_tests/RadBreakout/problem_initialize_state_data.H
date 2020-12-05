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

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real dx_sub = dx[0] / static_cast<Real>(problem::nsub);

    Real rhowind0 = problem::rhowind1 *
        std::pow(problem::rwind1/problem::rwind0, 2);
    Real rlast = model::profile(0).r(model::npts-1);
    Real rholast = model::profile(0).state(model::dens, model::npts-1);
    Real rbase = rlast * problem::rbasefac;

    Real Twind0 = Twind1 * std::pow(problem::rwind1/problem::rwind0, Tindex);
    Tlast = model_T(npts_model)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           xcl = problo(1) + delta(1) * dble(i)

           vtot = 0.e0_rt
           rho = 0.e0_rt
           T = 0.e0_rt
           u = 0.e0_rt
           Ye = 0.e0_rt
           Abar = 0.e0_rt

           do ii=0,nsub-1
              xx = xcl + (dble(ii)+0.5e0_rt) * dx_sub
              vsub = xx**2
              vtot = vtot + vsub
              if (xx .ge. model_r(npts_model)) then
                 if (xx .ge. rwind0 ) then
                    rho_tmp = rhowind1 * (rwind1/xx)**2
                 else
                    rho_tmp = rholast * (rhowind0/rholast)** &
                         ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
                 end if
                 rho = rho + rho_tmp * vsub

                 if (xx .ge. rwind0 ) then
                    T_tmp = Twind1 * (rwind1/xx)**Tindex
                 else
                    T_tmp = Tlast * (Twind0/Tlast)** &
                         ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
                 end if
                 T = T + vsub * T_tmp

                 Ye = Ye + model_Ye(npts_model) * vsub
                 Abar = Abar + model_Abar(npts_model) * vsub

                 u = u + 0.e0_rt

              else if (xx .le. model_r(1)) then
                 rho = rho + model_rho(1) * vsub
                 T = T + model_T(1) * vsub
                 u = u + 0.e0_rt
                 Ye = Ye + model_Ye(1) * vsub
                 Abar = Abar + model_Abar(1) * vsub
              else
                 rho = rho + interpolate(xx,npts_model,model_r,model_rho) * vsub
                 T   = T   + interpolate(xx,npts_model,model_r,model_T  ) * vsub
                 u   = u   + interpolate(xx,npts_model,model_r,model_v  ) * vsub
                 Ye  = Ye  + interpolate(xx,npts_model,model_r,model_Ye ) * vsub
                 Abar=Abar + interpolate(xx,npts_model,model_r,model_Abar) * vsub
              end if
           end do

           rho = rho / vtot
           T = T / vtot
           u = u / vtot
           Ye = Ye / vtot
           Abar = Abar / vtot

           invmu = (1.e0_rt+Abar*Ye)/Abar

           state(i,j,k,URHO)  = rho
           state(i,j,k,UTEMP)  = T
           state(i,j,k,UMX)   = rho * u
           state(i,j,k,UMY:UMZ) = 0.e0_rt

           ! set the composition to be all in the first species
           state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
           state(i,j,k,UFS  ) = state(i,j,k,URHO)
           state(i,j,k,UFX) = Ye*rho
           state(i,j,k,UFX+1) = invmu*rho

           ! set the internal energy via the EOS
           rhoInv = 1.e0_rt / state(i,j,k,URHO)
           eos_state % rho = state(i,j,k,URHO)
           eos_state % T   = state(i,j,k,UTEMP)
           eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
           eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                0.5e0_rt*sum(state(i,j,k,UMX:UMZ)**2)/state(i,j,k,URHO)
        end do
     end do
  end do

end subroutine ca_initdata
