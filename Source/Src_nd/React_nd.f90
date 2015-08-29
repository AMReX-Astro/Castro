
  subroutine ca_react_state(lo,hi, &
                            state,s_lo,s_hi, &
                            time,dt_react)

      use eos_module
      use network           , only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, small_dens, small_temp, allow_negative_energy
      use castro_burner_module
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: s_lo(3), s_hi(3)
      double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
      double precision :: time, dt_react

      integer          :: i,j,k,n
      double precision :: rho, rhoInv, u, v, w, ke, e_in, e_out, temp
      double precision :: delta_rho_e
      double precision :: x_in(nspec+naux), x_out(nspec+naux)

      type (eos_t) :: eos_state

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Define T from e
               rho           = state(i,j,k,URHO)
               rhoInv        = ONE / rho
               u             = state(i,j,k,UMX)*rhoInv
               v             = state(i,j,k,UMY)*rhoInv
               w             = state(i,j,k,UMZ)*rhoInv
               ke            = HALF * (u**2 + v**2 + w**2)
               temp          = state(i,j,k,UTEMP)
               x_in(1:nspec) = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
               x_in(nspec+1:nspec+naux)  = state(i,j,k,UFX:UFX+naux-1) * rhoInv

               e_in          = state(i,j,k,UEINT) * rhoInv

               eos_state % T   = temp
               eos_state % rho = rho
               eos_state % e   = e_in
               eos_state % xn  = x_in(1:nspec)
               eos_state % aux = x_in(nspec+1:nspec+naux)

               ! Use this call to define T
               
               call eos(eos_input_re, eos_state)

               temp = eos_state % T
               e_in = eos_state % e

               call burner(rho, temp, x_in, e_in, dt_react, time, x_out, e_out)

               ! Make sure that species emerge in the proper range: [0,1]
               do n = 1, nspec
                  x_out(n) = max(min(x_out(n),ONE),ZERO)
               end do

               delta_rho_e = rho * e_out - rho * e_in

               state(i,j,k,UEINT)           = state(i,j,k,UEINT) + delta_rho_e
               state(i,j,k,UEDEN)           = state(i,j,k,UEDEN) + delta_rho_e
               state(i,j,k,UFS:UFS+nspec-1) = rho * x_out(1:nspec)
               state(i,j,k,UFX:UFX+naux-1)  = rho * x_out(nspec+1:nspec+naux)
   
               ! Now update the temperature to match the new internal energy

               eos_state % rho = rho
               eos_state % e   = e_out
               eos_state % xn  = x_out(1:nspec)
               eos_state % aux = x_out(nspec+1:nspec+naux)

               call eos(eos_input_re, eos_state)

               state(i,j,k,UTEMP)           = eos_state % T

            end do
         end do
      end do

  end subroutine ca_react_state
