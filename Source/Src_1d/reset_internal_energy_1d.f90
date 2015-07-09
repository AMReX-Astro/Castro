! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine reset_internal_e(u,u_l1,u_h1,lo,hi,verbose)

      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, small_temp, allow_negative_energy, &
                                     dual_energy_eta2, dual_energy_update_E_from_e
      use bl_constants_module

      implicit none

      integer          :: lo(1), hi(1), verbose
      integer          :: u_l1,u_h1
      double precision :: u(u_l1:u_h1,NVAR)

      ! Local variables
      integer          :: i
      double precision :: Up, ke, rho_eint, eint_new

      type (eos_t) :: eos_state

      ! Reset internal energy if negative.
      if (allow_negative_energy .eq. 0) then
         do i = lo(1),hi(1)

            Up = u(i,UMX) / u(i,URHO)
            ke = HALF * (Up**2)
   
            rho_eint = u(i,UEDEN) - u(i,URHO) * ke
   
            ! Reset (rho e) if e is greater than 0.01% of E.
            if (rho_eint .gt. ZERO .and. rho_eint / u(i,UEDEN) .gt. dual_energy_eta2) then
   
                u(i,UEINT) = rho_eint
   
            ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
            else if (u(i,UEINT) .gt. ZERO .and. dual_energy_update_E_from_e) then

               u(i,UEDEN) = u(i,UEINT) + u(i,URHO) * ke

            ! If not resetting and little e is negative ...
            else if (u(i,UEINT) .le. ZERO) then

               eos_state % rho = u(i,URHO)
               eos_state % T   = small_temp
               eos_state % xn  = u(i,UFS:UFS+nspec-1) / u(i,URHO)   
               eos_state % aux = u(i,UFX:UFX+naux-1) / u(i,URHO)
               eos_state % loc = (/ i, -99, -99 /)
   
               call eos(eos_input_rt, eos_state)

               eint_new = eos_state % e

               if (verbose .gt. 0) then
                  print *,'   '
                  print *,'>>> Warning: Castro_1d::reset_internal_e ',i
                  print *,'>>> ... resetting neg. e from EOS using small_temp'
                  print *,'>>> ... from ',u(i,UEINT)/u(i,URHO),' to ', eint_new
                  print *,'    '
               end if
   
               u(i,UEINT) = u(i,URHO) *  eint_new
               u(i,UEDEN) = u(i,URHO) * (eint_new + ke)
   
            end if
         enddo

      ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
      else

         do i = lo(1),hi(1)

            Up = u(i,UMX) / u(i,URHO)
            ke = HALF * (Up**2)

            u(i,UEINT) = u(i,UEDEN) - u(i,URHO) * ke

         enddo

      endif
 
    end subroutine reset_internal_e

