! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine reset_internal_e(u,u_l1,u_h1,lo,hi,verbose)

      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, small_temp, allow_negative_energy

      implicit none

      integer          :: lo(1), hi(1), verbose
      integer          :: u_l1,u_h1
      double precision :: u(u_l1:u_h1,NVAR)

      ! Local variables
      integer          :: i
      integer          :: pt_index(1)
      double precision :: Up, ke, rho_eint, eint_new, x_in(1:nspec+naux), dummy_pres

      ! Reset internal energy if negative.
      if (allow_negative_energy .eq. 0) then
         do i = lo(1),hi(1)

            Up = u(i,UMX) / u(i,URHO)
            ke = 0.5d0 * (Up**2)
   
            rho_eint = u(i,UEDEN) - u(i,URHO) * ke
   
            ! Reset (rho e) if e is greater than 0.01% of E.
            if (rho_eint .gt. 0.d0 .and. rho_eint / u(i,UEDEN) .gt. 1.d-4) then
   
                u(i,UEINT) = rho_eint
   
            ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
            else if (u(i,UEINT) .gt. 0.d0) then

               u(i,UEDEN) = u(i,UEINT) + u(i,URHO) * ke

            ! If not resetting and little e is negative ...
            else if (u(i,UEINT) .le. 0.d0) then
   
               x_in(1:nspec) = u(i,UFS:UFS+nspec-1) / u(i,URHO)
               if (naux > 0) &
                 x_in(nspec+1:nspec+naux)  = u(i,UFX:UFX+naux -1) / u(i,URHO)
   
               pt_index(1) = i
               call eos_given_RTX(eint_new, dummy_pres, u(i,URHO), small_temp, x_in, pt_index=pt_index)
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
            ke = 0.5d0 * (Up**2)

            u(i,UEINT) = u(i,UEDEN) - u(i,URHO) * ke

         enddo

      endif
 
    end subroutine reset_internal_e

