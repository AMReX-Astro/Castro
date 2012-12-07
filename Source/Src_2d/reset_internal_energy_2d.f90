! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine reset_internal_e(u,u_l1,u_l2,u_h1,u_h2,lo,hi,verbose)

      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UFX, &
                                     small_temp, allow_negative_energy

      implicit none

      integer          :: lo(2), hi(2), verbose
      integer          :: u_l1,u_l2,u_h1,u_h2
      double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)

      ! Local variables
      integer          :: i,j
      integer          :: pt_index(2)
      double precision :: Up, Vp, ke
      double precision :: rho_eint, eint_new
      double precision :: x_in(1:nspec+naux), dummy_pres

      ! Reset internal energy if negative.
      if (allow_negative_energy .eq. 0) then

         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           Up = u(i,j,UMX) / u(i,j,URHO)
           Vp = u(i,j,UMY) / u(i,j,URHO)

           ke = 0.5d0 * u(i,j,URHO) * (Up**2 + Vp**2)

           rho_eint = u(i,j,UEDEN) - ke

           ! Reset (rho e) if e is greater than 0.01% of E.
           if (rho_eint .gt. 0.d0 .and. rho_eint / u(i,j,UEDEN) .gt. 1.d-4) then

               u(i,j,UEINT) = rho_eint

           ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
           else if (u(i,j,UEINT) .gt. 0.d0) then

               u(i,j,UEDEN) = u(i,j,UEINT) + ke

           ! If not resetting and little e is negative ...
           else if (u(i,j,UEINT) .le. 0.d0) then

              x_in(1:nspec) = u(i,j,UFS:UFS+nspec-1) / u(i,j,URHO)
              if (naux > 0) &
                x_in(nspec+1:nspec+naux)  = u(i,j,UFX:UFX+naux -1) / u(i,j,URHO)

              pt_index(1) = i
              pt_index(2) = j

              call eos_given_RTX(eint_new, dummy_pres, u(i,j,URHO), small_temp, x_in, pt_index=pt_index)
              if (verbose .gt. 0) then
                 print *,'   '
                 print *,'>>> Warning: Castro_2d::reset_internal_energy  ',i,j 
                 print *,'>>> ... resetting neg. e from EOS using small_temp'
                 print *,'>>> ... from ',u(i,j,UEINT)/u(i,j,URHO),' to ', eint_new
                 print *,'    '
              end if

              u(i,j,UEINT) = u(i,j,URHO) * eint_new
              u(i,j,UEDEN) = u(i,j,URHO) * eint_new + ke

           end if

         enddo
         enddo

      ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
      else

         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

           Up = u(i,j,UMX) / u(i,j,URHO)
           Vp = u(i,j,UMY) / u(i,j,URHO)
           ke = 0.5d0 * u(i,j,URHO) * (Up**2 + Vp**2)

           u(i,j,UEINT) = u(i,j,UEDEN) - ke

         enddo
         enddo

      endif

    end subroutine reset_internal_e
