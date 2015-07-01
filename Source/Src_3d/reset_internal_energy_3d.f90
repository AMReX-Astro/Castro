      subroutine reset_internal_e(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi,verbose)

      use eos_module 
      use eos_type_module
      use network, only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, &
                                     small_temp, allow_negative_energy, &
                                     dual_energy_eta2, dual_energy_update_E_from_e
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3), verbose
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
      double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)

      ! Local variables
      integer          :: i,j,k
      integer          :: pt_index(3)
      double precision :: Up, Vp, Wp, ke, rho_eint, eint_new

      type (eos_t) :: eos_state(1)

      ! Reset internal energy
      if (allow_negative_energy .eq. 0) then

         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

              Up = u(i,j,k,UMX) / u(i,j,k,URHO)
              Vp = u(i,j,k,UMY) / u(i,j,k,URHO)
              Wp = u(i,j,k,UMZ) / u(i,j,k,URHO)
              ke = HALF * (Up**2 + Vp**2 + Wp**2)

              rho_eint = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

              ! Reset (e from e) if it's greater than eta * E.
              if (rho_eint .gt. ZERO .and. rho_eint / u(i,j,k,UEDEN) .gt. dual_energy_eta2) then

                  u(i,j,k,UEINT) = rho_eint

              ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
              else if (u(i,j,k,UEINT) .gt. ZERO .and. dual_energy_update_E_from_e) then

                 u(i,j,k,UEDEN) = u(i,j,k,UEINT) + u(i,j,k,URHO) * ke

              ! If not resetting and little e is negative ...
              else if (u(i,j,k,UEINT) .le. ZERO) then

                 eos_state(1) % rho = u(i,j,k,URHO)
                 eos_state(1) % T   = small_temp 
                 eos_state(1) % e   = ZERO
                 eos_state(1) % xn  = u(i,j,k,UFS:UFS+nspec-1) / u(i,j,k,URHO)
                 eos_state(1) % aux = u(i,j,k,UFX:UFX+naux-1) / u(i,j,k,URHO)

                 call eos(eos_input_rt, eos_state, .false., pt_index = pt_index)

                 eint_new = eos_state(1) % e

                 if (verbose .gt. 0) then
                    print *,'   '
                    print *,'>>> Warning: Castro_3d::reset_internal_energy  ',i,j,k
                    print *,'>>> ... resetting neg. e from EOS using small_temp'
                    print *,'>>> ... from ',u(i,j,k,UEINT)/u(i,j,k,URHO),' to ', eint_new
                    print *,'    '
                 end if

                 u(i,j,k,UEINT) = u(i,j,k,URHO) *  eint_new
                 u(i,j,k,UEDEN) = u(i,j,k,UEINT) + u(i,j,k,URHO) * ke

              end if
         enddo
         enddo
         enddo

      ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
      else

         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)

              Up = u(i,j,k,UMX) / u(i,j,k,URHO)
              Vp = u(i,j,k,UMY) / u(i,j,k,URHO)
              Wp = u(i,j,k,UMZ) / u(i,j,k,URHO)
              ke = HALF * (Up**2 + Vp**2 + Wp**2)

              u(i,j,k,UEINT) = u(i,j,k,UEDEN) - u(i,j,k,URHO) * ke

         enddo
         enddo
         enddo

      endif
      
      end subroutine reset_internal_e
