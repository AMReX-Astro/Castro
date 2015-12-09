      subroutine compute_temp(lo,hi,state,s_lo,s_hi)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, allow_negative_energy
      use bl_constants_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_lo(3),s_hi(3)
      double precision, intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

      integer          :: i,j,k
      double precision :: rhoInv

      type (eos_t) :: eos_state

      ! First check the inputs for validity.
      
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               if (state(i,j,k,URHO) <= ZERO) then
                  print *,'   '
                  print *,'>>> Error: Castro_3d::compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,URHO)
                  print *,'    '
                  call bl_error("Error:: compute_temp_nd.f90")
               end if

               if (allow_negative_energy .eq. 0 .and. state(i,j,k,UEINT) <= ZERO) then
                  print *,'   '
                  print *,'>>> Warning: Castro_3d::compute_temp ',i,j,k
                  print *,'>>> ... negative (rho e) ',state(i,j,k,UEINT)
                  print *,'   '
                  call bl_error("Error:: compute_temp_nd.f90")
               end if
               
            enddo
         enddo
      enddo

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
               rhoInv = ONE / state(i,j,k,URHO)
               
               eos_state % rho = state(i,j,k,URHO)
               eos_state % T   = state(i,j,k,UTEMP) ! Initial guess for the EOS
               eos_state % e   = state(i,j,k,UEINT) * rhoInv
               eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

               call eos(eos_input_re, eos_state)

               state(i,j,k,UTEMP) = eos_state % T

               ! In case we've floored, or otherwise allowed the energy to change, update the energy accordingly.

               state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + (state(i,j,k,URHO) * eos_state % e - state(i,j,k,UEINT))
               state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e

            enddo
         enddo
      enddo
      
      end subroutine compute_temp
