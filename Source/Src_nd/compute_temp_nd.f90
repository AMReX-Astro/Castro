      subroutine compute_temp(lo,hi,state,s_lo,s_hi)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, UMX, UMY, UMZ, allow_negative_energy
      use bl_constants_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_lo(3),s_hi(3)
      double precision, intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

      integer          :: i,j,k
      double precision :: rhoInv

      type (eos_t_3D) :: eos_state

      eos_state = eos_t_3D(lo,hi)

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        if (state(i,j,k,URHO) <= ZERO) then
           print *,'   '
           print *,'>>> Error: Castro_3d::compute_temp ',i,j,k
           print *,'>>> ... negative density ',state(i,j,k,URHO)
           print *,'    '
           call bl_error("Error:: Castro_3d.f90 :: compute_temp")
        end if
      enddo
      enddo
      enddo

      if (allow_negative_energy.eq.0) then
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (state(i,j,k,UEINT) <= ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_3d::compute_temp ',i,j,k
                print *,'>>> ... (rho e) is negative '
                call bl_error("Error:: Castro_3d.f90 :: compute_temp")
            end if
         enddo
         enddo
         enddo
      end if

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rhoInv = ONE / state(i,j,k,URHO)
               
               eos_state % rho(i,j,k) = state(i,j,k,URHO)
               eos_state % T(i,j,k)   = state(i,j,k,UTEMP) ! Initial guess for the EOS
               eos_state % e(i,j,k)   = state(i,j,k,UEINT) * rhoInv

               eos_state % xn(i,j,k,:)  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
               eos_state % aux(i,j,k,:) = state(i,j,k,UFX:UFX+naux-1) * rhoInv
            enddo
         enddo
      enddo

      call eos(eos_input_re, eos_state)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,UTEMP) = eos_state % T(i,j,k)

               ! Reset energy in case we floored
               state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e(i,j,k)
               state(i,j,k,UEDEN) = state(i,j,k,UEINT) &
                                  + HALF * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

            enddo
         enddo
      enddo
      
      end subroutine compute_temp
