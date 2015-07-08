      subroutine compute_temp(lo,hi,state,state_l1,state_l2,state_l3, &
                              state_h1,state_h2,state_h3)

      use network, only : nspec, naux
      use eos_module
      use eos_type_module
      use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, UMX, UMY, UMZ, allow_negative_energy
      use bl_constants_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: state_l1,state_l2,state_l3
      integer         , intent(in   ) :: state_h1,state_h2,state_h3
      double precision, intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,&
                                               state_l3:state_h3,NVAR)

      integer          :: i,j,k,n
      double precision :: rhoInv

      type (eos_t) :: eos_state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

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
               eos_state(i,j,k) % rho = state(i,j,k,URHO)
               eos_state(i,j,k) % T   = state(i,j,k,UTEMP) ! Initial guess for the EOS
               eos_state(i,j,k) % e   = state(i,j,k,UEINT) / state(i,j,k,URHO)

               do n = 1, nspec
                  eos_state(i,j,k) % xn(n)  = state(i,j,k,UFS+n-1) / state(i,j,k,URHO)
               enddo
               do n = 1, naux
                  eos_state(i,j,k) % aux(n) = state(i,j,k,UFX+n-1) / state(i,j,k,URHO)
               enddo

               eos_state(i,j,k) % loc = (/ i, j, k /)
            enddo
         enddo
      enddo

      call eos(eos_input_re, eos_state)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               state(i,j,k,UTEMP) = eos_state(i,j,k) % T

               ! Reset energy in case we floored
               state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state(i,j,k) % e
               state(i,j,k,UEDEN) = state(i,j,k,UEINT) &
                                  + HALF * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

            enddo
         enddo
      enddo
      
      end subroutine compute_temp
