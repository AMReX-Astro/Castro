      subroutine compute_temp(lo,hi,state,state_l1,state_l2,state_h1,state_h2)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UEINT, UEDEN, UMX, UMY, UTEMP, UFS, UFX, &
                                     allow_negative_energy, small_temp
      use bl_constants_module

      implicit none
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: state_l1,state_h1,state_l2,state_h2
      double precision, intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

      integer          :: i,j
      double precision :: rhoInv

      type (eos_t) :: eos_state

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        if (state(i,j,URHO) <= ZERO) then
           print *,'   '
           print *,'>>> Error: Castro_2d::compute_temp ',i,j
           print *,'>>> ... negative density in compute_temp',i,j,state(i,j,URHO)
           call bl_error("Error:: Castro_2d.f90 :: compute_temp")
        end if
      enddo
      enddo

      if (allow_negative_energy.eq.0) then
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,UEINT) <= ZERO) then
                   print *,'   '
                   print *,'>>> Warning: Castro_2d::compute_temp ',i,j
                   print *,'>>> ... (rho e) is negative '
                   call bl_error("Error:: Castro_2d.f90 :: compute_temp")
               end if
            enddo
         enddo
      end if

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            rhoInv = ONE / state(i,j,URHO)

            eos_state % rho = state(i,j,URHO)
            eos_state % e   = state(i,j,UEINT) * rhoInv
            eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
            eos_state % aux = state(i,j,UFX:UFX+naux-1) * rhoInv

            eos_state % T   = state(i,j,UTEMP) ! Initial guess for EOS
            eos_state % loc = (/ i, j, -99 /)

            call eos(eos_input_re, eos_state)

            state(i,j,UTEMP) = eos_state % T

            ! Reset energy in case we floored

            state(i,j,UEDEN) = state(i,j,UEDEN) + state(i,j,URHO) * eos_state % e - state(i,j,UEINT)
            state(i,j,UEINT) = state(i,j,URHO) * eos_state % e

         enddo
      enddo

      end subroutine compute_temp
