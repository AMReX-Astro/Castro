      subroutine compute_temp(lo,hi,state,state_l1,state_h1)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UEINT, UEDEN, UTEMP, UFS, UFX, &
                                     small_temp, allow_negative_energy
      use bl_constants_module

      implicit none
      integer         , intent(in   ) :: lo(1),hi(1)
      integer         , intent(in   ) :: state_l1,state_h1
      double precision, intent(inout) :: state(state_l1:state_h1,NVAR)

      integer          :: i
      double precision :: eint,xn(nspec+naux)

      type (eos_t) :: eos_state

      do i = lo(1),hi(1)
        if (state(i,URHO) <= ZERO) then
           print *,'   '
           print *,'>>> Error: Castro_1d::compute_temp ',i
           print *,'>>> ... negative density ',state(i,URHO)
           call bl_error("Error:: Castro_1d.f90 :: compute_temp")
        end if
      enddo

      if (allow_negative_energy.eq.0) then
         do i = lo(1),hi(1)
            if (state(i,UEINT) <= ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_1d::compute_temp ',i
                print *,'>>> ... (rho e) is negative '
                call bl_error("Error:: Castro_1d.f90 :: compute_temp")
            end if
         end do
      end if

      do i = lo(1), hi(1)

         eos_state % rho = state(i,URHO)
         eos_state % e   = state(i,UEINT) / state(i,URHO)
         eos_state % xn  = state(i,UFS:UFS+nspec-1) / state(i,URHO)
         eos_state % aux = state(i,UFX:UFX+naux-1) / state(i,URHO)

         ! initial guess for iterations
         eos_state % T = state(i,UTEMP) 

      enddo

      call eos(eos_input_re, eos_state)

      do i = lo(1), hi(1)

         state(i,UTEMP) = eos_state % T

         ! In case we've floored, or otherwise allowed the energy to change, update the energy accordingly.                                         

         state(i,UEDEN) = state(i,UEDEN) + (state(i,URHO) * eos_state % e(i) - state(i,UEINT))
         state(i,UEINT) = state(i,URHO) * eos_state % e(i)

      enddo

      end subroutine compute_temp
