      subroutine compute_temp(lo,hi,state,state_l1,state_h1)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UEINT, UTEMP, UFS, UFX, &
                                     small_temp, allow_negative_energy

      implicit none
      integer         , intent(in   ) :: lo(1),hi(1)
      integer         , intent(in   ) :: state_l1,state_h1
      double precision, intent(inout) :: state(state_l1:state_h1,NVAR)

      integer          :: i
      integer          :: pt_index(1)
      double precision :: eint,xn(nspec+naux)
      double precision :: dummy_gam,dummy_pres,dummy_c,dummy_dpdr,dummy_dpde

      do i = lo(1),hi(1)
        if (state(i,URHO) <= 0.d0) then
           print *,'   '
           print *,'>>> Error: Castro_1d::compute_temp ',i
           print *,'>>> ... negative density ',state(i,URHO)
           call bl_error("Error:: Castro_1d.f90 :: compute_temp")
        end if
      enddo

      if (allow_negative_energy.eq.0) then
         do i = lo(1),hi(1)
            if (state(i,UEINT) <= 0.d0) then
                print *,'   '
                print *,'>>> Warning: Castro_1d::compute_temp ',i
                print *,'>>> ... (rho e) is negative '
                call bl_error("Error:: Castro_1d.f90 :: compute_temp")
            end if
         end do
      end if

      do i = lo(1),hi(1)

         xn(1:nspec)  = state(i,UFS:UFS+nspec-1) / state(i,URHO)
         if (naux > 0) &
           xn(nspec+1:nspec+naux) = state(i,UFX:UFX+naux-1) / state(i,URHO)

         eint = state(i,UEINT) / state(i,URHO)

         pt_index(1) = i

         call eos_given_ReX(dummy_gam, dummy_pres , dummy_c, state(i,UTEMP), &
                            dummy_dpdr, dummy_dpde, state(i,URHO), eint, xn, pt_index=pt_index)

      enddo

      end subroutine compute_temp
