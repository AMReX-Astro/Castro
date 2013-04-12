      subroutine compute_temp(lo,hi,state,state_l1,state_l2,state_h1,state_h2)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, UFX, &
                                     allow_negative_energy, small_temp

      implicit none
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: state_l1,state_h1,state_l2,state_h2
      double precision, intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

      integer          :: i,j
      integer          :: pt_index(2)
      double precision :: eint,xn(nspec+naux)
      double precision :: dummy_gam,dummy_pres,dummy_c,dummy_dpdr,dummy_dpde

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        if (state(i,j,URHO) < 0.d0) then
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
               if (state(i,j,UEINT) <= 0.d0) then
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

            xn(1:nspec)  = state(i,j,UFS:UFS+nspec-1) / state(i,j,URHO)
            if (naux > 0) &
              xn(nspec+1:nspec+naux) = state(i,j,UFX:UFX+naux-1) / state(i,j,URHO)

            eint = state(i,j,UEINT) / state(i,j,URHO)
   
            pt_index(1) = i
            pt_index(2) = j
            call eos_given_ReX(dummy_gam, dummy_pres , dummy_c, state(i,j,UTEMP), &
                               dummy_dpdr, dummy_dpde, state(i,j,URHO), eint, xn, pt_index=pt_index)
         enddo
      enddo

      end subroutine compute_temp
