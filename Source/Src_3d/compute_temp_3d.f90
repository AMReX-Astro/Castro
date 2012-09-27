      subroutine compute_temp(lo,hi,state,state_l1,state_l2,state_l3, &
                              state_h1,state_h2,state_h3)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, &
                                     UFS, UFX, allow_negative_energy

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: state_l1,state_l2,state_l3
      integer         , intent(in   ) :: state_h1,state_h2,state_h3
      double precision, intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,&
                                               state_l3:state_h3,NVAR)

      integer          :: i,j,k
      double precision :: rhoInv,eint,xn(nspec+naux)
      double precision :: dummy_gam,dummy_pres,dummy_c,dummy_dpdr,dummy_dpde
      integer          :: pt_index(3)

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        if (state(i,j,k,URHO) <= 0.d0) then
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
            if (state(i,j,k,UEINT) <= 0.d0) then
                print *,'   '
                print *,'>>> Warning: Castro_3d::compute_temp ',i,j,k
                print *,'>>> ... (rho e) is negative '
                call bl_error("Error:: Castro_3d.f90 :: compute_temp")
            end if
         enddo
         enddo
         enddo
      end if


      !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,xn,eint,pt_index,dummy_gam,dummy_pres,dummy_c,dummy_dpdr,dummy_dpde)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         rhoInv = 1.d0 / state(i,j,k,URHO)

         xn(1:nspec)  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
         if (naux > 0) &
           xn(nspec+1:nspec+naux) = state(i,j,k,UFX:UFX+naux -1) * rhoInv

         eint = state(i,j,k,UEINT) / state(i,j,k,URHO)

         pt_index(1) = i
         pt_index(2) = j
         pt_index(3) = k
         call eos_given_ReX(dummy_gam, dummy_pres , dummy_c, state(i,j,k,UTEMP), &
                            dummy_dpdr, dummy_dpde, state(i,j,k,URHO), eint, xn, pt_index=pt_index)

      enddo
      enddo
      enddo
      !$OMP END PARALLEL DO

      end subroutine compute_temp
