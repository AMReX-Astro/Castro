! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_compute_c_v(lo, hi, &
           cv, cv_l1, cv_l2, cv_h1, cv_h2, &
           ye, ye_l1, ye_l2, ye_h1, ye_h2, &
           temp, temp_l1, temp_l2, temp_h1, temp_h2, &
           state, state_l1, state_l2, state_h1, state_h2)

        use eos_module
        use network, only : nspec, naux
        use meth_params_module, only : NVAR, URHO, UFS, UFX

        implicit none
        integer, intent(in)           :: lo(2), hi(2)
        integer, intent(in)           :: cv_l1, cv_l2, cv_h1, cv_h2
        integer, intent(in)           :: ye_l1, ye_l2, ye_h1, ye_h2
        integer, intent(in)           :: temp_l1, temp_l2, temp_h1, temp_h2
        integer, intent(in)           :: state_l1, state_l2, state_h1, state_h2
        double precision, intent(out) :: cv(cv_l1:cv_h1,cv_l2:cv_h2)
        double precision, intent(in)  :: ye(ye_l1:ye_h1,ye_l2:ye_h2)
        double precision, intent(in)  :: temp(temp_l1:temp_h1,temp_l2:temp_h2)
        double precision, intent(in)  :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

        integer           :: i, j
        double precision :: rho, rhoInv
        double precision :: xn(nspec+naux)

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              rho = state(i,j,URHO)
              rhoInv = 1.d0 / rho
              xn(1:nspec) = state(i,j,UFS:UFS+nspec-1) * rhoInv

              if (naux > 0) then
                 xn(nspec+1:nspec+naux)  = state(i,j,UFX:UFX+naux-1) * rhoInv
                 xn(nspec+1) = ye(i,j)
              end if

              call eos_get_cv(cv(i,j), rho, temp(i,j), xn)

           enddo
        enddo

      end subroutine ca_compute_c_v

! :::
! ::: ------------------------------------------------------------------
! :::
! :::

      subroutine ca_get_rhoe(lo, hi, &
           rhoe, rhoe_l1, rhoe_l2, rhoe_h1, rhoe_h2, &
           temp, temp_l1, temp_l2, temp_h1, temp_h2, &
           ye, ye_l1, ye_l2, ye_h1, ye_h2, &
           state, state_l1, state_l2, state_h1, state_h2)

        use eos_module
        use network, only : nspec, naux
        use meth_params_module, only : NVAR, URHO, UMX, UMY, &
             UFS, UFX, small_temp, allow_negative_energy

        implicit none
        integer         , intent(in) :: lo(2), hi(2)
        integer         , intent(in) :: rhoe_l1, rhoe_l2, rhoe_h1, rhoe_h2
        integer         , intent(in) :: temp_l1, temp_l2, temp_h1, temp_h2
        integer         , intent(in) :: ye_l1, ye_l2, ye_h1, ye_h2
        integer         , intent(in) :: state_l1, state_l2, state_h1, state_h2
        double precision, intent(in) :: temp(temp_l1:temp_h1,temp_l2:temp_h2)
        double precision, intent(in) :: ye(ye_l1:ye_h1,ye_l2:ye_h2)
        double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
        double precision, intent(inout) :: rhoe(rhoe_l1:rhoe_h1,rhoe_l2:rhoe_h2)

        integer          :: i, j
        double precision :: dummy_pres
        double precision :: rho, rhoInv
        double precision :: xn(nspec+naux)

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              rho = state(i,j,URHO)
              rhoInv = 1.d0 / rho
              xn(1:nspec) = state(i,j,UFS:UFS+nspec-1) * rhoInv

              if (naux > 0) then
                 xn(nspec+1:nspec+naux)  = state(i,j,UFX:UFX+naux-1) * rhoInv
                 xn(nspec+1) = ye(i,j)
              end if

              call eos_given_RTX(rhoe(i,j), dummy_pres, rho, temp(i,j), xn)

              rhoe(i,j) = rho * rhoe(i,j)

           enddo
        enddo
      end subroutine ca_get_rhoe

! :::
! ::: ------------------------------------------------------------------
! :::
! :::

! temp enters as rhoe
      subroutine ca_compute_temp_for_rad(lo, hi, &
           temp, temp_l1, temp_l2, temp_h1, temp_h2, &
           ye, ye_l1, ye_l2, ye_h1, ye_h2, &
           tempGuess, tg_l1, tg_l2, tg_h1, tg_h2, &
           state, state_l1, state_l2, state_h1, state_h2)

        use network, only: nspec, naux
        use eos_module
        use meth_params_module, only : NVAR, URHO, UMX, UMY, UFS, UFX, &
             small_temp, allow_negative_energy

        implicit none
        integer         , intent(in   ) :: lo(2), hi(2)
        integer         , intent(in   ) :: temp_l1, temp_l2, temp_h1, temp_h2
        integer         , intent(in   ) :: ye_l1, ye_l2, ye_h1, ye_h2
        integer         , intent(in   ) :: tg_l1, tg_l2, tg_h1, tg_h2
        integer         , intent(in   ) :: state_l1, state_l2, state_h1, state_h2
        double precision, intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
        double precision, intent(in   ) :: ye(ye_l1:ye_h1,ye_l2:ye_h2)
        double precision, intent(in   ) :: tempGuess(tg_l1:tg_h1,tg_l2:tg_h2)
        double precision, intent(inout) :: temp(temp_l1:temp_h1,temp_l2:temp_h2)

        integer          :: i, j
        double precision :: u, v
        double precision :: rhoInv, e, xn(nspec+naux)
        double precision :: dummy_gam, dummy_pres, dummy_c, dummy_dpdr, dummy_dpde

        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              rhoInv = 1.d0 / state(i,j,URHO)
              u = state(i,j,UMX) * rhoInv
              v = state(i,j,UMY) * rhoInv
!              e = temp(i,j)*rhoInv - 0.5d0*(u**2+v**2)
! temp enters as rhoe, not rhoE
              e = temp(i,j)*rhoInv
              xn(1:nspec) = state(i,j,UFS:UFS+nspec-1)*rhoInv

              if (naux > 0) then
                 xn(nspec+1:nspec+naux)  = state(i,j,UFX:UFX+naux-1) * rhoInv
                 xn(nspec+1) = ye(i,j)
              end if

              if(allow_negative_energy.eq.0 .and. e.le.0.d0) then
                 temp(i,j) = small_temp
              else
                 ! set initial guess of temperature
                 temp(i,j) = tempGuess(i,j)

                 call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, temp(i,j), &
                      dummy_dpdr, dummy_dpde, state(i,j,URHO), e, xn)

              endif

              if(temp(i,j).lt.0.d0) then
                 print *,'negative temp in compute_temp_for_rad ', temp(i,j)
                 call bl_error("Error:: Compute_cv_2d.f90 :: ca_compute_temp_for_rad")
              endif

           enddo
        enddo
      end subroutine ca_compute_temp_for_rad


      subroutine ca_compute_temp_given_reye(lo, hi, &
           temp , temp_l1, temp_l2, temp_h1, temp_h2, &
           rhoe ,   re_l1,   re_l2,   re_h1,   re_h2, &
           ye   ,   ye_l1,   ye_l2,   ye_h1,   ye_h2, &
           state,state_l1,state_l2,state_h1,state_h2)

        use network, only: nspec, naux
        use eos_module
        use meth_params_module, only : NVAR, URHO, UMX, UMY, UFS, UFX, &
             small_temp, allow_negative_energy

        implicit none
        integer         , intent(in   ) :: lo(2), hi(2)
        integer         , intent(in   ) ::  temp_l1, temp_l2, temp_h1, temp_h2
        integer         , intent(in   ) ::    re_l1,   re_l2,   re_h1,   re_h2
        integer         , intent(in   ) ::    ye_l1,   ye_l2,   ye_h1,   ye_h2
        integer         , intent(in   ) :: state_l1,state_l2,state_h1,state_h2
        double precision, intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
        double precision, intent(in   ) ::  rhoe(   re_l1:   re_h1,   re_l2:   re_h2)
        double precision, intent(in   ) ::    ye(   ye_l1:   ye_h1,   ye_l2:   ye_h2)
        double precision, intent(inout) ::  temp( temp_l1: temp_h1, temp_l2: temp_h2)

        integer          :: i, j, n
        double precision :: rhoInv, e, xn(nspec+naux)
        double precision :: dummy_gam, dummy_pres, dummy_c, dummy_dpdr, dummy_dpde

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rhoInv = 1.d0 / state(i,j,URHO)
           e = rhoe(i,j)*rhoInv 
           xn(1:nspec) = state(i,j,UFS:UFS+nspec-1)*rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,j,UFX:UFX+naux-1) * rhoInv
              xn(nspec+1) = ye(i,j)
           end if

           if(allow_negative_energy.eq.0 .and. e.le.0.d0) then
              temp(i,j) = small_temp
           else

              call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, temp(i,j), &
                   dummy_dpdr, dummy_dpde, state(i,j,URHO), e, xn)

           endif

           if(temp(i,j).lt.0.d0) then
              print*,'negative temp in compute_temp_given_reye ', temp(i,j)
              call bl_error("Error:: Compute_cv_2d.f90 :: ca_compute_temp_given_reye")
           endif

        enddo
        enddo
      end subroutine ca_compute_temp_given_reye


      subroutine ca_compute_reye_given_ty(lo, hi, &
           rhoe, re_l1, re_l2, re_h1, re_h2, &
           rhoY, rY_l1, rY_l2, rY_h1, rY_h2, &
           temp, temp_l1, temp_l2, temp_h1, temp_h2, &
           ye, ye_l1, ye_l2, ye_h1, ye_h2, &
           state, state_l1, state_l2, state_h1, state_h2)

        use network, only: nspec, naux
        use eos_module, only : eos_given_RTX
        use meth_params_module, only : NVAR, URHO, UFS, UFX

        implicit none
        integer         , intent(in   ) :: lo(2), hi(2)
        integer         , intent(in   ) :: re_l1, re_h1, re_l2, re_h2
        integer         , intent(in   ) :: rY_l1, rY_h1, rY_l2, rY_h2
        integer         , intent(in   ) :: temp_l1, temp_h1, temp_l2, temp_h2
        integer         , intent(in   ) :: ye_l1, ye_h1, ye_l2, ye_h2
        integer         , intent(in   ) :: state_l1, state_h1, state_l2, state_h2
        double precision, intent(in   ) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
        double precision, intent(out  ) :: rhoe(re_l1:re_h1,re_l2:re_h2)
        double precision, intent(inout) :: rhoY(rY_l1:rY_h1,rY_l2:rY_h2)
        double precision, intent(in   ) :: ye(ye_l1:ye_h1,ye_l2:ye_h2)
        double precision, intent(in   ) :: temp(temp_l1:temp_h1,temp_l2:temp_h2)

        integer          :: i, j
        double precision :: dummy_pres
        double precision :: rho, rhoInv
        double precision :: xn(nspec+naux)

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rho = state(i,j,URHO)
           rhoInv = 1.d0 / rho
           xn(1:nspec) = state(i,j,UFS:UFS+nspec-1)*rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,j,UFX:UFX+naux-1) * rhoInv
              xn(nspec+1) = ye(i,j)
              rhoY(i,j) = rho*ye(i,j)
           end if

           call eos_given_RTX(rhoe(i,j), dummy_pres, rho, temp(i,j), xn)

           rhoe(i,j) = rho * rhoe(i,j)

        enddo
        enddo
      end subroutine ca_compute_reye_given_ty



subroutine ca_compute_temp_given_rhoe(lo,hi,  &
     temp,  temp_l1, temp_l2, temp_h1, temp_h2, &
     state,state_l1,state_l2,state_h1,state_h2)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  implicit none
  integer         , intent(in) :: lo(2),hi(2)
  integer         , intent(in) :: temp_l1, temp_l2, temp_h1, temp_h2, &
                                 state_l1,state_l2,state_h1,state_h2
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision, intent(inout) :: temp(temp_l1:temp_h1,temp_l2:temp_h2) ! temp contains rhoe as input

  integer :: i, j
  integer          :: pt_index(2)
  double precision :: eint,xn(nspec+naux)
  double precision :: dummy_gam,dummy_pres,dummy_c,dummy_dpdr,dummy_dpde

  do j = lo(2),hi(2)
  do i = lo(1),hi(1)
     
     xn(1:nspec)  = state(i,j,UFS:UFS+nspec-1) / state(i,j,URHO)
     if (naux > 0) &
          xn(nspec+1:nspec+naux) = state(i,j,UFX:UFX+naux-1) / state(i,j,URHO)
     
     eint = temp(i,j) / state(i,j,URHO) 
     
     pt_index(1) = i
     pt_index(2) = j
     
     call eos_given_ReX(dummy_gam, dummy_pres , dummy_c, temp(i,j), &
          dummy_dpdr, dummy_dpde, state(i,j,URHO), eint, xn, pt_index)
     
  enddo
  enddo

end subroutine ca_compute_temp_given_rhoe


subroutine ca_compute_temp_given_cv(lo,hi,  &
     temp,  temp_l1, temp_l2, temp_h1, temp_h2, &
     state,state_l1,state_l2,state_h1,state_h2, &
     const_c_v, c_v_exp_m, c_v_exp_n)

  use meth_params_module, only : NVAR, URHO

  implicit none
  integer         , intent(in) :: lo(2),hi(2)
  integer         , intent(in) :: temp_l1, temp_l2, temp_h1, temp_h2, &
                                  state_l1,state_l2,state_h1,state_h2
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision, intent(inout) :: temp(temp_l1:temp_h1,temp_l2:temp_h2) ! temp contains rho
  double precision, intent(in) :: const_c_v, c_v_exp_m, c_v_exp_n

  integer :: i, j
  double precision :: ex, alpha, rhoal, teff

  ex = 1.d0 / (1.d0 - c_v_exp_n)
  
  do j=lo(2), hi(2)
     do i=lo(1), hi(1)
        if (c_v_exp_m .eq. 0.d0) then
           alpha = const_c_v
        else
           alpha = const_c_v * state(i,j,URHO) ** c_v_exp_m
        endif
        rhoal = state(i,j,URHO) * alpha + 1.d-50
        if (c_v_exp_n .eq. 0.d0) then
           temp(i,j) = temp(i,j) / rhoal
        else
           teff = max(temp(i,j), 1.d-50)
           temp(i,j) = ((1.d0 - c_v_exp_n) * teff / rhoal) ** ex
        endif
     end do
  end do

end subroutine ca_compute_temp_given_cv


subroutine reset_eint_compute_temp( lo, hi, &
     state, s_l1, s_l2, s_h1, s_h2, &
     resetei, rela, abso )

  use network, only: nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UMX, UMY, UTEMP, UFS, UFX, &
       allow_negative_energy

  implicit none
  integer, intent(in) :: lo(2), hi(2), resetei
  integer, intent(in) :: s_l1, s_h1, s_l2, s_h2
  double precision, intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,NVAR)
  double precision, intent(inout) :: rela, abso

  integer :: i, j
  double precision :: xn(nspec+naux)
  double precision :: rhoInv,  ek, e1, e2, T1, T2, diff, rdiff
  double precision :: dummy_gam, dummy_pres, dummy_c, dummy_dpdr, dummy_dpde

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rhoInv = 1.d0 / state(i,j,URHO)
     xn(1:nspec) = state(i,j,UFS:UFS+nspec-1) * rhoInv
     
     if (naux > 0) &
          xn(nspec+1:nspec+naux) = state(i,j,UFX:UFX+naux-1) * rhoInv
     
     ek = 0.5d0 * (state(i,j,UMX)**2 + state(i,j,UMY)**2) &
          * rhoInv

     e1 = state(i,j,UEINT) * rhoInv
     T1 = state(i,j,UTEMP)
     call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, T1, &
          dummy_dpdr, dummy_dpde, state(i,j,URHO), e1, xn)

     if (resetei .eq. 0) then

        state(i,j,UEDEN) = state(i,j,UEINT) + ek
        state(i,j,UTEMP) = T1

     else if (allow_negative_energy.eq.0 .and. state(i,j,UEINT).le.1.d-4*ek) then

        state(i,j,UEDEN) = state(i,j,UEINT) + ek
        state(i,j,UTEMP) = T1

     else

        e2 = (state(i,j,UEDEN) - ek) * rhoInv
        T2 = T1
        call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, T2, &
             dummy_dpdr, dummy_dpde, state(i,j,URHO), e2, xn)

        diff = abs(T1-T2)
        rdiff = diff/T1
        
        abso = max(diff, abso)
        rela = max(rdiff, rela)

        if (rdiff .gt. 0.1d0) then
           state(i,j,UEDEN) = state(i,j,UEINT) + ek
           state(i,j,UTEMP) = T1
        else
           state(i,j,UEINT) = state(i,j,UEDEN) - ek
           state(i,j,UTEMP) = T2
        end if
           
     end if
     
  end do
  end do

end subroutine reset_eint_compute_temp

