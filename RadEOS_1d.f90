! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_compute_c_v(lo, hi, &
           cv, cv_l1, cv_h1, &
           ye, ye_l1, ye_h1, &
           temp, temp_l1, temp_h1, &
           state, state_l1, state_h1)

        use eos_module
        use network, only : nspec, naux
        use meth_params_module, only : NVAR, URHO, UFS, UFX

        implicit none
        integer, intent(in)           :: lo(1), hi(1)
        integer, intent(in)           :: cv_l1, cv_h1
        integer, intent(in)           :: ye_l1, ye_h1
        integer, intent(in)           :: temp_l1, temp_h1
        integer, intent(in)           :: state_l1, state_h1
        double precision, intent(out) :: cv(cv_l1:cv_h1)
        double precision, intent(in)  :: ye(ye_l1:ye_h1)
        double precision, intent(in)  :: temp(temp_l1:temp_h1)
        double precision, intent(in)  :: state(state_l1:state_h1,NVAR)

        integer           :: i
        double precision :: rho, rhoInv
        double precision :: xn(nspec+naux)

        do i = lo(1), hi(1)

           rho = state(i,URHO)
           rhoInv = 1.d0 / rho
           xn(1:nspec) = state(i,UFS:UFS+nspec-1) * rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,UFX:UFX+naux-1)*rhoInv
              xn(nspec+1) = ye(i)
           end if

           call eos_get_cv(cv(i), rho, temp(i), xn)

        enddo

      end subroutine ca_compute_c_v

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_get_rhoe(lo, hi, &
           rhoe, rhoe_l1, rhoe_h1, &
           temp, temp_l1, temp_h1, &
           ye, ye_l1, ye_h1, &
           state, state_l1, state_h1)

        use eos_module
        use network, only : nspec, naux
        use meth_params_module, only : NVAR, URHO, UMX, UMY, &
             UFS, UFX, small_temp, allow_negative_energy

        implicit none
        integer         , intent(in) :: lo(1), hi(1)
        integer         , intent(in) :: rhoe_l1, rhoe_h1
        integer         , intent(in) :: temp_l1, temp_h1
        integer         , intent(in) :: ye_l1, ye_h1
        integer         , intent(in) :: state_l1, state_h1
        double precision, intent(in) :: temp(temp_l1:temp_h1)
        double precision, intent(in) :: ye(ye_l1:ye_h1)
        double precision, intent(in) :: state(state_l1:state_h1,NVAR)
        double precision, intent(inout) :: rhoe(rhoe_l1:rhoe_h1)

        integer          :: i
        double precision :: dummy_pres
        double precision :: rho, rhoInv
        double precision :: xn(nspec+naux)

        do i = lo(1), hi(1)

           rho = state(i,URHO)
           rhoInv = 1.d0 / rho
           xn(1:nspec) = state(i,UFS:UFS+nspec-1) * rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,UFX:UFX+naux-1)*rhoInv
              xn(nspec+1) = ye(i)
           end if

           call eos_given_RTX(rhoe(i), dummy_pres, rho, temp(i), xn)

           rhoe(i) = rho * rhoe(i)

        enddo
      end subroutine ca_get_rhoe

! :::
! ::: ------------------------------------------------------------------
! :::

      ! temp enters as rhoe
      subroutine ca_compute_temp_for_rad(lo, hi, &
           temp, temp_l1, temp_h1, &
           ye, ye_l1, ye_h1, &
           tempGuess, tg_l1, tg_h1, &
           state, state_l1, state_h1)

        use network, only: nspec, naux
        use eos_module
        use meth_params_module, only : NVAR, URHO, UMX, UMY, UFS, UFX, &
             small_temp, allow_negative_energy

        implicit none
        integer         , intent(in   ) :: lo(1), hi(1)
        integer         , intent(in   ) :: temp_l1, temp_h1
        integer         , intent(in   ) :: ye_l1, ye_h1
        integer         , intent(in   ) :: tg_l1, tg_h1
        integer         , intent(in   ) :: state_l1, state_h1
        double precision, intent(in   ) :: state(state_l1:state_h1,NVAR)
        double precision, intent(in   ) :: ye(ye_l1:ye_h1)
        double precision, intent(in   ) :: tempGuess(tg_l1:tg_h1)
        double precision, intent(inout) :: temp(temp_l1:temp_h1)

        integer          :: i, n
        double precision :: u, v
        double precision :: rhoInv, e, xn(nspec+naux)
        double precision :: dummy_gam, dummy_pres, dummy_c, dummy_dpdr, dummy_dpde

        do i = lo(1), hi(1)

           rhoInv = 1.d0 / state(i,URHO)
           u = state(i,UMX) * rhoInv
!           e = temp(i)*rhoInv - 0.5d0*(u**2)
           e = temp(i)*rhoInv 
           xn(1:nspec) = state(i,UFS:UFS+nspec-1)*rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,UFX:UFX+naux-1)*rhoInv
              xn(nspec+1) = ye(i)
           end if

           if(allow_negative_energy.eq.0 .and. e.le.0.d0) then
              temp(i) = small_temp
           else
              ! set initial guess of temperature
              temp(i) = tempGuess(i)

              call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, temp(i), &
                   dummy_dpdr, dummy_dpde, state(i,URHO), e, xn)

           endif

           if(temp(i).lt.0.d0) then
              print*,'negative temp in compute_temp_for_rad ', temp(i)
              call bl_error("Error:: Compute_cv_1d.f90 :: ca_compute_temp_for_rad")
           endif

        enddo
      end subroutine ca_compute_temp_for_rad


      subroutine ca_compute_temp_given_reye(lo, hi, &
           temp, temp_l1, temp_h1, &
           rhoe, re_l1, re_h1, &
           ye, ye_l1, ye_h1, &
           state, state_l1, state_h1)

        use network, only: nspec, naux
        use eos_module
        use meth_params_module, only : NVAR, URHO, UMX, UMY, UFS, UFX, &
             small_temp, allow_negative_energy

        implicit none
        integer         , intent(in   ) :: lo(1), hi(1)
        integer         , intent(in   ) :: temp_l1, temp_h1
        integer         , intent(in   ) :: re_l1, re_h1
        integer         , intent(in   ) :: ye_l1, ye_h1
        integer         , intent(in   ) :: state_l1, state_h1
        double precision, intent(in   ) :: state(state_l1:state_h1,NVAR)
        double precision, intent(in   ) :: rhoe(re_l1:re_h1)
        double precision, intent(in   ) :: ye(ye_l1:ye_h1)
        double precision, intent(inout) :: temp(temp_l1:temp_h1)

        integer          :: i, n
        double precision :: rhoInv, e, xn(nspec+naux)
        double precision :: dummy_gam, dummy_pres, dummy_c, dummy_dpdr, dummy_dpde

        do i = lo(1), hi(1)

           rhoInv = 1.d0 / state(i,URHO)
           e = rhoe(i)*rhoInv 
           xn(1:nspec) = state(i,UFS:UFS+nspec-1)*rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,UFX:UFX+naux-1)*rhoInv
              xn(nspec+1) = ye(i)
           end if

           if(allow_negative_energy.eq.0 .and. e.le.0.d0) then
              temp(i) = small_temp
           else

              call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, temp(i), &
                   dummy_dpdr, dummy_dpde, state(i,URHO), e, xn)

           endif

           if(temp(i).lt.0.d0) then
              print*,'negative temp in compute_temp_given_reye ', temp(i)
              call bl_error("Error:: Compute_cv_1d.f90 :: ca_compute_temp_given_reye")
           endif

        enddo
      end subroutine ca_compute_temp_given_reye


      subroutine ca_compute_reye_given_ty(lo, hi, &
           rhoe, re_l1, re_h1, &
           rhoY, rY_l1, rY_h1, &
           temp, temp_l1, temp_h1, &
           ye, ye_l1, ye_h1, &
           state, state_l1, state_h1)

        use network, only: nspec, naux
        use eos_module, only : eos_given_RTX
        use meth_params_module, only : NVAR, URHO, UFS, UFX

        implicit none
        integer         , intent(in   ) :: lo(1), hi(1)
        integer         , intent(in   ) :: re_l1, re_h1
        integer         , intent(in   ) :: rY_l1, rY_h1
        integer         , intent(in   ) :: temp_l1, temp_h1
        integer         , intent(in   ) :: ye_l1, ye_h1
        integer         , intent(in   ) :: state_l1, state_h1
        double precision, intent(in   ) :: state(state_l1:state_h1,NVAR)
        double precision, intent(out  ) :: rhoe(re_l1:re_h1)
        double precision, intent(inout) :: rhoY(rY_l1:rY_h1)
        double precision, intent(in   ) :: ye(ye_l1:ye_h1)
        double precision, intent(in   ) :: temp(temp_l1:temp_h1)

        integer          :: i
        double precision :: dummy_pres
        double precision :: rho, rhoInv
        double precision :: xn(nspec+naux)

        do i = lo(1), hi(1)

           rho = state(i,URHO)
           rhoInv = 1.d0 / rho
           xn(1:nspec) = state(i,UFS:UFS+nspec-1)*rhoInv

           if (naux > 0) then
              xn(nspec+1:nspec+naux)  = state(i,UFX:UFX+naux-1)*rhoInv
              xn(nspec+1) = ye(i)
              rhoY(i) = rho*ye(i)
           end if

           call eos_given_RTX(rhoe(i), dummy_pres, rho, temp(i), xn)

           rhoe(i) = rho * rhoe(i)

        enddo
      end subroutine ca_compute_reye_given_ty


subroutine ca_compute_temp_given_rhoe(lo,hi,  &
     temp,  temp_l1, temp_h1, &
     state,state_l1,state_h1)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  implicit none
  integer         , intent(in) :: lo(1),hi(1)
  integer         , intent(in) :: temp_l1,temp_h1, state_l1,state_h1
  double precision, intent(in) :: state(state_l1:state_h1,NVAR)
  double precision, intent(inout) :: temp(temp_l1:temp_h1) ! temp contains rhoe as input

  integer :: i
  integer          :: pt_index(1)
  double precision :: eint,xn(nspec+naux)
  double precision :: dummy_gam,dummy_pres,dummy_c,dummy_dpdr,dummy_dpde

  do i = lo(1),hi(1)
     
     xn(1:nspec)  = state(i,UFS:UFS+nspec-1) / state(i,URHO)
     if (naux > 0) &
          xn(nspec+1:nspec+naux) = state(i,UFX:UFX+naux-1) / state(i,URHO)
     
     eint = temp(i) / state(i,URHO)
     
     pt_index(1) = i
     
     call eos_given_ReX(dummy_gam, dummy_pres , dummy_c, temp(i), &
          dummy_dpdr, dummy_dpde, state(i,URHO), eint, xn, pt_index)

  enddo

end subroutine ca_compute_temp_given_rhoe

subroutine ca_compute_temp_given_cv(lo,hi,  &
     temp,  temp_l1, temp_h1, &
     state,state_l1,state_h1, &
     const_c_v, c_v_exp_m, c_v_exp_n)

  use meth_params_module, only : NVAR, URHO

  implicit none
  integer         , intent(in) :: lo(1),hi(1)
  integer         , intent(in) :: temp_l1,temp_h1, state_l1,state_h1
  double precision, intent(in) :: state(state_l1:state_h1,NVAR)
  double precision, intent(inout) :: temp(temp_l1:temp_h1) ! temp contains rhoe as input
  double precision, intent(in) :: const_c_v, c_v_exp_m, c_v_exp_n

  integer :: i
  double precision :: ex, alpha, rhoal, teff

  ex = 1.d0 / (1.d0 - c_v_exp_n)

  do i=lo(1), hi(1)
     if (c_v_exp_m .eq. 0.d0) then
        alpha = const_c_v
     else
        alpha = const_c_v * state(i,URHO) ** c_v_exp_m
     endif
     rhoal = state(i,URHO) * alpha + 1.d-50
     if (c_v_exp_n .eq. 0.d0) then
        temp(i) = temp(i) / rhoal
     else
        teff = max(temp(i), 1.d-50)
        temp(i) = ((1.d0 - c_v_exp_n) * teff / rhoal) ** ex
     endif
  end do

end subroutine ca_compute_temp_given_cv


subroutine reset_eint_compute_temp( lo, hi, &
     state, s_l1, s_h1, &
     resetei, rela, abso )

  use network, only: nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, UEDEN, UEINT, URHO, UMX, UTEMP, UFS, UFX, &
       allow_negative_energy

  implicit none
  integer, intent(in) :: lo(1), hi(1), resetei
  integer, intent(in) :: s_l1, s_h1
  double precision, intent(inout) :: state(s_l1:s_h1,NVAR)
  double precision, intent(inout) :: rela, abso

  integer :: i
  double precision :: xn(nspec+naux)
  double precision :: rhoInv,  ek, e1, e2, T1, T2, diff, rdiff
  double precision :: dummy_gam, dummy_pres, dummy_c, dummy_dpdr, dummy_dpde

  do i = lo(1), hi(1)
     rhoInv = 1.d0 / state(i,URHO)
     xn(1:nspec) = state(i,UFS:UFS+nspec-1) * rhoInv
     
     if (naux > 0) &
          xn(nspec+1:nspec+naux) = state(i,UFX:UFX+naux-1) * rhoInv
     
     ek = 0.5d0 * state(i,UMX)**2 * rhoInv

     e1 = state(i,UEINT) * rhoInv
     T1 = state(i,UTEMP)
     call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, T1, &
          dummy_dpdr, dummy_dpde, state(i,URHO), e1, xn)

     if (resetei .eq. 0) then

        state(i,UEDEN) = state(i,UEINT) + ek
        state(i,UTEMP) = T1        

     else if (allow_negative_energy.eq.0 .and. state(i,UEINT).le.1.d-4*ek) then

        state(i,UEDEN) = state(i,UEINT) + ek
        state(i,UTEMP) = T1        

     else 
        e2 = (state(i,UEDEN) - ek) * rhoInv
        T2 = T1
        call eos_given_ReX(dummy_gam, dummy_pres, dummy_c, T2, &
             dummy_dpdr, dummy_dpde, state(i,URHO), e2, xn)

        diff = abs(T1-T2)
        rdiff = diff/T1
        
        abso = max(diff, abso)
        rela = max(rdiff, rela)

        if (rdiff .gt. 0.1d0) then
           state(i,UEDEN) = state(i,UEINT) + ek
           state(i,UTEMP) = T1
        else
           state(i,UEINT) = state(i,UEDEN) - ek
           state(i,UTEMP) = T2
        end if
     end if
     
  end do

end subroutine reset_eint_compute_temp

