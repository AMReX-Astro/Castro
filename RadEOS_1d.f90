
subroutine ca_compute_c_v(lo, hi, &
     cv, cv_l1, cv_h1, &
     temp, temp_l1, temp_h1, &
     state, state_l1, state_h1)
  
  use eos_module
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS, UFX
  
  implicit none
  integer, intent(in)           :: lo(1), hi(1)
  integer, intent(in)           :: cv_l1, cv_h1
  integer, intent(in)           :: temp_l1, temp_h1
  integer, intent(in)           :: state_l1, state_h1
  double precision, intent(out) :: cv(cv_l1:cv_h1)
  double precision, intent(in)  :: temp(temp_l1:temp_h1)
  double precision, intent(in)  :: state(state_l1:state_h1,NVAR)
  
  integer           :: i
  double precision :: rhoInv
  type(eos_t) :: eos_state
  
  do i = lo(1), hi(1)
     
     rhoInv = 1.d0 / state(i,URHO)
     eos_state % rho = state(i,URHO)
     eos_state % T = temp(i)
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = state(i,UFX:UFX+naux -1) * rhoInv

     call eos(eos_input_rt, eos_state)

     cv(i) = eos_state % cv
     
  enddo
  
end subroutine ca_compute_c_v


subroutine ca_get_rhoe(lo, hi, &
     rhoe, rhoe_l1, rhoe_h1, &
     temp, temp_l1, temp_h1, &
     state, state_l1, state_h1)
  
  use eos_module
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS, UFX
  
  implicit none
  integer         , intent(in) :: lo(1), hi(1)
  integer         , intent(in) :: rhoe_l1, rhoe_h1
  integer         , intent(in) :: temp_l1, temp_h1
  integer         , intent(in) :: state_l1, state_h1
  double precision, intent(in) :: temp(temp_l1:temp_h1)
  double precision, intent(in) :: state(state_l1:state_h1,NVAR)
  double precision, intent(inout) :: rhoe(rhoe_l1:rhoe_h1)
  
  integer          :: i
  double precision :: rhoInv
  type(eos_t) :: eos_state  

  do i = lo(1), hi(1)

     rhoInv = 1.d0 / state(i,URHO)
     eos_state % rho = state(i,URHO)
     eos_state % T   = temp(i)
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = state(i,UFX:UFX+naux -1) * rhoInv

     call eos(eos_input_rt, eos_state)

     rhoe(i) = eos_state % rho * eos_state % e
          
  enddo
end subroutine ca_get_rhoe


subroutine ca_compute_temp_given_rhoe(lo,hi,  &
     temp,  temp_l1, temp_h1, &
     state,state_l1,state_h1)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UFS, UFX, UTEMP

  implicit none
  integer         , intent(in) :: lo(1),hi(1)
  integer         , intent(in) :: temp_l1,temp_h1, state_l1,state_h1
  double precision, intent(in) :: state(state_l1:state_h1,NVAR)
  double precision, intent(inout) :: temp(temp_l1:temp_h1) ! temp contains rhoe as input

  integer :: i
  double precision :: rhoInv
  type (eos_t) :: eos_state

  do i = lo(1),hi(1)

     rhoInv = 1.d0 / state(i,URHO)
     eos_state % rho = state(i,URHO)
     eos_state % T   = state(i,UTEMP)
     eos_state % e   = temp(i)*rhoInv 
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = state(i,UFX:UFX+naux-1) * rhoInv

     call eos(eos_input_re, eos_state)

     temp(i) = eos_state % T

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
  double precision :: rhoInv,  ek, e1, e2, T1, T2, diff, rdiff
  type(eos_t) :: eos_state

  do i = lo(1), hi(1)

     rhoInv = 1.d0 / state(i,URHO)
     eos_state % rho = state(i,URHO)
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = state(i,UFX:UFX+naux -1) * rhoInv
     
     ek = 0.5d0 * state(i,UMX)**2 * rhoInv

     e1 = state(i,UEINT) * rhoInv
     eos_state % e = e1
     eos_state % T = state(i,UTEMP)

     call eos(eos_input_re, eos_state)
     
     T1 = eos_state % T

     if (resetei .eq. 0) then

        state(i,UEDEN) = state(i,UEINT) + ek
        state(i,UTEMP) = T1        

     else if (allow_negative_energy.eq.0 .and. state(i,UEINT).le.1.d-4*ek) then

        state(i,UEDEN) = state(i,UEINT) + ek
        state(i,UTEMP) = T1        

     else 
        e2 = (state(i,UEDEN) - ek) * rhoInv
        eos_state % e = e2
        eos_state % T = T1

        call eos(eos_input_re, eos_state)

        T2 = eos_state % T

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following routined are used by NEUTRINO only.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ca_compute_temp_given_reye(lo, hi, &
     temp, temp_l1, temp_h1, &
     rhoe, re_l1, re_h1, &
     ye, ye_l1, ye_h1, &
     state, state_l1, state_h1)
  
  use network, only: nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UFS, UFX, &
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

  integer          :: i
  double precision :: rhoInv
  type (eos_t) :: eos_state
  
  do i = lo(1), hi(1)

     if(allow_negative_energy.eq.0 .and. rhoe(i).le.0.d0) then
        temp(i) = small_temp
     else

        rhoInv = 1.d0 / state(i,URHO)
        eos_state % rho = state(i,URHO)
        ! set initial guess of temperature
        eos_state % T = temp(i)
        eos_state % e = rhoe(i)*rhoInv 
        eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
        if (naux > 0) then
           eos_state % aux = ye(i)
        end if

        call eos(eos_input_re, eos_state)

        temp(i) = eos_state % T

        if(temp(i).lt.0.d0) then
           print*,'negative temp in compute_temp_given_reye ', temp(i)
           call bl_error("Error:: ca_compute_temp_given_reye")
        endif
     
     end if

  enddo
end subroutine ca_compute_temp_given_reye


subroutine ca_compute_reye_given_ty(lo, hi, &
     rhoe, re_l1, re_h1, &
     rhoY, rY_l1, rY_h1, &
     temp, temp_l1, temp_h1, &
     ye, ye_l1, ye_h1, &
     state, state_l1, state_h1)
  
  use network, only: nspec, naux
  use eos_module
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
  double precision :: rhoInv
  type (eos_t) :: eos_state

  do i = lo(1), hi(1)

     rhoInv = 1.d0 / state(i,URHO)
     eos_state % rho = state(i,URHO)
     eos_state % T = temp(i)
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv

     if (naux > 0) then
        eos_state % aux = ye(i)
        rhoY(i) = state(i,URHO)*ye(i)        
     end if

     call eos(eos_input_rt, eos_state)

     rhoe(i) = eos_state % rho * eos_state % e

  enddo
end subroutine ca_compute_reye_given_ty

