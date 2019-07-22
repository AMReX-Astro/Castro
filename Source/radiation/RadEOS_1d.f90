
subroutine ca_compute_c_v(lo, hi, &
                          cv, cv_l1, cv_h1, &
                          temp, temp_l1, temp_h1, &
                          state, state_l1, state_h1) bind(C, name="ca_compute_c_v")

  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in)           :: lo(1), hi(1)
  integer, intent(in)           :: cv_l1, cv_h1
  integer, intent(in)           :: temp_l1, temp_h1
  integer, intent(in)           :: state_l1, state_h1
  real(rt)                      :: cv(cv_l1:cv_h1)
  real(rt)        , intent(in)  :: temp(temp_l1:temp_h1)
  real(rt)        , intent(in)  :: state(state_l1:state_h1,NVAR)

  integer           :: i
  real(rt)         :: rhoInv
  type(eos_t) :: eos_state

  do i = lo(1), hi(1)

     rhoInv = 1.e0_rt / state(i,URHO)
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
                       state, state_l1, state_h1) bind(C, name="ca_get_rhoe")

  use eos_module, only : eos
  use eos_type_module, only: eos_t, eos_input_rt
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer         , intent(in) :: lo(1), hi(1)
  integer         , intent(in) :: rhoe_l1, rhoe_h1
  integer         , intent(in) :: temp_l1, temp_h1
  integer         , intent(in) :: state_l1, state_h1
  real(rt)        , intent(in) :: temp(temp_l1:temp_h1)
  real(rt)        , intent(in) :: state(state_l1:state_h1,NVAR)
  real(rt)                     :: rhoe(rhoe_l1:rhoe_h1)

  integer          :: i
  real(rt)         :: rhoInv
  type(eos_t) :: eos_state

  do i = lo(1), hi(1)

     rhoInv = 1.e0_rt / state(i,URHO)
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
                                      state,state_l1,state_h1) bind(C, name="ca_compute_temp_given_rhoe")

  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only: eos_t, eos_input_re
  use meth_params_module, only : NVAR, URHO, UFS, UFX, UTEMP, small_temp

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer         , intent(in) :: lo(1),hi(1)
  integer         , intent(in) :: temp_l1,temp_h1, state_l1,state_h1
  real(rt)        , intent(in) :: state(state_l1:state_h1,NVAR)
  real(rt)                     :: temp(temp_l1:temp_h1) ! temp contains rhoe as input

  integer :: i
  real(rt)         :: rhoInv
  type (eos_t) :: eos_state

  do i = lo(1),hi(1)
     if (temp(i) .le. 0.e0_rt) then
        temp(i) = small_temp
     else
        rhoInv = 1.e0_rt / state(i,URHO)
        eos_state % rho = state(i,URHO)
        eos_state % T   = state(i,UTEMP)
        eos_state % e   = temp(i)*rhoInv
        eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = state(i,UFX:UFX+naux-1) * rhoInv

        call eos(eos_input_re, eos_state)
        temp(i) = eos_state % T
     end if
  enddo

end subroutine ca_compute_temp_given_rhoe


subroutine ca_compute_temp_given_cv(lo,hi,  &
                                    temp,  temp_l1, temp_h1, &
                                    state,state_l1,state_h1, &
                                    const_c_v, c_v_exp_m, c_v_exp_n) bind(C, name="ca_compute_temp_given_cv")

  use meth_params_module, only : NVAR, URHO

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer         , intent(in) :: lo(1),hi(1)
  integer         , intent(in) :: temp_l1,temp_h1, state_l1,state_h1
  real(rt)        , intent(in) :: state(state_l1:state_h1,NVAR)
  real(rt)                     :: temp(temp_l1:temp_h1) ! temp contains rhoe as input
  real(rt)        , intent(in) :: const_c_v, c_v_exp_m, c_v_exp_n

  integer :: i
  real(rt)         :: ex, alpha, rhoal, teff

  ex = 1.e0_rt / (1.e0_rt - c_v_exp_n)

  do i=lo(1), hi(1)
     if (c_v_exp_m .eq. 0.e0_rt) then
        alpha = const_c_v
     else
        alpha = const_c_v * state(i,URHO) ** c_v_exp_m
     endif
     rhoal = state(i,URHO) * alpha + 1.e-50_rt
     if (c_v_exp_n .eq. 0.e0_rt) then
        temp(i) = temp(i) / rhoal
     else
        teff = max(temp(i), 1.e-50_rt)
        temp(i) = ((1.e0_rt - c_v_exp_n) * teff / rhoal) ** ex
     endif
  end do

end subroutine ca_compute_temp_given_cv



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following routined are used by NEUTRINO only.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ca_compute_temp_given_reye(lo, hi, &
     temp, temp_l1, temp_h1, &
     rhoe, re_l1, re_h1, &
     ye, ye_l1, ye_h1, &
     state, state_l1, state_h1)

  use network, only: nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_re
  use meth_params_module, only : NVAR, URHO, UFS, UFX, &
       small_temp

  use castro_error_module, only : castro_error
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer         , intent(in) :: lo(1), hi(1)
  integer         , intent(in) :: temp_l1, temp_h1
  integer         , intent(in) :: re_l1, re_h1
  integer         , intent(in) :: ye_l1, ye_h1
  integer         , intent(in) :: state_l1, state_h1
  real(rt)        , intent(in) :: state(state_l1:state_h1,NVAR)
  real(rt)        , intent(in) :: rhoe(re_l1:re_h1)
  real(rt)        , intent(in) :: ye(ye_l1:ye_h1)
  real(rt)                     :: temp(temp_l1:temp_h1)

  integer          :: i
  real(rt)         :: rhoInv
  type (eos_t) :: eos_state

  do i = lo(1), hi(1)

     if (rhoe(i) .le. 0.e0_rt) then
        temp(i) = small_temp
     else

        rhoInv = 1.e0_rt / state(i,URHO)
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

        if(temp(i).lt.0.e0_rt) then
           print*,'negative temp in compute_temp_given_reye ', temp(i)
           call castro_error("Error:: ca_compute_temp_given_reye")
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
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer         , intent(in) :: lo(1), hi(1)
  integer         , intent(in) :: re_l1, re_h1
  integer         , intent(in) :: rY_l1, rY_h1
  integer         , intent(in) :: temp_l1, temp_h1
  integer         , intent(in) :: ye_l1, ye_h1
  integer         , intent(in) :: state_l1, state_h1
  real(rt)        , intent(in) :: state(state_l1:state_h1,NVAR)
  real(rt)                     :: rhoe(re_l1:re_h1)
  real(rt)                     :: rhoY(rY_l1:rY_h1)
  real(rt)        , intent(in) :: ye(ye_l1:ye_h1)
  real(rt)        , intent(in) :: temp(temp_l1:temp_h1)

  integer          :: i
  real(rt)         :: rhoInv
  type (eos_t) :: eos_state

  do i = lo(1), hi(1)

     rhoInv = 1.e0_rt / state(i,URHO)
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
