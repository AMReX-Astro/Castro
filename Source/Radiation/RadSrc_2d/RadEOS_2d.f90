
subroutine ca_compute_c_v(lo, hi, &
                          cv, cv_l1, cv_l2, cv_h1, cv_h2, &
                          temp, temp_l1, temp_l2, temp_h1, temp_h2, &
                          state, state_l1, state_l2, state_h1, state_h2) &
                          bind(C, name="ca_compute_c_v")

  use eos_module
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  implicit none
  integer, intent(in)          :: lo(2), hi(2)
  integer, intent(in)          :: cv_l1, cv_l2, cv_h1, cv_h2
  integer, intent(in)          :: temp_l1, temp_l2, temp_h1, temp_h2
  integer, intent(in)          :: state_l1, state_l2, state_h1, state_h2
  double precision             :: cv(cv_l1:cv_h1,cv_l2:cv_h2)
  double precision, intent(in) :: temp(temp_l1:temp_h1,temp_l2:temp_h2)
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  integer           :: i, j
  double precision :: rhoInv
  type(eos_t) :: eos_state

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rhoInv = 1.d0 / state(i,j,URHO)
        eos_state % rho = state(i,j,URHO)
        eos_state % T   =  temp(i,j)
        eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = state(i,j,UFX:UFX+naux -1) * rhoInv

        call eos(eos_input_rt, eos_state)

        cv(i,j) = eos_state % cv

     enddo
  enddo

end subroutine ca_compute_c_v


subroutine ca_get_rhoe(lo, hi, &
                       rhoe, rhoe_l1, rhoe_l2, rhoe_h1, rhoe_h2, &
                       temp, temp_l1, temp_l2, temp_h1, temp_h2, &
                       state, state_l1, state_l2, state_h1, state_h2) &
                       bind(C, name="ca_get_rhoe")

  use eos_module
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  implicit none
  integer         , intent(in) :: lo(2), hi(2)
  integer         , intent(in) :: rhoe_l1, rhoe_l2, rhoe_h1, rhoe_h2
  integer         , intent(in) :: temp_l1, temp_l2, temp_h1, temp_h2
  integer         , intent(in) :: state_l1, state_l2, state_h1, state_h2
  double precision, intent(in) :: temp(temp_l1:temp_h1,temp_l2:temp_h2)
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision             :: rhoe(rhoe_l1:rhoe_h1,rhoe_l2:rhoe_h2)

  integer          :: i, j
  double precision :: rhoInv
  type(eos_t) :: eos_state  

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rhoInv = 1.d0 / state(i,j,URHO)
        eos_state % rho = state(i,j,URHO)
        eos_state % T   =  temp(i,j)
        eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = state(i,j,UFX:UFX+naux -1) * rhoInv

        call eos(eos_input_rt, eos_state)

        rhoe(i,j) = eos_state % rho * eos_state % e

     enddo
  enddo
end subroutine ca_get_rhoe


subroutine ca_compute_temp_given_rhoe(lo,hi,  &
                                      temp,  temp_l1, temp_l2, temp_h1, temp_h2, &
                                      state,state_l1,state_l2,state_h1,state_h2) &
                                      bind(C, name="ca_compute_temp_given_rhoe")

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UTEMP, UFS, UFX, &
                                 small_temp, allow_negative_energy

  implicit none
  integer         , intent(in) :: lo(2),hi(2)
  integer         , intent(in) :: temp_l1, temp_l2, temp_h1, temp_h2, &
       state_l1,state_l2,state_h1,state_h2
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision :: temp(temp_l1:temp_h1,temp_l2:temp_h2) ! temp contains rhoe as input

  integer :: i, j
  double precision :: rhoInv
  type (eos_t) :: eos_state

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        if (allow_negative_energy.eq.0 .and. temp(i,j).le.0.d0) then
           temp(i,j) = small_temp
        else
           rhoInv = 1.d0 / state(i,j,URHO)
           eos_state % rho = state(i,j,URHO)
           eos_state % T   = state(i,j,UTEMP)
           eos_state % e   =  temp(i,j)*rhoInv 
           eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
           eos_state % aux = state(i,j,UFX:UFX+naux -1) * rhoInv
           
           call eos(eos_input_re, eos_state)
           temp(i,j) = eos_state % T
        end if
     enddo
  enddo

end subroutine ca_compute_temp_given_rhoe


subroutine ca_compute_temp_given_cv(lo,hi,  &
                                    temp,  temp_l1, temp_l2, temp_h1, temp_h2, &
                                    state,state_l1,state_l2,state_h1,state_h2, &
                                    const_c_v, c_v_exp_m, c_v_exp_n) &
                                    bind(C, name="ca_compute_temp_given_cv")

  use meth_params_module, only : NVAR, URHO

  implicit none
  integer         , intent(in) :: lo(2),hi(2)
  integer         , intent(in) :: temp_l1, temp_l2, temp_h1, temp_h2, &
       state_l1,state_l2,state_h1,state_h2
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision :: temp(temp_l1:temp_h1,temp_l2:temp_h2) ! temp contains rhoe
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
       allow_negative_energy, small_temp

  implicit none
  integer, intent(in) :: lo(2), hi(2), resetei
  integer, intent(in) :: s_l1, s_h1, s_l2, s_h2
  double precision :: state(s_l1:s_h1,s_l2:s_h2,NVAR)
  double precision, intent(inout) :: rela, abso

  integer :: i, j
  double precision :: rhoInv,  ek, e1, e2, T1, T2, diff, rdiff
  type(eos_t) :: eos_state

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rhoInv = 1.d0 / state(i,j,URHO)
        eos_state % rho = state(i,j,URHO)
        eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = state(i,j,UFX:UFX+naux -1) * rhoInv

        ek = 0.5d0 * (state(i,j,UMX)**2 + state(i,j,UMY)**2) * rhoInv

        e1 = state(i,j,UEINT) * rhoInv
        eos_state % e = e1
        eos_state % T = state(i,j,UTEMP)

        if (allow_negative_energy.eq.0 .and. eos_state%e.lt.0.d0) then
           e2 = (state(i,j,UEDEN) - ek) * rhoInv
           if (e2 .gt. 0.d0) then
              eos_state % e = e2
              e1 = e2
              state(i,j,UEINT) = eos_state%e * eos_state%rho
           else
              eos_state % T = small_temp
              call eos(eos_input_rt, eos_state)
              e1 = eos_state%e
              state(i,j,UEINT) = eos_state%e * eos_state%rho
           end if
        end if

        call eos(eos_input_re, eos_state)

        T1 = eos_state % T

        if (resetei .eq. 0) then

           state(i,j,UEDEN) = state(i,j,UEINT) + ek
           state(i,j,UTEMP) = T1

        else if (allow_negative_energy.eq.0 .and. state(i,j,UEINT).le.1.d-4*ek) then

           state(i,j,UEDEN) = state(i,j,UEINT) + ek
           state(i,j,UTEMP) = T1

        else
           e2 = (state(i,j,UEDEN) - ek) * rhoInv
           if (allow_negative_energy.eq.1 .or. e2.gt.0.d0) then
              eos_state % e = e2
              eos_state % T = T1
              call eos(eos_input_re, eos_state)
              T2 = eos_state % T
           else
              T2 = T1
              state(i,j,UEDEN) = state(i,j,UEINT) + ek
           end if

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following routined are used by NEUTRINO only.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  integer         , intent(in) :: lo(2), hi(2)
  integer         , intent(in) ::  temp_l1, temp_l2, temp_h1, temp_h2
  integer         , intent(in) ::    re_l1,   re_l2,   re_h1,   re_h2
  integer         , intent(in) ::    ye_l1,   ye_l2,   ye_h1,   ye_h2
  integer         , intent(in) :: state_l1,state_l2,state_h1,state_h2
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision, intent(in) ::  rhoe(   re_l1:   re_h1,   re_l2:   re_h2)
  double precision, intent(in) ::    ye(   ye_l1:   ye_h1,   ye_l2:   ye_h2)
  double precision             ::  temp( temp_l1: temp_h1, temp_l2: temp_h2)

  integer          :: i, j
  double precision :: rhoInv
  type (eos_t) :: eos_state

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        if(allow_negative_energy.eq.0 .and. rhoe(i,j).le.0.d0) then
           temp(i,j) = small_temp
        else

           rhoInv = 1.d0 / state(i,j,URHO)
           eos_state % rho = state(i,j,URHO)
           ! set initial guess of temperature
           eos_state % T = temp(i,j)
           eos_state % e = rhoe(i,j)*rhoInv 
           eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv
           if (naux > 0) then
              eos_state % aux = ye(i,j)
           end if

           call eos(eos_input_re, eos_state)

           temp(i,j) = eos_state % T

           if(temp(i,j).lt.0.d0) then
              print*,'negative temp in compute_temp_given_reye ', temp(i,j)
              call bl_error("Error :: ca_compute_temp_given_reye")
           endif

        end if

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
  use eos_module
  use meth_params_module, only : NVAR, URHO, UFS, UFX

  implicit none
  integer         , intent(in) :: lo(2), hi(2)
  integer         , intent(in) :: re_l1, re_h1, re_l2, re_h2
  integer         , intent(in) :: rY_l1, rY_h1, rY_l2, rY_h2
  integer         , intent(in) :: temp_l1, temp_h1, temp_l2, temp_h2
  integer         , intent(in) :: ye_l1, ye_h1, ye_l2, ye_h2
  integer         , intent(in) :: state_l1, state_h1, state_l2, state_h2
  double precision, intent(in) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision             :: rhoe(re_l1:re_h1,re_l2:re_h2)
  double precision             :: rhoY(rY_l1:rY_h1,rY_l2:rY_h2)
  double precision, intent(in) :: ye(ye_l1:ye_h1,ye_l2:ye_h2)
  double precision, intent(in) :: temp(temp_l1:temp_h1,temp_l2:temp_h2)

  integer          :: i, j
  double precision :: rhoInv
  type (eos_t) :: eos_state

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rhoInv = 1.d0 / state(i,j,URHO)
        eos_state % rho = state(i,j,URHO)
        eos_state % T   =  temp(i,j)
        eos_state % xn  = state(i,j,UFS:UFS+nspec-1) * rhoInv

        if (naux > 0) then
           eos_state % aux = ye(i,j)
           rhoY(i,j) = state(i,j,URHO)*ye(i,j)        
        end if

        call eos(eos_input_rt, eos_state)

        rhoe(i,j) = eos_state % rho * eos_state % e

     enddo
  enddo
end subroutine ca_compute_reye_given_ty

