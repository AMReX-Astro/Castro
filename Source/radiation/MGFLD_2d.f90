subroutine ca_state_update( lo, hi, &
     state,state_l1,state_l2,state_h1,state_h2, &
     rhoe,  rhoe_l1, rhoe_l2, rhoe_h1, rhoe_h2, &
     temp,  temp_l1, temp_l2, temp_h1, temp_h2, &
     msk ,   msk_l1,  msk_l2,  msk_h1,  msk_h2, &
     derat, dTrat)

  use meth_params_module, only : NVAR, UEDEN, UEINT, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2) 
  integer, intent(in) :: state_l1, state_h1, state_l2, state_h2
  integer, intent(in) ::  rhoe_l1,  rhoe_h1,  rhoe_l2,  rhoe_h2
  integer, intent(in) ::  temp_l1,  temp_h1,  temp_l2,  temp_h2
  integer, intent(in) ::   msk_l1,   msk_h1,   msk_l2,   msk_h2
  real(rt)        , intent(in) :: rhoe( rhoe_l1: rhoe_h1, rhoe_l2: rhoe_h2)
  real(rt)        , intent(in) :: temp( temp_l1: temp_h1, temp_l2: temp_h2)
  real(rt)        , intent(in) ::  msk(  msk_l1:  msk_h1,  msk_l2:  msk_h2)
  real(rt)                     ::state(state_l1:state_h1,state_l2:state_h2,NVAR)
  real(rt)        , intent(inout) :: derat, dTrat

  integer :: i, j
  real(rt)         :: ei, ek, Told

  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     ei = state(i,j,UEINT)
     derat = max(derat, abs((rhoe(i,j) - ei)*msk(i,j)/ max(ei, 1.e-50_rt)))
     ek = state(i,j,UEDEN) - state(i,j,UEINT)
     state(i,j,UEINT) = rhoe(i,j)
     state(i,j,UEDEN) = rhoe(i,j) + ek

     Told = state(i,j,UTEMP);
     dTrat = max(dTrat, abs((temp(i,j)-Told)*msk(i,j)/ max(Told, 1.e-50_rt)))
     state(i,j,UTEMP) = temp(i,j)
  end do
  end do

end subroutine ca_state_update



subroutine ca_flux_face2center( lo, hi, &
     t, t_l1, t_l2, t_h1, t_h2, &
     f, f_l1, f_l2, f_h1, f_h2, &
     x, x_l1, x_h1, nt, idim, iflx) bind(C, name="ca_flux_face2center")

  use rad_params_module, only : ngroups
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in)::t_l1,t_h1,t_l2,t_h2
  integer,intent(in)::f_l1,f_h1,f_l2,f_h2
  integer,intent(in)::x_l1,x_h1
  integer,intent(in) :: nt, idim, iflx
  real(rt)                   ::t(t_l1:t_h1,t_l2:t_h2,0:nt-1)
  real(rt)        ,intent(in)::f(f_l1:f_h1,f_l2:f_h2)
  real(rt)        ,intent(in)::x(x_l1:x_h1)

  integer it, i, j

  it = idim*ngroups + iflx

  if (idim .eq. 0) then
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           t(i,j,it) = (f(i,j)/(x(i)+1.e-50_rt) + f(i+1,j)/x(i+1)) * 0.5e0_rt
        end do
     end do
  else 
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           t(i,j,it) = (f(i,j)/x(i) + f(i,j+1)/x(i)) * 0.5e0_rt
        end do
     end do
  end if

end subroutine ca_flux_face2center
