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



subroutine ca_ncupdate_matter( lo, hi,  &
     Tp_n,Tp_n_l1,Tp_n_l2,Tp_n_h1,Tp_n_h2,  &
     Er_n,Er_n_l1,Er_n_l2,Er_n_h1,Er_n_h2,  &
     re_s,re_s_l1,re_s_l2,re_s_h1,re_s_h2,  &
     re_2,re_2_l1,re_2_l2,re_2_h1,re_2_h2,  &
     etTz,etTz_l1,etTz_l2,etTz_h1,etTz_h2,  &
      kpp, kpp_l1, kpp_l2, kpp_h1, kpp_h2,  &
       jg,  jg_l1,  jg_l2,  jg_h1,  jg_h2,  &
     dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(2),hi(2)
  integer,intent(in)::Tp_n_l1,Tp_n_h1,Tp_n_l2,Tp_n_h2
  integer,intent(in)::Er_n_l1,Er_n_h1,Er_n_l2,Er_n_h2
  integer,intent(in)::re_s_l1,re_s_h1,re_s_l2,re_s_h2
  integer,intent(in)::re_2_l1,re_2_h1,re_2_l2,re_2_h2
  integer,intent(in)::etTz_l1,etTz_h1,etTz_l2,etTz_h2
  integer,intent(in):: kpp_l1, kpp_h1, kpp_l2, kpp_h2
  integer,intent(in)::  jg_l1,  jg_h1,  jg_l2,  jg_h2
  real(rt)                   ::Tp_n(Tp_n_l1:Tp_n_h1,Tp_n_l2:Tp_n_h2)
  real(rt)        ,intent(in)::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,0:ngroups-1)
  real(rt)        ,intent(in)::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2)
  real(rt)        ,intent(in)::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2)
  real(rt)        ,intent(in)::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2)
  real(rt)        ,intent(in):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2,0:ngroups-1)
  real(rt)        ,intent(in)::  jg(  jg_l1:  jg_h1,  jg_l2:  jg_h2,0:ngroups-1)
  real(rt)        ,intent(in) :: dt

   integer :: i,j,g
   real(rt)         :: cdt1, cpT, scrch_re
   real(rt)         :: dTemp
   real(rt)        , parameter :: fac = 0.01e0_rt

   cdt1 = 1.e0_rt / (clight * dt)
   do j = lo(2), hi(2)
   do i = lo(1), hi(1)

      cpT = 0.e0_rt
      do g = 0, ngroups-1
         cpT = cpT + kpp(i,j,g)*Er_n(i,j,g) - jg(i,j,g)
      end do

      scrch_re = cpT - (re_s(i,j) - re_2(i,j)) * cdt1

      dTemp = etTz(i,j)*scrch_re

      if (abs(dTemp/(Tp_n(i,j)+1.e-50_rt)) > fac) then
         dTemp = sign(fac*Tp_n(i,j), dTemp)
      end if

     Tp_n(i,j) = Tp_n(i,j) + dTemp

  end do
  end do

end subroutine ca_ncupdate_matter



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

! =======================================================================
! used by the hyperbolic solver

subroutine ca_spalpha( lo, hi, &
     spa, spa_l1, spa_l2, spa_h1, spa_h2, &
     lmx, lmx_l1, lmx_l2, lmx_h1, lmx_h2, &
     lmy, lmy_l1, lmy_l2, lmy_h1, lmy_h2, &
     igroup) bind(C, name="ca_spalpha")

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDalpha
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: spa_l1, spa_h1, spa_l2, spa_h2
  integer, intent(in) :: lmx_l1, lmx_h1, lmx_l2, lmx_h2
  integer, intent(in) :: lmy_l1, lmy_h1, lmy_l2, lmy_h2
  integer, intent(in) :: igroup
  real(rt)                     :: spa(spa_l1:spa_h1,spa_l2:spa_h2)
  real(rt)        , intent(in) :: lmx(lmx_l1:lmx_h1,lmx_l2:lmx_h2,0:ngroups-1)
  real(rt)        , intent(in) :: lmy(lmy_l1:lmy_h1,lmy_l2:lmy_h2,0:ngroups-1)
  integer :: i,j
  real(rt)         :: lam

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     if ( i.eq.spa_l1 .or. i.eq.spa_h1 .or.  &
          j.eq.spa_l2 .or. j.eq.spa_h2 ) then
        lam = 0.25e0_rt*(lmx(i,j,igroup) + lmx(i+1,j  ,igroup)  &
             +        lmy(i,j,igroup) + lmy(i  ,j+1,igroup))
        spa(i,j) = FLDalpha(lam)
     end if
  end do
  end do

end subroutine ca_spalpha
