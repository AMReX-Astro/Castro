subroutine ca_state_update( lo, hi, &
     state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
     rhoe,  rhoe_l1, rhoe_l2, rhoe_l3, rhoe_h1, rhoe_h2, rhoe_h3, &
     temp,  temp_l1, temp_l2, temp_l3, temp_h1, temp_h2, temp_h3, &
     msk ,   msk_l1,  msk_l2,  msk_l3,  msk_h1,  msk_h2,  msk_h3, &
     derat, dTrat)

  use meth_params_module, only : NVAR, UEDEN, UEINT, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3) 
  integer, intent(in) :: state_l1, state_h1, state_l2, state_h2, state_l3, state_h3
  integer, intent(in) ::  rhoe_l1,  rhoe_h1,  rhoe_l2,  rhoe_h2,  rhoe_l3,  rhoe_h3
  integer, intent(in) ::  temp_l1,  temp_h1,  temp_l2,  temp_h2,  temp_l3,  temp_h3
  integer, intent(in) ::   msk_l1,   msk_h1,   msk_l2,   msk_h2,   msk_l3,   msk_h3
  real(rt)        ,intent(in):: rhoe( rhoe_l1: rhoe_h1, rhoe_l2: rhoe_h2, rhoe_l3: rhoe_h3)
  real(rt)        ,intent(in):: temp( temp_l1: temp_h1, temp_l2: temp_h2, temp_l3: temp_h3)
  real(rt)        ,intent(in)::  msk(  msk_l1:  msk_h1,  msk_l2:  msk_h2,  msk_l3:  msk_h3)
  real(rt)                   ::state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  real(rt)        , intent(inout) :: derat, dTrat

  integer :: i, j, k
  real(rt)         :: ei, ek, Told

  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     ei = state(i,j,k,UEINT)
     derat = max(derat, abs((rhoe(i,j,k) - ei)*msk(i,j,k)/ max(ei, 1.e-50_rt)))
     ek = state(i,j,k,UEDEN) - state(i,j,k,UEINT)
     state(i,j,k,UEINT) = rhoe(i,j,k)
     state(i,j,k,UEDEN) = rhoe(i,j,k) + ek

     Told = state(i,j,k,UTEMP);
     dTrat = max(dTrat, abs((temp(i,j,k)-Told)*msk(i,j,k)/ max(Told, 1.e-50_rt)))
     state(i,j,k,UTEMP) = temp(i,j,k)
  end do
  end do
  end do

end subroutine ca_state_update


subroutine ca_flux_face2center( lo, hi, &
     t, t_l1, t_l2, t_l3, t_h1, t_h2, t_h3, &
     f, f_l1, f_l2, f_l3, f_h1, f_h2, f_h3, &
     x, x_l1, x_h1, nt, idim, iflx) bind(C, name="ca_flux_face2center")

  use rad_params_module, only : ngroups
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in)::t_l1,t_h1,t_l2,t_h2,t_l3,t_h3
  integer,intent(in)::f_l1,f_h1,f_l2,f_h2,f_l3,f_h3
  integer,intent(in)::x_l1,x_h1
  integer,intent(in) :: nt, idim, iflx
  real(rt)                   ::t(t_l1:t_h1,t_l2:t_h2,t_l3:t_h3,0:nt-1)
  real(rt)        ,intent(in)::f(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)
  real(rt)        ,intent(in)::x(x_l1:x_h1)

  integer it, i, j, k

  it = idim*ngroups + iflx

  if (idim .eq. 0) then
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              t(i,j,k,it) = (f(i,j,k) + f(i+1,j,k)) * 0.5e0_rt
           end do
        end do
     end do
  else if (idim .eq. 1) then
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              t(i,j,k,it) = (f(i,j,k) + f(i,j+1,k)) * 0.5e0_rt
           end do
        end do
     end do
  else
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              t(i,j,k,it) = (f(i,j,k) + f(i,j,k+1)) * 0.5e0_rt
           end do
        end do
     end do
  end if

end subroutine ca_flux_face2center

! =======================================================================
! used by the hyperbolic solver

subroutine ca_spalpha( lo, hi, &
     spa, spa_l1, spa_l2, spa_l3, spa_h1, spa_h2, spa_h3, &
     lmx, lmx_l1, lmx_l2, lmx_l3, lmx_h1, lmx_h2, lmx_h3, &
     lmy, lmy_l1, lmy_l2, lmy_l3, lmy_h1, lmy_h2, lmy_h3, &
     lmz, lmz_l1, lmz_l2, lmz_l3, lmz_h1, lmz_h2, lmz_h3, &
     igroup) bind(C, name="ca_spalpha")

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDalpha
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: spa_l1, spa_h1, spa_l2, spa_h2, spa_l3, spa_h3
  integer, intent(in) :: lmx_l1, lmx_h1, lmx_l2, lmx_h2, lmx_l3, lmx_h3
  integer, intent(in) :: lmy_l1, lmy_h1, lmy_l2, lmy_h2, lmy_l3, lmy_h3
  integer, intent(in) :: lmz_l1, lmz_h1, lmz_l2, lmz_h2, lmz_l3, lmz_h3
  integer, intent(in) :: igroup
  real(rt)                     :: spa(spa_l1:spa_h1,spa_l2:spa_h2,spa_l3:spa_h3)
  real(rt)        , intent(in) :: lmx(lmx_l1:lmx_h1,lmx_l2:lmx_h2,lmx_l3:lmx_h3,0:ngroups-1)
  real(rt)        , intent(in) :: lmy(lmy_l1:lmy_h1,lmy_l2:lmy_h2,lmy_l3:lmy_h3,0:ngroups-1)
  real(rt)        , intent(in) :: lmz(lmz_l1:lmz_h1,lmz_l2:lmz_h2,lmz_l3:lmz_h3,0:ngroups-1)
  integer :: i,j,k
  real(rt)         :: lam

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     if ( i.eq.spa_l1 .or. i.eq.spa_h1 .or.  &
          j.eq.spa_l2 .or. j.eq.spa_h2 .or.  &
          k.eq.spa_l3 .or. k.eq.spa_h3 ) then
        lam = (lmx(i,j,k,igroup) + lmx(i+1,j  ,k  ,igroup) &
             + lmy(i,j,k,igroup) + lmy(i  ,j+1,k  ,igroup) &
             + lmz(i,j,k,igroup) + lmz(i  ,j  ,k+1,igroup)) / 6.e0_rt
        spa(i,j,k) = FLDalpha(lam)
     end if
  end do
  end do
  end do

end subroutine ca_spalpha

