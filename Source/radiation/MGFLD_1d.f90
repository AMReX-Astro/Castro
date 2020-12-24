subroutine ca_state_update( lo, hi, &
     state, state_l1, state_h1,  &
     rhoe,   rhoe_l1,  rhoe_h1,  &
     temp,   temp_l1,  temp_h1,  &
     msk ,    msk_l1,   msk_h1,  &
     derat, dTrat)

  use meth_params_module, only : NVAR, UEDEN, UEINT, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1) 
  integer, intent(in) :: state_l1, state_h1
  integer, intent(in) ::  rhoe_l1,  rhoe_h1
  integer, intent(in) ::  temp_l1,  temp_h1
  integer, intent(in) ::   msk_l1,   msk_h1
  real(rt)        , intent(in)   :: rhoe( rhoe_l1: rhoe_h1)
  real(rt)        , intent(in)   :: temp( temp_l1: temp_h1)
  real(rt)        , intent(in)   ::  msk(  msk_l1:  msk_h1)
  real(rt)                       ::state(state_l1:state_h1, NVAR)
  real(rt)        , intent(inout) :: derat, dTrat

  integer :: i
  real(rt)         :: ei, ek, Told

  do i=lo(1), hi(1)
     ei = state(i,UEINT)
     derat = max(derat, abs((rhoe(i) - ei)*msk(i)/ (ei + 1.e-50_rt)))
     ek = state(i,UEDEN) - state(i,UEINT)
     state(i,UEINT) = rhoe(i)
     state(i,UEDEN) = rhoe(i) + ek

     Told = state(i,UTEMP)
     dTrat = max(dTrat, abs((temp(i)-Told)*msk(i)/ (Told + 1.e-50_rt)))
     state(i,UTEMP) = temp(i)
  end do

end subroutine ca_state_update



subroutine ca_flux_face2center( lo, hi, &
     t, t_l1, t_h1, &
     f, f_l1, f_h1, &
     x, x_l1, x_h1, &
     nt, idim, it) bind(C, name="ca_flux_face2center")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in)::t_l1,t_h1
  integer,intent(in)::f_l1,f_h1
  integer,intent(in)::x_l1,x_h1
  integer,intent(in) :: nt, idim, it
  real(rt)                   ::t(t_l1:t_h1,0:nt-1)
  real(rt)        ,intent(in)::f(f_l1:f_h1)
  real(rt)        ,intent(in)::x(x_l1:x_h1)

  integer i

  do i=lo(1),hi(1)
     t(i,it) = (f(i)/(x(i)+1.e-50_rt) + f(i+1)/x(i+1)) * 0.5e0_rt
  end do

end subroutine ca_flux_face2center

! =======================================================================
! used by the hyperbolic solver

subroutine ca_spalpha( lo, hi, &
     spa, spa_l1, spa_h1, &
     lam, lam_l1, lam_h1, &
     igroup) bind(C, name="ca_spalpha")

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDalpha
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: spa_l1, spa_h1
  integer, intent(in) :: lam_l1, lam_h1
  integer, intent(in) :: igroup
  real(rt)                     :: spa(spa_l1:spa_h1)
  real(rt)        , intent(in) :: lam(lam_l1:lam_h1, 0:ngroups-1)
  integer :: i

  i = spa_l1
  if (i.ge.lo(1) .and. i.le.hi(1)) &
       spa(i) = FLDalpha(lam(i,igroup))

  i = spa_h1
  if (i.ge.lo(1) .and. i.le.hi(1)) &
       spa(i) = FLDalpha(lam(i+1,igroup))
end subroutine ca_spalpha


