! begin photon routine

subroutine ca_accel_acoe( lo, hi,  &
     eta1,eta1_l1,eta1_h1,   &
     spc , spc_l1, spc_h1,   &
     kap , kap_l1, kap_h1,   &
     aco , aco_l1, aco_h1,   &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: eta1_l1,eta1_h1
  integer, intent(in) ::  spc_l1, spc_h1
  integer, intent(in) ::  kap_l1, kap_h1
  integer, intent(in) ::  aco_l1, aco_h1
  real(rt)        , intent(in ) :: eta1(eta1_l1:eta1_h1)
  real(rt)        , intent(in ) :: spc ( spc_l1: spc_h1,0:ngroups-1)
  real(rt)        , intent(in ) :: kap ( kap_l1: kap_h1,0:ngroups-1)
  real(rt)                      :: aco ( aco_l1: aco_h1)
  real(rt)        , intent(in) :: dt, tau

  integer :: i
  real(rt)         :: kbar, H1, dt1

  dt1 = (1.e0_rt+tau)/dt

  do i = lo(1), hi(1)
     kbar = sum(spc(i,:) * kap(i,:))
     H1 = eta1(i)
     aco(i) = H1*kbar*clight + dt1
  end do

end subroutine ca_accel_acoe


subroutine ca_accel_rhs( lo, hi, &
     Ern , Ern_l1, Ern_h1,  &
     Erl , Erl_l1, Erl_h1,  &
     kap , kap_l1, kap_h1,   &
     etaT,etaT_l1,etaT_h1,   &
     rhs , rhs_l1, rhs_h1,   &
     dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Ern_l1, Ern_h1
  integer, intent(in) :: Erl_l1, Erl_h1
  integer, intent(in) :: kap_l1, kap_h1
  integer, intent(in) ::etaT_l1,etaT_h1
  integer, intent(in) :: rhs_l1, rhs_h1
  real(rt)        , intent(in ) ::Ern ( Ern_l1: Ern_h1,0:ngroups-1)
  real(rt)        , intent(in ) ::Erl ( Erl_l1: Erl_h1,0:ngroups-1)
  real(rt)        , intent(in ) :: kap( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        , intent(in ) ::etaT(etaT_l1:etaT_h1)
  real(rt)                      :: rhs(rhs_l1:rhs_h1)
  real(rt)        , intent(in) :: dt

  integer :: i 
  real(rt)         :: rt_term, H

  do i = lo(1), hi(1)
     rt_term = sum(kap(i,:)*(Ern(i,:)-Erl(i,:)))
     H = etaT(i)
     rhs(i) = clight*H*rt_term
  end do

end subroutine ca_accel_rhs


subroutine ca_accel_spec(lo, hi, &
     kap , kap_l1, kap_h1,  &
     mugT,mugT_l1,mugT_h1,  &
     spec,spec_l1,spec_h1,  &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in) :: lo(1), hi(1)
  integer,intent(in):: kap_l1, kap_h1
  integer,intent(in)::mugT_l1,mugT_h1
  integer,intent(in)::spec_l1,spec_h1
  real(rt)        ,intent(in )::kap ( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        ,intent(in )::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)                    ::spec(spec_l1:spec_h1,0:ngroups-1)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i 
  real(rt)         :: cdt1, sumeps
  real(rt)        ,dimension(0:ngroups-1):: epsilon, kapt

  cdt1 = 1.e0_rt/(clight*dt)

  do i = lo(1), hi(1)
     kapt = kap(i,:) + (1.e0_rt+tau)*cdt1
     epsilon = mugT(i,:) / kapt
     sumeps = sum(epsilon)
     if (sumeps .eq. 0.e0_rt) then
        spec(i,:) = 0.e0_rt
     else
        spec(i,:) = epsilon / sumeps
     end if
  end do
  
end subroutine ca_accel_spec


subroutine ca_compute_rhs( lo, hi, &
     rhs , rhs_l1, rhs_h1, &
     jg  ,  jg_l1,  jg_h1, &
     mugT,mugT_l1,mugT_h1, &
     cpT , cpT_l1, cpT_h1, &
     etaT,etaT_l1,etaT_h1, &
     Er2 , Er2_l1, Er2_h1, &
     re2 , re2_l1, re2_h1, &
     Ers , Ers_l1, Ers_h1, &
     res , res_l1, res_h1, &
     r, dt, igroup, tau) bind(C, name="ca_compute_rhs")

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1) 
  integer,intent(in):: rhs_l1, rhs_h1
  integer,intent(in)::  jg_l1,  jg_h1
  integer,intent(in)::mugT_l1,mugT_h1
  integer,intent(in):: cpT_l1, cpT_h1
  integer,intent(in)::etaT_l1,etaT_h1
  integer,intent(in):: Er2_l1, Er2_h1
  integer,intent(in):: re2_l1, re2_h1
  integer,intent(in):: Ers_l1, Ers_h1
  integer,intent(in):: res_l1, res_h1
  real(rt)                    ::rhs ( rhs_l1: rhs_h1)
  real(rt)        ,intent(in )::jg  (  jg_l1:  jg_h1,0:ngroups-1)
  real(rt)        ,intent(in )::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in )::cpT ( cpT_l1: cpT_h1)
  real(rt)        ,intent(in )::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in )::Er2 ( Er2_l1: Er2_h1,0:ngroups-1)
  real(rt)        ,intent(in )::re2 ( re2_l1: re2_h1)
  real(rt)        ,intent(in )::Ers ( Ers_l1: Ers_h1,0:ngroups-1)
  real(rt)        ,intent(in )::res ( res_l1: res_h1)
  real(rt)        ,intent(in) ::   r(lo(1):hi(1))
  real(rt)        ,intent(in) :: dt, tau
  integer, intent(in) :: igroup

  integer :: i
  real(rt)         :: Hg, dt1

  dt1 = 1.e0_rt/dt
  do i=lo(1),hi(1)
     Hg = mugT(i,igroup) * etaT(i)

     rhs(i) = clight*(jg(i,igroup) + Hg*cpT(i))  &
          + dt1 * (Er2(i,igroup) - Hg*(res(i)-re2(i)) &
          &        + tau*Ers(i,igroup))

     rhs(i) = r(i) * rhs(i)
   end do

end subroutine ca_compute_rhs


subroutine ca_compute_rhs_so( lo, hi, & ! MG Su-Olson
     rhs , rhs_l1, rhs_h1, &
     jg  ,  jg_l1,  jg_h1, &
     mugT,mugT_l1,mugT_h1, &
     cpT , cpT_l1, cpT_h1, &
     etaT,etaT_l1,etaT_h1, &
     Er2 , Er2_l1, Er2_h1, &
     re2 , re2_l1, re2_h1, &
     res , res_l1, res_h1, &
     x, t, dt, igroup) bind(C, name="ca_compute_rhs_so")

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1) 
  integer,intent(in):: rhs_l1, rhs_h1
  integer,intent(in)::  jg_l1,  jg_h1
  integer,intent(in)::mugT_l1,mugT_h1
  integer,intent(in):: cpT_l1, cpT_h1
  integer,intent(in)::etaT_l1,etaT_h1
  integer,intent(in):: Er2_l1, Er2_h1
  integer,intent(in):: re2_l1, re2_h1
  integer,intent(in):: res_l1, res_h1
  real(rt)                    ::rhs ( rhs_l1: rhs_h1)
  real(rt)        ,intent(in )::jg  (  jg_l1:  jg_h1,0:ngroups-1)
  real(rt)        ,intent(in )::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in )::cpT ( cpT_l1: cpT_h1)
  real(rt)        ,intent(in )::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in )::Er2 ( Er2_l1: Er2_h1,0:ngroups-1)
  real(rt)        ,intent(in )::re2 ( re2_l1: re2_h1)
  real(rt)        ,intent(in )::res ( res_l1: res_h1)
  real(rt)        ,intent(in) :: x(lo(1):hi(1))
  real(rt)        ,intent(in) :: t, dt
  integer, intent(in) :: igroup

  real(rt)        , parameter :: x0 = 0.5e0_rt
  real(rt)        , parameter :: t0 = 3.3356409519815202e-10_rt 
  real(rt)        , parameter :: qn = 1.134074546528399e20_rt 

  integer :: i
  real(rt)         :: Hg

  do i=lo(1),hi(1)
     Hg = mugT(i,igroup)*etaT(i)
     rhs(i) = clight*jg(i,igroup) + clight*cpt(i)*Hg &
          + (Er2(i,igroup) - (res(i)-re2(i))*Hg) / dt
     if (t .le. t0 .and. abs(x(i)) .le. x0) then
        rhs(i) = rhs(i) + qn ! (qn / dt) * dt
     end if
   end do

end subroutine ca_compute_rhs_so


subroutine ca_local_accel( lo, hi,  &
     Ern , Ern_l1, Ern_h1,  &
     Erl , Erl_l1, Erl_h1,  &
     kap , kap_l1, kap_h1,  &
     etaT,etaT_l1,etaT_h1,  &
     mugT,mugT_l1,mugT_h1,  &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: Ern_l1, Ern_h1
  integer,intent(in):: Erl_l1, Erl_h1
  integer,intent(in):: kap_l1, kap_h1
  integer,intent(in)::etaT_l1,etaT_h1
  integer,intent(in)::mugT_l1,mugT_h1
  real(rt)                      ::Ern ( Ern_l1: Ern_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::Erl ( Erl_l1: Erl_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::kap ( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in   )::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i 
  real(rt)         :: cdt1, rt_term, p
  real(rt)        ,dimension(0:ngroups-1)::Hg, epsilon, kapt, kk

  cdt1 = 1.e0_rt/(clight*dt)

  do i = lo(1), hi(1)
     rt_term = sum(kap(i,:)*(Ern(i,:)-Erl(i,:)))

     Hg = mugT(i,:)*etaT(i)

     kapt = kap(i,:) + (1.e0_rt+tau)*cdt1
     kk = kap(i,:) / kapt

     p = 1.e0_rt-sum(Hg*kk)
     epsilon = (Hg * rt_term) / (kapt*p + 1.e-50_rt)

     Ern(i,:) = Ern(i,:) + epsilon
  end do

end subroutine ca_local_accel


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



subroutine ca_ncupdate_matter( lo, hi,  &
     Tp_n, Tp_n_l1, Tp_n_h1,  &
     Er_n, Er_n_l1, Er_n_h1,  &
     re_s, re_s_l1, re_s_h1,  &
     re_2, re_2_l1, re_2_h1,  &
     etTz, etTz_l1, etTz_h1,  &
      kpp,  kpp_l1,  kpp_h1,  &
       jg,   jg_l1,   jg_h1,  &
     dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(1),hi(1)
  integer,intent(in)::Tp_n_l1, Tp_n_h1
  integer,intent(in)::Er_n_l1, Er_n_h1
  integer,intent(in)::re_s_l1, re_s_h1
  integer,intent(in)::re_2_l1, re_2_h1
  integer,intent(in)::etTz_l1, etTz_h1
  integer,intent(in):: kpp_l1,  kpp_h1
  integer,intent(in)::  jg_l1,   jg_h1
  real(rt)                   ::Tp_n(Tp_n_l1:Tp_n_h1)
  real(rt)        ,intent(in)::Er_n(Er_n_l1:Er_n_h1,0:ngroups-1)
  real(rt)        ,intent(in)::re_s(re_s_l1:re_s_h1)
  real(rt)        ,intent(in)::re_2(re_2_l1:re_2_h1)
  real(rt)        ,intent(in)::etTz(etTz_l1:etTz_h1)
  real(rt)        ,intent(in):: kpp( kpp_l1: kpp_h1,0:ngroups-1)
  real(rt)        ,intent(in)::  jg(  jg_l1:  jg_h1,0:ngroups-1)
  real(rt)        ,intent(in) :: dt

   integer :: i,g
   real(rt)         :: cdt1, cpT, scrch_re
   real(rt)         :: dTemp
   real(rt)        , parameter :: fac = 0.01e0_rt

   cdt1 = 1.e0_rt / (clight * dt)
   do i = lo(1), hi(1)

      cpT = 0.e0_rt
      do g = 0, ngroups-1
         cpT = cpT + kpp(i,g)*Er_n(i,g) - jg(i,g)
      end do

      scrch_re = cpT - (re_s(i) - re_2(i)) * cdt1

      dTemp = etTz(i)*scrch_re

      if (abs(dTemp/(Tp_n(i)+1.e-50_rt)) > fac) then
         dTemp = sign(fac*Tp_n(i), dTemp)
      end if

     Tp_n(i) = Tp_n(i) + dTemp

  end do

end subroutine ca_ncupdate_matter



! end photon routines
! ========================================================================

! other routines that work for both photon and neutrinos

subroutine ca_accel_ccoe( lo, hi, &
     bcgr, bcgr_l1, bcgr_h1,  &
     spec, spec_l1, spec_h1,  &
     ccoe, ccoe_l1, ccoe_h1,  &
     dx, idim, igroup)

  use rad_params_module, only : ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: bcgr_l1, bcgr_h1
  integer, intent(in) :: spec_l1, spec_h1
  integer, intent(in) :: ccoe_l1, ccoe_h1
  real(rt)        , intent(in) :: bcgr(bcgr_l1:bcgr_h1)
  real(rt)        , intent(in) :: spec(spec_l1:spec_h1, 0:ngroups-1)
  real(rt)                     :: ccoe(ccoe_l1:ccoe_h1, 0:1)
  real(rt)        , intent(in) :: dx(1)
  integer, intent(in) :: idim, igroup

  integer :: i
  real(rt)         :: grad_spec, foo

  do i = lo(1), hi(1)
     grad_spec = (spec(i,igroup) - spec(i-1,igroup)) / dx(1)
     foo = - 0.5e0_rt * bcgr(i) * grad_spec
     ccoe(i,0) = ccoe(i,0) + foo
     ccoe(i,1) = ccoe(i,1) + foo
  end do

end subroutine ca_accel_ccoe


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

subroutine ca_rhstoer(lo, hi, rhs, rhs_l1, rhs_h1, r, dt)
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer,intent(in):: lo(1), hi(1), rhs_l1, rhs_h1
  real(rt)                   ::rhs ( rhs_l1: rhs_h1)
  real(rt)        ,intent(in)::   r(lo(1):hi(1))
  real(rt)        ,intent(in):: dt
  integer :: i
  do i=lo(1),hi(1)
     rhs(i) = rhs(i)*dt/r(i)
  end do
end subroutine ca_rhstoer

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


