! begin photon routine

subroutine ca_accel_acoe( lo, hi,  &
     eta1,eta1_l1,eta1_l2,eta1_l3,eta1_h1,eta1_h2,eta1_h3, &
     spc , spc_l1, spc_l2, spc_l3, spc_h1, spc_h2, spc_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     aco , aco_l1, aco_l2, aco_l3, aco_h1, aco_h2, aco_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: eta1_l1,eta1_h1,eta1_l2,eta1_h2,eta1_l3,eta1_h3
  integer, intent(in) ::  spc_l1, spc_h1, spc_l2, spc_h2, spc_l3, spc_h3
  integer, intent(in) ::  kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer, intent(in) ::  aco_l1, aco_h1, aco_l2, aco_h2, aco_l3, aco_h3
  real(rt)        , intent(in)::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2,eta1_l3:eta1_h3)
  real(rt)        , intent(in)::spc ( spc_l1: spc_h1, spc_l2: spc_h2, spc_l3: spc_h3,0:ngroups-1)
  real(rt)        , intent(in)::kap ( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,0:ngroups-1)
  real(rt)                    ::aco ( aco_l1: aco_h1, aco_l2: aco_h2, aco_l3: aco_h3)
  real(rt)        , intent(in)::dt, tau

  integer :: i, j, k
  real(rt)         :: kbar, H1, dt1

  dt1 = (1.e0_rt+tau)/dt

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     kbar = sum(spc(i,j,k,:) * kap(i,j,k,:))
     H1 = eta1(i,j,k)
     aco(i,j,k) = H1*kbar*clight + dt1
  end do
  end do
  end do

end subroutine ca_accel_acoe


subroutine ca_accel_rhs( lo, hi, &
     Ern , Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3, &
     Erl , Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     rhs , rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3, &
     dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer, intent(in) :: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer, intent(in) :: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer, intent(in) ::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer, intent(in) :: rhs_l1, rhs_h1, rhs_l2, rhs_h2, rhs_l3, rhs_h3
  real(rt)        ,intent(in)::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,0:ngroups-1)
  real(rt)        ,intent(in)::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,0:ngroups-1)
  real(rt)        ,intent(in):: kap( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,0:ngroups-1)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  real(rt)                   :: rhs( rhs_l1: rhs_h1, rhs_l2: rhs_h2, rhs_l3: rhs_h3)
  real(rt)        ,intent(in) :: dt

  integer :: i, j, k
  real(rt)         :: rt_term, H

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rt_term = sum(kap(i,j,k,:)*(Ern(i,j,k,:)-Erl(i,j,k,:)))
     H = etaT(i,j,k)
     rhs(i,j,k) = clight*H*rt_term
  end do
  end do
  end do

end subroutine ca_accel_rhs


subroutine ca_accel_spec(lo, hi, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  integer,intent(in)::spec_l1,spec_h1,spec_l2,spec_h2,spec_l3,spec_h3
  real(rt)        ,intent(in)::kap ( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,0:ngroups-1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,0:ngroups-1)
  real(rt)                   ::spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,0:ngroups-1)
  real(rt)        ,intent(in):: dt, tau

  integer :: i, j, k
  real(rt)         :: cdt1, sumeps
  real(rt)        ,dimension(0:ngroups-1):: epsilon, kapt

  cdt1 = 1.e0_rt/(clight*dt)

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     kapt = kap(i,j,k,:) + (1.e0_rt+tau)*cdt1
     epsilon = mugT(i,j,k,:) / kapt
     sumeps = sum(epsilon)
     if (sumeps .eq. 0.e0_rt) then
        spec(i,j,k,:) = 0.e0_rt
     else
        spec(i,j,k,:) = epsilon / sumeps
     end if
  end do
  end do
  end do
  
end subroutine ca_accel_spec


subroutine ca_compute_rhs( lo, hi, &
     rhs , rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3, &
     jg  ,  jg_l1,  jg_l2,  jg_l3,  jg_h1,  jg_h2,  jg_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     cpT , cpT_l1, cpT_l2, cpT_l3, cpT_h1, cpT_h2, cpT_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     Er2 , Er2_l1, Er2_l2, Er2_l3, Er2_h1, Er2_h2, Er2_h3, &
     re2 , re2_l1, re2_l2, re2_l3, re2_h1, re2_h2, re2_h3, &
     Ers , Ers_l1, Ers_l2, Ers_l3, Ers_h1, Ers_h2, Ers_h3, &
     res , res_l1, res_l2, res_l3, res_h1, res_h2, res_h3, &
     r, dt, igroup, tau) bind(C, name="ca_compute_rhs")

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(3), hi(3) 
  integer,intent(in):: rhs_l1, rhs_h1, rhs_l2, rhs_h2, rhs_l3, rhs_h3
  integer,intent(in)::  jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  integer,intent(in):: cpT_l1, cpT_h1, cpT_l2, cpT_h2, cpT_l3, cpT_h3
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer,intent(in):: Er2_l1, Er2_h1, Er2_l2, Er2_h2, Er2_l3, Er2_h3
  integer,intent(in):: re2_l1, re2_h1, re2_l2, re2_h2, re2_l3, re2_h3
  integer,intent(in):: Ers_l1, Ers_h1, Ers_l2, Ers_h2, Ers_l3, Ers_h3
  integer,intent(in):: res_l1, res_h1, res_l2, res_h2, res_l3, res_h3
  real(rt)                   ::rhs ( rhs_l1: rhs_h1, rhs_l2: rhs_h2, rhs_l3: rhs_h3)
  real(rt)        ,intent(in)::jg  (  jg_l1:  jg_h1,  jg_l2:  jg_h2,  jg_l3:  jg_h3,0:ngroups-1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::cpT ( cpT_l1: cpT_h1, cpT_l2: cpT_h2, cpT_l3: cpT_h3)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  real(rt)        ,intent(in)::Er2 ( Er2_l1: Er2_h1, Er2_l2: Er2_h2, Er2_l3: Er2_h3,0:ngroups-1)
  real(rt)        ,intent(in)::re2 ( re2_l1: re2_h1, re2_l2: re2_h2, re2_l3: re2_h3)
  real(rt)        ,intent(in)::Ers ( Ers_l1: Ers_h1, Ers_l2: Ers_h2, Ers_l3: Ers_h3,0:ngroups-1)
  real(rt)        ,intent(in)::res ( res_l1: res_h1, res_l2: res_h2, res_l3: res_h3)
  real(rt)        ,intent(in):: r(lo(1):hi(1))
  real(rt)        ,intent(in):: dt, tau
  integer, intent(in) :: igroup

  integer :: i, j, k
  real(rt)         :: Hg, dt1

  dt1 = 1.e0_rt/dt
  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     Hg = mugT(i,j,k,igroup) * etaT(i,j,k)
     rhs(i,j,k) = clight*(jg(i,j,k,igroup) + Hg*cpT(i,j,k))  &
          + dt1 * (Er2(i,j,k,igroup) - Hg*(res(i,j,k)-re2(i,j,k)) &
          &        + tau*Ers(i,j,k,igroup))
   end do
   end do
   end do

end subroutine ca_compute_rhs


subroutine ca_compute_rhs_so( lo, hi, & ! MG Su-Olson
     rhs , rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3, &
     jg  ,  jg_l1,  jg_l2,  jg_l3,  jg_h1,  jg_h2,  jg_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     cpt , cpt_l1, cpt_l2, cpt_l3, cpt_h1, cpt_h2, cpt_h3, &
     eta , eta_l1, eta_l2, eta_l3, eta_h1, eta_h2, eta_h3, &
     Er2 , Er2_l1, Er2_l2, Er2_l3, Er2_h1, Er2_h2, Er2_h3, &
     re2 , re2_l1, re2_l2, re2_l3, re2_h1, re2_h2, re2_h3, &
     res , res_l1, res_l2, res_l3, res_h1, res_h2, res_h3, &
     x, t, dt, igroup) bind(C, name="ca_compute_rhs_so")

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(3), hi(3) 
  integer,intent(in):: rhs_l1, rhs_h1, rhs_l2, rhs_h2, rhs_l3, rhs_h3
  integer,intent(in)::  jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  integer,intent(in):: cpt_l1, cpt_h1, cpt_l2, cpt_h2, cpt_l3, cpt_h3
  integer,intent(in):: eta_l1, eta_h1, eta_l2, eta_h2, eta_l3, eta_h3
  integer,intent(in):: Er2_l1, Er2_h1, Er2_l2, Er2_h2, Er2_l3, Er2_h3
  integer,intent(in):: re2_l1, re2_h1, re2_l2, re2_h2, re2_l3, re2_h3
  integer,intent(in):: res_l1, res_h1, res_l2, res_h2, res_l3, res_h3
  real(rt)                  ::rhs ( rhs_l1: rhs_h1, rhs_l2: rhs_h2, rhs_l3: rhs_h3)
  real(rt)        ,intent(in)::jg  (  jg_l1:  jg_h1,  jg_l2:  jg_h2,  jg_l3:  jg_h3,0:ngroups-1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::cpt ( cpt_l1: cpt_h1, cpt_l2: cpt_h2, cpt_l3: cpt_h3)
  real(rt)        ,intent(in)::eta ( eta_l1: eta_h1, eta_l2: eta_h2, eta_l3: eta_h3)
  real(rt)        ,intent(in)::Er2 ( Er2_l1: Er2_h1, Er2_l2: Er2_h2, Er2_l3: Er2_h3,0:ngroups-1)
  real(rt)        ,intent(in)::re2 ( re2_l1: re2_h1, re2_l2: re2_h2, re2_l3: re2_h3)
  real(rt)        ,intent(in)::res ( res_l1: res_h1, res_l2: res_h2, res_l3: res_h3)
  real(rt)        ,intent(in):: x(lo(1):hi(1))
  real(rt)        ,intent(in):: t, dt
  integer, intent(in) :: igroup

  real(rt)        , parameter :: x0 = 0.5e0_rt
  real(rt)        , parameter :: t0 = 3.3356409519815202e-10_rt 
  real(rt)        , parameter :: qn = 1.134074546528399e20_rt 

  integer :: i, j, k
  real(rt)         :: Hg

  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     Hg = mugT(i,j,k,igroup)*eta(i,j,k)
     rhs(i,j,k) = clight*jg(i,j,k,igroup) + clight*cpt(i,j,k)*Hg &
          + (Er2(i,j,k,igroup) - (res(i,j,k)-re2(i,j,k))*Hg) / dt
     if (t .le. t0 .and. abs(x(i)) .le. x0) then
        rhs(i,j,k) = rhs(i,j,k) + qn ! (qn / dt) * dt
     end if
   end do
   end do
   end do

end subroutine ca_compute_rhs_so


subroutine ca_local_accel( lo, hi,  &
     Ern , Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3, &
     Erl , Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  real(rt)                   ::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,0:ngroups-1)
  real(rt)        ,intent(in)::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,0:ngroups-1)
  real(rt)        ,intent(in)::kap ( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,0:ngroups-1)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,0:ngroups-1)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i, j, k
  real(rt)         :: cdt1, rt_term, p
  real(rt)        ,dimension(0:ngroups-1)::Hg, epsilon, kapt, kk

  cdt1 = 1.e0_rt/(clight*dt)

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rt_term = sum(kap(i,j,k,:)*(Ern(i,j,k,:)-Erl(i,j,k,:)))

     Hg = mugT(i,j,k,:)*etaT(i,j,k)

     kapt = kap(i,j,k,:) + (1.e0_rt+tau)*cdt1
     kk = kap(i,j,k,:) / kapt

     p = 1.e0_rt-sum(Hg*kk)
     epsilon = (Hg * rt_term) / (kapt*p + 1.e-50_rt)

     Ern(i,j,k,:) = Ern(i,j,k,:) + epsilon
  end do
  end do
  end do

end subroutine ca_local_accel


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


subroutine ca_ncupdate_matter( lo, hi,  &
     Tp_n,Tp_n_l1,Tp_n_l2,Tp_n_l3,Tp_n_h1,Tp_n_h2,Tp_n_h3,  &
     Er_n,Er_n_l1,Er_n_l2,Er_n_l3,Er_n_h1,Er_n_h2,Er_n_h3,  &
     re_s,re_s_l1,re_s_l2,re_s_l3,re_s_h1,re_s_h2,re_s_h3,  &
     re_2,re_2_l1,re_2_l2,re_2_l3,re_2_h1,re_2_h2,re_2_h3,  &
     etTz,etTz_l1,etTz_l2,etTz_l3,etTz_h1,etTz_h2,etTz_h3,  &
      kpp, kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3,  &
       jg,  jg_l1,  jg_l2,  jg_l3,  jg_h1,  jg_h2,  jg_h3,  &
     dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(3),hi(3)
  integer,intent(in)::Tp_n_l1,Tp_n_h1,Tp_n_l2,Tp_n_h2,Tp_n_l3,Tp_n_h3
  integer,intent(in)::Er_n_l1,Er_n_h1,Er_n_l2,Er_n_h2,Er_n_l3,Er_n_h3
  integer,intent(in)::re_s_l1,re_s_h1,re_s_l2,re_s_h2,re_s_l3,re_s_h3
  integer,intent(in)::re_2_l1,re_2_h1,re_2_l2,re_2_h2,re_2_l3,re_2_h3
  integer,intent(in)::etTz_l1,etTz_h1,etTz_l2,etTz_h2,etTz_l3,etTz_h3
  integer,intent(in):: kpp_l1, kpp_h1, kpp_l2, kpp_h2, kpp_l3, kpp_h3
  integer,intent(in)::  jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  real(rt)                   ::Tp_n(Tp_n_l1:Tp_n_h1,Tp_n_l2:Tp_n_h2,Tp_n_l3:Tp_n_h3)
  real(rt)        ,intent(in)::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,Er_n_l3:Er_n_h3,0:ngroups-1)
  real(rt)        ,intent(in)::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2,re_s_l3:re_s_h3)
  real(rt)        ,intent(in)::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2,re_2_l3:re_2_h3)
  real(rt)        ,intent(in)::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2,etTz_l3:etTz_h3)
  real(rt)        ,intent(in):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,0:ngroups-1)
  real(rt)        ,intent(in)::  jg(  jg_l1:  jg_h1,  jg_l2:  jg_h2,  jg_l3:  jg_h3,0:ngroups-1)
  real(rt)        ,intent(in) :: dt

   integer :: i,j,k,g
   real(rt)         :: cdt1, cpT, scrch_re
   real(rt)         :: dTemp
   real(rt)        , parameter :: fac = 0.01e0_rt

   cdt1 = 1.e0_rt / (clight * dt)
   do k = lo(3), hi(3)
   do j = lo(2), hi(2)
   do i = lo(1), hi(1)

      cpT = 0.e0_rt
      do g = 0, ngroups-1
         cpT = cpT + kpp(i,j,k,g)*Er_n(i,j,k,g) - jg(i,j,k,g)
      end do

      scrch_re = cpT - (re_s(i,j,k) - re_2(i,j,k)) * cdt1

      dTemp = etTz(i,j,k)*scrch_re

      if (abs(dTemp/(Tp_n(i,j,k)+1.e-50_rt)) > fac) then
         dTemp = sign(fac*Tp_n(i,j,k), dTemp)
      end if

     Tp_n(i,j,k) = Tp_n(i,j,k) + dTemp

  end do
  end do
  end do

end subroutine ca_ncupdate_matter



! end photon routines
! ========================================================================

! other routines that work for both photon and neutrinos

subroutine ca_accel_ccoe( lo, hi, &
     bcgr,bcgr_l1,bcgr_l2,bcgr_l3,bcgr_h1,bcgr_h2,bcgr_h3, &
     spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3, &
     ccoe,ccoe_l1,ccoe_l2,ccoe_l3,ccoe_h1,ccoe_h2,ccoe_h3, &
     dx, idim, igroup)

  use rad_params_module, only : ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: bcgr_l1, bcgr_h1, bcgr_l2, bcgr_h2, bcgr_l3, bcgr_h3
  integer, intent(in) :: spec_l1, spec_h1, spec_l2, spec_h2, spec_l3, spec_h3
  integer, intent(in) :: ccoe_l1, ccoe_h1, ccoe_l2, ccoe_h2, ccoe_l3, ccoe_h3
  real(rt)        ,intent(in)::bcgr(bcgr_l1:bcgr_h1,bcgr_l2:bcgr_h2,bcgr_l3:bcgr_h3)
  real(rt)        ,intent(in)::spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,0:ngroups-1)
  real(rt)                   ::ccoe(ccoe_l1:ccoe_h1,ccoe_l2:ccoe_h2,ccoe_l3:ccoe_h3,0:1)
  real(rt)        , intent(in) :: dx(3)
  integer, intent(in) :: idim, igroup

  integer :: i, j, k, ioff, joff, koff
  real(rt)         :: grad_spec, foo, h1

  if (idim .eq. 0) then
     ioff = 1
     joff = 0
     koff = 0
     h1 = 1.e0_rt/dx(1)
  else if (idim .eq. 1) then
     ioff = 0
     joff = 1
     koff = 0
     h1 = 1.e0_rt/dx(2)
  else
     ioff = 0
     joff = 0
     koff = 1
     h1 = 1.e0_rt/dx(3)
  end if

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     grad_spec = (spec(i,j,k,igroup) - spec(i-ioff,j-joff,k-koff,igroup)) * h1
     foo = - 0.5e0_rt * bcgr(i,j,k) * grad_spec
     ccoe(i,j,k,0) = ccoe(i,j,k,0) + foo
     ccoe(i,j,k,1) = ccoe(i,j,k,1) + foo
  end do
  end do
  end do

end subroutine ca_accel_ccoe


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

subroutine ca_rhstoer( lo, hi, &
     rhs, rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3, &
     r, dt)
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: rhs_l1, rhs_h1, rhs_l2, rhs_h2, rhs_l3, rhs_h3
  real(rt)         ::rhs ( rhs_l1: rhs_h1, rhs_l2: rhs_h2, rhs_l3: rhs_h3)
  real(rt)        ,intent(in):: r(lo(1):hi(1))
  real(rt)        ,intent(in):: dt
  integer :: i, j, k
  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     rhs(i,j,k) = rhs(i,j,k)*dt
  end do
  end do
  end do
end subroutine ca_rhstoer

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

