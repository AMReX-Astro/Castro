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


subroutine ca_check_conv( lo, hi, &
     ren,ren_l1,ren_l2,ren_l3,ren_h1,ren_h2,ren_h3, &
     res,res_l1,res_l2,res_l3,res_h1,res_h2,res_h3, &
     re2,re2_l1,re2_l2,re2_l3,re2_h1,re2_h2,re2_h3, &
     ern,ern_l1,ern_l2,ern_l3,ern_h1,ern_h2,ern_h3, &
     Tmn,Tmn_l1,Tmn_l2,Tmn_l3,Tmn_h1,Tmn_h2,Tmn_h3, &
     Tms,Tms_l1,Tms_l2,Tms_l3,Tms_h1,Tms_h2,Tms_h3, &
     rho,rho_l1,rho_l2,rho_l3,rho_h1,rho_h2,rho_h3, &
     kap,kap_l1,kap_l2,kap_l3,kap_h1,kap_h2,kap_h3, &
     jg , jg_l1, jg_l2, jg_l3, jg_h1, jg_h2, jg_h3, &
     deT,deT_l1,deT_l2,deT_l3,deT_h1,deT_h2,deT_h3, &
     rel_re, abs_re, & 
     rel_FT, abs_FT, rel_T, abs_T, &
     dt)
  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(3),hi(3)
  integer,intent(in)::ren_l1, ren_h1, ren_l2, ren_h2, ren_l3, ren_h3
  integer,intent(in)::res_l1, res_h1, res_l2, res_h2, res_l3, res_h3
  integer,intent(in)::re2_l1, re2_h1, re2_l2, re2_h2, re2_l3, re2_h3
  integer,intent(in)::ern_l1, ern_h1, ern_l2, ern_h2, ern_l3, ern_h3
  integer,intent(in)::Tmn_l1, Tmn_h1, Tmn_l2, Tmn_h2, Tmn_l3, Tmn_h3
  integer,intent(in)::Tms_l1, Tms_h1, Tms_l2, Tms_h2, Tms_l3, Tms_h3
  integer,intent(in)::rho_l1, rho_h1, rho_l2, rho_h2, rho_l3, rho_h3
  integer,intent(in)::kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in):: jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  integer,intent(in)::deT_l1, deT_h1, deT_l2, deT_h2, deT_l3, deT_h3
  real(rt)        ,intent(in)::ren(ren_l1:ren_h1,ren_l2:ren_h2,ren_l3:ren_h3)
  real(rt)        ,intent(in)::res(res_l1:res_h1,res_l2:res_h2,res_l3:res_h3)
  real(rt)        ,intent(in)::re2(re2_l1:re2_h1,re2_l2:re2_h2,re2_l3:re2_h3)
  real(rt)        ,intent(in)::ern(ern_l1:ern_h1,ern_l2:ern_h2,ern_l3:ern_h3,0:ngroups-1)
  real(rt)        ,intent(in)::Tmn(Tmn_l1:Tmn_h1,Tmn_l2:Tmn_h2,Tmn_l3:Tmn_h3)
  real(rt)        ,intent(in)::Tms(Tms_l1:Tms_h1,Tms_l2:Tms_h2,Tms_l3:Tms_h3)
  real(rt)        ,intent(in)::rho(rho_l1:rho_h1,rho_l2:rho_h2,rho_l3:rho_h3)
  real(rt)        ,intent(in)::kap(kap_l1:kap_h1,kap_l2:kap_h2,kap_l3:kap_h3,0:ngroups-1)
  real(rt)        ,intent(in):: jg( jg_l1: jg_h1, jg_l2: jg_h2, jg_l3: jg_h3,0:ngroups-1)
  real(rt)        ,intent(in)::deT(deT_l1:deT_h1,deT_l2:deT_h2,deT_l3:deT_h3)
  real(rt)        ,intent(inout)::rel_re, abs_re 
  real(rt)        ,intent(inout)::rel_FT, abs_FT,rel_T, abs_T
  real(rt)        ,intent(in) :: dt

  integer :: i, j, k
  real(rt)         :: chg, relchg, FT, cdt, FTdenom, dTe

  cdt = clight*dt

  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
  do i=lo(1),hi(1)
     chg = abs(ren(i,j,k) - res(i,j,k))
     relchg = abs(chg/(ren(i,j,k)+1.e-50_rt))
     rel_re = max(rel_re,relchg)
     abs_re = max(abs_re,chg)

     chg = abs(Tmn(i,j,k) - Tms(i,j,k))
     relchg = abs(chg/(Tmn(i,j,k)+1.e-50_rt))
     rel_T = max(rel_T,relchg)
     abs_T = max(abs_T,chg)

     FT = abs((ren(i,j,k)-re2(i,j,k)) - cdt*sum(kap(i,j,k,:)*Ern(i,j,k,:)-jg(i,j,k,:)))

     dTe = Tmn(i,j,k)
     FTdenom = rho(i,j,k)*abs(deT(i,j,k)*dTe)
!     FTdenom = max(abs(ren(i,j,k)-re2(i,j,k)), abs(ren(i,j,k)*1.e-15_rt))

     rel_FT = max(rel_FT, FT/(FTdenom+1.e-50_rt))
     abs_FT = max(abs_FT, FT)
  end do
  end do
  end do

end subroutine ca_check_conv


subroutine ca_check_conv_er( lo, hi, &
     Ern,  Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3, &
     Erl,  Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3, &
     kap,  kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     etTz,etTz_l1,etTz_l2,etTz_l3,etTz_h1,etTz_h2,etTz_h3, &
     temp,temp_l1,temp_l2,temp_l3,temp_h1,temp_h2,temp_h3, &
     rela, abso, errr, dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in)::temp_l1,temp_h1,temp_l2,temp_h2,temp_l3,temp_h3
  integer,intent(in)::etTz_l1,etTz_h1,etTz_l2,etTz_h2,etTz_l3,etTz_h3
  real(rt)        ,intent(in):: Ern( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,0:ngroups-1)
  real(rt)        ,intent(in):: Erl( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,0:ngroups-1)
  real(rt)        ,intent(in):: kap( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,0:ngroups-1)
  real(rt)        ,intent(in)::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2,etTz_l3:etTz_h3)
  real(rt)        ,intent(in)::temp(temp_l1:temp_h1,temp_l2:temp_h2,temp_l3:temp_h3)
  real(rt)        , intent(inout) :: rela, abso, errr
  real(rt)        , intent(in) :: dt

  integer :: i, j, k, g
  real(rt)         :: chg, tot, cdt, der, kde, err_T, err

  cdt = clight * dt
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     chg = 0.e0_rt
     tot = 0.e0_rt
     kde = 0.e0_rt
     do g=0,ngroups-1
        der = Ern(i,j,k,g)-Erl(i,j,k,g)
        chg = chg + abs(der)
        tot = tot + abs(Ern(i,j,k,g))
        kde = kde + kap(i,j,k,g)*der
     end do
     abso = max(abso, chg)
     rela = max(rela, chg / (tot + 1.e-50_rt))

     err_T =  etTz(i,j,k)*kde
     err = abs(err_T/(temp(i,j,k)+1.e-50_rt))
     errr = max(errr, err)
  end do
  end do
  end do

end subroutine ca_check_conv_er


subroutine ca_compute_coupt( lo, hi,  &
     cpt, cpt_l1, cpt_l2, cpt_l3, cpt_h1, cpt_h2, cpt_h3, &
     kpp, kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     eg ,  eg_l1,  eg_l2,  eg_l3,  eg_h1,  eg_h2,  eg_h3, &
     jg ,  jg_l1,  jg_l2,  jg_l3,  jg_h1,  jg_h2,  jg_h3)
  
  use rad_params_module, only : ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: cpt_l1, cpt_h1, cpt_l2, cpt_h2, cpt_l3, cpt_h3 
  integer, intent(in) :: kpp_l1, kpp_h1, kpp_l2, kpp_h2, kpp_l3, kpp_h3 
  integer, intent(in) ::  eg_l1,  eg_h1,  eg_l2,  eg_h2,  eg_l3,  eg_h3
  integer, intent(in) ::  jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  real(rt)                   ::cpt(cpt_l1:cpt_h1,cpt_l2:cpt_h2,cpt_l3:cpt_h3)
  real(rt)        ,intent(in)::kpp(kpp_l1:kpp_h1,kpp_l2:kpp_h2,kpp_l3:kpp_h3,0:ngroups-1)
  real(rt)        ,intent(in):: eg( eg_l1: eg_h1, eg_l2: eg_h2, eg_l3: eg_h3,0:ngroups-1)
  real(rt)        ,intent(in):: jg( jg_l1: jg_h1, jg_l2: jg_h2, jg_l3: jg_h3,0:ngroups-1)

  integer :: i, j, k, g

  cpt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.e0_rt

  do g=0, ngroups-1
     do k=lo(3),hi(3)
     do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        cpt(i,j,k) = cpt(i,j,k) + (kpp(i,j,k,g) * eg(i,j,k,g) - jg(i,j,k,g))
     end do
     end do
     end do
  end do

end subroutine ca_compute_coupt


subroutine ca_compute_etat( lo, hi, &
     & etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     & etTz,etTz_l1,etTz_l2,etTz_l3,etTz_h1,etTz_h2,etTz_h3, &
     & eta1,eta1_l1,eta1_l2,eta1_l3,eta1_h1,eta1_h2,eta1_h3, &
     & djdT,djdT_l1,djdT_l2,djdT_l3,djdT_h1,djdT_h2,djdT_h3, &
     & dkdT,dkdT_l1,dkdT_l2,dkdT_l3,dkdT_h1,dkdT_h2,dkdT_h3, &
     & dedT,dedT_l1,dedT_l2,dedT_l3,dedT_h1,dedT_h2,dedT_h3, &
     & Ers , Ers_l1, Ers_l2, Ers_l3, Ers_h1, Ers_h2, Ers_h3, &
     & rho , rho_l1, rho_l2, rho_l3, rho_h1, rho_h2, rho_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: etaT_l1, etaT_h1, etaT_l2, etaT_h2, etaT_l3, etaT_h3
  integer, intent(in) :: etTz_l1, etTz_h1, etTz_l2, etTz_h2, etTz_l3, etTz_h3
  integer, intent(in) :: eta1_l1, eta1_h1, eta1_l2, eta1_h2, eta1_l3, eta1_h3
  integer, intent(in) :: djdT_l1, djdT_h1, djdT_l2, djdT_h2, djdT_l3, djdT_h3
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2, dkdT_l3, dkdT_h3
  integer, intent(in) :: dedT_l1, dedT_h1, dedT_l2, dedT_h2, dedT_l3, dedT_h3
  integer, intent(in) ::  Ers_l1,  Ers_h1,  Ers_l2,  Ers_h2,  Ers_l3,  Ers_h3
  integer, intent(in) ::  rho_l1,  rho_h1,  rho_l2,  rho_h2,  rho_l3,  rho_h3
  real(rt)                   ::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  real(rt)                   ::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2,etTz_l3:etTz_h3)
  real(rt)                   ::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2,eta1_l3:eta1_h3)
  real(rt)                   ::djdT(djdT_l1:djdT_h1,djdT_l2:djdT_h2,djdT_l3:djdT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,dkdT_l3:dkdT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::dedT(dedT_l1:dedT_h1,dedT_l2:dedT_h2,dedT_l3:dedT_h3)
  real(rt)        ,intent(in)::Ers ( Ers_l1: Ers_h1, Ers_l2: Ers_h2, Ers_l3: Ers_h3,0:ngroups-1)
  real(rt)        ,intent(in)::rho ( rho_l1: rho_h1, rho_l2: rho_h2, rho_l3: rho_h3)
  real(rt)        ,intent(in):: dt, tau

  integer :: i, j, k
  real(rt)         :: cdt, sigma
  real(rt)         :: dZdT(0:ngroups-1), sumdZdT, foo, bar

  sigma = 1.e0_rt + tau
  cdt = clight * dt

  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     dZdT = djdT(i,j,k,:) - dkdT(i,j,k,:)*Ers(i,j,k,:)
     sumdZdT = sum(dZdT)
     if (sumdZdT .eq. 0.e0_rt) then
        sumdZdT = 1.e-50_rt
     end if
     foo = cdt * sumdZdT
     bar = sigma*rho(i,j,k)*dedT(i,j,k)
     etaT(i,j,k) = foo / (foo + bar)
     etTz(i,j,k) = etaT(i,j,k) / sumdZdT
     eta1(i,j,k) = bar / (foo + bar)
     djdT(i,j,k,:) = dZdT / sumdZdT
  end do
  end do
  end do

end subroutine ca_compute_etat


subroutine ca_compute_emissivity( lo, hi, &
     jg  ,  jg_l1,  jg_l2,  jg_l3,  jg_h1,  jg_h2,  jg_h3, &
     djdT,djdT_l1,djdT_l2,djdT_l3,djdT_h1,djdT_h2,djdT_h3, &
        T,   T_l1,   T_l2,   T_l3,   T_h1,   T_h2,   T_h3, &
      kap, kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     dkdT,dkdT_l1,dkdT_l2,dkdT_l3,dkdT_h1,dkdT_h2,dkdT_h3, &
     pfc, use_WiensLaw, integrate_Planck, Tf)

  use rad_params_module, only : ngroups, nugroup, dnugroup, xnu,  &
       pi, clight, hplanck, kboltz, arad
  use blackbody_module, only : BdBdTIndefInteg

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in)  :: lo(3), hi(3) 
  integer, intent(in) ::   jg_l1,   jg_h1,   jg_l2,   jg_h2,   jg_l3,   jg_h3
  integer, intent(in) :: djdT_l1, djdT_h1, djdT_l2, djdT_h2, djdT_l3, djdT_h3
  integer, intent(in) ::    T_l1,    T_h1,    T_l2,    T_h2,    T_l3,    T_h3
  integer, intent(in) ::  kap_l1,  kap_h1,  kap_l2,  kap_h2,  kap_l3,  kap_h3
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2, dkdT_l3, dkdT_h3
  real(rt)                   ::jg  (  jg_l1:  jg_h1,  jg_l2:  jg_h2,  jg_l3:  jg_h3,0:ngroups-1)
  real(rt)                   ::djdT(djdT_l1:djdT_h1,djdT_l2:djdT_h2,djdT_l3:djdT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::   T(   T_l1:   T_h1,   T_l2:   T_h2,   T_l3:   T_h3)
  real(rt)        ,intent(in):: kap( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,0:ngroups-1)
  real(rt)        ,intent(in)::dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,dkdT_l3:dkdT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::pfc(0:ngroups-1)
  integer, intent(in) :: use_WiensLaw, integrate_Planck
  real(rt)        , intent(in) :: Tf

  integer :: i, j, k, g
  real(rt)         :: dBdT, Bg
  real(rt)         :: Teff, nu, num, nup, hoverk
  real(rt)         :: cB, Tfix
  real(rt)         :: B0, B1, dBdT0, dBdT1
  real(rt)         :: dnu, nubar, expnubar, cdBdT
  real(rt)         :: xnu_full(0:ngroups)

  if (ngroups .eq. 1) then

     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        Bg = arad*T(i,j,k)**4
        dBdT = 4.e0_rt*arad*T(i,j,k)**3
        g = 0
        jg(i,j,k,g) = Bg*kap(i,j,k,g)
        djdT(i,j,k,g) = dkdT(i,j,k,g)*Bg + dBdT*kap(i,j,k,g)
     end do     
     end do     
     end do     

  else if (pfc(0) > 0.e0_rt) then  ! a special case for picket-fence model in Su-Olson test
     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        Bg = arad*T(i,j,k)**4
        dBdT = 4.e0_rt*arad*T(i,j,k)**3
        do g=0, ngroups-1
           jg(i,j,k,g) = pfc(g) * Bg*kap(i,j,k,g)
           djdT(i,j,k,g) = pfc(g) * (dkdT(i,j,k,g)*Bg + dBdT*kap(i,j,k,g))
        end do
     end do
     end do
     end do
  else if (use_WiensLaw > 0) then

     hoverk = hplanck/kboltz
     cB = 8.*pi*kboltz/clight**3
     
     do g=0, ngroups-1
        nu = nugroup(g)
        num = xnu(g)
        nup = xnu(g+1)
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (Tf < 0.e0_rt) then
              Tfix = T(i,j,k)
           else
              Tfix = Tf
           end if
           dBdT = cB * nu**3 * (exp(-hoverk*num/Tfix) - exp(-hoverk*nup/Tfix))
           Bg = dBdT * T(i,j,k)
           jg(i,j,k,g) = Bg*kap(i,j,k,g)
           djdT(i,j,k,g) = dkdT(i,j,k,g)*Bg + dBdT*kap(i,j,k,g)
        end do
        end do
        end do
     end do

  else if (integrate_Planck > 0) then

     xnu_full = xnu(0:ngroups)
     xnu_full(0) = 0.e0_rt
     xnu_full(ngroups) = max(xnu(ngroups), 1.e25_rt)

     do k=lo(3), hi(3)
     do j=lo(2), hi(2)
     do i=lo(1), hi(1)
        Teff = max(T(i,j,k), 1.e-50_rt)
        call BdBdTIndefInteg(Teff, xnu_full(0), B1, dBdT1)
        do g=0, ngroups-1
           B0 = B1
           dBdT0 = dBdT1
           call BdBdTIndefInteg(Teff, xnu_full(g+1), B1, dBdT1)
           Bg = B1 - B0
           dBdT = dBdT1 - dBdT0

           jg(i,j,k,g) = Bg*kap(i,j,k,g)
           djdT(i,j,k,g) = dkdT(i,j,k,g)*Bg + dBdT*kap(i,j,k,g)
        end do
     end do
     end do
     end do

  else

     cB = 8.e0_rt*pi*hplanck / clight**3
     cdBdT = 8.e0_rt*pi*hplanck**2 / (kboltz*clight**3)

     do g=0, ngroups-1
        nu = nugroup(g)
        dnu = dnugroup(g)

        do k=lo(3), hi(3)
        do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           Teff = max(T(i,j,k), 1.e-50_rt)
           nubar = hplanck * nu / (kboltz * Teff)
           if (nubar > 100.e0_rt) then
              Bg = 0.e0_rt
              dBdT = 0.e0_rt
           else if (nubar < 1.e-15_rt) then
              Bg = 0.e0_rt
              dBdT = 0.e0_rt           
           else
              expnubar = exp(nubar)
              Bg = cB * nu**3 / (expnubar - 1.e0_rt) * dnu
              dBdT = cdBdT * nu**4 / Teff**2 * expnubar / (expnubar-1.e0_rt)**2 * dnu
           end if

           jg(i,j,k,g) = Bg*kap(i,j,k,g)
           djdT(i,j,k,g) = dkdT(i,j,k,g)*Bg + dBdT*kap(i,j,k,g)
        end do
        end do
        end do
     end do
     
  end if

end subroutine ca_compute_emissivity


subroutine ca_compute_kappas(lo, hi, &
     stt,stt_l1,stt_l2,stt_l3,stt_h1,stt_h2,stt_h3, &
     T  ,  T_l1,  T_l2,  T_l3,  T_h1,  T_h2,  T_h3, &
     kpp,kpp_l1,kpp_l2,kpp_l3,kpp_h1,kpp_h2,kpp_h3, &
     kpr,kpr_l1,kpr_l2,kpr_l3,kpr_h1,kpr_h2,kpr_h3, &
     kpT,kpT_l1,kpT_l2,kpT_l3,kpT_h1,kpT_h2,kpT_h3, &
     do_stme, use_dkdT, &
     c_kpp, kpp_m, kpp_n, kpp_p, &
     c_kpr, kpr_m, kpr_n, kpr_p, &
     c_sct, sct_m, sct_n, sct_p, &
     Tfloor)

  use rad_params_module, only : ngroups, nugroup
  use fundamental_constants_module, only : hplanck, k_B
  use meth_params_module, only : NVAR, URHO

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3) 
  integer, intent(in) :: stt_l1, stt_h1, stt_l2, stt_h2, stt_l3, stt_h3 
  integer, intent(in) ::   T_l1,   T_h1,   T_l2,   T_h2,   T_l3,   T_h3
  integer, intent(in) :: kpp_l1, kpp_h1, kpp_l2, kpp_h2, kpp_l3, kpp_h3 
  integer, intent(in) :: kpr_l1, kpr_h1, kpr_l2, kpr_h2, kpr_l3, kpr_h3
  integer, intent(in) :: kpT_l1, kpT_h1, kpT_l2, kpT_h2, kpT_l3, kpT_h3
  integer, intent(in) :: do_stme, use_dkdT
  real(rt)        ,intent(in)::stt(stt_l1:stt_h1,stt_l2:stt_h2,stt_l3:stt_h3,NVAR)
  real(rt)        ,intent(in)::  T(  T_l1:  T_h1,  T_l2:  T_h2,  T_l3:  T_h3)
  real(rt)                   ::kpp(kpp_l1:kpp_h1,kpp_l2:kpp_h2,kpp_l3:kpp_h3,0:ngroups-1)
  real(rt)                   ::kpr(kpr_l1:kpr_h1,kpr_l2:kpr_h2,kpr_l3:kpr_h3,0:ngroups-1)
  real(rt)                   ::kpT(kpT_l1:kpT_h1,kpT_l2:kpT_h2,kpT_l3:kpT_h3,0:ngroups-1)
  real(rt)        , intent(in)   :: c_kpp, kpp_m, kpp_n, kpp_p 
  real(rt)        , intent(in)   :: c_kpr, kpr_m, kpr_n, kpr_p 
  real(rt)        , intent(in)   :: c_sct, sct_m, sct_n, sct_p 
  real(rt)        , intent(in)   :: Tfloor

  integer :: i, j, k, g
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, nup_kpp, nup_kpr, nup_sct, sct, foo, hnuoverkt, exptmp

  do g=0, ngroups-1
     nup_kpp = nugroup(g)**kpp_p
     nup_kpr = nugroup(g)**kpr_p
     nup_sct = nugroup(g)**sct_p

     do k=lo(3), hi(3)
     do j=lo(2), hi(2)
     do i=lo(1), hi(1)

        Teff = max(T(i,j,k), tiny)
        Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
        kpp(i,j,k,g) = c_kpp * (stt(i,j,k,URHO) ** kpp_m) * (Teff ** (-kpp_n)) * nup_kpp
        if (do_stme > 0) then
           foo = kpp(i,j,k,g)
           hnuoverkt = hplanck*nugroup(g)/(k_B*Teff)
           exptmp = exp(-hnuoverkt)
           kpp(i,j,k,g) = foo * (1.e0_rt - exptmp)
           kpT(i,j,k,g) = foo*(-kpp_n/Teff)*(1.e0_rt-exptmp) - foo*exptmp*(hnuoverkt/Teff)
        else
           kpT(i,j,k,g) = kpp(i,j,k,g) * (-kpp_n/Teff)           
        end if

        if (use_dkdT.eq.0) then
           kpT(i,j,k,g) = 0.e0_rt
        end if

        if (c_kpr < 0.0e0_rt) then
           sct       = c_sct * (stt(i,j,k,URHO) ** sct_m) * (Teff ** (-sct_n)) * nup_sct
           kpr(i,j,k,g) = kpp(i,j,k,g) + sct
        else
           kpr(i,j,k,g) = c_kpr * (stt(i,j,k,URHO) ** kpr_m) * (Teff ** (-kpr_n)) * nup_kpr
        end if
           
     end do
     end do
     end do
  end do

end subroutine ca_compute_kappas


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


subroutine ca_update_matter( lo, hi,  &
     re_n,re_n_l1,re_n_l2,re_n_l3,re_n_h1,re_n_h2,re_n_h3, &
     Er_n,Er_n_l1,Er_n_l2,Er_n_l3,Er_n_h1,Er_n_h2,Er_n_h3, &
     Er_l,Er_l_l1,Er_l_l2,Er_l_l3,Er_l_h1,Er_l_h2,Er_l_h3, &
     re_s,re_s_l1,re_s_l2,re_s_l3,re_s_h1,re_s_h2,re_s_h3, &
     re_2,re_2_l1,re_2_l2,re_2_l3,re_2_h1,re_2_h2,re_2_h3, &
     eta1,eta1_l1,eta1_l2,eta1_l3,eta1_h1,eta1_h2,eta1_h3, &
      cpt, cpt_l1, cpt_l2, cpt_l3, cpt_h1, cpt_h2, cpt_h3, &
      kpp, kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     Snew,Snew_l1,Snew_l2,Snew_l3,Snew_h1,Snew_h2,Snew_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight
  use meth_params_module, only : NVAR

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(3),hi(3)
  integer,intent(in)::re_n_l1, re_n_h1, re_n_l2, re_n_h2, re_n_l3, re_n_h3
  integer,intent(in)::Er_n_l1, Er_n_h1, Er_n_l2, Er_n_h2, Er_n_l3, Er_n_h3
  integer,intent(in)::Er_l_l1, Er_l_h1, Er_l_l2, Er_l_h2, Er_l_l3, Er_l_h3
  integer,intent(in)::re_s_l1, re_s_h1, re_s_l2, re_s_h2, re_s_l3, re_s_h3
  integer,intent(in)::re_2_l1, re_2_h1, re_2_l2, re_2_h2, re_2_l3, re_2_h3
  integer,intent(in)::eta1_l1, eta1_h1, eta1_l2, eta1_h2, eta1_l3, eta1_h3
  integer,intent(in):: cpt_l1,  cpt_h1,  cpt_l2,  cpt_h2,  cpt_l3,  cpt_h3
  integer,intent(in):: kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2,  kpp_l3,  kpp_h3
  integer,intent(in)::mugT_l1, mugT_h1, mugT_l2, mugT_h2, mugT_l3, mugT_h3
  integer,intent(in)::Snew_l1, Snew_h1, Snew_l2, Snew_h2, Snew_l3, Snew_h3
  real(rt)                   ::re_n(re_n_l1:re_n_h1,re_n_l2:re_n_h2,re_n_l3:re_n_h3)
  real(rt)        ,intent(in)::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,Er_n_l3:Er_n_h3,0:ngroups-1)
  real(rt)        ,intent(in)::Er_l(Er_l_l1:Er_l_h1,Er_l_l2:Er_l_h2,Er_l_l3:Er_l_h3,0:ngroups-1)
  real(rt)        ,intent(in)::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2,re_s_l3:re_s_h3)
  real(rt)        ,intent(in)::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2,re_2_l3:re_2_h3)
  real(rt)        ,intent(in)::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2,eta1_l3:eta1_h3)
  real(rt)        ,intent(in):: cpt( cpt_l1: cpt_h1, cpt_l2: cpt_h2, cpt_l3: cpt_h3)
  real(rt)        ,intent(in):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,0:ngroups-1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,0:ngroups-1)
  real(rt)        ,intent(in)::Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,Snew_l3:Snew_h3,NVAR)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i,j,k
  real(rt)         :: cdt, H1, dkEE, chg

  cdt = clight * dt
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     H1 = eta1(i,j,k)

     dkEE = sum(kpp(i,j,k,:)*(Er_n(i,j,k,:)-Er_l(i,j,k,:)))

     chg = cdt*dkEE + H1*((re_2(i,j,k)-re_s(i,j,k)) + cdt*cpt(i,j,k))

     re_n(i,j,k) = re_s(i,j,k) + chg

     re_n(i,j,k) = (re_n(i,j,k) + tau*re_s(i,j,k)) / (1.e0_rt+tau)

     ! temperature will be updated after exiting this subroutine
  end do
  end do
  end do

end subroutine ca_update_matter


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


subroutine ca_opacs( lo, hi,  &
     Snew,Snew_l1,Snew_l2,Snew_l3,Snew_h1,Snew_h2,Snew_h3, &
     T   ,   T_l1,   T_l2,   T_l3,   T_h1,   T_h2,   T_h3, &
     Ts  ,  Ts_l1,  Ts_l2,  Ts_l3,  Ts_h1,  Ts_h2,  Ts_h3, &
     kpp , kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     kpr , kpr_l1, kpr_l2, kpr_l3, kpr_h1, kpr_h2, kpr_h3, &
     dkdT,dkdT_l1,dkdT_l2,dkdT_l3,dkdT_h1,dkdT_h2,dkdT_h3, &
     use_dkdT, validStar, lag_opac) 

  use rad_params_module, only : ngroups, nugroup
  use opacity_table_module, only : get_opacities
  use network, only : naux
  use meth_params_module, only : NVAR, URHO, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: Snew_l1, Snew_h1, Snew_l2, Snew_h2, Snew_l3, Snew_h3 
  integer, intent(in) ::    T_l1,    T_h1,    T_l2,    T_h2,    T_l3,    T_h3
  integer, intent(in) ::   Ts_l1,   Ts_h1,   Ts_l2,   Ts_h2,   Ts_l3,   Ts_h3
  integer, intent(in) ::  kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2,  kpp_l3,  kpp_h3 
  integer, intent(in) ::  kpr_l1,  kpr_h1,  kpr_l2,  kpr_h2,  kpr_l3,  kpr_h3
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2, dkdT_l3, dkdT_h3
  real(rt)        ,intent(in)::Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,Snew_l3:Snew_h3,NVAR)
  real(rt)        ,intent(in)::T   (   T_l1:   T_h1,   T_l2:   T_h2,   T_l3:   T_h3)
  real(rt)        ,intent(in)::Ts  (  Ts_l1:  Ts_h1,  Ts_l2:  Ts_h2,  Ts_l3:  Ts_h3)
  real(rt)                   ::kpp ( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,0:ngroups-1)
  real(rt)                   ::kpr ( kpr_l1: kpr_h1, kpr_l2: kpr_h2, kpr_l3: kpr_h3,0:ngroups-1)
  real(rt)                   ::dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,dkdT_l3:dkdT_h3,0:ngroups-1)
  integer, intent(in) :: use_dkdT, validStar, lag_opac

  integer :: i, j, k, g
  real(rt)         :: kp, kr, nu, rho, temp, Ye
  real(rt)         :: kp1, kr1
  real(rt)         :: kp2, kr2
  real(rt)         :: dT
  logical :: comp_kp, comp_kr
  real(rt)        , parameter :: fac = 0.5e0_rt, minfrac = 1.e-8_rt

  if (lag_opac .eq. 1) then
     dkdT(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.e0_rt
     return
  end if

  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     
     rho = Snew(i,j,k,URHO)
     temp = T(i,j,k)
     if (naux > 0) then
        Ye = Snew(i,j,k,UFX)
     else
        Ye = 0.e0_rt
     end if

     if (validStar > 0) then
        dT = fac*abs(Ts(i,j,k) - T(i,j,k))
        dT = max(dT, minfrac*T(i,j,k))
     else
        dT = T(i,j,k) * 1.e-3_rt + 1.e-50_rt
     end if

     do g=0, ngroups-1
        
        nu = nugroup(g)

        comp_kp = .true.
        comp_kr = .true.
        
        call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)
        kpp(i,j,k,g) = kp
        kpr(i,j,k,g) = kr 

        if (use_dkdT .eq. 0) then        
           dkdT(i,j,k,g) = 0.e0_rt
        else

           comp_kp = .true.
           comp_kr = .false.

           call get_opacities(kp1, kr1, rho, temp-dT, Ye, nu, comp_kp, comp_kr)
           call get_opacities(kp2, kr2, rho, temp+dT, Ye, nu, comp_kp, comp_kr)

           dkdT(i,j,k,g) = (kp2-kp1)/(2.e0_rt*dT)

        end if

     end do
  end do
  end do
  end do

end subroutine ca_opacs



subroutine ca_compute_planck( lo, hi,  &
     kpp , kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     stat,stat_l1,stat_l2,stat_l3,stat_h1,stat_h2,stat_h3 ) bind(C, name="ca_compute_planck")

  use rad_params_module, only : ngroups, nugroup
  use opacity_table_module, only : get_opacities
  use network, only : naux
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3
  integer, intent(in) :: stat_l1,stat_l2,stat_l3,stat_h1,stat_h2,stat_h3
  real(rt)                   ::kpp ( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,0:ngroups-1)
  real(rt)        ,intent(in)::stat(stat_l1:stat_h1,stat_l2:stat_h2,stat_l3:stat_h3,NVAR)

  integer :: i, j, k, g
  real(rt)         :: kp, kr, nu, rho, temp, Ye
  logical, parameter :: comp_kp = .true. 
  logical, parameter :: comp_kr = .false.

  do g=0, ngroups-1

     nu = nugroup(g)

     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rho = stat(i,j,k,URHO)
        temp = stat(i,j,k,UTEMP)
        if (naux > 0) then
           Ye = stat(i,j,k,UFX)
        else
           Ye = 0.e0_rt
        end if

        call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

        kpp(i,j,k,g) = kp

     end do
     end do
     end do
  end do

end subroutine ca_compute_planck


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


subroutine ca_compute_powerlaw_kappa_s( lo, hi,  &
     & kappa, k_l1, k_l2, k_l3, k_h1, k_h2, k_h3, &
     &     u, u_l1, u_l2, u_l3, u_h1, u_h2, u_h3, &
     & kappa0, m, n, p, s0, sm, sn, sp, &
     & Tfloor, kfloor) bind(C, name="ca_compute_powerlaw_kappa_s")

  use rad_params_module, only : ngroups, nugroup
  use meth_params_module, only : NVAR, URHO, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: k_l1, k_l2, k_l3, k_h1, k_h2, k_h3, &
                         u_l1, u_l2, u_l3, u_h1, u_h2, u_h3
  real(rt)                     :: kappa(k_l1:k_h1, k_l2:k_h2, k_l3:k_h3, 0:ngroups-1)
  real(rt)        , intent(in) ::     u(u_l1:u_h1, u_l2:u_h2, u_l3:u_h3, NVAR)
  real(rt)        , intent(in) :: kappa0, m, n, p, Tfloor, kfloor
  real(rt)        , intent(in) :: s0, sm, sn, sp

  integer :: i, j, k, g
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, kf, sct, nup, nusp

  if (  m.eq.0.e0_rt .and.  n.eq.0.e0_rt .and.  p.eq.0.e0_rt .and. &
       sm.eq.0.e0_rt .and. sn.eq.0.e0_rt .and. sp.eq.0.e0_rt ) then
     kappa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = kappa0 + s0
  else
     do g = 0, ngroups-1
        nup = nugroup(g)**p
        nusp = nugroup(g)**sp
        do k       = lo(3), hi(3)
           do j    = lo(2), hi(2)
              do i = lo(1), hi(1)
                 Teff = max(u(i,j,k,UTEMP), tiny)
                 Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
                 kf = kappa0 * (u(i,j,k,URHO) ** m) * (Teff ** (-n)) * nup
                 sct = s0 * (u(i,j,k,URHO) ** sm) * (Teff ** (-sn)) * nusp
                 kappa(i,j,k,g) = max(kf+sct, kfloor)
              end do
           end do
        end do
     end do
  end if

end subroutine ca_compute_powerlaw_kappa_s
   

subroutine ca_compute_powerlaw_kappa( lo, hi,  &
     kappa, k_l1, k_l2, k_l3, k_h1, k_h2, k_h3, &
     u    , u_l1, u_l2, u_l3, u_h1, u_h2, u_h3, &
     kappa0, m, n, p, Tfloor, kfloor) bind(C, name="ca_compute_powerlaw_kappa")

  use rad_params_module, only : ngroups, nugroup
  use meth_params_module, only : NVAR, URHO, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: k_l1, k_l2, k_l3, k_h1, k_h2, k_h3, &
                         u_l1, u_l2, u_l3, u_h1, u_h2, u_h3
  real(rt)                     :: kappa(k_l1:k_h1, k_l2:k_h2, k_l3:k_h3, 0:ngroups-1)
  real(rt)        , intent(in) ::     u(u_l1:u_h1, u_l2:u_h2, u_l3:u_h3, NVAR)
  real(rt)        , intent(in) :: kappa0, m, n, p, Tfloor, kfloor

  integer :: i, j, k, g
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, kf, nup

  if (  m.eq.0.e0_rt .and.  n.eq.0.e0_rt .and.  p.eq.0.e0_rt ) then
     kappa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = kappa0
  else
     do g = 0, ngroups-1
        nup = nugroup(g)**p

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           Teff = max(u(i,j,k,UTEMP), tiny)
           Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
           kf = kappa0 * (u(i,j,k,URHO) ** m) * (Teff ** (-n)) * nup
           kappa(i,j,k,g) = max(kf, kfloor)
        end do
        end do
        end do
     end do
  end if

end subroutine ca_compute_powerlaw_kappa


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

