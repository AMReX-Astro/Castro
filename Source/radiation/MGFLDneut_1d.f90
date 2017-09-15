subroutine ca_accel_acoe_neut( lo, hi,  &
     eta1,eta1_l1,eta1_h1,   &
     theT,theT_l1,theT_h1,   &
     theY,theY_l1,theY_h1,   &
     spc , spc_l1, spc_h1,   &
     kap , kap_l1, kap_h1,   &
     aco , aco_l1, aco_h1,   &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: eta1_l1,eta1_h1
  integer, intent(in) :: theT_l1,theT_h1
  integer, intent(in) :: theY_l1,theY_h1
  integer, intent(in) ::  spc_l1, spc_h1
  integer, intent(in) ::  kap_l1, kap_h1
  integer, intent(in) ::  aco_l1, aco_h1
  real(rt)        , intent(in) :: eta1(eta1_l1:eta1_h1)
  real(rt)        , intent(in) :: theT(theT_l1:theT_h1)
  real(rt)        , intent(in) :: theY(theY_l1:theY_h1)
  real(rt)        , intent(in) :: spc ( spc_l1: spc_h1,0:ngroups-1)
  real(rt)        , intent(in) :: kap ( kap_l1: kap_h1,0:ngroups-1)
  real(rt)                     :: aco ( aco_l1: aco_h1)
  real(rt)        , intent(in) :: dt, tau

  integer :: i, g
  real(rt)         :: kbar, kybar, H1, Theta, dt1, foo

  dt1 = (1.e0_rt+tau)/dt

  do i = lo(1), hi(1)
     kbar = 0.e0_rt
     kybar = 0.e0_rt
     do g=0, ngroups-1
        foo = spc(i,g) * kap(i,g)
        kbar = kbar + foo
        kybar = kybar + foo * erg2rhoYe(g)
     end do

     H1 = eta1(i)
     Theta = theY(i) - theT(i)

     aco(i) = (H1*kbar - Theta*kybar) * clight + dt1
  end do

end subroutine ca_accel_acoe_neut


subroutine ca_accel_rhs_neut( lo, hi, &
     Ern , Ern_l1, Ern_h1,  &
     Erl , Erl_l1, Erl_h1,  &
     kap , kap_l1, kap_h1,   &
     etaT,etaT_l1,etaT_h1,   &
     etaY,etaY_l1,etaY_h1,   &
     theT,theT_l1,theT_h1,   &
     theY,theY_l1,theY_h1,   &
     rhs , rhs_l1, rhs_h1,   &
     dt)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Ern_l1, Ern_h1
  integer, intent(in) :: Erl_l1, Erl_h1
  integer, intent(in) :: kap_l1, kap_h1
  integer, intent(in) ::etaT_l1,etaT_h1
  integer, intent(in) ::etaY_l1,etaY_h1
  integer, intent(in) ::theT_l1,theT_h1
  integer, intent(in) ::theY_l1,theY_h1
  integer, intent(in) :: rhs_l1, rhs_h1
  real(rt)        , intent(in) ::Ern ( Ern_l1: Ern_h1,0:ngroups-1)
  real(rt)        , intent(in) ::Erl ( Erl_l1: Erl_h1,0:ngroups-1)
  real(rt)        , intent(in) :: kap( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        , intent(in) ::etaT(etaT_l1:etaT_h1)
  real(rt)        , intent(in) ::etaY(etaY_l1:etaY_h1)
  real(rt)        , intent(in) ::theT(theT_l1:theT_h1)
  real(rt)        , intent(in) ::theY(theY_l1:theY_h1)
  real(rt)                     :: rhs(rhs_l1:rhs_h1)
  real(rt)        , intent(in) :: dt

  integer :: i, g
  real(rt)         :: rtt, ry, H, Theta, foo

  do i = lo(1), hi(1)
     rtt = 0.e0_rt
     ry = 0.e0_rt
     do g=0,ngroups-1
        foo = kap(i,g)*(Ern(i,g)-Erl(i,g))
        rtt = rtt + foo
        ry = ry + foo * erg2rhoYe(g)
     end do

     H = etaT(i) - etaY(i)
     Theta = theY(i) - theT(i)

     rhs(i) = clight*(H*rtt + Theta*ry)
  end do

end subroutine ca_accel_rhs_neut


subroutine ca_accel_spec_neut( lo, hi, &
     Ern , Ern_l1, Ern_h1,  &
     Erl , Erl_l1, Erl_h1,  &
     kap , kap_l1, kap_h1,  &
     etaT,etaT_l1,etaT_h1,  &
     etaY,etaY_l1,etaY_h1,  &
     theT,theT_l1,theT_h1,  &
     theY,theY_l1,theY_h1,  &
     mugT,mugT_l1,mugT_h1,  &
     mugY,mugY_l1,mugY_h1,  &
     spec,spec_l1,spec_h1,  &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in) :: lo(1), hi(1)
  integer,intent(in):: Ern_l1, Ern_h1
  integer,intent(in):: Erl_l1, Erl_h1
  integer,intent(in):: kap_l1, kap_h1
  integer,intent(in)::etaT_l1,etaT_h1
  integer,intent(in)::etaY_l1,etaY_h1
  integer,intent(in)::theT_l1,theT_h1
  integer,intent(in)::theY_l1,theY_h1
  integer,intent(in)::mugT_l1,mugT_h1
  integer,intent(in)::mugY_l1,mugY_h1
  integer,intent(in)::spec_l1,spec_h1
  real(rt)        ,intent(in)::Ern ( Ern_l1: Ern_h1,0:ngroups-1)
  real(rt)        ,intent(in)::Erl ( Erl_l1: Erl_h1,0:ngroups-1)
  real(rt)        ,intent(in)::kap ( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in)::etaY(etaY_l1:etaY_h1)
  real(rt)        ,intent(in)::theT(theT_l1:theT_h1)
  real(rt)        ,intent(in)::theY(theY_l1:theY_h1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in)::mugY(mugY_l1:mugY_h1,0:ngroups-1)
  real(rt)                   ::spec(spec_l1:spec_h1,0:ngroups-1)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i , g
  real(rt)         :: cdt1, rtt, ry, p, q, r, s, foo, sumeps
  real(rt)        ,dimension(0:ngroups-1)::Hg, Tg, epsilon, kapt, kk

  cdt1 = 1.e0_rt/(clight*dt)

  do i = lo(1), hi(1)

     rtt = 0.e0_rt
     ry = 0.e0_rt
     do g=0,ngroups-1
        foo = kap(i,g)*(Ern(i,g)-Erl(i,g))
        rtt = rtt + foo
        ry = ry + foo * erg2rhoYe(g)
     end do

     Hg = mugT(i,:)*etaT(i) - mugY(i,:)*etaY(i)
     Tg = -mugT(i,:)*theT(i) + mugY(i,:)*theY(i)

     kapt = kap(i,:) + (1.e0_rt+tau)*cdt1
     kk = kap(i,:) / kapt

     p = 1.e0_rt - sum(Hg*kk)
     q = 1.e0_rt - sum(Tg*erg2rhoYe*kk)
     r = sum(Hg*erg2rhoYe*kk)
     s = sum(Tg*kk)

     epsilon = ((r*Tg + q*Hg) * rtt + (s*Hg + p*Tg) * ry) / kapt
     
     sumeps = sum(epsilon)
     if (sumeps .eq. 0.e0_rt) then
        sumeps = 1.e-50_rt
     end if

     spec(i,:) = epsilon / sumeps
  end do
  
end subroutine ca_accel_spec_neut


subroutine ca_check_conv_neut( lo, hi, &
     ren, ren_l1, ren_h1,  &
     res, res_l1, res_h1,  &
     re2, re2_l1, re2_h1,  &
     ern, ern_l1, ern_h1,  &
     Tmn, Tmn_l1, Tmn_h1,  &
     Tms, Tms_l1, Tms_h1,  &
     rYn, rYn_l1, rYn_h1,  &
     rYs, rYs_l1, rYs_h1,  &
     rY2, rY2_l1, rY2_h1,  &
     rho, rho_l1, rho_h1,  &
     kap, kap_l1, kap_h1,  &
     jg ,  jg_l1,  jg_h1,  &
     deT, deT_l1, deT_h1,  &
     deY, deY_l1, deY_h1,  &
     rel_re, abs_re, &
     rel_FT, abs_FT, rel_T, abs_T, &
     rel_FY, abs_FY, rel_Y, abs_Y, &
     dt)
  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(1),hi(1)
  integer,intent(in)::ren_l1, ren_h1
  integer,intent(in)::res_l1, res_h1
  integer,intent(in)::re2_l1, re2_h1
  integer,intent(in)::ern_l1, ern_h1
  integer,intent(in)::Tmn_l1, Tmn_h1
  integer,intent(in)::Tms_l1, Tms_h1
  integer,intent(in)::rYn_l1, rYn_h1
  integer,intent(in)::rYs_l1, rYs_h1
  integer,intent(in)::rY2_l1, rY2_h1
  integer,intent(in)::rho_l1, rho_h1
  integer,intent(in)::kap_l1, kap_h1
  integer,intent(in):: jg_l1,  jg_h1
  integer,intent(in)::deT_l1, deT_h1
  integer,intent(in)::deY_l1, deY_h1
  real(rt)        ,intent(in   )::ren(ren_l1:ren_h1)
  real(rt)        ,intent(in   )::res(res_l1:res_h1)
  real(rt)        ,intent(in   )::re2(re2_l1:re2_h1)
  real(rt)        ,intent(in   )::ern(ern_l1:ern_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::Tmn(Tmn_l1:Tmn_h1)
  real(rt)        ,intent(in   )::Tms(Tms_l1:Tms_h1)
  real(rt)        ,intent(in   )::rYn(rYn_l1:rYn_h1)
  real(rt)        ,intent(in   )::rYs(rYs_l1:rYs_h1)
  real(rt)        ,intent(in   )::rY2(rY2_l1:rY2_h1)
  real(rt)        ,intent(in   )::rho(rho_l1:rho_h1)
  real(rt)        ,intent(in   )::kap(kap_l1:kap_h1,0:ngroups-1)
  real(rt)        ,intent(in   ):: jg( jg_l1: jg_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::deT(deT_l1:deT_h1)
  real(rt)        ,intent(in   )::deY(deY_l1:deY_h1)
  real(rt)        ,intent(inout)::rel_re, abs_re
  real(rt)        ,intent(inout)::rel_FT, abs_FT,rel_T, abs_T
  real(rt)        ,intent(inout)::rel_FY, abs_FY,rel_Y, abs_Y
  real(rt)        ,intent(in) :: dt

  integer :: i
  real(rt)         :: chg, relchg, FT, FY, cdt, FTdenom, FYdenom, dTe, dYe

  cdt = clight*dt

  do i=lo(1),hi(1)
     chg = abs(ren(i) - res(i))
     relchg = abs(chg/(ren(i)+1.e-50_rt))
     rel_re = max(rel_re,relchg)
     abs_re = max(abs_re,chg)

     chg = abs(Tmn(i) - Tms(i))
     relchg = abs(chg/(Tmn(i)+1.e-50_rt))
     rel_T = max(rel_T,relchg)
     abs_T = max(abs_T,chg)

     chg = abs(rYn(i) - rYs(i))
     relchg = abs(chg/(rYn(i)+1.e-50_rt))
     rel_Y = max(rel_Y,relchg)
     abs_Y = max(abs_Y,chg/rho(i))

     FT = abs((ren(i)-re2(i)) - cdt*sum(kap(i,:)*Ern(i,:)-jg(i,:)))
     FY = abs((rYn(i)-rY2(i)) - cdt*sum(erg2rhoYe*(kap(i,:)*Ern(i,:)-jg(i,:))))

     dTe = Tmn(i)
     dYe = rYn(i)/rho(i)
     FTdenom = rho(i)*(abs(deT(i)*dTe)+abs(deY(i)*dYe))
     FYdenom = rYn(i)

     rel_FT = max(rel_FT, FT/(FTdenom+1.e-50_rt))
     abs_FT = max(abs_FT, FT)

     rel_FY = max(rel_FY, FY/(FYdenom+1.e-50_rt))
     abs_FY = max(abs_FY, FY)
  end do

end subroutine ca_check_conv_neut


subroutine ca_check_conv_er_neut( lo, hi, &
     Ern, Ern_l1, Ern_h1, &
     Erl, Erl_l1, Erl_h1, &
     kap, kap_l1, kap_h1, &
     etTz, etTz_l1, etTz_h1,  &
     etYz, etYz_l1, etYz_h1,  &
     thTz, thTz_l1, thTz_h1,  &
     thYz, thYz_l1, thYz_h1,  &
     temp, temp_l1, temp_h1,  &
     Ye  ,   Ye_l1,   Ye_h1,  &
     rela, abso, errr, dt)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Ern_l1, Ern_h1
  integer, intent(in) :: Erl_l1, Erl_h1
  integer, intent(in) :: kap_l1, kap_h1
  integer,intent(in)::temp_l1, temp_h1
  integer,intent(in)::  Ye_l1,   Ye_h1
  integer,intent(in)::etTz_l1, etTz_h1
  integer,intent(in)::etYz_l1, etYz_h1
  integer,intent(in)::thTz_l1, thTz_h1
  integer,intent(in)::thYz_l1, thYz_h1
  real(rt)        , intent(in) :: Ern(Ern_l1:Ern_h1,0:ngroups-1)
  real(rt)        , intent(in) :: Erl(Erl_l1:Erl_h1,0:ngroups-1)
  real(rt)        , intent(in) :: kap(kap_l1:kap_h1,0:ngroups-1)
  real(rt)        ,intent(in )::etTz(etTz_l1:etTz_h1)
  real(rt)        ,intent(in )::etYz(etYz_l1:etYz_h1)
  real(rt)        ,intent(in )::thTz(thTz_l1:thTz_h1)
  real(rt)        ,intent(in )::thYz(thYz_l1:thYz_h1)
  real(rt)        ,intent(in )::temp(temp_l1:temp_h1)
  real(rt)        ,intent(in )::Ye  (  Ye_l1:  Ye_h1)
  real(rt)        , intent(inout) :: rela, abso, errr
  real(rt)        , intent(in) :: dt

  integer :: i, g
  real(rt)         :: chg, tot, cdt, der, kdeT, kdeY, err_T, err_Y, err

  cdt = clight * dt
  do i = lo(1), hi(1)
     chg = 0.e0_rt
     tot = 0.e0_rt
     kdeT = 0.e0_rt
     kdeY = 0.e0_rt
     do g=0,ngroups-1
        der = Ern(i,g)-Erl(i,g)
        chg = chg + abs(der)
        tot = tot + abs(Ern(i,g))
        kdeT = kdeT + kap(i,g)*der
        kdeY = kdeY + erg2rhoYe(g)*kap(i,g)*der
     end do
     abso = max(abso, chg)
     rela = max(rela, chg / (tot + 1.e-50_rt))

     err_T =  etTz(i)*kdeT - thTz(i)*kdeY 
     err_Y = -etYz(i)*kdeT + thYz(i)*kdeY 
     err = max(abs(err_T/(temp(i)+1.e-50_rt)), abs(err_Y/(Ye(i)+1.e-50_rt)))
     errr = max(errr, err)
  end do

end subroutine ca_check_conv_er_neut


subroutine ca_compute_coupty( lo, hi, &
     cpt, cpt_l1, cpt_h1, &
     cpy, cpy_l1, cpy_h1, &
     kpp, kpp_l1, kpp_h1, &
     eg ,  eg_l1,  eg_h1, &
     jg ,  jg_l1,  jg_h1)
  
  use rad_params_module, only : ngroups, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: cpt_l1, cpt_h1 
  integer, intent(in) :: cpy_l1, cpy_h1 
  integer, intent(in) :: kpp_l1, kpp_h1 
  integer, intent(in) ::  eg_l1,  eg_h1
  integer, intent(in) ::  jg_l1,  jg_h1
  real(rt)                     :: cpt(cpt_l1:cpt_h1)
  real(rt)                     :: cpy(cpy_l1:cpy_h1)
  real(rt)        , intent(in) :: kpp(kpp_l1:kpp_h1,0:ngroups-1)
  real(rt)        , intent(in) ::  eg( eg_l1: eg_h1,0:ngroups-1)
  real(rt)        , intent(in) ::  jg( jg_l1: jg_h1,0:ngroups-1)

  integer :: i, g
  real(rt)         :: foo

  cpt(lo(1):hi(1)) = 0.e0_rt
  cpy(lo(1):hi(1)) = 0.e0_rt

  do g=0, ngroups-1
     do i=lo(1),hi(1)
        foo = kpp(i,g) * eg(i,g) - jg(i,g)
        cpt(i) = cpt(i) + foo
        cpy(i) = cpy(i) + erg2rhoYe(g) * foo
     end do
  end do

end subroutine ca_compute_coupty


subroutine ca_compute_dedx( lo, hi,  &
     S   ,   S_l1,   S_h1, &
     T   ,   T_l1,   T_h1, &
     Ye  ,  Ye_l1,  Ye_h1, &
     Ts  ,  Ts_l1,  Ts_h1, &
     Yes , Yes_l1, Yes_h1, &
     dedT,dedT_l1,dedT_h1, &
     dedY,dedY_l1,dedY_h1, &
     validStar)

  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos, Tmin=>table_Tmin, Ymin=>table_Yemin
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::    S_l1,    S_h1 
  integer, intent(in) ::    T_l1,    T_h1
  integer, intent(in) ::   Ye_l1,   Ye_h1
  integer, intent(in) ::   Ts_l1,   Ts_h1
  integer, intent(in) ::  Yes_l1,  Yes_h1
  integer, intent(in) :: dedT_l1, dedT_h1
  integer, intent(in) :: dedY_l1, dedY_h1
  real(rt)        , intent(in) :: S   (   S_l1:   S_h1,NVAR)
  real(rt)        , intent(in) :: T   (   T_l1:   T_h1)
  real(rt)        , intent(in) :: Ye  (  Ye_l1:  Ye_h1)
  real(rt)        , intent(in) :: Ts  (  Ts_l1:  Ts_h1)
  real(rt)        , intent(in) :: Yes ( Yes_l1: Yes_h1)
  real(rt)                     :: dedT(dedT_l1:dedT_h1)
  real(rt)                     :: dedY(dedY_l1:dedY_h1)
  integer, intent(in) :: validStar

  integer :: i
  real(rt)         :: rhoinv, e1, e2
  type(eos_t) :: eos_state
  real(rt)         :: dT, dYe, T1, T2, Ye1, Ye2
  real(rt)        , parameter :: fac = 0.5e0_rt, minfrac = 1.e-8_rt

  do i=lo(1),hi(1)

     rhoinv = 1.e0_rt/S(i,URHO)
     eos_state % rho = S(i,URHO)
     eos_state % xn  = S(i,UFS:UFS+nspec-1) * rhoinv
     eos_state % aux = ye(i)

     if (validStar > 0) then
        dT = fac*abs(Ts(i) - T(i))
        dT = max(dT, minfrac*T(i))
        dYe = fac*abs(Yes(i) - Ye(i))
        dYe = max(dYe, minfrac*Ye(i))
     else
        dT = T(i) * 1.e-3_rt + 1.e-50_rt
        dYe = 1.e-4_rt
     end if

     T1 = T(i) - dT
     if (T1 < Tmin*(1.e0_rt+1.e-6_rt)) then
        T1 = T(i)
     end if
     T2 = T(i) + dT

     Ye1 = Ye(i) - dYe
     if (Ye1 < Ymin*(1.e0_rt+1.e-6_rt)) then
        Ye1 = Ye(i)
     end if
     Ye2 = Ye(i) + dYe

     eos_state % T = T1
     call eos(eos_input_rt, eos_state)
     e1 = eos_state % e

     eos_state % T = T2
     call eos(eos_input_rt, eos_state)
     e2 = eos_state % e

     dedT(i) = (e2-e1) / (T2-T1)

     eos_state % T = T(i)

     eos_state % aux = Ye1
     call eos(eos_input_rt, eos_state)
     e1 = eos_state % e

     eos_state % aux = Ye2
     call eos(eos_input_rt, eos_state)
     e2 = eos_state % e
     
     dedY(i) = (e2-e1) / (Ye2-Ye1)

  end do

end subroutine ca_compute_dedx


subroutine ca_compute_eta_the( lo, hi, &
     & etaT, etaT_l1, etaT_h1, &
     & etTz, etTz_l1, etTz_h1, &
     & etaY, etaY_l1, etaY_h1, &
     & etYz, etYz_l1, etYz_h1, &
     & eta1, eta1_l1, eta1_h1, &
     & theT, theT_l1, theT_h1, &
     & thTz, thTz_l1, thTz_h1, &
     & theY, theY_l1, theY_h1, &
     & thYz, thYz_l1, thYz_h1, &
     & the1, the1_l1, the1_h1, &
     & djdT, djdT_l1, djdT_h1, &
     & djdY, djdY_l1, djdY_h1, &
     & dkdT, dkdT_l1, dkdT_h1, &
     & dkdY, dkdY_l1, dkdY_h1, &
     & dedT, dedT_l1, dedT_h1, &
     & dedY, dedY_l1, dedY_h1, &
     & Ers ,  Ers_l1,  Ers_h1, &
     & rho ,  rho_l1,  rho_h1, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: etaT_l1, etaT_h1
  integer, intent(in) :: etTz_l1, etTz_h1
  integer, intent(in) :: etaY_l1, etaY_h1
  integer, intent(in) :: etYz_l1, etYz_h1
  integer, intent(in) :: eta1_l1, eta1_h1
  integer, intent(in) :: theT_l1, theT_h1
  integer, intent(in) :: thTz_l1, thTz_h1
  integer, intent(in) :: theY_l1, theY_h1
  integer, intent(in) :: thYz_l1, thYz_h1
  integer, intent(in) :: the1_l1, the1_h1
  integer, intent(in) :: djdT_l1, djdT_h1
  integer, intent(in) :: djdY_l1, djdY_h1
  integer, intent(in) :: dkdT_l1, dkdT_h1
  integer, intent(in) :: dkdY_l1, dkdY_h1
  integer, intent(in) :: dedT_l1, dedT_h1
  integer, intent(in) :: dedY_l1, dedY_h1
  integer, intent(in) ::  Ers_l1,  Ers_h1
  integer, intent(in) ::  rho_l1,  rho_h1
  real(rt)                     :: etaT(etaT_l1:etaT_h1)
  real(rt)                     :: etTz(etTz_l1:etTz_h1)
  real(rt)                     :: etaY(etaY_l1:etaY_h1)
  real(rt)                     :: etYz(etYz_l1:etYz_h1)
  real(rt)                     :: eta1(eta1_l1:eta1_h1)
  real(rt)                     :: theT(theT_l1:theT_h1)
  real(rt)                     :: thTz(thTz_l1:thTz_h1)
  real(rt)                     :: theY(theY_l1:theY_h1)
  real(rt)                     :: thYz(thYz_l1:thYz_h1)
  real(rt)                     :: the1(the1_l1:the1_h1)
  real(rt)                     :: djdT(djdT_l1:djdT_h1,0:ngroups-1)
  real(rt)                     :: djdY(djdY_l1:djdY_h1,0:ngroups-1)
  real(rt)        , intent(in) :: dkdT(dkdT_l1:dkdT_h1,0:ngroups-1)
  real(rt)        , intent(in) :: dkdY(dkdY_l1:dkdY_h1,0:ngroups-1)
  real(rt)        , intent(in) :: dedT(dedT_l1:dedT_h1)
  real(rt)        , intent(in) :: dedY(dedY_l1:dedY_h1)
  real(rt)        , intent(in) :: Ers ( Ers_l1: Ers_h1,0:ngroups-1)
  real(rt)        , intent(in) :: rho ( rho_l1: rho_h1)
  real(rt)        , intent(in) :: dt, tau

  integer :: i
  real(rt)         :: cdt, det, et, ey, tt, ty, sigma
  real(rt)         :: dZdT(0:ngroups-1), dZdY(0:ngroups-1)
  real(rt)         :: sumdZdT, sumdZdY, fooT, fooY, barT, barY

  sigma = 1.e0_rt + tau
  cdt = clight * dt

  do i = lo(1), hi(1)
     dZdT = djdT(i,:) - dkdT(i,:)*Ers(i,:)
     dZdY = djdY(i,:) - dkdY(i,:)*Ers(i,:)

     sumdZdT = sum(dZdT)
     sumdZdY = sum(dZdY)

     if (sumdZdT .eq. 0.e0_rt) then
        sumdZdT = 1.e-50_rt
     end if

     if (sumdZdY .eq. 0.e0_rt) then
        sumdZdY = 1.e-50_rt
     end if

     fooT = cdt * sumdZdT
     fooY = cdt * sumdZdY
     barT = sigma*rho(i)*dedT(i)
     barY = sigma*rho(i)*dedY(i)

     et = sigma*rho(i) + cdt * sum(erg2rhoYe*dZdY)
     ey = cdt * sum(erg2rhoYe*dZdT)
     tt = barY + fooY
     ty = barT + fooT
     det = ty*et - tt*ey
     etaT(i) = (et * fooT) / det
     etaY(i) = (ey * fooY) / det
     theT(i) = (tt * fooT) / det
     theY(i) = (ty * fooY) / det

     etTz(i) = etaT(i) / sumdZdT
     etYz(i) = etaY(i) / sumdZdY
     thTz(i) = theT(i) / sumdZdT
     thYz(i) = theY(i) / sumdZdY

     eta1(i) = (et*barT - ey*barY) / det
     the1(i) = (sigma*rho(i) * ty) / det

     djdT(i,:) = dZdT / sumdZdT
     djdY(i,:) = dZdY / sumdZdY
  end do

end subroutine ca_compute_eta_the


subroutine ca_compute_rhs_neut( lo, hi, &
     rhs , rhs_l1, rhs_h1, &
     jg  ,  jg_l1,  jg_h1, &
     mugT,mugT_l1,mugT_h1, &
     mugY,mugY_l1,mugY_h1, &
     cpT , cpT_l1, cpT_h1, &
     cpY , cpY_l1, cpY_h1, &
     etaT,etaT_l1,etaT_h1, &
     etay,etaY_l1,etaY_h1, &
     theT,theT_l1,theT_h1, &
     they,theY_l1,theY_h1, &
     Er2 , Er2_l1, Er2_h1, &
     re2 , re2_l1, re2_h1, &
     rY2 , rY2_l1, rY2_h1, &
     Ers , Ers_l1, Ers_h1, &
     res , res_l1, res_h1, &
     rYs , rYs_l1, rYs_h1, &
     r, dt, igroup, tau) bind(C, name="ca_compute_rhs_neut")

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: rhs_l1, rhs_h1
  integer,intent(in)::  jg_l1,  jg_h1
  integer,intent(in)::mugT_l1,mugT_h1
  integer,intent(in)::mugY_l1,mugY_h1
  integer,intent(in):: cpT_l1, cpT_h1
  integer,intent(in):: cpY_l1, cpY_h1
  integer,intent(in)::etaT_l1,etaT_h1
  integer,intent(in)::etaY_l1,etaY_h1
  integer,intent(in)::theT_l1,theT_h1
  integer,intent(in)::theY_l1,theY_h1
  integer,intent(in):: Er2_l1, Er2_h1
  integer,intent(in):: re2_l1, re2_h1
  integer,intent(in):: rY2_l1, rY2_h1
  integer,intent(in):: Ers_l1, Ers_h1
  integer,intent(in):: res_l1, res_h1
  integer,intent(in):: rYs_l1, rYs_h1
  real(rt)                   ::rhs ( rhs_l1: rhs_h1)
  real(rt)        ,intent(in)::jg  (  jg_l1:  jg_h1,0:ngroups-1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in)::mugY(mugY_l1:mugY_h1,0:ngroups-1)
  real(rt)        ,intent(in)::cpT ( cpT_l1: cpT_h1)
  real(rt)        ,intent(in)::cpY ( cpY_l1: cpY_h1)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in)::etaY(etaY_l1:etaY_h1)
  real(rt)        ,intent(in)::theT(theT_l1:theT_h1)
  real(rt)        ,intent(in)::theY(theY_l1:theY_h1)
  real(rt)        ,intent(in)::Er2 ( Er2_l1: Er2_h1,0:ngroups-1)
  real(rt)        ,intent(in)::re2 ( re2_l1: re2_h1)
  real(rt)        ,intent(in)::rY2 ( rY2_l1: rY2_h1)
  real(rt)        ,intent(in)::Ers ( Ers_l1: Ers_h1,0:ngroups-1)
  real(rt)        ,intent(in)::res ( res_l1: res_h1)
  real(rt)        ,intent(in)::rYs ( rYs_l1: rYs_h1)
  real(rt)        ,intent(in) ::   r(lo(1):hi(1))
  real(rt)        ,intent(in) :: dt, tau
  integer, intent(in) :: igroup

  integer :: i
  real(rt)         :: Hg, thetag, dt1

  dt1 = 1.e0_rt/dt
  do i=lo(1),hi(1)
     Hg = mugT(i,igroup)*etaT(i) - mugY(i,igroup)*etaY(i)
     thetag = mugY(i,igroup)*theY(i) - mugT(i,igroup)*theT(i)

     rhs(i) = clight*(jg(i,igroup) + Hg*cpT(i) + thetag*cpY(i)) &
          + dt1 * (Er2(i,igroup) - Hg*(res(i)-re2(i))  &
          &                      - thetag*(rYs(i)-rY2(i)) &
          &        + tau*Ers(i,igroup))

     rhs(i) = r(i) * rhs(i)
   end do

end subroutine ca_compute_rhs_neut


subroutine ca_local_accel_neut( lo, hi,  &
     Ern , Ern_l1, Ern_h1,  &
     Erl , Erl_l1, Erl_h1,  &
     kap , kap_l1, kap_h1,  &
     etaT,etaT_l1,etaT_h1,  &
     etaY,etaY_l1,etaY_h1,  &
     theT,theT_l1,theT_h1,  &
     theY,theY_l1,theY_h1,  &
     mugT,mugT_l1,mugT_h1,  &
     mugY,mugY_l1,mugY_h1,  &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: Ern_l1, Ern_h1
  integer,intent(in):: Erl_l1, Erl_h1
  integer,intent(in):: kap_l1, kap_h1
  integer,intent(in)::etaT_l1,etaT_h1
  integer,intent(in)::etaY_l1,etaY_h1
  integer,intent(in)::theT_l1,theT_h1
  integer,intent(in)::theY_l1,theY_h1
  integer,intent(in)::mugT_l1,mugT_h1
  integer,intent(in)::mugY_l1,mugY_h1
  real(rt)                   ::Ern ( Ern_l1: Ern_h1,0:ngroups-1)
  real(rt)        ,intent(in)::Erl ( Erl_l1: Erl_h1,0:ngroups-1)
  real(rt)        ,intent(in)::kap ( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in)::etaY(etaY_l1:etaY_h1)
  real(rt)        ,intent(in)::theT(theT_l1:theT_h1)
  real(rt)        ,intent(in)::theY(theY_l1:theY_h1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in)::mugY(mugY_l1:mugY_h1,0:ngroups-1)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i 
  real(rt)         :: cdt1, rtt, ry, p, q, r, s 
  real(rt)        ,dimension(0:ngroups-1)::Hg, Tg, epsilon, kapt, kk

  cdt1 = 1.e0_rt/(clight*dt)

  do i = lo(1), hi(1)
     rtt = sum(kap(i,:)*(Ern(i,:)-Erl(i,:)))
     ry = sum(erg2rhoYe(:)*kap(i,:)*(Ern(i,:)-Erl(i,:)))

     Hg = mugT(i,:)*etaT(i) - mugY(i,:)*etaY(i)
     Tg = -mugT(i,:)*theT(i) + mugY(i,:)*theY(i)

     kapt = kap(i,:) + (1.e0_rt+tau)*cdt1
     kk = kap(i,:) / kapt

     p = 1.e0_rt-sum(Hg*kk)
     q = 1.e0_rt-sum(Tg*erg2rhoYe*kk)
     r = sum(Hg*erg2rhoYe*kk)
     s = sum(Tg*kk)

     epsilon = ((r*Tg + q*Hg) * rtt + (s*Hg + p*Tg) * ry)  &
          / (kapt*(p*q-r*s) + 1.e-50_rt)

     Ern(i,:) = Ern(i,:) + epsilon
  end do
  
end subroutine ca_local_accel_neut


subroutine ca_opac_emis_neut( lo, hi,  &
     Snew,Snew_l1,Snew_h1, &
     T   ,   T_l1,   T_h1, &
     Ye  ,  Ye_l1,  Ye_h1, &
     Ts  ,  Ts_l1,  Ts_h1, &
     Yes , Yes_l1, Yes_h1, &
     kpp , kpp_l1, kpp_h1, &
     kpr , kpr_l1, kpr_h1, &
     jg  ,   j_l1,   j_h1, &
     djdT,djdT_l1,djdT_h1, &
     djdY,djdY_l1,djdY_h1, &
     dkdT,dkdT_l1,dkdT_h1, &
     dkdY,dkdY_l1,dkdY_h1, &
     use_dkdT, validStar, lag_opac) 

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Snew_l1, Snew_h1 
  integer, intent(in) ::    T_l1,    T_h1
  integer, intent(in) ::   Ye_l1,   Ye_h1
  integer, intent(in) ::   Ts_l1,   Ts_h1
  integer, intent(in) ::  Yes_l1,  Yes_h1
  integer, intent(in) ::  kpp_l1,  kpp_h1 
  integer, intent(in) ::  kpr_l1,  kpr_h1
  integer, intent(in) ::    j_l1,    j_h1
  integer, intent(in) :: djdT_l1, djdT_h1
  integer, intent(in) :: djdY_l1, djdY_h1
  integer, intent(in) :: dkdT_l1, dkdT_h1
  integer, intent(in) :: dkdY_l1, dkdY_h1
  real(rt)        , intent(in) :: Snew(Snew_l1:Snew_h1,NVAR)
  real(rt)        , intent(in) :: T   (   T_l1:   T_h1)
  real(rt)        , intent(in) :: Ye  (  Ye_l1:  Ye_h1)
  real(rt)        , intent(in) :: Ts  (  Ts_l1:  Ts_h1)
  real(rt)        , intent(in) :: Yes ( Yes_l1: Yes_h1)
  real(rt)                     :: kpp ( kpp_l1: kpp_h1,0:ngroups-1)
  real(rt)                     :: kpr ( kpr_l1: kpr_h1,0:ngroups-1)
  real(rt)                     :: jg  (   j_l1:   j_h1,0:ngroups-1)
  real(rt)                     :: djdT(djdT_l1:djdT_h1,0:ngroups-1)
  real(rt)                     :: djdY(djdY_l1:djdY_h1,0:ngroups-1)
  real(rt)                     :: dkdT(dkdT_l1:dkdT_h1,0:ngroups-1)
  real(rt)                     :: dkdY(dkdY_l1:dkdY_h1,0:ngroups-1)
  integer, intent(in) :: use_dkdT, validStar, lag_opac

  integer :: i, g, inu
  real(rt)         :: ab, sc, delta, eta, er, der, rho, temp
  real(rt)         :: ab1, sc1, delta1, eta1
  real(rt)         :: ab2, sc2, delta2, eta2 
  real(rt)         :: dT, dYe
  real(rt)         :: Bg, Bg1, Bg2
  logical :: comp_ab, comp_sc, comp_eta
  real(rt)        , parameter :: fac = 0.5e0_rt, minfrac = 1.e-8_rt

  do i=lo(1), hi(1)
     
     rho = Snew(i,URHO)
     temp = T(i)

     if (validStar > 0) then
        dT = fac*abs(Ts(i) - T(i))
        dT = max(dT, minfrac*T(i))
        dYe = fac*abs(Yes(i) - Ye(i))
        dYe = max(dYe, minfrac*Ye(i))
     else
        dT = T(i) * 1.e-3_rt + 1.e-50_rt
        dYe = 1.e-4_rt
     end if

     do g=0, ngroups-1
        
        call prep_opacity(g, inu, er, der)

        if (lag_opac .eq. 1) then

           dkdT(i,g) = 0.e0_rt
           dkdY(i,g) = 0.e0_rt              

           comp_ab = .true.
           comp_sc = .false.
           comp_eta = .true.

           call get_opacity_emissivity(ab, sc, delta, eta, &
                rho, Ye(i), temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg = eta * der / ab
           jg(i,g) = Bg * kpp(i,g)

           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i)-dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg1 = eta1 * der / ab1

           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i)+dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg2 = eta2 * der / ab2

           djdY(i,g) = (Bg2-Bg1)*kpp(i,g)/(2.e0_rt*dye)

           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i), temp-dT, er, inu, comp_ab, comp_sc, comp_eta)
           Bg1 = eta1 * der / ab1
           
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i), temp+dT, er, inu, comp_ab, comp_sc, comp_eta)
           Bg2 = eta2 * der / ab2

           djdT(i,g) = (Bg2-Bg1)*kpp(i,g)/(2.e0_rt*dT)
           
        else

           comp_ab = .true.
           comp_sc = .true.
           comp_eta = .true.
           
           call get_opacity_emissivity(ab, sc, delta, eta, &
                rho, Ye(i), temp, er, inu, comp_ab, comp_sc, comp_eta)
           kpp(i,g) = ab
           kpr(i,g) = ab + sc * (1.e0_rt - delta/3.e0_rt)
           jg(i,g) = eta * der 
           
           comp_ab = .true.
           comp_sc = .false.
           comp_eta = .true.
           
           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i)-dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i)+dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           djdY(i,g) = (eta2-eta1)*der/(2.e0_rt*dye)
           dkdY(i,g) = (ab2-ab1)/(2.e0_rt*dye)
           if (use_dkdT .eq. 0) then
              djdY(i,g) = (eta2/ab2 - eta1/ab1)*der/(2.e0_rt*dye) * ab
              dkdY(i,g) = 0.e0_rt
           end if
           
           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i), temp-dT, er, inu, comp_ab, comp_sc, comp_eta)
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i), temp+dT, er, inu, comp_ab, comp_sc, comp_eta)
           djdT(i,g) = (eta2-eta1)*der/(2.e0_rt*dT)
           dkdT(i,g) = (ab2-ab1)/(2.e0_rt*dT)
           if (use_dkdT .eq. 0) then        
              djdT(i,g) = (eta2/ab2 - eta1/ab1)*der/(2.e0_rt*dT) * ab
              dkdT(i,g) = 0.e0_rt
           end if

        end if

     end do
  end do

end subroutine ca_opac_emis_neut


subroutine ca_state_update_neut( lo, hi, &
     state, state_l1, state_h1,  &
     rhoe,   rhoe_l1,  rhoe_h1,  &
     Ye  ,     Ye_l1,    Ye_h1,  &
     temp,   temp_l1,  temp_h1,  &
     msk ,    msk_l1,   msk_h1,  &
     derat, dTrat, dye)

  use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1) 
  integer, intent(in) :: state_l1, state_h1
  integer, intent(in) ::  rhoe_l1,  rhoe_h1
  integer, intent(in) ::    Ye_l1,    Ye_h1
  integer, intent(in) ::  temp_l1,  temp_h1
  integer, intent(in) ::   msk_l1,   msk_h1
  real(rt)        , intent(in) :: rhoe( rhoe_l1: rhoe_h1)
  real(rt)        , intent(in) ::   Ye(   Ye_l1:   Ye_h1)
  real(rt)        , intent(in) :: temp( temp_l1: temp_h1)
  real(rt)        , intent(in) ::  msk(  msk_l1:  msk_h1)
  real(rt)                     ::state(state_l1:state_h1, NVAR)
  real(rt)        , intent(inout) :: derat, dTrat, dye

  integer :: i
  real(rt)         :: ei, ek, Told, Yeold

  do i=lo(1), hi(1)
     ei = state(i,UEINT)
     derat = max(derat, abs((rhoe(i) - ei)*msk(i)/ (ei + 1.e-50_rt)))
     ek = state(i,UEDEN) - state(i,UEINT)
     state(i,UEINT) = rhoe(i)
     state(i,UEDEN) = rhoe(i) + ek

     Told = state(i,UTEMP);
     dTrat = max(dTrat, abs((temp(i)-Told)*msk(i)/ (Told + 1.e-50_rt)))
     state(i,UTEMP) = temp(i)

     Yeold = state(i,UFX) / state(i,URHO)
     dye = max(dye, abs((Ye(i)-Yeold)*msk(i)))
     state(i,UFX) = state(i,URHO) * Ye(i)
  end do

end subroutine ca_state_update_neut


subroutine ca_update_matter_neut( lo, hi,  &
     re_n, re_n_l1, re_n_h1,  &
     rY_n, rY_n_l1, rY_n_h1,  &
     Ye_n, Ye_n_l1, Ye_n_h1,  &
     Er_n, Er_n_l1, Er_n_h1,  &
     Er_l, Er_l_l1, Er_l_h1,  &
     re_s, re_s_l1, re_s_h1,  &
     rY_s, rY_s_l1, rY_s_h1,  &
     re_2, re_2_l1, re_2_h1,  &
     rY_2, rY_2_l1, rY_2_h1,  &
     etaT, etaT_l1, etaT_h1,  &
     etaY, etaY_l1, etaY_h1,  &
     eta1, eta1_l1, eta1_h1,  &
     theT, theT_l1, theT_h1,  &
     theY, theY_l1, theY_h1,  &
     the1, the1_l1, the1_h1,  &
      cpt,  cpt_l1,  cpt_h1,  &
      cpy,  cpy_l1,  cpy_h1,  &
      kpp,  kpp_l1,  kpp_h1,  &
     mugT, mugT_l1, mugT_h1,  &
     mugY, mugY_l1, mugY_h1,  &
     Snew, Snew_l1, Snew_h1,  &
     dt, tau)

  use rad_params_module, only : ngroups, erg2rhoYe, clight
  use meth_params_module, only : NVAR, URHO

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(1),hi(1)
  integer,intent(in)::re_n_l1, re_n_h1
  integer,intent(in)::rY_n_l1, rY_n_h1
  integer,intent(in)::Ye_n_l1, Ye_n_h1
  integer,intent(in)::Er_n_l1, Er_n_h1
  integer,intent(in)::Er_l_l1, Er_l_h1
  integer,intent(in)::re_s_l1, re_s_h1
  integer,intent(in)::rY_s_l1, rY_s_h1
  integer,intent(in)::re_2_l1, re_2_h1
  integer,intent(in)::rY_2_l1, rY_2_h1
  integer,intent(in)::etaT_l1, etaT_h1
  integer,intent(in)::etaY_l1, etaY_h1
  integer,intent(in)::eta1_l1, eta1_h1
  integer,intent(in)::theT_l1, theT_h1
  integer,intent(in)::theY_l1, theY_h1
  integer,intent(in)::the1_l1, the1_h1
  integer,intent(in):: cpt_l1,  cpt_h1
  integer,intent(in):: cpy_l1,  cpy_h1
  integer,intent(in):: kpp_l1,  kpp_h1
  integer,intent(in)::mugT_l1, mugT_h1
  integer,intent(in)::mugY_l1, mugY_h1
  integer,intent(in)::Snew_l1, Snew_h1
  real(rt)                   ::re_n(re_n_l1:re_n_h1)
  real(rt)                   ::rY_n(rY_n_l1:rY_n_h1)
  real(rt)                   ::Ye_n(Ye_n_l1:Ye_n_h1)
  real(rt)        ,intent(in)::Er_n(Er_n_l1:Er_n_h1,0:ngroups-1)
  real(rt)        ,intent(in)::Er_l(Er_l_l1:Er_l_h1,0:ngroups-1)
  real(rt)        ,intent(in)::re_s(re_s_l1:re_s_h1)
  real(rt)        ,intent(in)::rY_s(rY_s_l1:rY_s_h1)
  real(rt)        ,intent(in)::re_2(re_2_l1:re_2_h1)
  real(rt)        ,intent(in)::rY_2(rY_2_l1:rY_2_h1)
  real(rt)        ,intent(in)::etaT(etaT_l1:etaT_h1)
  real(rt)        ,intent(in)::etaY(etaY_l1:etaY_h1)
  real(rt)        ,intent(in)::eta1(eta1_l1:eta1_h1)
  real(rt)        ,intent(in)::theT(theT_l1:theT_h1)
  real(rt)        ,intent(in)::theY(theY_l1:theY_h1)
  real(rt)        ,intent(in)::the1(the1_l1:the1_h1)
  real(rt)        ,intent(in):: cpt( cpt_l1: cpt_h1)
  real(rt)        ,intent(in):: cpy( cpy_l1: cpy_h1)
  real(rt)        ,intent(in):: kpp( kpp_l1: kpp_h1,0:ngroups-1)
  real(rt)        ,intent(in)::mugT(mugT_l1:mugT_h1,0:ngroups-1)
  real(rt)        ,intent(in)::mugY(mugY_l1:mugY_h1,0:ngroups-1)
  real(rt)        ,intent(in)::Snew(Snew_l1:Snew_h1,NVAR)
  real(rt)        ,intent(in) :: dt, tau

  integer :: i,g
  real(rt)         :: cdt, H1, Theta, Thbar1, Hbar 
  real(rt)         :: dkEE, dkEEY, foo, chg, chgY

  cdt = clight * dt
  do i = lo(1), hi(1)

     H1 = eta1(i)
     Thbar1 = the1(i)

     Theta = theY(i) - theT(i)
     Hbar = sum(erg2rhoYe*(mugT(i,:)*etaT(i) - mugY(i,:)*etaY(i)))

     dkEE = 0.e0_rt
     dkEEY = 0.e0_rt
     do g=0, ngroups-1
        foo = kpp(i,g)*(Er_n(i,g)-Er_l(i,g))
        dkEE = dkEE + foo
        dkEEY = dkEEY + foo * erg2rhoYe(g)
     end do

     chg = cdt*dkEE + H1*((re_2(i)-re_s(i)) + cdt*cpt(i)) &
          - Theta*((rY_2(i)-rY_s(i)) + cdt*cpy(i))

     chgY = cdt*dkEEY + Thbar1*((rY_2(i)-rY_s(i)) + cdt*cpy(i)) &
          - Hbar*((re_2(i)-re_s(i)) + cdt*cpt(i))

     re_n(i) = re_s(i) + chg
     rY_n(i) = rY_s(i) + chgY

     re_n(i) = (re_n(i) + tau*re_s(i)) / (1.e0_rt+tau)
     rY_n(i) = (rY_n(i) + tau*rY_s(i)) / (1.e0_rt+tau)

     Ye_n(i) = rY_n(i) / Snew(i,URHO)

     ! temperature will be updated after exiting this subroutine
  end do

end subroutine ca_update_matter_neut


subroutine ca_ncupdate_matter_neut( lo, hi,  &
     Tp_n, Tp_n_l1, Tp_n_h1,  &
     Ye_n, Ye_n_l1, Ye_n_h1,  &
     Er_n, Er_n_l1, Er_n_h1,  &
     re_s, re_s_l1, re_s_h1,  &
     rY_s, rY_s_l1, rY_s_h1,  &
     re_2, re_2_l1, re_2_h1,  &
     rY_2, rY_2_l1, rY_2_h1,  &
     etTz, etTz_l1, etTz_h1,  &
     etYz, etYz_l1, etYz_h1,  &
     thTz, thTz_l1, thTz_h1,  &
     thYz, thYz_l1, thYz_h1,  &
      kpp,  kpp_l1,  kpp_h1,  &
       jg,   jg_l1,   jg_h1,  &
     dt)

  use rad_params_module, only : ngroups, erg2rhoYe, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(1),hi(1)
  integer,intent(in)::Tp_n_l1, Tp_n_h1
  integer,intent(in)::Ye_n_l1, Ye_n_h1
  integer,intent(in)::Er_n_l1, Er_n_h1
  integer,intent(in)::re_s_l1, re_s_h1
  integer,intent(in)::rY_s_l1, rY_s_h1
  integer,intent(in)::re_2_l1, re_2_h1
  integer,intent(in)::rY_2_l1, rY_2_h1
  integer,intent(in)::etTz_l1, etTz_h1
  integer,intent(in)::etYz_l1, etYz_h1
  integer,intent(in)::thTz_l1, thTz_h1
  integer,intent(in)::thYz_l1, thYz_h1
  integer,intent(in):: kpp_l1,  kpp_h1
  integer,intent(in)::  jg_l1,   jg_h1
  real(rt)                   ::Tp_n(Tp_n_l1:Tp_n_h1)
  real(rt)                   ::Ye_n(Ye_n_l1:Ye_n_h1)
  real(rt)        ,intent(in)::Er_n(Er_n_l1:Er_n_h1,0:ngroups-1)
  real(rt)        ,intent(in)::re_s(re_s_l1:re_s_h1)
  real(rt)        ,intent(in)::rY_s(rY_s_l1:rY_s_h1)
  real(rt)        ,intent(in)::re_2(re_2_l1:re_2_h1)
  real(rt)        ,intent(in)::rY_2(rY_2_l1:rY_2_h1)
  real(rt)        ,intent(in)::etTz(etTz_l1:etTz_h1)
  real(rt)        ,intent(in)::etYz(etYz_l1:etYz_h1)
  real(rt)        ,intent(in)::thTz(thTz_l1:thTz_h1)
  real(rt)        ,intent(in)::thYz(thYz_l1:thYz_h1)
  real(rt)        ,intent(in):: kpp( kpp_l1: kpp_h1,0:ngroups-1)
  real(rt)        ,intent(in)::  jg(  jg_l1:  jg_h1,0:ngroups-1)
  real(rt)        ,intent(in) :: dt

   integer :: i,g
   real(rt)         :: cdt1, cpT, cpY, foo, scrch_re, scrch_rY
   real(rt)         :: dTemp, dYe
   real(rt)        , parameter :: fac = 0.01e0_rt

   cdt1 = 1.e0_rt / (clight * dt)
   do i = lo(1), hi(1)

      cpT = 0.e0_rt
      cpY = 0.e0_rt
      do g = 0, ngroups-1
         foo = kpp(i,g)*Er_n(i,g) - jg(i,g)
         cpT = cpT + foo
         cpY = cpY + erg2rhoYe(g)*foo
      end do

      scrch_re = cpT - (re_s(i) - re_2(i)) * cdt1
      scrch_rY = cpY - (rY_s(i) - rY_2(i)) * cdt1

      dTemp = etTz(i)*scrch_re - thTz(i)*scrch_rY
      dYe   = -etYz(i)*scrch_re + thYz(i)*scrch_rY

      if (abs(dTemp/(Tp_n(i)+1.e-50_rt)) > fac) then
         dTemp = sign(fac*Tp_n(i), dTemp)
      end if

      if (abs(dYe/(Ye_n(i)+1.e-50_rt)) > fac) then
         dYe = sign(fac*Ye_n(i), dYe)
      end if

     Tp_n(i) = Tp_n(i) + dTemp
     Ye_n(i) = Ye_n(i) + dYe

  end do

end subroutine ca_ncupdate_matter_neut

subroutine ca_compute_rosseland_neut( lo, hi,  &
     kpr , kpr_l1, kpr_h1, &
     stat,stat_l1,stat_h1 )

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  kpr_l1,  kpr_h1
  integer, intent(in) :: stat_l1, stat_h1 
  real(rt)                     :: kpr ( kpr_l1: kpr_h1,0:ngroups-1)
  real(rt)        , intent(in) :: stat(stat_l1:stat_h1,NVAR)

  integer :: i, g, inu
  real(rt)         :: ab, sc, delta, eta, er, der, rho, temp, Ye
  logical :: comp_ab, comp_sc, comp_eta

  comp_ab = .true.
  comp_sc = .true.
  comp_eta = .false.

  do g=0, ngroups-1

     call prep_opacity(g, inu, er, der)

     do i = lo(1), hi(1)

        rho = stat(i,URHO)
        temp = stat(i,UTEMP)
        Ye = stat(i,UFX) / rho

        call get_opacity_emissivity(ab, sc, delta, eta, &
             rho, Ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

        kpr(i,g) = ab + sc * (1.e0_rt - delta/3.e0_rt)

     end do
  end do

end subroutine ca_compute_rosseland_neut


subroutine ca_compute_planck_neut( lo, hi,  &
     kpp , kpp_l1, kpp_h1, &
     stat,stat_l1,stat_h1 )

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  kpp_l1,  kpp_h1
  integer, intent(in) :: stat_l1, stat_h1 
  real(rt)                     :: kpp ( kpp_l1: kpp_h1,0:ngroups-1)
  real(rt)        , intent(in) :: stat(stat_l1:stat_h1,NVAR)

  integer :: i, g, inu
  real(rt)         :: ab, sc, delta, eta, er, der, rho, temp, Ye
  logical :: comp_ab, comp_sc, comp_eta

  comp_ab = .true.
  comp_sc = .false.
  comp_eta = .false.

  do g=0, ngroups-1

     call prep_opacity(g, inu, er, der)

     do i = lo(1), hi(1)

        rho = stat(i,URHO)
        temp = stat(i,UTEMP)
        Ye = stat(i,UFX) / rho

        call get_opacity_emissivity(ab, sc, delta, eta, &
             rho, Ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

        kpp(i,g) = ab 

     end do
  end do

end subroutine ca_compute_planck_neut

