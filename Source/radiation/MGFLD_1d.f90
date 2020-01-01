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


subroutine ca_check_conv( lo, hi, &
     ren, ren_l1, ren_h1,  &
     res, res_l1, res_h1,  &
     re2, re2_l1, re2_h1,  &
     ern, ern_l1, ern_h1,  &
     Tmn, Tmn_l1, Tmn_h1,  &
     Tms, Tms_l1, Tms_h1,  &
     rho, rho_l1, rho_h1,  &
     kap, kap_l1, kap_h1,  &
     jg ,  jg_l1,  jg_h1,  &
     deT, deT_l1, deT_h1,  &
     rel_re, abs_re, & 
     rel_FT, abs_FT, rel_T, abs_T, &
     dt)
  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in)::lo(1),hi(1)
  integer,intent(in)::ren_l1, ren_h1
  integer,intent(in)::res_l1, res_h1
  integer,intent(in)::re2_l1, re2_h1
  integer,intent(in)::ern_l1, ern_h1
  integer,intent(in)::Tmn_l1, Tmn_h1
  integer,intent(in)::Tms_l1, Tms_h1
  integer,intent(in)::rho_l1, rho_h1
  integer,intent(in)::kap_l1, kap_h1
  integer,intent(in):: jg_l1,  jg_h1
  integer,intent(in)::deT_l1, deT_h1
  real(rt)        ,intent(in   )::ren(ren_l1:ren_h1)
  real(rt)        ,intent(in   )::res(res_l1:res_h1)
  real(rt)        ,intent(in   )::re2(re2_l1:re2_h1)
  real(rt)        ,intent(in   )::ern(ern_l1:ern_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::Tmn(Tmn_l1:Tmn_h1)
  real(rt)        ,intent(in   )::Tms(Tms_l1:Tms_h1)
  real(rt)        ,intent(in   )::rho(rho_l1:rho_h1)
  real(rt)        ,intent(in   )::kap(kap_l1:kap_h1,0:ngroups-1)
  real(rt)        ,intent(in   ):: jg( jg_l1: jg_h1,0:ngroups-1)
  real(rt)        ,intent(in   )::deT(deT_l1:deT_h1)
  real(rt)        ,intent(inout)::rel_re, abs_re 
  real(rt)        ,intent(inout)::rel_FT, abs_FT,rel_T, abs_T
  real(rt)        ,intent(in) :: dt

  integer :: i
  real(rt)         :: chg, relchg, FT, cdt, FTdenom, dTe

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

     FT = abs((ren(i)-re2(i)) - cdt*sum(kap(i,:)*Ern(i,:)-jg(i,:)))

     dTe = Tmn(i)
     FTdenom = rho(i)*abs(deT(i)*dTe)
!     FTdenom = max(abs(ren(i)-re2(i)), abs(ren(i)*1.e-15_rt))

     rel_FT = max(rel_FT, FT/(FTdenom+1.e-50_rt))
     abs_FT = max(abs_FT, FT)
  end do

end subroutine ca_check_conv


subroutine ca_check_conv_er( lo, hi, &
     Ern, Ern_l1, Ern_h1, &
     Erl, Erl_l1, Erl_h1, &
     kap, kap_l1, kap_h1, &
     etTz, etTz_l1, etTz_h1,  &
     temp, temp_l1, temp_h1,  &
     rela, abso, errr, dt)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Ern_l1, Ern_h1
  integer, intent(in) :: Erl_l1, Erl_h1
  integer, intent(in) :: kap_l1, kap_h1
  integer,intent(in)::temp_l1, temp_h1
  integer,intent(in)::etTz_l1, etTz_h1
  real(rt)        , intent(in) :: Ern(Ern_l1:Ern_h1,0:ngroups-1)
  real(rt)        , intent(in) :: Erl(Erl_l1:Erl_h1,0:ngroups-1)
  real(rt)        , intent(in) :: kap(kap_l1:kap_h1,0:ngroups-1)
  real(rt)        ,intent(in )::etTz(etTz_l1:etTz_h1)
  real(rt)        ,intent(in )::temp(temp_l1:temp_h1)
  real(rt)        , intent(inout) :: rela, abso, errr
  real(rt)        , intent(in) :: dt

  integer :: i, g
  real(rt)         :: chg, tot, cdt, der, kde, err_T, err

  cdt = clight * dt
  do i = lo(1), hi(1)
     chg = 0.e0_rt
     tot = 0.e0_rt
     kde = 0.e0_rt
     do g=0,ngroups-1
        der = Ern(i,g)-Erl(i,g)
        chg = chg + abs(der)
        tot = tot + abs(Ern(i,g))
        kde = kde + kap(i,g)*der
     end do
     abso = max(abso, chg)
     rela = max(rela, chg / (tot + 1.e-50_rt))

     err_T =  etTz(i)*kde
     err = abs(err_T/(temp(i)+1.e-50_rt))
     errr = max(errr, err)
  end do

end subroutine ca_check_conv_er


subroutine ca_compute_coupt( lo, hi, &
     cpt, cpt_l1, cpt_h1, &
     kpp, kpp_l1, kpp_h1, &
     eg ,  eg_l1,  eg_h1, &
     jg ,  jg_l1,  jg_h1)
  
  use rad_params_module, only : ngroups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: cpt_l1, cpt_h1 
  integer, intent(in) :: kpp_l1, kpp_h1 
  integer, intent(in) ::  eg_l1,  eg_h1
  integer, intent(in) ::  jg_l1,  jg_h1
  real(rt)                     :: cpt(cpt_l1:cpt_h1)
  real(rt)        , intent(in) :: kpp(kpp_l1:kpp_h1,0:ngroups-1)
  real(rt)        , intent(in) ::  eg( eg_l1: eg_h1,0:ngroups-1)
  real(rt)        , intent(in) ::  jg( jg_l1: jg_h1,0:ngroups-1)

  integer :: i, g

  cpt(lo(1):hi(1)) = 0.e0_rt

  do g=0, ngroups-1
     do i=lo(1),hi(1)
        cpt(i) = cpt(i) + (kpp(i,g) * eg(i,g) - jg(i,g))
     end do
  end do

end subroutine ca_compute_coupt


subroutine ca_compute_etat( lo, hi, &
     & etaT, etaT_l1, etaT_h1, &
     & etTz, etTz_l1, etTz_h1, &
     & eta1, eta1_l1, eta1_h1, &
     & djdT, djdT_l1, djdT_h1, &
     & dkdT, dkdT_l1, dkdT_h1, &
     & dedT, dedT_l1, dedT_h1, &
     & Ers ,  Ers_l1,  Ers_h1, &
     & rho ,  rho_l1,  rho_h1, &
     dt, tau)

  use rad_params_module, only : ngroups, clight

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: etaT_l1, etaT_h1
  integer, intent(in) :: etTz_l1, etTz_h1
  integer, intent(in) :: eta1_l1, eta1_h1
  integer, intent(in) :: djdT_l1, djdT_h1
  integer, intent(in) :: dkdT_l1, dkdT_h1
  integer, intent(in) :: dedT_l1, dedT_h1
  integer, intent(in) ::  Ers_l1,  Ers_h1
  integer, intent(in) ::  rho_l1,  rho_h1
  real(rt)                       :: etaT(etaT_l1:etaT_h1)
  real(rt)                       :: etTz(etTz_l1:etTz_h1)
  real(rt)                       :: eta1(eta1_l1:eta1_h1)
  real(rt)                       :: djdT(djdT_l1:djdT_h1,0:ngroups-1)
  real(rt)        , intent(in )  :: dkdT(dkdT_l1:dkdT_h1,0:ngroups-1)
  real(rt)        , intent(in )  :: dedT(dedT_l1:dedT_h1)
  real(rt)        , intent(in )  :: Ers ( Ers_l1: Ers_h1,0:ngroups-1)
  real(rt)        , intent(in )  :: rho ( rho_l1: rho_h1)
  real(rt)        , intent(in) :: dt, tau

  integer :: i
  real(rt)         :: cdt, sigma
  real(rt)         :: dZdT(0:ngroups-1), sumdZdT, foo, bar

  sigma = 1.e0_rt + tau
  cdt = clight * dt

  do i = lo(1), hi(1)
     dZdT = djdT(i,:) - dkdT(i,:)*Ers(i,:)
     sumdZdT = sum(dZdT)
     if (sumdZdT .eq. 0.e0_rt) then
        sumdZdT = 1.e-50_rt
     end if
     foo = cdt * sumdZdT
     bar = sigma*rho(i)*dedT(i)
     etaT(i) = foo / (foo + bar)
     etTz(i) = etaT(i) / sumdZdT
     eta1(i) = bar / (foo + bar)
     djdT(i,:) = dZdT / sumdZdT
  end do

end subroutine ca_compute_etat


subroutine ca_compute_emissivity( lo, hi, &
     jg  ,   jg_l1,   jg_h1, &
     djdT, djdT_l1, djdT_h1, &
        T,    T_l1,    T_h1, &
      kap,  kap_l1,  kap_h1, &
     dkdT, dkdT_l1, dkdT_h1, &
     pfc, use_WiensLaw, integrate_Planck, Tf)

  use rad_params_module, only : ngroups, nugroup, dnugroup, xnu,  &
       pi, clight, hplanck, kboltz, arad
  use blackbody_module, only : BdBdTIndefInteg

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in)  :: lo(1), hi(1) 
  integer, intent(in) ::   jg_l1,   jg_h1
  integer, intent(in) :: djdT_l1, djdT_h1
  integer, intent(in) ::    T_l1,    T_h1
  integer, intent(in) ::  kap_l1,  kap_h1
  integer, intent(in) :: dkdT_l1, dkdT_h1
  real(rt)                      :: jg  (  jg_l1:  jg_h1,0:ngroups-1)
  real(rt)                      :: djdT(djdT_l1:djdT_h1,0:ngroups-1)
  real(rt)        , intent(in ) ::    T(   T_l1:   T_h1)
  real(rt)        , intent(in ) ::  kap( kap_l1: kap_h1,0:ngroups-1)
  real(rt)        , intent(in ) :: dkdT(dkdT_l1:dkdT_h1,0:ngroups-1)
  real(rt)        , intent(in) :: pfc(0:ngroups-1)
  integer, intent(in) :: use_WiensLaw, integrate_Planck
  real(rt)        , intent(in) :: Tf

  integer :: i, g
  real(rt)         :: dBdT, Bg
  real(rt)         :: Teff, nu, num, nup, hoverk
  real(rt)         :: cB, Tfix
  real(rt)         :: B0, B1, dBdT0, dBdT1
  real(rt)         :: dnu, nubar, expnubar, cdBdT
  real(rt)         :: xnu_full(0:ngroups)

  if (ngroups .eq. 1) then

     do i = lo(1), hi(1)
        Bg = arad*T(i)**4
        dBdT = 4.e0_rt*arad*T(i)**3
        g = 0
        jg(i,g) = Bg*kap(i,g)
        djdT(i,g) = dkdT(i,g)*Bg + dBdT*kap(i,g)
     end do     

  else if (pfc(0) > 0.e0_rt) then  ! a special case for picket-fence model in Su-Olson test
     do i = lo(1), hi(1)
        Bg = arad*T(i)**4
        dBdT = 4.e0_rt*arad*T(i)**3
        do g=0, ngroups-1
           jg(i,g) = pfc(g) * Bg*kap(i,g)
           djdT(i,g) = pfc(g) * (dkdT(i,g)*Bg + dBdT*kap(i,g))
        end do
     end do
  else if (use_WiensLaw > 0) then

     hoverk = hplanck/kboltz
     cB = 8.*pi*kboltz/clight**3
     
     do g=0, ngroups-1
        nu = nugroup(g)
        num = xnu(g)
        nup = xnu(g+1)
        do i = lo(1), hi(1)
           if (Tf < 0.e0_rt) then
              Tfix = T(i)
           else
              Tfix = Tf
           end if
           dBdT = cB * nu**3 * (exp(-hoverk*num/Tfix) - exp(-hoverk*nup/Tfix))
           Bg = dBdT * T(i)
           jg(i,g) = Bg*kap(i,g)
           djdT(i,g) = dkdT(i,g)*Bg + dBdT*kap(i,g)
        end do
     end do

  else if (integrate_Planck > 0) then

     xnu_full = xnu(0:ngroups)
     xnu_full(0) = 0.e0_rt
     xnu_full(ngroups) = max(xnu(ngroups), 1.e25_rt)

     do i=lo(1), hi(1)
        Teff = max(T(i), 1.e-50_rt)
        call BdBdTIndefInteg(Teff, xnu_full(0), B1, dBdT1)
        do g=0, ngroups-1
           B0 = B1
           dBdT0 = dBdT1
           call BdBdTIndefInteg(Teff, xnu_full(g+1), B1, dBdT1)
           Bg = B1 - B0
           dBdT = dBdT1 - dBdT0

           jg(i,g) = Bg*kap(i,g)
           djdT(i,g) = dkdT(i,g)*Bg + dBdT*kap(i,g)
        end do
     end do

  else

     cB = 8.e0_rt*pi*hplanck / clight**3
     cdBdT = 8.e0_rt*pi*hplanck**2 / (kboltz*clight**3)

     do g=0, ngroups-1
        nu = nugroup(g)
        dnu = dnugroup(g)
        do i=lo(1), hi(1)
           Teff = max(T(i), 1.e-50_rt)
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

           jg(i,g) = Bg*kap(i,g)
           djdT(i,g) = dkdT(i,g)*Bg + dBdT*kap(i,g)
        end do
     end do
     
  end if

end subroutine ca_compute_emissivity


subroutine ca_compute_kappas(lo, hi, &
     stt, stt_l1, stt_h1, &
     T  ,   T_l1,   T_h1, &
     kpp, kpp_l1, kpp_h1, &
     kpr, kpr_l1, kpr_h1, &
     kpT, kpT_l1, kpT_h1, &
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

  integer, intent(in)  :: lo(1), hi(1) 
  integer, intent(in) :: stt_l1, stt_h1 
  integer, intent(in) ::   T_l1,   T_h1
  integer, intent(in) :: kpp_l1, kpp_h1 
  integer, intent(in) :: kpr_l1, kpr_h1
  integer, intent(in) :: kpT_l1, kpT_h1
  integer, intent(in) :: do_stme, use_dkdT
  real(rt)        , intent(in)   :: stt(stt_l1:stt_h1,NVAR)
  real(rt)        , intent(in)   ::   T(  T_l1:  T_h1)
  real(rt)                       :: kpp(kpp_l1:kpp_h1, 0:ngroups-1)
  real(rt)                       :: kpr(kpr_l1:kpr_h1, 0:ngroups-1)
  real(rt)                       :: kpT(kpT_l1:kpT_h1, 0:ngroups-1)
  real(rt)        , intent(in)   :: c_kpp, kpp_m, kpp_n, kpp_p 
  real(rt)        , intent(in)   :: c_kpr, kpr_m, kpr_n, kpr_p 
  real(rt)        , intent(in)   :: c_sct, sct_m, sct_n, sct_p 
  real(rt)        , intent(in)   :: Tfloor

  integer :: i, g
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, nup_kpp, nup_kpr, nup_sct, sct, foo, hnuoverkt, exptmp

  do g=0, ngroups-1
     nup_kpp = nugroup(g)**kpp_p
     nup_kpr = nugroup(g)**kpr_p
     nup_sct = nugroup(g)**sct_p
     do i=lo(1), hi(1)
        Teff = max(T(i), tiny)
        Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
        kpp(i,g) = c_kpp * (stt(i,URHO) ** kpp_m) * (Teff ** (-kpp_n)) * nup_kpp
        if (do_stme > 0) then
           foo = kpp(i,g)
           hnuoverkt = hplanck*nugroup(g)/(k_B*Teff)
           exptmp = exp(-hnuoverkt)
           kpp(i,g) = foo * (1.e0_rt - exptmp)
           kpT(i,g) = foo*(-kpp_n/Teff)*(1.e0_rt-exptmp) - foo*exptmp*(hnuoverkt/Teff)
        else
           kpT(i,g) = kpp(i,g) * (-kpp_n/Teff)           
        end if

        if (use_dkdT.eq.0) then
           kpT(i,g) = 0.e0_rt
        end if

        if (c_kpr < 0.0e0_rt) then
           sct       = c_sct * (stt(i,URHO) ** sct_m) * (Teff ** (-sct_n)) * nup_sct
           kpr(i,g) = kpp(i,g) + sct
        else
           kpr(i,g) = c_kpr * (stt(i,URHO) ** kpr_m) * (Teff ** (-kpr_n)) * nup_kpr
        end if
     end do
  end do

end subroutine ca_compute_kappas


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


subroutine ca_opacs( lo, hi,  &
     Snew,Snew_l1,Snew_h1, &
     T   ,   T_l1,   T_h1, &
     Ts  ,  Ts_l1,  Ts_h1, &
     kpp , kpp_l1, kpp_h1, &
     kpr , kpr_l1, kpr_h1, &
     dkdT,dkdT_l1,dkdT_h1, &
     use_dkdT, validStar, lag_opac) 

  use rad_params_module, only : ngroups, nugroup
  use opacity_table_module, only : get_opacities
  use network, only : naux
  use meth_params_module, only : NVAR, URHO, UFX

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: Snew_l1, Snew_h1 
  integer, intent(in) ::    T_l1,    T_h1
  integer, intent(in) ::   Ts_l1,   Ts_h1
  integer, intent(in) ::  kpp_l1,  kpp_h1 
  integer, intent(in) ::  kpr_l1,  kpr_h1
  integer, intent(in) :: dkdT_l1, dkdT_h1
  real(rt)        , intent(in ) :: Snew(Snew_l1:Snew_h1,NVAR)
  real(rt)        , intent(in ) :: T   (   T_l1:   T_h1)
  real(rt)        , intent(in ) :: Ts  (  Ts_l1:  Ts_h1)
  real(rt)         :: kpp ( kpp_l1: kpp_h1,0:ngroups-1)
  real(rt)         :: kpr ( kpr_l1: kpr_h1,0:ngroups-1)
  real(rt)         :: dkdT(dkdT_l1:dkdT_h1,0:ngroups-1)
  integer, intent(in) :: use_dkdT, validStar, lag_opac

  integer :: i, g
  real(rt)         :: kp, kr, nu, rho, temp, Ye
  real(rt)         :: kp1, kr1
  real(rt)         :: kp2, kr2
  real(rt)         :: dT
  logical :: comp_kp, comp_kr
  real(rt)        , parameter :: fac = 0.5e0_rt, minfrac = 1.e-8_rt

  if (lag_opac .eq. 1) then
     dkdT(lo(1):hi(1),:) = 0.e0_rt
     return
  end if

  do i=lo(1), hi(1)
     
     rho = Snew(i,URHO)
     temp = T(i)
     if (naux > 0) then
        Ye = Snew(i,UFX)
     else
        Ye = 0.e0_rt
     end if

     if (validStar > 0) then
        dT = fac*abs(Ts(i) - T(i))
        dT = max(dT, minfrac*T(i))
     else
        dT = T(i) * 1.e-3_rt + 1.e-50_rt
     end if

     do g=0, ngroups-1
        
        nu = nugroup(g)

        comp_kp = .true.
        comp_kr = .true.
        
        call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)
        kpp(i,g) = kp
        kpr(i,g) = kr 

        if (use_dkdT .eq. 0) then        
           dkdT(i,g) = 0.e0_rt
        else

           comp_kp = .true.
           comp_kr = .false.

           call get_opacities(kp1, kr1, rho, temp-dT, Ye, nu, comp_kp, comp_kr)
           call get_opacities(kp2, kr2, rho, temp+dT, Ye, nu, comp_kp, comp_kr)

           dkdT(i,g) = (kp2-kp1)/(2.e0_rt*dT)

        end if

     end do
  end do

end subroutine ca_opacs


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

subroutine ca_compute_powerlaw_kappa_s( lo, hi, &
     kappa, k_l1, k_h1, &
     u    , u_l1, u_h1, & 
     kappa0, m, n, p,  &
     s0, sm, sn, sp,  &
     Tfloor, kfloor) bind(C, name="ca_compute_powerlaw_kappa_s")

  use rad_params_module, only : ngroups, nugroup
  use meth_params_module, only : NVAR, URHO, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: k_l1, k_h1, &
                         u_l1, u_h1
  real(rt)                     :: kappa(k_l1:k_h1, 0:ngroups-1) 
  real(rt)        , intent(in) ::     u(u_l1:u_h1, NVAR)
  real(rt)        , intent(in) :: kappa0, m, n, p, Tfloor, kfloor
  real(rt)        , intent(in) :: s0, sm, sn, sp

  integer :: i, g
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, kf, sct, nup, nusp

  if (  m.eq.0.e0_rt .and.  n.eq.0.e0_rt .and.  p.eq.0.e0_rt .and. &
       sm.eq.0.e0_rt .and. sn.eq.0.e0_rt .and. sp.eq.0.e0_rt ) then
     kappa(lo(1):hi(1),:) = kappa0 + s0
  else
     do g = 0, ngroups-1
        nup = nugroup(g)**p
        nusp = nugroup(g)**sp
        do i = lo(1), hi(1)
           Teff = max(u(i,UTEMP), tiny)
           Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
           kf = kappa0 * (u(i,URHO) ** m) * (Teff ** (-n)) * nup
           sct = s0 * (u(i,URHO) ** sm) * (Teff ** (-sn)) * nusp
           kappa(i,g) = max(kf+sct, kfloor)
        end do
     end do
  end if
  
end subroutine ca_compute_powerlaw_kappa_s

subroutine ca_compute_powerlaw_kappa( lo, hi,  &
     kappa, k_l1, k_h1, &
     u    , u_l1, u_h1, & 
     kappa0, m, n, p,  &
     Tfloor, kfloor) bind(C, name="ca_compute_powerlaw_kappa")

  use rad_params_module, only : ngroups, nugroup
  use meth_params_module, only : NVAR, URHO, UTEMP

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: k_l1, k_h1, &
                         u_l1, u_h1
  real(rt)                     :: kappa(k_l1:k_h1, 0:ngroups-1) 
  real(rt)        , intent(in) ::     u(u_l1:u_h1, NVAR)
  real(rt)        , intent(in) :: kappa0, m, n, p, Tfloor, kfloor

  integer :: i, g
  real(rt)        , parameter :: tiny = 1.0e-50_rt
  real(rt)         :: Teff, kf, nup

  if (  m.eq.0.e0_rt .and.  n.eq.0.e0_rt .and.  p.eq.0.e0_rt ) then
     kappa(lo(1):hi(1),:) = kappa0
  else
     do g = 0, ngroups-1
        nup = nugroup(g)**p
        do i = lo(1), hi(1)
           Teff = max(u(i,UTEMP), tiny)
           Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
           kf = kappa0 * (u(i,URHO) ** m) * (Teff ** (-n)) * nup
           kappa(i,g) = max(kf, kfloor)
        end do
     end do
  end if
  
end subroutine ca_compute_powerlaw_kappa


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


