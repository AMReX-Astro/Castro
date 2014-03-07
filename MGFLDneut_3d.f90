
subroutine ca_accel_acoe_neut(   &
     eta1,eta1_l1,eta1_l2,eta1_l3,eta1_h1,eta1_h2,eta1_h3, &
     theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3, &
     theY,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3, &
     spc , spc_l1, spc_l2, spc_l3, spc_h1, spc_h2, spc_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     aco , aco_l1, aco_l2, aco_l3, aco_h1, aco_h2, aco_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: eta1_l1,eta1_h1,eta1_l2,eta1_h2,eta1_l3,eta1_h3
  integer, intent(in) :: theT_l1,theT_h1,theT_l2,theT_h2,theT_l3,theT_h3
  integer, intent(in) :: theY_l1,theY_h1,theY_l2,theY_h2,theY_l3,theY_h3
  integer, intent(in) ::  spc_l1, spc_h1, spc_l2, spc_h2, spc_l3, spc_h3
  integer, intent(in) ::  kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer, intent(in) ::  aco_l1, aco_h1, aco_l2, aco_h2, aco_l3, aco_h3
  double precision,intent(in )::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2,eta1_l3:eta1_h3)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(in )::spc ( spc_l1: spc_h1, spc_l2: spc_h2, spc_l3: spc_h3,&
       0:ngroups-1)
  double precision,intent(in )::kap ( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,&
       0:ngroups-1)
  double precision,intent(out)::aco ( aco_l1: aco_h1, aco_l2: aco_h2, aco_l3: aco_h3)
  double precision,intent(in) :: dt, tau

  integer :: i, j, k, g
  double precision :: kbar, kybar, H1, Theta, dt1, foo

  dt1 = (1.d0+tau)/dt

  !$omp parallel do private(i,j,k,g,kbar,kybar,H1,Theta,foo)
  do k = aco_l3, aco_h3
  do j = aco_l2, aco_h2
  do i = aco_l1, aco_h1
     kbar = 0.d0
     kybar = 0.d0
     do g=0, ngroups-1
        foo = spc(i,j,k,g) * kap(i,j,k,g)
        kbar = kbar + foo
        kybar = kybar + foo * erg2rhoYe(g)
     end do

     H1 = eta1(i,j,k)
     Theta = theY(i,j,k) - theT(i,j,k)

     aco(i,j,k) = (H1*kbar - Theta*kybar) * clight + dt1
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_accel_acoe_neut


subroutine ca_accel_rhs_neut(  &
     Ern , Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3,  &
     Erl , Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3,  &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3,  &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3,  &
     etaY,etaY_l1,etaY_l2,etaY_l3,etaY_h1,etaY_h2,etaY_h3,  &
     theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3,  &
     theY,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3,  &
     rhs , rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3,  &
     dt)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer, intent(in) :: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer, intent(in) :: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer, intent(in) ::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer, intent(in) ::etaY_l1,etaY_h1,etaY_l2,etaY_h2,etaY_l3,etaY_h3
  integer, intent(in) ::theT_l1,theT_h1,theT_l2,theT_h2,theT_l3,theT_h3
  integer, intent(in) ::theY_l1,theY_h1,theY_l2,theY_h2,theY_l3,theY_h3
  integer, intent(in) :: rhs_l1, rhs_h1, rhs_l2, rhs_h2, rhs_l3, rhs_h3
  double precision,intent(in )::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,&
       0:ngroups-1)
  double precision,intent(in )::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,&
       0:ngroups-1)
  double precision,intent(in ):: kap( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,&
       0:ngroups-1)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2,etaY_l3:etaY_h3)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(out):: rhs( rhs_l1: rhs_h1, rhs_l2: rhs_h2, rhs_l3: rhs_h3)
  double precision,intent(in) :: dt

  integer :: i, j, k, g
  double precision :: rt, ry, H, Theta, foo

  !$omp parallel do private(i,j,k,g,rt, ry, H, Theta, foo)
  do k = rhs_l3, rhs_h3
  do j = rhs_l2, rhs_h2
  do i = rhs_l1, rhs_h1
     rt = 0.d0
     ry = 0.d0
     do g=0,ngroups-1
        foo = kap(i,j,k,g)*(Ern(i,j,k,g)-Erl(i,j,k,g))
        rt = rt + foo
        ry = ry + foo * erg2rhoYe(g)
     end do

     H = etaT(i,j,k) - etaY(i,j,k)
     Theta = theY(i,j,k) - theT(i,j,k)

     rhs(i,j,k) = clight*(H*rt + Theta*ry)
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_accel_rhs_neut


subroutine ca_accel_spec_neut( lo, hi, &
     Ern , Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3, &
     Erl , Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     etaY,etaY_l1,etaY_l2,etaY_l3,etaY_h1,etaY_h2,etaY_h3, &
     theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3, &
     theY,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     mugY,mugY_l1,mugY_l2,mugY_l3,mugY_h1,mugY_h2,mugY_h3, &
     spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer,intent(in)::etaY_l1,etaY_h1,etaY_l2,etaY_h2,etaY_l3,etaY_h3
  integer,intent(in)::theT_l1,theT_h1,theT_l2,theT_h2,theT_l3,theT_h3
  integer,intent(in)::theY_l1,theY_h1,theY_l2,theY_h2,theY_l3,theY_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  integer,intent(in)::mugY_l1,mugY_h1,mugY_l2,mugY_h2,mugY_l3,mugY_h3
  integer,intent(in)::spec_l1,spec_h1,spec_l2,spec_h2,spec_l3,spec_h3
  double precision,intent(in )::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,&
       0:ngroups-1)
  double precision,intent(in )::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,&
       0:ngroups-1)
  double precision,intent(in )::kap ( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,&
       0:ngroups-1)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2,etaY_l3:etaY_h3)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(in )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,&
       0:ngroups-1)
  double precision,intent(in )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,mugY_l3:mugY_h3,&
       0:ngroups-1)
  double precision,intent(out)::spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,&
       0:ngroups-1)
  double precision,intent(in) :: dt, tau

  integer :: i, j, k, g
  double precision :: cdt1, rt, ry, p, q, r, s, foo, sumeps
  double precision,dimension(0:ngroups-1)::Hg, Tg, epsilon, kapt, kk

  cdt1 = 1.d0/(clight*dt)

  !$omp parallel do private(i,j,k,g,rt,ry,p,q,r,s,foo,sumeps) &
  !$omp private (Hg,Tg,epsilon,kapt,kk)
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     rt = 0.d0
     ry = 0.d0
     do g=0,ngroups-1
        foo = kap(i,j,k,g)*(Ern(i,j,k,g)-Erl(i,j,k,g))
        rt = rt + foo
        ry = ry + foo * erg2rhoYe(g)
     end do

     Hg = mugT(i,j,k,:)*etaT(i,j,k) - mugY(i,j,k,:)*etaY(i,j,k)
     Tg = -mugT(i,j,k,:)*theT(i,j,k) + mugY(i,j,k,:)*theY(i,j,k)

     kapt = kap(i,j,k,:) + (1.d0+tau)*cdt1
     kk = kap(i,j,k,:) / kapt

     p = 1.d0 - sum(Hg*kk)
     q = 1.d0 - sum(Tg*erg2rhoYe*kk)
     r = sum(Hg*erg2rhoYe*kk)
     s = sum(Tg*kk)

     epsilon = ((r*Tg + q*Hg) * rt + (s*Hg + p*Tg) * ry) / kapt

     sumeps = sum(epsilon)
     if (sumeps .eq. 0.d0) then
        sumeps = 1.d-50
     end if
     
     spec(i,j,k,:) = epsilon / sumeps
  end do
  end do
  end do
  !$omp end parallel do
  
end subroutine ca_accel_spec_neut


subroutine ca_check_conv_neut( lo, hi, &
     ren,ren_l1,ren_l2,ren_l3,ren_h1,ren_h2,ren_h3, &
     res,res_l1,res_l2,res_l3,res_h1,res_h2,res_h3, &
     re2,re2_l1,re2_l2,re2_l3,re2_h1,re2_h2,re2_h3, &
     ern,ern_l1,ern_l2,ern_l3,ern_h1,ern_h2,ern_h3, &
     Tmn,Tmn_l1,Tmn_l2,Tmn_l3,Tmn_h1,Tmn_h2,Tmn_h3, &
     Tms,Tms_l1,Tms_l2,Tms_l3,Tms_h1,Tms_h2,Tms_h3, &
     rYn,rYn_l1,rYn_l2,rYn_l3,rYn_h1,rYn_h2,rYn_h3, &
     rYs,rYs_l1,rYs_l2,rYs_l3,rYs_h1,rYs_h2,rYs_h3, &
     rY2,rY2_l1,rY2_l2,rY2_l3,rY2_h1,rY2_h2,rY2_h3, &
     rho,rho_l1,rho_l2,rho_l3,rho_h1,rho_h2,rho_h3, &
     kap,kap_l1,kap_l2,kap_l3,kap_h1,kap_h2,kap_h3, &
     jg , jg_l1, jg_l2, jg_l3, jg_h1, jg_h2, jg_h3, &
     deT,deT_l1,deT_l2,deT_l3,deT_h1,deT_h2,deT_h3, &
     deY,deY_l1,deY_l2,deY_l3,deY_h1,deY_h2,deY_h3, &
     rel_re, abs_re, &
     rel_FT, abs_FT, rel_T, abs_T, &
     rel_FY, abs_FY, rel_Y, abs_Y, &
     dt)
  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer,intent(in)::lo(3),hi(3)
  integer,intent(in)::ren_l1, ren_h1, ren_l2, ren_h2, ren_l3, ren_h3
  integer,intent(in)::res_l1, res_h1, res_l2, res_h2, res_l3, res_h3
  integer,intent(in)::re2_l1, re2_h1, re2_l2, re2_h2, re2_l3, re2_h3
  integer,intent(in)::ern_l1, ern_h1, ern_l2, ern_h2, ern_l3, ern_h3
  integer,intent(in)::Tmn_l1, Tmn_h1, Tmn_l2, Tmn_h2, Tmn_l3, Tmn_h3
  integer,intent(in)::Tms_l1, Tms_h1, Tms_l2, Tms_h2, Tms_l3, Tms_h3
  integer,intent(in)::rYn_l1, rYn_h1, rYn_l2, rYn_h2, rYn_l3, rYn_h3
  integer,intent(in)::rYs_l1, rYs_h1, rYs_l2, rYs_h2, rYs_l3, rYs_h3
  integer,intent(in)::rY2_l1, rY2_h1, rY2_l2, rY2_h2, rY2_l3, rY2_h3
  integer,intent(in)::rho_l1, rho_h1, rho_l2, rho_h2, rho_l3, rho_h3
  integer,intent(in)::kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in):: jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  integer,intent(in)::deT_l1, deT_h1, deT_l2, deT_h2, deT_l3, deT_h3
  integer,intent(in)::deY_l1, deY_h1, deY_l2, deY_h2, deY_l3, deY_h3
  double precision,intent(in)::ren(ren_l1:ren_h1,ren_l2:ren_h2,ren_l3:ren_h3)
  double precision,intent(in)::res(res_l1:res_h1,res_l2:res_h2,res_l3:res_h3)
  double precision,intent(in)::re2(re2_l1:re2_h1,re2_l2:re2_h2,re2_l3:re2_h3)
  double precision,intent(in)::ern(ern_l1:ern_h1,ern_l2:ern_h2,ern_l3:ern_h3,0:ngroups-1)
  double precision,intent(in)::Tmn(Tmn_l1:Tmn_h1,Tmn_l2:Tmn_h2,Tmn_l3:Tmn_h3)
  double precision,intent(in)::Tms(Tms_l1:Tms_h1,Tms_l2:Tms_h2,Tms_l3:Tms_h3)
  double precision,intent(in)::rYn(rYn_l1:rYn_h1,rYn_l2:rYn_h2,rYn_l3:rYn_h3)
  double precision,intent(in)::rYs(rYs_l1:rYs_h1,rYs_l2:rYs_h2,rYs_l3:rYs_h3)
  double precision,intent(in)::rY2(rY2_l1:rY2_h1,rY2_l2:rY2_h2,rY2_l3:rY2_h3)
  double precision,intent(in)::rho(rho_l1:rho_h1,rho_l2:rho_h2,rho_l3:rho_h3)
  double precision,intent(in)::kap(kap_l1:kap_h1,kap_l2:kap_h2,kap_l3:kap_h3,0:ngroups-1)
  double precision,intent(in):: jg( jg_l1: jg_h1, jg_l2: jg_h2, jg_l3: jg_h3,0:ngroups-1)
  double precision,intent(in)::deT(deT_l1:deT_h1,deT_l2:deT_h2,deT_l3:deT_h3)
  double precision,intent(in)::deY(deY_l1:deY_h1,deY_l2:deY_h2,deY_l3:deY_h3)
  double precision,intent(inout)::rel_re, abs_re
  double precision,intent(inout)::rel_FT, abs_FT,rel_T, abs_T
  double precision,intent(inout)::rel_FY, abs_FY,rel_Y, abs_Y
  double precision,intent(in) :: dt

  integer :: i, j, k
  double precision :: chg, relchg, FT, FY, cdt, FTdenom, FYdenom, dTe, dYe

  cdt = clight*dt

  !$omp parallel do private(i,j,k,chg,relchg,FT,FY,FTdenom,FYdenom,dTe,dYe) &
  !$omp reduction(max:rel_re,abs_re,rel_FT,abs_FT,rel_T,abs_T,rel_FY,abs_FY,rel_Y,abs_Y)
  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
  do i=lo(1),hi(1)
     chg = abs(ren(i,j,k) - res(i,j,k))
     relchg = abs(chg/(ren(i,j,k)+1.d-50))
     rel_re = max(rel_re,relchg)
     abs_re = max(abs_re,chg)

     chg = abs(Tmn(i,j,k) - Tms(i,j,k))
     relchg = abs(chg/(Tmn(i,j,k)+1.d-50))
     rel_T = max(rel_T,relchg)
     abs_T = max(abs_T,chg)

     chg = abs(rYn(i,j,k) - rYs(i,j,k))
     relchg = abs(chg/(rYn(i,j,k)+1.d-50))
     rel_Y = max(rel_Y,relchg)
     abs_Y = max(abs_Y,chg/rho(i,j,k))

     FT = abs((ren(i,j,k)-re2(i,j,k)) - cdt*sum(kap(i,j,k,:)*Ern(i,j,k,:)-jg(i,j,k,:)))
     FY = abs((rYn(i,j,k)-rY2(i,j,k)) - cdt*sum(erg2rhoYe*(kap(i,j,k,:)*Ern(i,j,k,:)-jg(i,j,k,:))))

     dTe = Tmn(i,j,k)
     dYe = rYn(i,j,k)/rho(i,j,k)
     FTdenom = rho(i,j,k)*(abs(deT(i,j,k)*dTe)+abs(deY(i,j,k)*dYe))
     FYdenom = rYn(i,j,k)

     rel_FT = max(rel_FT, FT/(FTdenom+1.d-50))
     abs_FT = max(abs_FT, FT)

     rel_FY = max(rel_FY, FY/(FYdenom+1.d-50))
     abs_FY = max(abs_FY, FY)
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_check_conv_neut


subroutine ca_check_conv_er_neut( lo, hi,  &
     Ern , Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3, &
     Erl , Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     etTz,etTz_l1,etTz_l2,etTz_l3,etTz_h1,etTz_h2,etTz_h3, &
     etYz,etYz_l1,etYz_l2,etYz_l3,etYz_h1,etYz_h2,etYz_h3, &
     thTz,thTz_l1,thTz_l2,thTz_l3,thTz_h1,thTz_h2,thTz_h3, &
     thYz,thYz_l1,thYz_l2,thYz_l3,thYz_h1,thYz_h2,thYz_h3, &
     temp,temp_l1,temp_l2,temp_l3,temp_h1,temp_h2,temp_h3, &
     Ye  ,  Ye_l1,  Ye_l2,  Ye_l3,  Ye_h1,  Ye_h2,  Ye_h3, &
     rela, abso, errr, dt)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in)::temp_l1,temp_h1,temp_l2,temp_h2,temp_l3,temp_h3
  integer,intent(in)::  Ye_l1,  Ye_h1,  Ye_l2,  Ye_h2,  Ye_l3,  Ye_h3
  integer,intent(in)::etTz_l1,etTz_h1,etTz_l2,etTz_h2,etTz_l3,etTz_h3
  integer,intent(in)::etYz_l1,etYz_h1,etYz_l2,etYz_h2,etYz_l3,etYz_h3
  integer,intent(in)::thTz_l1,thTz_h1,thTz_l2,thTz_h2,thTz_l3,thTz_h3
  integer,intent(in)::thYz_l1,thYz_h1,thYz_l2,thYz_h2,thYz_l3,thYz_h3
  double precision,intent(in):: Ern( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,&
       0:ngroups-1)
  double precision,intent(in):: Erl( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,&
       0:ngroups-1)
  double precision,intent(in):: kap( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,&
       0:ngroups-1)
  double precision,intent(in)::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2,etTz_l3:etTz_h3)
  double precision,intent(in)::etYz(etYz_l1:etYz_h1,etYz_l2:etYz_h2,etYz_l3:etYz_h3)
  double precision,intent(in)::thTz(thTz_l1:thTz_h1,thTz_l2:thTz_h2,thTz_l3:thTz_h3)
  double precision,intent(in)::thYz(thYz_l1:thYz_h1,thYz_l2:thYz_h2,thYz_l3:thYz_h3)
  double precision,intent(in)::temp(temp_l1:temp_h1,temp_l2:temp_h2,temp_l3:temp_h3)
  double precision,intent(in)::Ye  (  Ye_l1:  Ye_h1,  Ye_l2:  Ye_h2,  Ye_l3:  Ye_h3)
  double precision, intent(inout) :: rela, abso, errr
  double precision, intent(in) :: dt

  integer :: i, j, k, g
  double precision :: chg, tot, cdt, der, kdeT, kdeY, err_T, err_Y, err

  cdt = clight * dt
  !$omp parallel do private(i,j,k,g,chg,tot,der,kdeT,kdeY,err_T,err_Y,err) &
  !$omp reduction(max:rela, abso, errr)
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     chg = 0.d0
     tot = 0.d0
     kdeT = 0.d0
     kdeY = 0.d0
     do g=0,ngroups-1
        der = Ern(i,j,k,g)-Erl(i,j,k,g)
        chg = chg + abs(der)
        tot = tot + abs(Ern(i,j,k,g))
        kdeT = kdeT + kap(i,j,k,g)*der
        kdeY = kdeY + erg2rhoYe(g)*kap(i,j,k,g)*der
     end do
     abso = max(abso, chg)
     rela = max(rela, chg / (tot + 1.d-50))

     err_T =  etTz(i,j,k)*kdeT - thTz(i,j,k)*kdeY 
     err_Y = -etYz(i,j,k)*kdeT + thYz(i,j,k)*kdeY 
     err = max(abs(err_T/(temp(i,j,k)+1.d-50)), abs(err_Y/(Ye(i,j,k)+1.d-50)))
     errr = max(errr, err)
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_check_conv_er_neut


subroutine ca_compute_coupty(  &
     cpt,cpt_l1,cpt_l2,cpt_l3,cpt_h1,cpt_h2,cpt_h3, &
     cpy,cpy_l1,cpy_l2,cpy_l3,cpy_h1,cpy_h2,cpy_h3, &
     kpp,kpp_l1,kpp_l2,kpp_l3,kpp_h1,kpp_h2,kpp_h3, &
     eg , eg_l1, eg_l2, eg_l3, eg_h1, eg_h2, eg_h3, &
     jg , jg_l1, jg_l2, jg_l3, jg_h1, jg_h2, jg_h3)
  
  use rad_params_module, only : ngroups, erg2rhoYe

  implicit none

  integer, intent(in) :: cpt_l1, cpt_h1, cpt_l2, cpt_h2, cpt_l3, cpt_h3
  integer, intent(in) :: cpy_l1, cpy_h1, cpy_l2, cpy_h2, cpy_l3, cpy_h3
  integer, intent(in) :: kpp_l1, kpp_h1, kpp_l2, kpp_h2, kpp_l3, kpp_h3
  integer, intent(in) ::  eg_l1,  eg_h1,  eg_l2,  eg_h2,  eg_l3,  eg_h3
  integer, intent(in) ::  jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  double precision, intent(out):: cpt(cpt_l1:cpt_h1,cpt_l2:cpt_h2,cpt_l3:cpt_h3)
  double precision, intent(out):: cpy(cpy_l1:cpy_h1,cpy_l2:cpy_h2,cpy_l3:cpy_h3)
  double precision, intent(in) :: kpp(kpp_l1:kpp_h1,kpp_l2:kpp_h2,kpp_l3:kpp_h3,0:ngroups-1)
  double precision, intent(in) ::  eg( eg_l1: eg_h1, eg_l2: eg_h2, eg_l3: eg_h3,0:ngroups-1)
  double precision, intent(in) ::  jg( jg_l1: jg_h1, jg_l2: jg_h2, jg_l3: jg_h3,0:ngroups-1)

  integer :: i, j, k, g
  double precision :: foo

  !$omp parallel private(i,j,k,g,foo)

  !$omp do
  do k=cpt_l3, cpt_h3
     do j=cpt_l2, cpt_h2
        do i=cpt_l1, cpt_h1
           cpt(i,j,k) = 0.d0
           cpy(i,j,k) = 0.d0
        end do
     end do
  end do
  !$omp end do nowait

  do g=0, ngroups-1
     !$omp do
     do k=cpt_l3, cpt_h3
     do j=cpt_l2, cpt_h2
     do i=cpt_l1, cpt_h1
        foo = kpp(i,j,k,g) * eg(i,j,k,g) - jg(i,j,k,g)
        cpt(i,j,k) = cpt(i,j,k) + foo
        cpy(i,j,k) = cpy(i,j,k) + erg2rhoYe(g) * foo
     end do
     end do
     end do
     !$omp end do nowait
  end do

  !$omp end parallel

end subroutine ca_compute_coupty


subroutine ca_compute_dedx(  &
     S   ,   S_l1,   S_l2,   S_l3,   S_h1,   S_h2,   S_h3, &
     T   ,   T_l1,   T_l2,   T_l3,   T_h1,   T_h2,   T_h3, &
     Ye  ,  Ye_l1,  Ye_l2,  Ye_l3,  Ye_h1,  Ye_h2,  Ye_h3, &
     Ts  ,  Ts_l1,  Ts_l2,  Ts_l3,  Ts_h1,  Ts_h2,  Ts_h3, &
     Yes , Yes_l1, Yes_l2, Yes_l3, Yes_h1, Yes_h2, Yes_h3, &
     dedT,dedT_l1,dedT_l2,dedT_l3,dedT_h1,dedT_h2,dedT_h3, &
     dedY,dedY_l1,dedY_l2,dedY_l3,dedY_h1,dedY_h2,dedY_h3, &
     validStar)

  use eos_module, only : eos, eos_t, eos_input_rt, Tmin=>table_Tmin, Ymin=>table_Yemin
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UFS

  implicit none

  integer, intent(in) ::    S_l1,    S_h1,    S_l2,    S_h2,    S_l3,    S_h3 
  integer, intent(in) ::    T_l1,    T_h1,    T_l2,    T_h2,    T_l3,    T_h3
  integer, intent(in) ::   Ye_l1,   Ye_h1,   Ye_l2,   Ye_h2,   Ye_l3,   Ye_h3
  integer, intent(in) ::   Ts_l1,   Ts_h1,   Ts_l2,   Ts_h2,   Ts_l3,   Ts_h3
  integer, intent(in) ::  Yes_l1,  Yes_h1,  Yes_l2,  Yes_h2,  Yes_l3,  Yes_h3
  integer, intent(in) :: dedT_l1, dedT_h1, dedT_l2, dedT_h2, dedT_l3, dedT_h3
  integer, intent(in) :: dedY_l1, dedY_h1, dedY_l2, dedY_h2, dedY_l3, dedY_h3
  double precision, intent(in ) :: S   (   S_l1:   S_h1,   S_l2:   S_h2,   S_l3:   S_h3,NVAR)
  double precision, intent(in ) :: T   (   T_l1:   T_h1,   T_l2:   T_h2,   T_l3:   T_h3)
  double precision, intent(in ) :: Ye  (  Ye_l1:  Ye_h1,  Ye_l2:  Ye_h2,  Ye_l3:  Ye_h3)
  double precision, intent(in ) :: Ts  (  Ts_l1:  Ts_h1,  Ts_l2:  Ts_h2,  Ts_l3:  Ts_h3)
  double precision, intent(in ) :: Yes ( Yes_l1: Yes_h1, Yes_l2: Yes_h2, Yes_l3: Yes_h3)
  double precision, intent(out) :: dedT(dedT_l1:dedT_h1,dedT_l2:dedT_h2,dedT_l3:dedT_h3)
  double precision, intent(out) :: dedY(dedY_l1:dedY_h1,dedY_l2:dedY_h2,dedY_l3:dedY_h3)
  integer, intent(in) :: validStar

  integer :: i, j, k
  double precision :: rhoinv, e1, e2
  type(eos_t) :: eos_state
  double precision :: dT, dYe, T1, T2, Ye1, Ye2
  double precision, parameter :: fac = 0.5d0, minfrac = 1.d-8

  !$OMP PARALLEL DO PRIVATE(i,j,k,rhoinv,eos_state,dT,dYe,T1,T2,Ye1,Ye2,e1,e2)
  do k=dedT_l3, dedT_h3
  do j=dedT_l2, dedT_h2
  do i=dedT_l1, dedT_h1

     rhoinv = 1.d0/S(i,j,k,URHO)
     eos_state % rho = S(i,j,k,URHO)
     eos_state % xn  = S(i,j,k,UFS:UFS+nspec-1) * rhoinv
     eos_state % aux = ye(i,j,k)

     if (validStar > 0) then
        dT = fac*abs(Ts(i,j,k) - T(i,j,k))
        dT = max(dT, minfrac*T(i,j,k))
        dYe = fac*abs(Yes(i,j,k) - Ye(i,j,k))
        dYe = max(dYe, minfrac*Ye(i,j,k))
     else
        dT = T(i,j,k) * 1.d-3 + 1.d-50
        dYe = 1.d-4
     end if

     T1 = T(i,j,k) - dT
     if (T1 < Tmin*(1.d0+1.d-6)) then
        T1 = T(i,j,k)
     end if
     T2 = T(i,j,k) + dT

     Ye1 = Ye(i,j,k) - dYe
     if (Ye1 < Ymin*(1.d0+1.d-6)) then
        Ye1 = Ye(i,j,k)
     end if
     Ye2 = Ye(i,j,k) + dYe

     eos_state % T = T1
     call eos(eos_input_rt, eos_state)
     e1 = eos_state % e

     eos_state % T = T2
     call eos(eos_input_rt, eos_state)
     e2 = eos_state % e

     dedT(i,j,k) = (e2-e1)/(T2-T1)

     eos_state % T = T(i,j,k)

     eos_state % aux = Ye1
     call eos(eos_input_rt, eos_state)
     e1 = eos_state % e

     eos_state % aux = Ye2
     call eos(eos_input_rt, eos_state)
     e2 = eos_state % e

     dedY(i,j,k) = (e2-e1)/(Ye2-Ye1)

  end do
  end do
  end do
  !$OMP END PARALLEL DO 

end subroutine ca_compute_dedx


subroutine ca_compute_eta_the( lo, hi, &
     & etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     & etTz,etTz_l1,etTz_l2,etTz_l3,etTz_h1,etTz_h2,etTz_h3, &
     & etaY,etaY_l1,etaY_l2,etaY_l3,etaY_h1,etaY_h2,etaY_h3, &
     & etYz,etYz_l1,etYz_l2,etYz_l3,etYz_h1,etYz_h2,etYz_h3, &
     & eta1,eta1_l1,eta1_l2,eta1_l3,eta1_h1,eta1_h2,eta1_h3, &
     & theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3, &
     & thTz,thTz_l1,thTz_l2,thTz_l3,thTz_h1,thTz_h2,thTz_h3, &
     & theY,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3, &
     & thYz,thYz_l1,thYz_l2,thYz_l3,thYz_h1,thYz_h2,thYz_h3, &
     & the1,the1_l1,the1_l2,the1_l3,the1_h1,the1_h2,the1_h3, &
     & djdT,djdT_l1,djdT_l2,djdT_l3,djdT_h1,djdT_h2,djdT_h3, &
     & djdY,djdY_l1,djdY_l2,djdY_l3,djdY_h1,djdY_h2,djdY_h3, &
     & dkdT,dkdT_l1,dkdT_l2,dkdT_l3,dkdT_h1,dkdT_h2,dkdT_h3, &
     & dkdY,dkdY_l1,dkdY_l2,dkdY_l3,dkdY_h1,dkdY_h2,dkdY_h3, &
     & dedT,dedT_l1,dedT_l2,dedT_l3,dedT_h1,dedT_h2,dedT_h3, &
     & dedY,dedY_l1,dedY_l2,dedY_l3,dedY_h1,dedY_h2,dedY_h3, &
     & Ers , Ers_l1, Ers_l2, Ers_l3, Ers_h1, Ers_h2, Ers_h3, &
     & rho , rho_l1, rho_l2, rho_l3, rho_h1, rho_h2, rho_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: etaT_l1, etaT_h1, etaT_l2, etaT_h2, etaT_l3, etaT_h3
  integer, intent(in) :: etTz_l1, etTz_h1, etTz_l2, etTz_h2, etTz_l3, etTz_h3
  integer, intent(in) :: etaY_l1, etaY_h1, etaY_l2, etaY_h2, etaY_l3, etaY_h3
  integer, intent(in) :: etYz_l1, etYz_h1, etYz_l2, etYz_h2, etYz_l3, etYz_h3
  integer, intent(in) :: eta1_l1, eta1_h1, eta1_l2, eta1_h2, eta1_l3, eta1_h3
  integer, intent(in) :: theT_l1, theT_h1, theT_l2, theT_h2, theT_l3, theT_h3
  integer, intent(in) :: thTz_l1, thTz_h1, thTz_l2, thTz_h2, thTz_l3, thTz_h3
  integer, intent(in) :: theY_l1, theY_h1, theY_l2, theY_h2, theY_l3, theY_h3
  integer, intent(in) :: thYz_l1, thYz_h1, thYz_l2, thYz_h2, thYz_l3, thYz_h3
  integer, intent(in) :: the1_l1, the1_h1, the1_l2, the1_h2, the1_l3, the1_h3
  integer, intent(in) :: djdT_l1, djdT_h1, djdT_l2, djdT_h2, djdT_l3, djdT_h3
  integer, intent(in) :: djdY_l1, djdY_h1, djdY_l2, djdY_h2, djdY_l3, djdY_h3
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2, dkdT_l3, dkdT_h3
  integer, intent(in) :: dkdY_l1, dkdY_h1, dkdY_l2, dkdY_h2, dkdY_l3, dkdY_h3
  integer, intent(in) :: dedT_l1, dedT_h1, dedT_l2, dedT_h2, dedT_l3, dedT_h3
  integer, intent(in) :: dedY_l1, dedY_h1, dedY_l2, dedY_h2, dedY_l3, dedY_h3
  integer, intent(in) ::  Ers_l1,  Ers_h1,  Ers_l2,  Ers_h2,  Ers_l3,  Ers_h3
  integer, intent(in) ::  rho_l1,  rho_h1,  rho_l2,  rho_h2,  rho_l3,  rho_h3
  double precision,intent(out)  ::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  double precision,intent(out)  ::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2,etTz_l3:etTz_h3)
  double precision,intent(out)  ::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2,etaY_l3:etaY_h3)
  double precision,intent(out)  ::etYz(etYz_l1:etYz_h1,etYz_l2:etYz_h2,etYz_l3:etYz_h3)
  double precision,intent(out)  ::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2,eta1_l3:eta1_h3)
  double precision,intent(out)  ::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(out)  ::thTz(thTz_l1:thTz_h1,thTz_l2:thTz_h2,thTz_l3:thTz_h3)
  double precision,intent(out)  ::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(out)  ::thYz(thYz_l1:thYz_h1,thYz_l2:thYz_h2,thYz_l3:thYz_h3)
  double precision,intent(out)  ::the1(the1_l1:the1_h1,the1_l2:the1_h2,the1_l3:the1_h3)
  double precision,intent(inout)::djdT(djdT_l1:djdT_h1,djdT_l2:djdT_h2,djdT_l3:djdT_h3,&
       0:ngroups-1)
  double precision,intent(inout)::djdY(djdY_l1:djdY_h1,djdY_l2:djdY_h2,djdY_l3:djdY_h3,&
       0:ngroups-1)
  double precision,intent(in )  ::dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,dkdT_l3:dkdT_h3,&
       0:ngroups-1)
  double precision,intent(in )  ::dkdY(dkdY_l1:dkdY_h1,dkdY_l2:dkdY_h2,dkdY_l3:dkdY_h3,&
       0:ngroups-1)
  double precision,intent(in )  ::dedT(dedT_l1:dedT_h1,dedT_l2:dedT_h2,dedT_l3:dedT_h3)
  double precision,intent(in )  ::dedY(dedY_l1:dedY_h1,dedY_l2:dedY_h2,dedY_l3:dedY_h3)
  double precision,intent(in )  ::Ers ( Ers_l1: Ers_h1, Ers_l2: Ers_h2, Ers_l3: Ers_h3,&
       0:ngroups-1)
  double precision,intent(in )  ::rho ( rho_l1: rho_h1, rho_l2: rho_h2, rho_l3: rho_h3)
  double precision,intent(in) :: dt, tau

  integer :: i, j, k
  double precision :: cdt, det, et, ey, tt, ty, sigma
  double precision :: dZdT(0:ngroups-1), dZdY(0:ngroups-1)
  double precision :: sumdZdT, sumdZdY, fooT, fooY, barT, barY

  sigma = 1.d0 + tau
  cdt = clight * dt

  !$omp parallel do private(i,j,k,det,et,ey,tt,ty,dZdT,dZdY,sumdZdT,sumdZdY) &
  !$omp private(fooT,fooY,barT,barY)
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     dZdT = djdT(i,j,k,:) - dkdT(i,j,k,:)*Ers(i,j,k,:)
     dZdY = djdY(i,j,k,:) - dkdY(i,j,k,:)*Ers(i,j,k,:)

     sumdZdT = sum(dZdT)
     sumdZdY = sum(dZdY)

     if (sumdZdT .eq. 0.d0) then
        sumdZdT = 1.d-50
     end if

     if (sumdZdY .eq. 0.d0) then
        sumdZdY = 1.d-50
     end if

     fooT = cdt * sumdZdT
     fooY = cdt * sumdZdY
     barT = sigma*rho(i,j,k)*dedT(i,j,k)
     barY = sigma*rho(i,j,k)*dedY(i,j,k)

     et = sigma*rho(i,j,k) + cdt * sum(erg2rhoYe*dZdY)
     ey = cdt * sum(erg2rhoYe*dZdT)
     tt = barY + fooY
     ty = barT + fooT
     det = ty*et - tt*ey
     etaT(i,j,k) = (et * fooT) / det
     etaY(i,j,k) = (ey * fooY) / det
     theT(i,j,k) = (tt * fooT) / det
     theY(i,j,k) = (ty * fooY) / det

     etTz(i,j,k) = etaT(i,j,k) / sumdZdT
     etYz(i,j,k) = etaY(i,j,k) / sumdZdY
     thTz(i,j,k) = theT(i,j,k) / sumdZdT
     thYz(i,j,k) = theY(i,j,k) / sumdZdY

     eta1(i,j,k) = (et*barT - ey*barY) / det
     the1(i,j,k) = (sigma*rho(i,j,k) * ty) / det

     djdT(i,j,k,:) = dZdT / sumdZdT
     djdY(i,j,k,:) = dZdY / sumdZdY
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_compute_eta_the


subroutine ca_compute_rhs_neut(  &
     rhs , rhs_l1, rhs_l2, rhs_l3, rhs_h1, rhs_h2, rhs_h3, &
     jg  ,  jg_l1,  jg_l2,  jg_l3,  jg_h1,  jg_h2,  jg_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     mugY,mugY_l1,mugY_l2,mugY_l3,mugY_h1,mugY_h2,mugY_h3, &
     cpT , cpT_l1, cpT_l2, cpT_l3, cpT_h1, cpT_h2, cpT_h3, &
     cpY , cpY_l1, cpY_l2, cpY_l3, cpY_h1, cpY_h2, cpY_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     etay,etaY_l1,etaY_l2,etaY_l3,etaY_h1,etaY_h2,etaY_h3, &
     theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3, &
     they,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3, &
     Er2 , Er2_l1, Er2_l2, Er2_l3, Er2_h1, Er2_h2, Er2_h3, &
     re2 , re2_l1, re2_l2, re2_l3, re2_h1, re2_h2, re2_h3, &
     rY2 , rY2_l1, rY2_l2, rY2_l3, rY2_h1, rY2_h2, rY2_h3, &
     Ers , Ers_l1, Ers_l2, Ers_l3, Ers_h1, Ers_h2, Ers_h3, &
     res , res_l1, res_l2, res_l3, res_h1, res_h2, res_h3, &
     rYs , rYs_l1, rYs_l2, rYs_l3, rYs_h1, rYs_h2, rYs_h3, &
     r, dt, igroup, tau)

  use rad_params_module, only : ngroups, clight

  implicit none

  integer,intent(in):: rhs_l1, rhs_h1, rhs_l2, rhs_h2, rhs_l3, rhs_h3
  integer,intent(in)::  jg_l1,  jg_h1,  jg_l2,  jg_h2,  jg_l3,  jg_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  integer,intent(in)::mugY_l1,mugY_h1,mugY_l2,mugY_h2,mugY_l3,mugY_h3
  integer,intent(in):: cpT_l1, cpT_h1, cpT_l2, cpT_h2, cpT_l3, cpT_h3
  integer,intent(in):: cpY_l1, cpY_h1, cpY_l2, cpY_h2, cpY_l3, cpY_h3
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer,intent(in)::etaY_l1,etaY_h1,etaY_l2,etaY_h2,etaY_l3,etaY_h3
  integer,intent(in)::theT_l1,theT_h1,theT_l2,theT_h2,theT_l3,theT_h3
  integer,intent(in)::theY_l1,theY_h1,theY_l2,theY_h2,theY_l3,theY_h3
  integer,intent(in):: Er2_l1, Er2_h1, Er2_l2, Er2_h2, Er2_l3, Er2_h3
  integer,intent(in):: re2_l1, re2_h1, re2_l2, re2_h2, re2_l3, re2_h3
  integer,intent(in):: rY2_l1, rY2_h1, rY2_l2, rY2_h2, rY2_l3, rY2_h3
  integer,intent(in):: Ers_l1, Ers_h1, Ers_l2, Ers_h2, Ers_l3, Ers_h3
  integer,intent(in):: res_l1, res_h1, res_l2, res_h2, res_l3, res_h3
  integer,intent(in):: rYs_l1, rYs_h1, rYs_l2, rYs_h2, rYs_l3, rYs_h3
  double precision,intent(out)::rhs ( rhs_l1: rhs_h1, rhs_l2: rhs_h2, rhs_l3: rhs_h3)
  double precision,intent(in )::jg  (  jg_l1:  jg_h1,  jg_l2:  jg_h2,  jg_l3:  jg_h3,&
       0:ngroups-1)
  double precision,intent(in )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,&
       0:ngroups-1)
  double precision,intent(in )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,mugY_l3:mugY_h3,&
       0:ngroups-1)
  double precision,intent(in )::cpT ( cpT_l1: cpT_h1, cpT_l2: cpT_h2, cpT_l3: cpT_h3)
  double precision,intent(in )::cpY ( cpY_l1: cpY_h1, cpY_l2: cpY_h2, cpY_l3: cpY_h3)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2,etaY_l3:etaY_h3)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(in )::Er2 ( Er2_l1: Er2_h1, Er2_l2: Er2_h2, Er2_l3: Er2_h3,&
       0:ngroups-1)
  double precision,intent(in )::re2 ( re2_l1: re2_h1, re2_l2: re2_h2, re2_l3: re2_h3)
  double precision,intent(in )::rY2 ( rY2_l1: rY2_h1, rY2_l2: rY2_h2, rY2_l3: rY2_h3)
  double precision,intent(in )::Ers ( Ers_l1: Ers_h1, Ers_l2: Ers_h2, Ers_l3: Ers_h3,&
       0:ngroups-1)
  double precision,intent(in )::res ( res_l1: res_h1, res_l2: res_h2, res_l3: res_h3)
  double precision,intent(in )::rYs ( rYs_l1: rYs_h1, rYs_l2: rYs_h2, rYs_l3: rYs_h3)
  double precision,intent(in) ::   r( rhs_l1: rhs_h1)
  double precision,intent(in) :: dt, tau             
  integer, intent(in) :: igroup

  integer :: i, j, k
  double precision :: Hg, thetag, dt1

  dt1 = 1.d0/dt
  !$omp parallel do private(i,j,k,Hg,thetag)
  do k=rhs_l3, rhs_h3
  do j=rhs_l2, rhs_h2
  do i=rhs_l1, rhs_h1
     Hg = mugT(i,j,k,igroup)*etaT(i,j,k) - mugY(i,j,k,igroup)*etaY(i,j,k)
     thetag = mugY(i,j,k,igroup)*theY(i,j,k) - mugT(i,j,k,igroup)*theT(i,j,k)

     rhs(i,j,k) = clight*(jg(i,j,k,igroup) + Hg*cpT(i,j,k) + thetag*cpY(i,j,k)) &
          + dt1 * (Er2(i,j,k,igroup) - Hg*(res(i,j,k)-re2(i,j,k))  &
          &                        - thetag*(rYs(i,j,k)-rY2(i,j,k)) &
          &        + tau*Ers(i,j,k,igroup))
   end do
   end do
   end do
   !$omp end parallel do

end subroutine ca_compute_rhs_neut


subroutine ca_local_accel_neut( lo, hi,  &
     Ern , Ern_l1, Ern_l2, Ern_l3, Ern_h1, Ern_h2, Ern_h3, &
     Erl , Erl_l1, Erl_l2, Erl_l3, Erl_h1, Erl_h2, Erl_h3, &
     kap , kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     etaY,etaY_l1,etaY_l2,etaY_l3,etaY_h1,etaY_h2,etaY_h3, &
     theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3, &
     theY,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     mugY,mugY_l1,mugY_l2,mugY_l3,mugY_h1,mugY_h2,mugY_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, clight, erg2rhoYe

  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: Ern_l1, Ern_h1, Ern_l2, Ern_h2, Ern_l3, Ern_h3
  integer,intent(in):: Erl_l1, Erl_h1, Erl_l2, Erl_h2, Erl_l3, Erl_h3
  integer,intent(in):: kap_l1, kap_h1, kap_l2, kap_h2, kap_l3, kap_h3
  integer,intent(in)::etaT_l1,etaT_h1,etaT_l2,etaT_h2,etaT_l3,etaT_h3
  integer,intent(in)::etaY_l1,etaY_h1,etaY_l2,etaY_h2,etaY_l3,etaY_h3
  integer,intent(in)::theT_l1,theT_h1,theT_l2,theT_h2,theT_l3,theT_h3
  integer,intent(in)::theY_l1,theY_h1,theY_l2,theY_h2,theY_l3,theY_h3
  integer,intent(in)::mugT_l1,mugT_h1,mugT_l2,mugT_h2,mugT_l3,mugT_h3
  integer,intent(in)::mugY_l1,mugY_h1,mugY_l2,mugY_h2,mugY_l3,mugY_h3
  double precision,intent(inout)::Ern ( Ern_l1: Ern_h1, Ern_l2: Ern_h2, Ern_l3: Ern_h3,&
       0:ngroups-1)
  double precision,intent(in   )::Erl ( Erl_l1: Erl_h1, Erl_l2: Erl_h2, Erl_l3: Erl_h3,&
       0:ngroups-1)
  double precision,intent(in   )::kap ( kap_l1: kap_h1, kap_l2: kap_h2, kap_l3: kap_h3,&
       0:ngroups-1)
  double precision,intent(in   )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  double precision,intent(in   )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2,etaY_l3:etaY_h3)
  double precision,intent(in   )::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(in   )::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(in   )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,&
       0:ngroups-1)
  double precision,intent(in   )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,mugY_l3:mugY_h3,&
       0:ngroups-1)
  double precision,intent(in) :: dt, tau

  integer :: i, j, k
  double precision :: cdt1, rt, ry, p, q, r, s 
  double precision,dimension(0:ngroups-1)::Hg, Tg, epsilon, kapt, kk

  cdt1 = 1.d0/(clight*dt)

  !$omp parallel do private(i,j,k,rt,ry,p,q,r,s ,Hg,Tg,epsilon,kapt,kk)
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rt = sum(kap(i,j,k,:)*(Ern(i,j,k,:)-Erl(i,j,k,:)))
     ry = sum(erg2rhoYe(:)*kap(i,j,k,:)*(Ern(i,j,k,:)-Erl(i,j,k,:)))

     Hg = mugT(i,j,k,:)*etaT(i,j,k) - mugY(i,j,k,:)*etaY(i,j,k)
     Tg = -mugT(i,j,k,:)*theT(i,j,k) + mugY(i,j,k,:)*theY(i,j,k)

     kapt = kap(i,j,k,:) + (1.d0+tau)*cdt1
     kk = kap(i,j,k,:) / kapt

     p = 1.d0-sum(Hg*kk)
     q = 1.d0-sum(Tg*erg2rhoYe*kk)
     r = sum(Hg*erg2rhoYe*kk)
     s = sum(Tg*kk)

     epsilon = ((r*Tg + q*Hg) * rt + (s*Hg + p*Tg) * ry)  &
          / (kapt*(p*q-r*s) + 1.d-50)

     Ern(i,j,k,:) = Ern(i,j,k,:) + epsilon
  end do
  end do
  end do
  !$omp end parallel do
  
end subroutine ca_local_accel_neut


subroutine ca_opac_emis_neut( lo, hi,  &
     Snew,Snew_l1,Snew_l2,Snew_l3,Snew_h1,Snew_h2,Snew_h3, &
     T   ,   T_l1,   T_l2,   T_l3,   T_h1,   T_h2,   T_h3, &
     Ye  ,  Ye_l1,  Ye_l2,  Ye_l3,  Ye_h1,  Ye_h2,  Ye_h3, &
     Ts  ,  Ts_l1,  Ts_l2,  Ts_l3,  Ts_h1,  Ts_h2,  Ts_h3, &
     Yes , Yes_l1, Yes_l2, Yes_l3, Yes_h1, Yes_h2, Yes_h3, &
     kpp , kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     kpr , kpr_l1, kpr_l2, kpr_l3, kpr_h1, kpr_h2, kpr_h3, &
     jg  ,   j_l1,   j_l2,   j_l3,   j_h1,   j_h2,   j_h3, &
     djdT,djdT_l1,djdT_l2,djdT_l3,djdT_h1,djdT_h2,djdT_h3, &
     djdY,djdY_l1,djdY_l2,djdY_l3,djdY_h1,djdY_h2,djdY_h3, &
     dkdT,dkdT_l1,dkdT_l2,dkdT_l3,dkdT_h1,dkdT_h2,dkdT_h3, &
     dkdY,dkdY_l1,dkdY_l2,dkdY_l3,dkdY_h1,dkdY_h2,dkdY_h3, &
     use_dkdT, validStar, lag_opac) 

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: Snew_l1, Snew_h1, Snew_l2, Snew_h2, Snew_l3, Snew_h3 
  integer, intent(in) ::    T_l1,    T_h1,    T_l2,    T_h2,    T_l3,    T_h3
  integer, intent(in) ::   Ye_l1,   Ye_h1,   Ye_l2,   Ye_h2,   Ye_l3,   Ye_h3
  integer, intent(in) ::   Ts_l1,   Ts_h1,   Ts_l2,   Ts_h2,   Ts_l3,   Ts_h3
  integer, intent(in) ::  Yes_l1,  Yes_h1,  Yes_l2,  Yes_h2,  Yes_l3,  Yes_h3
  integer, intent(in) ::  kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2,  kpp_l3,  kpp_h3 
  integer, intent(in) ::  kpr_l1,  kpr_h1,  kpr_l2,  kpr_h2,  kpr_l3,  kpr_h3
  integer, intent(in) ::    j_l1,    j_h1,    j_l2,    j_h2,    j_l3,    j_h3
  integer, intent(in) :: djdT_l1, djdT_h1, djdT_l2, djdT_h2, djdT_l3, djdT_h3
  integer, intent(in) :: djdY_l1, djdY_h1, djdY_l2, djdY_h2, djdY_l3, djdY_h3
  integer, intent(in) :: dkdT_l1, dkdT_h1, dkdT_l2, dkdT_h2, dkdT_l3, dkdT_h3
  integer, intent(in) :: dkdY_l1, dkdY_h1, dkdY_l2, dkdY_h2, dkdY_l3, dkdY_h3
  double precision, intent(in ) :: Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,Snew_l3:Snew_h3,NVAR)
  double precision, intent(in ) :: T   (   T_l1:   T_h1,   T_l2:   T_h2,   T_l3:   T_h3)
  double precision, intent(in ) :: Ye  (  Ye_l1:  Ye_h1,  Ye_l2:  Ye_h2,  Ye_l3:  Ye_h3)
  double precision, intent(in ) :: Ts  (  Ts_l1:  Ts_h1,  Ts_l2:  Ts_h2,  Ts_l3:  Ts_h3)
  double precision, intent(in ) :: Yes ( Yes_l1: Yes_h1, Yes_l2: Yes_h2, Yes_l3: Yes_h3)
  double precision, intent(inout) :: kpp ( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,&
       0:ngroups-1)
  double precision, intent(inout) :: kpr ( kpr_l1: kpr_h1, kpr_l2: kpr_h2, kpr_l3: kpr_h3,&
       0:ngroups-1)
  double precision, intent(out  ) :: jg  (   j_l1:   j_h1,   j_l2:   j_h2,   j_l3:   j_h3,&
       0:ngroups-1)
  double precision, intent(out  ) :: djdT(djdT_l1:djdT_h1,djdT_l2:djdT_h2,djdT_l3:djdT_h3,&
       0:ngroups-1)
  double precision, intent(out  ) :: djdY(djdY_l1:djdY_h1,djdY_l2:djdY_h2,djdY_l3:djdY_h3,&
       0:ngroups-1)
  double precision, intent(inout) :: dkdT(dkdT_l1:dkdT_h1,dkdT_l2:dkdT_h2,dkdT_l3:dkdT_h3,&
       0:ngroups-1)
  double precision, intent(inout) :: dkdY(dkdY_l1:dkdY_h1,dkdY_l2:dkdY_h2,dkdY_l3:dkdY_h3,&
       0:ngroups-1)
  integer, intent(in) :: use_dkdT, validStar, lag_opac

  integer :: i, j, k, g, inu
  double precision :: ab, sc, delta, eta, er, der, rho, temp
  double precision :: ab1, sc1, delta1, eta1
  double precision :: ab2, sc2, delta2, eta2 
  double precision :: dT, dYe
  double precision :: Bg, Bg1, Bg2
  logical :: comp_ab, comp_sc, comp_eta
  double precision, parameter :: fac = 0.5d0, minfrac = 1.d-8

  !$OMP PARALLEL DO PRIVATE(i,j,k,g,rho,temp,dT,dYe,inu,er,der,comp_ab,comp_sc,comp_eta) &
  !$OMP PRIVATE(ab ,sc ,delta ,eta ,Bg ) &
  !$OMP PRIVATE(ab1,sc1,delta1,eta1,Bg1) &
  !$OMP PRIVATE(ab2,sc2,delta2,eta2,Bg2)
  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     
     rho = Snew(i,j,k,URHO)
     temp = T(i,j,k)

     if (validStar > 0) then
        dT = fac*abs(Ts(i,j,k) - T(i,j,k))
        dT = max(dT, minfrac*T(i,j,k))
        dYe = fac*abs(Yes(i,j,k) - Ye(i,j,k))
        dYe = max(dYe, minfrac*Ye(i,j,k))
     else
        dT = T(i,j,k) * 1.d-3 + 1.d-50
        dYe = 1.d-4
     end if

     do g=0, ngroups-1
        
        call prep_opacity(g, inu, er, der)

        if (lag_opac .eq. 1) then

           dkdT(i,j,k,g) = 0.d0
           dkdY(i,j,k,g) = 0.d0              

           comp_ab = .true.
           comp_sc = .false.
           comp_eta = .true.

           call get_opacity_emissivity(ab, sc, delta, eta, &
                rho, Ye(i,j,k), temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg = eta * der / ab
           jg(i,j,k,g) = Bg * kpp(i,j,k,g)

           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j,k)-dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg1 = eta1 * der / ab1

           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j,k)+dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           Bg2 = eta2 * der / ab2

           djdY(i,j,k,g) = (Bg2-Bg1)*kpp(i,j,k,g)/(2.d0*dye)

           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j,k), temp-dT, er, inu, comp_ab, comp_sc, comp_eta)
           Bg1 = eta1 * der / ab1
           
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j,k), temp+dT, er, inu, comp_ab, comp_sc, comp_eta)
           Bg2 = eta2 * der / ab2

           djdT(i,j,k,g) = (Bg2-Bg1)*kpp(i,j,k,g)/(2.d0*dT)
           
        else

           comp_ab = .true.
           comp_sc = .true.
           comp_eta = .true.
           
           call get_opacity_emissivity(ab, sc, delta, eta, &
                rho, Ye(i,j,k), temp, er, inu, comp_ab, comp_sc, comp_eta)
           kpp(i,j,k,g) = ab
           kpr(i,j,k,g) = ab + sc * (1.d0 - delta/3.d0)
           jg(i,j,k,g) = eta * der 
           
           comp_ab = .true.
           comp_sc = .false.
           comp_eta = .true.
           
           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j,k)-dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j,k)+dye, temp, er, inu, comp_ab, comp_sc, comp_eta)
           djdY(i,j,k,g) = (eta2-eta1)*der/(2.d0*dye)
           dkdY(i,j,k,g) = (ab2-ab1)/(2.d0*dye)
           if (use_dkdT .eq. 0) then
              djdY(i,j,k,g) = (eta2/ab2 - eta1/ab1)*der/(2.d0*dye) * ab
              dkdY(i,j,k,g) = 0.d0
           end if
           
           call get_opacity_emissivity(ab1, sc1, delta1, eta1, &
                rho, Ye(i,j,k), temp-dT, er, inu, comp_ab, comp_sc, comp_eta)
           call get_opacity_emissivity(ab2, sc2, delta2, eta2, &
                rho, Ye(i,j,k), temp+dT, er, inu, comp_ab, comp_sc, comp_eta)
           djdT(i,j,k,g) = (eta2-eta1)*der/(2.d0*dT)
           dkdT(i,j,k,g) = (ab2-ab1)/(2.d0*dT)
           if (use_dkdT .eq. 0) then        
              djdT(i,j,k,g) = (eta2/ab2 - eta1/ab1)*der/(2.d0*dT) * ab
              dkdT(i,j,k,g) = 0.d0
           end if

        end if

     end do
  end do
  end do
  end do
  !$OMP END PARALLEL DO 

end subroutine ca_opac_emis_neut


subroutine ca_state_update_neut( lo, hi, &
     state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
     rhoe,  rhoe_l1, rhoe_l2, rhoe_l3, rhoe_h1, rhoe_h2, rhoe_h3, &
     Ye  ,    Ye_l1,   Ye_l2,   Ye_l3,   Ye_h1,   Ye_h2,   Ye_h3, &
     temp,  temp_l1, temp_l2, temp_l3, temp_h1, temp_h2, temp_h3, &
     msk ,   msk_l1,  msk_l2,  msk_l3,  msk_h1,  msk_h2,  msk_h3, &
     derat, dTrat, dye)

  use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, UTEMP, UFS, UFX

  implicit none

  integer, intent(in) :: lo(3), hi(3) 
  integer, intent(in) :: state_l1, state_h1, state_l2, state_h2, state_l3, state_h3
  integer, intent(in) ::  rhoe_l1,  rhoe_h1,  rhoe_l2,  rhoe_h2,  rhoe_l3,  rhoe_h3
  integer, intent(in) ::    Ye_l1,    Ye_h1,    Ye_l2,    Ye_h2,    Ye_l3,    Ye_h3
  integer, intent(in) ::  temp_l1,  temp_h1,  temp_l2,  temp_h2,  temp_l3,  temp_h3
  integer, intent(in) ::   msk_l1,   msk_h1,   msk_l2,   msk_h2,   msk_l3,   msk_h3
  double precision,intent(in)   :: rhoe( rhoe_l1: rhoe_h1, rhoe_l2: rhoe_h2, rhoe_l3: rhoe_h3)
  double precision,intent(in)   ::   Ye(   Ye_l1:   Ye_h1,   Ye_l2:   Ye_h2,   Ye_l3:   Ye_h3)
  double precision,intent(in)   :: temp( temp_l1: temp_h1, temp_l2: temp_h2, temp_l3: temp_h3)
  double precision,intent(in)   ::  msk(  msk_l1:  msk_h1,  msk_l2:  msk_h2,  msk_l3:  msk_h3)
  double precision,intent(inout)::state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3&
       ,NVAR)
  double precision, intent(inout) :: derat, dTrat, dye

  integer :: i, j, k
  double precision :: ei, ek, Told, Yeold

  !$omp parallel do private(i,j,k,ei,ek,Told,Yeold) reduction(max:derat,dTrat,dye)
  do k=lo(3), hi(3)
  do j=lo(2), hi(2)
  do i=lo(1), hi(1)
     ei = state(i,j,k,UEINT)
     derat = max(derat, abs((rhoe(i,j,k) - ei)*msk(i,j,k)/ max(ei, 1.d-50)))
     ek = state(i,j,k,UEDEN) - state(i,j,k,UEINT)
     state(i,j,k,UEINT) = rhoe(i,j,k)
     state(i,j,k,UEDEN) = rhoe(i,j,k) + ek

     Told = state(i,j,k,UTEMP)
     dTrat = max(dTrat, abs((temp(i,j,k)-Told)*msk(i,j,k)/ max(Told, 1.d-50)))
     state(i,j,k,UTEMP) = temp(i,j,k)

     Yeold = state(i,j,k,UFX) / state(i,j,k,URHO)
     dye = max(dye, abs((Ye(i,j,k)-Yeold)*msk(i,j,k)))
     state(i,j,k,UFX) = state(i,j,k,URHO) * Ye(i,j,k)
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_state_update_neut


subroutine ca_update_matter_neut( lo, hi,  &
     re_n,re_n_l1,re_n_l2,re_n_l3,re_n_h1,re_n_h2,re_n_h3, &
     rY_n,rY_n_l1,rY_n_l2,rY_n_l3,rY_n_h1,rY_n_h2,rY_n_h3, &
     Ye_n,Ye_n_l1,Ye_n_l2,Ye_n_l3,Ye_n_h1,Ye_n_h2,Ye_n_h3, &
     Er_n,Er_n_l1,Er_n_l2,Er_n_l3,Er_n_h1,Er_n_h2,Er_n_h3, &
     Er_l,Er_l_l1,Er_l_l2,Er_l_l3,Er_l_h1,Er_l_h2,Er_l_h3, &
     re_s,re_s_l1,re_s_l2,re_s_l3,re_s_h1,re_s_h2,re_s_h3, &
     rY_s,rY_s_l1,rY_s_l2,rY_s_l3,rY_s_h1,rY_s_h2,rY_s_h3, &
     re_2,re_2_l1,re_2_l2,re_2_l3,re_2_h1,re_2_h2,re_2_h3, &
     rY_2,rY_2_l1,rY_2_l2,rY_2_l3,rY_2_h1,rY_2_h2,rY_2_h3, &
     etaT,etaT_l1,etaT_l2,etaT_l3,etaT_h1,etaT_h2,etaT_h3, &
     etaY,etaY_l1,etaY_l2,etaY_l3,etaY_h1,etaY_h2,etaY_h3, &
     eta1,eta1_l1,eta1_l2,eta1_l3,eta1_h1,eta1_h2,eta1_h3, &
     theT,theT_l1,theT_l2,theT_l3,theT_h1,theT_h2,theT_h3, &
     theY,theY_l1,theY_l2,theY_l3,theY_h1,theY_h2,theY_h3, &
     the1,the1_l1,the1_l2,the1_l3,the1_h1,the1_h2,the1_h3, &
      cpt, cpt_l1, cpt_l2, cpt_l3, cpt_h1, cpt_h2, cpt_h3, &
      cpy, cpy_l1, cpy_l2, cpy_l3, cpy_h1, cpy_h2, cpy_h3, &
      kpp, kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     mugT,mugT_l1,mugT_l2,mugT_l3,mugT_h1,mugT_h2,mugT_h3, &
     mugY,mugY_l1,mugY_l2,mugY_l3,mugY_h1,mugY_h2,mugY_h3, &
     Snew,Snew_l1,Snew_l2,Snew_l3,Snew_h1,Snew_h2,Snew_h3, &
     dt, tau)

  use rad_params_module, only : ngroups, erg2rhoYe, clight
  use meth_params_module, only : NVAR, URHO

  implicit none

  integer,intent(in)::lo(3),hi(3)
  integer,intent(in)::re_n_l1, re_n_h1, re_n_l2, re_n_h2, re_n_l3, re_n_h3
  integer,intent(in)::rY_n_l1, rY_n_h1, rY_n_l2, rY_n_h2, rY_n_l3, rY_n_h3
  integer,intent(in)::Ye_n_l1, Ye_n_h1, Ye_n_l2, Ye_n_h2, Ye_n_l3, Ye_n_h3
  integer,intent(in)::Er_n_l1, Er_n_h1, Er_n_l2, Er_n_h2, Er_n_l3, Er_n_h3
  integer,intent(in)::Er_l_l1, Er_l_h1, Er_l_l2, Er_l_h2, Er_l_l3, Er_l_h3
  integer,intent(in)::re_s_l1, re_s_h1, re_s_l2, re_s_h2, re_s_l3, re_s_h3
  integer,intent(in)::rY_s_l1, rY_s_h1, rY_s_l2, rY_s_h2, rY_s_l3, rY_s_h3
  integer,intent(in)::re_2_l1, re_2_h1, re_2_l2, re_2_h2, re_2_l3, re_2_h3
  integer,intent(in)::rY_2_l1, rY_2_h1, rY_2_l2, rY_2_h2, rY_2_l3, rY_2_h3
  integer,intent(in)::etaT_l1, etaT_h1, etaT_l2, etaT_h2, etaT_l3, etaT_h3
  integer,intent(in)::etaY_l1, etaY_h1, etaY_l2, etaY_h2, etaY_l3, etaY_h3
  integer,intent(in)::eta1_l1, eta1_h1, eta1_l2, eta1_h2, eta1_l3, eta1_h3
  integer,intent(in)::theT_l1, theT_h1, theT_l2, theT_h2, theT_l3, theT_h3
  integer,intent(in)::theY_l1, theY_h1, theY_l2, theY_h2, theY_l3, theY_h3
  integer,intent(in)::the1_l1, the1_h1, the1_l2, the1_h2, the1_l3, the1_h3
  integer,intent(in):: cpt_l1,  cpt_h1,  cpt_l2,  cpt_h2,  cpt_l3,  cpt_h3
  integer,intent(in):: cpy_l1,  cpy_h1,  cpy_l2,  cpy_h2,  cpy_l3,  cpy_h3
  integer,intent(in):: kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2,  kpp_l3,  kpp_h3
  integer,intent(in)::mugT_l1, mugT_h1, mugT_l2, mugT_h2, mugT_l3, mugT_h3
  integer,intent(in)::mugY_l1, mugY_h1, mugY_l2, mugY_h2, mugY_l3, mugY_h3
  integer,intent(in)::Snew_l1, Snew_h1, Snew_l2, Snew_h2, Snew_l3, Snew_h3
  double precision,intent(out)::re_n(re_n_l1:re_n_h1,re_n_l2:re_n_h2,re_n_l3:re_n_h3)
  double precision,intent(out)::rY_n(rY_n_l1:rY_n_h1,rY_n_l2:rY_n_h2,rY_n_l3:rY_n_h3)
  double precision,intent(out)::Ye_n(Ye_n_l1:Ye_n_h1,Ye_n_l2:Ye_n_h2,Ye_n_l3:Ye_n_h3)
  double precision,intent(in )::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,Er_n_l3:Er_n_h3,&
       0:ngroups-1)
  double precision,intent(in )::Er_l(Er_l_l1:Er_l_h1,Er_l_l2:Er_l_h2,Er_l_l3:Er_l_h3,&
       0:ngroups-1)
  double precision,intent(in )::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2,re_s_l3:re_s_h3)
  double precision,intent(in )::rY_s(rY_s_l1:rY_s_h1,rY_s_l2:rY_s_h2,rY_s_l3:rY_s_h3)
  double precision,intent(in )::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2,re_2_l3:re_2_h3)
  double precision,intent(in )::rY_2(rY_2_l1:rY_2_h1,rY_2_l2:rY_2_h2,rY_2_l3:rY_2_h3)
  double precision,intent(in )::etaT(etaT_l1:etaT_h1,etaT_l2:etaT_h2,etaT_l3:etaT_h3)
  double precision,intent(in )::etaY(etaY_l1:etaY_h1,etaY_l2:etaY_h2,etaY_l3:etaY_h3)
  double precision,intent(in )::eta1(eta1_l1:eta1_h1,eta1_l2:eta1_h2,eta1_l3:eta1_h3)
  double precision,intent(in )::theT(theT_l1:theT_h1,theT_l2:theT_h2,theT_l3:theT_h3)
  double precision,intent(in )::theY(theY_l1:theY_h1,theY_l2:theY_h2,theY_l3:theY_h3)
  double precision,intent(in )::the1(the1_l1:the1_h1,the1_l2:the1_h2,the1_l3:the1_h3)
  double precision,intent(in ):: cpt( cpt_l1: cpt_h1, cpt_l2: cpt_h2, cpt_l3: cpt_h3)
  double precision,intent(in ):: cpy( cpy_l1: cpy_h1, cpy_l2: cpy_h2, cpy_l3: cpy_h3)
  double precision,intent(in ):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,&
       0:ngroups-1)
  double precision,intent(in )::mugT(mugT_l1:mugT_h1,mugT_l2:mugT_h2,mugT_l3:mugT_h3,&
       0:ngroups-1)
  double precision,intent(in )::mugY(mugY_l1:mugY_h1,mugY_l2:mugY_h2,mugY_l3:mugY_h3,&
       0:ngroups-1)
  double precision,intent(in )::Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,Snew_l3:Snew_h3,NVAR)
  double precision,intent(in) :: dt, tau

  integer :: i,j,k,g
  double precision :: cdt, H1, Theta, Thbar1, Hbar 
  double precision :: dkEE, dkEEY, foo, chg, chgY

  cdt = clight * dt
  !$omp parallel do private(i,j,k,g,H1,Theta,Thbar1,Hbar,dkEE,dkEEY,foo,chg,chgY)
  do k = lo(3), hi(3)
  do j = lo(2), hi(2)
  do i = lo(1), hi(1)

     H1 = eta1(i,j,k)
     Thbar1 = the1(i,j,k)

     Theta = theY(i,j,k) - theT(i,j,k)
     Hbar = sum(erg2rhoYe*(mugT(i,j,k,:)*etaT(i,j,k) - mugY(i,j,k,:)*etaY(i,j,k)))

     dkEE = 0.d0
     dkEEY = 0.d0
     do g=0, ngroups-1
        foo = kpp(i,j,k,g)*(Er_n(i,j,k,g)-Er_l(i,j,k,g))
        dkEE = dkEE + foo
        dkEEY = dkEEY + foo * erg2rhoYe(g)
     end do

     chg = cdt*dkEE + H1*((re_2(i,j,k)-re_s(i,j,k)) + cdt*cpt(i,j,k)) &
          - Theta*((rY_2(i,j,k)-rY_s(i,j,k)) + cdt*cpy(i,j,k))

     chgY = cdt*dkEEY + Thbar1*((rY_2(i,j,k)-rY_s(i,j,k)) + cdt*cpy(i,j,k)) &
          - Hbar*((re_2(i,j,k)-re_s(i,j,k)) + cdt*cpt(i,j,k))

     re_n(i,j,k) = re_s(i,j,k) + chg
     rY_n(i,j,k) = rY_s(i,j,k) + chgY

     re_n(i,j,k) = (re_n(i,j,k) + tau*re_s(i,j,k)) / (1.d0+tau)
     rY_n(i,j,k) = (rY_n(i,j,k) + tau*rY_s(i,j,k)) / (1.d0+tau)

     Ye_n(i,j,k) = rY_n(i,j,k) / Snew(i,j,k,URHO)

     ! temperature will be updated after exiting this subroutine
  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_update_matter_neut


subroutine ca_ncupdate_matter_neut( lo, hi,  &
     Tp_n, Tp_n_l1, Tp_n_l2, Tp_n_l3, Tp_n_h1, Tp_n_h2, Tp_n_h3,  &
     Ye_n, Ye_n_l1, Ye_n_l2, Ye_n_l3, Ye_n_h1, Ye_n_h2, Ye_n_h3,  &
     Er_n, Er_n_l1, Er_n_l2, Er_n_l3, Er_n_h1, Er_n_h2, Er_n_h3,  &
     re_s, re_s_l1, re_s_l2, re_s_l3, re_s_h1, re_s_h2, re_s_h3,  &
     rY_s, rY_s_l1, rY_s_l2, rY_s_l3, rY_s_h1, rY_s_h2, rY_s_h3,  &
     re_2, re_2_l1, re_2_l2, re_2_l3, re_2_h1, re_2_h2, re_2_h3,  &
     rY_2, rY_2_l1, rY_2_l2, rY_2_l3, rY_2_h1, rY_2_h2, rY_2_h3,  &
     etTz, etTz_l1, etTz_l2, etTz_l3, etTz_h1, etTz_h2, etTz_h3,  &
     etYz, etYz_l1, etYz_l2, etYz_l3, etYz_h1, etYz_h2, etYz_h3,  &
     thTz, thTz_l1, thTz_l2, thTz_l3, thTz_h1, thTz_h2, thTz_h3,  &
     thYz, thYz_l1, thYz_l2, thYz_l3, thYz_h1, thYz_h2, thYz_h3,  &
      kpp,  kpp_l1,  kpp_l2,  kpp_l3,  kpp_h1,  kpp_h2,  kpp_h3,  &
       jg,   jg_l1,   jg_l2,   jg_l3,   jg_h1,   jg_h2,   jg_h3,  &
     dt)

  use rad_params_module, only : ngroups, erg2rhoYe, clight
  use meth_params_module, only : NVAR, URHO

  implicit none

  integer,intent(in)::lo(3),hi(3)
  integer,intent(in)::Tp_n_l1, Tp_n_h1, Tp_n_l2, Tp_n_h2, Tp_n_l3, Tp_n_h3
  integer,intent(in)::Ye_n_l1, Ye_n_h1, Ye_n_l2, Ye_n_h2, Ye_n_l3, Ye_n_h3
  integer,intent(in)::Er_n_l1, Er_n_h1, Er_n_l2, Er_n_h2, Er_n_l3, Er_n_h3
  integer,intent(in)::re_s_l1, re_s_h1, re_s_l2, re_s_h2, re_s_l3, re_s_h3
  integer,intent(in)::rY_s_l1, rY_s_h1, rY_s_l2, rY_s_h2, rY_s_l3, rY_s_h3
  integer,intent(in)::re_2_l1, re_2_h1, re_2_l2, re_2_h2, re_2_l3, re_2_h3
  integer,intent(in)::rY_2_l1, rY_2_h1, rY_2_l2, rY_2_h2, rY_2_l3, rY_2_h3
  integer,intent(in)::etTz_l1, etTz_h1, etTz_l2, etTz_h2, etTz_l3, etTz_h3
  integer,intent(in)::etYz_l1, etYz_h1, etYz_l2, etYz_h2, etYz_l3, etYz_h3
  integer,intent(in)::thTz_l1, thTz_h1, thTz_l2, thTz_h2, thTz_l3, thTz_h3
  integer,intent(in)::thYz_l1, thYz_h1, thYz_l2, thYz_h2, thYz_l3, thYz_h3
  integer,intent(in):: kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2,  kpp_l3,  kpp_h3
  integer,intent(in)::  jg_l1,   jg_h1,   jg_l2,   jg_h2,   jg_l3,   jg_h3
  double precision,intent(inout)::Tp_n(Tp_n_l1:Tp_n_h1,Tp_n_l2:Tp_n_h2,Tp_n_l3:Tp_n_h3)
  double precision,intent(inout)::Ye_n(Ye_n_l1:Ye_n_h1,Ye_n_l2:Ye_n_h2,Ye_n_l3:Ye_n_h3)
  double precision,intent(in   )::Er_n(Er_n_l1:Er_n_h1,Er_n_l2:Er_n_h2,Er_n_l3:Er_n_h3,&
       0:ngroups-1)
  double precision,intent(in   )::re_s(re_s_l1:re_s_h1,re_s_l2:re_s_h2,re_s_l3:re_s_h3)
  double precision,intent(in   )::rY_s(rY_s_l1:rY_s_h1,rY_s_l2:rY_s_h2,rY_s_l3:rY_s_h3)
  double precision,intent(in   )::re_2(re_2_l1:re_2_h1,re_2_l2:re_2_h2,re_2_l3:re_2_h3)
  double precision,intent(in   )::rY_2(rY_2_l1:rY_2_h1,rY_2_l2:rY_2_h2,rY_2_l3:rY_2_h3)
  double precision,intent(in   )::etTz(etTz_l1:etTz_h1,etTz_l2:etTz_h2,etTz_l3:etTz_h3)
  double precision,intent(in   )::etYz(etYz_l1:etYz_h1,etYz_l2:etYz_h2,etYz_l3:etYz_h3)
  double precision,intent(in   )::thTz(thTz_l1:thTz_h1,thTz_l2:thTz_h2,thTz_l3:thTz_h3)
  double precision,intent(in   )::thYz(thYz_l1:thYz_h1,thYz_l2:thYz_h2,thYz_l3:thYz_h3)
  double precision,intent(in   ):: kpp( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,&
       0:ngroups-1)
  double precision,intent(in   )::  jg(  jg_l1:  jg_h1,  jg_l2:  jg_h2,  jg_l3:  jg_h3,&
       0:ngroups-1)
  double precision,intent(in) :: dt

   integer :: i,j,k,g
   double precision :: cdt1, cpT, cpY, foo, scrch_re, scrch_rY
   double precision :: dTemp, dYe
   double precision, parameter :: fac = 0.01d0

   cdt1 = 1.d0 / (clight * dt)
   !$omp parallel do private(i,j,k,g,cpT,cpY,foo,scrch_re,scrch_rY,dTemp,dYe)
   do k = lo(3), hi(3)
   do j = lo(2), hi(2)
   do i = lo(1), hi(1)

      cpT = 0.d0
      cpY = 0.d0
      do g = 0, ngroups-1
         foo = kpp(i,j,k,g)*Er_n(i,j,k,g) - jg(i,j,k,g)
         cpT = cpT + foo
         cpY = cpY + erg2rhoYe(g)*foo
      end do

      scrch_re = cpT - (re_s(i,j,k) - re_2(i,j,k)) * cdt1
      scrch_rY = cpY - (rY_s(i,j,k) - rY_2(i,j,k)) * cdt1

      dTemp = etTz(i,j,k)*scrch_re - thTz(i,j,k)*scrch_rY
      dYe   = -etYz(i,j,k)*scrch_re + thYz(i,j,k)*scrch_rY

      if (abs(dTemp/(Tp_n(i,j,k)+1.d-50)) > fac) then
         dTemp = sign(fac*Tp_n(i,j,k), dTemp)
      end if

      if (abs(dYe/(Ye_n(i,j,k)+1.d-50)) > fac) then
         dYe = sign(fac*Ye_n(i,j,k), dYe)
      end if

     Tp_n(i,j,k) = Tp_n(i,j,k) + dTemp
     Ye_n(i,j,k) = Ye_n(i,j,k) + dYe

  end do
  end do
  end do
  !$omp end parallel do

end subroutine ca_ncupdate_matter_neut


subroutine ca_compute_rosseland_neut( lo, hi, &
     kpr , kpr_l1, kpr_l2, kpr_l3, kpr_h1, kpr_h2, kpr_h3, &
     stat,stat_l1,stat_l2,stat_l3,stat_h1,stat_h2,stat_h3 )

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  kpr_l1,  kpr_h1,  kpr_l2,  kpr_h2,  kpr_l3,  kpr_h3
  integer, intent(in) :: stat_l1, stat_h1, stat_l2, stat_h2, stat_l3, stat_h3 
  double precision,intent(out)::kpr ( kpr_l1: kpr_h1, kpr_l2: kpr_h2, kpr_l3: kpr_h3,&
       0:ngroups-1)
  double precision,intent(in )::stat(stat_l1:stat_h1,stat_l2:stat_h2,stat_l3:stat_h3,NVAR)

  integer :: i, j, k, g, inu
  double precision :: ab, sc, delta, eta, er, der, rho, temp, Ye
  logical :: comp_ab, comp_sc, comp_eta

  comp_ab = .true.
  comp_sc = .true.
  comp_eta = .false.

  !$OMP PARALLEL DO SHARED(comp_ab,comp_sc,comp_eta) &
  !$OMP PRIVATE(g,inu,er,der,i,j,k,rho,temp,Ye,ab,sc,delta,eta)
  do g=0, ngroups-1

     call prep_opacity(g, inu, er, der)

     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rho = stat(i,j,k,URHO)
        temp = stat(i,j,k,UTEMP)
        Ye = stat(i,j,k,UFX) / rho

        call get_opacity_emissivity(ab, sc, delta, eta, &
             rho, Ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

        kpr(i,j,k,g) = ab + sc * (1.d0 - delta/3.d0)

     end do
     end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine ca_compute_rosseland_neut


subroutine ca_compute_planck_neut( lo, hi,  &
     kpp , kpp_l1, kpp_l2, kpp_l3, kpp_h1, kpp_h2, kpp_h3, &
     stat,stat_l1,stat_l2,stat_l3,stat_h1,stat_h2,stat_h3 )

  use rad_params_module, only : ngroups
  use opacity_table_module, only : prep_opacity, get_opacity_emissivity
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) ::  kpp_l1,  kpp_h1,  kpp_l2,  kpp_h2,  kpp_l3,  kpp_h3
  integer, intent(in) :: stat_l1, stat_h1, stat_l2, stat_h2, stat_l3, stat_h3 
  double precision, intent(out) :: kpp ( kpp_l1: kpp_h1, kpp_l2: kpp_h2, kpp_l3: kpp_h3,&
       0:ngroups-1)
  double precision, intent(in ) :: stat(stat_l1:stat_h1,stat_l2:stat_h2,stat_l3:stat_h3,NVAR)

  integer :: i, j, k, g, inu
  double precision :: ab, sc, delta, eta, er, der, rho, temp, Ye
  logical :: comp_ab, comp_sc, comp_eta
  
  comp_ab = .true.
  comp_sc = .false.
  comp_eta = .false.

  !$OMP PARALLEL DO SHARED(comp_ab,comp_sc,comp_eta) &
  !$OMP PRIVATE(g,inu,er,der,i,j,k,rho,temp,Ye,ab,sc,delta,eta)
  do g=0, ngroups-1

     call prep_opacity(g, inu, er, der)

     do k = lo(3), hi(3)
     do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        rho = stat(i,j,k,URHO)
        temp = stat(i,j,k,UTEMP)
        Ye = stat(i,j,k,UFX) / rho

        call get_opacity_emissivity(ab, sc, delta, eta, &
             rho, Ye, temp, er, inu, comp_ab, comp_sc, comp_eta)

        kpp(i,j,k,g) = ab 

     end do
     end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine ca_compute_planck_neut
