
subroutine ca_er_com2lab(lo, hi, &
     Snew,  S_l1,  S_l2,  S_h1,  S_h2, &
     Ecom, Ec_l1, Ec_l2, Ec_h1, Ec_h2, &
     F,     F_l1,  F_l2,  F_h1,  F_h2, iflx, nflx, & 
     Elab, El_l1, El_l2, El_h1, El_h2, ier, npv)
  use meth_params_module, only : NVAR, URHO, UMX, UMY
  use rad_params_module, only : ngroups, clight, ng0, ng1, dlognu
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) ::  S_l1,  S_l2,  S_h1,  S_h2
  integer, intent(in) :: Ec_l1, Ec_l2, Ec_h1, Ec_h2
  integer, intent(in) ::  F_l1,  F_l2,  F_h1,  F_h2, iflx, nflx
  integer, intent(in) :: El_l1, El_l2, El_h1, El_h2, ier, npv
  real(rt)        ,intent(in)   ::Snew( S_l1: S_h1, S_l2: S_h2,NVAR)
  real(rt)        ,intent(in)   ::Ecom(Ec_l1:Ec_h1,Ec_l2:Ec_h2,0:ngroups-1)
  real(rt)        ,intent(in)   ::F   ( F_l1: F_h1, F_l2: F_h2,0:nflx-1)
  real(rt)        ,intent(inout)::Elab(El_l1:El_h1,El_l2:El_h2,0:npv-1)

  integer :: i, j, g, ifx, ify
  real(rt)         :: rhoInv, c2, vxc2, vyc2
  real(rt)         :: nufnux(-1:ngroups), nufnuy(-1:ngroups)
  real(rt)         :: dlognuInv(0:ngroups-1)

  ifx = iflx
  ify = iflx + ngroups

  c2 = 1.e0_rt/clight**2

  if (ngroups > 1) dlognuInv = 1.e0_rt/dlognu

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rhoInv = 1.e0_rt/Snew(i,j,URHO)
     vxc2 = Snew(i,j,UMX)*rhoInv*c2
     vyc2 = Snew(i,j,UMY)*rhoInv*c2
     
     do g = 0, ngroups-1
        Elab(i,j,g+ier) = Ecom(i,j,g) + 2.e0_rt*(vxc2*F(i,j,ifx+g) &
             &                              + vyc2*F(i,j,ify+g))
     end do
     
     if (ngroups > 1) then
        do g=0,ngroups-1
           nufnux(g) = F(i,j,ifx+g)*dlognuInv(g)
           nufnuy(g) = F(i,j,ify+g)*dlognuInv(g)
        end do
        nufnux(-1) = -nufnux(0)
        nufnuy(-1) = -nufnuy(0)
        nufnux(ngroups) = -nufnux(ngroups-1)
        nufnuy(ngroups) = -nufnuy(ngroups-1)              
        do g=0,ngroups-1
           Elab(i,j,g+ier) = Elab(i,j,g+ier) - vxc2*0.5e0_rt*(nufnux(g+1)-nufnux(g-1)) &
                - vyc2*0.5e0_rt*(nufnuy(g+1)-nufnuy(g-1))
        end do
     end if
  end do
  end do

end subroutine ca_er_com2lab


subroutine ca_compute_fcc(lo, hi, &
     lamx, lamx_l1, lamx_l2, lamx_h1, lamx_h2, &
     lamy, lamy_l1, lamy_l2, lamy_h1, lamy_h2, nlam, &
     Eddf, Eddf_l1, Eddf_l2, Eddf_h1, Eddf_h2)
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: lamx_l1, lamx_l2, lamx_h1, lamx_h2
  integer, intent(in) :: lamy_l1, lamy_l2, lamy_h1, lamy_h2, nlam
  integer, intent(in) :: Eddf_l1, Eddf_l2, Eddf_h1, Eddf_h2
  real(rt)        ,intent(in   )::lamx(lamx_l1:lamx_h1,lamx_l2:lamx_h2,0:nlam-1)
  real(rt)        ,intent(in   )::lamy(lamy_l1:lamy_h1,lamy_l2:lamy_h2,0:nlam-1)
  real(rt)        ,intent(inout)::Eddf(Eddf_l1:Eddf_h1,Eddf_l2:Eddf_h2,0:ngroups-1)

  integer :: i, j, g, ilam
  real(rt)         :: lamcc

  do g=0,ngroups-1
     ilam = min(g,nlam-1)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           lamcc = 0.25e0_rt*(lamx(i,j,ilam)+lamx(i+1,j,ilam)+lamy(i,j,ilam)+lamy(i,j+1,ilam))
           Eddf(i,j,g) = Edd_factor(lamcc)
        end do
     end do
  end do
end subroutine ca_compute_fcc


subroutine ca_transform_flux (lo, hi, flag, &
     Snew,  S_l1,  S_l2,  S_h1,  S_h2, &
     f,     f_l1,  f_l2,  f_h1,  f_h2, &
     Er,   Er_l1, Er_l2, Er_h1, Er_h2, &
     Fi,   Fi_l1, Fi_l2, Fi_h1, Fi_h2, ifi, nfi, & 
     Fo,   Fo_l1, Fo_l2, Fo_h1, Fo_h2, ifo, nfo)
  use meth_params_module, only : NVAR, URHO, UMX, UMY
  use rad_params_module, only : ngroups, ng0, ng1, dlognu
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  real(rt)        , intent(in) :: flag
  integer, intent(in) ::  S_l1,  S_l2,  S_h1,  S_h2
  integer, intent(in) ::  f_l1,  f_l2,  f_h1,  f_h2
  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2
  integer, intent(in) :: Fi_l1, Fi_l2, Fi_h1, Fi_h2, ifi, nfi
  integer, intent(in) :: Fo_l1, Fo_l2, Fo_h1, Fo_h2, ifo, nfo
  real(rt)        ,intent(in   )::Snew( S_l1: S_h1, S_l2: S_h2,NVAR)
  real(rt)        ,intent(in   )::   f( f_l1: f_h1, f_l2: f_h2,0:ngroups-1)
  real(rt)        ,intent(in   )::  Er(Er_l1:Er_h1,Er_l2:Er_h2,0:ngroups-1)
  real(rt)        ,intent(in   )::  Fi(Fi_l1:Fi_h1,Fi_l2:Fi_h2,0:nfi-1)
  real(rt)        ,intent(inout)::  Fo(Fo_l1:Fo_h1,Fo_l2:Fo_h2,0:nfo-1)

  integer :: i, j, g, ifix, ifiy, ifox, ifoy
  real(rt)         :: rhoInv,  vx, vy, f1, f2, nx, ny, foo, vdotn
  real(rt)         :: nuvpnux(-1:ngroups), nuvpnuy(-1:ngroups)
  real(rt)         :: dlognuInv(0:ngroups-1)
  real(rt)         :: vdotpx(0:ngroups-1), vdotpy(0:ngroups-1)

  ifix = ifi
  ifiy = ifi + ngroups

  ifox = ifo
  ifoy = ifo + ngroups

  if (ngroups > 1) dlognuInv = 1.e0_rt/dlognu

  do j = lo(2), hi(2)
  do i = lo(1), hi(1)
     rhoInv = 1.e0_rt/Snew(i,j,URHO)
     vx = Snew(i,j,UMX)*rhoInv*flag
     vy = Snew(i,j,UMY)*rhoInv*flag

     do g = 0, ngroups-1
        f1 = (1.e0_rt-f(i,j,g))
        f2 = (3.e0_rt*f(i,j,g)-1.e0_rt)
        foo = 1.e0_rt/sqrt(Fi(i,j,ifix+g)**2+Fi(i,j,ifiy+g)**2+1.e-50_rt)
        nx = Fi(i,j,ifix+g)*foo
        ny = Fi(i,j,ifiy+g)*foo
        vdotn = vx*nx+vy*ny
        vdotpx(g) = 0.5e0_rt*Er(i,j,g)*(f1*vx + f2*vdotn*nx)
        vdotpy(g) = 0.5e0_rt*Er(i,j,g)*(f1*vy + f2*vdotn*ny)
        Fo(i,j,ifox+g) = Fi(i,j,ifix+g) + vx*Er(i,j,g) + vdotpx(g)
        Fo(i,j,ifoy+g) = Fi(i,j,ifiy+g) + vy*Er(i,j,g) + vdotpy(g)
     end do

     if (ngroups > 1) then
        do g=0,ngroups-1
           nuvpnux(g) = vdotpx(g)*dlognuInv(g)
           nuvpnuy(g) = vdotpy(g)*dlognuInv(g)
        end do
        nuvpnux(-1) = -nuvpnux(0)
        nuvpnuy(-1) = -nuvpnuy(0)
        nuvpnux(ngroups) = -nuvpnux(ngroups-1)
        nuvpnuy(ngroups) = -nuvpnuy(ngroups-1)              
        do g=0,ngroups-1
           Fo(i,j,ifox+g) = Fo(i,j,ifox+g) - 0.5e0_rt*(nuvpnux(g+1)-nuvpnux(g-1))
           Fo(i,j,ifoy+g) = Fo(i,j,ifoy+g) - 0.5e0_rt*(nuvpnuy(g+1)-nuvpnuy(g-1))
        end do
     end if
  end do
  end do

end subroutine ca_transform_flux
