
subroutine ca_er_com2lab(lo, hi, &
     Snew,  S_l1,  S_h1, &
     Ecom, Ec_l1, Ec_h1, &
     F,     F_l1,  F_h1, iflx, nflx, & 
     Elab, El_l1, El_h1, ier, npv)
  use meth_params_module, only : NVAR, URHO, UMX
  use rad_params_module, only : ngroups, clight, nnuspec, ng0, ng1, dlognu
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) ::  S_l1,  S_h1
  integer, intent(in) :: Ec_l1, Ec_h1
  integer, intent(in) ::  F_l1,  F_h1, iflx, nflx
  integer, intent(in) :: El_l1, El_h1, ier, npv
  double precision,intent(in)   ::Snew( S_l1: S_h1,NVAR)
  double precision,intent(in)   ::Ecom(Ec_l1:Ec_h1,0:ngroups-1)
  double precision,intent(in)   ::F   ( F_l1: F_h1,0:nflx-1)
  double precision,intent(inout)::Elab(El_l1:El_h1,0:npv-1)

  integer :: i, g, ifx
  double precision :: rhoInv, c2, vxc2
  double precision :: nufnux(-1:ngroups)
  double precision :: dlognuInv(0:ngroups-1)

  ifx = iflx

  c2 = 1.d0/clight**2

  if (ngroups > 1) dlognuInv = 1.d0/dlognu

  do i = lo(1), hi(1)
     rhoInv = 1.d0/Snew(i,URHO)
     vxc2 = Snew(i,UMX)*rhoInv*c2
     
     do g = 0, ngroups-1
        Elab(i,g+ier) = Ecom(i,g) + 2.d0*vxc2*F(i,ifx+g)
     end do
     
     if (ngroups > 1) then
        if (nnuspec .eq. 0) then
           
           do g=0,ngroups-1
              nufnux(g) = F(i,ifx+g)*dlognuInv(g)
           end do
           nufnux(-1) = -nufnux(0)
           nufnux(ngroups) = -nufnux(ngroups-1)
           do g=0,ngroups-1
              Elab(i,g+ier) = Elab(i,g+ier) - vxc2*0.5d0*(nufnux(g+1)-nufnux(g-1))
           end do
           
        else
           
           do g=0,ng0-1
              nufnux(g) = F(i,ifx+g)*dlognuInv(g)
           end do
           nufnux(-1) = -nufnux(0)
           nufnux(ng0) = -nufnux(ng0-1)
           do g=0,ng0-1
              Elab(i,g+ier) = Elab(i,g+ier) - vxc2*0.5d0*(nufnux(g+1)-nufnux(g-1))
           end do
           
           if (nnuspec >= 2) then
              
              do g=ng0,ng0+ng1-1
                 nufnux(g) = F(i,ifx+g)*dlognuInv(g)
              end do
              nufnux(ng0-1) = -nufnux(ng0)
              nufnux(ng0+ng1) = -nufnux(ng0+ng1-1)
              do g=ng0,ng0+ng1-1
                 Elab(i,g+ier) = Elab(i,g+ier) - vxc2*0.5d0*(nufnux(g+1)-nufnux(g-1))
              end do
              
           end if
           
           if (nnuspec == 3) then
              
              do g=ng0+ng1,ngroups-1
                 nufnux(g) = F(i,ifx+g)*dlognuInv(g)
              end do
              nufnux(ng0+ng1-1) = -nufnux(ng0+ng1)
              nufnux(ngroups) = -nufnux(ngroups-1)
              do g=ng0+ng1,ngroups-1
                 Elab(i,g+ier) = Elab(i,g+ier) - vxc2*0.5d0*(nufnux(g+1)-nufnux(g-1))
              end do
              
           end if

        end if
     end if
  end do

end subroutine ca_er_com2lab


subroutine ca_compute_fcc(lo, hi, &
     lamx, lamx_l1, lamx_h1, nlam, &
     Eddf, Eddf_l1, Eddf_h1)
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
  implicit none
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: lamx_l1, lamx_h1, nlam
  integer, intent(in) :: Eddf_l1, Eddf_h1
  double precision,intent(in   )::lamx(lamx_l1:lamx_h1,0:nlam-1)
  double precision,intent(inout)::Eddf(Eddf_l1:Eddf_h1,0:ngroups-1)

  integer :: i, g, ilam
  double precision :: lamcc

  do g=0,ngroups-1
     ilam = min(g,nlam-1)
     do i = lo(1), hi(1)
        lamcc = 0.5d0*(lamx(i,ilam)+lamx(i+1,ilam))
        Eddf(i,g) = Edd_factor(lamcc)
     end do
  end do
end subroutine ca_compute_fcc


subroutine ca_transform_flux (lo, hi, flag, &
     Snew,  S_l1,  S_h1, &
     f,     f_l1,  f_h1, &
     Er,   Er_l1, Er_h1, &
     Fi,   Fi_l1, Fi_h1, ifi, nfi, & 
     Fo,   Fo_l1, Fo_h1, ifo, nfo)
  use meth_params_module, only : NVAR, URHO, UMX
  use rad_params_module, only : ngroups, nnuspec, ng0, ng1, dlognu
  implicit none

  integer, intent(in) :: lo(1), hi(1)
  double precision, intent(in) :: flag
  integer, intent(in) ::  S_l1,  S_h1
  integer, intent(in) ::  f_l1,  f_h1
  integer, intent(in) :: Er_l1, Er_h1
  integer, intent(in) :: Fi_l1, Fi_h1, ifi, nfi
  integer, intent(in) :: Fo_l1, Fo_h1, ifo, nfo
  double precision,intent(in   )::Snew( S_l1: S_h1,NVAR)
  double precision,intent(in   )::   f( f_l1: f_h1,0:ngroups-1)
  double precision,intent(in   )::  Er(Er_l1:Er_h1,0:ngroups-1)
  double precision,intent(in   )::  Fi(Fi_l1:Fi_h1,0:nfi-1)
  double precision,intent(inout)::  Fo(Fo_l1:Fo_h1,0:nfo-1)

  integer :: i, g, ifix, ifox
  double precision :: rhoInv,  vx, f1, f2, nx, foo, vdotn
  double precision :: nuvpnux(-1:ngroups)
  double precision :: dlognuInv(0:ngroups-1)
  double precision :: vdotpx(0:ngroups-1)

  ifix = ifi
  ifox = ifo

  if (ngroups > 1) dlognuInv = 1.d0/dlognu

  do i = lo(1), hi(1)
     rhoInv = 1.d0/Snew(i,URHO)
     vx = Snew(i,UMX)*rhoInv*flag

     do g = 0, ngroups-1
        f1 = (1.d0-f(i,g))
        f2 = (3.d0*f(i,g)-1.d0)
        foo = 1.d0/abs(Fi(i,ifix+g)+1.d-50)
        nx = Fi(i,ifix+g)*foo
        vdotn = vx*nx
        vdotpx(g) = 0.5d0*Er(i,g)*(f1*vx + f2*vdotn*nx)
        Fo(i,ifox+g) = Fi(i,ifix+g) + vx*Er(i,g) + vdotpx(g)
     end do

     if (ngroups > 1) then
        if (nnuspec .eq. 0) then
           
           do g=0,ngroups-1
              nuvpnux(g) = vdotpx(g)*dlognuInv(g)
           end do
           nuvpnux(-1) = -nuvpnux(0)
           nuvpnux(ngroups) = -nuvpnux(ngroups-1)
           do g=0,ngroups-1
              Fo(i,ifox+g) = Fo(i,ifox+g) - 0.5d0*(nuvpnux(g+1)-nuvpnux(g-1))
           end do

        else
           
           do g=0,ng0-1
              nuvpnux(g) = vdotpx(g)*dlognuInv(g)
           end do
           nuvpnux(-1) = -nuvpnux(0)
           nuvpnux(ng0) = -nuvpnux(ng0-1)
           do g=0,ng0-1
              Fo(i,ifox+g) = Fo(i,ifox+g) - 0.5d0*(nuvpnux(g+1)-nuvpnux(g-1))
           end do

           if (nnuspec >= 2) then

              do g=ng0,ng0+ng1-1
                 nuvpnux(g) = vdotpx(g)*dlognuInv(g)
              end do
              nuvpnux(ng0-1) = -nuvpnux(ng0)
              nuvpnux(ng0+ng1) = -nuvpnux(ng0+ng1-1)
              do g=ng0,ng0+ng1-1
                 Fo(i,ifox+g) = Fo(i,ifox+g) - 0.5d0*(nuvpnux(g+1)-nuvpnux(g-1))
              end do

           end if

           if (nnuspec == 3) then

              do g=ng0+ng1,ngroups-1
                 nuvpnux(g) = vdotpx(g)*dlognuInv(g)
              end do
              nuvpnux(ng0+ng1-1) = -nuvpnux(ng0+ng1)
              nuvpnux(ngroups) = -nuvpnux(ngroups-1)
              do g=ng0+ng1,ngroups-1
                 Fo(i,ifox+g) = Fo(i,ifox+g) - 0.5d0*(nuvpnux(g+1)-nuvpnux(g-1))
              end do

           end if

        end if
     end if
  end do

end subroutine ca_transform_flux
