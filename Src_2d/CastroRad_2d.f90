subroutine ca_umdrv_rad(is_finest_level,time,&
                        lo,hi,domlo,domhi,&
                        uin     ,     uin_l1,     uin_l2,     uin_h1,     uin_h2,&
                        uout    ,    uout_l1,    uout_l2,    uout_h1,    uout_h2,&
                        Erin    ,    Erin_l1,    Erin_l2,    Erin_h1,    Erin_h2,&
                        lam     ,     lam_l1,     lam_l2,     lam_h1,     lam_h2,&
                        Erout   ,   Erout_l1,   Erout_l2,   Erout_h1,   Erout_h2,&
                        ugdx    ,    ugdx_l1,    ugdx_l2,    ugdx_h1,    ugdx_h2, &
                        ugdy    ,    ugdy_l1,    ugdy_l2,    ugdy_h1,    ugdy_h2, &
                        src     ,     src_l1,     src_l2,     src_h1,     src_h2, &
                        grav    ,      gv_l1,      gv_l2,      gv_h1,      gv_h2, &
                        delta,dt,&
                        flux1   ,   flux1_l1,   flux1_l2,   flux1_h1,   flux1_h2, &
                        flux2   ,   flux2_l1,   flux2_l2,   flux2_h1,   flux2_h2, &
                        radflux1,radflux1_l1,radflux1_l2,radflux1_h1,radflux1_h2, &
                        radflux2,radflux2_l1,radflux2_l2,radflux2_h1,radflux2_h2, &
                        area1   ,   area1_l1,   area1_l2,   area1_h1,   area1_h2, &
                        area2   ,   area2_l1,   area2_l2,   area2_h1,   area2_h2, &
                        dloga   ,   dloga_l1,   dloga_l2,   dloga_h1,   dloga_h2, &
                        vol     ,     vol_l1,     vol_l2,     vol_h1,     vol_h2, &
                        courno,verbose, nstep_fsp)

  use meth_params_module, only : QVAR, NVAR, NHYP, do_sponge, normalize_species
  use rad_params_module, only : ngroups
  use radhydro_params_module, only : QRADVAR
  use advection_module, only : enforce_minimum_density, normalize_new_species, divu
  use rad_advection_module, only : umeth2d_rad, ctoprim_rad, consup_rad
  use sponge_module, only : sponge

  implicit none

  integer nstep_fsp
  integer is_finest_level
  integer lo(2),hi(2),verbose
  integer domlo(2),domhi(2)
  integer      uin_l1,     uin_l2,     uin_h1,     uin_h2
  integer     Erin_l1,    Erin_l2,    Erin_h1,    Erin_h2
  integer      lam_l1,     lam_l2,     lam_h1,     lam_h2
  integer     uout_l1,    uout_l2,    uout_h1,    uout_h2
  integer    Erout_l1,   Erout_l2,   Erout_h1,   Erout_h2
  integer     ugdx_l1,    ugdx_l2,    ugdx_h1,    ugdx_h2
  integer     ugdy_l1,    ugdy_l2,    ugdy_h1,    ugdy_h2
  integer    flux1_l1,   flux1_l2,   flux1_h1,   flux1_h2
  integer    flux2_l1,   flux2_l2,   flux2_h1,   flux2_h2
  integer radflux1_l1,radflux1_l2,radflux1_h1,radflux1_h2
  integer radflux2_l1,radflux2_l2,radflux2_h1,radflux2_h2
  integer    area1_l1,   area1_l2,   area1_h1,   area1_h2
  integer    area2_l1,   area2_l2,   area2_h1,   area2_h2
  integer    dloga_l1,   dloga_l2,   dloga_h1,   dloga_h2
  integer      vol_l1,     vol_l2,     vol_h1,     vol_h2
  integer      src_l1,     src_l2,     src_h1,     src_h2
  integer       gv_l1,      gv_l2,      gv_h1,      gv_h2

  double precision uin     (     uin_l1:     uin_h1,     uin_l2:     uin_h2,NVAR)
  double precision uout    (    uout_l1:    uout_h1,    uout_l2:    uout_h2,NVAR)
  double precision Erin    (    Erin_l1:    Erin_h1,    Erin_l2:    Erin_h2,0:ngroups-1)
  double precision lam     (     lam_l1:     lam_h1,     lam_l2:     lam_h2,0:ngroups-1)
  double precision Erout   (   Erout_l1:   Erout_h1,   Erout_l2:   Erout_h2,0:ngroups-1)
  double precision ugdx    (    ugdx_l1:    ugdx_h1,    ugdx_l2:    ugdx_h2)
  double precision ugdy    (    ugdy_l1:    ugdy_h1,    ugdy_l2:    ugdy_h2)
  double precision src     (     src_l1:     src_h1,     src_l2:     src_h2,NVAR)
  double precision grav    (      gv_l1:      gv_h1,      gv_l2:      gv_h2,2)
  double precision flux1   (   flux1_l1:   flux1_h1,   flux1_l2:   flux1_h2,NVAR)
  double precision flux2   (   flux2_l1:   flux2_h1,   flux2_l2:   flux2_h2,NVAR)
  double precision radflux1(radflux1_l1:radflux1_h1,radflux1_l2:radflux1_h2,0:ngroups-1)
  double precision radflux2(radflux2_l1:radflux2_h1,radflux2_l2:radflux2_h2,0:ngroups-1)
  double precision area1   (   area1_l1:   area1_h1,   area1_l2:   area1_h2)
  double precision area2   (   area2_l1:   area2_h1,   area2_l2:   area2_h2)
  double precision dloga   (   dloga_l1:   dloga_h1,   dloga_l2:   dloga_h2)
  double precision vol     (     vol_l1:     vol_h1,     vol_l2:     vol_h2)
  double precision delta(2),dt,time,courno
  
  !     Automatic arrays for workspace
  double precision, allocatable:: q(:,:,:)
  double precision, allocatable:: gamc(:,:)
  double precision, allocatable:: gamcg(:,:)
  double precision, allocatable:: flatn(:,:)
  double precision, allocatable:: c(:,:)
  double precision, allocatable:: cg(:,:)
  double precision, allocatable:: csml(:,:)
  double precision, allocatable:: div(:,:)
  double precision, allocatable:: pgdx(:,:)
  double precision, allocatable:: pgdy(:,:)
  double precision, allocatable:: ergdx(:,:,:)
  double precision, allocatable:: ergdy(:,:,:)
  double precision, allocatable:: lamgdx(:,:,:)
  double precision, allocatable:: lamgdy(:,:,:)
  double precision, allocatable:: srcQ(:,:,:)
  double precision, allocatable:: pdivu(:,:)

  double precision, allocatable :: uy_xfc(:,:), ux_yfc(:,:)

  integer ngq,ngf,iflaten
  integer q_l1, q_l2, q_h1, q_h2
  double precision dx,dy,mass_added,eint_added,eden_added
  double precision E_added_sponge,xmom_added_sponge,ymom_added_sponge

  ngq = NHYP
  ngf = 1
  iflaten = 1

  q_l1 = lo(1)-NHYP
  q_l2 = lo(2)-NHYP
  q_h1 = hi(1)+NHYP
  q_h2 = hi(2)+NHYP
  
  allocate(     q(q_l1:q_h1,q_l2:q_h2,QRADVAR))
  allocate(  gamc(q_l1:q_h1,q_l2:q_h2))
  allocate( gamcg(q_l1:q_h1,q_l2:q_h2))
  allocate( flatn(q_l1:q_h1,q_l2:q_h2))
  allocate(     c(q_l1:q_h1,q_l2:q_h2))
  allocate(    cg(q_l1:q_h1,q_l2:q_h2))
  allocate(  csml(q_l1:q_h1,q_l2:q_h2))

  allocate(  srcQ(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,QVAR))

  allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))
  allocate( pdivu(lo(1)  :hi(1)  ,lo(2)  :hi(2)))
  allocate(  pgdx(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1))
  allocate(  pgdy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1))
  allocate( ergdx(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1,0:ngroups-1))
  allocate( ergdy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1,0:ngroups-1))
  allocate(lamgdx(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1,0:ngroups-1))
  allocate(lamgdy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1,0:ngroups-1))

  allocate(uy_xfc(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2))
  allocate(ux_yfc(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2))
 
  dx = delta(1)
  dy = delta(2)
  
  !     Translate to primitive variables, compute sound speeds
  !     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
  !       and set to correspond to coordinates of (lo:hi)
  call ctoprim_rad(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2, &
       Erin,Erin_l1,Erin_l2,Erin_h1,Erin_h2, &
       lam,lam_l1,lam_l2,lam_h1,lam_h2, &
       q,c,cg,gamc,gamcg,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
       src,src_l1,src_l2,src_h1,src_h2, &
       srcQ,lo(1)-1,lo(2)-1,hi(1)+1,hi(2)+1, &
       courno,dx,dy,dt,ngq,ngf,iflaten)

!     Compute hyperbolic fluxes using unsplit Godunov
  call umeth2d_rad(q,c,cg,gamc,gamcg,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
       lam, lam_l1,lam_l2,lam_h1,lam_h2, &
       srcQ,lo(1)-1,lo(2)-1,hi(1)+1,hi(2)+1, &
       grav,gv_l1,gv_l2,gv_h1,gv_h2, &
       lo(1),lo(2),hi(1),hi(2),dx,dy,dt, &
       flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
       flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
       radflux1,radflux1_l1,radflux1_l2,radflux1_h1,radflux1_h2, &
       radflux2,radflux2_l1,radflux2_l2,radflux2_h1,radflux2_h2, &
       pgdx,  lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
       pgdy,  lo(1)-1, lo(2), hi(1)+1, hi(2)+1, &
       ergdx, lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
       ergdy, lo(1)-1, lo(2), hi(1)+1, hi(2)+1, &
       lamgdx, lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
       lamgdy, lo(1)-1, lo(2), hi(1)+1, hi(2)+1, &
       ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
       ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
       area1, area1_l1, area1_l2, area1_h1, area1_h2, &
       area2, area2_l1, area2_l2, area2_h1, area2_h2, &
       pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
       uy_xfc, ux_yfc, &
       dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2)

  !     Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2, &
       delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

  !     Conservative update
  call consup_rad(uin,  uin_l1,  uin_l2,  uin_h1,  uin_h2, &
       uout,  uout_l1, uout_l2, uout_h1, uout_h2, &
       Erin,Erin_l1,Erin_l2,Erin_h1,Erin_h2, &
       Erout,Erout_l1,Erout_l2,Erout_h1,Erout_h2, &
       pgdx,   lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
       pgdy, lo(1)-1,   lo(2), hi(1)+1, hi(2)+1, &
       ergdx,  lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
       ergdy,lo(1)-1,   lo(2), hi(1)+1, hi(2)+1, &
       lamgdx,  lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
       lamgdy,lo(1)-1,   lo(2), hi(1)+1, hi(2)+1, &
       ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
       ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
       src,    src_l1,  src_l2,  src_h1,  src_h2, &
       grav,    gv_l1,   gv_l2,   gv_h1,   gv_h2, &
       flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
       flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
       radflux1,radflux1_l1,radflux1_l2,radflux1_h1,radflux1_h2, &
       radflux2,radflux2_l1,radflux2_l2,radflux2_h1,radflux2_h2, &
       area1,area1_l1,area1_l2,area1_h1,area1_h2, &
       area2,area2_l1,area2_l2,area2_h1,area2_h2, &
       vol,    vol_l1,  vol_l2,  vol_h1,  vol_h2, &
       div,pdivu, uy_xfc, ux_yfc, &
       lo,hi,dx,dy,dt, nstep_fsp)

  ! Enforce the density >= small_dens.
  mass_added = 0.d0
  eint_added = 0.d0
  eden_added = 0.d0
  call enforce_minimum_density( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
       uout,uout_l1,uout_l2,uout_h1,uout_h2,&
       lo,hi,mass_added,eint_added,eden_added,verbose)
  
  ! Enforce the species >= 0
  call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)
  
  ! Normalize the species 
  if (normalize_species .eq. 1) &
       call normalize_new_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)
  
  if (do_sponge .eq. 1) then
     E_added_sponge = 0.d0
     xmom_added_sponge = 0.d0
     ymom_added_sponge = 0.d0
     call sponge(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi, &
          time,dt,dx,dy,domlo,domhi, &
          E_added_sponge,xmom_added_sponge,ymom_added_sponge)
  end if
  
  deallocate(q,gamc,gamcg,flatn,c,cg,csml,div,pgdx,pgdy,ergdx,ergdy)
  deallocate(lamgdx,lamgdy,srcQ,pdivu,uy_xfc, ux_yfc)

end subroutine ca_umdrv_rad

! This subroutine cannot be tiled
subroutine ca_compute_lamborder(Er, Er_l1, Er_l2, Er_h1, Er_h2, &
     kap, kap_l1, kap_l2, kap_h1, kap_h2, &
     lam, lam_l1, lam_l2, lam_h1, lam_h2, &
     dx, ngrow, limiter, filter_T, S)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda
  use filter_module

  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2, kap_l1, kap_l2, kap_h1, kap_h2, &
       lam_l1, lam_l2, lam_h1, lam_h2
  integer, intent(in) :: ngrow, limiter, filter_T, S
  double precision, intent(in) :: dx(2)
  double precision, intent(in) :: kap(kap_l1:kap_h1, kap_l2:kap_h2)
  double precision, intent(in) :: Er(Er_l1:Er_h1, Er_l2:Er_h2, 0:ngroups-1)
  double precision, intent(out) :: lam(lam_l1:lam_h1, lam_l2:lam_h2, 0:ngroups-1)

  integer :: i, j, reg_l1, reg_l2, reg_h1, reg_h2, g
  double precision :: r, r1, r2

  double precision, allocatable :: lamfil(:,:)

  lam = -1.d50

  reg_l1 = lam_l1 + ngrow
  reg_l2 = lam_l2 + ngrow
  reg_h1 = lam_h1 - ngrow
  reg_h2 = lam_h2 - ngrow

  if (filter_T .gt. 0) then
     allocate(lamfil(reg_l1:reg_h1,lam_l2:lam_h2))
  end if

  do g = 0, ngroups-1

  do j=lam_l2, lam_h2
     do i=lam_l1, lam_h1
        if (Er(i,j,g) .eq. -1.d0) then
           cycle
        end if

        if (Er(i-1,j,g) .eq. -1.d0) then
           r1 = (Er(i+1,j,g) - Er(i,j,g)) / (dx(1))
        else if (Er(i+1,j,g) .eq. -1.d0) then
           r1 = (Er(i,j,g) - Er(i-1,j,g)) / (dx(1))
        else
           r1 = (Er(i+1,j,g) - Er(i-1,j,g)) / (2.d0*dx(1))
        end if

        if (Er(i,j-1,g) .eq. -1.d0) then
           r2 = (Er(i,j+1,g) - Er(i,j,g)) / (dx(2))
        else if (Er(i,j+1,g) .eq. -1.d0) then
           r2 = (Er(i,j,g) - Er(i,j-1,g)) / (dx(2))
        else
           r2 = (Er(i,j+1,g) - Er(i,j-1,g)) / (2.d0*dx(2))
        end if

        r = sqrt(r1**2 + r2**2)
        r = r / (kap(i,j) * max(Er(i,j,g), 1.d-50))
        
        lam(i,j,g) = FLDlambda(r, limiter)
     end do
  end do

  ! filter
  if (filter_T .eq. 1) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.d0) then
           lamfil(:,j) = -1.d-50
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff1(0) * lam(i,j,g) &
                &      + ff1(1) * (lam(i-1,j,g)+lam(i+1,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.d0) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff1b, lam(i:i+1,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.d0) then
            i = reg_h1
            lamfil(i,j) = dot_product(ff1b(1:0:-1), lam(i-1:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff1(0) * lamfil(i,j) &
                &     + ff1(1) * (lamfil(i,j-1)+lamfil(i,j+1))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.d0) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff1b, lamfil(i,j:j+1))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.d0) then
        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff1b(1:0:-1), lamfil(i,j-1:j))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

  else if (filter_T .eq. 2) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.d0) then
           lamfil(:,j) = -1.d-50
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff2(0,S) * lam(i,j,g) &
                &      + ff2(1,S) * (lam(i-1,j,g)+lam(i+1,j,g)) &
                &      + ff2(2,S) * (lam(i-2,j,g)+lam(i+2,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.d0) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff2b0, lam(i:i+2,j,g))

            i = reg_l1 + 1
            lamfil(i,j) = dot_product(ff2b1, lam(i-1:i+2,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.d0) then
            i = reg_h1 - 1
            lamfil(i,j) = dot_product(ff2b1(2:-1:-1), lam(i-2:i+1,j,g))

            i = reg_h1
            lamfil(i,j) = dot_product(ff2b0(2:0:-1), lam(i-2:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff2(0,S) * lamfil(i,j) &
                &     + ff2(1,S) * (lamfil(i,j-1)+lamfil(i,j+1)) &
                &     + ff2(2,S) * (lamfil(i,j-2)+lamfil(i,j+2))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.d0) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b0, lamfil(i,j:j+2))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_l2 + 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b1, lamfil(i,j-1:j+2))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.d0) then
        j = reg_h2 - 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b1(2:-1:-1), lamfil(i,j-2:j+1))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff2b0(2:0:-1), lamfil(i,j-2:j))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

  else if (filter_T .eq. 3) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.d0) then
           lamfil(:,j) = -1.d-50
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff3(0,S) * lam(i,j,g) &
                &      + ff3(1,S) * (lam(i-1,j,g)+lam(i+1,j,g)) &
                &      + ff3(2,S) * (lam(i-2,j,g)+lam(i+2,j,g)) &
                &      + ff3(3,S) * (lam(i-3,j,g)+lam(i+3,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.d0) then
            i = reg_l1
            lamfil(i,j) = dot_product(ff3b0, lam(i:i+3,j,g))

            i = reg_l1 + 1
            lamfil(i,j) = dot_product(ff3b1, lam(i-1:i+3,j,g))

            i = reg_l1 + 2
            lamfil(i,j) = dot_product(ff3b2, lam(i-2:i+3,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.d0) then
            i = reg_h1 - 2
            lamfil(i,j) = dot_product(ff3b2(3:-2:-1), lam(i-3:i+2,j,g))

            i = reg_h1 - 1
            lamfil(i,j) = dot_product(ff3b1(3:-1:-1), lam(i-3:i+1,j,g))

            i = reg_h1
            lamfil(i,j) = dot_product(ff3b0(3:0:-1), lam(i-3:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff3(0,S) * lamfil(i,j) &
                &     + ff3(1,S) * (lamfil(i,j-1)+lamfil(i,j+1)) &
                &     + ff3(2,S) * (lamfil(i,j-2)+lamfil(i,j+2)) &
                &     + ff3(3,S) * (lamfil(i,j-3)+lamfil(i,j+3))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.d0) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b0, lamfil(i,j:j+3))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_l2 + 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b1, lamfil(i,j-1:j+3))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_l2 + 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b2, lamfil(i,j-2:j+3))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.d0) then
        j = reg_h2 - 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b2(3:-2:-1), lamfil(i,j-3:j+2))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_h2 - 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b1(3:-1:-1), lamfil(i,j-3:j+1))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff3b0(3:0:-1), lamfil(i,j-3:j))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

  else if (filter_T .eq. 4) then

     do j=lam_l2, lam_h2
        if (Er(reg_l1,j,g) .eq. -1.d0) then
           lamfil(:,j) = -1.d-50
           cycle
        endif

        do i=reg_l1, reg_h1
           lamfil(i,j) = ff4(0,S) * lam(i,j,g) &
                &      + ff4(1,S) * (lam(i-1,j,g)+lam(i+1,j,g)) &
                &      + ff4(2,S) * (lam(i-2,j,g)+lam(i+2,j,g)) &
                &      + ff4(3,S) * (lam(i-3,j,g)+lam(i+3,j,g)) &
                &      + ff4(4,S) * (lam(i-4,j,g)+lam(i+4,j,g))
        end do

        if (Er(reg_l1-1,j,g) .eq. -1.d0) then
            i = reg_l1 
            lamfil(i,j) = dot_product(ff4b0, lam(i:i+4,j,g))

            i = reg_l1 + 1
            lamfil(i,j) = dot_product(ff4b1, lam(i-1:i+4,j,g))

            i = reg_l1 + 2
            lamfil(i,j) = dot_product(ff4b2, lam(i-2:i+4,j,g))

            i = reg_l1 + 3
            lamfil(i,j) = dot_product(ff4b3, lam(i-3:i+4,j,g))
         end if

         if (Er(reg_h1+1,j,g) .eq. -1.d0) then
            i = reg_h1 - 3
            lamfil(i,j) = dot_product(ff4b3(4:-3:-1), lam(i-4:i+3,j,g))

            i = reg_h1 - 2
            lamfil(i,j) = dot_product(ff4b2(4:-2:-1), lam(i-4:i+2,j,g))

            i = reg_h1 - 1
            lamfil(i,j) = dot_product(ff4b1(4:-1:-1), lam(i-4:i+1,j,g))

            i = reg_h1
            lamfil(i,j) = dot_product(ff4b0(4:0:-1), lam(i-4:i,j,g))
         end if
     end do

     do j=reg_l2, reg_h2
        do i=reg_l1, reg_h1
           lam(i,j,g) = ff4(0,S) * lamfil(i,j) &
                &     + ff4(1,S) * (lamfil(i,j-1)+lamfil(i,j+1)) &
                &     + ff4(2,S) * (lamfil(i,j-2)+lamfil(i,j+2)) &
                &     + ff4(3,S) * (lamfil(i,j-3)+lamfil(i,j+3)) &
                &     + ff4(4,S) * (lamfil(i,j-4)+lamfil(i,j+4))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end do

     if (Er(reg_l1,reg_l2-1,g) .eq. -1.d0) then
        j = reg_l2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b0, lamfil(i,j:j+4))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_l2 + 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b1, lamfil(i,j-1:j+4))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_l2 + 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b2, lamfil(i,j-2:j+4))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_l2 + 3
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b3, lamfil(i,j-3:j+4))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

     if (Er(reg_l1,reg_h2+1,g) .eq. -1.d0) then
        j = reg_h2 - 3
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b3(4:-3:-1), lamfil(i,j-4:j+3))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_h2 - 2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b2(4:-2:-1), lamfil(i,j-4:j+2))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_h2 - 1
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b1(4:-1:-1), lamfil(i,j-4:j+1))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do

        j = reg_h2
        do i=reg_l1,reg_h1
           lam(i,j,g) = dot_product(ff4b0(4:0:-1), lamfil(i,j-4:j))
           lam(i,j,g) = min(1.d0/3.d0, max(1.d-25, lam(i,j,g)))
        end do
     end if

  end if

  ! lo-x lo-y 
  do j=lam_l2,reg_l2-1
     do i=lam_l1,reg_l1-1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(reg_l1,reg_l2,g)
        end if
     end do
  end do

  ! reg-x lo-y 
  do j=lam_l2,reg_l2-1
     do i=reg_l1,reg_h1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(i,reg_l2,g)
        end if
     end do
  end do
  
  ! hi-x lo-y 
  do j=lam_l2,reg_l2-1
     do i=reg_h1+1,lam_h1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(reg_h1,reg_l2,g)
        end if
     end do
  end do

  ! lo-x reg-y
  do j=reg_l2,reg_h2
     do i=lam_l1,reg_l1-1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(reg_l1,j,g)
        end if
     end do
  end do

  ! hi-x reg-y
  do j=reg_l2,reg_h2
     do i=reg_h1+1,lam_h1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(reg_h1,j,g)
        end if
     end do
  end do

  ! lo-x hi-y 
  do j=reg_h2+1,lam_h2
     do i=lam_l1,reg_l1-1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(reg_l1,reg_h2,g)
        end if
     end do
  end do

  ! reg-x hi-y 
  do j=reg_h2+1,lam_h2
     do i=reg_l1,reg_h1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(i,reg_h2,g)
        end if
     end do
  end do

  ! hi-x hi-y
  do j=reg_h2+1,lam_h2
     do i=reg_h1+1,lam_h1
        if (Er(i,j,g).eq.-1.d0) then
           lam(i,j,g) = lam(reg_h1,reg_h2,g)
        end if
     end do
  end do

  end do

  if (filter_T .gt. 0) then
     deallocate(lamfil)
  end if

  return
end subroutine ca_compute_lamborder


subroutine ca_get_v_dcf( lo, hi, &
     er ,  er_l1,  er_l2,  er_h1,  er_h2, &
     s  ,   s_l1,   s_l2,   s_h1,   s_h2, &
     T  ,   T_l1,   T_l2,   T_h1,   T_h2, &
     c_v, c_v_l1, c_v_l2, c_v_h1, c_v_h2, &
     kr ,  kr_l1,  kr_l2,  kr_h1,  kr_h2, &
     kp ,  kp_l1,  kp_l2,  kp_h1,  kp_h2, &
     kp2, kp2_l1, kp2_l2, kp2_h1, kp2_h2, &
     dtemp, dtime, sigma, c, &
     v  ,   v_l1,   v_l2,   v_h1,   v_h2, &
     dcf, dcf_l1, dcf_l2, dcf_h1, dcf_h2)

  use meth_params_module, only : NVAR, URHO, UMX, UMY

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: er_l1,er_l2,er_h1,er_h2,s_l1,s_l2,s_h1,s_h2
  integer, intent(in) :: T_l1,T_l2,T_h1,T_h2,c_v_l1,c_v_l2,c_v_h1,c_v_h2
  integer, intent(in) :: kr_l1,kr_l2,kr_h1,kr_h2
  integer, intent(in) :: kp_l1,kp_l2,kp_h1,kp_h2,kp2_l1,kp2_l2,kp2_h1,kp2_h2
  integer, intent(in) :: v_l1,v_l2,v_h1,v_h2,dcf_l1,dcf_l2,dcf_h1,dcf_h2
  double precision, intent(in)  ::  er( er_l1: er_h1,  er_l2: er_h2)
  double precision, intent(in)  ::   s(  s_l1:  s_h1,   s_l2:  s_h2, NVAR)
  double precision, intent(in)  ::   T(  T_l1:  T_h1,   T_l2:  T_h2)
  double precision, intent(in)  :: c_v(c_v_l1:c_v_h1, c_v_l2:c_v_h2)
  double precision, intent(in ) ::  kr( kr_l1: kr_h1,  kr_l2: kr_h2)
  double precision, intent(in ) ::  kp( kp_l1: kp_h1,  kp_l2: kp_h2)
  double precision, intent(in ) :: kp2(kp2_l1:kp2_h1, kp2_l2:kp2_h2)
  double precision, intent(in) :: dtemp, dtime, sigma, c
  double precision              ::   v(  v_l1:  v_h1,   v_l2:  v_h2, 2)
  double precision              :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2)

  integer :: i, j
  double precision :: etainv, fac0, fac2, alpha, frc

  fac0 = 4.d0 * sigma * dtime / dtemp
  fac2 = c * dtime / dtemp

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        v(i,j,1) = s(i,j,UMX)/s(i,j,URHO)
        v(i,j,2) = s(i,j,UMY)/s(i,j,URHO)
        
        alpha = fac0 * (kp2(i,j) * (T(i,j) + dtemp) ** 4    &
             -          kp (i,j) * (T(i,j)        ) ** 4)   &
             -  fac2 * (kp2(i,j) - kp(i,j)) * er(i,j)
        
        frc = s(i,j,URHO) * c_v(i,j) + 1.0d-50
        etainv = frc / (alpha + frc)
        
        dcf(i,j) = 2.d0 * etainv * (kp(i,j) / kr(i,j))
     end do
  end do

end subroutine ca_get_v_dcf


subroutine ca_compute_dcoefs( lo, hi, &
     d  ,   d_l1,   d_l2,   d_h1,   d_h2, &
     lam, lam_l1, lam_l2, lam_h1, lam_h2, &
     v ,    v_l1,   v_l2,   v_h1,   v_h2, &
     dcf, dcf_l1, dcf_l2, dcf_h1, dcf_h2, &
     r, idir)

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: d_l1, d_l2, d_h1, d_h2, &
       & lam_l1, lam_l2, lam_h1, lam_h2, &
       &   v_l1,   v_l2,   v_h1,   v_h2, &
       & dcf_l1, dcf_l2, dcf_h1, dcf_h2, &
       idir

  double precision              ::   d(  d_l1:  d_h1,   d_l2:  d_h2)
  double precision, intent(in)  :: lam(lam_l1:lam_h1, lam_l2:lam_h2)
  double precision, intent(in)  ::   v(  v_l1:  v_h1,   v_l2:  v_h2, 2)
  double precision, intent(in)  :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2)
  double precision, intent(in)  ::   r( lo(1): hi(1))

  integer :: i, j

  if (idir.eq.0) then
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (v(i-1,j,1) + v(i,j,1) .gt. 0.d0) then
              d(i,j) = dcf(i-1,j) * v(i-1,j,1) * lam(i,j)
           else if (v(i-1,j,1) + v(i,j,1) .lt. 0.d0) then
              d(i,j) = dcf(i,j) * v(i,j,1) * lam(i,j)
           else
              d(i,j) = 0.0
           end if
           d(i,j) = d(i,j) * r(i)
        end do
     end do
  else
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (v(i,j-1,2) + v(i,j,2) .gt. 0.d0) then
              d(i,j) = dcf(i,j-1) * v(i,j-1,2) * lam(i,j)
           else if (v(i,j-1,2) + v(i,j,2) .lt. 0.d0) then
              d(i,j) = dcf(i,j) * v(i,j,2) * lam(i,j)
           else
              d(i,j) = 0.0
           end if
           d(i,j) = d(i,j) * r(i)
        end do
     end do
  end if

end subroutine ca_compute_dcoefs


subroutine ca_update_dcf(lo, hi, &
     dcf, dcf_l1, dcf_l2, dcf_h1, dcf_h2, &
     etainv, eti_l1, eti_l2, eti_h1, eti_h2, &
     kp, kp_l1, kp_l2, kp_h1, kp_h2, kr, kr_l1, kr_l2, kr_h1, kr_h2)

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: dcf_l1, dcf_l2, dcf_h1, dcf_h2, eti_l1, eti_l2, eti_h1, eti_h2, &
       kp_l1, kp_l2, kp_h1, kp_h2, kr_l1, kr_l2, kr_h1, kr_h2
  double precision, intent(in) :: etainv(eti_l1:eti_h1, eti_l2:eti_h2)
  double precision, intent(in) :: kp(kp_l1:kp_h1, kp_l2:kp_h2)
  double precision, intent(in) :: kr(kr_l1:kr_h1, kr_l2:kr_h2)
  double precision             :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2)

  integer :: i, j

  do j=lo(2), hi(2)
     do i=lo(1), hi(1)
        dcf(i,j) = 2.d0 * etainv(i,j) * (kp(i,j)/kr(i,j))
     end do
  end do

end subroutine ca_update_dcf


subroutine ca_set_dterm_face( lo, hi, &
     Er, Er_l1, Er_l2, Er_h1, Er_h2, &
     dc, dc_l1, dc_l2, dc_h1, dc_h2, &
     dtf, dtf_l1, dtf_l2, dtf_h1, dtf_h2, dx, idir)
  implicit none
  
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2,  &
       dc_l1, dc_l2, dc_h1, dc_h2, dtf_l1, dtf_l2, dtf_h1, dtf_h2, idir
  double precision, intent(in) :: dx(2)
  double precision, intent(in) :: Er(Er_l1:Er_h1,Er_l2:Er_h2)
  double precision, intent(in) :: dc(dc_l1:dc_h1,dc_l2:dc_h2)
  double precision             :: dtf(dtf_l1:dtf_h1,dtf_l2:dtf_h2)
  integer :: i, j

  if (idir .eq. 0) then
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           dtf(i,j) = (Er(i,j) - Er(i-1,j)) / dx(1) * dc(i,j)
        end do
     end do
  else
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           dtf(i,j) = (Er(i,j) - Er(i,j-1)) / dx(2) * dc(i,j)
        end do
     end do
  end if

end subroutine ca_set_dterm_face


subroutine ca_face2center( lo, hi, &
     scomp, dcomp, ncomp, nf, nc, &
     foox, foox_l1, foox_l2, foox_h1, foox_h2, &
     fooy, fooy_l1, fooy_l2, fooy_h1, fooy_h2, &
     fooc, fooc_l1, fooc_l2, fooc_h1, fooc_h2)

  implicit none

  integer, intent(in) :: lo(2), hi(2), scomp,dcomp,ncomp,nf,nc
  integer, intent(in) :: foox_l1, foox_l2, foox_h1, foox_h2
  integer, intent(in) :: fooy_l1, fooy_l2, fooy_h1, fooy_h2
  integer, intent(in) :: fooc_l1, fooc_l2, fooc_h1, fooc_h2
  double precision, intent(in)  :: foox(foox_l1:foox_h1,foox_l2:foox_h2,0:nf-1)
  double precision, intent(in)  :: fooy(fooy_l1:fooy_h1,fooy_l2:fooy_h2,0:nf-1)
  double precision              :: fooc(fooc_l1:fooc_h1,fooc_l2:fooc_h2,0:nc-1)

  integer :: i,j,n

  do n = 0, ncomp-1
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           fooc(i,j,dcomp+n) = (foox(i,j,scomp+n) + foox(i+1,j,scomp+n) &
                &             + fooy(i,j,scomp+n) + fooy(i,j+1,scomp+n)) * 0.25d0
        end do
     end do
  end do

end subroutine ca_face2center


! no tiling
subroutine ca_correct_dterm(  & 
     dfx, dfx_l1, dfx_l2, dfx_h1, dfx_h2, &
     dfy, dfy_l1, dfy_l2, dfy_h1, dfy_h2, &
     re, rc)

  implicit none

  integer, intent(in) :: dfx_l1, dfx_l2, dfx_h1, dfx_h2
  integer, intent(in) :: dfy_l1, dfy_l2, dfy_h1, dfy_h2
  double precision, intent(inout) :: dfx(dfx_l1:dfx_h1,dfx_l2:dfx_h2)
  double precision, intent(inout) :: dfy(dfy_l1:dfy_h1,dfy_l2:dfy_h2)
  double precision, intent(in) :: re(dfx_l1:dfx_h1), rc(dfy_l1:dfy_h1)

  integer :: i, j

  do j=dfx_l2, dfx_h2
     do i=dfx_l1, dfx_h1
        dfx(i,j) = dfx(i,j) / (re(i) + 1.d-50)
     end do
  end do

  do j=dfy_l2, dfy_h2
     do i=dfy_l1, dfy_h1
        dfy(i,j) = dfy(i,j) / rc(i)
     end do
  end do

end subroutine ca_correct_dterm


subroutine ca_estdt_rad(u,u_l1,u_l2,u_h1,u_h2, &
     gpr,gpr_l1,gpr_l2,gpr_h1,gpr_h2, &
     lo,hi,dx,dt)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
       allow_negative_energy
  
  implicit none
  
  integer          :: u_l1,u_l2,u_h1,u_h2
  integer          :: gpr_l1,gpr_l2,gpr_h1,gpr_h2
  integer          :: lo(2), hi(2)
  double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
  double precision :: gpr(gpr_l1:gpr_h1,gpr_l2:gpr_h2)
  double precision :: dx(2),dt

  double precision :: rhoInv,ux,uy,dt1,dt2,c
  integer          :: i,j
  type(eos_t) :: eos_state
  
  !    Translate to primitive variables, compute sound speed (call eos), get dtmax
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)

        rhoInv = 1.d0 / u(i,j,URHO)

        eos_state % rho = u(i,j,URHO)
        eos_state % T   = u(i,j,UTEMP)
        eos_state % e   = u(i,j,UEINT)*rhoInv
        eos_state % xn  = u(i,j,UFS:UFS+nspec-1) * rhoInv
        eos_state % aux = u(i,j,UFX:UFX+naux -1) * rhoInv

        if (eos_state % e .gt. 0.d0 .or. allow_negative_energy.eq.1) then
           call eos(eos_input_re, eos_state)
           c = eos_state % cs
        else
           c = 0.d0
        end if

        c = sqrt(c**2 + gpr(i,j)*rhoInv)

        ux = u(i,j,UMX)*rhoInv
        uy = u(i,j,UMY)*rhoInv
                
        dt1 = dx(1)/(c + abs(ux))
        dt2 = dx(2)/(c + abs(uy))
        dt = min(dt,dt1,dt2)
     enddo
  enddo

end subroutine ca_estdt_rad


! this is tiling safe
subroutine ca_est_gpr0(Er, Er_l1, Er_l2, Er_h1, Er_h2, &
     gPr, gPr_l1, gPr_l2, gPr_h1, gPr_h2)

  use rad_params_module, only : ngroups

  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2, gpr_l1, gpr_l2, gpr_h1, gpr_h2
  double precision, intent(in) :: Er(Er_l1:Er_h1,Er_l2:Er_h2, 0:ngroups-1)
  double precision, intent(out) :: gPr(gPr_l1:gPr_h1,gPr_l2:gPr_h2)

  integer :: i, j, g

  gPr = 0.d0

  do g = 0, ngroups-1
     do j = gPr_l2, gPr_h2
        do i = gPr_l1, gPr_h1
           gPr(i,j) = gPr(i,j) + 4.d0/9.d0*Er(i,j,g)
        end do
     end do
  end do

end subroutine ca_est_gpr0


! this is tiling safe
subroutine ca_est_gpr2(kap, kap_l1, kap_l2, kap_h1, kap_h2, &
     Er, Er_l1, Er_l2, Er_h1, Er_h2, &
     gPr, gPr_l1, gPr_l2, gPr_h1, gPr_h2, vlo, vhi, dx, limiter, comoving)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda, Edd_factor

  implicit none

  integer, intent(in) :: kap_l1, kap_l2, kap_h1, kap_h2, &
       Er_l1, Er_l2, Er_h1, Er_h2, gPr_l1, gPr_l2, gPr_h1, gPr_h2
  integer, intent(in) :: vlo(2), vhi(2)  ! the region with valid Er
  integer, intent(in) :: limiter, comoving
  double precision, intent(in) :: dx(2)
  double precision, intent(in) :: kap(kap_l1:kap_h1,kap_l2:kap_h2, 0:ngroups-1), &
       Er(Er_l1:Er_h1,Er_l2:Er_h2, 0:ngroups-1)
  double precision, intent(out) :: gPr(gPr_l1:gPr_h1,gPr_l2:gPr_h2)

  integer :: i, j, g
  double precision :: gE(gPr_l1:gPr_h1,gPr_l2:gPr_h2)
  double precision :: lam, gE1, gE2, r, f, gamr 
  integer :: im, ip, jm, jp
  double precision :: xm, xp, ym, yp

  if (gPr_l1-1 .ge. vlo(1)) then
     im = 1
     xm = 2.d0
  else
     im = 0
     xm = 1.d0
  end if

  if (gPr_h1+1 .le. vhi(1)) then
     ip = 1
     xp = 2.d0
  else
     ip = 0
     xp = 1.d0
  end if
  
  if (gPr_l2-1 .ge. vlo(2)) then
     jm = 1
     ym = 2.d0
  else
     jm = 0
     ym = 1.d0
  end if

  if (gPr_h2+1 .le. vhi(2)) then
     jp = 1
     yp = 2.d0
  else
     jp = 0
     yp = 1.d0
  end if

  gPr = 0.0d0

  do g = 0, ngroups-1

     do j = gPr_l2+1, gPr_h2-1
        do i = gPr_l1+1, gPr_h1-1
           gE1 = (Er(i+1,j,g) - Er(i-1,j,g)) / (2.d0*dx(1))
           gE2 = (Er(i,j+1,g) - Er(i,j-1,g)) / (2.d0*dx(2))
           gE(i,j) = sqrt(gE1**2 + gE2**2)
        end do
     end do

     ! lo-x lo-y corner
     i = gPr_l1
     j = gPr_l2
     gE1 = (Er(i+1,j  ,g) - Er(i-im,j   ,g)) / (xm*dx(1))
     gE2 = (Er(i  ,j+1,g) - Er(i   ,j-jm,g)) / (ym*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)

     ! med-x lo-y side
     j = gPr_l2
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j  ,g) - Er(i-1,j   ,g)) / (2.d0*dx(1))
        gE2 = (Er(i  ,j+1,g) - Er(i  ,j-jm,g)) / (  ym*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do
     
     ! hi-x lo-y corner
     i = gPr_h1
     j = gPr_l2
     gE1 = (Er(i+ip,j  ,g) - Er(i-1,j   ,g)) / (xp*dx(1))
     gE2 = (Er(i   ,j+1,g) - Er(i  ,j-jm,g)) / (ym*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)
  
     ! lo-x med-y side
     i = gPr_l1
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i+1,j  ,g) - Er(i-im,j  ,g)) / (  xm*dx(1))
        gE2 = (Er(i  ,j+1,g) - Er(i   ,j-1,g)) / (2.d0*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do
     
     ! hi-x med-y side
     i = gPr_h1
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i+ip,j  ,g) - Er(i-1,j  ,g)) / (  xp*dx(1))
        gE2 = (Er(i   ,j+1,g) - Er(i  ,j-1,g)) / (2.d0*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do
     
     ! lo-x hi-y corner
     i = gPr_l1
     j = gPr_h2
     gE1 = (Er(i+1,j   ,g) - Er(i-im,j  ,g)) / (xm*dx(1))
     gE2 = (Er(i  ,j+jp,g) - Er(i   ,j-1,g)) / (yp*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)
     
     ! med-x hi-y side
     j = gPr_h2
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j   ,g) - Er(i-1,j,g)) / (2.d0*dx(1))
        gE2 = (Er(i  ,j+jp,g) - Er(i,j-1,g)) / (  yp*dx(2))
        gE(i,j) = sqrt(gE1**2 + gE2**2)
     end do
     
     ! hi-x hi-y corner
     i = gPr_h1
     j = gPr_h2
     gE1 = (Er(i+ip,j   ,g) - Er(i-1,j  ,g)) / (xp*dx(1))
     gE2 = (Er(i   ,j+jp,g) - Er(i  ,j-1,g)) / (yp*dx(2))
     gE(i,j) = sqrt(gE1**2 + gE2**2)
     
     do j = gPr_l2, gPr_h2
        do i = gPr_l1, gPr_h1
           r = gE(i,j) / (kap(i,j,g) * max(Er(i,j,g), 1.d-50))
           lam = FLDlambda(r, limiter)
           if (comoving .eq. 1) then
              f = Edd_factor(lam)
              gamr = (3.d0-f)/2.d0
           else
              gamr = lam + 1.d0
           end if
           gPr(i,j) = gPr(i,j) + lam * gamr * Er(i,j,g)
        end do
     end do
     
  end do

end subroutine ca_est_gpr2

