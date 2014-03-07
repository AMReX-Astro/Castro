subroutine ca_umdrv_rad(is_finest_level,time,lo,hi,domlo,domhi, &
     uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
     uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
     Erin,Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3, &
     lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
     Erout,Erout_l1,Erout_l2,Erout_l3,Erout_h1,Erout_h2,Erout_h3, &
     ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
     ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
     ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
     src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
     grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
     delta,dt, &
     flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
     flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
     flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
     radflux1,radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3, &
     radflux2,radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3, &
     radflux3,radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3, &
     area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
     area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
     area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
     vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
     courno,verbose, &
     nstep_fsp, fspace_advection_type, comoving)
  
  use meth_params_module, only : QVAR, QU, QPRES, NVAR, NHYP, URHO, use_colglaz, do_sponge, &
       normalize_species
  use rad_params_module, only : ngroups
  use radhydro_params_module, only : init_radhydro_params, QRADVAR
  use advection_module, only : enforce_minimum_density, normalize_new_species, divu
  use rad_advection_module, only : umeth3d_rad, ctoprim_rad, consup_rad
  use sponge_module, only : sponge

  implicit none

  integer nstep_fsp, fspace_advection_type, comoving
  integer is_finest_level
  integer lo(3),hi(3),verbose
  integer domlo(3),domhi(3)
  integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
  integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
  integer Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3
  integer Erout_l1,Erout_l2,Erout_l3,Erout_h1,Erout_h2,Erout_h3
  integer lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3
  integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
  integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
  integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
  integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
  integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
  integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
  integer radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3
  integer radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3
  integer radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3
  integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
  integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
  integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
  integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
  integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
  integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
  double precision   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,     uin_l3:uin_h3,  NVAR)
  double precision  uout( uout_l1:uout_h1,  uout_l2:uout_h2,   uout_l3:uout_h3, NVAR)
  double precision Erin ( Erin_l1: Erin_h1, Erin_l2: Erin_h2, Erin_l3: Erin_h3,0:ngroups-1)
  double precision Erout(Erout_l1:Erout_h1,Erout_l2:Erout_h2,Erout_l3:Erout_h3,0:ngroups-1)
  double precision lam  (  lam_l1:  lam_h1,  lam_l2:  lam_h2,  lam_l3:  lam_h3,0:ngroups-1)
  double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
  double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
  double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
  double precision   src(  src_l1:src_h1,    src_l2:src_h2,     src_l3:src_h3,  NVAR)
  double precision  grav( gv_l1:gv_h1,  gv_l2:gv_h2,   gv_l3:gv_h3,    3)
  double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
  double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
  double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)
  double precision radflux1(radflux1_l1:radflux1_h1,radflux1_l2:radflux1_h2, &
                            radflux1_l3:radflux1_h3,0:ngroups-1)  
  double precision radflux2(radflux2_l1:radflux2_h1,radflux2_l2:radflux2_h2, &
                            radflux2_l3:radflux2_h3,0:ngroups-1)  
  double precision radflux3(radflux3_l1:radflux3_h1,radflux3_l2:radflux3_h2, &
                            radflux3_l3:radflux3_h3,0:ngroups-1)  
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2, area1_l3:area1_h3)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2, area2_l3:area2_h3)
  double precision area3(area3_l1:area3_h1,area3_l2:area3_h2, area3_l3:area3_h3)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2, vol_l3:vol_h3)
  double precision delta(3),dt,time,courno
  
  !     Automatic arrays for workspace
  double precision, allocatable:: q(:,:,:,:)
  double precision, allocatable:: gamc(:,:,:)
  double precision, allocatable:: gamcg(:,:,:)
  double precision, allocatable:: flatn(:,:,:)
  double precision, allocatable:: c(:,:,:)
  double precision, allocatable:: cg(:,:,:)
  double precision, allocatable:: csml(:,:,:)
  double precision, allocatable:: div(:,:,:)
  double precision, allocatable:: pdivu(:,:,:)
  double precision, allocatable:: srcQ(:,:,:,:)
  double precision, allocatable:: ergdx(:,:,:,:)
  double precision, allocatable:: ergdy(:,:,:,:)
  double precision, allocatable:: ergdz(:,:,:,:)
  double precision, allocatable:: lmgdx(:,:,:,:)
  double precision, allocatable:: lmgdy(:,:,:,:)
  double precision, allocatable:: lmgdz(:,:,:,:)

  double precision, allocatable ::uy_xfc(:,:,:)
  double precision, allocatable ::uz_xfc(:,:,:)
  double precision, allocatable ::ux_yfc(:,:,:)
  double precision, allocatable ::uz_yfc(:,:,:)
  double precision, allocatable ::ux_zfc(:,:,:)
  double precision, allocatable ::uy_zfc(:,:,:)

  integer ngq,ngf,iflaten
  double precision dx,dy,dz, mass_added,eint_added,eden_added

  logical, save :: first_call = .true.

  if (first_call) then
     call init_radhydro_params(fspace_advection_type, comoving)
     first_call = .false.
  end if
  
  allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,QRADVAR))
  allocate(  gamc(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate( gamcg(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate( flatn(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate(     c(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate(    cg(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate(  csml(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate(   div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))
  
  allocate( pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

  allocate(  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR))

  allocate( ergdx(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3,0:ngroups-1))
  allocate( ergdy(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3,0:ngroups-1))
  allocate( ergdz(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3,0:ngroups-1))

  allocate( lmgdx(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3,0:ngroups-1))
  allocate( lmgdy(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3,0:ngroups-1))
  allocate( lmgdz(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3,0:ngroups-1))

  allocate (uy_xfc(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3))
  allocate (uz_xfc(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3))
  allocate (ux_yfc(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3))
  allocate (uz_yfc(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3))
  allocate (ux_zfc(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3))
  allocate (uy_zfc(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3))
 
  dx = delta(1)
  dy = delta(2)
  dz = delta(3)

  ngq = NHYP
  ngf = 1
  iflaten = 1
  
  !     Translate to primitive variables, compute sound speeds
  !     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
  !       and set to correspond to coordinates of (lo:hi)
  call ctoprim_rad(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       Erin,Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3, &
       lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       q,c,cg,gamc,gamcg,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       src,srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
       courno,dx,dy,dz,dt,ngq,ngf,iflaten)

!     Compute hyperbolic fluxes using unsplit Godunov
  call umeth3d_rad(q,c,cg,gamc,gamcg,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       lam,lam_l1,lam_l2,lam_l3,lam_h1,lam_h2,lam_h3, &
       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
       grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
       lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       radflux1,radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3, &
       radflux2,radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3, &
       radflux3,radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3, &
       ugdnvx_out, ergdx, lmgdx, &
       ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
       ugdnvy_out, ergdy, lmgdy, &
       ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
       ugdnvz_out, ergdz, lmgdz, &
       ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
       pdivu, uy_xfc, uz_xfc, ux_yfc, uz_yfc, ux_zfc, uy_zfc)

  !     Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo,hi,q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       dx,dy,dz,div,lo(1),lo(2),lo(3),hi(1)+1,hi(2)+1,hi(3)+1)

  !     Conservative update
  call consup_rad(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
       Erin,Erin_l1,Erin_l2,Erin_l3,Erin_h1,Erin_h2,Erin_h3, &
       Erout,Erout_l1,Erout_l2,Erout_l3,Erout_h1,Erout_h2,Erout_h3, &
       ugdnvx_out, ergdx, lmgdx, &
       ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
       ugdnvy_out, ergdy, lmgdy, &
       ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
       ugdnvz_out, ergdz, lmgdz, &
       ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
       src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
       grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       radflux1,radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3, &
       radflux2,radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3, &
       radflux3,radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3, &
       area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
       area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
       area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
       vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
       div,pdivu, uy_xfc, uz_xfc, ux_yfc, uz_yfc, ux_zfc, uy_zfc, &
       lo,hi,dx,dy,dz,dt, nstep_fsp)

  ! Enforce the density >= small_dens.
  mass_added = 0.d0
  eint_added = 0.d0
  eden_added = 0.d0
  call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
       uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
       lo,hi,mass_added,eint_added,eden_added,verbose)

  ! Enforce species >= 0
  call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
       uout_h1,uout_h2,uout_h3,lo,hi)
 
  ! Re-normalize the species
  if (normalize_species .eq. 1) then
     call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
          lo,hi)
  end if
      
  ! Impose sponge
  if (do_sponge .eq. 1) then
     call sponge(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,lo,hi, &
          time,dt,dx,dy,dz,domlo,domhi)
  end if
  
  deallocate(q,gamc,gamcg,flatn,c,cg,csml,div,srcQ,pdivu,ergdx,ergdy,ergdz,lmgdx,lmgdy,lmgdz)
  deallocate(uy_xfc,uz_xfc,ux_yfc,uz_yfc,ux_zfc,uy_zfc)

end subroutine ca_umdrv_rad

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine ca_compute_lamborder(Er, Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
     kap, kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     lam, lam_l1, lam_l2, lam_l3, lam_h1, lam_h2, lam_h3, &
     dx, ngrow, limiter, filter_T, S)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda
  use filter_module

  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
       kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
       lam_l1, lam_l2, lam_l3, lam_h1, lam_h2, lam_h3
  integer, intent(in) :: ngrow, limiter, filter_T, S
  double precision, intent(in) :: dx(3)
  double precision, intent(in) :: kap(kap_l1:kap_h1, kap_l2:kap_h2, kap_l3:kap_h3, 0:ngroups-1)
  double precision, intent(in) :: Er(Er_l1:Er_h1, Er_l2:Er_h2, Er_l3:Er_h3, 0:ngroups-1)
  double precision, intent(out) :: lam(lam_l1:lam_h1, lam_l2:lam_h2, lam_l3:lam_h3, 0:ngroups-1)

  integer :: i, j, k, reg_l1, reg_l2, reg_l3, reg_h1, reg_h2, reg_h3, g
  double precision :: r, r1, r2, r3

  double precision, allocatable :: lamfil(:,:,:)

  reg_l1 = lam_l1 + ngrow
  reg_l2 = lam_l2 + ngrow
  reg_l3 = lam_l3 + ngrow
  reg_h1 = lam_h1 - ngrow
  reg_h2 = lam_h2 - ngrow
  reg_h3 = lam_h3 - ngrow

  !$omp parallel if (ngroups>1) private(g,i,j,k,r,r1,r2,r3,lamfil)

  if (filter_T .gt. 0) then
     allocate(lamfil(reg_l1:reg_h1,lam_l2:lam_h2,lam_l3:lam_h3))
  end if

  !$omp do
  do g = 0, ngroups-1

  do k=lam_l3, lam_h3
     do j=lam_l2, lam_h2
        do i=lam_l1, lam_h1

           lam(i,j,k,g) = -1.d50

           if (Er(i,j,k,g) .eq. -1.d0) then
              cycle
           end if
           
           if (Er(i-1,j,k,g) .eq. -1.d0) then
              r1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
           else if (Er(i+1,j,k,g) .eq. -1.d0) then
              r1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
           else
              r1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
           end if

           if (Er(i,j-1,k,g) .eq. -1.d0) then
              r2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
           else if (Er(i,j+1,k,g) .eq. -1.d0) then
              r2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
           else
              r2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
           end if

           if (Er(i,j,k-1,g) .eq. -1.d0) then
              r3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / dx(3)
           else if (Er(i,j,k+1,g) .eq. -1.d0) then
              r3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / dx(3)
           else
              r3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
           end if

           r = sqrt(r1**2 + r2**2 + r3**2)
           r = r / (kap(i,j,k,g) * max(Er(i,j,k,g), 1.d-50))

           lam(i,j,k,g) = FLDlambda(r, limiter)
        end do
     end do
  end do

  ! filter
  if (filter_T .eq. 1) then
     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.d0) then
              lamfil(:,j,k) = -1.d-50
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff1(0) * lam(i,j,k,g) &
                   &        + ff1(1) * (lam(i-1,j,k,g)+lam(i+1,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.d0) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff1b, lam(i:i+1,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.d0) then
              i = reg_h1
              lamfil(i,j,k) = dot_product(ff1b(1:0:-1), lam(i-1:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.d0) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff1(0) * lamfil(i,j,k) &
                   &       + ff1(1) * (lamfil(i,j-1,k)+lamfil(i,j+1,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.d0) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff1b, lamfil(i,j:j+1,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.d0) then
           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff1b(1:0:-1), lamfil(i,j-1:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff1(0) * lam(i,j,k,g) &
                   &        + ff1(1) * (lam(i,j,k-1,g)+lam(i,j,k+1,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.d0) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff1b, lam(i,j,k:k+1,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.d0) then
        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff1b(1:0:-1), lam(i,j,k-1:k,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  else if (filter_T .eq. 2) then

     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.d0) then
              lamfil(:,j,k) = -1.d-50
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff2(0,S) * lam(i,j,k,g) &
                   &        + ff2(1,S) * (lam(i-1,j,k,g)+lam(i+1,j,k,g)) &
                   &        + ff2(2,S) * (lam(i-2,j,k,g)+lam(i+2,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.d0) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff2b0, lam(i:i+2,j,k,g))

              i = reg_l1 + 1
              lamfil(i,j,k) = dot_product(ff2b1, lam(i-1:i+2,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.d0) then
              i = reg_h1 - 1
              lamfil(i,j,k) = dot_product(ff2b1(2:-1:-1), lam(i-2:i+1,j,k,g))

              i = reg_h1 
              lamfil(i,j,k) = dot_product(ff2b0(2:0:-1), lam(i-2:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.d0) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff2(0,S) * lamfil(i,j,k) &
                   &       + ff2(1,S) * (lamfil(i,j-1,k)+lamfil(i,j+1,k)) &
                   &       + ff2(2,S) * (lamfil(i,j-2,k)+lamfil(i,j+2,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.d0) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b0, lamfil(i,j:j+2,k))
           end do

           j = reg_l2 + 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b1, lamfil(i,j-1:j+2,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.d0) then
           j = reg_h2 - 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b1(2:-1:-1), lamfil(i,j-2:j+1,k))
           end do

           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff2b0(2:0:-1), lamfil(i,j-2:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff2(0,S) * lam(i,j,k,g) &
                   &        + ff2(1,S) * (lam(i,j,k-1,g)+lam(i,j,k+1,g)) &
                   &        + ff2(2,S) * (lam(i,j,k-2,g)+lam(i,j,k+2,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.d0) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b0, lam(i,j,k:k+2,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b1, lam(i,j,k-1:k+2,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.d0) then
        k = reg_h3 - 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b1(2:-1:-1), lam(i,j,k-2:k+1,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff2b0(2:0:-1), lam(i,j,k-2:k,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  else if (filter_T .eq. 3) then

     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.d0) then
              lamfil(:,j,k) = -1.d-50
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff3(0,S) * lam(i,j,k,g) &
                   &        + ff3(1,S) * (lam(i-1,j,k,g)+lam(i+1,j,k,g)) &
                   &        + ff3(2,S) * (lam(i-2,j,k,g)+lam(i+2,j,k,g)) &
                   &        + ff3(3,S) * (lam(i-3,j,k,g)+lam(i+3,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.d0) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff3b0, lam(i:i+3,j,k,g))

              i = reg_l1 + 1
              lamfil(i,j,k) = dot_product(ff3b1, lam(i-1:i+3,j,k,g))

              i = reg_l1 + 2
              lamfil(i,j,k) = dot_product(ff3b2, lam(i-2:i+3,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.d0) then
              i = reg_h1 - 2
              lamfil(i,j,k) = dot_product(ff3b2(3:-2:-1), lam(i-3:i+2,j,k,g))

              i = reg_h1 - 1
              lamfil(i,j,k) = dot_product(ff3b1(3:-1:-1), lam(i-3:i+1,j,k,g))

              i = reg_h1 
              lamfil(i,j,k) = dot_product(ff3b0(3:0:-1), lam(i-3:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.d0) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff3(0,S) * lamfil(i,j,k) &
                   &       + ff3(1,S) * (lamfil(i,j-1,k)+lamfil(i,j+1,k)) &
                   &       + ff3(2,S) * (lamfil(i,j-2,k)+lamfil(i,j+2,k)) &
                   &       + ff3(3,S) * (lamfil(i,j-3,k)+lamfil(i,j+3,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.d0) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b0, lamfil(i,j:j+3,k))
           end do

           j = reg_l2 + 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b1, lamfil(i,j-1:j+3,k))
           end do

           j = reg_l2 + 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b2, lamfil(i,j-2:j+3,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.d0) then
           j = reg_h2 - 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b2(3:-2:-1), lamfil(i,j-3:j+2,k))
           end do

           j = reg_h2 - 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b1(3:-1:-1), lamfil(i,j-3:j+1,k))
           end do

           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff3b0(3:0:-1), lamfil(i,j-3:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff3(0,S) * lam(i,j,k,g) &
                   &        + ff3(1,S) * (lam(i,j,k-1,g)+lam(i,j,k+1,g)) &
                   &        + ff3(2,S) * (lam(i,j,k-2,g)+lam(i,j,k+2,g)) &
                   &        + ff3(3,S) * (lam(i,j,k-3,g)+lam(i,j,k+3,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.d0) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b0, lam(i,j,k:k+3,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b1, lam(i,j,k-1:k+3,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b2, lam(i,j,k-2:k+3,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.d0) then
        k = reg_h3 - 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b2(3:-2:-1), lam(i,j,k-3:k+2,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3 - 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b1(3:-1:-1), lam(i,j,k-3:k+1,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff3b0(3:0:-1), lam(i,j,k-3:k,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  else if (filter_T .eq. 4) then

     do k=lam_l3, lam_h3
        do j=lam_l2, lam_h2
           if (Er(reg_l1,j,k,g) .eq. -1.d0) then
              lamfil(:,j,k) = -1.d-50
              cycle
           end if

           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff4(0,S) * lam(i,j,k,g) &
                   &        + ff4(1,S) * (lam(i-1,j,k,g)+lam(i+1,j,k,g)) &
                   &        + ff4(2,S) * (lam(i-2,j,k,g)+lam(i+2,j,k,g)) &
                   &        + ff4(3,S) * (lam(i-3,j,k,g)+lam(i+3,j,k,g)) &
                   &        + ff4(4,S) * (lam(i-4,j,k,g)+lam(i+4,j,k,g))
           end do

           if (Er(reg_l1-1,j,k,g) .eq. -1.d0) then
              i = reg_l1
              lamfil(i,j,k) = dot_product(ff4b0, lam(i:i+4,j,k,g))

              i = reg_l1 + 1
              lamfil(i,j,k) = dot_product(ff4b1, lam(i-1:i+4,j,k,g))

              i = reg_l1 + 2
              lamfil(i,j,k) = dot_product(ff4b2, lam(i-2:i+4,j,k,g))

              i = reg_l1 + 3
              lamfil(i,j,k) = dot_product(ff4b3, lam(i-3:i+4,j,k,g))
           end if

           if (Er(reg_h1+1,j,k,g) .eq. -1.d0) then
              i = reg_h1 - 3 
              lamfil(i,j,k) = dot_product(ff4b3(4:-3:-1), lam(i-4:i+3,j,k,g))

              i = reg_h1 - 2
              lamfil(i,j,k) = dot_product(ff4b2(4:-2:-1), lam(i-4:i+2,j,k,g))

              i = reg_h1 - 1
              lamfil(i,j,k) = dot_product(ff4b1(4:-1:-1), lam(i-4:i+1,j,k,g))

              i = reg_h1 
              lamfil(i,j,k) = dot_product(ff4b0(4:0:-1), lam(i-4:i,j,k,g))
           end if
        end do
     end do

     do k=lam_l3, lam_h3
        if (Er(reg_l1,reg_l2,k,g) .eq. -1.d0) then
           cycle
        end if

        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lam(i,j,k,g) = ff4(0,S) * lamfil(i,j,k) &
                   &       + ff4(1,S) * (lamfil(i,j-1,k)+lamfil(i,j+1,k)) &
                   &       + ff4(2,S) * (lamfil(i,j-2,k)+lamfil(i,j+2,k)) &
                   &       + ff4(3,S) * (lamfil(i,j-3,k)+lamfil(i,j+3,k)) &
                   &       + ff4(4,S) * (lamfil(i,j-4,k)+lamfil(i,j+4,k))
           end do
        end do

        if (Er(reg_l1,reg_l2-1,k,g) .eq. -1.d0) then
           j = reg_l2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b0, lamfil(i,j:j+4,k))
           end do

           j = reg_l2 + 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b1, lamfil(i,j-1:j+4,k))
           end do

           j = reg_l2 + 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b2, lamfil(i,j-2:j+4,k))
           end do

           j = reg_l2 + 3
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b3, lamfil(i,j-3:j+4,k))
           end do
        end if

        if (Er(reg_l1,reg_h2+1,k,g) .eq. -1.d0) then
           j = reg_h2 - 3
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b3(4:-3:-1), lamfil(i,j-4:j+3,k))
           end do

           j = reg_h2 - 2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b2(4:-2:-1), lamfil(i,j-4:j+2,k))
           end do

           j = reg_h2 - 1
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b1(4:-1:-1), lamfil(i,j-4:j+1,k))
           end do

           j = reg_h2
           do i=reg_l1,reg_h1
              lam(i,j,k,g) = dot_product(ff4b0(4:0:-1), lamfil(i,j-4:j,k))
           end do
        end if
     end do

     do k=reg_l3, reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = ff4(0,S) * lam(i,j,k,g) &
                   &        + ff4(1,S) * (lam(i,j,k-1,g)+lam(i,j,k+1,g)) &
                   &        + ff4(2,S) * (lam(i,j,k-2,g)+lam(i,j,k+2,g)) &
                   &        + ff4(3,S) * (lam(i,j,k-3,g)+lam(i,j,k+3,g)) &
                   &        + ff4(4,S) * (lam(i,j,k-4,g)+lam(i,j,k+4,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end do

     if (Er(reg_l1,reg_l2,reg_l3-1,g) .eq. -1.d0) then
        k = reg_l3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b0, lam(i,j,k:k+4,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b1, lam(i,j,k-1:k+4,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b2, lam(i,j,k-2:k+4,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_l3 + 3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b3, lam(i,j,k-3:k+4,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     if (Er(reg_l1,reg_l2,reg_h3+1,g) .eq. -1.d0) then
        k = reg_h3 - 3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b3(4:-3:-1), lam(i,j,k-4:k+3,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3 - 2
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b2(4:-2:-1), lam(i,j,k-4:k+2,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3 - 1
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b1(4:-1:-1), lam(i,j,k-4:k+1,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do

        k = reg_h3
        do j=reg_l2, reg_h2
           do i=reg_l1, reg_h1
              lamfil(i,j,k) = dot_product(ff4b0(4:0:-1), lam(i,j,k-4:k,g))
              lamfil(i,j,k) = min(1.d0/3.d0, max(1.d-25, lamfil(i,j,k)))
           end do
        end do
     end if

     lam(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3,g) = &
          lamfil(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)

  end if

  ! lo-x lo-y lo-z
  do k=lam_l3,reg_l3-1
     do j=lam_l2,reg_l2-1
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,reg_l2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! reg-x lo-y lo-z
  do k=lam_l3,reg_l3-1
     do j=lam_l2,reg_l2-1
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,reg_l2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! hi-x lo-y lo-z
  do k=lam_l3,reg_l3-1
     do j=lam_l2,reg_l2-1
        do i=reg_h1+1,lam_h1 
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,reg_l2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! lo-x reg-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_l2,reg_h2
        do i=lam_l1,reg_l1-1
           lam(i,j,k,g) = lam(reg_l1,j,reg_l3,g)
        end do
     end do
  end do

  ! reg-x reg-y lo-z side
  do k=lam_l3,reg_l3-1
     do j=reg_l2,reg_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,j,reg_l3,g)
           end if
        end do
     end do
  end do

  ! hi-x reg-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_l2,reg_h2
        do i=reg_h1+1,lam_h1 
           lam(i,j,k,g) = lam(reg_h1,j,reg_l3,g)
        end do
     end do
  end do

  ! lo-x hi-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_h2+1,lam_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,reg_h2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! reg-x hi-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_h2+1,lam_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,reg_h2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! hi-x hi-y lo-z
  do k=lam_l3,reg_l3-1
     do j=reg_h2+1,lam_h2
        do i=reg_h1+1,lam_h1 
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,reg_h2,reg_l3,g)
           end if
        end do
     end do
  end do

  ! lo-x lo-y reg-z
  do k=reg_l3,reg_h3
     do j=lam_l2,reg_l2-1
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,reg_l2,k,g)
           end if
        end do
     end do
  end do

  ! reg-x lo-y reg-z
  do k=reg_l3,reg_h3
     do j=lam_l2,reg_l2-1
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,reg_l2,k,g)
           end if
        end do
     end do
  end do

  ! hi-x lo-y reg-z
  do k=reg_l3,reg_h3
     do j=lam_l2,reg_l2-1
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,reg_l2,k,g)
           end if
        end do
     end do
  end do

  ! lo-x reg-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_l2,reg_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,j,k,g)
           end if
        end do
     end do
  end do

  ! hi-x reg-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_l2,reg_h2
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,j,k,g)
           end if
        end do
     end do
  end do

  ! lo-x hi-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_h2+1,lam_h2 
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,reg_h2,k,g)
           end if
        end do
     end do
  end do

  ! reg-x hi-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_h2+1,lam_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,reg_h2,k,g)
           end if
        end do
     end do
  end do

  ! hi-x hi-y reg-z
  do k=reg_l3,reg_h3
     do j=reg_h2+1,lam_h2
        do i=reg_h1+1,lam_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,reg_h2,k,g)
           end if
        end do
     end do
  end do

  ! lo-x lo-y hi-z
  do k=reg_h3+1,lam_h3
     do j=lam_l2,reg_l2-1
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,reg_l2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! reg-x lo-y hi-z
  do k=reg_h3+1,lam_h3
     do j=lam_l2,reg_l2-1
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,reg_l2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! hi-x lo-y hi-z
  do k=reg_h3+1,lam_h3
     do j=lam_l2,reg_l2-1
        do i=reg_h1+1,lam_h1 
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,reg_l2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! lo-x reg-y hi-z
  do k=reg_h3+1,lam_h3 
     do j=reg_l2,reg_h2
        do i=lam_l1,reg_l1-1
           lam(i,j,k,g) = lam(reg_l1,j,reg_h3,g)
        end do
     end do
  end do

  ! reg-x reg-y hi-z 
  do k=reg_h3+1,lam_h3
     do j=reg_l2,reg_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,j,reg_h3,g)
           end if
        end do
     end do
  end do

  ! hi-x reg-y hi-z
  do k=reg_h3+1,lam_h3 
     do j=reg_l2,reg_h2
        do i=reg_h1+1,lam_h1
           lam(i,j,k,g) = lam(reg_h1,j,reg_h3,g)
        end do
     end do
  end do

  ! lo-x hi-y hi-z
  do k=reg_h3+1,lam_h3 
     do j=reg_h2+1,lam_h2
        do i=lam_l1,reg_l1-1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_l1,reg_h2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! reg-x hi-y hi-z
  do k=reg_h3+1,lam_h3
     do j=reg_h2+1,lam_h2
        do i=reg_l1,reg_h1
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(i,reg_h2,reg_h3,g)
           end if
        end do
     end do
  end do

  ! hi-x hi-y hi-z
  do k=reg_h3+1,lam_h3 
     do j=reg_h2+1,lam_h2
        do i=reg_h1+1,lam_h1 
           if (Er(i,j,k,g).eq.-1.d0) then
              lam(i,j,k,g) = lam(reg_h1,reg_h2,reg_h3,g)
           end if
        end do
     end do
  end do

  end do ! do g=
  !$omp end do

  if (filter_T .gt. 0) then
     deallocate(lamfil)
  end if

  !$omp end parallel

  return
end subroutine ca_compute_lamborder


subroutine ca_get_v_dcf( &
     er ,  er_l1,  er_l2,  er_l3,  er_h1,  er_h2,  er_h3, &
     s  ,   s_l1,   s_l2,   s_l3,   s_h1,   s_h2,   s_h3, &
     T  ,   T_l1,   T_l2,   T_l3,   T_h1,   T_h2,   T_h3, &
     c_v, c_v_l1, c_v_l2, c_v_l3, c_v_h1, c_v_h2, c_v_h3, &
     kr ,  kr_l1,  kr_l2,  kr_l3,  kr_h1,  kr_h2,  kr_h3, &
     kp ,  kp_l1,  kp_l2,  kp_l3,  kp_h1,  kp_h2,  kp_h3, &
     kp2, kp2_l1, kp2_l2, kp2_l3, kp2_h1, kp2_h2, kp2_h3, &
     dtemp, dtime, sigma, c, &
     v  ,   v_l1,   v_l2,   v_l3,   v_h1,   v_h2,   v_h3, &
     dcf, dcf_l1, dcf_l2, dcf_l3, dcf_h1, dcf_h2, dcf_h3)

  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ

  implicit none

  integer, intent(in) :: er_l1,er_l2,er_l3,er_h1,er_h2,er_h3
  integer, intent(in) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
  integer, intent(in) :: T_l1,T_l2,T_l3,T_h1,T_h2,T_h3
  integer, intent(in) :: c_v_l1,c_v_l2,c_v_l3,c_v_h1,c_v_h2,c_v_h3
  integer, intent(in) :: kr_l1,kr_l2,kr_l3,kr_h1,kr_h2,kr_h3
  integer, intent(in) :: kp_l1,kp_l2,kp_l3,kp_h1,kp_h2,kp_h3
  integer, intent(in) :: kp2_l1,kp2_l2,kp2_l3,kp2_h1,kp2_h2,kp2_h3
  integer, intent(in) :: v_l1,v_l2,v_l3,v_h1,v_h2,v_h3
  integer, intent(in) :: dcf_l1,dcf_l2,dcf_l3,dcf_h1,dcf_h2,dcf_h3
  double precision, intent(in)  ::  er( er_l1: er_h1,  er_l2: er_h2,  er_l3: er_h3)
  double precision, intent(in)  ::   s(  s_l1:  s_h1,   s_l2:  s_h2,   s_l3:  s_h3, NVAR)
  double precision, intent(in)  ::   T(  T_l1:  T_h1,   T_l2:  T_h2,   T_l3:  T_h3)
  double precision, intent(in)  :: c_v(c_v_l1:c_v_h1, c_v_l2:c_v_h2, c_v_l3:c_v_h3)
  double precision, intent(in ) ::  kr( kr_l1: kr_h1,  kr_l2: kr_h2,  kr_l3: kr_h3)
  double precision, intent(in ) ::  kp( kp_l1: kp_h1,  kp_l2: kp_h2,  kp_l3: kp_h3)
  double precision, intent(in ) :: kp2(kp2_l1:kp2_h1, kp2_l2:kp2_h2, kp2_l3:kp2_h3)
  double precision, intent(in) :: dtemp, dtime, sigma, c
  double precision, intent(out) ::   v(  v_l1:  v_h1,   v_l2:  v_h2,   v_l3:  v_h3, 3)
  double precision, intent(out) :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2, dcf_l3:dcf_h3)

  integer :: i, j, k
  double precision :: etainv, fac0, fac2, alpha, frc

  fac0 = 4.d0 * sigma * dtime / dtemp
  fac2 = c * dtime / dtemp

  do k=v_l3, v_h3
     do j=v_l2, v_h2
        do i=v_l1, v_h1
           v(i,j,k,1) = s(i,j,k,UMX)/s(i,j,k,URHO)
           v(i,j,k,2) = s(i,j,k,UMY)/s(i,j,k,URHO)
           v(i,j,k,3) = s(i,j,k,UMZ)/s(i,j,k,URHO)
           
           alpha = fac0 * (kp2(i,j,k) * (T(i,j,k) + dtemp) ** 4    &
                -          kp (i,j,k) * (T(i,j,k)        ) ** 4)   &
                -  fac2 * (kp2(i,j,k) - kp(i,j,k)) * er(i,j,k)
           
           frc = s(i,j,k,URHO) * c_v(i,j,k) + 1.0d-50
           etainv = frc / (alpha + frc)
           
           dcf(i,j,k) = 2.d0 * etainv * (kp(i,j,k) / kr(i,j,k))
        end do
     end do
  end do

end subroutine ca_get_v_dcf


subroutine ca_compute_dcoefs(d, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
     lam, lam_l1, lam_l2, lam_l3, lam_h1, lam_h2, lam_h3, &
     v ,    v_l1,   v_l2,   v_l3,   v_h1,   v_h2,   v_h3, &
     dcf, dcf_l1, dcf_l2, dcf_l3, dcf_h1, dcf_h2, dcf_h3, &
     r, idir)

  implicit none

  integer, intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
       & lam_l1, lam_l2, lam_l3, lam_h1, lam_h2, lam_h3, &
       &   v_l1,   v_l2,   v_l3,   v_h1,   v_h2,   v_h3, &
       & dcf_l1, dcf_l2, dcf_l3, dcf_h1, dcf_h2, dcf_h3, &
       idir

  double precision, intent(out) ::   d(  d_l1:  d_h1,   d_l2:  d_h2,   d_l3:  d_h3)
  double precision, intent(in)  :: lam(lam_l1:lam_h1, lam_l2:lam_h2, lam_l3:lam_h3)
  double precision, intent(in)  ::   v(  v_l1:  v_h1,   v_l2:  v_h2,   v_l3:  v_h3, 3)
  double precision, intent(in)  :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2, dcf_l3:dcf_h3)
  double precision, intent(in)  ::   r(  d_l1:  d_h1)

  integer :: i, j, k

  if (idir.eq.0) then
     do k = d_l3, d_h3
        do j = d_l2, d_h2
           do i = d_l1, d_h1
              if (v(i-1,j,k,1) + v(i,j,k,1) .gt. 0.d0) then
                 d(i,j,k) = dcf(i-1,j,k) * v(i-1,j,k,1) * lam(i,j,k)
              else if (v(i-1,j,k,1) + v(i,j,k,1) .lt. 0.d0) then
                 d(i,j,k) = dcf(i,j,k) * v(i,j,k,1) * lam(i,j,k)
              else
                 d(i,j,k) = 0.0
              end if
           end do
        end do
     end do
  else if (idir.eq.1) then
     do k = d_l3, d_h3
        do j = d_l2, d_h2
           do i = d_l1, d_h1
              if (v(i,j-1,k,2) + v(i,j,k,2) .gt. 0.d0) then
                 d(i,j,k) = dcf(i,j-1,k) * v(i,j-1,k,2) * lam(i,j,k)
              else if (v(i,j-1,k,2) + v(i,j,k,2) .lt. 0.d0) then
                 d(i,j,k) = dcf(i,j,k) * v(i,j,k,2) * lam(i,j,k)
              else
                 d(i,j,k) = 0.0
              end if
           end do
        end do
     end do
  else
     do k = d_l3, d_h3
        do j = d_l2, d_h2
           do i = d_l1, d_h1
              if (v(i,j,k-1,3) + v(i,j,k,3) .gt. 0.d0) then
                 d(i,j,k) = dcf(i,j,k-1) * v(i,j,k-1,3) * lam(i,j,k)
              else if (v(i,j,k-1,3) + v(i,j,k,3) .lt. 0.d0) then
                 d(i,j,k) = dcf(i,j,k) * v(i,j,k,3) * lam(i,j,k)
              else
                 d(i,j,k) = 0.0
              end if
           end do
        end do
     end do
  end if

end subroutine ca_compute_dcoefs


subroutine ca_update_dcf(dcf, dcf_l1, dcf_l2, dcf_l3, dcf_h1, dcf_h2, dcf_h3, &
     etainv, eti_l1, eti_l2, eti_l3, eti_h1, eti_h2, eti_h3, &
     kp, kp_l1, kp_l2, kp_l3, kp_h1, kp_h2, kp_h3, &
     kr, kr_l1, kr_l2, kr_l3, kr_h1, kr_h2, kr_h3)

  implicit none

  integer, intent(in) :: dcf_l1, dcf_l2, dcf_l3, dcf_h1, dcf_h2, dcf_h3
  integer, intent(in) :: eti_l1, eti_l2, eti_l3, eti_h1, eti_h2, eti_h3
  integer, intent(in) :: kp_l1, kp_l2, kp_l3, kp_h1, kp_h2, kp_h3
  integer, intent(in) :: kr_l1, kr_l2, kr_l3, kr_h1, kr_h2, kr_h3
  double precision, intent(in) :: etainv(eti_l1:eti_h1, eti_l2:eti_h2, eti_l3:eti_h3)
  double precision, intent(in) :: kp(kp_l1:kp_h1, kp_l2:kp_h2, kp_l3:kp_h3)
  double precision, intent(in) :: kr(kr_l1:kr_h1, kr_l2:kr_h2, kr_l3:kr_h3)
  double precision, intent(out) :: dcf(dcf_l1:dcf_h1, dcf_l2:dcf_h2, dcf_l3:dcf_h3)

  integer :: i, j, k

  do k=dcf_l3+1, dcf_h3-1
     do j=dcf_l2+1, dcf_h2-1
        do i=dcf_l1+1, dcf_h1-1
           dcf(i,j,k) = 2.d0 * etainv(i,j,k) * (kp(i,j,k)/kr(i,j,k))
        end do
     end do
  end do

end subroutine ca_update_dcf


subroutine ca_set_dterm_face(Er, Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
     dc, dc_l1, dc_l2, dc_l3, dc_h1, dc_h2, dc_h3, &
     dtf, dtf_l1, dtf_l2, dtf_l3, dtf_h1, dtf_h2, dtf_h3, dx, idir)
  implicit none
  
  integer, intent(in) :: Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
       dc_l1, dc_l2, dc_l3, dc_h1, dc_h2, dc_h3, &
       dtf_l1, dtf_l2, dtf_l3, dtf_h1, dtf_h2, dtf_h3, idir
  double precision, intent(in) :: dx(3)
  double precision, intent(in) :: Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3)
  double precision, intent(in) :: dc(dc_l1:dc_h1,dc_l2:dc_h2,dc_l3:dc_h3)
  double precision, intent(out) :: dtf(dtf_l1:dtf_h1,dtf_l2:dtf_h2,dtf_l3:dtf_h3)
  integer :: i, j, k

  if (idir .eq. 0) then
     do k = dtf_l3, dtf_h3
        do j = dtf_l2, dtf_h2
           do i = dtf_l1, dtf_h1
              dtf(i,j,k) = (Er(i,j,k) - Er(i-1,j,k)) / dx(1) * dc(i,j,k)
           end do
        end do
     end do
  else if (idir .eq. 1) then
     do k = dtf_l3, dtf_h3
        do j = dtf_l2, dtf_h2
           do i = dtf_l1, dtf_h1
              dtf(i,j,k) = (Er(i,j,k) - Er(i,j-1,k)) / dx(2) * dc(i,j,k)
           end do
        end do
     end do
  else
     do k = dtf_l3, dtf_h3
        do j = dtf_l2, dtf_h2
           do i = dtf_l1, dtf_h1
              dtf(i,j,k) = (Er(i,j,k) - Er(i,j,k-1)) / dx(2) * dc(i,j,k)
           end do
        end do
     end do
  end if

end subroutine ca_set_dterm_face


subroutine ca_face2center(  &
     foox, foox_l1, foox_l2, foox_l3, foox_h1, foox_h2, foox_h3, &
     fooy, fooy_l1, fooy_l2, fooy_l3, fooy_h1, fooy_h2, fooy_h3, &
     fooz, fooz_l1, fooz_l2, fooz_l3, fooz_h1, fooz_h2, fooz_h3, &
     fooc, fooc_l1, fooc_l2, fooc_l3, fooc_h1, fooc_h2, fooc_h3)

  implicit none

  integer, intent(in) :: foox_l1, foox_l2, foox_l3, foox_h1, foox_h2, foox_h3
  integer, intent(in) :: fooy_l1, fooy_l2, fooy_l3, fooy_h1, fooy_h2, fooy_h3
  integer, intent(in) :: fooz_l1, fooz_l2, fooz_l3, fooz_h1, fooz_h2, fooz_h3
  integer, intent(in) :: fooc_l1, fooc_l2, fooc_l3, fooc_h1, fooc_h2, fooc_h3
  double precision, intent(in)  :: foox(foox_l1:foox_h1,foox_l2:foox_h2,foox_l3:foox_h3)
  double precision, intent(in)  :: fooy(fooy_l1:fooy_h1,fooy_l2:fooy_h2,fooy_l3:fooy_h3)
  double precision, intent(in)  :: fooz(fooz_l1:fooz_h1,fooz_l2:fooz_h2,fooz_l3:fooz_h3)
  double precision, intent(out) :: fooc(fooc_l1:fooc_h1,fooc_l2:fooc_h2,fooc_l3:fooc_h3)

  integer :: i,j,k

  !$omp parallel do private(i,j,k)
  do k=fooc_l3,fooc_h3
     do j=fooc_l2,fooc_h2
        do i=fooc_l1,fooc_h1
           fooc(i,j,k) = (foox(i,j,k) + foox(i+1,j,k) &
                &       + fooy(i,j,k) + fooy(i,j+1,k) &
                &       + fooz(i,j,k) + fooz(i,j,k+1) )/6.d0
        end do
     end do
  end do
  !$omp end parallel do

end subroutine ca_face2center


subroutine ca_correct_dterm(  & 
     dfx, dfx_l1, dfx_l2, dfx_l3, dfx_h1, dfx_h2, dfx_h3, &
     dfy, dfy_l1, dfy_l2, dfy_l3, dfy_h1, dfy_h2, dfy_h3, &
     dfz, dfz_l1, dfz_l2, dfz_l3, dfz_h1, dfz_h2, dfz_h3, &
     re, rc)

  implicit none

  integer, intent(in) :: dfx_l1, dfx_l2, dfx_l3, dfx_h1, dfx_h2, dfx_h3
  integer, intent(in) :: dfy_l1, dfy_l2, dfy_l3, dfy_h1, dfy_h2, dfy_h3
  integer, intent(in) :: dfz_l1, dfz_l2, dfz_l3, dfz_h1, dfz_h2, dfz_h3
  double precision, intent(inout) :: dfx(dfx_l1:dfx_h1,dfx_l2:dfx_h2,dfx_l3:dfx_h3)
  double precision, intent(inout) :: dfy(dfy_l1:dfy_h1,dfy_l2:dfy_h2,dfy_l3:dfy_h3)
  double precision, intent(inout) :: dfz(dfz_l1:dfz_h1,dfz_l2:dfz_h2,dfz_l3:dfz_h3)
  double precision, intent(in) :: re(1), rc(1)

end subroutine ca_correct_dterm


subroutine ca_estdt_rad(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3, &
     gpr,gpr_l1,gpr_l2,gpr_l3,gpr_h1,gpr_h2,gpr_h3, &
     lo,hi,dx,dt)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, UFS, &
       UFX, allow_negative_energy
  
  implicit none
  
  integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
  integer          :: gpr_l1,gpr_l2,gpr_l3,gpr_h1,gpr_h2,gpr_h3
  integer          :: lo(3), hi(3)
  double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
  double precision :: gpr(gpr_l1:gpr_h1,gpr_l2:gpr_h2,gpr_l3:gpr_h3)
  double precision :: dx(3), dt
  
  double precision :: rhoInv,ux,uy,uz,dt1,dt2,dt3,c
  integer          :: i,j,k
  type(eos_t) :: eos_state

  ! Translate to primitive variables, compute sound speed (call eos)
  !$omp parallel do private(rhoInv,ux,uy,uz,dt1,dt2,dt3,i,j,k,eos_state,c) &
  !$omp reduction(min:dt)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           
           rhoInv = 1.d0/u(i,j,k,URHO)

           eos_state % rho = u(i,j,k,URHO)
           eos_state % T   = u(i,j,k,UTEMP)
           eos_state % e   = u(i,j,k,UEINT)*rhoInv
           eos_state % xn  = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
           eos_state % aux = u(i,j,k,UFX:UFX+naux -1) * rhoInv

           ! Protect against negative e
           if (eos_state%e .gt. 0.d0 .or. allow_negative_energy.eq.1) then
              call eos(eos_input_re, eos_state)
              c = eos_state % cs
           else
              c = 0.d0
           end if
           
           c = sqrt(c**2 + gpr(i,j,k)*rhoInv)

           ux = u(i,j,k,UMX)*rhoInv
           uy = u(i,j,k,UMY)*rhoInv
           uz = u(i,j,k,UMZ)*rhoInv

           dt1 = dx(1)/(c + abs(ux))
           dt2 = dx(2)/(c + abs(uy))
           dt3 = dx(3)/(c + abs(uz))
           dt = min(dt,dt1,dt2,dt3)
           
        enddo
     enddo
  enddo
  !$omp end parallel do
  
end subroutine ca_estdt_rad


subroutine ca_est_gpr0(Er, Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
     gPr, gPr_l1, gPr_l2, gPr_l3, gPr_h1, gPr_h2, gPr_h3)

  use rad_params_module, only : ngroups

  implicit none

  integer, intent(in) :: Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3
  integer, intent(in) :: gpr_l1, gpr_l2, gpr_l3, gpr_h1, gpr_h2, gpr_h3
  double precision, intent(in) :: Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3, 0:ngroups-1)
  double precision, intent(out) :: gPr(gPr_l1:gPr_h1,gPr_l2:gPr_h2,gPr_l3:gPr_h3)

  integer :: i, j, k, ig

  !$omp parallel private(i,j,k,ig)

  !$omp do
  do k = gPr_l3, gPr_h3
     do j = gPr_l2, gPr_h2
        do i = gPr_l1, gPr_h1
           gPr(i,j,k) = 0.d0
        end do
     end do
  end do
  !$omp end do nowait

  do ig = 0, ngroups-1
     !$omp do
     do k = gPr_l3, gPr_h3
        do j = gPr_l2, gPr_h2
           do i = gPr_l1, gPr_h1
              gPr(i,j,k) = gPr(i,j,k) + 4.d0/9.d0*Er(i,j,k,ig)
           end do
        end do
     end do
     !$omp end do nowait
  end do

  !$omp end parallel

end subroutine ca_est_gpr0


subroutine ca_est_gpr2(kap, kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
     Er, Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
     gPr, gPr_l1, gPr_l2, gPr_l3, gPr_h1, gPr_h2, gPr_h3, dx, limiter, comoving)

  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : FLDlambda, Edd_factor

  implicit none

  integer, intent(in) :: kap_l1, kap_l2, kap_l3, kap_h1, kap_h2, kap_h3, &
       Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
       gPr_l1, gPr_l2, gPr_l3, gPr_h1, gPr_h2, gPr_h3
  integer, intent(in) :: limiter, comoving
  double precision, intent(in) :: dx(3)
  double precision, intent(in) :: kap(kap_l1:kap_h1,kap_l2:kap_h2,kap_l3:kap_h3,0:ngroups-1), &
       Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3,0:ngroups-1)
  double precision, intent(out) :: gPr(gPr_l1:gPr_h1,gPr_l2:gPr_h2,gPr_l3:gPr_h3)

  integer :: i, j, k, g
  double precision, allocatable :: gE(:,:,:)
  double precision :: lam, gE1, gE2, gE3, r, f, gamr 

  allocate(gE(gPr_l1:gPr_h1,gPr_l2:gPr_h2,gPr_l3:gPr_h3))

  gPr = 0.d0

  do g = 0, ngroups-1

     do k = gPr_l3+1, gPr_h3-1
        do j = gPr_l2+1, gPr_h2-1
           do i = gPr_l1+1, gPr_h1-1
              gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
              gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
              gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
              gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
           end do
        end do
     end do
     
     ! lo-x lo-y lo-z
     i = gPr_l1
     j = gPr_l2
     k = gPr_l3
     gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
     gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
     gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     
     ! med-x lo-y lo-z
     j = gPr_l2
     k = gPr_l3
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! hi-x lo-y lo-z
     i = gPr_h1
     j = gPr_l2
     k = gPr_l3
     gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
     gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
     gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     
     ! lo-x med-y lo-z
     i = gPr_l1
     k = gPr_l3
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! med-x med-y lo-z side
     k = gPr_l3
     do j = gPr_l2+1, gPr_h2-1
        do i = gPr_l1+1, gPr_h1-1
           gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
           gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
           gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
           gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
        end do
     end do
     
     ! hi-x med-y lo-z
     i = gPr_h1
     k = gPr_l3
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! lo-x hi-y lo-z
     i = gPr_l1
     j = gPr_h2
     k = gPr_l3
     gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
     gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
     gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     
     ! med-x hi-y lo-z
     j = gPr_h2
     k = gPr_l3
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
        gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! hi-x hi-y lo-z
     i = gPr_h1
     j = gPr_h2
     k = gPr_l3
     gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
     gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
     gE3 = (Er(i,j,k+1,g) - Er(i,j,k,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)

     ! lo-x lo-y med-z
     i = gPr_l1
     j = gPr_l2
     do k = gPr_l3+1, gPr_h3-1
        gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! med-x lo-y med-z
     j = gPr_l2
     do k = gPr_l3+1, gPr_h3-1
        do i = gPr_l1+1, gPr_h1-1
           gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
           gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
           gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
           gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
        end do
     end do
     
     ! hi-x lo-y med-z
     i = gPr_h1
     j = gPr_l2
     do k = gPr_l3+1, gPr_h3-1
        gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! lo-x med-y med-z
     i = gPr_l1
     do k = gPr_l3+1, gPr_h3-1
        do j = gPr_l2+1, gPr_h2-1
           gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
           gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
           gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
           gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
        end do
     end do
     
     ! hi-x med-y med-z
     i = gPr_h1
     do k = gPr_l3+1, gPr_h3-1
        do j = gPr_l2+1, gPr_h2-1
           gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
           gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
           gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
           gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
        end do
     end do
     
     ! lo-x hi-y med-z
     i = gPr_l1
     j = gPr_h2
     do k = gPr_l3+1, gPr_h3-1
        gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
        gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! med-x hi-y med-z
     j = gPr_h2
     do k = gPr_l3+1, gPr_h3-1
        do i = gPr_l1+1, gPr_h1-1
           gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
           gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
           gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
           gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
        end do
     end do
     
     ! hi-x hi-y med-z
     i = gPr_h1
     j = gPr_h2
     do k = gPr_l3+1, gPr_h3-1
        gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
        gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
        gE3 = (Er(i,j,k+1,g) - Er(i,j,k-1,g)) / (2.d0*dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! lo-x lo-y hi-z
     i = gPr_l1
     j = gPr_l2
     k = gPr_h3
     gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
     gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
     gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)

     ! med-x lo-y hi-z
     j = gPr_l2
     k = gPr_h3
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
        gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! hi-x lo-y hi-z
     i = gPr_h1
     j = gPr_l2
     k = gPr_h3
     gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
     gE2 = (Er(i,j+1,k,g) - Er(i,j,k,g)) / (dx(2))
     gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     
     ! lo-x med-y hi-z
     i = gPr_l1
     k = gPr_h3
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
        gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! med-x med-y hi-z 
     k = gPr_h3
     do j = gPr_l2+1, gPr_h2-1
        do i = gPr_l1+1, gPr_h1-1
           gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
           gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
           gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
           gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
        end do
     end do
     
     ! hi-x med-y hi-z
     i = gPr_h1
     k = gPr_h3
     do j = gPr_l2+1, gPr_h2-1
        gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
        gE2 = (Er(i,j+1,k,g) - Er(i,j-1,k,g)) / (2.d0*dx(2))
        gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! lo-x hi-y hi-z
     i = gPr_l1
     j = gPr_h2
     k = gPr_h3
     gE1 = (Er(i+1,j,k,g) - Er(i,j,k,g)) / (dx(1))
     gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
     gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     
     ! med-x hi-y hi-z
     j = gPr_h2
     k = gPr_h3
     do i = gPr_l1+1, gPr_h1-1
        gE1 = (Er(i+1,j,k,g) - Er(i-1,j,k,g)) / (2.d0*dx(1))
        gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
        gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
        gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     end do
     
     ! hi-x hi-y hi-z
     i = gPr_h1
     j = gPr_h2
     k = gPr_h3
     gE1 = (Er(i,j,k,g) - Er(i-1,j,k,g)) / (dx(1))
     gE2 = (Er(i,j,k,g) - Er(i,j-1,k,g)) / (dx(2))
     gE3 = (Er(i,j,k,g) - Er(i,j,k-1,g)) / (dx(3))
     gE(i,j,k) = sqrt(gE1**2 + gE2**2 + gE3**2)
     
     do k = gPr_l3, gPr_h3
        do j = gPr_l2, gPr_h2
           do i = gPr_l1, gPr_h1
              r = gE(i,j,k) / (kap(i,j,k,g) * max(Er(i,j,k,g), 1.d-50))
              lam = FLDlambda(r, limiter)
              if (comoving .eq. 1) then
                 f = Edd_factor(lam)
                 gamr = (3.d0-f)/2.d0
              else
                 gamr = lam + 1.d0
              end if
              gPr(i,j,k) = gPr(i,j,k) +  lam * gamr * Er(i,j,k,g)
           end do
        end do
     end do
     
  end do

  deallocate(gE)

end subroutine ca_est_gpr2

