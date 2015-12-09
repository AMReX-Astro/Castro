subroutine ca_umdrv(is_finest_level,time,lo,hi,domlo,domhi, &
                    uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                    uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                    ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                    ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                    ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                    src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                    delta,dt, &
                    flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                    flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                    flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                    area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                    area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                    area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                    vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                    courno,verbose,mass_added,eint_added,eden_added,&
                    xmom_added_flux,ymom_added_flux,zmom_added_flux,&
                    E_added_flux)

  use mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : QVAR, NVAR, NHYP, &
                                 normalize_species
  use advection_module, only : umeth3d, ctoprim, consup
  use advection_util_module, only : divu
  use advection_util_module, only : enforce_minimum_density, normalize_new_species

  implicit none

  integer is_finest_level
  integer lo(3),hi(3),verbose
  integer domlo(3),domhi(3)
  integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
  integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
  integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
  integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
  integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
  integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
  integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
  integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
  integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
  integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
  integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
  integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
  integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
  double precision   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,     uin_l3:uin_h3,  NVAR)
  double precision  uout( uout_l1:uout_h1,  uout_l2:uout_h2,   uout_l3:uout_h3, NVAR)
  double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
  double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
  double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
  double precision   src(  src_l1:src_h1,    src_l2:src_h2,     src_l3:src_h3,  NVAR)
  double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
  double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
  double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2, area1_l3:area1_h3)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2, area2_l3:area2_h3)
  double precision area3(area3_l1:area3_h1,area3_l2:area3_h2, area3_l3:area3_h3)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2, vol_l3:vol_h3)
  double precision delta(3),dt,time,courno,E_added_flux
  double precision mass_added,eint_added,eden_added
  double precision xmom_added_flux,ymom_added_flux,zmom_added_flux

  ! Automatic arrays for workspace
  double precision, pointer:: q(:,:,:,:)
  double precision, pointer:: gamc(:,:,:)
  double precision, pointer:: flatn(:,:,:)
  double precision, pointer:: c(:,:,:)
  double precision, pointer:: csml(:,:,:)
  double precision, pointer:: div(:,:,:)
  double precision, pointer:: pdivu(:,:,:)
  double precision, pointer:: srcQ(:,:,:,:)
  
  double precision dx,dy,dz
  integer ngq,ngf
  integer q_l1, q_l2, q_l3, q_h1, q_h2, q_h3

  ngq = NHYP
  ngf = 1
    
  q_l1 = lo(1)-NHYP
  q_l2 = lo(2)-NHYP
  q_l3 = lo(3)-NHYP
  q_h1 = hi(1)+NHYP
  q_h2 = hi(2)+NHYP
  q_h3 = hi(3)+NHYP

  call bl_allocate(     q, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3,1,QVAR)
  call bl_allocate(  gamc, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
  call bl_allocate( flatn, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
  call bl_allocate(     c, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
  call bl_allocate(  csml, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)

  call bl_allocate(   div, lo(1),hi(1)+1,lo(2),hi(2)+1,lo(3),hi(3)+1)  
  call bl_allocate( pdivu, lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)  )

  call bl_allocate(  srcQ, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3,1,QVAR)
  
  dx = delta(1)
  dy = delta(2)
  dz = delta(3)

  ! 1) Translate conserved variables (u) to primitive variables (q).
  ! 2) Compute sound speeds (c) and gamma (gamc).
  !    Note that (q,c,gamc,csml,flatn) are all dimensioned the same
  !    and set to correspond to coordinates of (lo:hi)
  ! 3) Translate source terms
  call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
               q,c,gamc,csml,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
               src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
               srcQ,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
               courno,dx,dy,dz,dt,ngq,ngf)

  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth3d(q,c,gamc,csml,flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
               srcQ,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
               lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
               flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
               flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
               flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
               ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
               ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
               ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
               pdivu, domlo, domhi)

  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo,hi,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
            dx,dy,dz,div,lo(1),lo(2),lo(3),hi(1)+1,hi(2)+1,hi(3)+1)

  ! Conservative update
  call consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
              uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
              src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
              flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
              flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
              flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
              area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
              area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
              area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
              vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
              div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux,&
              xmom_added_flux,ymom_added_flux,zmom_added_flux)

  ! Add the radiative cooling -- for SGS only.
  ! if (radiative_cooling_type.eq.2) then
  !    call post_step_radiative_cooling(lo,hi,dt, &
  !         uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
  ! endif

  ! Enforce the density >= small_dens.
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

  call bl_deallocate(     q)
  call bl_deallocate(  gamc)
  call bl_deallocate( flatn)
  call bl_deallocate(     c)
  call bl_deallocate(  csml)

  call bl_deallocate(   div)
  call bl_deallocate( pdivu)

  call bl_deallocate(  srcQ)

end subroutine ca_umdrv
