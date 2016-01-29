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
                    E_added_flux) bind(C)

  use mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : QVAR, NVAR, NHYP, NGDNV, &
                                 normalize_species, GDU, GDV, GDW
  use advection_module, only : umeth3d, ctoprim, consup
  use advection_util_module, only : divu, enforce_minimum_density, normalize_new_species
  use castro_util_3d_module, only : ca_enforce_nonnegative_species

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

  ! Edge-centered primitive variables (Riemann state)
  double precision, pointer:: q1(:,:,:,:)
  double precision, pointer:: q2(:,:,:,:)
  double precision, pointer:: q3(:,:,:,:)
  
  integer :: ngq, ngf
  integer :: uin_lo(3), uin_hi(3)
  integer :: uout_lo(3), uout_hi(3)
  integer :: flux1_lo(3), flux1_hi(3)
  integer :: flux2_lo(3), flux2_hi(3)
  integer :: flux3_lo(3), flux3_hi(3)  
  integer :: area1_lo(3), area1_hi(3)
  integer :: area2_lo(3), area2_hi(3)
  integer :: area3_lo(3), area3_hi(3)
  integer :: ugdnvx_lo(3), ugdnvx_hi(3)
  integer :: ugdnvy_lo(3), ugdnvy_hi(3)
  integer :: ugdnvz_lo(3), ugdnvz_hi(3)
  integer :: vol_lo(3), vol_hi(3)
  integer :: src_lo(3), src_hi(3)
  integer :: q_lo(3), q_hi(3)
  integer :: q1_lo(3), q1_hi(3), q2_lo(3), q2_hi(3), q3_lo(3), q3_hi(3)

  ngq = NHYP
  ngf = 1
    
  q_lo(:) = lo(:) - NHYP
  q_hi(:) = hi(:) + NHYP

  uin_lo = [ uin_l1, uin_l2, uin_l3 ]
  uin_hi = [ uin_h1, uin_h2, uin_h3 ]
  
  uout_lo = [ uout_l1, uout_l2, uout_l3 ]
  uout_hi = [ uout_h1, uout_h2, uout_h3 ]

  flux1_lo = [ flux1_l1, flux1_l2, flux1_l3 ]
  flux1_hi = [ flux1_h1, flux1_h2, flux1_h3 ]

  flux2_lo = [ flux2_l1, flux2_l2, flux2_l3 ]
  flux2_hi = [ flux2_h1, flux2_h2, flux2_h3 ]

  flux3_lo = [ flux3_l1, flux3_l2, flux3_l3 ]
  flux3_hi = [ flux3_h1, flux3_h2, flux3_h3 ]  
  
  area1_lo = [ area1_l1, area1_l2, area1_l3 ]
  area1_hi = [ area1_h1, area1_h2, area1_h3 ]

  area2_lo = [ area2_l1, area2_l2, area2_l3 ]
  area2_hi = [ area2_h1, area2_h2, area2_h3 ]

  area3_lo = [ area3_l1, area3_l2, area3_l3 ]
  area3_hi = [ area3_h1, area3_h2, area3_h3 ]

  vol_lo = [ vol_l1, vol_l2, vol_l3 ]
  vol_hi = [ vol_h1, vol_h2, vol_h3 ]

  src_lo = [ src_l1, src_l2, src_l3 ]
  src_hi = [ src_h1, src_h2, src_h3 ]  
  
  call bl_allocate(     q, q_lo, q_hi, QVAR)
  call bl_allocate(  gamc, q_lo, q_hi)
  call bl_allocate( flatn, q_lo, q_hi)
  call bl_allocate(     c, q_lo, q_hi)
  call bl_allocate(  csml, q_lo, q_hi)

  call bl_allocate(   div, lo(1), hi(1)+1, lo(2), hi(2)+1, lo(3), hi(3)+1)  
  call bl_allocate( pdivu, lo(1), hi(1)  , lo(2), hi(2)  , lo(3), hi(3)  )

  call bl_allocate(  srcQ, q_lo, q_hi, QVAR)

  ugdnvx_lo = [ugdnvx_l1, ugdnvx_l2, ugdnvx_l3]
  ugdnvx_hi = [ugdnvx_h1, ugdnvx_h2, ugdnvx_h3]
  ugdnvy_lo = [ugdnvy_l1, ugdnvy_l2, ugdnvy_l3]
  ugdnvy_hi = [ugdnvy_h1, ugdnvy_h2, ugdnvy_h3]
  ugdnvz_lo = [ugdnvz_l1, ugdnvz_l2, ugdnvz_l3]
  ugdnvz_hi = [ugdnvz_h1, ugdnvz_h2, ugdnvz_h3]

  q1_lo = ugdnvx_lo
  q1_hi = ugdnvx_hi
  q2_lo = ugdnvy_lo
  q2_hi = ugdnvy_hi
  q3_lo = ugdnvz_lo
  q3_hi = ugdnvz_hi
  
  call bl_allocate(q1, q1_lo, q1_hi, NGDNV)
  call bl_allocate(q2, q2_lo, q2_hi, NGDNV)
  call bl_allocate(q3, q3_lo, q3_hi, NGDNV)
  
  ! 1) Translate conserved variables (u) to primitive variables (q).
  ! 2) Compute sound speeds (c) and gamma (gamc).
  !    Note that (q,c,gamc,csml,flatn) are all dimensioned the same
  !    and set to correspond to coordinates of (lo:hi)
  ! 3) Translate source terms
  call ctoprim(lo,hi,uin,uin_lo,uin_hi, &
               q,c,gamc,csml,flatn,q_lo,q_hi, &
               src,src_lo,src_hi, &
               srcQ,q_lo,q_hi, &
               courno,delta,dt,ngq,ngf)

  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth3d(q,c,gamc,csml,flatn,q_lo,q_hi, &
               srcQ,q_lo,q_hi, &
               lo,hi,delta,dt, &
               flux1,flux1_lo,flux1_hi, &
               flux2,flux2_lo,flux2_hi, &
               flux3,flux3_lo,flux3_hi, &
               q1, q1_lo, q1_hi, &
               q2, q2_lo, q2_hi, &
               q3, q3_lo, q3_hi, &
               pdivu, domlo, domhi)

  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo,hi,q,q_lo,q_hi,delta,div,lo,hi+1)

  ! Conservative update
  call consup(uin ,  uin_lo , uin_hi, &
              uout, uout_lo, uout_hi, &
              src ,  src_lo,  src_hi, &
              flux1, flux1_lo, flux1_hi, &
              flux2, flux2_lo, flux2_hi, &
              flux3, flux3_lo, flux3_hi, &
              q1, q1_lo, q1_hi, &
              q2, q2_lo, q2_hi, &
              q3, q3_lo, q3_hi, &
              area1, area1_lo, area1_hi, &
              area2, area2_lo, area2_hi, &
              area3, area3_lo, area3_hi, &
              vol, vol_lo, vol_hi, &
              div,pdivu,lo,hi,delta,dt,E_added_flux, &
              xmom_added_flux,ymom_added_flux,zmom_added_flux, &
              verbose)

  ! Add the radiative cooling -- for SGS only.
  ! if (radiative_cooling_type.eq.2) then
  !    call post_step_radiative_cooling(lo,hi,dt, &
  !         uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
  ! endif

  ! Enforce the density >= small_dens.
  call enforce_minimum_density(uin,uin_lo,uin_hi,uout,uout_lo,uout_hi, &
                               lo,hi,mass_added,eint_added,eden_added,verbose)

  ! Enforce species >= 0
  call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                      uout_h1,uout_h2,uout_h3,lo,hi)
 
  ! Re-normalize the species
  if (normalize_species .eq. 1) then
     call normalize_new_species(uout,uout_lo,uout_hi,lo,hi)
  end if

  ! Copy data from the edge-centered state into ugdnv

  ugdnvx_out(:,:,:) = q1(:,:,:,GDU)
  ugdnvy_out(:,:,:) = q2(:,:,:,GDV)
  ugdnvz_out(:,:,:) = q3(:,:,:,GDW)
  
  call bl_deallocate(     q)
  call bl_deallocate(  gamc)
  call bl_deallocate( flatn)
  call bl_deallocate(     c)
  call bl_deallocate(  csml)

  call bl_deallocate(   div)
  call bl_deallocate( pdivu)

  call bl_deallocate(  srcQ)

  call bl_deallocate(    q1)
  call bl_deallocate(    q2)
  call bl_deallocate(    q3)
  
end subroutine ca_umdrv
