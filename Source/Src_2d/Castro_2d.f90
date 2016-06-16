subroutine ca_umdrv(is_finest_level,time,lo,hi,domlo,domhi, &
                    uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                    uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                    ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                    ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                    src,src_l1,src_l2,src_h1,src_h2, &
                    delta,dt, &
                    flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                    flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                    area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                    area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                    dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                    vol,vol_l1,vol_l2,vol_h1,vol_h2,&
                    courno,verbose,mass_added,eint_added,eden_added,frac_change, &
                    mass_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux, &
                    E_added_flux,mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                    eden_lost,xang_lost,yang_lost,zang_lost) bind(C, name="ca_umdrv")

  use meth_params_module, only : QVAR, NVAR, NHYP, ngdnv, GDU, GDV
  use advection_module, only : umeth2d, ctoprim, consup
  use advection_util_2d_module, only : divu
  use castro_util_module, only : ca_normalize_species
  use advection_util_module, only : enforce_minimum_density

  implicit none

  integer is_finest_level
  integer lo(2),hi(2),verbose
  integer domlo(2),domhi(2)
  integer uin_l1,uin_l2,uin_h1,uin_h2
  integer uout_l1,uout_l2,uout_h1,uout_h2
  integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
  integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
  integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
  integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
  integer area1_l1,area1_l2,area1_h1,area1_h2
  integer area2_l1,area2_l2,area2_h1,area2_h2
  integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
  integer vol_l1,vol_l2,vol_h1,vol_h2
  integer src_l1,src_l2,src_h1,src_h2

  double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
  double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
  double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
  double precision src(src_l1:src_h1,src_l2:src_h2,NVAR)
  double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
  double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
  double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
  double precision delta(2),dt,time,courno
  double precision E_added_flux, mass_added_flux
  double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
  double precision mass_added,eint_added,eden_added,frac_change
  double precision mass_lost,xmom_lost,ymom_lost,zmom_lost
  double precision eden_lost,xang_lost,yang_lost,zang_lost

  ! Automatic arrays for workspace
  double precision, allocatable :: q(:,:,:)
  double precision, allocatable :: gamc(:,:)
  double precision, allocatable :: flatn(:,:)
  double precision, allocatable :: c(:,:)
  double precision, allocatable :: csml(:,:)
  double precision, allocatable :: div(:,:)
  double precision, allocatable :: srcQ(:,:,:)
  double precision, allocatable :: pdivu(:,:)

  ! Edge-centered primitive variables (Riemann state)
  double precision, allocatable :: q1(:,:,:)
  double precision, allocatable :: q2(:,:,:)

  integer ngq,ngf
  double precision dx,dy

  integer q_l1, q_l2, q_h1, q_h2

  integer ::  uin_lo(2),  uin_hi(2)
  integer :: uout_lo(2), uout_hi(2)

  uin_lo  = [uin_l1, uin_l2]
  uin_hi  = [uin_h1, uin_h2]

  uout_lo = [uout_l1, uout_l2]
  uout_hi = [uout_h1, uout_h2]
  
  ngq = NHYP
  ngf = 1

  q_l1 = lo(1)-NHYP
  q_l2 = lo(2)-NHYP
  q_h1 = hi(1)+NHYP
  q_h2 = hi(2)+NHYP

  allocate(     q(q_l1:q_h1,q_l2:q_h2,QVAR))
  allocate(  gamc(q_l1:q_h1,q_l2:q_h2))
  allocate( flatn(q_l1:q_h1,q_l2:q_h2))
  allocate(     c(q_l1:q_h1,q_l2:q_h2))
  allocate(  csml(q_l1:q_h1,q_l2:q_h2))

  allocate(  srcQ(q_l1:q_h1,q_l2:q_h2,QVAR))

  allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))
  allocate( pdivu(lo(1)  :hi(1)  ,lo(2)  :hi(2)))

  allocate(q1(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ngdnv))
  allocate(q2(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ngdnv))

  dx = delta(1)
  dy = delta(2)

  ! Translate to primitive variables, compute sound speeds.  Note that
  ! (q,c,gamc,csml,flatn) are all dimensioned the same and set to
  ! correspond to coordinates of (lo:hi)
  call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2, &
               q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
               src,src_l1,src_l2,src_h1,src_h2, &
               srcQ,q_l1,q_l2,q_h1,q_h2, &
               courno,dx,dy,dt,ngq,ngf)
  
  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth2d(q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
               srcQ, q_l1,q_l2,q_h1,q_h2, &
               lo(1),lo(2),hi(1),hi(2),dx,dy,dt, &
               uout,uout_l1,uout_l2,uout_h1,uout_h2, &
               flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
               flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
               q1, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
               q2, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
               area1, area1_l1, area1_l2, area1_h1, area1_h2, &
               area2, area2_l1, area2_l2, area2_h1, area2_h2, &
               pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
               dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
               domlo, domhi)

  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  ! this is used for the artifical viscosity
  call divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2, &
            delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

  !     Conservative update
  call consup(uin,    uin_l1,  uin_l2,  uin_h1,  uin_h2, &
              uout,  uout_l1, uout_l2, uout_h1, uout_h2, &
              q1, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
              q2, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
              src,    src_l1,  src_l2,  src_h1,  src_h2, &
              flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
              flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
              area1,area1_l1,area1_l2,area1_h1,area1_h2, &
              area2,area2_l1,area2_l2,area2_h1,area2_h2, &
              vol,    vol_l1,  vol_l2,  vol_h1,  vol_h2, &
              div,pdivu,lo,hi,dx,dy,dt,mass_added_flux,E_added_flux, &
              xmom_added_flux,ymom_added_flux,zmom_added_flux, &
              mass_lost,xmom_lost,ymom_lost,zmom_lost, &
              eden_lost,xang_lost,yang_lost,zang_lost, &
              verbose)

  ! Enforce the density >= small_dens.
  call enforce_minimum_density(uin,[uin_lo(1), uin_lo(2), 0],[uin_hi(1), uin_hi(2), 0], &
                               uout,[uout_lo(1), uout_lo(2), 0],[uout_hi(1),uout_hi(2),0], &
                               [lo(1), lo(2), 0],[hi(1), hi(2), 0], &
                               mass_added,eint_added,eden_added,frac_change,verbose)

  ! Renormalize species mass fractions
  call ca_normalize_species(uout, &
                            [uout_lo(1), uout_lo(2), 0], [uout_hi(1), uout_hi(2), 0], &
                            [lo(1), lo(2), 0], [hi(1), hi(2), 0])

  ugdx(:,:) = q1(:,:,GDU)
  ugdy(:,:) = q2(:,:,GDV)

  deallocate(q,gamc,flatn,c,csml,div,q1,q2,srcQ,pdivu)

end subroutine ca_umdrv
