subroutine ca_umdrv(is_finest_level,time,&
     lo,hi,domlo,domhi,&
     uin,uin_l1,uin_h1,&
     uout,uout_l1,uout_h1,&
     update,updt_l1,updt_h1,&
     ugdnv,ugdnv_l1,ugdnv_h1,&
     src,src_l1,src_h1, &
     delta,dt,&
     flux,flux_l1,flux_h1,&
     area,area_l1,area_h1,&
     dloga,dloga_l1,dloga_h1,&
     vol,vol_l1,vol_h1,courno,verbose,&
     mass_added,eint_added,eden_added,frac_change,&
     mass_added_flux,xmom_added_flux,ymom_added_flux,zmom_added_flux,&
     E_added_flux,mass_lost,xmom_lost,ymom_lost,zmom_lost, &
     eden_lost,xang_lost,yang_lost,zang_lost) bind(C, name="ca_umdrv")

  use meth_params_module, only : QVAR, QU, QV, QW, QPRES, NVAR, NHYP, use_flattening
  use advection_module  , only : umeth1d, ctoprim, consup
  use bl_constants_module, only : ZERO, HALF, ONE
  use advection_util_module, only : compute_cfl
  use flatten_module, only : uflaten

  implicit none

  integer is_finest_level
  integer lo(1),hi(1),verbose
  integer domlo(1),domhi(1)
  integer uin_l1,uin_h1
  integer uout_l1,uout_h1
  integer updt_l1,updt_h1
  integer ugdnv_l1,ugdnv_h1
  integer flux_l1,flux_h1
  integer area_l1,area_h1
  integer dloga_l1,dloga_h1
  integer vol_l1,vol_h1
  integer src_l1,src_h1
  double precision   uin(  uin_l1:  uin_h1,NVAR)
  double precision  uout( uout_l1: uout_h1,NVAR)
  double precision update(updt_l1: updt_h1,NVAR)
  double precision ugdnv(ugdnv_l1:ugdnv_h1)
  double precision   src(  src_l1:  src_h1,NVAR)
  double precision  flux( flux_l1: flux_h1,NVAR)
  double precision  area( area_l1: area_h1     )
  double precision dloga(dloga_l1:dloga_h1     )
  double precision   vol(  vol_l1: vol_h1      )
  double precision delta(1),dt,time,courno

  !     Automatic arrays for workspace
  double precision, allocatable:: q(:,:)
  double precision, allocatable:: gamc(:)
  double precision, allocatable:: flatn(:)
  double precision, allocatable:: c(:)
  double precision, allocatable:: csml(:)
  double precision, allocatable:: div(:)
  double precision, allocatable:: pgdnv(:)
  double precision, allocatable:: srcQ(:,:)
  double precision, allocatable:: pdivu(:)

  double precision :: dx,E_added_flux,mass_added_flux
  double precision :: xmom_added_flux,ymom_added_flux,zmom_added_flux
  double precision :: mass_added, eint_added, eden_added, frac_change
  double precision :: mass_lost,xmom_lost,ymom_lost,zmom_lost
  double precision :: eden_lost,xang_lost,yang_lost,zang_lost
  integer i,ngf,ngq
  integer q_l1, q_h1

  integer ::  uin_lo(1),  uin_hi(1)
  integer :: uout_lo(1), uout_hi(1)

  integer :: lo_3D(3), hi_3D(3)
  integer :: q_lo_3D(3), q_hi_3D(3)
  double precision :: dx_3D(3)

  uin_lo  = [uin_l1]
  uin_hi  = [uin_h1]

  uout_lo = [uout_l1]
  uout_hi = [uout_h1]

  ngq = NHYP
  ngf = 1

  q_l1 = lo(1)-NHYP
  q_h1 = hi(1)+NHYP

  lo_3D   = [lo(1), 0, 0]
  hi_3D   = [hi(1), 0, 0]
  q_lo_3D = [q_l1, 0, 0]
  q_hi_3D = [q_h1, 0, 0]
  dx_3D   = [delta(1), ZERO, ZERO]

  allocate(     q(q_l1:q_h1,QVAR))
  allocate(     c(q_l1:q_h1))
  allocate(  gamc(q_l1:q_h1))
  allocate( flatn(q_l1:q_h1))
  allocate(  csml(q_l1:q_h1))

  allocate(  srcQ(q_l1:q_h1,QVAR))

  allocate(   div(lo(1):hi(1)+1))
  allocate( pdivu(lo(1):hi(1)  ))
  allocate( pgdnv(lo(1):hi(1)+1))

  dx = delta(1)

  !     Translate to primitive variables, compute sound speeds
  !     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
  !       and set to correspond to coordinates of (lo:hi)

  call ctoprim(lo,hi,uin,uin_l1,uin_h1, &
       q,c,gamc,csml,q_l1,q_h1, &
       src,src_l1,src_h1, &
       srcQ,q_l1,q_h1, &
       dx,dt,ngq)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo_3D, q_hi_3D, c, q_lo_3D, q_hi_3D, lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coef for slope calculations
  if (use_flattening == 1) then
     call uflaten((/ lo(1) - ngf, 0, 0 /), (/ hi(1) + ngf, 0, 0 /), &
                  q(:,QPRES), q(:,QU), q(:,QV), q(:,QW), &
                  flatn,(/ q_l1, 0, 0 /), (/ q_h1, 0, 0 /))
  else
     flatn = ONE
  endif

  call umeth1d(lo,hi,domlo,domhi, &
       q,c,gamc,csml,flatn,q_l1,q_h1, &
       srcQ, q_l1,q_h1, &
       lo(1),hi(1),dx,dt, &
       flux,flux_l1,flux_h1, &
       pgdnv,lo(1),hi(1)+1, &
       ugdnv,ugdnv_l1,ugdnv_h1, &
       dloga,dloga_l1,dloga_h1)

  ! Define p*divu
  do i = lo(1), hi(1)
     pdivu(i) = HALF * &
          (pgdnv(i+1)+pgdnv(i))*(ugdnv(i+1)*area(i+1)-ugdnv(i)*area(i)) / vol(i)
  end do

  ! Define divu on surroundingNodes(lo,hi)
  do i = lo(1),hi(1)+1
     div(i) = (q(i,QU)-q(i-1,QU)) / dx
  enddo

  !     Conservative update
  call consup(uin,uin_l1,uin_h1, &
       uout,uout_l1,uout_h1, &
       update,updt_l1,updt_h1, &
       pgdnv,lo(1),hi(1)+1, &
       src , src_l1, src_h1, &
       flux,flux_l1,flux_h1, &
       area,area_l1,area_h1, &
       vol , vol_l1, vol_h1, &
       div ,pdivu,lo,hi,dx,dt,mass_added_flux,E_added_flux, &
       xmom_added_flux,ymom_added_flux,zmom_added_flux, &
       mass_lost,xmom_lost,ymom_lost,zmom_lost, &
       eden_lost,xang_lost,yang_lost,zang_lost, &
       verbose)

  deallocate(q,c,gamc,flatn,csml,srcQ,div,pdivu,pgdnv)

end subroutine ca_umdrv
