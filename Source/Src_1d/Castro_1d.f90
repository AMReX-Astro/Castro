subroutine ca_umdrv(is_finest_level,time,&
     lo,hi,domlo,domhi,&
     uin,uin_l1,uin_h1,&
     uout,uout_l1,uout_h1,&
     q,q_l1,q_h1,&
     srcQ,srQ_l1,srQ_h1,&
     update,updt_l1,updt_h1,&
     ugdnv,ugdnv_l1,ugdnv_h1,&
     delta,dt,&
     flux,flux_l1,flux_h1,&
     area,area_l1,area_h1,&
     dloga,dloga_l1,dloga_h1,&
     vol,vol_l1,vol_h1,courno,verbose,&
     mass_added_flux,xmom_added_flux,ymom_added_flux,zmom_added_flux,&
     E_added_flux,mass_lost,xmom_lost,ymom_lost,zmom_lost, &
     eden_lost,xang_lost,yang_lost,zang_lost) bind(C, name="ca_umdrv")

  use meth_params_module, only : QVAR, QU, QV, QW, QPRES, NVAR, NHYP, use_flattening
  use advection_module  , only : umeth1d, consup
  use bl_constants_module, only : ZERO, HALF, ONE
  use advection_util_module, only : compute_cfl
  use flatten_module, only : uflaten

  implicit none

  integer is_finest_level
  integer lo(1),hi(1),verbose
  integer domlo(1),domhi(1)
  integer uin_l1,uin_h1
  integer uout_l1,uout_h1
  integer q_l1,q_h1
  integer srQ_l1,srQ_h1
  integer updt_l1,updt_h1
  integer ugdnv_l1,ugdnv_h1
  integer flux_l1,flux_h1
  integer area_l1,area_h1
  integer dloga_l1,dloga_h1
  integer vol_l1,vol_h1
  double precision   uin(  uin_l1:  uin_h1,NVAR)
  double precision  uout( uout_l1: uout_h1,NVAR)
  double precision     q(    q_l1:    q_h1,QVAR)
  double precision  srcQ(  srQ_l1:  srQ_h1,QVAR)
  double precision update(updt_l1: updt_h1,NVAR)
  double precision ugdnv(ugdnv_l1:ugdnv_h1)
  double precision  flux( flux_l1: flux_h1,NVAR)
  double precision  area( area_l1: area_h1     )
  double precision dloga(dloga_l1:dloga_h1     )
  double precision   vol(  vol_l1: vol_h1      )
  double precision delta(1),dt,time,courno

  !     Automatic arrays for workspace
  double precision, allocatable:: flatn(:)
  double precision, allocatable:: div(:)
  double precision, allocatable:: pgdnv(:)
  double precision, allocatable:: pdivu(:)

  double precision :: dx,E_added_flux,mass_added_flux
  double precision :: xmom_added_flux,ymom_added_flux,zmom_added_flux
  double precision :: mass_lost,xmom_lost,ymom_lost,zmom_lost
  double precision :: eden_lost,xang_lost,yang_lost,zang_lost
  integer i,ngf,ngq

  integer ::  uin_lo(1),  uin_hi(1)
  integer :: uout_lo(1), uout_hi(1)
  integer :: q_lo(1), q_hi(1)

  integer :: lo_3D(3), hi_3D(3)
  integer :: q_lo_3D(3), q_hi_3D(3)
  integer :: uin_lo_3D(3), uin_hi_3D(3)
  double precision :: dx_3D(3)

  uin_lo  = [uin_l1]
  uin_hi  = [uin_h1]

  uout_lo = [uout_l1]
  uout_hi = [uout_h1]

  ngq = NHYP
  ngf = 1

  q_lo = [q_l1]
  q_hi = [q_h1]

  lo_3D   = [lo(1), 0, 0]
  hi_3D   = [hi(1), 0, 0]
  q_lo_3D = [q_l1, 0, 0]
  q_hi_3D = [q_h1, 0, 0]
  uin_lo_3D = [uin_l1, 0, 0]
  uin_hi_3D = [uin_h1, 0, 0]
  dx_3D   = [delta(1), ZERO, ZERO]

  allocate( flatn(q_l1:q_h1))

  allocate(   div(lo(1):hi(1)+1))
  allocate( pdivu(lo(1):hi(1)  ))
  allocate( pgdnv(lo(1):hi(1)+1))

  dx = delta(1)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo_3D, q_hi_3D, lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (use_flattening == 1) then
     call uflaten((/ lo(1) - ngf, 0, 0 /), (/ hi(1) + ngf, 0, 0 /), &
                  q(:,QPRES), q(:,QU), q(:,QV), q(:,QW), &
                  flatn,(/ q_l1, 0, 0 /), (/ q_h1, 0, 0 /))
  else
     flatn = ONE
  endif

  call umeth1d(lo,hi,domlo,domhi, &
       q,flatn,q_l1,q_h1, &
       srcQ,srQ_l1,srQ_h1, &
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
       flux,flux_l1,flux_h1, &
       area,area_l1,area_h1, &
       vol , vol_l1, vol_h1, &
       div ,pdivu,lo,hi,dx,dt,mass_added_flux,E_added_flux, &
       xmom_added_flux,ymom_added_flux,zmom_added_flux, &
       mass_lost,xmom_lost,ymom_lost,zmom_lost, &
       eden_lost,xang_lost,yang_lost,zang_lost, &
       verbose)

  deallocate(flatn,div,pdivu,pgdnv)

end subroutine ca_umdrv
