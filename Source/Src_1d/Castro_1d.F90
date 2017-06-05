subroutine ca_ctu_update(is_finest_level, time, &
                    lo, hi, domlo, domhi, &
                    uin, uin_l1, uin_h1, &
                    uout, uout_l1, uout_h1, &
#ifdef RADIATION
                    Erin, Erin_l1, Erin_h1, &
                    Erout, Erout_l1, Erout_h1, &
#endif
                    q, q_l1, q_h1, &
                    qaux, qa_l1, qa_h1, &
                    srcQ, srQ_l1, srQ_h1, &
                    update, updt_l1, updt_h1, &
                    delta, dt, &
                    flux, flux_l1, flux_h1, &
#ifdef RADIATION
                    radflux, radflux_l1, radflux_h1,&
#endif
                    pradial, p_l1, p_h1, &
                    area, area_l1, area_h1, &
                    dloga, dloga_l1, dloga_h1, &
                    vol, vol_l1, vol_h1, &
                    courno, verbose, &
#ifdef RADIATION
                    nstep_fsp, &
#endif
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

  use meth_params_module, only : NQ, QVAR, QU, QV, QW, QPRES, &
#ifdef RADIATION
                                 QPTOT, &
#endif
                                 NQAUX, NVAR, NHYP, use_flattening, &
                                 NGDNV, GDU, GDPRES, first_order_hydro
  use bl_constants_module, only : ZERO, HALF, ONE
  use advection_util_module, only : compute_cfl
  use flatten_module, only : uflaten
  use prob_params_module, only : coord_type
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use flatten_module, only : rad_flaten
#endif
  use ctu_advection_module  , only : umeth1d, consup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: is_finest_level
  integer, intent(in) :: lo(1), hi(1), verbose
  integer, intent(in) :: domlo(1), domhi(1)
  integer, intent(in) :: uin_l1, uin_h1
  integer, intent(in) :: uout_l1, uout_h1
  integer, intent(in) :: q_l1, q_h1
  integer, intent(in) :: qa_l1, qa_h1
  integer, intent(in) :: srQ_l1,srQ_h1
  integer, intent(in) :: updt_l1, updt_h1
  integer, intent(in) :: flux_l1, flux_h1
#ifdef RADIATION
  integer, intent(inout) :: nstep_fsp
  integer, intent(in) :: Erin_l1, Erin_h1
  integer, intent(in) :: Erout_l1, Erout_h1
  integer, intent(in) :: radflux_l1, radflux_h1
#endif
  integer, intent(in) :: p_l1, p_h1
  integer, intent(in) :: area_l1, area_h1
  integer, intent(in) :: dloga_l1, dloga_h1
  integer, intent(in) :: vol_l1, vol_h1

  real(rt)        , intent(in) ::      uin(  uin_l1:  uin_h1,NVAR)
  real(rt)        , intent(inout) ::  uout( uout_l1: uout_h1,NVAR)
  real(rt)        , intent(inout) ::     q(    q_l1:    q_h1,NQ)
  real(rt)        , intent(in) ::     qaux(   qa_l1:   qa_h1,NQAUX)
  real(rt)        , intent(in) ::     srcQ(  srQ_l1:  srQ_h1,QVAR)
  real(rt)        , intent(inout) :: update(updt_l1: updt_h1,NVAR)
  real(rt)        , intent(inout) ::  flux( flux_l1: flux_h1,NVAR)
#ifdef RADIATION
  real(rt)        , intent(inout) :: Erout(Erout_l1:Erout_h1, 0:ngroups-1)
  real(rt)        , intent(inout) :: radflux(radflux_l1: radflux_h1, 0:ngroups-1)
  real(rt)        , intent(in) :: Erin( Erin_l1: Erin_h1, 0:ngroups-1)
#endif
  real(rt)        , intent(inout) :: pradial(  p_l1:   p_h1)
  real(rt)        , intent(in) :: area( area_l1: area_h1     )
  real(rt)        , intent(in) :: dloga(dloga_l1:dloga_h1     )
  real(rt)        , intent(in) ::   vol(  vol_l1: vol_h1      )
  real(rt)        , intent(in) :: delta(1), dt, time
  real(rt)        , intent(inout) :: courno

  real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
  real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

  ! Automatic arrays for workspace
  real(rt)        , allocatable:: flatn(:)
  real(rt)        , allocatable:: div(:)
  real(rt)        , allocatable:: pdivu(:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , allocatable :: q1(:,:)

  real(rt)         :: dx

  integer i,ngf, ngq

  integer ::  uin_lo(1),  uin_hi(1)
  integer :: uout_lo(1), uout_hi(1)
  integer :: q_lo(1), q_hi(1)

  integer :: lo_3D(3), hi_3D(3)
  integer :: q_lo_3D(3), q_hi_3D(3)
  integer :: uin_lo_3D(3), uin_hi_3D(3)
  real(rt)         :: dx_3D(3)

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
  allocate(    q1(flux_l1-1:flux_h1+1,NGDNV))

  dx = delta(1)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo_3D, q_hi_3D, &
                   qaux, [qa_l1, 0, 0], [qa_h1, 0, 0], &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  else if (use_flattening == 1) then
#ifdef RADIATION
     call rad_flaten([lo(1)-ngf, 0, 0], [hi(1)+ngf, 0, 0], &
                      q(:,qpres), q(:,qptot), &
                      q(:,QU), q(:,QV), q(:,QW), &
                      flatn, [q_l1, 0, 0], [q_h1, 0, 0])
#else
     call uflaten([lo(1) - ngf, 0, 0], [hi(1) + ngf, 0, 0], &
                  q(:,QPRES), q(:,QU), q(:,QV), q(:,QW), &
                  flatn, [q_l1, 0, 0], [q_h1, 0, 0])
#endif
  else
     flatn = ONE
  endif

  call umeth1d(lo, hi, domlo, domhi, &
               q, q_l1, q_h1, &
               flatn, &
               qaux, qa_l1, qa_h1, &
               srcQ, srQ_l1, srQ_h1, &
               lo(1), hi(1), dx, dt, &
               uout, uout_l1, uout_h1, &
               flux, flux_l1, flux_h1, &
#ifdef RADIATION
               radflux,radflux_l1,radflux_h1, &
#endif
               q1, flux_l1-1, flux_h1+1, &
               dloga, dloga_l1, dloga_h1)


  ! Define p*divu
  do i = lo(1), hi(1)
     pdivu(i) = HALF * &
          (q1(i+1,GDPRES) + q1(i,GDPRES))* &
          (q1(i+1,GDU)*area(i+1) - q1(i,GDU)*area(i)) / vol(i)
  end do

  ! Define divu on surroundingNodes(lo,hi)
  do i = lo(1), hi(1)+1
     div(i) = (q(i,QU) - q(i-1,QU)) / dx
  enddo


  ! Conservative update
  call consup(uin, uin_l1, uin_h1, &
              uout, uout_l1, uout_h1, &
              update, updt_l1, updt_h1, &
              q, q_l1, q_h1, &
              flux, flux_l1, flux_h1, &
              q1, flux_l1-1, flux_h1+1, &
#ifdef RADIATION
              Erin, Erin_l1, Erin_h1, &
              Erout, Erout_l1, Erout_h1, &
              radflux, radflux_l1, radflux_h1, &
              nstep_fsp, &
#endif
              area, area_l1, area_h1, &
              vol, vol_l1, vol_h1, &
              div, pdivu, lo, hi, dx, dt, &
              mass_lost, xmom_lost, ymom_lost, zmom_lost, &
              eden_lost, xang_lost, yang_lost, zang_lost, &
              verbose)

  if (coord_type .gt. 0) then
     pradial(lo(1):hi(1)+1) = q1(lo(1):hi(1)+1,GDPRES) * dt
  end if

  deallocate(flatn,div,pdivu,q1)

end subroutine ca_ctu_update
