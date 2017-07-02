subroutine ca_ctu_update(is_finest_level, time, &
                         lo, hi, domlo, domhi, &
                         uin, uin_lo, uin_hi, &
                         uout, uout_lo, uout_hi, &
#ifdef RADIATION
                         Erin, Erin_lo, Erin_hi, &
                         Erout, Erout_lo, Erout_hi, &
#endif
                         q, q_lo, q_hi, &
                         qaux, qa_lo, qa_hi, &
                         srcQ, srQ_lo, srQ_hi, &
                         update, updt_lo, updt_hi, &
                         delta, dt, &
                         flux, flux_lo, flux_hi, &
#ifdef RADIATION
                         radflux, radflux_lo, radflux_hi, &
#endif
                         pradial, p_lo, p_hi, &
                         area, area_lo, area_hi, &
                         dloga, dloga_lo, dloga_hi, &
                         vol, vol_lo, vol_hi, &
                         courno, verbose, &
#ifdef RADIATION
                         nstep_fsp, &
#endif
                         mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                         eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

  use meth_params_module, only : NQ, QVAR, QU, QPRES, &
#ifdef RADIATION
                                 QPTOT, &
#endif
                                 NQAUX, NVAR, NHYP, use_flattening, &
                                 NGDNV, GDU, GDPRES, first_order_hydro
  use bl_constants_module, only : ZERO, HALF, ONE
  use advection_util_module, only : compute_cfl
  use flatten_module, only : uflatten
  use prob_params_module, only : coord_type
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use flatten_module, only : rad_flatten
#endif
  use ctu_advection_module  , only : umeth1d, consup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: is_finest_level
  integer, intent(in) :: lo(1), hi(1), verbose
  integer, intent(in) :: domlo(1), domhi(1)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srQ_lo(3),srQ_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: flux_lo(3), flux_hi(3)
#ifdef RADIATION
  integer, intent(inout) :: nstep_fsp
  integer, intent(in) :: Erin_lo(3), Erin_hi(3)
  integer, intent(in) :: Erout_lo(3), Erout_hi(3)
  integer, intent(in) :: radflux_lo(3), radflux_hi(3)
#endif
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: area_lo(3), area_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) ::      uin(  uin_lo(1):  uin_hi(1),NVAR)
  real(rt)        , intent(inout) ::  uout( uout_lo(1): uout_hi(1),NVAR)
  real(rt)        , intent(inout) ::     q(    q_lo(1):    q_hi(1),NQ)
  real(rt)        , intent(in) ::     qaux(   qa_lo(1):   qa_hi(1),NQAUX)
  real(rt)        , intent(in) ::     srcQ(  srQ_lo(1):  srQ_hi(1),QVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1): updt_hi(1),NVAR)
  real(rt)        , intent(inout) ::  flux( flux_lo(1): flux_hi(1),NVAR)
#ifdef RADIATION
  real(rt)        , intent(inout) :: Erout(Erout_lo(1):Erout_hi(1), 0:ngroups-1)
  real(rt)        , intent(inout) :: radflux(radflux_lo(1): radflux_hi(1), 0:ngroups-1)
  real(rt)        , intent(in) :: Erin( Erin_lo(1): Erin_hi(1), 0:ngroups-1)
#endif
  real(rt)        , intent(inout) :: pradial(  p_lo(1):   p_hi(1))
  real(rt)        , intent(in) :: area( area_lo(1): area_hi(1)     )
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1)     )
  real(rt)        , intent(in) ::   vol(  vol_lo(1): vol_hi(1)      )
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

  integer i,ngf

  integer :: lo_3D(3), hi_3D(3)
  real(rt)         :: dx_3D(3)

  ngf = 1

  lo_3D   = [lo(1), 0, 0]
  hi_3D   = [hi(1), 0, 0]
  dx_3D   = [delta(1), ZERO, ZERO]

  allocate( flatn(q_lo(1):q_hi(1)))

  allocate(   div(lo(1):hi(1)+1))
  allocate( pdivu(lo(1):hi(1)  ))
  allocate(    q1(flux_lo(1)-1:flux_hi(1)+1,NGDNV))

  dx = delta(1)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  else if (use_flattening == 1) then
#ifdef RADIATION
     call rad_flatten([lo(1)-ngf, 0, 0], [hi(1)+ngf, 0, 0], &
                      q, flatn, q_lo, q_hi)
#else
     call uflatten([lo(1) - ngf, 0, 0], [hi(1) + ngf, 0, 0], &
                   q, flatn, q_lo, q_hi, QPRES)
#endif
  else
     flatn = ONE
  endif

  call umeth1d(lo, hi, domlo, domhi, &
               q, q_lo, q_hi, &
               flatn, &
               qaux, qa_lo, qa_hi, &
               srcQ, srQ_lo, srQ_hi, &
               lo(1), hi(1), dx, dt, &
               uout, uout_lo, uout_hi, &
               flux, flux_lo, flux_hi, &
#ifdef RADIATION
               radflux,radflux_lo, radflux_hi, &
#endif
               q1, flux_lo-1, flux_hi+1, &
               dloga, dloga_lo, dloga_hi)


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
  call consup(uin, uin_lo, uin_hi, &
              uout, uout_lo, uout_hi, &
              update, updt_lo, updt_hi, &
              q, q_lo, q_hi, &
              flux, flux_lo, flux_hi, &
              q1, flux_lo-1, flux_hi+1, &
#ifdef RADIATION
              Erin, Erin_lo, Erin_hi, &
              Erout, Erout_lo, Erout_hi, &
              radflux, radflux_lo, radflux_hi, &
              nstep_fsp, &
#endif
              area, area_lo, area_hi, &
              vol, vol_lo, vol_hi, &
              div, pdivu, lo, hi, dx, dt, &
              mass_lost, xmom_lost, ymom_lost, zmom_lost, &
              eden_lost, xang_lost, yang_lost, zang_lost, &
              verbose)

  if (coord_type .gt. 0) then
     pradial(lo(1):hi(1)+1) = q1(lo(1):hi(1)+1,GDPRES) * dt
  end if

  deallocate(flatn,div,pdivu,q1)

end subroutine ca_ctu_update
