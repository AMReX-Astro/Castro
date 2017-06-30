! advection routines in support of the CTU unsplit advection scheme

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
                         flux1, flux1_lo, flux1_hi, &
                         flux2, flux2_lo, flux2_hi, &
#ifdef RADIATION
                         radflux1, radflux1_lo, radflux1_hi, &
                         radflux2, radflux2_lo, radflux2_hi, &
#endif
                         pradial, p_lo, p_hi, &
                         area1, area1_lo, area1_hi, &
                         area2, area2_lo, area2_hi, &
                         dloga, dloga_lo, dloga_hi, &
                         vol, vol_lo, vol_hi, &
                         courno, verbose, &
#ifdef RADIATION
                         nstep_fsp, &
#endif
                         mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                         eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

  use meth_params_module, only : NQ, QVAR, NVAR, NHYP, NGDNV, GDPRES, UMX, &
#ifdef RADIATION
                                 QPTOT, &
#endif
                                 use_flattening, QU, QV, QW, QPRES, NQAUX, &
                                 first_order_hydro
  use advection_util_2d_module, only : divu
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, ONE
  use flatten_module, only : uflatten
  use prob_params_module, only : mom_flux_has_p
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use flatten_module, only : rad_flatten
#endif
  use ctu_advection_module, only : umeth2d, consup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

#ifdef RADIATION
  integer, intent(inout) :: nstep_fsp
#endif
  integer, intent(in) :: is_finest_level
  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
#ifdef RADIATION
  integer, intent(in) :: Erin_lo(3), Erin_hi(3)
  integer, intent(in) :: Erout_lo(3), Erout_hi(3)
#endif
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srQ_lo(3), srQ_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
#ifdef RADIATION
  integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
  integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),NVAR)
  real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),NVAR)
#ifdef RADIATION
  real(rt)        , intent(in) :: Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),0:ngroups-1)
  real(rt)        , intent(inout) :: Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),0:ngroups-1)
#endif
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
  real(rt)        , intent(in) :: srcQ(srQ_lo(1):srQ_hi(1),srQ_lo(2):srQ_hi(2),QVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),NVAR)
#ifdef RADIATION
  real(rt)        , intent(inout) :: radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),0:ngroups-1)
  real(rt)        , intent(inout) :: radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),0:ngroups-1)
#endif
  real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
  real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2))
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2))
  real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))
  real(rt)        , intent(in) :: delta(2), dt, time
  real(rt)        , intent(inout) :: courno

  real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
  real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

  ! Automatic arrays for workspace
  real(rt)        , allocatable :: flatn(:,:)
  real(rt)        , allocatable :: div(:,:)
  real(rt)        , allocatable :: pdivu(:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , allocatable :: q1(:,:,:)
  real(rt)        , allocatable :: q2(:,:,:)

  integer :: ngf
  real(rt)         dx,dy

  integer :: lo_3D(3), hi_3D(3)
  real(rt)         :: dx_3D(3)

  
  ngf = 1

  lo_3D   = [lo(1), lo(2), 0]
  hi_3D   = [hi(1), hi(2), 0]

  dx_3D   = [delta(1), delta(2), ZERO]

  allocate( flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2)))

  allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))
  allocate( pdivu(lo(1)  :hi(1)  ,lo(2)  :hi(2)))

  allocate(q1(flux1_lo(1)-1:flux1_hi(1)+1,flux1_lo(2)-1:flux1_hi(2)+1,NGDNV))
  allocate(q2(flux2_lo(1)-1:flux2_hi(1)+1,flux2_lo(2)-1:flux2_hi(2)+1,NGDNV))

  dx = delta(1)
  dy = delta(2)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  elseif (use_flattening == 1) then
#ifdef RADIATION
     call rad_flatten([lo(1)-ngf, lo(2)-ngf, 0], [hi(1)+ngf, hi(2)+ngf, 0], &
                      q, flatn, q_lo, q_hi)
#else
     call uflatten([lo(1) - ngf, lo(2) - ngf, 0], [hi(1) + ngf, hi(2) + ngf, 0], &
                   q, flatn, q_lo, q_hi, QPRES)
#endif
  else
     flatn = ONE
  endif

  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth2d(q, q_lo, q_hi, &
               flatn, &
               qaux, qa_lo, qa_hi, &
               srcQ, srQ_lo, srQ_hi, &
               lo, hi, dx, dy, dt, &
               uout, uout_lo, uout_hi, &
               flux1, flux1_lo, flux1_hi, &
               flux2, flux2_lo, flux2_hi, &
#ifdef RADIATION
               radflux1, radflux1_lo, radflux1_hi, &
               radflux2, radflux2_lo, radflux2_hi, &
#endif               
               q1, flux1_lo-1, flux1_hi+1, &
               q2, flux2_lo-1, flux2_hi+1, &
               area1, area1_lo, area1_hi, &
               area2, area2_lo, area2_hi, &
               pdivu, vol, vol_lo, vol_hi, &
               dloga, dloga_lo, dloga_hi, &
               domlo, domhi)


  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  ! this is used for the artifical viscosity
  call divu(lo,hi,q,q_lo(1),q_lo(2),q_hi(1),q_hi(2), &
            delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

  ! Conservative update
  call consup(uin,  uin_lo,  uin_hi, &
              q, q_lo, q_hi, &
              uout,  uout_lo, uout_hi, &
              update, updt_lo, updt_hi, &
              q1, flux1_lo-1, flux1_hi+1, &
              q2, flux2_lo-1, flux2_hi+1, &
              flux1, flux1_lo, flux1_hi, &
              flux2, flux2_lo, flux2_hi, &
#ifdef RADIATION
              Erin, Erin_lo, Erin_hi, &
              Erout, Erout_lo, Erout_hi, &
              radflux1, radflux1_lo, radflux1_hi, &
              radflux2, radflux2_lo, radflux2_hi, &
              nstep_fsp, &
#endif
              area1, area1_lo, area1_hi, &
              area2, area2_lo, area2_hi, &
              vol, vol_lo, vol_hi, &
              div, pdivu, lo, hi, dx, dy, dt, &
              mass_lost, xmom_lost, ymom_lost, zmom_lost, &
              eden_lost, xang_lost, yang_lost, zang_lost, &
              verbose)

  if (.not. mom_flux_has_p(1)%comp(UMX)) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2)) = q1(lo(1):hi(1)+1,lo(2):hi(2),GDPRES) * dt
  end if

  deallocate(flatn,div,q1,q2,pdivu)

end subroutine ca_ctu_update
