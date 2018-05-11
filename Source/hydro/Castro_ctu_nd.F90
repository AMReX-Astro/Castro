! advection routines in support of the CTU unsplit advection scheme

subroutine ca_ctu_update(lo, hi, is_finest_level, time, &
                         domlo, domhi, &
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
#if BL_SPACEDIM >= 2
                         flux2, flux2_lo, flux2_hi, &
#endif
#if BL_SPACEDIM == 3
                         flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                         radflux1, radflux1_lo, radflux1_hi, &
#if BL_SPACEDIM >= 2
                         radflux2, radflux2_lo, radflux2_hi, &
#endif
#if BL_SPACEDIM == 3
                         radflux3, radflux3_lo, radflux3_hi, &
#endif
#endif
                         area1, area1_lo, area1_hi, &
#if BL_SPACEDIM >= 2
                         area2, area2_lo, area2_hi, &
#endif
#if BL_SPACEDIM == 3
                         area3, area3_lo, area3_hi, &
#endif
#if BL_SPACEDIM <= 2
                         pradial, p_lo, p_hi, &
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vol, vol_lo, vol_hi, &
                         verbose, &
#ifdef RADIATION
                         nstep_fsp, &
#endif
                         mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                         eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

  use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
  use meth_params_module, only : NQ, QVAR, QPRES, NQAUX, NVAR, NHYP, NGDNV, UMX, GDPRES, &
#ifdef RADIATION
                                 QPTOT, &
#endif
                                 use_flattening, &
                                 first_order_hydro
  use advection_util_module, only : divu
  use amrex_constants_module, only : ZERO, ONE
  use flatten_module, only: uflatten
  use prob_params_module, only : mom_flux_has_p, dg, coord_type
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use flatten_module, only : rad_flatten
#endif
  use ctu_advection_module, only : umeth, consup

  use amrex_fort_module, only : rt => amrex_real
  implicit none

#ifdef RADIATION
  integer, intent(inout) :: nstep_fsp
#endif
  integer, intent(in) :: is_finest_level
  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) :: domlo(3), domhi(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
#ifdef RADIATION
  integer, intent(in) :: Erin_lo(3), Erin_hi(3)
  integer, intent(in) :: Erout_lo(3), Erout_hi(3)
#endif
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srQ_lo(3), srQ_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
#if BL_SPACEDIM >= 2
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
#endif
#if BL_SPACEDIM == 3
  integer, intent(in) :: flux3_lo(3), flux3_hi(3)
#endif
#ifdef RADIATION
  integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
#if BL_SPACEDIM >= 2
  integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
#if BL_SPACEDIM == 3
  integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
#endif
#endif
  integer, intent(in) :: area1_lo(3), area1_hi(3)
#if BL_SPACEDIM >= 2
  integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if BL_SPACEDIM == 3
  integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
  integer, intent(in) :: vol_lo(3), vol_hi(3)
#if BL_SPACEDIM <= 2
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif

  real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
#ifdef RADIATION
  real(rt)        , intent(in) :: Erin(Erin_lo(1):Erin_hi(1), Erin_lo(2):Erin_hi(2), Erin_lo(3):Erin_hi(3), 0:ngroups-1)
  real(rt)        , intent(inout) :: Erout(Erout_lo(1):Erout_hi(1), Erout_lo(2):Erout_hi(2), Erout_lo(3):Erout_hi(3), 0:ngroups-1)
#endif
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
  real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt)        , intent(in) :: srcQ(srQ_lo(1):srQ_hi(1), srQ_lo(2):srQ_hi(2), srQ_lo(3):srQ_hi(3), QVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
#if BL_SPACEDIM >= 2
  real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
#endif
#if BL_SPACEDIM == 3
  real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
#endif
#ifdef RADIATION
  real(rt)        , intent(inout) :: radflux1(radflux1_lo(1):radflux1_hi(1), radflux1_lo(2):radflux1_hi(2), &
                                              radflux1_lo(3):radflux1_hi(3), 0:ngroups-1)
#if BL_SPACEDIM >= 2
  real(rt)        , intent(inout) :: radflux2(radflux2_lo(1):radflux2_hi(1), radflux2_lo(2):radflux2_hi(2), &
                                              radflux2_lo(3):radflux2_hi(3), 0:ngroups-1)
#endif
#if BL_SPACEDIM == 3
  real(rt)        , intent(inout) :: radflux3(radflux3_lo(1):radflux3_hi(1), radflux3_lo(2):radflux3_hi(2), &
                                              radflux3_lo(3):radflux3_hi(3), 0:ngroups-1)
#endif
#endif
  real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if BL_SPACEDIM >= 2
  real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if BL_SPACEDIM == 3
  real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
  real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))

#if BL_SPACEDIM < 3
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
  real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
#endif

  real(rt)        , intent(in) :: delta(3), dt, time

  real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
  real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

  ! Automatic arrays for workspace
  real(rt)        , pointer:: flatn(:,:,:)
  real(rt)        , pointer:: div(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , pointer:: q1(:,:,:,:)
  real(rt)        , pointer:: q2(:,:,:,:)
  real(rt)        , pointer:: q3(:,:,:,:)

  integer :: ngf
  integer :: q1_lo(3), q1_hi(3), q2_lo(3), q2_hi(3), q3_lo(3), q3_hi(3)

  ngf = 1

  call amrex_allocate(   div, lo, hi+dg)

  q1_lo = flux1_lo - dg
  q1_hi = flux1_hi + dg
#if BL_SPACEDIM >= 2
  q2_lo = flux2_lo - dg
  q2_hi = flux2_hi + dg
#endif
#if BL_SPACEDIM == 3
  q3_lo = flux3_lo - dg
  q3_hi = flux3_hi + dg
#endif

  call amrex_allocate(q1, q1_lo, q1_hi, NGDNV)
#if BL_SPACEDIM >= 2
  call amrex_allocate(q2, q2_lo, q2_hi, NGDNV)
#endif
#if BL_SPACEDIM == 3
  call amrex_allocate(q3, q3_lo, q3_hi, NGDNV)
#endif

  ! Compute flattening coefficient for slope calculations.
  call amrex_allocate( flatn, q_lo, q_hi)

  if (first_order_hydro == 1) then
     flatn = ZERO
  elseif (use_flattening == 1) then
#ifdef RADIATION
     call rad_flatten(lo-dg*ngf, hi+dg*ngf, &
                      q, flatn, q_lo, q_hi)
#else
     call uflatten(lo-dg*ngf, hi+dg*ngf, &
                   q, flatn, q_lo, q_hi, QPRES)
#endif
  else
     flatn = ONE
  endif

  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth(q, q_lo, q_hi, &
             flatn, &
             qaux, qa_lo, qa_hi, &
             srcQ, srQ_lo, srQ_hi, &
             lo, hi, delta, dt, &
             uout, uout_lo, uout_hi, &
             flux1, flux1_lo, flux1_hi, &
#if BL_SPACEDIM >= 2
             flux2, flux2_lo, flux2_hi, &
#endif
#if BL_SPACEDIM == 3
             flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
             radflux1, radflux1_lo, radflux1_hi, &
#if BL_SPACEDIM >= 2
             radflux2, radflux2_lo, radflux2_hi, &
#endif
#if BL_SPACEDIM == 3
             radflux3, radflux3_lo, radflux3_hi, &
#endif
#endif
             q1, q1_lo, q1_hi, &
#if BL_SPACEDIM >= 2
             q2, q2_lo, q2_hi, &
#endif
#if BL_SPACEDIM == 3
             q3, q3_lo, q3_hi, &
#endif
              area1, area1_lo, area1_hi, &
#if BL_SPACEDIM >= 2
              area2, area2_lo, area2_hi, &
#endif
#if BL_SPACEDIM == 3
              area3, area3_lo, area3_hi, &
#endif
              vol, vol_lo, vol_hi, &
#if BL_SPACEDIM < 3
              dloga, dloga_lo, dloga_hi, &
#endif
              domlo, domhi)


  call amrex_deallocate( flatn)

  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo, hi, q, q_lo, q_hi, delta, div, lo, hi+dg)

  ! Conservative update
  call consup(uin, uin_lo, uin_hi, &
              q, q_lo, q_hi, &
              uout, uout_lo, uout_hi, &
              update, updt_lo, updt_hi, &
              flux1, flux1_lo, flux1_hi, &
#if BL_SPACEDIM >= 2
              flux2, flux2_lo, flux2_hi, &
#endif
#if BL_SPACEDIM == 3
              flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
              Erin, Erin_lo, Erin_hi, &
              Erout, Erout_lo, Erout_hi, &
              radflux1, radflux1_lo, radflux1_hi, &
#if BL_SPACEDIM >= 2
              radflux2, radflux2_lo, radflux2_hi, &
#endif
#if BL_SPACEDIM == 3
              radflux3, radflux3_lo, radflux3_hi, &
#endif
              nstep_fsp, &
#endif
              q1, q1_lo, q1_hi, &
#if BL_SPACEDIM >= 2
              q2, q2_lo, q2_hi, &
#endif
#if BL_SPACEDIM == 3
              q3, q3_lo, q3_hi, &
#endif
              area1, area1_lo, area1_hi, &
#if BL_SPACEDIM >= 2
              area2, area2_lo, area2_hi, &
#endif
#if BL_SPACEDIM == 3
              area3, area3_lo, area3_hi, &
#endif
              vol, vol_lo, vol_hi, &
              div, lo, hi, delta, dt, &
              mass_lost,xmom_lost,ymom_lost,zmom_lost, &
              eden_lost,xang_lost,yang_lost,zang_lost, &
              verbose)


#if BL_SPACEDIM == 1
  if (coord_type > 0) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
  end if
#endif
#if BL_SPACEDIM == 2
  if (.not. mom_flux_has_p(1)%comp(UMX)) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
  end if
#endif

  call amrex_deallocate(   div)

  call amrex_deallocate(    q1)
#if BL_SPACEDIM >= 2
  call amrex_deallocate(    q2)
#endif
#if BL_SPACEDIM == 3
  call amrex_deallocate(    q3)
#endif

end subroutine ca_ctu_update
