subroutine ca_umdrv(is_finest_level, time, &
                    lo, hi, domlo, domhi, &
                    uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                    uout, uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
#ifdef RADIATION
                    Erin, Erin_l1, Erin_l2, Erin_l3, Erin_h1, Erin_h2, Erin_h3, &
                    Erout, Erout_l1, Erout_l2, Erout_l3, Erout_h1, Erout_h2, Erout_h3, &
#endif
                    q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                    qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3, &
                    srcQ, srQ_l1, srQ_l2, srQ_l3, srQ_h1, srQ_h2, srQ_h3, &
                    update, updt_l1, updt_l2, updt_l3, updt_h1, updt_h2, updt_h3, &
                    delta, dt, &
                    flux1, flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3, &
                    flux2, flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3, &
                    flux3, flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3, &
#ifdef RADIATION
                    radflux1,radflux1_l1,radflux1_l2,radflux1_l3,radflux1_h1,radflux1_h2,radflux1_h3, &
                    radflux2,radflux2_l1,radflux2_l2,radflux2_l3,radflux2_h1,radflux2_h2,radflux2_h3, &
                    radflux3,radflux3_l1,radflux3_l2,radflux3_l3,radflux3_h1,radflux3_h2,radflux3_h3, &
#endif
                    area1, area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3, &
                    area2, area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3, &
                    area3, area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3, &
                    vol, vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3, &
                    courno, verbose, &
#ifdef RADIATION
                    nstep_fsp, &
#endif
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_umdrv")

  use mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, QU, QV, QW, QPRES, NQAUX, NVAR, NHYP, NGDNV, &
#ifdef RADIATION
                                 QPTOT, &
#endif
                                 use_flattening, &
                                 first_order_hydro
  use advection_util_3d_module, only : divu
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, ONE
  use flatten_module, only: uflaten
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use flatten_module, only : rad_flaten
#endif
  use advection_module, only : umeth3d, consup

  use bl_fort_module, only : rt => c_real
  implicit none

#ifdef RADIATION
  integer, intent(inout) :: nstep_fsp
#endif
  integer, intent(in) :: is_finest_level
  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) ::  domlo(3), domhi(3)
  integer, intent(in) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  integer, intent(in) :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
#ifdef RADIATION
  integer, intent(in) :: Erin_l1, Erin_l2, Erin_l3, Erin_h1, Erin_h2, Erin_h3
  integer, intent(in) :: Erout_l1, Erout_l2, Erout_l3, Erout_h1, Erout_h2, Erout_h3
#endif
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  integer, intent(in) :: srQ_l1, srQ_l2, srQ_l3, srQ_h1, srQ_h2, srQ_h3
  integer, intent(in) :: updt_l1, updt_l2, updt_l3, updt_h1, updt_h2, updt_h3
  integer, intent(in) :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
  integer, intent(in) :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
  integer, intent(in) :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
#ifdef RADIATION
  integer, intent(in) :: radflux1_l1, radflux1_l2, radflux1_l3, radflux1_h1, radflux1_h2, radflux1_h3
  integer, intent(in) :: radflux2_l1, radflux2_l2, radflux2_l3, radflux2_h1, radflux2_h2, radflux2_h3
  integer, intent(in) :: radflux3_l1, radflux3_l2, radflux3_l3, radflux3_h1, radflux3_h2, radflux3_h3
#endif
  integer, intent(in) :: area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3
  integer, intent(in) :: area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3
  integer, intent(in) :: area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3
  integer, intent(in) :: vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3

  real(rt)        , intent(in) :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
  real(rt)        , intent(inout) :: uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NVAR)
#ifdef RADIATION
  real(rt)        , intent(in) :: Erin(Erin_l1:Erin_h1, Erin_l2:Erin_h2, Erin_l3:Erin_h3, 0:ngroups-1)
  real(rt)        , intent(inout) :: Erout(Erout_l1:Erout_h1, Erout_l2:Erout_h2, Erout_l3:Erout_h3, 0:ngroups-1)
#endif
  real(rt)        , intent(inout) :: q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, NQ)
  real(rt)        , intent(inout) :: qaux(qa_l1:qa_h1, qa_l2:qa_h2, qa_l3:qa_h3, NQAUX)
  real(rt)        , intent(in) :: srcQ(srQ_l1:srQ_h1, srQ_l2:srQ_h2, srQ_l3:srQ_h3, QVAR)
  real(rt)        , intent(inout) :: update(updt_l1:updt_h1, updt_l2:updt_h2, updt_l3:updt_h3, NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR)
  real(rt)        , intent(inout) :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR)
#ifdef RADIATION
  real(rt)        , intent(inout) :: radflux1(radflux1_l1:radflux1_h1, radflux1_l2:radflux1_h2, &
                                              radflux1_l3:radflux1_h3, 0:ngroups-1)
  real(rt)        , intent(inout) :: radflux2(radflux2_l1:radflux2_h1, radflux2_l2:radflux2_h2, &
                                              radflux2_l3:radflux2_h3, 0:ngroups-1)
  real(rt)        , intent(inout) :: radflux3(radflux3_l1:radflux3_h1, radflux3_l2:radflux3_h2, &
                                              radflux3_l3:radflux3_h3, 0:ngroups-1)
#endif
  real(rt)        , intent(in) :: area1(area1_l1:area1_h1, area1_l2:area1_h2, area1_l3:area1_h3)
  real(rt)        , intent(in) :: area2(area2_l1:area2_h1, area2_l2:area2_h2, area2_l3:area2_h3)
  real(rt)        , intent(in) :: area3(area3_l1:area3_h1, area3_l2:area3_h2, area3_l3:area3_h3)
  real(rt)        , intent(in) :: vol(vol_l1:vol_h1, vol_l2:vol_h2, vol_l3:vol_h3)
  real(rt)        , intent(in) :: delta(3), dt, time
  real(rt)        , intent(inout) :: courno

  real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
  real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

  ! Automatic arrays for workspace
  real(rt)        , pointer:: flatn(:,:,:)
  real(rt)        , pointer:: div(:,:,:)
  real(rt)        , pointer:: pdivu(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , pointer:: q1(:,:,:,:)
  real(rt)        , pointer:: q2(:,:,:,:)
  real(rt)        , pointer:: q3(:,:,:,:)

  integer :: ngq, ngf
  integer :: uin_lo(3), uin_hi(3)
  integer :: uout_lo(3), uout_hi(3)
#ifdef RADIATION
  integer :: Erin_lo(3), Erin_hi(3)                                             
  integer :: Erout_lo(3), Erout_hi(3)                                           
#endif
  integer :: updt_lo(3), updt_hi(3)
  integer :: flux1_lo(3), flux1_hi(3)
  integer :: flux2_lo(3), flux2_hi(3)
  integer :: flux3_lo(3), flux3_hi(3)
#ifdef RADIATION
  integer :: radflux1_lo(3), radflux1_hi(3)
  integer :: radflux2_lo(3), radflux2_hi(3)
  integer :: radflux3_lo(3), radflux3_hi(3)
#endif
  integer :: area1_lo(3), area1_hi(3)
  integer :: area2_lo(3), area2_hi(3)
  integer :: area3_lo(3), area3_hi(3)
  integer :: vol_lo(3), vol_hi(3)
  integer :: q_lo(3), q_hi(3)
  integer :: qa_lo(3), qa_hi(3)
  integer :: srQ_lo(3), srQ_hi(3)
  integer :: q1_lo(3), q1_hi(3), q2_lo(3), q2_hi(3), q3_lo(3), q3_hi(3)

  ngq = NHYP
  ngf = 1

  q_lo = [ q_l1, q_l2, q_l3 ]
  q_hi = [ q_h1, q_h2, q_h3 ]

  qa_lo = [ qa_l1, qa_l2, qa_l3 ]
  qa_hi = [ qa_h1, qa_h2, qa_h3 ]

  srQ_lo = [ srQ_l1, srQ_l2, srQ_l3 ]
  srQ_hi = [ srQ_h1, srQ_h2, srQ_h3 ]

  uin_lo = [ uin_l1, uin_l2, uin_l3 ]
  uin_hi = [ uin_h1, uin_h2, uin_h3 ]

  uout_lo = [ uout_l1, uout_l2, uout_l3 ]
  uout_hi = [ uout_h1, uout_h2, uout_h3 ]

#ifdef RADIATION
  Erin_lo = [Erin_l1, Erin_l2, Erin_l3]                                         
  Erin_hi = [Erin_h1, Erin_h2, Erin_h3]                                         
                                                                                
  Erout_lo = [Erout_l1, Erout_l2, Erout_l3]                                     
  Erout_hi = [Erout_h1, Erout_h2, Erout_h3]     
#endif

  updt_lo = [ updt_l1, updt_l2, updt_l3 ]
  updt_hi = [ updt_h1, updt_h2, updt_h3 ]

  flux1_lo = [ flux1_l1, flux1_l2, flux1_l3 ]
  flux1_hi = [ flux1_h1, flux1_h2, flux1_h3 ]

  flux2_lo = [ flux2_l1, flux2_l2, flux2_l3 ]
  flux2_hi = [ flux2_h1, flux2_h2, flux2_h3 ]

  flux3_lo = [ flux3_l1, flux3_l2, flux3_l3 ]
  flux3_hi = [ flux3_h1, flux3_h2, flux3_h3 ]

#ifdef RADIATION
  radflux1_lo = [ radflux1_l1, radflux1_l2, radflux1_l3 ]
  radflux1_hi = [ radflux1_h1, radflux1_h2, radflux1_h3 ]

  radflux2_lo = [ radflux2_l1, radflux2_l2, radflux2_l3 ]
  radflux2_hi = [ radflux2_h1, radflux2_h2, radflux2_h3 ]

  radflux3_lo = [ radflux3_l1, radflux3_l2, radflux3_l3 ]
  radflux3_hi = [ radflux3_h1, radflux3_h2, radflux3_h3 ]
#endif

  area1_lo = [ area1_l1, area1_l2, area1_l3 ]
  area1_hi = [ area1_h1, area1_h2, area1_h3 ]

  area2_lo = [ area2_l1, area2_l2, area2_l3 ]
  area2_hi = [ area2_h1, area2_h2, area2_h3 ]

  area3_lo = [ area3_l1, area3_l2, area3_l3 ]
  area3_hi = [ area3_h1, area3_h2, area3_h3 ]

  vol_lo = [ vol_l1, vol_l2, vol_l3 ]
  vol_hi = [ vol_h1, vol_h2, vol_h3 ]

  call bl_allocate(   div, lo(1), hi(1)+1, lo(2), hi(2)+1, lo(3), hi(3)+1)
  call bl_allocate( pdivu, lo(1), hi(1)  , lo(2), hi(2)  , lo(3), hi(3)  )

  q1_lo = flux1_lo - 1
  q1_hi = flux1_hi + 1
  q2_lo = flux2_lo - 1
  q2_hi = flux2_hi + 1
  q3_lo = flux3_lo - 1
  q3_hi = flux3_hi + 1

  call bl_allocate(q1, q1_lo, q1_hi, NGDNV)
  call bl_allocate(q2, q2_lo, q2_hi, NGDNV)
  call bl_allocate(q3, q3_lo, q3_hi, NGDNV)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo, hi, dt, delta, courno)

  ! Compute flattening coefficient for slope calculations.
  call bl_allocate( flatn, q_lo, q_hi)

  if (first_order_hydro == 1) then
     flatn = ZERO
  elseif (use_flattening == 1) then
#ifdef RADIATION
     call rad_flaten([lo(1)-ngf, lo(2)-ngf, lo(3)-ngf], &
                     [hi(1)+ngf, hi(2)+ngf, hi(3)+ngf], &
                     q(:,:,:,qpres), q(:,:,:,qptot), &
                     q(:,:,:,QU), q(:,:,:,QV), q(:,:,:,QW), &
                     flatn, q_lo, q_hi)
#else
     call uflaten(lo - ngf, hi + ngf, &
                  q(:,:,:,QPRES), q(:,:,:,QU), q(:,:,:,QV), q(:,:,:,QW), &
                  flatn, q_lo, q_hi)
#endif
  else
     flatn = ONE
  endif

  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth3d(q, q_lo, q_hi, &
               flatn, &
               qaux, qa_lo, qa_hi, &
               srcQ, srQ_lo, srQ_hi, &
               lo, hi, delta, dt, &
               uout, uout_lo, uout_hi, &
               flux1, flux1_lo, flux1_hi, &
               flux2, flux2_lo, flux2_hi, &
               flux3, flux3_lo, flux3_hi, &
#ifdef RADIATION
               radflux1, radflux1_lo, radflux1_hi, &
               radflux2, radflux2_lo, radflux2_hi, &
               radflux3, radflux3_lo, radflux3_hi, &
#endif
               q1, q1_lo, q1_hi, &
               q2, q2_lo, q2_hi, &
               q3, q3_lo, q3_hi, &
               pdivu, domlo, domhi)


  call bl_deallocate( flatn)

  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo,hi,q,q_lo,q_hi,delta,div,lo,hi+1)

  ! Conservative update
  call consup(uin ,  uin_lo , uin_hi, &
              q, q_lo, q_hi, &
              uout, uout_lo, uout_hi, &
              update, updt_lo, updt_hi, &
              flux1, flux1_lo, flux1_hi, &
              flux2, flux2_lo, flux2_hi, &
              flux3, flux3_lo, flux3_hi, &
#ifdef RADIATION
              Erin, Erin_lo, Erin_hi, &
              Erout, Erout_lo, Erout_hi, &
              radflux1, radflux1_lo, radflux1_hi, &
              radflux2, radflux2_lo, radflux2_hi, &
              radflux3, radflux3_lo, radflux3_hi, &
              nstep_fsp, &
#endif
              q1, q1_lo, q1_hi, &
              q2, q2_lo, q2_hi, &
              q3, q3_lo, q3_hi, &
              area1, area1_lo, area1_hi, &
              area2, area2_lo, area2_hi, &
              area3, area3_lo, area3_hi, &
              vol, vol_lo, vol_hi, &
              div,pdivu,lo,hi,delta,dt, &
              mass_lost,xmom_lost,ymom_lost,zmom_lost, &
              eden_lost,xang_lost,yang_lost,zang_lost, &
              verbose)

  call bl_deallocate(   div)
  call bl_deallocate( pdivu)

  call bl_deallocate(    q1)
  call bl_deallocate(    q2)
  call bl_deallocate(    q3)

end subroutine ca_umdrv
