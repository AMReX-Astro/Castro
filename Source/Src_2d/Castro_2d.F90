subroutine ca_umdrv(is_finest_level, time, &
                    lo, hi, domlo, domhi, &
                    uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                    uout, uout_l1, uout_l2, uout_h1, uout_h2, &
#ifdef RADIATION
                    Erin, Erin_l1, Erin_l2, Erin_h1, Erin_h2, &
                    Erout, Erout_l1, Erout_l2, Erout_h1, Erout_h2, &
#endif
                    q, q_l1, q_l2, q_h1, q_h2, &
                    qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
                    srcQ, srQ_l1, srQ_l2, srQ_h1, srQ_h2, &
                    update, updt_l1, updt_l2, updt_h1, updt_h2, &
                    delta, dt, &
                    flux1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
                    flux2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
#ifdef RADIATION
                    radflux1, radflux1_l1, radflux1_l2, radflux1_h1, radflux1_h2, &
                    radflux2, radflux2_l1, radflux2_l2, radflux2_h1, radflux2_h2, &
#endif
                    pradial, p_l1, p_l2, p_h1, p_h2, &
                    area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                    area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                    dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                    vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                    courno, verbose, &
#ifdef RADIATION
                    nstep_fsp, &
#endif
                    mass_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux, &
                    E_added_flux, mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_umdrv")

  use meth_params_module, only : NQ, QVAR, NVAR, NHYP, NGDNV, GDPRES, &
#ifdef RADIATION
                                 QPTOT, &
#endif
                                 use_flattening, QU, QV, QW, QPRES, NQAUX, &
                                 first_order_hydro
  use advection_util_2d_module, only : divu
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, ONE
  use flatten_module, only : uflaten
  use prob_params_module, only : coord_type
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use flatten_module, only : rad_flaten
#endif
  use advection_module, only : umeth2d, consup

  implicit none

#ifdef RADIATION
  integer, intent(inout) :: nstep_fsp
#endif
  integer, intent(in) :: is_finest_level
  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
  integer, intent(in) :: uin_l1, uin_l2, uin_h1, uin_h2
#ifdef RADIATION
  integer, intent(in) :: Erin_l1, Erin_l2, Erin_h1, Erin_h2
  integer, intent(in) :: Erout_l1, Erout_l2, Erout_h1, Erout_h2
#endif
  integer, intent(in) :: uout_l1,uout_l2,uout_h1,uout_h2
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  integer, intent(in) :: srQ_l1, srQ_l2, srQ_h1, srQ_h2
  integer, intent(in) :: updt_l1,updt_l2,updt_h1,updt_h2
  integer, intent(in) :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
  integer, intent(in) :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
#ifdef RADIATION
  integer, intent(in) :: radflux1_l1,radflux1_l2,radflux1_h1,radflux1_h2
  integer, intent(in) :: radflux2_l1,radflux2_l2,radflux2_h1,radflux2_h2
#endif
  integer, intent(in) :: p_l1,p_l2,p_h1,p_h2
  integer, intent(in) :: area1_l1,area1_l2,area1_h1,area1_h2
  integer, intent(in) :: area2_l1,area2_l2,area2_h1,area2_h2
  integer, intent(in) :: dloga_l1,dloga_l2,dloga_h1,dloga_h2
  integer, intent(in) :: vol_l1,vol_l2,vol_h1,vol_h2

  double precision, intent(in) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  double precision, intent(inout) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
#ifdef RADIATION
  double precision, intent(in) :: Erin(Erin_l1:Erin_h1,Erin_l2:Erin_h2,0:ngroups-1)
  double precision, intent(inout) :: Erout(Erout_l1:Erout_h1,Erout_l2:Erout_h2,0:ngroups-1)
#endif
  double precision, intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,NQ)
  double precision, intent(in) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  double precision, intent(in) :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,QVAR)
  double precision, intent(inout) :: update(updt_l1:updt_h1,updt_l2:updt_h2,NVAR)
  double precision, intent(inout) :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
  double precision, intent(inout) :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
#ifdef RADIATION
  double precision, intent(inout) :: radflux1(radflux1_l1:radflux1_h1,radflux1_l2:radflux1_h2,0:ngroups-1)
  double precision, intent(inout) :: radflux2(radflux2_l1:radflux2_h1,radflux2_l2:radflux2_h2,0:ngroups-1)
#endif
  double precision, intent(inout) :: pradial(p_l1:p_h1,p_l2:p_h2)
  double precision, intent(in) :: area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision, intent(in) :: area2(area2_l1:area2_h1,area2_l2:area2_h2)
  double precision, intent(in) :: dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  double precision, intent(in) :: vol(vol_l1:vol_h1,vol_l2:vol_h2)
  double precision, intent(in) :: delta(2), dt, time
  double precision, intent(inout) :: courno

  double precision, intent(inout) :: E_added_flux, mass_added_flux
  double precision, intent(inout) :: xmom_added_flux, ymom_added_flux, zmom_added_flux
  double precision, intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
  double precision, intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

  ! Automatic arrays for workspace
  double precision, allocatable :: flatn(:,:)
  double precision, allocatable :: div(:,:)
  double precision, allocatable :: pdivu(:,:)

  ! Edge-centered primitive variables (Riemann state)
  double precision, allocatable :: q1(:,:,:)
  double precision, allocatable :: q2(:,:,:)

  integer ngq, ngf
  double precision dx,dy

  integer ::  uin_lo(2),  uin_hi(2)
  integer :: uout_lo(2), uout_hi(2)
  integer :: q_lo(2), q_hi(2)

  integer :: lo_3D(3), hi_3D(3)
  integer :: q_lo_3D(3), q_hi_3D(3)
  integer :: uin_lo_3D(3), uin_hi_3D(3)
  double precision :: dx_3D(3)

  uin_lo  = [uin_l1, uin_l2]
  uin_hi  = [uin_h1, uin_h2]

  uout_lo = [uout_l1, uout_l2]
  uout_hi = [uout_h1, uout_h2]
  
  ngq = NHYP
  ngf = 1

  q_lo = [q_l1, q_l2]
  q_hi = [q_h1, q_h2]

  lo_3D   = [lo(1), lo(2), 0]
  hi_3D   = [hi(1), hi(2), 0]

  q_lo_3D = [q_l1, q_l2, 0]
  q_hi_3D = [q_h1, q_h2, 0]

  uin_lo_3D = [uin_l1, uin_l2, 0]
  uin_hi_3D = [uin_h1, uin_h2, 0]

  dx_3D   = [delta(1), delta(2), ZERO]

  allocate( flatn(q_l1:q_h1,q_l2:q_h2))

  allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))
  allocate( pdivu(lo(1)  :hi(1)  ,lo(2)  :hi(2)))

  allocate(q1(flux1_l1-1:flux1_h1+1,flux1_l2-1:flux1_h2+1,NGDNV))
  allocate(q2(flux2_l1-1:flux2_h1+1,flux2_l2-1:flux2_h2+1,NGDNV))

  dx = delta(1)
  dy = delta(2)

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo_3D, q_hi_3D, &
                   qaux, [qa_l1, qa_l2, 0], [qa_h1, qa_h2, 0], &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  elseif (use_flattening == 1) then
#ifdef RADIATION
     call rad_flaten([lo(1)-ngf, lo(2)-ngf, 0], [hi(1)+ngf, hi(2)+ngf, 0], &
                     q(:,:,qpres), q(:,:,qptot), &
                     q(:,:,QU), q(:,:,QV), q(:,:,QW), &
                     flatn, [q_l1, q_l2, 0], [q_h1, q_h2, 0])
#else
     call uflaten([lo(1) - ngf, lo(2) - ngf, 0], [hi(1) + ngf, hi(2) + ngf, 0], &
                  q(:,:,QPRES), q(:,:,QU), q(:,:,QV), q(:,:,QW), &
                  flatn, [q_l1, q_l2, 0], [q_h1, q_h2, 0])
#endif
  else
     flatn = ONE
  endif

  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth2d(q, q_l1, q_l2, q_h1, q_h2, &
               flatn, &
               qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
               srcQ, srQ_l1, srQ_l2, srQ_h1, srQ_h2, &
               lo(1), lo(2), hi(1), hi(2), dx, dy, dt, &
               uout, uout_l1, uout_l2, uout_h1, uout_h2, &
               flux1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
               flux2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
#ifdef RADIATION
               radflux1,radflux1_l1,radflux1_l2,radflux1_h1,radflux1_h2, &
               radflux2,radflux2_l1,radflux2_l2,radflux2_h1,radflux2_h2, &
#endif               
               q1, flux1_l1-1, flux1_l2-1, flux1_h1+1, flux1_h2+1, &
               q2, flux2_l1-1, flux2_l2-1, flux2_h1+1, flux2_h2+1, &
               area1, area1_l1, area1_l2, area1_h1, area1_h2, &
               area2, area2_l1, area2_l2, area2_h1, area2_h2, &
               pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
               dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
               domlo, domhi)


  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  ! this is used for the artifical viscosity
  call divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2, &
            delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

  ! Conservative update
  call consup(uin,  uin_l1,  uin_l2,  uin_h1,  uin_h2, &
              q, q_l1, q_l2, q_h1, q_h2, &
              uout,  uout_l1, uout_l2, uout_h1, uout_h2, &
              update, updt_l1, updt_l2, updt_h1, updt_h2, &
              q1, flux1_l1-1, flux1_l2-1, flux1_h1+1, flux1_h2+1, &
              q2, flux2_l1-1, flux2_l2-1, flux2_h1+1, flux2_h2+1, &
              flux1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
              flux2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
#ifdef RADIATION
              Erin, Erin_l1, Erin_l2, Erin_h1, Erin_h2, &
              Erout, Erout_l1, Erout_l2, Erout_h1, Erout_h2, &
              radflux1, radflux1_l1, radflux1_l2, radflux1_h1, radflux1_h2, &
              radflux2, radflux2_l1, radflux2_l2, radflux2_h1, radflux2_h2, &
              nstep_fsp, &
#endif
              area1, area1_l1, area1_l2, area1_h1, area1_h2, &
              area2, area2_l1, area2_l2, area2_h1, area2_h2, &
              vol, vol_l1, vol_l2, vol_h1, vol_h2, &
              div, pdivu, lo, hi, dx, dy, dt, &
              mass_added_flux, E_added_flux, &
              xmom_added_flux, ymom_added_flux, zmom_added_flux, &
              mass_lost, xmom_lost, ymom_lost, zmom_lost, &
              eden_lost, xang_lost, yang_lost, zang_lost, &
              verbose)

  if (coord_type .eq. 1) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2)) = q1(lo(1):hi(1)+1,lo(2):hi(2),GDPRES) * dt
  end if

  deallocate(flatn,div,q1,q2,pdivu)

end subroutine ca_umdrv
