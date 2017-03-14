! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(time, &
                               lo, hi, domlo, domhi, &
                               uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                               uout, uout_l1, uout_l2, uout_h1, uout_h2, &
                               q, q_l1, q_l2, q_h1, q_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
                               srcQ, srQ_l1, srQ_l2, srQ_h1, srQ_h2, &
                               update, updt_l1, updt_l2, updt_h1, updt_h2, &
                               delta, dt, &
                               flux1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
                               flux2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
                               pradial, p_l1, p_l2, p_h1, p_h2, &
                               area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                               area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                               dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                               vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                               courno, verbose) bind(C, name="ca_ctu_update")

  use meth_params_module, only : NQ, QVAR, NVAR, NHYP, NGDNV, GDPRES, &
                                 use_flattening, QU, QV, QW, QPRES, NQAUX, &
                                 first_order_hydro
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, ONE
  use flatten_module, only : uflaten
  use prob_params_module, only : coord_type
  use ctu_advection_module, only : umeth2d, consup

  use bl_fort_module, only : rt => c_real
  implicit none

  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
  integer, intent(in) :: uin_l1, uin_l2, uin_h1, uin_h2
  integer, intent(in) :: uout_l1,uout_l2,uout_h1,uout_h2
  integer, intent(in) :: q_l1, q_l2, q_h1, q_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  integer, intent(in) :: srQ_l1, srQ_l2, srQ_h1, srQ_h2
  integer, intent(in) :: updt_l1,updt_l2,updt_h1,updt_h2
  integer, intent(in) :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
  integer, intent(in) :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
  integer, intent(in) :: p_l1,p_l2,p_h1,p_h2
  integer, intent(in) :: area1_l1,area1_l2,area1_h1,area1_h2
  integer, intent(in) :: area2_l1,area2_l2,area2_h1,area2_h2
  integer, intent(in) :: dloga_l1,dloga_l2,dloga_h1,dloga_h2
  integer, intent(in) :: vol_l1,vol_l2,vol_h1,vol_h2

  real(rt)        , intent(in) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  real(rt)        , intent(inout) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
  real(rt)        , intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,NQ)
  real(rt)        , intent(inout) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  real(rt)        , intent(in) :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,QVAR)
  real(rt)        , intent(inout) :: update(updt_l1:updt_h1,updt_l2:updt_h2,NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
  real(rt)        , intent(inout) :: pradial(p_l1:p_h1,p_l2:p_h2)
  real(rt)        , intent(in) :: area1(area1_l1:area1_h1,area1_l2:area1_h2)
  real(rt)        , intent(in) :: area2(area2_l1:area2_h1,area2_l2:area2_h2)
  real(rt)        , intent(in) :: dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  real(rt)        , intent(in) :: vol(vol_l1:vol_h1,vol_l2:vol_h2)
  real(rt)        , intent(in) :: delta(2), dt, time
  real(rt)        , intent(inout) :: courno

  ! Automatic arrays for workspace
  real(rt)        , allocatable :: flatn(:,:)
  real(rt)        , allocatable :: div(:,:)

  integer ngf
  real(rt)         dx,dy

  integer ::  uin_lo(2),  uin_hi(2)
  integer :: uout_lo(2), uout_hi(2)
  integer :: q_lo(2), q_hi(2)

  integer :: lo_3D(3), hi_3D(3)
  integer :: q_lo_3D(3), q_hi_3D(3)
  integer :: uin_lo_3D(3), uin_hi_3D(3)
  real(rt)         :: dx_3D(3)

  uin_lo  = [uin_l1, uin_l2]
  uin_hi  = [uin_h1, uin_h2]

  uout_lo = [uout_l1, uout_l2]
  uout_hi = [uout_h1, uout_h2]
  
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
     call uflaten([lo(1) - ngf, lo(2) - ngf, 0], [hi(1) + ngf, hi(2) + ngf, 0], &
                  q(:,:,QPRES), q(:,:,QU), q(:,:,QV), q(:,:,QW), &
                  flatn, [q_l1, q_l2, 0], [q_h1, q_h2, 0])
  else
     flatn = ONE
  endif

  ! Do PPM reconstruction


  ! Construct the interface states -- this is essentially just a
  ! reshuffling of interface states from zone-center indexing to
  ! edge-centered indexing


  ! Get the fluxes from the Riemann solver


  ! Compute the artifical viscosity


  ! Make the update for this state


  ! Store fluxes for flux correction
  if (coord_type .eq. 1) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2)) = q1(lo(1):hi(1)+1,lo(2):hi(2),GDPRES) * dt
  end if

  deallocate(flatn,div,q1,q2,pdivu)

end subroutine ca_mol_single_stage
