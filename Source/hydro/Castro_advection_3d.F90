module ctu_advection_module

  use bl_constants_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public umeth, consup

contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH     Compute hyperbolic fluxes using unsplit second
! ::: ::           order Godunov integrator.
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: qaux        => (const)  auxiliary hydro data
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: nz          => (const)  number of cells in Z direction
! ::: :: dx          => (const)  grid spacing in X, Y, Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
! ::: ----------------------------------------------------------------

  subroutine umeth(q, qd_lo, qd_hi, &
                   flatn, &
                   qaux, qa_lo, qa_hi, &
                   srcQ, src_lo, src_hi, &
                   lo, hi, dx, dt, &
                   uout, uout_lo, uout_hi, &
                   flux1, fd1_lo, fd1_hi, &
                   flux2, fd2_lo, fd2_hi, &
                   flux3, fd3_lo, fd3_hi, &
#ifdef RADIATION
                   rflux1, rfd1_lo, rfd1_hi, &
                   rflux2, rfd2_lo, rfd2_hi, &
                   rflux3, rfd3_lo, rfd3_hi, &
#endif
                   q1, q1_lo, q1_hi, &
                   q2, q2_lo, q2_hi, &
                   q3, q3_lo, q3_hi, &
                   area1, area1_lo, area1_hi, &
                   area2, area2_lo, area2_hi, &
                   area3, area3_lo, area3_hi, &
                   vol, vol_lo, vol_hi, &
                   domlo, domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NQ, NVAR, QPRES, QRHO, QU, QW, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   ppm_type, &
                                   use_pslope, ppm_temp_fix, &
                                   hybrid_riemann
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
    use trace_module, only : tracexy, tracez
    use transverse_module, only : transx1, transx2, transy1, transy2, transz, &
                                  transxy, transxz, transyz
    use ppm_module, only : ppm_reconstruct, ppm_int_profile
    use slope_module, only : uslope, pslope
    use actual_network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use riemann_module, only: cmpflx
    use bl_constants_module
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : tracexy_ppm_rad, tracez_ppm_rad
#else
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use advection_util_module, only : shock
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: fd1_lo(3), fd1_hi(3)
    integer, intent(in) :: fd2_lo(3), fd2_hi(3)
    integer, intent(in) :: fd3_lo(3), fd3_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
#ifdef RADIATION
    integer, intent(in) :: rfd1_lo(3), rfd1_hi(3)
    integer, intent(in) :: rfd2_lo(3), rfd2_hi(3)
    integer, intent(in) :: rfd3_lo(3), rfd3_hi(3)
#endif

    real(rt)        , intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)        , intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)        , intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    real(rt)        , intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)

    real(rt)        , intent(inout) ::  uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt)        , intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    real(rt)        , intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    real(rt)        , intent(inout) :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)
    real(rt)        , intent(inout) ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    real(rt)        , intent(inout) ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt)        , intent(inout) ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    real(rt)        , intent(in) :: dx(3), dt

#ifdef RADIATION
    real(rt)         rflux1(rfd1_lo(1):rfd1_hi(1),rfd1_lo(2):rfd1_hi(2),rfd1_lo(3):rfd1_hi(3),0:ngroups-1)
    real(rt)         rflux2(rfd2_lo(1):rfd2_hi(1),rfd2_lo(2):rfd2_hi(2),rfd2_lo(3):rfd2_hi(3),0:ngroups-1)
    real(rt)         rflux3(rfd3_lo(1):rfd3_hi(1),rfd3_lo(2):rfd3_hi(2),rfd3_lo(3):rfd3_hi(3),0:ngroups-1)
#endif

    real(rt)         :: dxinv, dyinv, dzinv
    real(rt)         :: dtdx, dtdy, dtdz, hdt
    real(rt)         :: cdtdx, cdtdy, cdtdz
    real(rt)         :: hdtdx, hdtdy, hdtdz

    integer :: km, kc, kt, k3d, n
    integer :: i, j, iwave, idim

    ! Left and right state arrays (edge centered, cell centered)
    real(rt)        , pointer :: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
    real(rt)        , pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
    real(rt)        , pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

    real(rt)        , pointer :: qmxy(:,:,:,:), qpxy(:,:,:,:)
    real(rt)        , pointer :: qmxz(:,:,:,:), qpxz(:,:,:,:)

    real(rt)        , pointer :: qmyx(:,:,:,:), qpyx(:,:,:,:)
    real(rt)        , pointer :: qmyz(:,:,:,:), qpyz(:,:,:,:)

    real(rt)        , pointer :: qmzx(:,:,:,:), qpzx(:,:,:,:)
    real(rt)        , pointer :: qmzy(:,:,:,:), qpzy(:,:,:,:)

    real(rt)        , pointer :: qxl(:,:,:,:), qxr(:,:,:,:)
    real(rt)        , pointer :: qyl(:,:,:,:), qyr(:,:,:,:)
    real(rt)        , pointer :: qzl(:,:,:,:), qzr(:,:,:,:)

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    real(rt)        , pointer ::  fx(:,:,:,:), fy(:,:,:,:), fz(:,:,:,:)

    real(rt)        , pointer :: fxy(:,:,:,:), fxz(:,:,:,:)
    real(rt)        , pointer :: fyx(:,:,:,:), fyz(:,:,:,:)
    real(rt)        , pointer :: fzx(:,:,:,:), fzy(:,:,:,:)

    real(rt)        , pointer :: qgdnvx(:,:,:,:), qgdnvxf(:,:,:,:), qgdnvtmpx(:,:,:,:)
    real(rt)        , pointer :: qgdnvy(:,:,:,:), qgdnvyf(:,:,:,:), qgdnvtmpy(:,:,:,:)
    real(rt)        , pointer :: qgdnvz(:,:,:,:), qgdnvzf(:,:,:,:), qgdnvtmpz1(:,:,:,:), qgdnvtmpz2(:,:,:,:)

#ifdef RADIATION
    real(rt)        , pointer ::  rfx(:,:,:,:), rfy(:,:,:,:),rfz(:,:,:,:)
    real(rt)        , pointer ::rfxy(:,:,:,:),rfxz(:,:,:,:)
    real(rt)        , pointer ::rfyx(:,:,:,:),rfyz(:,:,:,:)
    real(rt)        , pointer ::rfzx(:,:,:,:),rfzy(:,:,:,:)
#endif

    real(rt)        , pointer :: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    real(rt)        , pointer :: Ip_src(:,:,:,:,:,:), Im_src(:,:,:,:,:,:)
    real(rt)        , pointer :: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    real(rt)        , pointer :: shk(:,:,:)

    real(rt)        , pointer :: sxm(:,:,:), sym(:,:,:), szm(:,:,:)
    real(rt)        , pointer :: sxp(:,:,:), syp(:,:,:), szp(:,:,:)

    type (eos_t) :: eos_state

    integer :: qt_lo(3), qt_hi(3)
    integer :: It_lo(3), It_hi(3)
    integer :: shk_lo(3), shk_hi(3)
    integer :: fx_lo(3), fx_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: fz_lo(3), fz_hi(3)

    qt_lo = [lo(1) - 1, lo(2) - 1, 1]
    qt_hi = [hi(1) + 2, hi(2) + 2, 2]

    It_lo = [lo(1) - 1, lo(2) - 1, 1]
    It_hi = [hi(1) + 1, hi(2) + 1, 2]

    shk_lo(:) = lo(:) - 1
    shk_hi(:) = hi(:) + 1

    fx_lo = [lo(1)    , lo(2) - 1, 1]
    fx_hi = [hi(1) + 1, hi(2) + 1, 2]

    fy_lo = [lo(1) - 1, lo(2)    , 1]
    fy_hi = [hi(1) + 1, hi(2) + 1, 2]

    fz_lo = [lo(1) - 1, lo(2) - 1, 1]
    fz_hi = [hi(1) + 1, hi(2) + 1, 2]

    call bl_allocate (     qgdnvx, qt_lo, qt_hi, NGDNV)
    call bl_allocate (    qgdnvxf, qt_lo, qt_hi, NGDNV)
    call bl_allocate (  qgdnvtmpx, qt_lo, qt_hi, NGDNV)

    call bl_allocate (     qgdnvy, qt_lo, qt_hi, NGDNV)
    call bl_allocate (    qgdnvyf, qt_lo, qt_hi, NGDNV)
    call bl_allocate (  qgdnvtmpy, qt_lo, qt_hi, NGDNV)

    call bl_allocate (     qgdnvz, qt_lo, qt_hi, NGDNV)
    call bl_allocate (    qgdnvzf, qt_lo, qt_hi, NGDNV)
    call bl_allocate ( qgdnvtmpz1, qt_lo, qt_hi, NGDNV)
    call bl_allocate ( qgdnvtmpz2, qt_lo, qt_hi, NGDNV)

    call bl_allocate ( qxm, qt_lo, qt_hi, NQ)
    call bl_allocate ( qxp, qt_lo, qt_hi, NQ)

    call bl_allocate ( qmxy, qt_lo, qt_hi, NQ)
    call bl_allocate ( qpxy, qt_lo, qt_hi, NQ)

    call bl_allocate ( qmxz, qt_lo, qt_hi, NQ)
    call bl_allocate ( qpxz, qt_lo, qt_hi, NQ)

    call bl_allocate ( qym, qt_lo, qt_hi, NQ)
    call bl_allocate ( qyp, qt_lo, qt_hi, NQ)

    call bl_allocate ( qmyx, qt_lo, qt_hi, NQ)
    call bl_allocate ( qpyx, qt_lo, qt_hi, NQ)

    call bl_allocate ( qmyz, qt_lo, qt_hi, NQ)
    call bl_allocate ( qpyz, qt_lo, qt_hi, NQ)

    call bl_allocate ( qzm, qt_lo, qt_hi, NQ)
    call bl_allocate ( qzp, qt_lo, qt_hi, NQ)

    call bl_allocate ( qxl, qt_lo, qt_hi, NQ)
    call bl_allocate ( qxr, qt_lo, qt_hi, NQ)
    call bl_allocate ( qyl, qt_lo, qt_hi, NQ)
    call bl_allocate ( qyr, qt_lo, qt_hi, NQ)
    call bl_allocate ( qzl, qt_lo, qt_hi, NQ)
    call bl_allocate ( qzr, qt_lo, qt_hi, NQ)

    call bl_allocate ( qmzx, qt_lo, qt_hi, NQ)
    call bl_allocate ( qpzx, qt_lo, qt_hi, NQ)

    call bl_allocate ( qmzy, qt_lo, qt_hi, NQ)
    call bl_allocate ( qpzy, qt_lo, qt_hi, NQ)

    call bl_allocate ( fx, fx_lo, fx_hi, NVAR)
    call bl_allocate ( fy, fy_lo, fy_hi, NVAR)
    call bl_allocate ( fz, fz_lo, fz_hi, NVAR)

    call bl_allocate ( fxy, fx_lo, fx_hi, NVAR)
    call bl_allocate ( fxz, fx_lo, fx_hi, NVAR)

    call bl_allocate ( fyx, fy_lo, fy_hi, NVAR)
    call bl_allocate ( fyz, fy_lo, fy_hi, NVAR)

    call bl_allocate ( fzx, fz_lo, fz_hi, NVAR)
    call bl_allocate ( fzy, fz_lo, fz_hi, NVAR)

#ifdef RADIATION
    call bl_allocate (rfx, fx_lo(1), fx_hi(1), fx_lo(2), fx_hi(2), fx_lo(3), fx_hi(3), 0, ngroups-1)
    call bl_allocate (rfy, fy_lo(1), fy_hi(1), fy_lo(2), fy_hi(2), fy_lo(3), fy_hi(3), 0, ngroups-1)
    call bl_allocate (rfz, fz_lo(1), fz_hi(1), fz_lo(2), fz_hi(2), fz_lo(3), fz_hi(3), 0, ngroups-1)
    call bl_allocate (rfxy, fx_lo(1), fx_hi(1), fx_lo(2), fx_hi(2), fx_lo(3), fx_hi(3), 0, ngroups-1)
    call bl_allocate (rfxz, fx_lo(1), fx_hi(1), fx_lo(2), fx_hi(2), fx_lo(3), fx_hi(3), 0, ngroups-1)
    call bl_allocate (rfyx, fy_lo(1), fy_hi(1), fy_lo(2), fy_hi(2), fy_lo(3), fy_hi(3), 0, ngroups-1)
    call bl_allocate (rfyz, fy_lo(1), fy_hi(1), fy_lo(2), fy_hi(2), fy_lo(3), fy_hi(3), 0, ngroups-1)
    call bl_allocate (rfzx, fz_lo(1), fz_hi(1), fz_lo(2), fz_hi(2), fz_lo(3), fz_hi(3), 0, ngroups-1)
    call bl_allocate (rfzy, fz_lo(1), fz_hi(1), fz_lo(2), fz_hi(2), fz_lo(3), fz_hi(3), 0, ngroups-1)
#endif

    if (ppm_type .gt. 0) then
       ! x-index, y-index, z-index, dim, characteristics, variables
       call bl_allocate ( Ip, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,NQ)
       call bl_allocate ( Im, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,NQ)

       ! for source terms
       call bl_allocate ( Ip_src, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,NQ)
       call bl_allocate ( Im_src, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,NQ)

       ! for gamc -- needed for the reference state in eigenvectors
       call bl_allocate ( Ip_gc, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,1)
       call bl_allocate ( Im_gc, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,1)
    else
       call bl_allocate ( dqx, qt_lo, qt_hi, NQ)
       call bl_allocate ( dqy, qt_lo, qt_hi, NQ)
       call bl_allocate ( dqz, qt_lo, qt_hi, NQ)
    end if

    ! for the hybrid Riemann solver
    call bl_allocate(shk, shk_lo, shk_hi)

    ! Local constants
    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)
    dtdx = dt*dxinv
    dtdy = dt*dyinv
    dtdz = dt*dzinv
    hdt = HALF*dt
    hdtdx = HALF*dtdx
    hdtdy = HALF*dtdy
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD

#ifdef SHOCK_VAR
    uout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),USHK) = ZERO

    call shock(q,qd_lo,qd_hi,shk,shk_lo,shk_hi,lo,hi,dx)

    ! Store the shock data for future use in the burning step.

    do k3d = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,k3d,USHK) = shk(i,j,k3d)
          enddo
       enddo
    enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:,:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q,qd_lo,qd_hi,shk,shk_lo,shk_hi,lo,hi,dx)
    else
       shk(:,:,:) = ZERO
    endif
#endif

    ! We come into this routine with a 3-d box of data, but we operate
    ! on it locally by considering 2 planes that encompass all of the
    ! x, y indices of the original box, but each plane corresponds to
    ! a single z index.
    !
    ! In the notation below, k3d will always been the index into the
    ! original 3-d box.  kc will be the z-index in the local "planar"
    ! data and km will be the previously used index in the local
    ! planar data.
    !
    ! With each loop in the k direction, we will overwrite the old
    ! data in the planar arrays.


    ! Initialize kc (current k-level) and km (previous k-level)
    kc = 1
    km = 2


    call bl_allocate(sxm, It_lo, It_hi)
    call bl_allocate(sxp, It_lo, It_hi)
    call bl_allocate(sym, It_lo, It_hi)
    call bl_allocate(syp, It_lo, It_hi)
    call bl_allocate(szm, It_lo, It_hi)
    call bl_allocate(szp, It_lo, It_hi)

    do k3d = lo(3)-1, hi(3)+1

       ! Swap pointers to levels
       kt = km
       km = kc
       kc = kt

       if (ppm_type .gt. 0) then

          do n = 1, NQ
             call ppm_reconstruct(q(:,:,:,n  ), qd_lo, qd_hi, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

             call ppm_int_profile(q(:,:,:,n  ), qd_lo, qd_hi, &
                                  q(:,:,:,QU:QW), qd_lo, qd_hi, &
                                  qaux(:,:,:,QC), qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  Ip(:,:,:,:,:,n), Im(:,:,:,:,:,n), It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, dt, k3d, kc)
          end do

          ! source terms
          do n=1,QVAR
             call ppm_reconstruct(srcQ(:,:,:,n), src_lo, src_hi, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

             call ppm_int_profile(srcQ(:,:,:,n), src_lo, src_hi, &
                                  q(:,:,:,QU:QW), qd_lo, qd_hi, &
                                  qaux(:,:,:,QC), qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  Ip_src(:,:,:,:,:,n), Im_src(:,:,:,:,:,n), It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, dt, k3d, kc)
          enddo

          ! this probably doesn't support radiation
          if (ppm_temp_fix /= 1) then
             call ppm_reconstruct(qaux(:,:,:,QGAMC), qa_lo, qa_hi, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

             call ppm_int_profile(qaux(:,:,:,QGAMC), qa_lo, qa_hi, &
                                  q(:,:,:,QU:QW), qd_lo, qd_hi, &
                                  qaux(:,:,:,QC), qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  Ip_gc(:,:,:,:,:,1), Im_gc(:,:,:,:,:,1), It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, dt, k3d, kc)
          else

             do iwave = 1, 3
                do idim = 1, 3

                   do j = lo(2)-1, hi(2)+1
                      do i = lo(1)-1, hi(1)+1
                         eos_state % rho = Ip(i,j,kc,idim,iwave,QRHO)
                         eos_state % T   = Ip(i,j,kc,idim,iwave,QTEMP)

                         eos_state % xn  = Ip(i,j,kc,idim,iwave,QFS:QFS+nspec-1)
                         eos_state % aux = Ip(i,j,kc,idim,iwave,QFX:QFX+naux-1)

                         call eos(eos_input_rt, eos_state)

                         Ip(i,j,kc,idim,iwave,QPRES)  = eos_state % p
                         Ip(i,j,kc,idim,iwave,QREINT) = eos_state % e * Ip(i,j,kc,idim,iwave,QRHO)
                         Ip_gc(i,j,kc,idim,iwave,1)   = eos_state % gam1
                      enddo
                   enddo

                   do j = lo(2)-1, hi(2)+1
                      do i = lo(1)-1, hi(1)+1
                         eos_state % rho = Im(i,j,kc,idim,iwave,QRHO)
                         eos_state % T   = Im(i,j,kc,idim,iwave,QTEMP)

                         eos_state % xn  = Im(i,j,kc,idim,iwave,QFS:QFS+nspec-1)
                         eos_state % aux = Im(i,j,kc,idim,iwave,QFX:QFX+naux-1)

                         call eos(eos_input_rt, eos_state)

                         Im(i,j,kc,idim,iwave,QPRES)  = eos_state % p
                         Im(i,j,kc,idim,iwave,QREINT) = eos_state % e * Im(i,j,kc,idim,iwave,QRHO)
                         Im_gc(i,j,kc,idim,iwave,1)   = eos_state % gam1
                      enddo
                   enddo

                enddo
             enddo

          endif

          ! Compute U_x and U_y at kc (k3d)

#ifdef RADIATION
          call tracexy_ppm_rad(q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               Ip, Im, Ip_src, Im_src, It_lo, It_hi, &
                               qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                               lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                               dx, dt, kc, k3d)
#else
          call tracexy_ppm(q, qd_lo, qd_hi, &
                           qaux, qa_lo, qa_hi, &
                           Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, It_lo, It_hi, &
                           qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                           lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                           dx, dt, kc, k3d)
#endif

       else

#ifdef RADIATION
          call bl_error("ppm_type <=0 is not supported in with radiation")
#endif

          ! Compute all slopes at kc (k3d)
          call uslope(q, flatn, qd_lo, qd_hi, &
                      dqx, dqy, dqz, qt_lo, qt_hi, &
                      lo(1), lo(2), hi(1), hi(2), kc, k3d)

          if (use_pslope .eq. 1) &
               call pslope(q, flatn, qd_lo, qd_hi, &
                           dqx, dqy, dqz, qt_lo, qt_hi, &
                           srcQ, src_lo, src_hi, &
                           lo(1), lo(2), hi(1), hi(2), kc, k3d, dx)

          ! Compute U_x and U_y at kc (k3d)
          call tracexy(q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       dqx, dqy, qt_lo, qt_hi, &
                       qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                       lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                       dx, dt, kc, k3d)

       end if

       ! Compute \tilde{F}^x at kc (k3d)
       call cmpflx(qxm, qxp, qt_lo, qt_hi, &
                   fx, fx_lo, fx_hi, &
                   qgdnvx, qt_lo, qt_hi, &
#ifdef RADIATION
                   rfx, fx_lo, fx_hi, &
#endif
                   qaux, qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym, qyp, qt_lo, qt_hi, &
                   fy, fy_lo, fy_hi, &
                   qgdnvy, qt_lo, qt_hi, &
#ifdef RADIATION
                   rfy, fy_lo, fy_hi, &
#endif
                   qaux, qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, kc, kc, k3d, domlo, domhi)

       ! Compute U'^y_x at kc (k3d)
       call transy1(qxm, qmxy, qxp, qpxy, qt_lo, qt_hi, &
                    qaux, qa_lo, qa_hi, &
                    fy, &
#ifdef RADIATION
                    rfy, &
#endif
                    fy_lo, fy_hi, &
                    qgdnvy, qt_lo, qt_hi, &
                    cdtdy, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, k3d)

       ! Compute U'^x_y at kc (k3d)
       call transx1(qym, qmyx, qyp, qpyx, qt_lo, qt_hi, &
                    qaux, qa_lo, qa_hi, &
                    fx, &
#ifdef RADIATION
                    rfx, &
#endif
                    fx_lo, fx_hi, &
                    qgdnvx, qt_lo, qt_hi, &
                    cdtdx, lo(1), hi(1), lo(2)-1, hi(2)+1,kc,k3d)

       ! Compute F^{x|y} at kc (k3d)
       call cmpflx(qmxy, qpxy, qt_lo, qt_hi, &
                   fxy, fx_lo, fx_hi, &
                   qgdnvtmpx, qt_lo, qt_hi, &
#ifdef RADIATION
                   rfxy, fx_lo, fx_hi, &
#endif
                   qaux, qa_lo, qa_hi, &
                   shk,shk_lo,shk_hi, &
                   1, lo(1), hi(1)+1, lo(2), hi(2), kc, kc, k3d, domlo, domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx, qpyx, qt_lo, qt_hi, &
                   fyx, fy_lo, fy_hi, &
                   qgdnvtmpy, qt_lo, qt_hi, &
#ifdef RADIATION
                   rfyx, fy_lo, fy_hi, &
#endif
                   qaux, qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   2, lo(1), hi(1), lo(2), hi(2)+1, kc, kc, k3d, domlo, domhi)

       if (k3d.ge.lo(3)) then

          ! Compute U_z at kc (k3d)
          if (ppm_type .gt. 0) then
#ifdef RADIATION
             call tracez_ppm_rad(q, qd_lo, qd_hi, &
                                 qaux, qa_lo, qa_hi, &
                                 Ip, Im, Ip_src, Im_src, It_lo, It_hi, &
                                 qzm, qzp, qt_lo, qt_hi, &
                                 lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                                 dt, km, kc, k3d)

#else
             call tracez_ppm(q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, It_lo, It_hi, &
                             qzm, qzp, qt_lo, qt_hi, &
                             lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                             dt, km, kc, k3d)
#endif
          else
             ! we should not land here with radiation
             call tracez(q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         dqz, qt_lo, qt_hi, &
                         qzm, qzp, qt_lo, qt_hi, &
                         lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                         dx, dt, km, kc, k3d)
          end if

          ! Compute \tilde{F}^z at kc (k3d)
          call cmpflx(qzm,qzp,qt_lo,qt_hi, &
                      fz,fz_lo,fz_hi, &
                      qgdnvz,qt_lo,qt_hi, &
#ifdef RADIATION
                      rfz, fz_lo, fz_hi, &
#endif
                      qaux, qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1)-1, hi(1)+1, lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

          ! Compute U'^y_z at kc (k3d)
          call transy2(qzm, qmzy, qzp, qpzy, qt_lo, qt_hi, &
                       qaux, qa_lo, qa_hi, &
                       fy, &
#ifdef RADIATION
                       rfy, &
#endif
                       fy_lo, fy_hi, &
                       qgdnvy, qt_lo, qt_hi, &
                       cdtdy, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, km, k3d)

          ! Compute U'^x_z at kc (k3d)
          call transx2(qzm, qmzx, qzp, qpzx, qt_lo, qt_hi, &
                       qaux, qa_lo, qa_hi, &
                       fx, &
#ifdef RADIATION
                       rfx, &
#endif
                       fx_lo, fx_hi, &
                       qgdnvx, qt_lo, qt_hi, &
                       cdtdx, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, km, k3d)

          ! Compute F^{z|x} at kc (k3d)
          call cmpflx(qmzx, qpzx, qt_lo, qt_hi, &
                      fzx, fz_lo, fz_hi, &
                      qgdnvtmpz1, qt_lo, qt_hi, &
#ifdef RADIATION
                      rfzx, fz_lo, fz_hi, &
#endif
                      qaux, qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy, qpzy, qt_lo, qt_hi, &
                      fzy, fz_lo, fz_hi, &
                      qgdnvtmpz2, qt_lo, qt_hi, &
#ifdef RADIATION
                      rfzy, fz_lo, fz_hi, &
#endif
                      qaux, qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, kc, k3d, domlo, domhi)

          ! Compute U''_z at kc (k3d)
          call transxy(qzm, qzl, qzp, qzr, qt_lo, qt_hi, &
                       qaux, qa_lo, qa_hi, &
                       fxy, &
#ifdef RADIATION
                       rfxy, &
#endif
                       fx_lo, fx_hi, &
                       fyx,&
#ifdef RADIATION
                       rfyx, &
#endif
                       fy_lo, fy_hi, &
                       qgdnvtmpx, qt_lo, qt_hi, &
                       qgdnvtmpy, qt_lo, qt_hi, &
                       srcQ, src_lo, src_hi,&
                       hdt, hdtdx, hdtdy, lo(1), hi(1), lo(2), hi(2), kc, km, k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl, qzr, qt_lo, qt_hi, &
                      flux3, fd3_lo, fd3_hi, &
                      qgdnvzf, qt_lo, qt_hi, &
#ifdef RADIATION
                      rflux3, rfd3_lo, rfd3_hi, &
#endif
                      qaux, qa_lo, qa_hi, &
                      shk,shk_lo,shk_hi, &
                      3,lo(1),hi(1),lo(2),hi(2),kc,k3d,k3d,domlo,domhi)

          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                q3(i,j,k3d,:) = qgdnvzf(i,j,kc,:)
             end do
          end do

          if (k3d.gt.lo(3)) then

             ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
             call transz(qxm, qmxz, qxp, qpxz, qym, qmyz, qyp, qpyz, qt_lo, qt_hi, &
                         qaux, qa_lo, qa_hi, &
                         fz, &
#ifdef RADIATION
                         rfz, &
#endif
                         fz_lo, fz_hi, &
                         qgdnvz,qt_lo,qt_hi, &
                         cdtdz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,km,kc,k3d)

             ! Compute F^{x|z} at km (k3d-1)
             call cmpflx(qmxz, qpxz, qt_lo, qt_hi, &
                         fxz, fx_lo, fx_hi, &
                         qgdnvx, qt_lo, qt_hi, &
#ifdef RADIATION
                         rfxz, fx_lo, fx_hi, &
#endif
                         qaux, qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, km, km, k3d-1, domlo, domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz, qpyz, qt_lo, qt_hi, &
                         fyz, fy_lo, fy_hi, &
                         qgdnvy, qt_lo, qt_hi, &
#ifdef RADIATION
                         rfyz, fy_lo, fy_hi, &
#endif
                         qaux, qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, km, km, k3d-1, domlo, domhi)

             ! Compute U''_x at km (k3d-1)
             call transyz(qxm, qxl, qxp, qxr, qt_lo, qt_hi, &
                          qaux, qa_lo, qa_hi, &
                          fyz, &
#ifdef RADIATION
                          rfyz, &
#endif
                          fy_lo, fy_hi, &
                          fzy, &
#ifdef RADIATION
                          rfzy, &
#endif
                          fz_lo, fz_hi, &
                          qgdnvy, qt_lo, qt_hi, &
                          qgdnvtmpz2, qt_lo, qt_hi, &
                          srcQ, src_lo, src_hi, &
                          hdt, hdtdy, hdtdz, lo(1)-1, hi(1)+1, lo(2), hi(2), km, kc, k3d-1)

             ! Compute U''_y at km (k3d-1)
             call transxz(qym, qyl, qyp, qyr, qt_lo, qt_hi, &
                          qaux, qa_lo, qa_hi, &
                          fxz, &
#ifdef RADIATION
                          rfxz, &
#endif
                          fx_lo, fx_hi, &
                          fzx, &
#ifdef RADIATION
                          rfzx, &
#endif
                          fz_lo, fz_hi, &
                          qgdnvx, qt_lo, qt_hi, &
                          qgdnvtmpz1, qt_lo, qt_hi, &
                          srcQ, src_lo, src_hi, &
                          hdt, hdtdx, hdtdz, lo(1), hi(1), lo(2)-1, hi(2)+1, km, kc, k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl, qxr, qt_lo, qt_hi, &
                         flux1, fd1_lo, fd1_hi, &
                         qgdnvxf, qt_lo, qt_hi, &
#ifdef RADIATION
                         rflux1, rfd1_lo, rfd1_hi, &
#endif
                         qaux, qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         1, lo(1), hi(1)+1, lo(2), hi(2), km, k3d-1, k3d-1, domlo, domhi)

             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+2
                   q1(i,j,k3d-1,:) = qgdnvxf(i,j,km,:)
                end do
             end do

             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl, qyr, qt_lo, qt_hi, &
                         flux2, fd2_lo, fd2_hi, &
                         qgdnvyf, qt_lo, qt_hi, &
#ifdef RADIATION
                         rflux2, rfd2_lo, rfd2_hi, &
#endif
                         qaux, qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         2, lo(1), hi(1), lo(2), hi(2)+1, km, k3d-1, k3d-1, domlo, domhi)

             do j=lo(2)-1,hi(2)+2
                do i=lo(1)-1,hi(1)+1
                   q2(i,j,k3d-1,:) = qgdnvyf(i,j,km,:)
                end do
             end do

          end if
       end if
    enddo

    ! Deallocate arrays
    call bl_deallocate(sxm)
    call bl_deallocate(sxp)
    call bl_deallocate(sym)
    call bl_deallocate(syp)
    call bl_deallocate(szm)
    call bl_deallocate(szp)

    call bl_deallocate ( qgdnvx)
    call bl_deallocate ( qgdnvxf)
    call bl_deallocate ( qgdnvtmpx)

    call bl_deallocate ( qgdnvy)
    call bl_deallocate ( qgdnvyf)
    call bl_deallocate ( qgdnvtmpy)

    call bl_deallocate ( qgdnvz)
    call bl_deallocate ( qgdnvzf)
    call bl_deallocate ( qgdnvtmpz1)
    call bl_deallocate ( qgdnvtmpz2)

    call bl_deallocate ( qxm)
    call bl_deallocate ( qxp)

    call bl_deallocate ( qmxy)
    call bl_deallocate ( qpxy)

    call bl_deallocate ( qmxz)
    call bl_deallocate ( qpxz)

    call bl_deallocate ( qym)
    call bl_deallocate ( qyp)

    call bl_deallocate ( qmyx)
    call bl_deallocate ( qpyx)

    call bl_deallocate ( qmyz)
    call bl_deallocate ( qpyz)

    call bl_deallocate ( qzm)
    call bl_deallocate ( qzp)

    call bl_deallocate ( qxl)
    call bl_deallocate ( qxr)
    call bl_deallocate ( qyl)
    call bl_deallocate ( qyr)
    call bl_deallocate ( qzl)
    call bl_deallocate ( qzr)

    call bl_deallocate ( qmzx)
    call bl_deallocate ( qpzx)

    call bl_deallocate ( qmzy)
    call bl_deallocate ( qpzy)

    call bl_deallocate ( fx)
    call bl_deallocate ( fy)
    call bl_deallocate ( fz)

    call bl_deallocate ( fxy)
    call bl_deallocate ( fxz)

    call bl_deallocate ( fyx)
    call bl_deallocate ( fyz)

    call bl_deallocate ( fzx)
    call bl_deallocate ( fzy)

#ifdef RADIATION
    call bl_deallocate (rfx)
    call bl_deallocate (rfy)
    call bl_deallocate (rfz)

    call bl_deallocate(rfxy)
    call bl_deallocate(rfxz)

    call bl_deallocate(rfyx)
    call bl_deallocate(rfyz)

    call bl_deallocate(rfzx)
    call bl_deallocate(rfxy)
#endif

    if (ppm_type .gt. 0) then
       call bl_deallocate ( Ip)
       call bl_deallocate ( Im)

       call bl_deallocate ( Ip_src)
       call bl_deallocate ( Im_src)

       call bl_deallocate ( Ip_gc)
       call bl_deallocate ( Im_gc)
    else
       call bl_deallocate ( dqx)
       call bl_deallocate ( dqy)
       call bl_deallocate ( dqz)
    end if

    call bl_deallocate(shk)

  end subroutine umeth

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine consup(uin, uin_lo, uin_hi, &
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
                    qx, qx_lo, qx_hi, &
                    qy, qy_lo, qy_hi, &
                    qz, qz_lo, qz_hi, &
                    area1, area1_lo, area1_hi, &
                    area2, area2_lo, area2_hi, &
                    area3, area3_lo, area3_hi, &
                    vol,vol_lo,vol_hi, &
                    div, lo, hi, dx, dt, &
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost, &
                    verbose)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, NGDNV, NQ, &
#ifdef RADIATION
                                   fspace_type, comoving, &
                                   GDPRES, GDU, GDV, GDW, GDLAMS, GDERADS, &
#endif
                                   track_grid_losses, limit_fluxes_on_small_dens
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens, normalize_species_fluxes, calc_pdivu
    use castro_util_module, only : position, linear_to_angular_momentum
    use prob_params_module, only : domlo_level, domhi_level, center
    use amrinfo_module, only : amr_level
#ifdef RADIATION
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
#endif
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) ::  uout_lo(3),  uout_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)

#ifdef RADIATION
    integer, intent(in) :: Erout_lo(3), Erout_hi(3)
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
    integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
    integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
    integer, intent(inout) :: nstep_fsp
#endif

    integer, intent(in) :: verbose


    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)
    real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt)        , intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt)        , intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt)        , intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt)        , intent(in) :: div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(rt)        , intent(in) :: dx(3), dt

#ifdef RADIATION
    real(rt)          Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt)         Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),Erout_lo(3):Erout_hi(3),0:ngroups-1)
    real(rt)         radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),radflux1_lo(3):radflux1_hi(3),0:ngroups-1)
    real(rt)         radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),radflux2_lo(3):radflux2_hi(3),0:ngroups-1)
    real(rt)         radflux3(radflux3_lo(1):radflux3_hi(1),radflux3_lo(2):radflux3_hi(2),radflux3_lo(3):radflux3_hi(3),0:ngroups-1)
#endif

    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    real(rt)         :: div1, volinv
    integer          :: i, j, g, k, n
    integer          :: domlo(3), domhi(3)
    real(rt)         :: loc(3), ang_mom(3)

#ifdef RADIATION
    real(rt)        , dimension(0:ngroups-1) :: Erscale
    real(rt)        , dimension(0:ngroups-1) :: ustar, af
    real(rt)         :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp
    real(rt)         :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    real(rt)         :: Gf1E(3)
    real(rt)         :: ux, uy, uz, divu, lamc, Egdc
    real(rt)         :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    real(rt)         :: dprdx, dprdy, dprdz, ek1, ek2, dek
    real(rt)         :: urho_new
    real(rt)         :: umx_new1, umy_new1, umz_new1
    real(rt)         :: umx_new2, umy_new2, umz_new2
#endif
    real(rt)        , pointer:: pdivu(:,:,:)

#ifdef RADIATION
    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if
#endif

    call bl_allocate(pdivu, lo, hi)

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, area1_lo, area1_hi, &
                    qy, qy_lo, qy_hi, &
                    area2, area2_lo, area2_hi, &
                    qz, qz_lo, qz_hi, &
                    area3, area3_lo, area3_hi, &
                    vol, vol_lo, vol_hi, &
                    dx, pdivu, lo, hi)

    do n = 1, NVAR

       if ( n == UTEMP ) then

          flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
          flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
          flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO

#ifdef SHOCK_VAR
       else if ( n == USHK ) then

          flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
          flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
          flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

       else

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   div1 = FOURTH*(div(i,j,k) + div(i,j+1,k) + div(i,j,k+1) + div(i,j+1,k+1))
                   div1 = difmag*min(ZERO,div1)

                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx(1) * div1 * (uin(i,j,k,n)-uin(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                   div1 = difmag*min(ZERO,div1)

                   flux2(i,j,k,n) = flux2(i,j,k,n) + dx(2) * div1 * (uin(i,j,k,n)-uin(i,j-1,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag*min(ZERO,div1)

                   flux3(i,j,k,n) = flux3(i,j,k,n) + dx(3) * div1 * (uin(i,j,k,n)-uin(i,j,k-1,n))
                enddo
             enddo
          enddo

       endif

    enddo

#ifdef RADIATION
    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                div1 = FOURTH*(div(i,j,k) + div(i,j+1,k) + div(i,j,k+1) + div(i,j+1,k+1))
                div1 = difmag*min(ZERO, div1)

                radflux1(i,j,k,g) = radflux1(i,j,k,g) + dx(1)*div1*(Erin(i,j,k,g)-Erin(i-1,j,k,g))
             enddo
          enddo
       enddo
    enddo

    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                div1 = difmag*min(ZERO, div1)

                radflux2(i,j,k,g) = radflux2(i,j,k,g) + dx(2)*div1*(Erin(i,j,k,g)-Erin(i,j-1,k,g))
             enddo
          enddo
       enddo
    enddo

    do g=0,ngroups-1
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                div1 = difmag*min(ZERO, div1)

                radflux3(i,j,k,g) = radflux3(i,j,k,g) + dx(3)*div1*(Erin(i,j,k,g)-Erin(i,j,k-1,g))
             enddo
          enddo
       enddo
    enddo
#endif

    if (limit_fluxes_on_small_dens == 1) then
       call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
                                             q,q_lo,q_hi, &
                                             vol,vol_lo,vol_hi, &
                                             flux1,flux1_lo,flux1_hi, &
                                             area1,area1_lo,area1_hi, &
                                             flux2,flux2_lo,flux2_hi, &
                                             area2,area2_lo,area2_hi, &
                                             flux3,flux3_lo,flux3_hi, &
                                             area3,area3_lo,area3_hi, &
                                             lo,hi,dt,dx)

    endif

    call normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                  flux2,flux2_lo,flux2_hi, &
                                  flux3,flux3_lo,flux3_hi, &
                                  lo,hi)

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                       flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                       flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

#ifdef HYBRID_MOMENTUM
    call add_hybrid_advection_source(lo, hi, dt, &
                                     update, uout_lo, uout_hi, &
                                     qx, qx_lo, qx_hi, &
                                     qy, qy_lo, qy_hi, &
                                     qz, qz_lo, qz_hi)
#endif


#ifdef RADIATION
    ! radiation energy update.  For the moment, we actually update things
    ! fully here, instead of creating a source term for the update
    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Erout(i,j,k,g) = Erin(i,j,k,g) + dt * &
                     ( radflux1(i,j,k,g) * area1(i,j,k) - radflux1(i+1,j,k,g) * area1(i+1,j,k) + &
                       radflux2(i,j,k,g) * area2(i,j,k) - radflux2(i,j+1,k,g) * area2(i,j+1,k) + &
                       radflux3(i,j,k,g) * area3(i,j,k) - radflux3(i,j,k+1,g) * area3(i,j,k+1) ) / vol(i,j,k)
             enddo
          enddo
       enddo
    enddo

    ! add radiation force terms (the hydro pressure gradient is always
    ! included in the momentum flux in 3-d)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             dprdx = ZERO
             dprdy = ZERO
             dprdz = ZERO
             do g=0,ngroups-1
                lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + &
                        qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.e0_rt
                dprdx = dprdx + lamc*(qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx(1)
                dprdy = dprdy + lamc*(qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx(2)
                dprdz = dprdz + lamc*(qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx(3)
             end do

             ! we now want to compute the change in the kinetic energy -- we
             ! base this off of uout, since that has the source terms in it.
             ! But note that this does not have the fluxes (since we are
             ! using update)

             urho_new = uout(i,j,k,URHO) + dt * update(i,j,k,URHO)

             ! this update includes the hydro fluxes and grad{p} from hydro
             umx_new1 = uout(i,j,k,UMX) + dt * update(i,j,k,UMX)
             umy_new1 = uout(i,j,k,UMY) + dt * update(i,j,k,UMY)
             umz_new1 = uout(i,j,k,UMZ) + dt * update(i,j,k,UMZ)

             ek1 = (umx_new1**2 + umy_new1**2 + umz_new1**2) / (TWO*urho_new)

             ! now add the radiation pressure gradient
             update(i,j,k,UMX) = update(i,j,k,UMX) - dprdx
             update(i,j,k,UMY) = update(i,j,k,UMY) - dprdy
             update(i,j,k,UMZ) = update(i,j,k,UMZ) - dprdz

             umx_new2 = umx_new1 - dt * dprdx
             umy_new2 = umy_new1 - dt * dprdy
             umz_new2 = umz_new1 - dt * dprdz

             ek2 = (umx_new2**2 + umy_new2**2 + umz_new2**2) / (TWO*urho_new)

             dek = ek2 - ek1

             ! the update is a source / dt, so scale accordingly
             update(i,j,k,UEDEN) = update(i,j,k,UEDEN) + dek/dt

             if (.not. comoving) then ! mixed-frame (single group only)
                Erout(i,j,k,0) = Erout(i,j,k,0) - dek
             end if

          end do
       end do
    end do

    ! Add radiation source terms
    if (comoving) then
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                ux = HALF*(qx(i,j,k,GDU) + qx(i+1,j,k,GDU))
                uy = HALF*(qy(i,j,k,GDV) + qy(i,j+1,k,GDV))
                uz = HALF*(qz(i,j,k,GDW) + qz(i,j,k+1,GDW))

                dudx(1) = (qx(i+1,j,k,GDU) - qx(i,j,k,GDU))/dx(1)
                dudx(2) = (qx(i+1,j,k,GDV) - qx(i,j,k,GDV))/dx(1)
                dudx(3) = (qx(i+1,j,k,GDW) - qx(i,j,k,GDW))/dx(1)

                dudy(1) = (qy(i,j+1,k,GDU) - qy(i,j,k,GDU))/dx(2)
                dudy(2) = (qy(i,j+1,k,GDV) - qy(i,j,k,GDV))/dx(2)
                dudy(3) = (qy(i,j+1,k,GDW) - qy(i,j,k,GDW))/dx(2)

                dudz(1) = (qz(i,j,k+1,GDU) - qz(i,j,k,GDU))/dx(3)
                dudz(2) = (qz(i,j,k+1,GDV) - qz(i,j,k,GDV))/dx(3)
                dudz(3) = (qz(i,j,k+1,GDW) - qz(i,j,k,GDW))/dx(3)

                divu = dudx(1) + dudy(2) + dudz(3)

                ! Note that for single group, fspace_type is always 1
                do g=0, ngroups-1

                   nhat(1) = (qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx(1)
                   nhat(2) = (qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx(2)
                   nhat(3) = (qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx(3)

                   GnDotu(1) = dot_product(nhat, dudx)
                   GnDotu(2) = dot_product(nhat, dudy)
                   GnDotu(3) = dot_product(nhat, dudz)

                   nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.e-50_rt)

                   lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                           qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + &
                           qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.e0_rt
                   Eddf = Edd_factor(lamc)
                   f1 = (ONE-Eddf)*HALF
                   f2 = (3.e0_rt*Eddf-ONE)*HALF
                   af(g) = -(f1*divu + f2*nnColonDotGu)

                   if (fspace_type .eq. 1) then
                      Eddfxp = Edd_factor(qx(i+1,j  ,k  ,GDLAMS+g))
                      Eddfxm = Edd_factor(qx(i  ,j  ,k  ,GDLAMS+g))
                      Eddfyp = Edd_factor(qy(i  ,j+1,k  ,GDLAMS+g))
                      Eddfym = Edd_factor(qy(i  ,j  ,k  ,GDLAMS+g))
                      Eddfzp = Edd_factor(qz(i  ,j  ,k+1,GDLAMS+g))
                      Eddfzm = Edd_factor(qz(i  ,j  ,k  ,GDLAMS+g))

                      f1xp = HALF*(ONE-Eddfxp)
                      f1xm = HALF*(ONE-Eddfxm)
                      f1yp = HALF*(ONE-Eddfyp)
                      f1ym = HALF*(ONE-Eddfym)
                      f1zp = HALF*(ONE-Eddfzp)
                      f1zm = HALF*(ONE-Eddfzm)

                      Gf1E(1) = (f1xp*qx(i+1,j,k,GDERADS+g) - f1xm*qx(i,j,k,GDERADS+g)) / dx(1)
                      Gf1E(2) = (f1yp*qy(i,j+1,k,GDERADS+g) - f1ym*qy(i,j,k,GDERADS+g)) / dx(2)
                      Gf1E(3) = (f1zp*qz(i,j,k+1,GDERADS+g) - f1zm*qz(i,j,k,GDERADS+g)) / dx(3)

                      Egdc = (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) &
                           +  qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g) &
                           +  qz(i,j,k,GDERADS+g) + qz(i,j,k+1,GDERADS+g) ) / 6.e0_rt

                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)+uz*Gf1E(3)) &
                           - dt*f2*Egdc*nnColonDotGu
                   end if

                end do

                if (ngroups.gt.1) then
                   ustar = Erout(i,j,k,:) / Erscale
                   call advect_in_fspace(ustar, af, dt, nstep_fsp)
                   Erout(i,j,k,:) = ustar * Erscale
                end if
             enddo
          enddo
       enddo
    endif
#endif

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo

#ifdef RADIATION
    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                radflux1(i,j,k,g) = dt * radflux1(i,j,k,g) * area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                radflux2(i,j,k,g) = dt * radflux2(i,j,k,g) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do g = 0, ngroups-1
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                radflux3(i,j,k,g) = dt * radflux3(i,j,k,g) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif


    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)

       if (lo(3) .le. domlo(3) .and. hi(3) .ge. domlo(3)) then

          k = domlo(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccz=.false.)

                mass_lost = mass_lost - flux3(i,j,k,URHO)
                xmom_lost = xmom_lost - flux3(i,j,k,UMX)
                ymom_lost = ymom_lost - flux3(i,j,k,UMY)
                zmom_lost = zmom_lost - flux3(i,j,k,UMZ)
                eden_lost = eden_lost - flux3(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux3(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(3) .le. domhi(3) .and. hi(3) .ge. domhi(3)) then

          k = domhi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccz=.false.)

                mass_lost = mass_lost + flux3(i,j,k,URHO)
                xmom_lost = xmom_lost + flux3(i,j,k,UMX)
                ymom_lost = ymom_lost + flux3(i,j,k,UMY)
                zmom_lost = zmom_lost + flux3(i,j,k,UMZ)
                eden_lost = eden_lost + flux3(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux3(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

       if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

          j = domlo(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost - flux2(i,j,k,URHO)
                xmom_lost = xmom_lost - flux2(i,j,k,UMX)
                ymom_lost = ymom_lost - flux2(i,j,k,UMY)
                zmom_lost = zmom_lost - flux2(i,j,k,UMZ)
                eden_lost = eden_lost - flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

          j = domhi(2) + 1
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost + flux2(i,j,k,URHO)
                xmom_lost = xmom_lost + flux2(i,j,k,UMX)
                ymom_lost = ymom_lost + flux2(i,j,k,UMY)
                zmom_lost = zmom_lost + flux2(i,j,k,UMZ)
                eden_lost = eden_lost + flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost - flux1(i,j,k,URHO)
                xmom_lost = xmom_lost - flux1(i,j,k,UMX)
                ymom_lost = ymom_lost - flux1(i,j,k,UMY)
                zmom_lost = zmom_lost - flux1(i,j,k,UMZ)
                eden_lost = eden_lost - flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost + flux1(i,j,k,URHO)
                xmom_lost = xmom_lost + flux1(i,j,k,UMX)
                ymom_lost = ymom_lost + flux1(i,j,k,UMY)
                zmom_lost = zmom_lost + flux1(i,j,k,UMZ)
                eden_lost = eden_lost + flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

    endif

    call bl_deallocate(pdivu)

  end subroutine consup

end module ctu_advection_module
