module ctu_advection_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public umeth

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

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NQ, NVAR, QPRES, QRHO, QU, QW, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   ppm_type, &
                                   use_pslope, ppm_temp_fix, &
                                   hybrid_riemann
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
    use trace_module, only : tracexy, tracez
    use ppm_module, only : ppm_reconstruct, ppm_int_profile
    use slope_module, only : uslope, pslope
    use actual_network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use riemann_module, only: cmpflx
    use amrex_constants_module
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

    logical :: source_nonzero(QVAR)

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

    ! preprocess the sources -- we don't want to trace under a source that is empty
    if (ppm_type > 0) then
       do n = 1, QVAR
          if (minval(srcQ(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,n)) == ZERO .and. &
              maxval(srcQ(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,n)) == ZERO) then
             source_nonzero(n) = .false.
          else
             source_nonzero(n) = .true.
          endif
       enddo
    endif

    do k3d = lo(3)-1, hi(3)+1

       ! Swap pointers to levels
       kt = km
       km = kc
       kc = kt

       if (ppm_type > 0) then

          do n = 1, NQ
             call ppm_reconstruct(q, qd_lo, qd_hi, NQ, n, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

             call ppm_int_profile(q, qd_lo, qd_hi, NQ, n, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  Ip, Im, It_lo, It_hi, NQ, n, &
                                  lo(1), lo(2), hi(1), hi(2), dx, dt, k3d, kc)
          end do

          ! source terms
          do n = 1, QVAR
             if (source_nonzero(n)) then
                call ppm_reconstruct(srcQ, src_lo, src_hi, QVAR, n, &
                                     flatn, qd_lo, qd_hi, &
                                     sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                     lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

                call ppm_int_profile(srcQ, src_lo, src_hi, QVAR, n, &
                                     q, qd_lo, qd_hi, &
                                     qaux, qa_lo, qa_hi, &
                                     sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                     Ip_src, Im_src, It_lo, It_hi, QVAR, n, &
                                     lo(1), lo(2), hi(1), hi(2), dx, dt, k3d, kc)
             else
                Ip_src(It_lo(1):It_hi(1),It_lo(2):It_hi(2),kc,:,:,n) = ZERO
                Im_src(It_lo(1):It_hi(1),It_lo(2):It_hi(2),kc,:,:,n) = ZERO
             endif

          enddo

          ! this probably doesn't support radiation
          if (ppm_temp_fix /= 1) then
             call ppm_reconstruct(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

             call ppm_int_profile(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  Ip_gc, Im_gc, It_lo, It_hi, 1, 1, &
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
#ifndef AMREX_USE_CUDA
          call amrex_error("ppm_type <=0 is not supported in with radiation")
#endif          
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


  subroutine reset_edge_state_thermo(qedge, qd_lo, qd_hi, ii, jj, kk)

  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer, intent(in) :: ii, jj, kk
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    real(rt)        , intent(inout) :: qedge(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)

    logical :: reset
    type (eos_t) :: eos_state

    reset = .false.

    if (transverse_reset_rhoe == 1) then
       ! if we are still negative, then we need to reset
       if (qedge(ii,jj,kk,QREINT) < ZERO) then
          reset = .true.

          eos_state % rho = qedge(ii,jj,kk,QRHO)
          eos_state % T = small_temp
          eos_state % xn(:) = qedge(ii,jj,kk,QFS:QFS-1+nspec)
          eos_state % aux(:) = qedge(ii,jj,kk,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          qedge(ii,jj,kk,QREINT) = qedge(ii,jj,kk,QRHO)*eos_state % e
          qedge(ii,jj,kk,QPRES) = eos_state % p
       endif

    end if

    if (ppm_predict_gammae == 0 ) then

       if (transverse_use_eos == 1) then
          eos_state % rho = qedge(ii,jj,kk,QRHO)
          eos_state % e   = qedge(ii,jj,kk,QREINT) / qedge(ii,jj,kk,QRHO)
          eos_state % T   = small_temp
          eos_state % xn  = qedge(ii,jj,kk,QFS:QFS+nspec-1)
          eos_state % aux = qedge(ii,jj,kk,QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          qedge(ii,jj,kk,QREINT) = eos_state % e * eos_state % rho
          qedge(ii,jj,kk,QPRES) = max(eos_state % p, small_pres)
       end if

    else
       if (reset) then
          ! recompute the p edge state from this and (rho e), since we reset
          ! qreint  (actually, is this code even necessary?)
          qedge(ii,jj,kk,QPRES) = qedge(ii,jj,kk,QREINT)*(qedge(ii,jj,kk,QGAME)-ONE)
          qedge(ii,jj,kk,QPRES) = max(qedge(ii,jj,kk,QPRES), small_pres)
       end if
    end if

  end subroutine reset_edge_state_thermo


  !===========================================================================
  ! transx1
  !===========================================================================
  subroutine transx1(qym, qymo, qyp, qypo, qd_lo, qd_hi, &
                     qaux, qa_lo, qa_hi, &
                     fx, &
#ifdef RADIATION
                     rfx, &
#endif
                     fx_lo, fx_hi, &
                     qx, qx_lo, qx_hi, &
                     cdtdx, ilo, ihi, jlo, jhi, kc, k3d)

    ! Note that what we call ilo here is ilo = lo(1)
    ! Note that what we call ihi here is ihi = hi(1)
    ! Note that what we call jlo here is jlo = lo(2) - 1
    ! Note that what we call jhi here is jhi = hi(2) + 1


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer ilo,ihi,jlo,jhi,kc,k3d

#ifdef RADIATION
    real(rt)         rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
#endif

    real(rt)          qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qymo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qypo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)         fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt)         qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt)         cdtdx

    integer i, j, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         rrry, rrly
    real(rt)         rury, ruly
    real(rt)         rvry, rvly
    real(rt)         rwry, rwly
    real(rt)         ekenry, ekenly
    real(rt)         rery, rely
    real(rt)         rrnewry, rrnewly
    real(rt)         runewry, runewly
    real(rt)         rvnewry, rvnewly
    real(rt)         rwnewry, rwnewly
    real(rt)         renewry, renewly
    real(rt)         pnewry, pnewly
    real(rt)         rhoekenry, rhoekenly
    real(rt)         compn, compu
    real(rt)         pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state


    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transverse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             if (j >= jlo+1) then
                rr = qyp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qyp(i,j,kc,nqp) - compn
                qypo(i,j,kc,nqp) = compu/rrnew
             end if

             if (j <= jhi-1) then
                rr = qym(i,j+1,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qym(i,j+1,kc,nqp) - compn
                qymo(i,j+1,kc,nqp) = compu/rrnew
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to y-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          pggp  = qx(i+1,j,kc,GDPRES)
          pggm  = qx(i  ,j,kc,GDPRES)
          ugp  = qx(i+1,j,kc,GDU   )
          ugm  = qx(i  ,j,kc,GDU   )
          gegp = qx(i+1,j,kc,GDGAME)
          gegm = qx(i  ,j,kc,GDGAME)

#ifdef RADIATION
          lambda = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
          ugc = HALF*(ugp+ugm)
          ergp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,k3d,QGAMCG)
#else
          gamc = qaux(i,j,k3d,QGAMC)
#endif

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          luge = ugc * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          if (j >= jlo+1) then

             ! Convert to conservation form
             rrry = qyp(i,j,kc,QRHO)
             rury = rrry*qyp(i,j,kc,QU)
             rvry = rrry*qyp(i,j,kc,QV)
             rwry = rrry*qyp(i,j,kc,QW)
             ekenry = HALF*rrry* &
                  (qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
             rery = qyp(i,j,kc,QREINT) + ekenry
#ifdef RADIATION
             err  = qyp(i,j,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewry = rury - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
             rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
#ifdef RADIATION
             runewry = runewry + dmom
             renewry = renewry + dre
             ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewry < ZERO) then
                rrnewry = rrry
                runewry = rury
                rvnewry = rvry
                rwnewry = rwry
                renewry = rery
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             qypo(i,j,kc,QRHO) = rrnewry
             rhoinv = ONE/rrnewry
             qypo(i,j,kc,QU) = runewry*rhoinv
             qypo(i,j,kc,QV) = rvnewry*rhoinv
             qypo(i,j,kc,QW) = rwnewry*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
             qypo(i,j,kc,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe .eq. 1 .and. qypo(i,j,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qypo(i,j,kc,QREINT) = qyp(i,j,kc,QREINT) - &
                        cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qypo(i,j,kc,QPRES) = max(pnewry,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qypo(i,j,kc,QGAME) = qyp(i,j,kc,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qypo(i,j,kc,QPRES) = qypo(i,j,kc,QREINT)*(qypo(i,j,kc,QGAME)-ONE)
                   qypo(i,j,kc,QPRES) = max(qypo(i,j,kc,QPRES),small_pres)
                end if
             else
                qypo(i,j,kc,QPRES) = qyp(i,j,kc,QPRES)
                qypo(i,j,kc,QGAME) = qyp(i,j,kc,QGAME)
             endif

             call reset_edge_state_thermo(qypo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
             qypo(i,j,kc,qrad:qradhi) = ernewr(:)
             qypo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,kc,QPRES)
             qypo(i,j,kc,qreitot) = sum(qypo(i,j,kc,qrad:qradhi)) + qypo(i,j,kc,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------

          if (j <= jhi-1) then

             ! Convert to conservation form
             rrly = qym(i,j+1,kc,QRHO)
             ruly = rrly*qym(i,j+1,kc,QU)
             rvly = rrly*qym(i,j+1,kc,QV)
             rwly = rrly*qym(i,j+1,kc,QW)
             ekenly = HALF*rrly* &
                  (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
             rely = qym(i,j+1,kc,QREINT) + ekenly
#ifdef RADIATION
             erl  = qym(i,j+1,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
             rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
#ifdef RADIATION
             runewly = runewly + dmom
             renewly = renewly + dre
             ernewl  = erl(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewly < ZERO) then
                rrnewly = rrly
                runewly = ruly
                rvnewly = rvly
                rwnewly = rwly
                renewly = rely
#ifdef RADIATION
                ernewl  = erl(:)
#endif
                reset_state = .true.
             endif

             qymo(i,j+1,kc,QRHO) = rrnewly
             rhoinv = ONE/rrnewly
             qymo(i,j+1,kc,QU) = runewly*rhoinv
             qymo(i,j+1,kc,QV) = rvnewly*rhoinv
             qymo(i,j+1,kc,QW) = rwnewly*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
             qymo(i,j+1,kc,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary
                if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qymo(i,j+1,kc,QREINT) = qym(i,j+1,kc,QREINT) - &
                        cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qymo(i,j+1,kc,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qymo(i,j+1,kc,QGAME) = qym(i,j+1,kc,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qymo(i,j+1,kc,QPRES) = qymo(i,j+1,kc,QREINT)*(qymo(i,j+1,kc,QGAME)-ONE)
                   qymo(i,j+1,kc,QPRES) = max(qymo(i,j+1,kc,QPRES), small_pres)
                end if
             else
                qymo(i,j+1,kc,QPRES) = qym(i,j+1,kc,QPRES)
                qymo(i,j+1,kc,QGAME) = qym(i,j+1,kc,QGAME)
             endif

             call reset_edge_state_thermo(qymo, qd_lo, qd_hi, i, j+1, kc)

#ifdef RADIATION
             qymo(i,j+1,kc,qrad:qradhi) = ernewl(:)
             qymo(i,j+1,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j+1,kc,QPRES)
             qymo(i,j+1,kc,qreitot) = sum(qymo(i,j+1,kc,qrad:qradhi)) + qymo(i,j+1,kc,QREINT)
#endif
          endif
       enddo

    enddo

  end subroutine transx1


  !===========================================================================
  ! transx2
  !===========================================================================
  subroutine transx2(qzm, qzmo, qzp, qzpo, qd_lo, qd_hi, &
                     qaux, qa_lo, qa_hi, &
                     fx, &
#ifdef RADIATION
                     rfx, &
#endif
                     fx_lo, fx_hi, &
                     qx, qx_lo, qx_hi, &
                     cdtdx, ilo, ihi, jlo, jhi, kc, km, k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer ilo,ihi,jlo,jhi,kc,km,k3d

#ifdef RADIATION
    real(rt)         rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
#endif

    real(rt)          qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qzmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qzpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt)         fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt)         qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt)         cdtdx

    integer i, j, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         rrrz, rrlz
    real(rt)         rurz, rulz
    real(rt)         rvrz, rvlz
    real(rt)         rwrz, rwlz
    real(rt)         ekenrz, ekenlz
    real(rt)         rerz, relz
    real(rt)         rrnewrz, rrnewlz
    real(rt)         runewrz, runewlz
    real(rt)         rvnewrz, rvnewlz
    real(rt)         rwnewrz, rwnewlz
    real(rt)         renewrz, renewlz
    real(rt)         pnewrz, pnewlz
    real(rt)         rhoekenrz, rhoekenlz
    real(rt)         compn, compu
    real(rt)         pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nqp) - compn
             qzpo(i,j,kc,nqp) = compu/rrnew

             compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nqp) - compn
             qzmo(i,j,kc,nqp) = compu/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to z-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------

          pggp  = qx(i+1,j,kc,GDPRES)
          pggm  = qx(i  ,j,kc,GDPRES)
          ugp  = qx(i+1,j,kc,GDU   )
          ugm  = qx(i  ,j,kc,GDU   )
          gegp = qx(i+1,j,kc,GDGAME)
          gegm = qx(i  ,j,kc,GDGAME)
#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
          ugc = HALF*(ugp+ugm)
          ergp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,k3d,QGAMCG)
#else
          gamc = qaux(i,j,k3d,QGAMC)
#endif

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if
#endif


          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
#ifdef RADIATION
          err  = qzp(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
#ifdef RADIATION
          runewrz = runewrz + dmom
          renewrz = renewrz + dre
          ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewrz < ZERO) then
             rrnewrz = rrrz
             runewrz = rurz
             rvnewrz = rvrz
             rwnewrz = rwrz
             renewrz = rerz
#ifdef RADIATION
             ernewr  = err(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzpo(i,j,kc,QRHO) = rrnewrz
          rhoinv = ONE/rrnewrz
          qzpo(i,j,kc,QU) = runewrz*rhoinv
          qzpo(i,j,kc,QV) = rvnewrz*rhoinv
          qzpo(i,j,kc,QW) = rwnewrz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz

          if (.not. reset_state) then
             ! do the transverse terms for p, gamma, and rhoe, as necessary

             if (transverse_reset_rhoe == 1 .and. qzpo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzpo(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) - &
                     cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
             end if

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT)*(qzpo(i,j,kc,QGAME)-ONE)
                qzpo(i,j,kc,QPRES) = max(qzpo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzpo(i,j,kc,QPRES) = qzp(i,j,kc,QPRES)
             qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)
#endif

          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------

          pggp  = qx(i+1,j,km,GDPRES)
          pggm  = qx(i  ,j,km,GDPRES)
          ugp  = qx(i+1,j,km,GDU   )
          ugm  = qx(i  ,j,km,GDU   )
          gegp = qx(i+1,j,km,GDGAME)
          gegm = qx(i  ,j,km,GDGAME)
#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d-1,QLAMS:QLAMS+ngroups-1)
          ugc = HALF*(ugp+ugm)
          ergp = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,k3d-1,QGAMCG)
#else
          gamc = qaux(i,j,k3d-1,QGAMC)
#endif

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if
#endif

          ! Convert to conservation form
          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
#ifdef RADIATION
          erl  = qzm(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
          runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX) - fx(i,j,km,UMX))
          rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY) - fx(i,j,km,UMY))
          rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ) - fx(i,j,km,UMZ))
          renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN))
#ifdef RADIATION
          runewlz = runewlz + dmom
          renewlz = renewlz + dre
          ernewl  = erl(:) - cdtdx*(rfx(i+1,j,km,:) - rfx(i,j,km,:)) &
               +der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewlz < ZERO) then
             rrnewlz = rrlz
             runewlz = rulz
             rvnewlz = rvlz
             rwnewlz = rwlz
             renewlz = relz
#ifdef RADIATION
             ernewl  = erl(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzmo(i,j,kc,QRHO) = rrnewlz
          rhoinv = ONE/rrnewlz
          qzmo(i,j,kc,QU) = runewlz*rhoinv
          qzmo(i,j,kc,QV) = rvnewlz*rhoinv
          qzmo(i,j,kc,QW) = rwnewlz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz

          if (.not. reset_state) then
             ! do the transverse terms for p, gamma, and rhoe, as necessary

             if (transverse_reset_rhoe == 1 .and. qzmo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzmo(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) - &
                     cdtdx*(fx(i+1,j,km,UEINT) - fx(i,j,km,UEINT) + pav*du)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT)*(qzmo(i,j,kc,QGAME)-ONE)
                qzmo(i,j,kc,QPRES) = max(qzmo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzmo(i,j,kc,QPRES) = qzm(i,j,kc,QPRES)
             qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzmo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)
#endif

       enddo
    enddo

  end subroutine transx2


  !===========================================================================
  ! transy1
  !===========================================================================
  subroutine transy1(qxm, qxmo, qxp, qxpo, qd_lo, qd_hi, &
                     qaux, qa_lo, qa_hi, &
                     fy, &
#ifdef RADIATION
                     rfy, &
#endif
                     fy_lo, fy_hi, &
                     qy, qy_lo, qy_hi, &
                     cdtdy, ilo, ihi, jlo, jhi, kc, k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer ilo,ihi,jlo,jhi,kc,k3d

#ifdef RADIATION
    real(rt)         rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    real(rt)          qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qxmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qxpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)         fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt)         qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt)         cdtdy

    integer i, j, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         compn, compu
    real(rt)         rrrx, rrlx
    real(rt)         rurx, rulx
    real(rt)         rvrx, rvlx
    real(rt)         rwrx, rwlx
    real(rt)         ekenrx, ekenlx
    real(rt)         rerx, relx
    real(rt)         rrnewrx, rrnewlx
    real(rt)         runewrx, runewlx
    real(rt)         rvnewrx, rvnewlx
    real(rt)         rwnewrx, rwnewlx
    real(rt)         renewrx, renewlx
    real(rt)         pnewrx, pnewlx
    real(rt)         rhoekenrx, rhoekenlx
    real(rt)         pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             if (i >= ilo+1) then
                rr = qxp(i,j,kc,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                compu = rr*qxp(i,j,kc,nqp) - compn
                qxpo(i,j,kc,nqp) = compu/rrnew
             end if

             if (i <= ihi-1) then
                rr = qxm(i+1,j,kc,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                compu = rr*qxm(i+1,j,kc,nqp) - compn
                qxmo(i+1,j,kc,nqp) = compu/rrnew
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to x-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          pggp  = qy(i,j+1,kc,GDPRES)
          pggm  = qy(i,j  ,kc,GDPRES)
          ugp  = qy(i,j+1,kc,GDV   )
          ugm  = qy(i,j  ,kc,GDV   )
          gegp = qy(i,j+1,kc,GDGAME)
          gegm = qy(i,j  ,kc,GDGAME)
#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
          ugc = HALF*(ugp+ugm)
          ergp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,k3d,QGAMCG)
#else
          gamc = qaux(i,j,k3d,QGAMC)
#endif

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          luge = ugc * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------

          if (i >= ilo+1) then
             ! Convert to conservation form
             rrrx = qxp(i,j,kc,QRHO)
             rurx = rrrx*qxp(i,j,kc,QU)
             rvrx = rrrx*qxp(i,j,kc,QV)
             rwrx = rrrx*qxp(i,j,kc,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
                  + qxp(i,j,kc,QW)**2)
             rerx = qxp(i,j,kc,QREINT) + ekenrx
#ifdef RADIATION
             err  = qxp(i,j,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
             rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
             rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
             renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))
#ifdef RADIATION
             rvnewrx = rvnewrx + dmom
             renewrx = renewrx + dre
             ernewr = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewrx < ZERO) then
                rrnewrx = rrrx
                runewrx = rurx
                rvnewrx = rvrx
                rwnewrx = rwrx
                renewrx = rerx
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             qxpo(i,j,kc,QRHO) = rrnewrx
             rhoinv = ONE/rrnewrx
             qxpo(i,j,kc,QU) = runewrx*rhoinv
             qxpo(i,j,kc,QV) = rvnewrx*rhoinv
             qxpo(i,j,kc,QW) = rwnewrx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
             qxpo(i,j,kc,QREINT) = renewrx - rhoekenrx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxpo(i,j,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qxpo(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) - &
                        cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qxpo(i,j,kc,QPRES) = max(pnewrx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxpo(i,j,kc,QGAME) = qxp(i,j,kc,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxpo(i,j,kc,QPRES) = qxpo(i,j,kc,QREINT)*(qxpo(i,j,kc,QGAME)-ONE)
                   qxpo(i,j,kc,QPRES) = max(qxpo(i,j,kc,QPRES), small_pres)
                endif
             else
                qxpo(i,j,kc,QPRES) = qxp(i,j,kc,QPRES)
                qxpo(i,j,kc,QGAME) = qxp(i,j,kc,QGAME)
             endif

             call reset_edge_state_thermo(qxpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
             qxpo(i,j,kc,qrad:qradhi) = ernewr(:)
             qxpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,kc,QPRES)
             qxpo(i,j,kc,qreitot) = sum(qxpo(i,j,kc,qrad:qradhi)) + qxpo(i,j,kc,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i <= ihi-1) then
             ! Convert to conservation form
             rrlx = qxm(i+1,j,kc,QRHO)
             rulx = rrlx*qxm(i+1,j,kc,QU)
             rvlx = rrlx*qxm(i+1,j,kc,QV)
             rwlx = rrlx*qxm(i+1,j,kc,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
                  + qxm(i+1,j,kc,QW)**2)
             relx = qxm(i+1,j,kc,QREINT) + ekenlx
#ifdef RADIATION
             erl  = qxm(i+1,j,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
             rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
             rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
             renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN)- fy(i,j,kc,UEDEN))
#ifdef RADIATION
             rvnewlx = rvnewlx + dmom
             renewlx = renewlx + dre
             ernewl  = erl(:) - cdtdy*(rfy(i,j+1,kc,:)- rfy(i,j,kc,:)) + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewlx < ZERO) then
                rrnewlx = rrlx
                runewlx = rulx
                rvnewlx = rvlx
                rwnewlx = rwlx
                renewlx = relx
#ifdef RADIATION
                ernewl  = erl(:)
#endif
                reset_state = .true.
             endif

             qxmo(i+1,j,kc,QRHO) = rrnewlx
             rhoinv = ONE/rrnewlx
             qxmo(i+1,j,kc,QU) = runewlx*rhoinv
             qxmo(i+1,j,kc,QV) = rvnewlx*rhoinv
             qxmo(i+1,j,kc,QW) = rwnewlx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
             qxmo(i+1,j,kc,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) - &
                        cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qxmo(i+1,j,kc,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo(i+1,j,kc,QGAME) = qxm(i+1,j,kc,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo(i+1,j,kc,QPRES) = qxmo(i+1,j,kc,QREINT)*(qxmo(i+1,j,kc,QGAME)-ONE)
                   qxmo(i+1,j,kc,QPRES) = max(qxmo(i+1,j,kc,QPRES), small_pres)
                endif
             else
                qxmo(i+1,j,kc,QPRES) = qxm(i+1,j,kc,QPRES)
                qxmo(i+1,j,kc,QGAME) = qxm(i+1,j,kc,QGAME)
             endif

             call reset_edge_state_thermo(qxmo, qd_lo, qd_hi, i+1, j, kc)

#ifdef RADIATION
             qxmo(i+1,j,kc,qrad:qradhi) = ernewl(:)
             qxmo(i+1,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i+1,j,kc,QPRES)
             qxmo(i+1,j,kc,qreitot) = sum(qxmo(i+1,j,kc,qrad:qradhi)) + qxmo(i+1,j,kc,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transy1


  !===========================================================================
  ! transy2
  !===========================================================================
  subroutine transy2(qzm, qzmo, qzp, qzpo, qd_lo, qd_hi, &
                     qaux, qa_lo, qa_hi, &
                     fy, &
#ifdef RADIATION
                     rfy, &
#endif
                     fy_lo, fy_hi, &
                     qy, qy_lo, qy_hi, &
                     cdtdy, ilo, ihi, jlo, jhi, kc, km, k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer ilo,ihi,jlo,jhi,kc,km,k3d

#ifdef RADIATION
    real(rt)         rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    real(rt)          qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qzmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qzpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)         fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt)         qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt)         cdtdy

    integer i, j, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         compn, compu
    real(rt)         rrrz, rrlz
    real(rt)         rurz, rulz
    real(rt)         rvrz, rvlz
    real(rt)         rwrz, rwlz
    real(rt)         ekenrz, ekenlz
    real(rt)         rerz, relz
    real(rt)         rrnewrz, rrnewlz
    real(rt)         runewrz, runewlz
    real(rt)         rvnewrz, rvnewlz
    real(rt)         rwnewrz, rwnewlz
    real(rt)         renewrz, renewlz
    real(rt)         pnewrz, pnewlz
    real(rt)         rhoekenrz, rhoekenlz
    real(rt)         pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nqp) - compn
             qzpo(i,j,kc,nqp) = compu/rrnew

             compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nqp) - compn
             qzmo(i,j,kc,nqp) = compu/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to z-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          !-------------------------------------------------------------------
          ! qzpo states
          !-------------------------------------------------------------------

          pggp  = qy(i,j+1,kc,GDPRES)
          pggm  = qy(i,j  ,kc,GDPRES)
          ugp  = qy(i,j+1,kc,GDV   )
          ugm  = qy(i,j  ,kc,GDV   )
          gegp = qy(i,j+1,kc,GDGAME)
          gegm = qy(i,j  ,kc,GDGAME)
#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
          ugc = HALF*(ugp+ugm)
          ergp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          gamc = qaux(i,j,k3d,QGAMCG)
#else
          gamc = qaux(i,j,k3d,QGAMC)
#endif


#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
               + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
#ifdef RADIATION
          err  = qzp(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))
#ifdef RADIATION
          rvnewrz = rvnewrz + dmom
          renewrz = renewrz + dre
          ernewr  = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewrz < ZERO) then
             rrnewrz = rrrz
             runewrz = rurz
             rvnewrz = rvrz
             rwnewrz = rwrz
             renewrz = rerz
#ifdef RADIATION
             ernewr  = err(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzpo(i,j,kc,QRHO) = rrnewrz
          rhoinv = ONE/rrnewrz
          qzpo(i,j,kc,QU) = runewrz*rhoinv
          qzpo(i,j,kc,QV) = rvnewrz*rhoinv
          qzpo(i,j,kc,QW) = rwnewrz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qzpo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzpo(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) - &
                     cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT)*(qzpo(i,j,kc,QGAME)-ONE)
                qzpo(i,j,kc,QPRES) = max(qzpo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzpo(i,j,kc,QPRES) = qzp(i,j,kc,QPRES)
             qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)
#endif

          !-------------------------------------------------------------------
          ! qzmo states
          !-------------------------------------------------------------------

          pggp  = qy(i,j+1,km,GDPRES)
          pggm  = qy(i,j  ,km,GDPRES)
          ugp  = qy(i,j+1,km,GDV   )
          ugm  = qy(i,j  ,km,GDV   )
          gegp = qy(i,j+1,km,GDGAME)
          gegm = qy(i,j  ,km,GDGAME)
#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d-1,QLAMS:QLAMS+ngroups-1)
          ugc = HALF*(ugp+ugm)
          ergp = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,k3d-1,QGAMCG)
#else
          gamc = qaux(i,j,k3d-1,QGAMC)
#endif

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          ! Convert to conservation form
          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
               + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
#ifdef RADIATION
          erl  = qzm(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
          runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX) - fy(i,j,km,UMX))
          rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY) - fy(i,j,km,UMY))
          rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ) - fy(i,j,km,UMZ))
          renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN)- fy(i,j,km,UEDEN))
#ifdef RADIATION
          rvnewlz = rvnewlz + dmom
          renewlz = renewlz + dre
          ernewl  = erl(:) - cdtdy*(rfy(i,j+1,km,:)- rfy(i,j,km,:)) &
               + der
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewlz < ZERO) then
             rrnewlz = rrlz
             runewlz = rulz
             rvnewlz = rvlz
             rwnewlz = rwlz
             renewlz = relz
#ifdef RADIATION
             ernewl  = erl(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzmo(i,j,kc,QRHO) = rrnewlz
          rhoinv = ONE/rrnewlz
          qzmo(i,j,kc,QU) = runewlz*rhoinv
          qzmo(i,j,kc,QV) = rvnewlz*rhoinv
          qzmo(i,j,kc,QW) = rwnewlz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qzmo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzmo(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) - &
                     cdtdy*(fy(i,j+1,km,UEINT) - fy(i,j,km,UEINT) + pav*du)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT)*(qzmo(i,j,kc,QGAME)-ONE)
                qzmo(i,j,kc,QPRES) = max(qzmo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzmo(i,j,kc,QPRES) = qzm(i,j,kc,QPRES)
             qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzmo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)
#endif

       enddo
    enddo

  end subroutine transy2


  !===========================================================================
  ! transz
  !===========================================================================
  subroutine transz(qxm,qxmo,qxp,qxpo,qym,qymo,qyp,qypo,qd_lo,qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    fz, &
#ifdef RADIATION
                    rfz, &
#endif
                    fz_lo,fz_hi, &
                    qz,qz_lo,qz_hi, &
                    cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fz_lo(3),fz_hi(3)
    integer :: qz_lo(3),qz_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    real(rt)         rfz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    real(rt)          qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qxmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qxpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qymo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qypo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)         fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt)         qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt)         cdtdz

    integer n, nqp, i, j, ipassive

    real(rt)         rrnew, rr
    real(rt)         compn, compu
    real(rt)         rrrx, rrry, rrlx, rrly
    real(rt)         rurx, rury, rulx, ruly
    real(rt)         rvrx, rvry, rvlx, rvly
    real(rt)         rwrx, rwry, rwlx, rwly
    real(rt)         ekenrx, ekenry, ekenlx, ekenly
    real(rt)         rerx, rery, relx, rely
    real(rt)         rrnewrx, rrnewry, rrnewlx, rrnewly
    real(rt)         runewrx, runewry, runewlx, runewly
    real(rt)         rvnewrx, rvnewry, rvnewlx, rvnewly
    real(rt)         rwnewrx, rwnewry, rwnewlx, rwnewly
    real(rt)         renewrx, renewry, renewlx, renewly
    real(rt)         pnewrx, pnewry, pnewlx, pnewly
    real(rt)         rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
    real(rt)         pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, luge, lamge, &
         ergp, errx, ernewrx, erry, ernewry, ergm, erlx, ernewlx, erly, ernewly
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

             if (i >= ilo+1) then
                rr = qxp(i,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qxp(i,j,km,nqp) - compn
                qxpo(i,j,km,nqp) = compu/rrnew
             end if

             if (j >= jlo+1) then
                rr = qyp(i,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qyp(i,j,km,nqp) - compn
                qypo(i,j,km,nqp) = compu/rrnew
             end if

             if (i <= ihi-1) then
                rr = qxm(i+1,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qxm(i+1,j,km,nqp) - compn
                qxmo(i+1,j,km,nqp) = compu/rrnew
             end if

             if (j <= jhi-1) then
                rr = qym(i,j+1,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qym(i,j+1,km,nqp) - compn
                qymo(i,j+1,km,nqp) = compu/rrnew
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add transverse flux difference in the z-direction to the x- and
          ! y-states for the fluid variables
          !-------------------------------------------------------------------

          pggp  = qz(i,j,kc,GDPRES)
          pggm  = qz(i,j,km,GDPRES)
          ugp  = qz(i,j,kc,GDW   )
          ugm  = qz(i,j,km,GDW   )
          gegp = qz(i,j,kc,GDGAME)
          gegm = qz(i,j,km,GDGAME)
#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d-1,QLAMS:QLAMS+ngroups-1)
          ergp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,k3d-1,QGAMCG)
#else
          gamc = qaux(i,j,k3d-1,QGAMC)
#endif

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmz = - cdtdz*sum(lamge)
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdz*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdz*HALF*(ugp+ugm)*(ergp(g)-ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdz*HALF*(ergp(g)+ergm(g))*(ugm-ugp)*f1
             end do
          else ! mixed frame
             der(:) = cdtdz * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------

          if (i >= ilo+1) then
             ! Convert to conservation form
             rrrx = qxp(i,j,km,QRHO)
             rurx = rrrx*qxp(i,j,km,QU)
             rvrx = rrrx*qxp(i,j,km,QV)
             rwrx = rrrx*qxp(i,j,km,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,km,QU)**2 + qxp(i,j,km,QV)**2 &
                  + qxp(i,j,km,QW)**2)
             rerx = qxp(i,j,km,QREINT) + ekenrx
#ifdef RADIATION
             errx = qxp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewrx = rrrx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewrx = rurx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewrx = rvrx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewrx = rwrx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewrx = rerx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewrx = rwnewrx + dmz
             renewrx = renewrx + dre
             ernewrx = errx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewrx < ZERO) then
                rrnewrx = rrrx
                runewrx = rurx
                rvnewrx = rvrx
                rwnewrx = rwrx
                renewrx = rerx
#ifdef RADIATION
                ernewrx = errx(:)
#endif
                reset_state = .true.
             end if

             qxpo(i,j,km,QRHO) = rrnewrx
             qxpo(i,j,km,QU) = runewrx/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QV) = rvnewrx/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QW) = rwnewrx/qxpo(i,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QREINT) = renewrx - rhoekenrx

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qxpo(i,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxpo(i,j,km,QREINT) = qxp(i,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qxpo(i,j,km,QPRES) = max(pnewrx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxpo(i,j,km,QGAME) = qxp(i,j,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxpo(i,j,km,QPRES) = qxpo(i,j,km,QREINT)*(qxpo(i,j,km,QGAME)-ONE)
                   qxpo(i,j,km,QPRES) = max(qxpo(i,j,km,QPRES), small_pres)
                endif
             else
                qxpo(i,j,km,QPRES) = qxp(i,j,km,QPRES)
                qxpo(i,j,km,QGAME) = qxp(i,j,km,QGAME)
             endif

             call reset_edge_state_thermo(qxpo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qxpo(i,j,km,qrad:qradhi) = ernewrx(:)
             qxpo(i,j,km,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,km,QPRES)
             qxpo(i,j,km,qreitot) = sum(qxpo(i,j,km,qrad:qradhi)) + qxpo(i,j,km,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          if (j >= jlo+1) then
             ! Convert to conservation form
             rrry = qyp(i,j,km,QRHO)
             rury = rrry*qyp(i,j,km,QU)
             rvry = rrry*qyp(i,j,km,QV)
             rwry = rrry*qyp(i,j,km,QW)
             ekenry = HALF*rrry*(qyp(i,j,km,QU)**2 + qyp(i,j,km,QV)**2 &
                  + qyp(i,j,km,QW)**2)
             rery = qyp(i,j,km,QREINT) + ekenry
#ifdef RADIATION
             erry = qyp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewry = rrry - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewry = rury - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewry = rvry - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewry = rwry - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewry = rery - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewry = rwnewry + dmz
             renewry = renewry + dre
             ernewry = erry(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewry < ZERO) then
                rrnewry = rrry
                runewry = rury
                rvnewry = rvry
                rwnewry = rwry
                renewry = rery
#ifdef RADIATION
                ernewry = erry(:)
#endif
                reset_state = .true.
             end if

             qypo(i,j,km,QRHO) = rrnewry
             qypo(i,j,km,QU) = runewry/qypo(i,j,km,QRHO)
             qypo(i,j,km,QV) = rvnewry/qypo(i,j,km,QRHO)
             qypo(i,j,km,QW) = rwnewry/qypo(i,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,km,QRHO)
             qypo(i,j,km,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qypo(i,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qypo(i,j,km,QREINT) = qyp(i,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qypo(i,j,km,QPRES) = max(pnewry,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qypo(i,j,km,QGAME) = qyp(i,j,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qypo(i,j,km,QPRES) = qypo(i,j,km,QREINT)*(qypo(i,j,km,QGAME)-ONE)
                   qypo(i,j,km,QPRES) = max(qypo(i,j,km,QPRES), small_pres)
                endif
             else
                qypo(i,j,km,QPRES) = qyp(i,j,km,QPRES)
                qypo(i,j,km,QGAME) = qyp(i,j,km,QGAME)
             endif

             call reset_edge_state_thermo(qypo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qypo(i,j,km,qrad:qradhi) = ernewry(:)
             qypo(i,j,km,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,km,QPRES)
             qypo(i,j,km,qreitot) = sum(qypo(i,j,km,qrad:qradhi)) + qypo(i,j,km,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i <= ihi-1) then

             ! Convert to conservation form
             rrlx = qxm(i+1,j,km,QRHO)
             rulx = rrlx*qxm(i+1,j,km,QU)
             rvlx = rrlx*qxm(i+1,j,km,QV)
             rwlx = rrlx*qxm(i+1,j,km,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,km,QU)**2 + qxm(i+1,j,km,QV)**2 &
                  + qxm(i+1,j,km,QW)**2)
             relx = qxm(i+1,j,km,QREINT) + ekenlx
#ifdef RADIATION
             erlx = qxm(i+1,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlx = rrlx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewlx = rulx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewlx = rvlx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewlx = rwlx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewlx = relx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewlx = rwnewlx + dmz
             renewlx = renewlx + dre
             ernewlx = erlx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewlx < ZERO) then
                rrnewlx = rrlx
                runewlx = rulx
                rvnewlx = rvlx
                rwnewlx = rwlx
                renewlx = relx
#ifdef RADIATION
                ernewlx = erlx(:)
#endif
                reset_state = .true.
             end if

             qxmo(i+1,j,km,QRHO) = rrnewlx
             qxmo(i+1,j,km,QU) = runewlx/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QV) = rvnewlx/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QW) = rwnewlx/qxmo(i+1,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo(i+1,j,km,QREINT) = qxm(i+1,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qxmo(i+1,j,km,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo(i+1,j,km,QGAME) = qxm(i+1,j,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo(i+1,j,km,QPRES) = qxmo(i+1,j,km,QREINT)*(qxmo(i+1,j,km,QGAME)-ONE)
                   qxmo(i+1,j,km,QPRES) = max(qxmo(i+1,j,km,QPRES), small_pres)
                end if
             else
                qxmo(i+1,j,km,QPRES) = qxm(i+1,j,km,QPRES)
                qxmo(i+1,j,km,QGAME) = qxm(i+1,j,km,QGAME)
             endif

             call reset_edge_state_thermo(qxmo, qd_lo, qd_hi, i+1, j, km)

#ifdef RADIATION
             qxmo(i+1,j,km,qrad:qradhi) = ernewlx(:)
             qxmo(i+1,j,km,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i+1,j,km,QPRES)
             qxmo(i+1,j,km,qreitot) = sum(qxmo(i+1,j,km,qrad:qradhi)) + qxmo(i+1,j,km,QREINT)
#endif

          endif


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------

          if (j <= jhi-1) then
             ! Convert to conservation form
             rrly = qym(i,j+1,km,QRHO)
             ruly = rrly*qym(i,j+1,km,QU)
             rvly = rrly*qym(i,j+1,km,QV)
             rwly = rrly*qym(i,j+1,km,QW)
             ekenly = HALF*rrly*(qym(i,j+1,km,QU)**2 + qym(i,j+1,km,QV)**2 &
                  + qym(i,j+1,km,QW)**2)
             rely = qym(i,j+1,km,QREINT) + ekenly
#ifdef RADIATION
             erly = qym(i,j+1,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewly = rrly - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewly = ruly - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewly = rvly - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewly = rwly - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewly = rely - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewly = rwnewly + dmz
             renewly = renewly + dre
             ernewly = erly(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewly < ZERO) then
                rrnewly = rrly
                runewly = ruly
                rvnewly = rvly
                rwnewly = rwly
                renewly = rely
#ifdef RADIATION
                ernewly = erly(:)
#endif
                reset_state = .true.
             endif

             qymo(i,j+1,km,QRHO) = rrnewly
             qymo(i,j+1,km,QU) = runewly/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QV) = rvnewly/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QW) = rwnewly/qymo(i,j+1,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qymo(i,j+1,km,QREINT) = qym(i,j+1,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j+1,km,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qymo(i,j+1,km,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qymo(i,j+1,km,QGAME) = qym(i,j+1,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qymo(i,j+1,km,QPRES) = qymo(i,j+1,km,QREINT)*(qymo(i,j+1,km,QGAME)-ONE)
                   qymo(i,j+1,km,QPRES) = max(qymo(i,j+1,km,QPRES), small_pres)
                endif
             else
                qymo(i,j+1,km,QPRES) = qym(i,j+1,km,QPRES)
                qymo(i,j+1,km,QGAME) = qym(i,j+1,km,QGAME)
             endif

             call reset_edge_state_thermo(qymo, qd_lo, qd_hi, i, j+1, km)

#ifdef RADIATION
             qymo(i,j+1,km,qrad:qradhi) = ernewly(:)
             qymo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j+1,km,QPRES)
             qymo(i,j+1,km,qreitot) = sum(qymo(i,j+1,km,qrad:qradhi)) + qymo(i,j+1,km,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transz


  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy(qm,qmo,qp,qpo,qd_lo,qd_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxy, &
#ifdef RADIATION
                     rfxy, &
#endif
                     fx_lo,fx_hi, &
                     fyx, &
#ifdef RADIATION
                     rfyx, &
#endif
                     fy_lo,fy_hi, &
                     qx,qx_lo,qx_hi, &
                     qy,qy_lo,qy_hi, &
                     srcQ,src_lo,src_hi, &
                     hdt,cdtdx,cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer :: src_lo(3),src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    real(rt)         rfxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    real(rt)         rfyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    real(rt)          qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)         fxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt)         fyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt)          qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt)          qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    real(rt)         hdt,cdtdx,cdtdy

    integer i, j, n, nqp, ipassive

    real(rt)         rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt)         rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt)         rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt)         rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt)         pnewr, pnewl
    real(rt)         pggxp, pggxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt)         pggyp, pggym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt)         uxav, gexav, dgex, uyav, geyav, dgey
    real(rt)         pggxpm, pggxmm, ugxpm, ugxmm, gegxpm, gegxmm, duxpm, pxavm, duxm, pxnewm, gexnewm
    real(rt)         pggypm, pggymm, ugypm, ugymm, gegypm, gegymm, duypm, pyavm, duym, pynewm, geynewm
    real(rt)         uxavm, gexavm, dgexm, uyavm, geyavm, dgeym
    real(rt)         compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt)         :: dmx, dmy, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lamc, lamm, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym, ergxpm, ergypm, ergxmm, ergymm
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             rrr = qp(i,j,kc,QRHO)
             rrl = qm(i,j,kc,QRHO)

             compr = rrr*qp(i,j,kc,nqp)
             compl = rrl*qm(i,j,kc,nqp)

             rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                          - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
             rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                          - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                            - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
             compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                            - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

             qpo(i,j,kc,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k3d  ,nqp)
             qmo(i,j,kc,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nqp)

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse xy and yx differences to the z-states for the
          ! fluid variables
          !-------------------------------------------------------------------

          pggxp  = qx(i+1,j,kc,GDPRES)
          pggxm  = qx(i  ,j,kc,GDPRES)
          ugxp  = qx(i+1,j,kc,GDU   )
          ugxm  = qx(i  ,j,kc,GDU   )
          gegxp = qx(i+1,j,kc,GDGAME)
          gegxm = qx(i  ,j,kc,GDGAME)

          pggyp  = qy(i,j+1,kc,GDPRES)
          pggym  = qy(i,j  ,kc,GDPRES)
          ugyp  = qy(i,j+1,kc,GDV   )
          ugym  = qy(i,j  ,kc,GDV   )
          gegyp = qy(i,j+1,kc,GDGAME)
          gegym = qy(i,j  ,kc,GDGAME)

#ifdef RADIATION
          lamc(:) = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
          ergxp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergxm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)
          ergyp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergym = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)
#endif

          pggxpm  = qx(i+1,j,km,GDPRES)
          pggxmm  = qx(i  ,j,km,GDPRES)
          ugxpm  = qx(i+1,j,km,GDU   )
          ugxmm  = qx(i  ,j,km,GDU   )
          gegxpm = qx(i+1,j,km,GDGAME)
          gegxmm = qx(i  ,j,km,GDGAME)

          pggypm  = qy(i,j+1,km,GDPRES)
          pggymm  = qy(i,j  ,km,GDPRES)
          ugypm  = qy(i,j+1,km,GDV   )
          ugymm  = qy(i,j  ,km,GDV   )
          gegypm = qy(i,j+1,km,GDGAME)
          gegymm = qy(i,j  ,km,GDGAME)

#ifdef RADIATION
          lamm(:) = qaux(i,j,k3d-1,QLAMS:QLAMS+ngroups-1)
          ergxpm = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergxmm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)
          ergypm = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergymm = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)
#endif

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = HALF*(pggxp+pggxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
#ifdef RADIATION
          pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k3d,QGAMCG) - ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k3d,QGAMCG))*dux - uxav*dgex )
#else
          pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k3d,QGAMC) - ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k3d,QGAMC))*dux - uxav*dgex )
#endif


          duxpm = pggxpm*ugxpm - pggxmm*ugxmm
          pxavm = HALF*(pggxpm+pggxmm)
          uxavm = HALF*(ugxpm+ugxmm)
          gexavm = HALF*(gegxpm+gegxmm)
          duxm = ugxpm-ugxmm
          dgexm = gegxpm-gegxmm
#ifdef RADIATION
          pxnewm = cdtdx*(duxpm + pxavm*duxm*(qaux(i,j,k3d-1,QGAMCG) - ONE))
          gexnewm = cdtdx*( (gexavm-ONE)*(gexavm - qaux(i,j,k3d-1,QGAMCG))*duxm - uxavm*dgexm )
#else
          pxnewm = cdtdx*(duxpm + pxavm*duxm*(qaux(i,j,k3d-1,QGAMC) - ONE))
          gexnewm = cdtdx*( (gexavm-ONE)*(gexavm - qaux(i,j,k3d-1,QGAMC))*duxm - uxavm*dgexm )
#endif
          duyp = pggyp*ugyp - pggym*ugym
          pyav = HALF*(pggyp+pggym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
#ifdef RADIATION
          pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k3d,QGAMCG) - ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k3d,QGAMCG))*duy - uyav*dgey )
#else
          pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k3d,QGAMC) - ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k3d,QGAMC))*duy - uyav*dgey )
#endif
          duypm = pggypm*ugypm - pggymm*ugymm
          pyavm = HALF*(pggypm+pggymm)
          uyavm = HALF*(ugypm+ugymm)
          geyavm = HALF*(gegypm+gegymm)
          duym = ugypm-ugymm
          dgeym = gegypm-gegymm
#ifdef RADIATION
          pynewm = cdtdy*(duypm + pyavm*duym*(qaux(i,j,k3d-1,QGAMCG) - ONE))
          geynewm = cdtdy*( (geyavm-ONE)*(geyavm - qaux(i,j,k3d-1,QGAMCG))*duym - uyavm*dgeym )
#else
          pynewm = cdtdy*(duypm + pyavm*duym*(qaux(i,j,k3d-1,QGAMC) - ONE))
          geynewm = cdtdy*( (geyavm-ONE)*(geyavm - qaux(i,j,k3d-1,QGAMC))*duym - uyavm*dgeym )
#endif

          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------

          ! Convert to conservation form
          rrr = qp(i,j,kc,QRHO)
          rur = rrr*qp(i,j,kc,QU)
          rvr = rrr*qp(i,j,kc,QV)
          rwr = rrr*qp(i,j,kc,QW)
          ekenr = HALF*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
               qp(i,j,kc,QW)**2)
          rer = qp(i,j,kc,QREINT) + ekenr
#ifdef RADIATION
          err = qp(i,j,kc,qrad:qradhi)

          lgex = lamc(:) * (ergxp(:)-ergxm(:))
          lgey = lamc(:) * (ergyp(:)-ergym(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          lugex = HALF*(ugxp+ugxm) * lgex(:)
          lugey = HALF*(ugyp+ugym) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if
#endif

          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                       - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
          runewr = rur - cdtdx*(fxy(i+1,j,kc,UMX) - fxy(i,j,kc,UMX)) &
                       - cdtdy*(fyx(i,j+1,kc,UMX) - fyx(i,j,kc,UMX))
          rvnewr = rvr - cdtdx*(fxy(i+1,j,kc,UMY) - fxy(i,j,kc,UMY)) &
                       - cdtdy*(fyx(i,j+1,kc,UMY) - fyx(i,j,kc,UMY))
          rwnewr = rwr - cdtdx*(fxy(i+1,j,kc,UMZ) - fxy(i,j,kc,UMZ)) &
                       - cdtdy*(fyx(i,j+1,kc,UMZ) - fyx(i,j,kc,UMZ))
          renewr = rer - cdtdx*(fxy(i+1,j,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
                       - cdtdy*(fyx(i,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))
#ifdef RADIATION
          runewr = runewr + dmx
          rvnewr = rvnewr + dmy
          renewr = renewr + dre
          ernewr = err(:) - cdtdx*(rfxy(i+1,j,kc,:) - rfxy(i,j,kc,:)) &
               &          - cdtdy*(rfyx(i,j+1,kc,:) - rfyx(i,j,kc,:))  &
               &          + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
             rrnewr = rrr
             runewr = rur
             rvnewr = rvr
             rwnewr = rwr
             renewr = rer
#ifdef RADIATION
             ernewr = err(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qpo(i,j,kc,QRHO  ) = rrnewr        
          qpo(i,j,kc,QU    ) = runewr/rrnewr
          qpo(i,j,kc,QV    ) = rvnewr/rrnewr
          qpo(i,j,kc,QW    ) = rwnewr/rrnewr

          ! for ppm_type > 0 we already added the piecewise parabolic traced
          ! source terms to the normal edge states.
          if (ppm_type == 0) then
             qpo(i,j,kc,QRHO  ) = qpo(i,j,kc,QRHO  ) + hdt*srcQ(i,j,k3d,QRHO)
             qpo(i,j,kc,QU:QW) = qpo(i,j,kc,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
          endif

          ! note: we run the risk of (rho e) being negative here
          rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
          qpo(i,j,kc,QREINT) = renewr - rhoekenr
          if (ppm_type == 0) then
             qpo(i,j,kc,QREINT) = qpo(i,j,kc,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
          endif

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qpo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qpo(i,j,kc,QREINT) = qp(i,j,kc,QREINT) &
                     - cdtdx*(fxy(i+1,j,kc,UEINT) - fxy(i,j,kc,UEINT) + pxav*dux) &
                     - cdtdy*(fyx(i,j+1,kc,UEINT) - fyx(i,j,kc,UEINT) + pyav*duy) 
             endif

             if (ppm_type == 0) then
                qpo(i,j,kc,QREINT) = qpo(i,j,kc,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
                qpo(i,j,kc,QPRES) = pnewr
                if (ppm_type == 0) then
                   qpo(i,j,kc,QPRES) = qpo(i,j,kc,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                endif
             else
                ! Update gammae with its transverse terms
                qpo(i,j,kc,QGAME) = qp(i,j,kc,QGAME) + gexnew + geynew

                ! and compute the p edge state from this and (rho e)
                qpo(i,j,kc,QPRES) = qpo(i,j,kc,QREINT)*(qpo(i,j,kc,QGAME)-ONE)
             endif
          else
             qpo(i,j,kc,QPRES) = qp(i,j,kc,QPRES)
             if (ppm_type == 0) then
                qpo(i,j,kc,QPRES) = qpo(i,j,kc,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
             endif
             qpo(i,j,kc,QGAME) = qp(i,j,kc,QGAME)
          endif

          qpo(i,j,kc,QPRES) = max(qpo(i,j,kc,QPRES), small_pres)

          call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qpo(i,j,kc,qptot  ) = sum(lamc(:)*ernewr(:)) + qpo(i,j,kc,QPRES)
          qpo(i,j,kc,qreitot) = sum(qpo(i,j,kc,qrad:qradhi)) + qpo(i,j,kc,QREINT)
#endif

          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------

          ! Convert to conservation form
          rrl = qm(i,j,kc,QRHO)
          rul = rrl*qm(i,j,kc,QU)
          rvl = rrl*qm(i,j,kc,QV)
          rwl = rrl*qm(i,j,kc,QW)
          ekenl = HALF*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
               qm(i,j,kc,QW)**2)
          rel = qm(i,j,kc,QREINT) + ekenl
#ifdef RADIATION
          erl = qm(i,j,kc,qrad:qradhi)

          lgex = lamm(:) * (ergxpm(:)-ergxmm(:))
          lgey = lamm(:) * (ergypm(:)-ergymm(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          lugex = HALF*(ugxpm+ugxmm) * lgex(:)
          lugey = HALF*(ugypm+ugymm) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ugxpm+ugxmm)*(ergxpm(g)-ergxmm(g)) &
                     +       cdtdy*HALF*(ugypm+ugymm)*(ergypm(g)-ergymm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ergxpm(g)+ergxmm(g))*(ugxmm-ugxpm) &
                     +       cdtdy*HALF*(ergypm(g)+ergymm(g))*(ugymm-ugypm) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if
#endif

          ! Add transverse predictor
          rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                       - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxy(i+1,j,km,UMX) - fxy(i,j,km,UMX)) &
                       - cdtdy*(fyx(i,j+1,km,UMX) - fyx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxy(i+1,j,km,UMY) - fxy(i,j,km,UMY)) &
                       - cdtdy*(fyx(i,j+1,km,UMY) - fyx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxy(i+1,j,km,UMZ) - fxy(i,j,km,UMZ)) &
                       - cdtdy*(fyx(i,j+1,km,UMZ) - fyx(i,j,km,UMZ))
          renewl = rel - cdtdx*(fxy(i+1,j,km,UEDEN) - fxy(i,j,km,UEDEN)) &
                       - cdtdy*(fyx(i,j+1,km,UEDEN) - fyx(i,j,km,UEDEN))
#ifdef RADIATION
          runewl = runewl + dmx
          rvnewl = rvnewl + dmy
          renewl = renewl + dre
          ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,km,:) - rfxy(i,j,km,:)) &
               &          - cdtdy*(rfyx(i  ,j+1,km,:) - rfyx(i,j,km,:)) &
               &          + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
             rrnewl = rrl
             runewl = rul
             rvnewl = rvl
             rwnewl = rwl
             renewl = rel
#ifdef RADIATION
             ernewl = erl(:)
#endif
             reset_state = .true.
          endif

          qmo(i,j,kc,QRHO  ) = rrnewl       
          qmo(i,j,kc,QU    ) = runewl/rrnewl
          qmo(i,j,kc,QV    ) = rvnewl/rrnewl
          qmo(i,j,kc,QW    ) = rwnewl/rrnewl

          ! for ppm_type > 0 we already added the piecewise parabolic traced
          ! source terms to the normal edge states.
          if (ppm_type == 0) then
             qmo(i,j,kc,QRHO  ) = qmo(i,j,kc,QRHO  ) + hdt*srcQ(i,j,k3d-1,QRHO)
             qmo(i,j,kc,QU:QW) = qmo(i,j,kc,QU:QW) + hdt * srcQ(i,j,k3d-1,QU:QW)
          endif

          ! note: we run the risk of (rho e) being negative here
          rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
          qmo(i,j,kc,QREINT) = renewl - rhoekenl
          if (ppm_type == 0) then
             qmo(i,j,kc,QREINT) = qmo(i,j,kc,QREINT) + hdt*srcQ(i,j,k3d-1,QREINT)
          endif

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qmo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qmo(i,j,kc,QREINT) = qm(i,j,kc,QREINT) &
                     - cdtdx*(fxy(i+1,j,km,UEINT) - fxy(i,j,km,UEINT) + pxavm*duxm) &
                     - cdtdy*(fyx(i,j+1,km,UEINT) - fyx(i,j,km,UEINT) + pyavm*duym)
                if (ppm_type == 0) then
                   qmo(i,j,kc,QREINT) = qmo(i,j,kc,QREINT) + hdt*srcQ(i,j,k3d-1,QREINT)
                endif
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm
                qmo(i,j,kc,QPRES) = pnewl
                if (ppm_type == 0) then
                   qmo(i,j,kc,QPRES) = qmo(i,j,kc,QPRES) + hdt*srcQ(i,j,k3d-1,QPRES)
                endif
             else
                ! Update gammae with its transverse terms
                qmo(i,j,kc,QGAME) = qm(i,j,kc,QGAME) + gexnewm + geynewm

                ! and compute the p edge state from this and (rho e)
                qmo(i,j,kc,QPRES) = qmo(i,j,kc,QREINT)*(qmo(i,j,kc,QGAME)-ONE)
             endif
          else
             qmo(i,j,kc,QPRES) = qm(i,j,kc,QPRES)
             if (ppm_type == 0) then
                qmo(i,j,kc,QPRES) = qmo(i,j,kc,QPRES) + hdt*srcQ(i,j,k3d-1,QPRES)
             endif
             qmo(i,j,kc,QGAME) = qm(i,j,kc,QGAME)
          endif

          qmo(i,j,kc,QPRES) = max(qmo(i,j,kc,QPRES), small_pres)

          call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qmo(i,j,kc,qptot  ) = sum(lamm(:)*ernewl(:)) + qmo(i,j,kc,QPRES)
          qmo(i,j,kc,qreitot) = sum(qmo(i,j,kc,qrad:qradhi)) + qmo(i,j,kc,QREINT)
#endif

       enddo
    enddo

  end subroutine transxy


  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz(qm, qmo, qp, qpo, qd_lo, qd_hi, &
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
                     qx, qx_lo, qx_hi, &
                     qz, qz_lo, qz_hi, &
                     srcQ, src_lo, src_hi, &
                     hdt, cdtdx, cdtdz, ilo, ihi, jlo, jhi, km, kc, k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: fz_lo(3),fz_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer :: qz_lo(3),qz_hi(3)
    integer :: src_lo(3),src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    real(rt)         rfxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    real(rt)         rfzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    real(rt)          qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)          qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)         fxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt)         fzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt)          qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt)          qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    real(rt)         hdt,cdtdx,cdtdz

    integer i, j, n, nqp, ipassive

    real(rt)         rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt)         rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt)         rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt)         rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt)         pnewr, pnewl
    real(rt)         pggxp, pggxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt)         pggzp, pggzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt)         uxav, gexav, dgex, uzav, gezav, dgez
    real(rt)         compr, compl, compnr, compnl, drr, dcompn

#ifdef RADIATION
    real(rt)         :: dmx, dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             drr    = - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                      - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             dcompn = - cdtdx*(fxz(i+1,j,km,n   ) - fxz(i,j,km,n)) &
                      - cdtdz*(fzx(i  ,j,kc,n   ) - fzx(i,j,km,n))

             if (j >= jlo+1) then
                rrr = qp(i,j,km,QRHO)
                compr = rrr*qp(i,j,km,nqp)

                rrnewr = rrr + drr
                compnr = compr + dcompn

                qpo(i,j  ,km,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nqp)
             end if

             if (j <= jhi-1) then
                rrl = qm(i,j+1,km,QRHO)
                compl = rrl*qm(i,j+1,km,nqp)

                rrnewl = rrl + drr
                compnl = compl + dcompn

                qmo(i,j+1,km,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nqp)
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse xz and zx differences to the y-states for the
          ! fluid variables
          !-------------------------------------------------------------------

#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
#endif

          pggxp  = qx(i+1,j,km,GDPRES)
          pggxm  = qx(i  ,j,km,GDPRES)
          ugxp  = qx(i+1,j,km,GDU   )
          ugxm  = qx(i  ,j,km,GDU   )
          gegxp = qx(i+1,j,km,GDGAME)
          gegxm = qx(i  ,j,km,GDGAME)
#ifdef RADIATION
          ergxp = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergxm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          pggzp  = qz(i,j,kc,GDPRES)
          pggzm  = qz(i,j,km,GDPRES)
          ugzp  = qz(i,j,kc,GDW   )
          ugzm  = qz(i,j,km,GDW   )
          gegzp = qz(i,j,kc,GDGAME)
          gegzm = qz(i,j,km,GDGAME)
#ifdef RADIATION
          ergzp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergzm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = HALF*(pggxp+pggxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
#ifdef RADIATION
          pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k3d,QGAMCG) - ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k3d,QGAMCG))*dux - uxav*dgex )
#else
          pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k3d,QGAMC) - ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k3d,QGAMC))*dux - uxav*dgex )
#endif

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = HALF*(pggzp+pggzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
#ifdef RADIATION
          pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k3d,QGAMCG) - ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k3d,QGAMCG))*duz - uzav*dgez )
#else
          pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k3d,QGAMC) - ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k3d,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
          lgex = lambda(:) * (ergxp(:)-ergxm(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmx = - cdtdx*sum(lgex)
          dmz = - cdtdz*sum(lgez)
          lugex = HALF*(ugxp+ugxm) * lgex(:)
          lugez = HALF*(ugzp+ugzm) * lgez(:)
          dre = -cdtdx * sum(lugex) - cdtdz * sum(lugez)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdx*lugex + cdtdz*lugez
          end if
#endif

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          if (j >= jlo+1) then
             ! Convert to conservation form
             rrr = qp(i,j,km,QRHO)
             rur = rrr*qp(i,j,km,QU)
             rvr = rrr*qp(i,j,km,QV)
             rwr = rrr*qp(i,j,km,QW)
             ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
             rer = qp(i,j,km,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                          - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
             runewr = rur - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                          - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
             rvnewr = rvr - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                          - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
             rwnewr = rwr - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                          - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
             renewr = rer - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                          - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))
#ifdef RADIATION
             runewr = runewr + dmx
             rwnewr = rwnewr + dmz
             renewr = renewr + dre
             ernewr = err(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
                             - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             qpo(i,j,km,QRHO  ) = rrnewr       
             qpo(i,j,km,QU    ) = runewr/rrnewr
             qpo(i,j,km,QV    ) = rvnewr/rrnewr
             qpo(i,j,km,QW    ) = rwnewr/rrnewr

             ! for ppm_type > 0 we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_type == 0) then
                qpo(i,j,km,QRHO  ) = qpo(i,j,km,QRHO  ) + hdt*srcQ(i,j,k3d,QRHO)
                qpo(i,j,km,QU:QW) = qpo(i,j,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,km,QREINT) = renewr - rhoekenr
             if (ppm_type == 0) then
                qpo(i,j,km,QREINT) = qpo(i,j,km,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,km,QREINT) <= ZERO) then
                   qpo(i,j,km,QREINT) = qp(i,j,km,QREINT) &
                        - cdtdx*(fxz(i+1,j,km,UEINT) - fxz(i,j,km,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i  ,j,kc,UEINT) - fzx(i,j,km,UEINT) + pzav*duz)
                   if (ppm_type == 0) then
                      qpo(i,j,km,QREINT) = qpo(i,j,km,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
                   endif
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,km,QPRES) - pxnew - pznew
                   qpo(i,j,km,QPRES) = pnewr
                   if (ppm_type == 0) then
                      qpo(i,j,km,QPRES) = qpo(i,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                   endif
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,km,QGAME) = qp(i,j,km,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT)*(qpo(i,j,km,QGAME)-ONE)
                endif
             else
                qpo(i,j,km,QPRES) = qp(i,j,km,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,km,QPRES) = qpo(i,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                endif
                qpo(i,j,km,QGAME) = qp(i,j,km,QGAME)
             endif

             qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qpo(i,j,km,qrad:qradhi) = ernewr(:)
             qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
             qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)
#endif
          end if


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------

          if (j <= jhi-1) then
             ! Convert to conservation form
             rrl = qm(i,j+1,km,QRHO)
             rul = rrl*qm(i,j+1,km,QU)
             rvl = rrl*qm(i,j+1,km,QV)
             rwl = rrl*qm(i,j+1,km,QW)
             ekenl = HALF*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
             rel = qm(i,j+1,km,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i,j+1,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                          - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
             runewl = rul - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                          - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
             rvnewl = rvl - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                          - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
             rwnewl = rwl - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                          - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
             renewl = rel - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                          - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
                             - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
#ifdef RADIATION
                ernewl = erl(:)
#endif
                reset_state = .true.
             endif

             qmo(i,j+1,km,QRHO  ) = rrnewl        
             qmo(i,j+1,km,QU    ) = runewl/rrnewl
             qmo(i,j+1,km,QV    ) = rvnewl/rrnewl
             qmo(i,j+1,km,QW    ) = rwnewl/rrnewl

             ! for ppm_type > 0 we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_type == 0) then
                qmo(i,j+1,km,QRHO  ) = qmo(i,j+1,km,QRHO  ) + hdt*srcQ(i,j,k3d,QRHO)
                qmo(i,j+1,km,QU:QW) = qmo(i,j+1,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j+1,km,QREINT) = renewl - rhoekenl
             if (ppm_type == 0) then
                qmo(i,j+1,km,QREINT) = qmo(i,j+1,km,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j+1,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j+1,km,QREINT) = qm(i,j+1,km,QREINT) &
                        - cdtdx*(fxz(i+1,j,km,UEINT) - fxz(i,j,km,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i,j,kc,UEINT) - fzx(i,j,km,UEINT) + pzav*duz)
                   if (ppm_type == 0) then
                      qmo(i,j+1,km,QREINT) = qmo(i,j+1,km,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
                   endif
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew
                   qmo(i,j+1,km,QPRES) = pnewl
                   if (ppm_type == 0) then
                      qmo(i,j+1,km,QPRES) = qmo(i,j+1,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                   endif
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j+1,km,QGAME) = qm(i,j+1,km,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j+1,km,QPRES) = qmo(i,j+1,km,QREINT)*(qmo(i,j+1,km,QGAME)-ONE)
                endif
             else
                qmo(i,j+1,km,QPRES) = qm(i,j+1,km,QPRES)
                if (ppm_type == 0) then
                   qmo(i,j+1,km,QPRES) = qmo(i,j+1,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                endif
                qmo(i,j+1,km,QGAME) = qm(i,j+1,km,QGAME)
             endif

             qmo(i,j+1,km,QPRES) = max(qmo(i,j+1,km,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i, j+1, km)

#ifdef RADIATION
             qmo(i,j+1,km,qrad:qradhi) = ernewl(:)
             qmo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j+1,km,QPRES)
             qmo(i,j+1,km,qreitot) = sum(qmo(i,j+1,km,qrad:qradhi)) + qmo(i,j+1,km,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transxz


  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz(qm, qmo, qp, qpo, qd_lo, qd_hi, &
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
                     qy, qy_lo, qy_hi, &
                     qz, qz_lo, qz_hi, &
                     srcQ, src_lo, src_hi, &
                     hdt, cdtdy, cdtdz, ilo, ihi, jlo, jhi, km, kc, k3d)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer :: qd_lo(3),qd_hi(3)
    integer :: qa_lo(3),qa_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: fz_lo(3),fz_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer :: qz_lo(3),qz_hi(3)
    integer :: src_lo(3),src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    real(rt)         rfyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    real(rt)         rfzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    real(rt)         qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)         fyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt)         fzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt)          qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt)          qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    real(rt)         hdt,cdtdy,cdtdz

    integer i, j, n, nqp, ipassive

    real(rt)         rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt)         rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt)         rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt)         rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt)         pnewr, pnewl
    real(rt)         pggyp, pggym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt)         pggzp, pggzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt)         uyav, geyav, dgey, uzav, gezav, dgez
    real(rt)         compr, compl, compnr, compnl
    real(rt)         drr, dcompn

#ifdef RADIATION
    real(rt)         :: dmy, dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             drr    = - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                      - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             dcompn = - cdtdy*(fyz(i,j+1,km,n   ) - fyz(i,j,km,n)) &
                      - cdtdz*(fzy(i,j  ,kc,n   ) - fzy(i,j,km,n))

             if (i >= ilo+1) then
                rrr = qp(i,j,km,QRHO)
                compr = rrr*qp(i,j,km,nqp)

                rrnewr = rrr +drr
                compnr = compr +dcompn

                qpo(i  ,j,km,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nqp)
             end if

             if (i <= ihi-1) then
                rrl = qm(i+1,j,km,QRHO)
                compl = rrl*qm(i+1,j,km,nqp)

                rrnewl = rrl + drr
                compnl = compl +dcompn

                qmo(i+1,j,km,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nqp)
             end if
          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse yz and zy differences to the x-states for the
          ! fluid variables
          !-------------------------------------------------------------------

#ifdef RADIATION
          lambda(:) = qaux(i,j,k3d,QLAMS:QLAMS+ngroups-1)
#endif

          pggyp  = qy(i,j+1,km,GDPRES)
          pggym  = qy(i,j  ,km,GDPRES)
          ugyp  = qy(i,j+1,km,GDV   )
          ugym  = qy(i,j  ,km,GDV   )
          gegyp = qy(i,j+1,km,GDGAME)
          gegym = qy(i,j  ,km,GDGAME)
#ifdef RADIATION
          ergyp = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergym = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)
#endif

          pggzp  = qz(i,j,kc,GDPRES)
          pggzm  = qz(i,j,km,GDPRES)
          ugzp  = qz(i,j,kc,GDW   )
          ugzm  = qz(i,j,km,GDW   )
          gegzp = qz(i,j,kc,GDGAME)
          gegzm = qz(i,j,km,GDGAME)
#ifdef RADIATION
          ergzp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergzm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          duyp = pggyp*ugyp - pggym*ugym
          pyav = HALF*(pggyp+pggym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
#ifdef RADIATION
          pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k3d,QGAMCG) - ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k3d,QGAMCG))*duy - uyav*dgey )
#else
          pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k3d,QGAMC) - ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k3d,QGAMC))*duy - uyav*dgey )
#endif

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = HALF*(pggzp+pggzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
#ifdef RADIATION
          pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k3d,QGAMCG) - ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k3d,QGAMCG))*duz - uzav*dgez )
#else
          pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k3d,QGAMC) - ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k3d,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
          lgey = lambda(:) * (ergyp(:)-ergym(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmy = - cdtdy*sum(lgey)
          dmz = - cdtdz*sum(lgez)
          lugey = HALF*(ugyp+ugym) * lgey(:)
          lugez = HALF*(ugzp+ugzm) * lgez(:)
          dre = -cdtdy*sum(lugey) - cdtdz*sum(lugez)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) &
                     +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) &
                     +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdy*lugey + cdtdz*lugez
          end if
#endif

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------

          if (i >= ilo+1) then
             ! Convert to conservation form
             rrr = qp(i,j,km,QRHO)
             rur = rrr*qp(i,j,km,QU)
             rvr = rrr*qp(i,j,km,QV)
             rwr = rrr*qp(i,j,km,QW)
             ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
                  qp(i,j,km,QW)**2)
             rer = qp(i,j,km,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                          - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
             runewr = rur - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                          - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
             rvnewr = rvr - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                          - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
             rwnewr = rwr - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                          - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
             renewr = rer - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                          - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))
#ifdef RADIATION
             rvnewr = rvnewr + dmy
             rwnewr = rwnewr + dmz
             renewr = renewr + dre
             ernewr = err(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
                             - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             end if

             qpo(i,j,km,QRHO  ) = rrnewr        
             qpo(i,j,km,QU    ) = runewr/rrnewr
             qpo(i,j,km,QV    ) = rvnewr/rrnewr
             qpo(i,j,km,QW    ) = rwnewr/rrnewr

             ! for ppm_type > 0 we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_type == 0) then
                qpo(i,j,km,QRHO  ) = qpo(i,j,km,QRHO  ) + hdt*srcQ(i,j,k3d,QRHO)
                qpo(i,j,km,QU:QW) = qpo(i,j,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,km,QREINT) = renewr - rhoekenr
             if (ppm_type == 0) then
                qpo(i,j,km,QREINT) = qpo(i,j,km,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qpo(i,j,km,QREINT) = qp(i,j,km,QREINT) &
                        - cdtdy*(fyz(i,j+1,km,UEINT) - fyz(i,j,km,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,kc,UEINT) - fzy(i,j,km,UEINT) + pzav*duz)
                   if (ppm_type == 0) then
                      qpo(i,j,km,QREINT) = qpo(i,j,km,QREINT) + hdt*srcQ(i,j,k3d,QREINT)
                   endif
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,km,QPRES) - pynew - pznew
                   qpo(i,j,km,QPRES) = pnewr
                   if (ppm_type == 0) then
                      qpo(i,j,km,QPRES) = qpo(i,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                   endif
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,km,QGAME) = qp(i,j,km,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT)*(qpo(i,j,km,QGAME)-ONE)
                end if
             else
                qpo(i,j,km,QPRES) = qp(i,j,km,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,km,QPRES) = qpo(i,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                endif
                qpo(i,j,km,QGAME) = qp(i,j,km,QGAME)
             endif

             qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qpo(i,j,km,qrad:qradhi) = ernewr(:)
             qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
             qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)
#endif

          endif

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i <=ihi-1) then
             ! Convert to conservation form
             rrl = qm(i+1,j,km,QRHO)
             rul = rrl*qm(i+1,j,km,QU)
             rvl = rrl*qm(i+1,j,km,QV)
             rwl = rrl*qm(i+1,j,km,QW)
             ekenl = HALF*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
                  qm(i+1,j,km,QW)**2)
             rel = qm(i+1,j,km,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i+1,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                          - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
             runewl = rul - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                          - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
             rvnewl = rvl - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                          - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
             rwnewl = rwl - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                          - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
             renewl = rel - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                          - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))
#ifdef RADIATION
             rvnewl = rvnewl + dmy
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
                  &          - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
                  &          + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
#ifdef RADIATION
                ernewl = erl(:)
#endif
                reset_state = .true.
             endif

             qmo(i+1,j,km,QRHO   ) = rrnewl       
             qmo(i+1,j,km,QU     ) = runewl/rrnewl
             qmo(i+1,j,km,QV     ) = rvnewl/rrnewl
             qmo(i+1,j,km,QW     ) = rwnewl/rrnewl

             ! for ppm_type > 0 we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_type == 0) then
                qmo(i+1,j,km,QRHO   ) = qmo(i+1,j,km,QRHO   ) + hdt*srcQ(i,j,k3d,QRHO)
                qmo(i+1,j,km,QU:QW) = qmo(i+1,j,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i+1,j,km,QREINT ) = renewl - rhoekenl
             if (ppm_type == 0) then
                qmo(i+1,j,km,QREINT ) = qmo(i+1,j,km,QREINT ) + hdt*srcQ(i,j,k3d,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i+1,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i+1,j,km,QREINT ) = qm(i+1,j,km,QREINT) &
                        - cdtdy*(fyz(i,j+1,km,UEINT) - fyz(i,j,km,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,kc,UEINT) - fzy(i,j,km,UEINT) + pzav*duz)
                   if (ppm_type == 0) then
                      qmo(i+1,j,km,QREINT ) = qmo(i+1,j,km,QREINT ) + hdt*srcQ(i,j,k3d,QREINT)
                   endif
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i+1,j,km,QPRES) - pynew - pznew
                   qmo(i+1,j,km,QPRES  ) = pnewl
                   if (ppm_type == 0) then
                      qmo(i+1,j,km,QPRES  ) = qmo(i+1,j,km,QPRES  ) + hdt*srcQ(i,j,k3d,QPRES)
                   endif
                else
                   ! Update gammae with its transverse terms
                   qmo(i+1,j,km,QGAME) = qm(i+1,j,km,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i+1,j,km,QPRES) = qmo(i+1,j,km,QREINT)*(qmo(i+1,j,km,QGAME)-ONE)
                end if
             else
                qmo(i+1,j,km,QPRES  ) = qm(i+1,j,km,QPRES)
                if (ppm_type == 0) then 
                   qmo(i+1,j,km,QPRES  ) = qmo(i+1,j,km,QPRES  ) + hdt*srcQ(i,j,k3d,QPRES)
                endif
                qmo(i+1,j,km,QGAME) = qm(i+1,j,km,QGAME)
             endif

             qmo(i+1,j,km,QPRES) = max(qmo(i+1,j,km,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i+1, j, km)

#ifdef RADIATION
             qmo(i+1,j,km,qrad:qradhi) = ernewl(:)
             qmo(i+1,j,km,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i+1,j,km,QPRES)
             qmo(i+1,j,km,qreitot) = sum(qmo(i+1,j,km,qrad:qradhi)) + qmo(i+1,j,km,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transyz

end module ctu_advection_module
