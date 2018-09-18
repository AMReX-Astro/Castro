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
    use transverse_module, only : transx1, transx2, transy1, transy2, transz, &
                                  transxy, transxz, transyz
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

end module ctu_advection_module
