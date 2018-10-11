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
    use prob_params_module, only : dg

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

    integer :: n
    integer :: i, j, k, iwave, idim

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

    qt_lo = lo(:) - dg(:)
    qt_hi = hi(:) + 2*dg(:)

    It_lo = lo(:) - dg(:)
    It_hi = hi(:) + dg(:)

    shk_lo(:) = lo(:) - dg(:)
    shk_hi(:) = hi(:) + dg(:)

    fx_lo = [lo(1)    , lo(2) - 1, lo(3) - 1]
    fx_hi = [hi(1) + 1, hi(2) + 1, hi(3) + 1]

    fy_lo = [lo(1) - 1, lo(2)    , lo(3) - 1]
    fy_hi = [hi(1) + 1, hi(2) + 1, hi(3) + 1]

    fz_lo = [lo(1) - 1, lo(2) - 1, lo(3)    ]
    fz_hi = [hi(1) + 1, hi(2) + 1, hi(3) + 1]

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

    call shock(q, qd_lo, qd_hi, &
               shk, shk_lo, shk_hi, &
               lo, hi, dx)

    ! Store the shock data for future use in the burning step.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,k,USHK) = shk(i,j,k)
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
       call shock(q, qd_lo, qd_hi, &
                  shk, shk_lo, shk_hi, &
                  lo, hi, dx)
    else
       shk(:,:,:) = ZERO
    endif
#endif


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

    if (ppm_type > 0) then

       do n = 1, NQ
          call ppm_reconstruct(q, qd_lo, qd_hi, NQ, n, &
                               flatn, qd_lo, qd_hi, &
                               sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                               lo, hi, dx)

          call ppm_int_profile(q, qd_lo, qd_hi, NQ, n, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                               Ip, Im, It_lo, It_hi, NQ, n, &
                               lo, hi, dx, dt)
       end do

       ! source terms
       do n = 1, QVAR
          if (source_nonzero(n)) then
             call ppm_reconstruct(srcQ, src_lo, src_hi, QVAR, n, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  lo, hi, dx)

             call ppm_int_profile(srcQ, src_lo, src_hi, QVAR, n, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                                  Ip_src, Im_src, It_lo, It_hi, QVAR, n, &
                                  lo, hi, dx, dt)
          else
             Ip_src(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),:,:,n) = ZERO
             Im_src(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),:,:,n) = ZERO
          endif

       enddo

       ! this probably doesn't support radiation
       if (ppm_temp_fix /= 1) then
          call ppm_reconstruct(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               flatn, qd_lo, qd_hi, &
                               sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                               lo, hi, dx)

          call ppm_int_profile(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                               Ip_gc, Im_gc, It_lo, It_hi, 1, 1, &
                               lo, hi, dx, dt)
       else

          do iwave = 1, 3
             do idim = 1, 3

                do k = lo(3)-1, hi(3)+1
                   do j = lo(2)-1, hi(2)+1
                      do i = lo(1)-1, hi(1)+1

                         eos_state % rho = Ip(i,j,k,idim,iwave,QRHO)
                         eos_state % T   = Ip(i,j,k,idim,iwave,QTEMP)

                         eos_state % xn  = Ip(i,j,k,idim,iwave,QFS:QFS+nspec-1)
                         eos_state % aux = Ip(i,j,k,idim,iwave,QFX:QFX+naux-1)

                         call eos(eos_input_rt, eos_state)

                         Ip(i,j,k,idim,iwave,QPRES)  = eos_state % p
                         Ip(i,j,k,idim,iwave,QREINT) = eos_state % e * Ip(i,j,k,idim,iwave,QRHO)
                         Ip_gc(i,j,k,idim,iwave,1)   = eos_state % gam1
                      end do
                   end do
                end do

                do k = lo(3)-1, hi(3)+1
                   do j = lo(2)-1, hi(2)+1
                      do i = lo(1)-1, hi(1)+1
                         eos_state % rho = Im(i,j,k,idim,iwave,QRHO)
                         eos_state % T   = Im(i,j,k,idim,iwave,QTEMP)

                         eos_state % xn  = Im(i,j,k,idim,iwave,QFS:QFS+nspec-1)
                         eos_state % aux = Im(i,j,k,idim,iwave,QFX:QFX+naux-1)

                         call eos(eos_input_rt, eos_state)

                         Im(i,j,k,idim,iwave,QPRES)  = eos_state % p
                         Im(i,j,k,idim,iwave,QREINT) = eos_state % e * Im(i,j,k,idim,iwave,QRHO)
                         Im_gc(i,j,k,idim,iwave,1)   = eos_state % gam1
                      end do
                   end do
                end do

             end do
          end do

       end if

       ! Compute U_x and U_y at kc (k3d)

#ifdef RADIATION
       call tracexy_ppm_rad(q, qd_lo, qd_hi, &
                            qaux, qa_lo, qa_hi, &
                            Ip, Im, Ip_src, Im_src, It_lo, It_hi, &
                            qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                            lo, hi, domlo, domhi, &
                            dx, dt)

       call tracez_ppm_rad(q, qd_lo, qd_hi, &
                           qaux, qa_lo, qa_hi, &
                           Ip, Im, Ip_src, Im_src, It_lo, It_hi, &
                           qzm, qzp, qt_lo, qt_hi, &
                           lo, hi, domlo, domhi, &
                           dt)

#else
       call tracexy_ppm(q, qd_lo, qd_hi, &
                        qaux, qa_lo, qa_hi, &
                        Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, It_lo, It_hi, &
                        qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                        lo, hi, domlo, domhi, &
                        dx, dt)

       call tracez_ppm(q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, It_lo, It_hi, &
                       qzm, qzp, qt_lo, qt_hi, &
                       lo, hi, domlo, domhi, &
                       dt)
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
                   lo, hi)

       if (use_pslope .eq. 1) &
            call pslope(q, flatn, qd_lo, qd_hi, &
                        dqx, dqy, dqz, qt_lo, qt_hi, &
                        srcQ, src_lo, src_hi, &
                        lo, hi, dx)

       ! Compute U_x and U_y at kc (k3d)
       call tracexy(q, qd_lo, qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    dqx, dqy, qt_lo, qt_hi, &
                    qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                    lo, hi, domlo, domhi, &
                    dx, dt)

       ! we should not land here with radiation
       call tracez(q, qd_lo, qd_hi, &
                   qaux, qa_lo, qa_hi, &
                   dqz, qt_lo, qt_hi, &
                   qzm, qzp, qt_lo, qt_hi, &
                   lo, hi, domlo, domhi, &
                   dx, dt)

    end if  ! ppm test

    ! Compute \tilde{F}^x at kc (k3d)
    call cmpflx(qxm, qxp, qt_lo, qt_hi, &
                fx, fx_lo, fx_hi, &
                qgdnvx, qt_lo, qt_hi, &
#ifdef RADIATION
                rfx, fx_lo, fx_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! Compute \tilde{F}^y at kc (k3d)
    call cmpflx(qym, qyp, qt_lo, qt_hi, &
                fy, fy_lo, fy_hi, &
                qgdnvy, qt_lo, qt_hi, &
#ifdef RADIATION
                rfy, fy_lo, fy_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! Compute U'^y_x at kc (k3d)
    call transy1(qxm, qmxy, qxp, qpxy, qt_lo, qt_hi, &
                 qaux, qa_lo, qa_hi, &
                 fy, &
#ifdef RADIATION
                 rfy, &
#endif
                 fy_lo, fy_hi, &
                 qgdnvy, qt_lo, qt_hi, &
                 cdtdy, [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1])

    ! Compute U'^x_y at kc (k3d)
    call transx1(qym, qmyx, qyp, qpyx, qt_lo, qt_hi, &
                 qaux, qa_lo, qa_hi, &
                 fx, &
#ifdef RADIATION
                 rfx, &
#endif
                 fx_lo, fx_hi, &
                 qgdnvx, qt_lo, qt_hi, &
                 cdtdx, [lo(1), lo(2)-1, lo(3)-1], [hi(1), hi(2)+1, hi(3)+1])

    ! Compute F^{x|y} at kc (k3d)
    call cmpflx(qmxy, qpxy, qt_lo, qt_hi, &
                fxy, fx_lo, fx_hi, &
                qgdnvtmpx, qt_lo, qt_hi, &
#ifdef RADIATION
                rfxy, fx_lo, fx_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk,shk_lo,shk_hi, &
                1, [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1], &
                domlo, domhi)

       ! Compute F^{y|x} at kc (k3d)
    call cmpflx(qmyx, qpyx, qt_lo, qt_hi, &
                fyx, fy_lo, fy_hi, &
                qgdnvtmpy, qt_lo, qt_hi, &
#ifdef RADIATION
                rfyx, fy_lo, fy_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1], &
                domlo, domhi)


    ! Compute \tilde{F}^z at kc (k3d)
    call cmpflx(qzm,qzp,qt_lo,qt_hi, &
                fz,fz_lo,fz_hi, &
                qgdnvz,qt_lo,qt_hi, &
#ifdef RADIATION
                rfz, fz_lo, fz_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                3, [lo(1)-1, lo(2)-1, lo(3)], [ hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! Compute U'^y_z at kc (k3d)
    call transy2(qzm, qmzy, qzp, qpzy, qt_lo, qt_hi, &
                 qaux, qa_lo, qa_hi, &
                 fy, &
#ifdef RADIATION
                 rfy, &
#endif
                 fy_lo, fy_hi, &
                 qgdnvy, qt_lo, qt_hi, &
                 cdtdy, [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1])

    ! Compute U'^x_z at kc (k3d)
    call transx2(qzm, qmzx, qzp, qpzx, qt_lo, qt_hi, &
                 qaux, qa_lo, qa_hi, &
                 fx, &
#ifdef RADIATION
                 rfx, &
#endif
                 fx_lo, fx_hi, &
                 qgdnvx, qt_lo, qt_hi, &
                 cdtdx, [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1])

    ! Compute F^{z|x} at kc (k3d)
    call cmpflx(qmzx, qpzx, qt_lo, qt_hi, &
                fzx, fz_lo, fz_hi, &
                qgdnvtmpz1, qt_lo, qt_hi, &
#ifdef RADIATION
                rfzx, fz_lo, fz_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                3, [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! Compute F^{z|y} at kc (k3d)
    call cmpflx(qmzy, qpzy, qt_lo, qt_hi, &
                fzy, fz_lo, fz_hi, &
                qgdnvtmpz2, qt_lo, qt_hi, &
#ifdef RADIATION
                rfzy, fz_lo, fz_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                3, [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1], &
                domlo, domhi)

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
                 hdt, hdtdx, hdtdy, [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1])

    ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
    call cmpflx(qzl, qzr, qt_lo, qt_hi, &
                flux3, fd3_lo, fd3_hi, &
                qgdnvzf, qt_lo, qt_hi, &
#ifdef RADIATION
                rflux3, rfd3_lo, rfd3_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk,shk_lo,shk_hi, &
                3, [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1], &
                domlo, domhi)

    do k = lo(3), hi(3)+1
       do j=lo(2)-1,hi(2)+1
          do i=lo(1)-1,hi(1)+1
             q3(i,j,k,:) = qgdnvzf(i,j,k,:)
          end do
       end do
    end do

    ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
    call transz(qxm, qmxz, qxp, qpxz, qym, qmyz, qyp, qpyz, qt_lo, qt_hi, &
                qaux, qa_lo, qa_hi, &
                fz, &
#ifdef RADIATION
                rfz, &
#endif
                fz_lo, fz_hi, &
                qgdnvz,qt_lo,qt_hi, &
                cdtdz, [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1])

    ! Compute F^{x|z} at km (k3d-1)
    call cmpflx(qmxz, qpxz, qt_lo, qt_hi, &
                fxz, fx_lo, fx_hi, &
                qgdnvx, qt_lo, qt_hi, &
#ifdef RADIATION
                rfxz, fx_lo, fx_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! Compute F^{y|z} at km (k3d-1)
    call cmpflx(qmyz, qpyz, qt_lo, qt_hi, &
                fyz, fy_lo, fy_hi, &
                qgdnvy, qt_lo, qt_hi, &
#ifdef RADIATION
                rfyz, fy_lo, fy_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

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
                 hdt, hdtdy, hdtdz, [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1])

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
                 hdt, hdtdx, hdtdz, [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1])

    ! Compute F^x at km (k3d-1)
    call cmpflx(qxl, qxr, qt_lo, qt_hi, &
                flux1, fd1_lo, fd1_hi, &
                qgdnvxf, qt_lo, qt_hi, &
#ifdef RADIATION
                rflux1, rfd1_lo, rfd1_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], domlo, domhi)

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i=lo(1)-1, hi(1)+2
             q1(i,j,k,:) = qgdnvxf(i,j,k,:)
          end do
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
                2, [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], domlo, domhi)

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+2
          do i = lo(1)-1, hi(1)+1
             q2(i,j,k,:) = qgdnvyf(i,j,k,:)
          end do
       end do
    end do


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
                     cdtdx, lo, hi)

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


    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                compn = cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))

                if (j >= lo(2)+1) then
                   rr = qyp(i,j,k,QRHO)
                   rrnew = rr - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                   compu = rr*qyp(i,j,k,nqp) - compn
                   qypo(i,j,k,nqp) = compu/rrnew
                end if

                if (j <= hi(2)-1) then
                   rr = qym(i,j+1,k,QRHO)
                   rrnew = rr - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                   compu = rr*qym(i,j+1,k,nqp) - compn
                   qymo(i,j+1,k,nqp) = compu/rrnew
                end if

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse flux difference in the x-direction to y-states
             ! for the fluid variables
             !-------------------------------------------------------------------

             pggp  = qx(i+1,j,k,GDPRES)
             pggm  = qx(i  ,j,k,GDPRES)
             ugp  = qx(i+1,j,k,GDU   )
             ugm  = qx(i  ,j,k,GDU   )
             gegp = qx(i+1,j,k,GDGAME)
             gegm = qx(i  ,j,k,GDGAME)

#ifdef RADIATION
             lambda = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
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

             if (j >= lo(2)+1) then

                ! Convert to conservation form
                rrry = qyp(i,j,k,QRHO)
                rury = rrry*qyp(i,j,k,QU)
                rvry = rrry*qyp(i,j,k,QV)
                rwry = rrry*qyp(i,j,k,QW)
                ekenry = HALF*rrry* &
                     (qyp(i,j,k,QU)**2 + qyp(i,j,k,QV)**2 + qyp(i,j,k,QW)**2)
                rery = qyp(i,j,k,QREINT) + ekenry
#ifdef RADIATION
                err  = qyp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewry = rrry - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                runewry = rury - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
                rvnewry = rvry - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
                rwnewry = rwry - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
                renewry = rery - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
                runewry = runewry + dmom
                renewry = renewry + dre
                ernewr  = err(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) &
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

                qypo(i,j,k,QRHO) = rrnewry
                rhoinv = ONE/rrnewry
                qypo(i,j,k,QU) = runewry*rhoinv
                qypo(i,j,k,QV) = rvnewry*rhoinv
                qypo(i,j,k,QW) = rwnewry*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
                qypo(i,j,k,QREINT) = renewry - rhoekenry

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe .eq. 1 .and. qypo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating (rho e).
                      qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewry = qyp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                      qypo(i,j,k,QPRES) = max(pnewry, small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qypo(i,j,k,QPRES) = qypo(i,j,k,QREINT)*(qypo(i,j,k,QGAME)-ONE)
                      qypo(i,j,k,QPRES) = max(qypo(i,j,k,QPRES),small_pres)
                   end if
                else
                   qypo(i,j,k,QPRES) = qyp(i,j,k,QPRES)
                   qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qypo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
                qypo(i,j,k,qrad:qradhi) = ernewr(:)
                qypo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,k,QPRES)
                qypo(i,j,k,qreitot) = sum(qypo(i,j,k,qrad:qradhi)) + qypo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             if (j <= hi(2)-1) then

                ! Convert to conservation form
                rrly = qym(i,j+1,k,QRHO)
                ruly = rrly*qym(i,j+1,k,QU)
                rvly = rrly*qym(i,j+1,k,QV)
                rwly = rrly*qym(i,j+1,k,QW)
                ekenly = HALF*rrly* &
                     (qym(i,j+1,k,QU)**2 + qym(i,j+1,k,QV)**2 + qym(i,j+1,k,QW)**2)
                rely = qym(i,j+1,k,QREINT) + ekenly
#ifdef RADIATION
                erl  = qym(i,j+1,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewly = rrly - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                runewly = ruly - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
                rvnewly = rvly - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
                rwnewly = rwly - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
                renewly = rely - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
                runewly = runewly + dmom
                renewly = renewly + dre
                ernewl  = erl(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) &
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

                qymo(i,j+1,k,QRHO) = rrnewly
                rhoinv = ONE/rrnewly
                qymo(i,j+1,k,QU) = runewly*rhoinv
                qymo(i,j+1,k,QV) = rvnewly*rhoinv
                qymo(i,j+1,k,QW) = rwnewly*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
                qymo(i,j+1,k,QREINT) = renewly - rhoekenly

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary
                   if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qymo(i,j+1,k,QREINT) = qym(i,j+1,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewly = qym(i,j+1,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                      qymo(i,j+1,k,QPRES) = max(pnewly,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qymo(i,j+1,k,QPRES) = qymo(i,j+1,k,QREINT)*(qymo(i,j+1,k,QGAME)-ONE)
                      qymo(i,j+1,k,QPRES) = max(qymo(i,j+1,k,QPRES), small_pres)
                   end if
                else
                   qymo(i,j+1,k,QPRES) = qym(i,j+1,k,QPRES)
                   qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME)
                endif

                call reset_edge_state_thermo(qymo, qd_lo, qd_hi, i, j+1, k)

#ifdef RADIATION
                qymo(i,j+1,k,qrad:qradhi) = ernewl(:)
                qymo(i,j+1,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j+1,k,QPRES)
                qymo(i,j+1,k,qreitot) = sum(qymo(i,j+1,k,qrad:qradhi)) + qymo(i,j+1,k,QREINT)
#endif
             endif
          end do
       end do
    end do

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
                     cdtdx, lo, hi)


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


    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                compn = cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))

                rr = qzp(i,j,k,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = rr*qzp(i,j,k,nqp) - compn
                qzpo(i,j,k,nqp) = compu/rrnew

                compn = cdtdx*(fx(i+1,j,k-1,n) - fx(i,j,k-1,n))

                rr = qzm(i,j,k,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,k-1,URHO) - fx(i,j,k-1,URHO))
                compu = rr*qzm(i,j,k,nqp) - compn
                qzmo(i,j,k,nqp) = compu/rrnew

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse flux difference in the x-direction to z-states
             ! for the fluid variables
             !-------------------------------------------------------------------

             !-------------------------------------------------------------------
             ! qzpo state
             !-------------------------------------------------------------------

             pggp  = qx(i+1,j,k,GDPRES)
             pggm  = qx(i  ,j,k,GDPRES)
             ugp  = qx(i+1,j,k,GDU   )
             ugm  = qx(i  ,j,k,GDU   )
             gegp = qx(i+1,j,k,GDGAME)
             gegm = qx(i  ,j,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
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
             rrrz = qzp(i,j,k,QRHO)
             rurz = rrrz*qzp(i,j,k,QU)
             rvrz = rrrz*qzp(i,j,k,QV)
             rwrz = rrrz*qzp(i,j,k,QW)
             ekenrz = HALF*rrrz*(qzp(i,j,k,QU)**2 + qzp(i,j,k,QV)**2 + qzp(i,j,k,QW)**2)
             rerz = qzp(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
             err  = qzp(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewrz = rrrz - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
             runewrz = rurz - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
             rvnewrz = rvrz - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
             rwnewrz = rwrz - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
             renewrz = rerz - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
             runewrz = runewrz + dmom
             renewrz = renewrz + dre
             ernewr  = err(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) &
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
             qzpo(i,j,k,QRHO) = rrnewrz
             rhoinv = ONE/rrnewrz
             qzpo(i,j,k,QU) = runewrz*rhoinv
             qzpo(i,j,k,QV) = rvnewrz*rhoinv
             qzpo(i,j,k,QW) = rwnewrz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
             qzpo(i,j,k,QREINT) = renewrz - rhoekenrz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzpo(i,j,k,QREINT) = qzp(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrz = qzp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qzpo(i,j,k,QPRES) = max(pnewrz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzpo(i,j,k,QPRES) = qzpo(i,j,k,QREINT)*(qzpo(i,j,k,QGAME)-ONE)
                   qzpo(i,j,k,QPRES) = max(qzpo(i,j,k,QPRES), small_pres)
                endif
             else
                qzpo(i,j,k,QPRES) = qzp(i,j,k,QPRES)
                qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzpo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
             qzpo(i,j,k,qrad:qradhi) = ernewr(:)
             qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
             qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qradhi)) + qzpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             pggp  = qx(i+1,j,k-1,GDPRES)
             pggm  = qx(i  ,j,k-1,GDPRES)
             ugp  = qx(i+1,j,k-1,GDU   )
             ugm  = qx(i  ,j,k-1,GDU   )
             gegp = qx(i+1,j,k-1,GDGAME)
             gegm = qx(i  ,j,k-1,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j,k-1,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k-1,QGAMCG)
#else
             gamc = qaux(i,j,k-1,QGAMC)
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
             rrlz = qzm(i,j,k,QRHO)
             rulz = rrlz*qzm(i,j,k,QU)
             rvlz = rrlz*qzm(i,j,k,QV)
             rwlz = rrlz*qzm(i,j,k,QW)
             ekenlz = HALF*rrlz*(qzm(i,j,k,QU)**2 + qzm(i,j,k,QV)**2 + qzm(i,j,k,QW)**2)
             relz = qzm(i,j,k,QREINT) + ekenlz
#ifdef RADIATION
             erl  = qzm(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlz = rrlz - cdtdx*(fx(i+1,j,k-1,URHO) - fx(i,j,k-1,URHO))
             runewlz = rulz - cdtdx*(fx(i+1,j,k-1,UMX) - fx(i,j,k-1,UMX))
             rvnewlz = rvlz - cdtdx*(fx(i+1,j,k-1,UMY) - fx(i,j,k-1,UMY))
             rwnewlz = rwlz - cdtdx*(fx(i+1,j,k-1,UMZ) - fx(i,j,k-1,UMZ))
             renewlz = relz - cdtdx*(fx(i+1,j,k-1,UEDEN) - fx(i,j,k-1,UEDEN))
#ifdef RADIATION
             runewlz = runewlz + dmom
             renewlz = renewlz + dre
             ernewl  = erl(:) - cdtdx*(rfx(i+1,j,k-1,:) - rfx(i,j,k-1,:)) &
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
             qzmo(i,j,k,QRHO) = rrnewlz
             rhoinv = ONE/rrnewlz
             qzmo(i,j,k,QU) = runewlz*rhoinv
             qzmo(i,j,k,QV) = rvnewlz*rhoinv
             qzmo(i,j,k,QW) = rwnewlz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
             qzmo(i,j,k,QREINT) = renewlz - rhoekenlz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzmo(i,j,k,QREINT) = qzm(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k-1,UEINT) - fx(i,j,k-1,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   pnewlz = qzm(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qzmo(i,j,k,QPRES) = max(pnewlz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzmo(i,j,k,QPRES) = qzmo(i,j,k,QREINT)*(qzmo(i,j,k,QGAME)-ONE)
                   qzmo(i,j,k,QPRES) = max(qzmo(i,j,k,QPRES), small_pres)
                endif
             else
                qzmo(i,j,k,QPRES) = qzm(i,j,k,QPRES)
                qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzmo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
             qzmo(i,j,k,qrad:qradhi) = ernewl(:)
             qzmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k,QPRES)
             qzmo(i,j,k,qreitot) = sum(qzmo(i,j,k,qrad:qradhi)) + qzmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

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
                     cdtdy, lo, hi)


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


    integer, intent(in) :: qd_lo(3),qd_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fy_lo(3),fy_hi(3)
    integer, intent(in) :: qy_lo(3),qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                compn = cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))

                if (i >= lo(1)+1) then
                   rr = qxp(i,j,k,QRHO)
                   rrnew = rr - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                   compu = rr*qxp(i,j,k,nqp) - compn
                   qxpo(i,j,k,nqp) = compu/rrnew
                end if

                if (i <= hi(1)-1) then
                   rr = qxm(i+1,j,k,QRHO)
                   rrnew = rr - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                   compu = rr*qxm(i+1,j,k,nqp) - compn
                   qxmo(i+1,j,k,nqp) = compu/rrnew
                end if

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse flux difference in the y-direction to x-states
             ! for the fluid variables
             !-------------------------------------------------------------------

             pggp  = qy(i,j+1,k,GDPRES)
             pggm  = qy(i,j  ,k,GDPRES)
             ugp  = qy(i,j+1,k,GDV   )
             ugm  = qy(i,j  ,k,GDV   )
             gegp = qy(i,j+1,k,GDGAME)
             gegm = qy(i,j  ,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
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

             if (i >= lo(1)+1) then
                ! Convert to conservation form
                rrrx = qxp(i,j,k,QRHO)
                rurx = rrrx*qxp(i,j,k,QU)
                rvrx = rrrx*qxp(i,j,k,QV)
                rwrx = rrrx*qxp(i,j,k,QW)
                ekenrx = HALF*rrrx*(qxp(i,j,k,QU)**2 + qxp(i,j,k,QV)**2 &
                     + qxp(i,j,k,QW)**2)
                rerx = qxp(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
                err  = qxp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewrx = rrrx - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                runewrx = rurx - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
                rvnewrx = rvrx - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
                rwnewrx = rwrx - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
                renewrx = rerx - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewrx = rvnewrx + dmom
                renewrx = renewrx + dre
                ernewr = err(:) - cdtdy*(rfy(i,j+1,k,:) - rfy(i,j,k,:)) &
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

                qxpo(i,j,k,QRHO) = rrnewrx
                rhoinv = ONE/rrnewrx
                qxpo(i,j,k,QU) = runewrx*rhoinv
                qxpo(i,j,k,QV) = rvnewrx*rhoinv
                qxpo(i,j,k,QW) = rwnewrx*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
                qxpo(i,j,k,QREINT) = renewrx - rhoekenrx

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qxpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating (rho e).
                      qxpo(i,j,k,QREINT) = qxp(i,j,k,QREINT) - &
                           cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewrx = qxp(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                      qxpo(i,j,k,QPRES) = max(pnewrx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME) + &
                           cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxpo(i,j,k,QPRES) = qxpo(i,j,k,QREINT)*(qxpo(i,j,k,QGAME)-ONE)
                      qxpo(i,j,k,QPRES) = max(qxpo(i,j,k,QPRES), small_pres)
                   endif
                else
                   qxpo(i,j,k,QPRES) = qxp(i,j,k,QPRES)
                   qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qxpo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
                qxpo(i,j,k,qrad:qradhi) = ernewr(:)
                qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,k,QPRES)
                qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qradhi)) + qxpo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             if (i <= hi(1)-1) then
                ! Convert to conservation form
                rrlx = qxm(i+1,j,k,QRHO)
                rulx = rrlx*qxm(i+1,j,k,QU)
                rvlx = rrlx*qxm(i+1,j,k,QV)
                rwlx = rrlx*qxm(i+1,j,k,QW)
                ekenlx = HALF*rrlx*(qxm(i+1,j,k,QU)**2 + qxm(i+1,j,k,QV)**2 &
                     + qxm(i+1,j,k,QW)**2)
                relx = qxm(i+1,j,k,QREINT) + ekenlx
#ifdef RADIATION
                erl  = qxm(i+1,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewlx = rrlx - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                runewlx = rulx - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
                rvnewlx = rvlx - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
                rwnewlx = rwlx - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
                renewlx = relx - cdtdy*(fy(i,j+1,k,UEDEN)- fy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewlx = rvnewlx + dmom
                renewlx = renewlx + dre
                ernewl  = erl(:) + der(:)
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

                qxmo(i+1,j,k,QRHO) = rrnewlx
                rhoinv = ONE/rrnewlx
                qxmo(i+1,j,k,QU) = runewlx*rhoinv
                qxmo(i+1,j,k,QV) = rvnewlx*rhoinv
                qxmo(i+1,j,k,QW) = rwnewlx*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
                qxmo(i+1,j,k,QREINT) = renewlx - rhoekenlx

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qxmo(i+1,j,k,QREINT) = qxm(i+1,j,k,QREINT) - &
                           cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewlx = qxm(i+1,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                      qxmo(i+1,j,k,QPRES) = max(pnewlx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME) + &
                           cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxmo(i+1,j,k,QPRES) = qxmo(i+1,j,k,QREINT)*(qxmo(i+1,j,k,QGAME)-ONE)
                      qxmo(i+1,j,k,QPRES) = max(qxmo(i+1,j,k,QPRES), small_pres)
                   endif
                else
                   qxmo(i+1,j,k,QPRES) = qxm(i+1,j,k,QPRES)
                   qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qxmo, qd_lo, qd_hi, i+1, j, k)

#ifdef RADIATION
                qxmo(i+1,j,k,qrad:qradhi) = ernewl(:)
                qxmo(i+1,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i+1,j,k,QPRES)
                qxmo(i+1,j,k,qreitot) = sum(qxmo(i+1,j,k,qrad:qradhi)) + qxmo(i+1,j,k,QREINT)
#endif

             endif

          end do
       end do
    end do

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
                     cdtdy, lo, hi)


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


    integer, intent(in) :: qd_lo(3),qd_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fy_lo(3),fy_hi(3)
    integer, intent(in) :: qy_lo(3),qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                compn = cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))

                rr = qzp(i,j,k,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = rr*qzp(i,j,k,nqp) - compn
                qzpo(i,j,k,nqp) = compu/rrnew

                compn = cdtdy*(fy(i,j+1,k-1,n) - fy(i,j,k-1,n))

                rr = qzm(i,j,k,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,k-1,URHO) - fy(i,j,k-1,URHO))
                compu = rr*qzm(i,j,k,nqp) - compn
                qzmo(i,j,k,nqp) = compu/rrnew

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse flux difference in the y-direction to z-states
             ! for the fluid variables
             !-------------------------------------------------------------------

             !-------------------------------------------------------------------
             ! qzpo states
             !-------------------------------------------------------------------

             pggp  = qy(i,j+1,k,GDPRES)
             pggm  = qy(i,j  ,k,GDPRES)
             ugp  = qy(i,j+1,k,GDV   )
             ugm  = qy(i,j  ,k,GDV   )
             gegp = qy(i,j+1,k,GDGAME)
             gegm = qy(i,j  ,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
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
             rrrz = qzp(i,j,k,QRHO)
             rurz = rrrz*qzp(i,j,k,QU)
             rvrz = rrrz*qzp(i,j,k,QV)
             rwrz = rrrz*qzp(i,j,k,QW)
             ekenrz = HALF*rrrz*(qzp(i,j,k,QU)**2 + qzp(i,j,k,QV)**2 &
                  + qzp(i,j,k,QW)**2)
             rerz = qzp(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
             err  = qzp(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewrz = rrrz - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
             runewrz = rurz - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
             rvnewrz = rvrz - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
             rwnewrz = rwrz - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
             renewrz = rerz - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))
#ifdef RADIATION
             rvnewrz = rvnewrz + dmom
             renewrz = renewrz + dre
             ernewr  = err(:) - cdtdy*(rfy(i,j+1,k,:) - rfy(i,j,k,:)) &
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
             qzpo(i,j,k,QRHO) = rrnewrz
             rhoinv = ONE/rrnewrz
             qzpo(i,j,k,QU) = runewrz*rhoinv
             qzpo(i,j,k,QV) = rvnewrz*rhoinv
             qzpo(i,j,k,QW) = rwnewrz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
             qzpo(i,j,k,QREINT) = renewrz - rhoekenrz

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qzpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzpo(i,j,k,QREINT) = qzp(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrz = qzp(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qzpo(i,j,k,QPRES) = max(pnewrz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzpo(i,j,k,QPRES) = qzpo(i,j,k,QREINT)*(qzpo(i,j,k,QGAME)-ONE)
                   qzpo(i,j,k,QPRES) = max(qzpo(i,j,k,QPRES), small_pres)
                endif
             else
                qzpo(i,j,k,QPRES) = qzp(i,j,k,QPRES)
                qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzpo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
             qzpo(i,j,k,qrad:qradhi) = ernewr(:)
             qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
             qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qradhi)) + qzpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qzmo states
             !-------------------------------------------------------------------

             pggp  = qy(i,j+1,k-1,GDPRES)
             pggm  = qy(i,j  ,k-1,GDPRES)
             ugp  = qy(i,j+1,k-1,GDV   )
             ugm  = qy(i,j  ,k-1,GDV   )
             gegp = qy(i,j+1,k-1,GDGAME)
             gegm = qy(i,j  ,k-1,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i,j+1,k-1,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i,j  ,k-1,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k-1,QGAMCG)
#else
             gamc = qaux(i,j,k-1,QGAMC)
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
             rrlz = qzm(i,j,k,QRHO)
             rulz = rrlz*qzm(i,j,k,QU)
             rvlz = rrlz*qzm(i,j,k,QV)
             rwlz = rrlz*qzm(i,j,k,QW)
             ekenlz = HALF*rrlz*(qzm(i,j,k,QU)**2 + qzm(i,j,k,QV)**2 &
                  + qzm(i,j,k,QW)**2)
             relz = qzm(i,j,k,QREINT) + ekenlz
#ifdef RADIATION
             erl  = qzm(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlz = rrlz - cdtdy*(fy(i,j+1,k-1,URHO) - fy(i,j,k-1,URHO))
             runewlz = rulz - cdtdy*(fy(i,j+1,k-1,UMX) - fy(i,j,k-1,UMX))
             rvnewlz = rvlz - cdtdy*(fy(i,j+1,k-1,UMY) - fy(i,j,k-1,UMY))
             rwnewlz = rwlz - cdtdy*(fy(i,j+1,k-1,UMZ) - fy(i,j,k-1,UMZ))
             renewlz = relz - cdtdy*(fy(i,j+1,k-1,UEDEN)- fy(i,j,k-1,UEDEN))
#ifdef RADIATION
             rvnewlz = rvnewlz + dmom
             renewlz = renewlz + dre
             ernewl  = erl(:) - cdtdy*(rfy(i,j+1,k-1,:)- rfy(i,j,k-1,:)) &
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
             qzmo(i,j,k,QRHO) = rrnewlz
             rhoinv = ONE/rrnewlz
             qzmo(i,j,k,QU) = runewlz*rhoinv
             qzmo(i,j,k,QV) = rvnewlz*rhoinv
             qzmo(i,j,k,QW) = rwnewlz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
             qzmo(i,j,k,QREINT) = renewlz - rhoekenlz

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qzmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzmo(i,j,k,QREINT) = qzm(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k-1,UEINT) - fy(i,j,k-1,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qzmo(i,j,k,QPRES) = max(pnewlz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzmo(i,j,k,QPRES) = qzmo(i,j,k,QREINT)*(qzmo(i,j,k,QGAME)-ONE)
                   qzmo(i,j,k,QPRES) = max(qzmo(i,j,k,QPRES), small_pres)
                endif
             else
                qzmo(i,j,k,QPRES) = qzm(i,j,k,QPRES)
                qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzmo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
             qzmo(i,j,k,qrad:qradhi) = ernewl(:)
             qzmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k,QPRES)
             qzmo(i,j,k,qreitot) = sum(qzmo(i,j,k,qrad:qradhi)) + qzmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

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
                    cdtdz,lo, hi)


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


    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    integer, intent(in) :: qz_lo(3), qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer n, nqp, i, j, k, ipassive

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

        do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                compn = cdtdz*(fz(i,j,k,n) - fz(i,j,k-1,n))

                if (i >= lo(1)+1) then
                   rr = qxp(i,j,k-1,QRHO)
                   rrnew = rr - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                   compu = rr*qxp(i,j,k-1,nqp) - compn
                   qxpo(i,j,k-1,nqp) = compu/rrnew
                end if

                if (j >= lo(2)+1) then
                   rr = qyp(i,j,k-1,QRHO)
                   rrnew = rr - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                   compu = rr*qyp(i,j,k-1,nqp) - compn
                   qypo(i,j,k-1,nqp) = compu/rrnew
                end if

                if (i <= hi(1)-1) then
                   rr = qxm(i+1,j,k-1,QRHO)
                   rrnew = rr - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                   compu = rr*qxm(i+1,j,k-1,nqp) - compn
                   qxmo(i+1,j,k-1,nqp) = compu/rrnew
                end if

                if (j <= hi(2)-1) then
                   rr = qym(i,j+1,k-1,QRHO)
                   rrnew = rr - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                   compu = rr*qym(i,j+1,k-1,nqp) - compn
                   qymo(i,j+1,k-1,nqp) = compu/rrnew
                end if

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add transverse flux difference in the z-direction to the x- and
             ! y-states for the fluid variables
             !-------------------------------------------------------------------

             pggp  = qz(i,j,k,GDPRES)
             pggm  = qz(i,j,k-1,GDPRES)
             ugp  = qz(i,j,k,GDW   )
             ugm  = qz(i,j,k-1,GDW   )
             gegp = qz(i,j,k,GDGAME)
             gegm = qz(i,j,k-1,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             ergp = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
             ergm = qz(i,j,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             dup = pggp*ugp - pggm*ugm
             pav = HALF*(pggp+pggm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k-1,QGAMCG)
#else
             gamc = qaux(i,j,k-1,QGAMC)
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

             if (i >= lo(1)+1) then
                ! Convert to conservation form
                rrrx = qxp(i,j,k-1,QRHO)
                rurx = rrrx*qxp(i,j,k-1,QU)
                rvrx = rrrx*qxp(i,j,k-1,QV)
                rwrx = rrrx*qxp(i,j,k-1,QW)
                ekenrx = HALF*rrrx*(qxp(i,j,k-1,QU)**2 + qxp(i,j,k-1,QV)**2 &
                     + qxp(i,j,k-1,QW)**2)
                rerx = qxp(i,j,k-1,QREINT) + ekenrx
#ifdef RADIATION
                errx = qxp(i,j,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewrx = rrrx - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                runewrx = rurx - cdtdz*(fz(i,j,k,UMX) - fz(i,j,k-1,UMX))
                rvnewrx = rvrx - cdtdz*(fz(i,j,k,UMY) - fz(i,j,k-1,UMY))
                rwnewrx = rwrx - cdtdz*(fz(i,j,k,UMZ) - fz(i,j,k-1,UMZ))
                renewrx = rerx - cdtdz*(fz(i,j,k,UEDEN) - fz(i,j,k-1,UEDEN))
#ifdef RADIATION
                rwnewrx = rwnewrx + dmz
                renewrx = renewrx + dre
                ernewrx = errx(:) - cdtdz*(rfz(i,j,k,:) - rfz(i,j,k-1,:)) &
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

                qxpo(i,j,k-1,QRHO) = rrnewrx
                qxpo(i,j,k-1,QU) = runewrx/qxpo(i,j,k-1,QRHO)
                qxpo(i,j,k-1,QV) = rvnewrx/qxpo(i,j,k-1,QRHO)
                qxpo(i,j,k-1,QW) = rwnewrx/qxpo(i,j,k-1,QRHO)

                ! note: we run the risk of (rho e) being negative here
                rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,k-1,QRHO)
                qxpo(i,j,k-1,QREINT) = renewrx - rhoekenrx

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qxpo(i,j,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qxpo(i,j,k-1,QREINT) = qxp(i,j,k-1,QREINT) - &
                           cdtdz*(fz(i,j,k,UEINT) - fz(i,j,k-1,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewrx = qxp(i,j,k-1,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qxpo(i,j,k-1,QPRES) = max(pnewrx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxpo(i,j,k-1,QGAME) = qxp(i,j,k-1,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxpo(i,j,k-1,QPRES) = qxpo(i,j,k-1,QREINT)*(qxpo(i,j,k-1,QGAME)-ONE)
                      qxpo(i,j,k-1,QPRES) = max(qxpo(i,j,k-1,QPRES), small_pres)
                   endif
                else
                   qxpo(i,j,k-1,QPRES) = qxp(i,j,k-1,QPRES)
                   qxpo(i,j,k-1,QGAME) = qxp(i,j,k-1,QGAME)
                endif

                call reset_edge_state_thermo(qxpo, qd_lo, qd_hi, i, j, k-1)

#ifdef RADIATION
                qxpo(i,j,k-1,qrad:qradhi) = ernewrx(:)
                qxpo(i,j,k-1,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,k-1,QPRES)
                qxpo(i,j,k-1,qreitot) = sum(qxpo(i,j,k-1,qrad:qradhi)) + qxpo(i,j,k-1,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qypo state
             !-------------------------------------------------------------------

             if (j >= lo(2)+1) then
                ! Convert to conservation form
                rrry = qyp(i,j,k-1,QRHO)
                rury = rrry*qyp(i,j,k-1,QU)
                rvry = rrry*qyp(i,j,k-1,QV)
                rwry = rrry*qyp(i,j,k-1,QW)
                ekenry = HALF*rrry*(qyp(i,j,k-1,QU)**2 + qyp(i,j,k-1,QV)**2 &
                     + qyp(i,j,k-1,QW)**2)
                rery = qyp(i,j,k-1,QREINT) + ekenry
#ifdef RADIATION
                erry = qyp(i,j,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewry = rrry - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                runewry = rury - cdtdz*(fz(i,j,k,UMX) - fz(i,j,k-1,UMX))
                rvnewry = rvry - cdtdz*(fz(i,j,k,UMY) - fz(i,j,k-1,UMY))
                rwnewry = rwry - cdtdz*(fz(i,j,k,UMZ) - fz(i,j,k-1,UMZ))
                renewry = rery - cdtdz*(fz(i,j,k,UEDEN) - fz(i,j,k-1,UEDEN))
#ifdef RADIATION
                rwnewry = rwnewry + dmz
                renewry = renewry + dre
                ernewry = erry(:) - cdtdz*(rfz(i,j,k,:) - rfz(i,j,k-1,:)) &
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

                qypo(i,j,k-1,QRHO) = rrnewry
                qypo(i,j,k-1,QU) = runewry/qypo(i,j,k-1,QRHO)
                qypo(i,j,k-1,QV) = rvnewry/qypo(i,j,k-1,QRHO)
                qypo(i,j,k-1,QW) = rwnewry/qypo(i,j,k-1,QRHO)

                ! note: we run the risk of (rho e) being negative here
                rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,k-1,QRHO)
                qypo(i,j,k-1,QREINT) = renewry - rhoekenry

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qypo(i,j,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qypo(i,j,k-1,QREINT) = qyp(i,j,k-1,QREINT) - &
                           cdtdz*(fz(i,j,k,UEINT) - fz(i,j,k-1,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewry = qyp(i,j,k-1,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qypo(i,j,k-1,QPRES) = max(pnewry,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qypo(i,j,k-1,QGAME) = qyp(i,j,k-1,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qypo(i,j,k-1,QPRES) = qypo(i,j,k-1,QREINT)*(qypo(i,j,k-1,QGAME)-ONE)
                      qypo(i,j,k-1,QPRES) = max(qypo(i,j,k-1,QPRES), small_pres)
                   endif
                else
                   qypo(i,j,k-1,QPRES) = qyp(i,j,k-1,QPRES)
                   qypo(i,j,k-1,QGAME) = qyp(i,j,k-1,QGAME)
                endif

                call reset_edge_state_thermo(qypo, qd_lo, qd_hi, i, j, k-1)

#ifdef RADIATION
                qypo(i,j,k-1,qrad:qradhi) = ernewry(:)
                qypo(i,j,k-1,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,k-1,QPRES)
                qypo(i,j,k-1,qreitot) = sum(qypo(i,j,k-1,qrad:qradhi)) + qypo(i,j,k-1,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             if (i <= hi(1)-1) then

                ! Convert to conservation form
                rrlx = qxm(i+1,j,k-1,QRHO)
                rulx = rrlx*qxm(i+1,j,k-1,QU)
                rvlx = rrlx*qxm(i+1,j,k-1,QV)
                rwlx = rrlx*qxm(i+1,j,k-1,QW)
                ekenlx = HALF*rrlx*(qxm(i+1,j,k-1,QU)**2 + qxm(i+1,j,k-1,QV)**2 &
                     + qxm(i+1,j,k-1,QW)**2)
                relx = qxm(i+1,j,k-1,QREINT) + ekenlx
#ifdef RADIATION
                erlx = qxm(i+1,j,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewlx = rrlx - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                runewlx = rulx - cdtdz*(fz(i,j,k,UMX) - fz(i,j,k-1,UMX))
                rvnewlx = rvlx - cdtdz*(fz(i,j,k,UMY) - fz(i,j,k-1,UMY))
                rwnewlx = rwlx - cdtdz*(fz(i,j,k,UMZ) - fz(i,j,k-1,UMZ))
                renewlx = relx - cdtdz*(fz(i,j,k,UEDEN) - fz(i,j,k-1,UEDEN))
#ifdef RADIATION
                rwnewlx = rwnewlx + dmz
                renewlx = renewlx + dre
                ernewlx = erlx(:) - cdtdz*(rfz(i,j,k,:) - rfz(i,j,k-1,:)) &
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

                qxmo(i+1,j,k-1,QRHO) = rrnewlx
                qxmo(i+1,j,k-1,QU) = runewlx/qxmo(i+1,j,k-1,QRHO)
                qxmo(i+1,j,k-1,QV) = rvnewlx/qxmo(i+1,j,k-1,QRHO)
                qxmo(i+1,j,k-1,QW) = rwnewlx/qxmo(i+1,j,k-1,QRHO)

                ! note: we run the risk of (rho e) being negative here
                rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,k-1,QRHO)
                qxmo(i+1,j,k-1,QREINT) = renewlx - rhoekenlx

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qxmo(i+1,j,k-1,QREINT) = qxm(i+1,j,k-1,QREINT) - &
                           cdtdz*(fz(i,j,k,UEINT) - fz(i,j,k-1,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewlx = qxm(i+1,j,k-1,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qxmo(i+1,j,k-1,QPRES) = max(pnewlx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxmo(i+1,j,k-1,QGAME) = qxm(i+1,j,k-1,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxmo(i+1,j,k-1,QPRES) = qxmo(i+1,j,k-1,QREINT)*(qxmo(i+1,j,k-1,QGAME)-ONE)
                      qxmo(i+1,j,k-1,QPRES) = max(qxmo(i+1,j,k-1,QPRES), small_pres)
                   end if
                else
                   qxmo(i+1,j,k-1,QPRES) = qxm(i+1,j,k-1,QPRES)
                   qxmo(i+1,j,k-1,QGAME) = qxm(i+1,j,k-1,QGAME)
                endif

                call reset_edge_state_thermo(qxmo, qd_lo, qd_hi, i+1, j, k-1)

#ifdef RADIATION
                qxmo(i+1,j,k-1,qrad:qradhi) = ernewlx(:)
                qxmo(i+1,j,k-1,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i+1,j,k-1,QPRES)
                qxmo(i+1,j,k-1,qreitot) = sum(qxmo(i+1,j,k-1,qrad:qradhi)) + qxmo(i+1,j,k-1,QREINT)
#endif

             endif


             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             if (j <= hi(2)-1) then
                ! Convert to conservation form
                rrly = qym(i,j+1,k-1,QRHO)
                ruly = rrly*qym(i,j+1,k-1,QU)
                rvly = rrly*qym(i,j+1,k-1,QV)
                rwly = rrly*qym(i,j+1,k-1,QW)
                ekenly = HALF*rrly*(qym(i,j+1,k-1,QU)**2 + qym(i,j+1,k-1,QV)**2 &
                     + qym(i,j+1,k-1,QW)**2)
                rely = qym(i,j+1,k-1,QREINT) + ekenly
#ifdef RADIATION
                erly = qym(i,j+1,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewly = rrly - cdtdz*(fz(i,j,k,URHO) - fz(i,j,k-1,URHO))
                runewly = ruly - cdtdz*(fz(i,j,k,UMX) - fz(i,j,k-1,UMX))
                rvnewly = rvly - cdtdz*(fz(i,j,k,UMY) - fz(i,j,k-1,UMY))
                rwnewly = rwly - cdtdz*(fz(i,j,k,UMZ) - fz(i,j,k-1,UMZ))
                renewly = rely - cdtdz*(fz(i,j,k,UEDEN) - fz(i,j,k-1,UEDEN))
#ifdef RADIATION
                rwnewly = rwnewly + dmz
                renewly = renewly + dre
                ernewly = erly(:) - cdtdz*(rfz(i,j,k,:) - rfz(i,j,k-1,:)) &
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

                qymo(i,j+1,k-1,QRHO) = rrnewly
                qymo(i,j+1,k-1,QU) = runewly/qymo(i,j+1,k-1,QRHO)
                qymo(i,j+1,k-1,QV) = rvnewly/qymo(i,j+1,k-1,QRHO)
                qymo(i,j+1,k-1,QW) = rwnewly/qymo(i,j+1,k-1,QRHO)

                ! note: we run the risk of (rho e) being negative here
                rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,k-1,QRHO)
                qymo(i,j+1,k-1,QREINT) = renewly - rhoekenly

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qymo(i,j+1,k-1,QREINT) = qym(i,j+1,k-1,QREINT) - &
                           cdtdz*(fz(i,j,k,UEINT) - fz(i,j,k-1,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewly = qym(i,j+1,k-1,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qymo(i,j+1,k-1,QPRES) = max(pnewly,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qymo(i,j+1,k-1,QGAME) = qym(i,j+1,k-1,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qymo(i,j+1,k-1,QPRES) = qymo(i,j+1,k-1,QREINT)*(qymo(i,j+1,k-1,QGAME)-ONE)
                      qymo(i,j+1,k-1,QPRES) = max(qymo(i,j+1,k-1,QPRES), small_pres)
                   endif
                else
                   qymo(i,j+1,k-1,QPRES) = qym(i,j+1,k-1,QPRES)
                   qymo(i,j+1,k-1,QGAME) = qym(i,j+1,k-1,QGAME)
                endif

                call reset_edge_state_thermo(qymo, qd_lo, qd_hi, i, j+1, k-1)

#ifdef RADIATION
                qymo(i,j+1,k-1,qrad:qradhi) = ernewly(:)
                qymo(i,j+1,k-1,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j+1,k-1,QPRES)
                qymo(i,j+1,k-1,qreitot) = sum(qymo(i,j+1,k-1,qrad:qradhi)) + qymo(i,j+1,k-1,QREINT)
#endif

             endif

          end do
       end do
    end do
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
                     hdt,cdtdx,cdtdy,lo,hi)


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


    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rrr = qp(i,j,k,QRHO)
                rrl = qm(i,j,k,QRHO)

                compr = rrr*qp(i,j,k,nqp)
                compl = rrl*qm(i,j,k,nqp)

                rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                             - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))

                compnr = compr - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                               - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))
                compnl = compl - cdtdx*(fxy(i+1,j,k-1,n) - fxy(i,j,k-1,n)) &
                               - cdtdy*(fyx(i,j+1,k-1,n) - fyx(i,j,k-1,n))

                qpo(i,j,k,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k  ,nqp)
                qmo(i,j,k,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k-1,nqp)

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse xy and yx differences to the z-states for the
             ! fluid variables
             !-------------------------------------------------------------------

             pggxp  = qx(i+1,j,k,GDPRES)
             pggxm  = qx(i  ,j,k,GDPRES)
             ugxp  = qx(i+1,j,k,GDU   )
             ugxm  = qx(i  ,j,k,GDU   )
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i  ,j,k,GDGAME)

             pggyp  = qy(i,j+1,k,GDPRES)
             pggym  = qy(i,j  ,k,GDPRES)
             ugyp  = qy(i,j+1,k,GDV   )
             ugym  = qy(i,j  ,k,GDV   )
             gegyp = qy(i,j+1,k,GDGAME)
             gegym = qy(i,j  ,k,GDGAME)

#ifdef RADIATION
             lamc(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ergxp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
             ergyp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             pggxpm  = qx(i+1,j,k-1,GDPRES)
             pggxmm  = qx(i  ,j,k-1,GDPRES)
             ugxpm  = qx(i+1,j,k-1,GDU   )
             ugxmm  = qx(i  ,j,k-1,GDU   )
             gegxpm = qx(i+1,j,k-1,GDGAME)
             gegxmm = qx(i  ,j,k-1,GDGAME)

             pggypm  = qy(i,j+1,k-1,GDPRES)
             pggymm  = qy(i,j  ,k-1,GDPRES)
             ugypm  = qy(i,j+1,k-1,GDV   )
             ugymm  = qy(i,j  ,k-1,GDV   )
             gegypm = qy(i,j+1,k-1,GDGAME)
             gegymm = qy(i,j  ,k-1,GDGAME)

#ifdef RADIATION
             lamm(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             ergxpm = qx(i+1,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergxmm = qx(i  ,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergypm = qy(i,j+1,k-1,GDERADS:GDERADS-1+ngroups)
             ergymm = qy(i,j  ,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pggxp*ugxp - pggxm*ugxm
             pxav = HALF*(pggxp+pggxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )
#endif


             duxpm = pggxpm*ugxpm - pggxmm*ugxmm
             pxavm = HALF*(pggxpm+pggxmm)
             uxavm = HALF*(ugxpm+ugxmm)
             gexavm = HALF*(gegxpm+gegxmm)
             duxm = ugxpm-ugxmm
             dgexm = gegxpm-gegxmm
#ifdef RADIATION
             pxnewm = cdtdx*(duxpm + pxavm*duxm*(qaux(i,j,k-1,QGAMCG) - ONE))
             gexnewm = cdtdx*( (gexavm-ONE)*(gexavm - qaux(i,j,k-1,QGAMCG))*duxm - uxavm*dgexm )
#else
             pxnewm = cdtdx*(duxpm + pxavm*duxm*(qaux(i,j,k-1,QGAMC) - ONE))
             gexnewm = cdtdx*( (gexavm-ONE)*(gexavm - qaux(i,j,k-1,QGAMC))*duxm - uxavm*dgexm )
#endif
             duyp = pggyp*ugyp - pggym*ugym
             pyav = HALF*(pggyp+pggym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMC))*duy - uyav*dgey )
#endif
             duypm = pggypm*ugypm - pggymm*ugymm
             pyavm = HALF*(pggypm+pggymm)
             uyavm = HALF*(ugypm+ugymm)
             geyavm = HALF*(gegypm+gegymm)
             duym = ugypm-ugymm
             dgeym = gegypm-gegymm
#ifdef RADIATION
             pynewm = cdtdy*(duypm + pyavm*duym*(qaux(i,j,k-1,QGAMCG) - ONE))
             geynewm = cdtdy*( (geyavm-ONE)*(geyavm - qaux(i,j,k-1,QGAMCG))*duym - uyavm*dgeym )
#else
             pynewm = cdtdy*(duypm + pyavm*duym*(qaux(i,j,k-1,QGAMC) - ONE))
             geynewm = cdtdy*( (geyavm-ONE)*(geyavm - qaux(i,j,k-1,QGAMC))*duym - uyavm*dgeym )
#endif

             !-------------------------------------------------------------------
             ! qzpo state
             !-------------------------------------------------------------------

             ! Convert to conservation form
             rrr = qp(i,j,k,QRHO)
             rur = rrr*qp(i,j,k,QU)
             rvr = rrr*qp(i,j,k,QV)
             rwr = rrr*qp(i,j,k,QW)
             ekenr = HALF*rrr*(qp(i,j,k,QU)**2 + qp(i,j,k,QV)**2 + &
                  qp(i,j,k,QW)**2)
             rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,k,qrad:qradhi)

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
             rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                          - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
             runewr = rur - cdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                          - cdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
             rvnewr = rvr - cdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                          - cdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
             rwnewr = rwr - cdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
             renewr = rer - cdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))
#ifdef RADIATION
             runewr = runewr + dmx
             rvnewr = rvnewr + dmy
             renewr = renewr + dre
             ernewr = err(:) - cdtdx*(rfxy(i+1,j,k,:) - rfxy(i,j,k,:)) &
                             - cdtdy*(rfyx(i,j+1,k,:) - rfyx(i,j,k,:))  &
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

             ! Convert back to primitive form
             qpo(i,j,k,QRHO  ) = rrnewr
             qpo(i,j,k,QU    ) = runewr/rrnewr
             qpo(i,j,k,QV    ) = rvnewr/rrnewr
             qpo(i,j,k,QW    ) = rwnewr/rrnewr

             ! for ppm_type > 0 we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_type == 0) then
                qpo(i,j,k,QRHO  ) = qpo(i,j,k,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                qpo(i,j,k,QU:QW) = qpo(i,j,k,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr
             if (ppm_type == 0) then
                qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy)
                endif

                if (ppm_type == 0) then
                   qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pxnew - pynew
                   qpo(i,j,k,QPRES) = pnewr
                   if (ppm_type == 0) then
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geynew
                   
                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                endif
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                endif
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
             qpo(i,j,k,qrad:qradhi) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lamc(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             ! Convert to conservation form
             rrl = qm(i,j,k,QRHO)
             rul = rrl*qm(i,j,k,QU)
             rvl = rrl*qm(i,j,k,QV)
             rwl = rrl*qm(i,j,k,QW)
             ekenl = HALF*rrl*(qm(i,j,k,QU)**2 + qm(i,j,k,QV)**2 + &
                  qm(i,j,k,QW)**2)
             rel = qm(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i,j,k,qrad:qradhi)

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
             rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                          - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
             runewl = rul - cdtdx*(fxy(i+1,j,k-1,UMX) - fxy(i,j,k-1,UMX)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMX) - fyx(i,j,k-1,UMX))
             rvnewl = rvl - cdtdx*(fxy(i+1,j,k-1,UMY) - fxy(i,j,k-1,UMY)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMY) - fyx(i,j,k-1,UMY))
             rwnewl = rwl - cdtdx*(fxy(i+1,j,k-1,UMZ) - fxy(i,j,k-1,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMZ) - fyx(i,j,k-1,UMZ))
             renewl = rel - cdtdx*(fxy(i+1,j,k-1,UEDEN) - fxy(i,j,k-1,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k-1,UEDEN) - fyx(i,j,k-1,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rvnewl = rvnewl + dmy
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,k-1,:) - rfxy(i,j,k-1,:)) &
                             - cdtdy*(rfyx(i  ,j+1,k-1,:) - rfyx(i,j,k-1,:)) &
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

             qmo(i,j,k,QRHO  ) = rrnewl
             qmo(i,j,k,QU    ) = runewl/rrnewl
             qmo(i,j,k,QV    ) = rvnewl/rrnewl
             qmo(i,j,k,QW    ) = rwnewl/rrnewl

             ! for ppm_type > 0 we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_type == 0) then
                qmo(i,j,k,QRHO  ) = qmo(i,j,k,QRHO  ) + hdt*srcQ(i,j,k-1,QRHO)
                qmo(i,j,k,QU:QW) = qmo(i,j,k,QU:QW) + hdt * srcQ(i,j,k-1,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT) = renewl - rhoekenl
             if (ppm_type == 0) then
                qmo(i,j,k,QREINT) = qmo(i,j,k,QREINT) + hdt*srcQ(i,j,k-1,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k-1,UEINT) - fxy(i,j,k-1,UEINT) + pxavm*duxm) &
                        - cdtdy*(fyx(i,j+1,k-1,UEINT) - fyx(i,j,k-1,UEINT) + pyavm*duym)
                   if (ppm_type == 0) then
                      qmo(i,j,k,QREINT) = qmo(i,j,k,QREINT) + hdt*srcQ(i,j,k-1,QREINT)
                   endif
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j,k,QPRES) - pxnewm - pynewm
                   qmo(i,j,k,QPRES) = pnewl
                   if (ppm_type == 0) then
                      qmo(i,j,k,QPRES) = qmo(i,j,k,QPRES) + hdt*srcQ(i,j,k-1,QPRES)
                   endif
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + gexnewm + geynewm

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QREINT)*(qmo(i,j,k,QGAME)-ONE)
                endif
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                if (ppm_type == 0) then
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QPRES) + hdt*srcQ(i,j,k-1,QPRES)
                endif
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif
             
             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i, j, k)

#ifdef RADIATION
             qmo(i,j,k,qrad:qradhi) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lamm(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qradhi)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

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
                     hdt, cdtdx, cdtdz, lo, hi)


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


    integer, intent(in) :: qd_lo(3),qd_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fx_lo(3),fx_hi(3)
    integer, intent(in) :: fz_lo(3),fz_hi(3)
    integer, intent(in) :: qx_lo(3),qx_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                drr    = - cdtdx*(fxz(i+1,j,k-1,URHO) - fxz(i,j,k-1,URHO)) &
                         - cdtdz*(fzx(i  ,j,k,URHO) - fzx(i,j,k-1,URHO))
                dcompn = - cdtdx*(fxz(i+1,j,k-1,n   ) - fxz(i,j,k-1,n)) &
                         - cdtdz*(fzx(i  ,j,k,n   ) - fzx(i,j,k-1,n))

                if (j >= lo(2)+1) then
                   rrr = qp(i,j,k-1,QRHO)
                   compr = rrr*qp(i,j,k-1,nqp)

                   rrnewr = rrr + drr
                   compnr = compr + dcompn

                   qpo(i,j  ,k-1,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k,nqp)
                end if

                if (j <= hi(2)-1) then
                   rrl = qm(i,j+1,k-1,QRHO)
                   compl = rrl*qm(i,j+1,k-1,nqp)

                   rrnewl = rrl + drr
                   compnl = compl + dcompn

                   qmo(i,j+1,k-1,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k,nqp)
                end if

             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse xz and zx differences to the y-states for the
             ! fluid variables
             !-------------------------------------------------------------------

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
#endif

             pggxp  = qx(i+1,j,k-1,GDPRES)
             pggxm  = qx(i  ,j,k-1,GDPRES)
             ugxp  = qx(i+1,j,k-1,GDU   )
             ugxm  = qx(i  ,j,k-1,GDU   )
             gegxp = qx(i+1,j,k-1,GDGAME)
             gegxm = qx(i  ,j,k-1,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             pggzp  = qz(i,j,k,GDPRES)
             pggzm  = qz(i,j,k-1,GDPRES)
             ugzp  = qz(i,j,k,GDW   )
             ugzm  = qz(i,j,k-1,GDW   )
             gegzp = qz(i,j,k,GDGAME)
             gegzm = qz(i,j,k-1,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pggxp*ugxp - pggxm*ugxm
             pxav = HALF*(pggxp+pggxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )
#endif

             duzp = pggzp*ugzp - pggzm*ugzm
             pzav = HALF*(pggzp+pggzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMC))*duz - uzav*dgez )
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

             if (j >= lo(2)+1) then
                ! Convert to conservation form
                rrr = qp(i,j,k-1,QRHO)
                rur = rrr*qp(i,j,k-1,QU)
                rvr = rrr*qp(i,j,k-1,QV)
                rwr = rrr*qp(i,j,k-1,QW)
                ekenr = HALF*rrr*(qp(i,j,k-1,QU)**2 + qp(i,j,k-1,QV)**2 + qp(i,j,k-1,QW)**2)
                rer = qp(i,j,k-1,QREINT) + ekenr
#ifdef RADIATION
                err = qp(i,j,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewr = rrr - cdtdx*(fxz(i+1,j,k-1,URHO) - fxz(i,j,k-1,URHO)) &
                             - cdtdz*(fzx(i,j,k,URHO) - fzx(i,j,k-1,URHO))
                runewr = rur - cdtdx*(fxz(i+1,j,k-1,UMX) - fxz(i,j,k-1,UMX)) &
                             - cdtdz*(fzx(i,j,k,UMX) - fzx(i,j,k-1,UMX))
                rvnewr = rvr - cdtdx*(fxz(i+1,j,k-1,UMY) - fxz(i,j,k-1,UMY)) &
                             - cdtdz*(fzx(i,j,k,UMY) - fzx(i,j,k-1,UMY))
                rwnewr = rwr - cdtdx*(fxz(i+1,j,k-1,UMZ) - fxz(i,j,k-1,UMZ)) &
                             - cdtdz*(fzx(i,j,k,UMZ) - fzx(i,j,k-1,UMZ))
                renewr = rer - cdtdx*(fxz(i+1,j,k-1,UEDEN) - fxz(i,j,k-1,UEDEN)) &
                             - cdtdz*(fzx(i,j,k,UEDEN) - fzx(i,j,k-1,UEDEN))
#ifdef RADIATION
                runewr = runewr + dmx
                rwnewr = rwnewr + dmz
                renewr = renewr + dre
                ernewr = err(:) - cdtdx*(rfxz(i+1,j,k-1,:) - rfxz(i,j,k-1,:)) &
                                - cdtdz*(rfzx(i  ,j,k,:) - rfzx(i,j,k-1,:)) &
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

                qpo(i,j,k-1,QRHO  ) = rrnewr
                qpo(i,j,k-1,QU    ) = runewr/rrnewr
                qpo(i,j,k-1,QV    ) = rvnewr/rrnewr
                qpo(i,j,k-1,QW    ) = rwnewr/rrnewr
                
                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qpo(i,j,k-1,QRHO  ) = qpo(i,j,k-1,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qpo(i,j,k-1,QU:QW) = qpo(i,j,k-1,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                qpo(i,j,k-1,QREINT) = renewr - rhoekenr
                if (ppm_type == 0) then
                   qpo(i,j,k-1,QREINT) = qpo(i,j,k-1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qpo(i,j,k-1,QREINT) <= ZERO) then
                      qpo(i,j,k-1,QREINT) = qp(i,j,k-1,QREINT) &
                           - cdtdx*(fxz(i+1,j,k-1,UEINT) - fxz(i,j,k-1,UEINT) + pxav*dux) &
                           - cdtdz*(fzx(i  ,j,k,UEINT) - fzx(i,j,k-1,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qpo(i,j,k-1,QREINT) = qpo(i,j,k-1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewr = qp(i,j,k-1,QPRES) - pxnew - pznew
                      qpo(i,j,k-1,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,k-1,QPRES) = qpo(i,j,k-1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                   ! Update gammae with its transverse terms
                      qpo(i,j,k-1,QGAME) = qp(i,j,k-1,QGAME) + gexnew + geznew
                      
                      ! and compute the p edge state from this and (rho e)
                      qpo(i,j,k-1,QPRES) = qpo(i,j,k-1,QREINT)*(qpo(i,j,k-1,QGAME)-ONE)
                   endif
                else
                   qpo(i,j,k-1,QPRES) = qp(i,j,k-1,QPRES)
                   if (ppm_type == 0) then
                      qpo(i,j,k-1,QPRES) = qpo(i,j,k-1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qpo(i,j,k-1,QGAME) = qp(i,j,k-1,QGAME)
                endif

                qpo(i,j,k-1,QPRES) = max(qpo(i,j,k-1,QPRES), small_pres)

                call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, k-1)

#ifdef RADIATION
                qpo(i,j,k-1,qrad:qradhi) = ernewr(:)
                qpo(i,j,k-1,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k-1,QPRES)
                qpo(i,j,k-1,qreitot) = sum(qpo(i,j,k-1,qrad:qradhi)) + qpo(i,j,k-1,QREINT)
#endif
             end if


             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             if (j <= hi(2)-1) then
                ! Convert to conservation form
                rrl = qm(i,j+1,k-1,QRHO)
                rul = rrl*qm(i,j+1,k-1,QU)
                rvl = rrl*qm(i,j+1,k-1,QV)
                rwl = rrl*qm(i,j+1,k-1,QW)
                ekenl = HALF*rrl*(qm(i,j+1,k-1,QU)**2 + qm(i,j+1,k-1,QV)**2 + qm(i,j+1,k-1,QW)**2)
                rel = qm(i,j+1,k-1,QREINT) + ekenl
#ifdef RADIATION
                erl = qm(i,j+1,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewl = rrl - cdtdx*(fxz(i+1,j,k-1,URHO) - fxz(i,j,k-1,URHO)) &
                             - cdtdz*(fzx(i,j,k,URHO) - fzx(i,j,k-1,URHO))
                runewl = rul - cdtdx*(fxz(i+1,j,k-1,UMX) - fxz(i,j,k-1,UMX)) &
                             - cdtdz*(fzx(i,j,k,UMX) - fzx(i,j,k-1,UMX))
                rvnewl = rvl - cdtdx*(fxz(i+1,j,k-1,UMY) - fxz(i,j,k-1,UMY)) &
                             - cdtdz*(fzx(i,j,k,UMY) - fzx(i,j,k-1,UMY))
                rwnewl = rwl - cdtdx*(fxz(i+1,j,k-1,UMZ) - fxz(i,j,k-1,UMZ)) &
                             - cdtdz*(fzx(i,j,k,UMZ) - fzx(i,j,k-1,UMZ))
                renewl = rel - cdtdx*(fxz(i+1,j,k-1,UEDEN) - fxz(i,j,k-1,UEDEN)) &
                             - cdtdz*(fzx(i,j,k,UEDEN) - fzx(i,j,k-1,UEDEN))
#ifdef RADIATION
                runewl = runewl + dmx
                rwnewl = rwnewl + dmz
                renewl = renewl + dre
                ernewl = erl(:) - cdtdx*(rfxz(i+1,j,k-1,:) - rfxz(i,j,k-1,:)) &
                                - cdtdz*(rfzx(i  ,j,k,:) - rfzx(i,j,k-1,:)) &
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

                qmo(i,j+1,k-1,QRHO  ) = rrnewl
                qmo(i,j+1,k-1,QU    ) = runewl/rrnewl
                qmo(i,j+1,k-1,QV    ) = rvnewl/rrnewl
                qmo(i,j+1,k-1,QW    ) = rwnewl/rrnewl

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qmo(i,j+1,k-1,QRHO  ) = qmo(i,j+1,k-1,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qmo(i,j+1,k-1,QU:QW) = qmo(i,j+1,k-1,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
                qmo(i,j+1,k-1,QREINT) = renewl - rhoekenl
                if (ppm_type == 0) then
                   qmo(i,j+1,k-1,QREINT) = qmo(i,j+1,k-1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qmo(i,j+1,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qmo(i,j+1,k-1,QREINT) = qm(i,j+1,k-1,QREINT) &
                           - cdtdx*(fxz(i+1,j,k-1,UEINT) - fxz(i,j,k-1,UEINT) + pxav*dux) &
                           - cdtdz*(fzx(i,j,k,UEINT) - fzx(i,j,k-1,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qmo(i,j+1,k-1,QREINT) = qmo(i,j+1,k-1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewl = qm(i,j+1,k-1,QPRES) - pxnew - pznew
                      qmo(i,j+1,k-1,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i,j+1,k-1,QPRES) = qmo(i,j+1,k-1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qmo(i,j+1,k-1,QGAME) = qm(i,j+1,k-1,QGAME) + gexnew + geznew
                      
                      ! and compute the p edge state from this and (rho e)
                      qmo(i,j+1,k-1,QPRES) = qmo(i,j+1,k-1,QREINT)*(qmo(i,j+1,k-1,QGAME)-ONE)
                   endif
                else
                   qmo(i,j+1,k-1,QPRES) = qm(i,j+1,k-1,QPRES)
                   if (ppm_type == 0) then
                      qmo(i,j+1,k-1,QPRES) = qmo(i,j+1,k-1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qmo(i,j+1,k-1,QGAME) = qm(i,j+1,k-1,QGAME)
                endif
                
                qmo(i,j+1,k-1,QPRES) = max(qmo(i,j+1,k-1,QPRES), small_pres)

                call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i, j+1, k-1)

#ifdef RADIATION
                qmo(i,j+1,k-1,qrad:qradhi) = ernewl(:)
                qmo(i,j+1,k-1,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j+1,k-1,QPRES)
                qmo(i,j+1,k-1,qreitot) = sum(qmo(i,j+1,k-1,qrad:qradhi)) + qmo(i,j+1,k-1,QREINT)
#endif

             endif

          end do
       end do
    end do

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
                     hdt, cdtdy, cdtdz, lo, hi)


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


    integer, intent(in) :: qd_lo(3),qd_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fy_lo(3),fy_hi(3)
    integer, intent(in) :: fz_lo(3),fz_hi(3)
    integer, intent(in) :: qy_lo(3),qy_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

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

    integer i, j, k, n, nqp, ipassive

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

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                drr    = - cdtdy*(fyz(i,j+1,k-1,URHO) - fyz(i,j,k-1,URHO)) &
                         - cdtdz*(fzy(i,j  ,k,URHO) - fzy(i,j,k-1,URHO))
                dcompn = - cdtdy*(fyz(i,j+1,k-1,n   ) - fyz(i,j,k-1,n)) &
                         - cdtdz*(fzy(i,j  ,k,n   ) - fzy(i,j,k-1,n))

                if (i >= lo(1)+1) then
                   rrr = qp(i,j,k-1,QRHO)
                   compr = rrr*qp(i,j,k-1,nqp)

                   rrnewr = rrr +drr
                   compnr = compr +dcompn

                   qpo(i  ,j,k-1,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k,nqp)
                end if

                if (i <= hi(1)-1) then
                   rrl = qm(i+1,j,k-1,QRHO)
                   compl = rrl*qm(i+1,j,k-1,nqp)

                   rrnewl = rrl + drr
                   compnl = compl +dcompn

                   qmo(i+1,j,k-1,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k,nqp)
                end if
             end do
          end do
       end do
    end do

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! add the transverse yz and zy differences to the x-states for the
             ! fluid variables
             !-------------------------------------------------------------------

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
#endif

             pggyp  = qy(i,j+1,k-1,GDPRES)
             pggym  = qy(i,j  ,k-1,GDPRES)
             ugyp  = qy(i,j+1,k-1,GDV   )
             ugym  = qy(i,j  ,k-1,GDV   )
             gegyp = qy(i,j+1,k-1,GDGAME)
             gegym = qy(i,j  ,k-1,GDGAME)
#ifdef RADIATION
             ergyp = qy(i,j+1,k-1,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             pggzp  = qz(i,j,k,GDPRES)
             pggzm  = qz(i,j,k-1,GDPRES)
             ugzp  = qz(i,j,k,GDW   )
             ugzm  = qz(i,j,k-1,GDW   )
             gegzp = qz(i,j,k,GDGAME)
             gegzm = qz(i,j,k-1,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             duyp = pggyp*ugyp - pggym*ugym
             pyav = HALF*(pggyp+pggym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMC))*duy - uyav*dgey )
#endif

             duzp = pggzp*ugzp - pggzm*ugzm
             pzav = HALF*(pggzp+pggzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMC))*duz - uzav*dgez )
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

             if (i >= lo(1)+1) then
                ! Convert to conservation form
                rrr = qp(i,j,k-1,QRHO)
                rur = rrr*qp(i,j,k-1,QU)
                rvr = rrr*qp(i,j,k-1,QV)
                rwr = rrr*qp(i,j,k-1,QW)
                ekenr = HALF*rrr*(qp(i,j,k-1,QU)**2 + qp(i,j,k-1,QV)**2 + &
                     qp(i,j,k-1,QW)**2)
                rer = qp(i,j,k-1,QREINT) + ekenr
#ifdef RADIATION
                err = qp(i,j,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewr = rrr - cdtdy*(fyz(i,j+1,k-1,URHO) - fyz(i,j,k-1,URHO)) &
                             - cdtdz*(fzy(i,j,k,URHO) - fzy(i,j,k-1,URHO))
                runewr = rur - cdtdy*(fyz(i,j+1,k-1,UMX) - fyz(i,j,k-1,UMX)) &
                             - cdtdz*(fzy(i,j,k,UMX) - fzy(i,j,k-1,UMX))
                rvnewr = rvr - cdtdy*(fyz(i,j+1,k-1,UMY) - fyz(i,j,k-1,UMY)) &
                             - cdtdz*(fzy(i,j,k,UMY) - fzy(i,j,k-1,UMY))
                rwnewr = rwr - cdtdy*(fyz(i,j+1,k-1,UMZ) - fyz(i,j,k-1,UMZ)) &
                             - cdtdz*(fzy(i,j,k,UMZ) - fzy(i,j,k-1,UMZ))
                renewr = rer - cdtdy*(fyz(i,j+1,k-1,UEDEN) - fyz(i,j,k-1,UEDEN)) &
                             - cdtdz*(fzy(i,j,k,UEDEN) - fzy(i,j,k-1,UEDEN))
#ifdef RADIATION
                rvnewr = rvnewr + dmy
                rwnewr = rwnewr + dmz
                renewr = renewr + dre
                ernewr = err(:) - cdtdy*(rfyz(i,j+1,k-1,:) - rfyz(i,j,k-1,:)) &
                                - cdtdz*(rfzy(i,j  ,k,:) - rfzy(i,j,k-1,:)) &
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

                qpo(i,j,k-1,QRHO  ) = rrnewr
                qpo(i,j,k-1,QU    ) = runewr/rrnewr
                qpo(i,j,k-1,QV    ) = rvnewr/rrnewr
                qpo(i,j,k-1,QW    ) = rwnewr/rrnewr

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qpo(i,j,k-1,QRHO  ) = qpo(i,j,k-1,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qpo(i,j,k-1,QU:QW) = qpo(i,j,k-1,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                qpo(i,j,k-1,QREINT) = renewr - rhoekenr
                if (ppm_type == 0) then
                   qpo(i,j,k-1,QREINT) = qpo(i,j,k-1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qpo(i,j,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qpo(i,j,k-1,QREINT) = qp(i,j,k-1,QREINT) &
                           - cdtdy*(fyz(i,j+1,k-1,UEINT) - fyz(i,j,k-1,UEINT) + pyav*duy) &
                           - cdtdz*(fzy(i,j  ,k,UEINT) - fzy(i,j,k-1,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qpo(i,j,k-1,QREINT) = qpo(i,j,k-1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewr = qp(i,j,k-1,QPRES) - pynew - pznew
                      qpo(i,j,k-1,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,k-1,QPRES) = qpo(i,j,k-1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qpo(i,j,k-1,QGAME) = qp(i,j,k-1,QGAME) + geynew + geznew
                      
                      ! and compute the p edge state from this and (rho e)
                      qpo(i,j,k-1,QPRES) = qpo(i,j,k-1,QREINT)*(qpo(i,j,k-1,QGAME)-ONE)
                   end if
                else
                   qpo(i,j,k-1,QPRES) = qp(i,j,k-1,QPRES)
                   if (ppm_type == 0) then
                      qpo(i,j,k-1,QPRES) = qpo(i,j,k-1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qpo(i,j,k-1,QGAME) = qp(i,j,k-1,QGAME)
                endif

                qpo(i,j,k-1,QPRES) = max(qpo(i,j,k-1,QPRES), small_pres)

                call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, k-1)

#ifdef RADIATION
                qpo(i,j,k-1,qrad:qradhi) = ernewr(:)
                qpo(i,j,k-1,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k-1,QPRES)
                qpo(i,j,k-1,qreitot) = sum(qpo(i,j,k-1,qrad:qradhi)) + qpo(i,j,k-1,QREINT)
#endif

             endif

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             if (i <=hi(1)-1) then
                ! Convert to conservation form
                rrl = qm(i+1,j,k-1,QRHO)
                rul = rrl*qm(i+1,j,k-1,QU)
                rvl = rrl*qm(i+1,j,k-1,QV)
                rwl = rrl*qm(i+1,j,k-1,QW)
                ekenl = HALF*rrl*(qm(i+1,j,k-1,QU)**2 + qm(i+1,j,k-1,QV)**2 + &
                     qm(i+1,j,k-1,QW)**2)
                rel = qm(i+1,j,k-1,QREINT) + ekenl
#ifdef RADIATION
                erl = qm(i+1,j,k-1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewl = rrl - cdtdy*(fyz(i,j+1,k-1,URHO) - fyz(i,j,k-1,URHO)) &
                             - cdtdz*(fzy(i,j,k,URHO) - fzy(i,j,k-1,URHO))
                runewl = rul - cdtdy*(fyz(i,j+1,k-1,UMX) - fyz(i,j,k-1,UMX)) &
                             - cdtdz*(fzy(i,j,k,UMX) - fzy(i,j,k-1,UMX))
                rvnewl = rvl - cdtdy*(fyz(i,j+1,k-1,UMY) - fyz(i,j,k-1,UMY)) &
                             - cdtdz*(fzy(i,j,k,UMY) - fzy(i,j,k-1,UMY))
                rwnewl = rwl - cdtdy*(fyz(i,j+1,k-1,UMZ) - fyz(i,j,k-1,UMZ)) &
                             - cdtdz*(fzy(i,j,k,UMZ) - fzy(i,j,k-1,UMZ))
                renewl = rel - cdtdy*(fyz(i,j+1,k-1,UEDEN) - fyz(i,j,k-1,UEDEN)) &
                             - cdtdz*(fzy(i,j,k,UEDEN) - fzy(i,j,k-1,UEDEN))
#ifdef RADIATION
                rvnewl = rvnewl + dmy
                rwnewl = rwnewl + dmz
                renewl = renewl + dre
                ernewl = erl(:) - cdtdy*(rfyz(i,j+1,k-1,:) - rfyz(i,j,k-1,:)) &
                                - cdtdz*(rfzy(i,j  ,k,:) - rfzy(i,j,k-1,:)) &
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

                qmo(i+1,j,k-1,QRHO   ) = rrnewl
                qmo(i+1,j,k-1,QU     ) = runewl/rrnewl
                qmo(i+1,j,k-1,QV     ) = rvnewl/rrnewl
                qmo(i+1,j,k-1,QW     ) = rwnewl/rrnewl

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qmo(i+1,j,k-1,QRHO   ) = qmo(i+1,j,k-1,QRHO   ) + hdt*srcQ(i,j,k,QRHO)
                   qmo(i+1,j,k-1,QU:QW) = qmo(i+1,j,k-1,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
                qmo(i+1,j,k-1,QREINT ) = renewl - rhoekenl
                if (ppm_type == 0) then
                   qmo(i+1,j,k-1,QREINT ) = qmo(i+1,j,k-1,QREINT ) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qmo(i+1,j,k-1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qmo(i+1,j,k-1,QREINT ) = qm(i+1,j,k-1,QREINT) &
                           - cdtdy*(fyz(i,j+1,k-1,UEINT) - fyz(i,j,k-1,UEINT) + pyav*duy) &
                           - cdtdz*(fzy(i,j  ,k,UEINT) - fzy(i,j,k-1,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qmo(i+1,j,k-1,QREINT ) = qmo(i+1,j,k-1,QREINT ) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewl = qm(i+1,j,k-1,QPRES) - pynew - pznew
                      qmo(i+1,j,k-1,QPRES  ) = pnewl
                      if (ppm_type == 0) then
                         qmo(i+1,j,k-1,QPRES  ) = qmo(i+1,j,k-1,QPRES  ) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qmo(i+1,j,k-1,QGAME) = qm(i+1,j,k-1,QGAME) + geynew + geznew
                      
                      ! and compute the p edge state from this and (rho e)
                      qmo(i+1,j,k-1,QPRES) = qmo(i+1,j,k-1,QREINT)*(qmo(i+1,j,k-1,QGAME)-ONE)
                   end if
                else
                   qmo(i+1,j,k-1,QPRES  ) = qm(i+1,j,k-1,QPRES)
                   if (ppm_type == 0) then
                      qmo(i+1,j,k-1,QPRES  ) = qmo(i+1,j,k-1,QPRES  ) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qmo(i+1,j,k-1,QGAME) = qm(i+1,j,k-1,QGAME)
                endif

                qmo(i+1,j,k-1,QPRES) = max(qmo(i+1,j,k-1,QPRES), small_pres)

                call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i+1, j, k-1)

#ifdef RADIATION
                qmo(i+1,j,k-1,qrad:qradhi) = ernewl(:)
                qmo(i+1,j,k-1,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i+1,j,k-1,QPRES)
                qmo(i+1,j,k-1,qreitot) = sum(qmo(i+1,j,k-1,qrad:qradhi)) + qmo(i+1,j,k-1,QREINT)
#endif

             endif

          end do
       end do
    end do

  end subroutine transyz

end module ctu_advection_module
