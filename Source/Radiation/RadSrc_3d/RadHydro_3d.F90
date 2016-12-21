module rad_advection_module

  use bl_constants_module

  implicit none

  private

  public umeth3d_rad, consup_rad

contains

  ! ::: ---------------------------------------------------------------
  ! ::: :: UMETH3D     Compute hyperbolic fluxes using unsplit second
  ! ::: ::               order Godunov integrator.
  ! ::: ::
  ! ::: :: inputs/outputs
  ! ::: :: q           => (const)  input state, primitives
  ! ::: :: qaux        => (const)  auxillary hydro data
  ! ::: :: flatn       => (const)  flattening parameter
  ! ::: :: src         => (const)  source
  ! ::: :: nx          => (const)  number of cells in X direction
  ! ::: :: ny          => (const)  number of cells in Y direction
  ! ::: :: nz          => (const)  number of cells in Z direction
  ! ::: :: delta       => (const)  grid spacing in X, Y, Z direction
  ! ::: :: dt          => (const)  time stepsize
  ! ::: :: flux1      <=  (modify) flux in X direction on X edges
  ! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
  ! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
  ! ::: ----------------------------------------------------------------

  subroutine umeth3d_rad(q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         lam,lam_lo,lam_hi, &
                         flatn, &
                         srcQ, src_lo, src_hi, &
                         lo, hi, dx, dt, &
                         flux1, fd1_lo, fd1_hi, &
                         flux2, fd2_lo, fd2_hi, &
                         flux3, fd3_lo, fd3_hi, &
                         rflux1,rfd1_lo, rfd1_hi, &
                         rflux2,rfd2_lo, rfd2_hi, &
                         rflux3,rfd3_lo, rfd3_hi, &
                         q1, q1_lo, q1_hi, &
                         q2, q2_lo, q2_hi, &
                         q3, q3_lo, q3_hi, &                         
                         pdivu, domlo, domhi)

    use meth_params_module, only : QVAR, NVAR, QU, ppm_type, hybrid_riemann, &
                                   QRADVAR, NQ, &
                                   GDPRES, GDU, GDV, GDW, GDERADS, GDLAMS, &
                                   ppm_trace_sources, QC, QCG, QCSML, QGAMC, QGAMCG, NQAUX, NGDNV

    use trace_ppm_rad_module, only : tracexy_ppm_rad, tracez_ppm_rad
    use transverse_module
    use ppm_module
    use rad_params_module, only : ngroups
    use riemann_module, only : cmpflx, shock
    use mempool_module, only : bl_allocate, bl_deallocate
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    implicit none

    integer :: lam_lo(3), lam_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer :: rfd1_lo(3), rfd1_hi(3)
    integer :: rfd2_lo(3), rfd2_hi(3)
    integer :: rfd3_lo(3), rfd3_hi(3)
    integer :: qd_lo(3), qd_hi(3)
    integer :: qa_lo(3), qa_hi(3)
    integer :: q1_lo(3), q1_hi(3)
    integer :: q2_lo(3), q2_hi(3)
    integer :: q3_lo(3), q3_hi(3)
    integer :: lo(3), hi(3)
    integer :: domlo(3), domhi(3)

    integer :: fd1_lo(3), fd1_hi(3)
    integer :: fd2_lo(3), fd2_hi(3)
    integer :: fd3_lo(3), fd3_hi(3)

    double precision     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    double precision flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    double precision flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    double precision flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)
    double precision ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    double precision ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    double precision ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    double precision dx(3), dt

    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)

    double precision rflux1(rfd1_lo(1):rfd1_hi(1),rfd1_lo(2):rfd1_hi(2),rfd1_lo(3):rfd1_hi(3),0:ngroups-1)
    double precision rflux2(rfd2_lo(1):rfd2_hi(1),rfd2_lo(2):rfd2_hi(2),rfd2_lo(3):rfd2_hi(3),0:ngroups-1)
    double precision rflux3(rfd3_lo(1):rfd3_hi(1),rfd3_lo(2):rfd3_hi(2),rfd3_lo(3):rfd3_hi(3),0:ngroups-1)


    ! Local variables

    integer km,kc,kt,k3d,n
    integer i,j

    double precision dtdx, dtdy, dtdz, hdt
    double precision cdtdx, cdtdy, cdtdz
    double precision hdtdx, hdtdy, hdtdz

    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer :: qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
    double precision, pointer :: qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)

    double precision, pointer :: qmxy(:,:,:,:),qpxy(:,:,:,:)
    double precision, pointer :: qmxz(:,:,:,:),qpxz(:,:,:,:)

    double precision, pointer :: qmyx(:,:,:,:),qpyx(:,:,:,:)
    double precision, pointer :: qmyz(:,:,:,:),qpyz(:,:,:,:)

    double precision, pointer :: qmzx(:,:,:,:),qpzx(:,:,:,:)
    double precision, pointer :: qmzy(:,:,:,:),qpzy(:,:,:,:)

    double precision, pointer :: qxl(:,:,:,:),qxr(:,:,:,:)
    double precision, pointer :: qyl(:,:,:,:),qyr(:,:,:,:)
    double precision, pointer :: qzl(:,:,:,:),qzr(:,:,:,:)

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    double precision, pointer ::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)
    double precision, pointer ::  rfx(:,:,:,:), rfy(:,:,:,:),rfz(:,:,:,:)

    double precision, pointer :: fxy(:,:,:,:), fxz(:,:,:,:)
    double precision, pointer :: fyx(:,:,:,:), fyz(:,:,:,:)
    double precision, pointer :: fzx(:,:,:,:), fzy(:,:,:,:)
    double precision, pointer ::rfxy(:,:,:,:),rfxz(:,:,:,:)
    double precision, pointer ::rfyx(:,:,:,:),rfyz(:,:,:,:)
    double precision, pointer ::rfzx(:,:,:,:),rfzy(:,:,:,:)

    double precision, pointer :: qgdnvx(:,:,:,:), qgdnvxf(:,:,:,:), qgdnvtmpx(:,:,:,:)
    double precision, pointer :: qgdnvy(:,:,:,:), qgdnvyf(:,:,:,:), qgdnvtmpy(:,:,:,:)
    double precision, pointer :: qgdnvz(:,:,:,:), qgdnvtmpz1(:,:,:,:), qgdnvtmpz2(:,:,:,:), qgdnvzf(:,:,:,:)

    double precision, pointer :: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, pointer :: Ip_src(:,:,:,:,:,:), Im_src(:,:,:,:,:,:)

    double precision, pointer :: shk(:,:,:)

    integer :: qt_lo(3), qt_hi(3)
    integer :: It_lo(3), It_hi(3)
    integer :: shk_lo(3), shk_hi(3)
    integer :: fx_lo(3), fx_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: fz_lo(3), fz_hi(3)

    qt_lo = [lo(1)-1, lo(2)-1, 1]
    qt_hi = [hi(1)+2, hi(2)+2, 2]

    It_lo = [lo(1)-1, lo(2)-1, 1]
    It_hi = [hi(1)+1, hi(2)+1, 2]

    shk_lo(:) = lo(:) - 1
    shk_hi(:) = hi(:) + 1

    fx_lo = [lo(1)    , lo(2) - 1, 1]
    fx_hi = [hi(1) + 1, hi(2) + 1, 2]

    fy_lo = [lo(1) - 1, lo(2)    , 1]
    fy_hi = [hi(1) + 1, hi(2) + 1, 2]

    fz_lo = [lo(1) - 1, lo(2) - 1, 1]
    fz_hi = [hi(1) + 1, hi(2) + 1, 2]

    call bl_allocate (    qgdnvx, qt_lo, qt_hi, ngdnv)
    call bl_allocate (   qgdnvxf, qt_lo, qt_hi, ngdnv)
    call bl_allocate ( qgdnvtmpx, qt_lo, qt_hi, ngdnv)

    call bl_allocate (    qgdnvy, qt_lo, qt_hi, ngdnv)
    call bl_allocate (   qgdnvyf, qt_lo, qt_hi, ngdnv)
    call bl_allocate ( qgdnvtmpy, qt_lo, qt_hi, ngdnv)

    call bl_allocate (    qgdnvz, qt_lo, qt_hi, ngdnv)
    call bl_allocate (   qgdnvzf, qt_lo, qt_hi, ngdnv)
    call bl_allocate (qgdnvtmpz1, qt_lo, qt_hi, ngdnv)
    call bl_allocate (qgdnvtmpz2, qt_lo, qt_hi, ngdnv)

    call bl_allocate ( qxm, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qxp, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qmxy, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qpxy, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qmxz, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qpxz, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qym, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qyp, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qmyx, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qpyx, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qmyz, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qpyz, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qzm, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qzp, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qxl, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qxr, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qyl, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qyr, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qzl, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qzr, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qmzx, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qpzx, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( qmzy, qt_lo, qt_hi, QRADVAR)
    call bl_allocate ( qpzy, qt_lo, qt_hi, QRADVAR)

    call bl_allocate ( fx, fx_lo, fx_hi, NVAR)
    call bl_allocate (rfx, fx_lo(1), fx_hi(1), fx_lo(2), fx_hi(2), fx_lo(3), fx_hi(3), 0, ngroups-1)

    call bl_allocate ( fy, fy_lo, fy_hi, NVAR)
    call bl_allocate (rfy, fy_lo(1), fy_hi(1), fy_lo(2), fy_hi(2), fy_lo(3), fy_hi(3), 0, ngroups-1)

    call bl_allocate ( fz, fz_lo, fz_hi, NVAR)
    call bl_allocate (rfz, fz_lo(1), fz_hi(1), fz_lo(2), fz_hi(2), fz_lo(3), fz_hi(3), 0, ngroups-1)

    call bl_allocate ( fxy, fx_lo, fx_hi, NVAR)
    call bl_allocate (rfxy, fx_lo(1), fx_hi(1), fx_lo(2), fx_hi(2), fx_lo(3), fx_hi(3), 0, ngroups-1)
    call bl_allocate ( fxz, fx_lo, fx_hi, NVAR)
    call bl_allocate (rfxz, fx_lo(1), fx_hi(1), fx_lo(2), fx_hi(2), fx_lo(3), fx_hi(3), 0, ngroups-1)

    call bl_allocate ( fyx, fy_lo, fy_hi, NVAR)
    call bl_allocate (rfyx, fy_lo(1), fy_hi(1), fy_lo(2), fy_hi(2), fy_lo(3), fy_hi(3), 0, ngroups-1)
    call bl_allocate ( fyz, fy_lo, fy_hi, NVAR)
    call bl_allocate (rfyz, fy_lo(1), fy_hi(1), fy_lo(2), fy_hi(2), fy_lo(3), fy_hi(3), 0, ngroups-1)

    call bl_allocate ( fzx, fz_lo, fz_hi, NVAR)
    call bl_allocate (rfzx, fz_lo(1), fz_hi(1), fz_lo(2), fz_hi(2), fz_lo(3), fz_hi(3), 0, ngroups-1)
    call bl_allocate ( fzy, fz_lo, fz_hi, NVAR)
    call bl_allocate (rfzy, fz_lo(1), fz_hi(1), fz_lo(2), fz_hi(2), fz_lo(3), fz_hi(3), 0, ngroups-1)

    ! x-index, y-index, z-index, dim, characteristics, variables
    call bl_allocate ( Ip, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 1, 3, 1, 3, 1, QRADVAR)
    call bl_allocate ( Im, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 1, 3, 1, 3, 1, QRADVAR)

    call bl_allocate ( Ip_src, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 1, 3, 1, 3, 1, QRADVAR)
    call bl_allocate ( Im_src, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 1, 3, 1, 3, 1, QRADVAR)

    ! for the hybrid Riemann solver
    call bl_allocate(shk, shk_lo, shk_hi)

    ! Local constants
    dtdx = dt/dx(1)
    dtdy = dt/dx(2)
    dtdz = dt/dx(3)
    hdt = HALF*dt
    hdtdx = HALF*dtdx
    hdtdy = HALF*dtdy
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD

    ! Initialize pdivu to zero
    pdivu(:,:,:) = ZERO


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

    do k3d = lo(3)-1, hi(3)+1

       ! Swap pointers to levels
       kt = km
       km = kc
       kc = kt

       if (ppm_type .le. 0) then
          call bl_error("ppm_type <=0 is not supported in umeth3d_rad")
       else

          do n=1,QRADVAR
             call ppm(q(:,:,:,n),qd_lo,qd_hi, &
                      q(:,:,:,QU:),qaux(:,:,:,QC),qd_lo,qd_hi, &
                      flatn,qd_lo,qd_hi, &
                      Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), It_lo, It_hi, &
                      lo(1),lo(2),hi(1),hi(2), dx, dt,k3d,kc)
          end do

          if (ppm_trace_sources .eq. 1) then
             do n=1,QVAR
                call ppm(srcQ(:,:,:,n),src_lo,src_hi, &
                         q(:,:,:,QU:),qaux(:,:,:,QC),qd_lo,qd_hi, &
                         flatn,qd_lo,qd_hi, &
                         Ip_src(:,:,:,:,:,n),Im_src(:,:,:,:,:,n),It_lo,It_hi, &
                         lo(1),lo(2),hi(1),hi(2), dx, dt,k3d,kc)
             enddo
          endif

          ! Compute U_x and U_y at kc (k3d)
          call tracexy_ppm_rad(lam, lam_lo, lam_hi, &
                               q, qaux(:,:,:,QC), qaux(:,:,:,QCG), flatn, qd_lo, qd_hi, &
                               Ip, Im, Ip_src, Im_src, &
                               qxm, qxp, qym, qyp, qt_lo, qt_hi, &
                               qaux(:,:,:,QGAMC), qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                               lo(1),lo(2),hi(1),hi(2),dt,kc,k3d)

       end if

       ! Compute \tilde{F}^x at kc (k3d)
       call cmpflx(qxm, qxp, qt_lo, qt_hi, &
                   fx, fx_lo, fx_hi, &
                   qgdnvx, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfx, fx_lo, fx_hi, &
                   qaux(:,:,:,QGAMCG), &
                   qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym, qyp, qt_lo, qt_hi, &
                   fy, fy_lo, fy_hi, &
                   qgdnvy, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfy, fy_lo, fy_hi, &
                   qaux(:,:,:,QGAMCG), &
                   qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, kc, kc, k3d, domlo, domhi)

       ! Compute U'^y_x at kc (k3d)
       call transy1(&
                    lam, lam_lo, lam_hi, &
                    qxm, qmxy, qxp, qpxy, qt_lo, qt_hi, &
                    fy, &
                    rfy, &
                    fy_lo, fy_hi, &
                    qgdnvy, qt_lo, qt_hi, &
                    qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                    cdtdy, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, k3d)

       ! Compute U'^x_y at kc (k3d)
       call transx1(&
                    lam, lam_lo, lam_hi, &
                    qym, qmyx, qyp, qpyx, qt_lo, qt_hi, &
                    fx, &
                    rfx, &
                    fx_lo, fx_hi, &
                    qgdnvx, qt_lo, qt_hi, &
                    qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                    cdtdx, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, k3d)

       ! Compute F^{x|y} at kc (k3d)
       call cmpflx(qmxy, qpxy, qt_lo, qt_hi, &
                   fxy, fx_lo, fx_hi, &
                   qgdnvtmpx, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfxy, fx_lo, fx_hi, &
                   qaux(:,:,:,QGAMCG), &
                   qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   1, lo(1), hi(1)+1, lo(2), hi(2), kc, kc, k3d, domlo, domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx, qpyx, qt_lo, qt_hi, &
                   fyx, fy_lo, fy_hi, &
                   qgdnvtmpy, qt_lo, qt_hi, &
                   lam, lam_lo, lam_hi, &
                   rfyx, fy_lo, fy_hi, &
                   qaux(:,:,:,QGAMCG), &
                   qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                   shk, shk_lo, shk_hi, &
                   2, lo(1), hi(1), lo(2), hi(2)+1, kc, kc, k3d, domlo, domhi)

       if (k3d .ge. lo(3)) then

          ! Compute U_z at kc (k3d)
          call tracez_ppm_rad(lam, lam_lo, lam_hi, &
                              q, qaux(:,:,:,QC), qaux(:,:,:,QCG), flatn, qd_lo, qd_hi, &
                              Ip, Im, Ip_src, Im_src, &
                              qzm, qzp, qt_lo, qt_hi, &
                              qaux(:,:,:,QGAMC), qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                              lo(1), lo(2), hi(1), hi(2), dt, km, kc, k3d)

          ! Compute \tilde{F}^z at kc (k3d)
          call cmpflx(qzm, qzp, qt_lo, qt_hi, &
                      fz, fz_lo, fz_hi, &
                      qgdnvz, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rfz, fz_lo, fz_hi, &
                      qaux(:,:,:,QGAMCG), &
                      qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1)-1, hi(1)+1, lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

          ! Compute U'^y_z at kc (k3d)
          call transy2(&
                       lam, lam_lo, lam_hi, &
                       qzm, qmzy, qzp, qpzy, qt_lo, qt_hi, &
                       fy, &
                       rfy, &
                       fy_lo, fy_hi, &
                       qgdnvy, qt_lo, qt_hi, &
                       qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                       cdtdy, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, km, k3d)

          ! Compute U'^x_z at kc (k3d)
          call transx2(&
                       lam, lam_lo, lam_hi, &
                       qzm, qmzx, qzp, qpzx, qt_lo, qt_hi, &
                       fx, &
                       rfx, &
                       fx_lo, fx_hi, &
                       qgdnvx, qt_lo, qt_hi, &
                       qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                       cdtdx, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, km, k3d)

          ! Compute F^{z|x} at kc (k3d)
          call cmpflx(qmzx, qpzx, qt_lo, qt_hi, &
                      fzx, fz_lo, fz_hi, &
                      qgdnvtmpz1, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rfzx, fz_lo, fz_hi, &
                      qaux(:,:,:,QGAMCG), &
                      qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1), hi(1), lo(2)-1, hi(2)+1, kc, kc, k3d, domlo, domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy, qpzy, qt_lo, qt_hi, &
                      fzy, fz_lo, fz_hi, &
                      qgdnvtmpz2, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rfzy, fz_lo, fz_hi, &
                      qaux(:,:,:,QGAMCG), &
                      qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3, lo(1)-1, hi(1)+1, lo(2), hi(2), kc, kc, k3d, domlo, domhi)

          ! Compute U''_z at kc (k3d)
          call transxy(&
                       lam, lam_lo, lam_hi, &
                       qzm, qzl, qzp, qzr, qt_lo, qt_hi, &
                       fxy, &
                       rfxy, &
                       fx_lo, fx_hi, &
                       fyx, &
                       rfyx, &
                       fy_lo, fy_hi, &
                       qgdnvtmpx, qt_lo, qt_hi, &
                       qgdnvtmpy, qt_lo, qt_hi, &
                       qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                       srcQ, src_lo, src_hi, &
                       hdt, hdtdx, hdtdy, lo(1), hi(1), lo(2), hi(2), kc, km, k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl, qzr, qt_lo, qt_hi, &
                      flux3, fd3_lo, fd3_hi, &
                      qgdnvzf, qt_lo, qt_hi, &
                      lam, lam_lo, lam_hi, &
                      rflux3, rfd3_lo, rfd3_hi, &
                      qaux(:,:,:,QGAMCG), &
                      qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                      shk, shk_lo, shk_hi, &
                      3,lo(1),hi(1),lo(2),hi(2),kc,k3d,k3d,domlo,domhi)

          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                q3(i,j,k3d,:) = qgdnvzf(i,j,kc,:)
             end do
          end do

          if (k3d .ge. lo(3)+1 .and. k3d .le. hi(3)+1) then
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        HALF*(qgdnvzf(i,j,kc,GDPRES) + qgdnvzf(i,j,km,GDPRES)) * &
                             (qgdnvzf(i,j,kc,GDW) - qgdnvzf(i,j,km,GDW))/dx(3)
                end do
             end do
          end if

          if (k3d .gt. lo(3)) then

             ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
             call transz(&
                         lam, lam_lo, lam_hi, &
                         qxm, qmxz, qxp, qpxz, qym, qmyz, qyp, qpyz, qt_lo, qt_hi, &
                         fz, &
                         rfz, &
                         fz_lo, fz_hi, &
                         qgdnvz, qt_lo, qt_hi, &
                         qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                         cdtdz, lo(1)-1, hi(1)+1, lo(2)-1, hi(2)+1, km, kc, k3d)

             ! Compute F^{x|z} at km (k3d-1)
             call cmpflx(qmxz, qpxz, qt_lo, qt_hi, &
                         fxz, fx_lo, fx_hi, &
                         qgdnvx, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rfxz, fx_lo, fx_hi, &
                         qaux(:,:,:,QGAMCG), &
                         qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, km, km, k3d-1, domlo, domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz, qpyz, qt_lo, qt_hi, &
                         fyz, fy_lo, fy_hi, &
                         qgdnvy, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rfyz, fy_lo, fy_hi, &
                         qaux(:,:,:,QGAMCG), &
                         qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, km, km, k3d-1, domlo, domhi)

             ! Compute U''_x at km (k3d-1)
             call transyz(&
                          lam, lam_lo, lam_hi, &
                          qxm, qxl, qxp, qxr, qt_lo, qt_hi, &
                          fyz, &
                          rfyz, &
                          fy_lo, fy_hi, &
                          fzy, &
                          rfzy, &
                          fz_lo, fz_hi, &
                          qgdnvy, qt_lo, qt_hi, &
                          qgdnvtmpz2, qt_lo, qt_hi, &
                          qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                          srcQ, src_lo, src_hi, &
                          hdt, hdtdy, hdtdz, lo(1)-1, hi(1)+1, lo(2), hi(2), km, kc, k3d-1)

             ! Compute U''_y at km (k3d-1)
             call transxz(&
                          lam, lam_lo, lam_hi, &
                          qym, qyl, qyp, qyr, qt_lo, qt_hi, &
                          fxz, &
                          rfxz, &
                          fx_lo, fx_hi, &
                          fzx, &
                          rfzx, &
                          fz_lo, fz_hi, &
                          qgdnvx, qt_lo, qt_hi, &
                          qgdnvtmpz1, qt_lo, qt_hi, &
                          qaux(:,:,:,QGAMCG), qa_lo, qa_hi, &
                          srcQ, src_lo, src_hi, &
                          hdt, hdtdx, hdtdz, lo(1), hi(1), lo(2)-1, hi(2)+1, km, kc, k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl, qxr, qt_lo, qt_hi, &
                         flux1, fd1_lo, fd1_hi, &
                         qgdnvxf, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rflux1, rfd1_lo, rfd1_hi, &
                         qaux(:,:,:,QGAMCG), &
                         qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         1,lo(1),hi(1)+1,lo(2),hi(2),km,k3d-1,k3d-1,domlo,domhi)

             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+2
                   q1(i,j,k3d-1,:) = qgdnvxf(i,j,km,:)
                end do
             end do


             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl, qyr, qt_lo, qt_hi, &
                         flux2,  fd2_lo, fd2_hi, &
                         qgdnvyf, qt_lo, qt_hi, &
                         lam, lam_lo, lam_hi, &
                         rflux2, rfd2_lo, rfd2_hi, &
                         qaux(:,:,:,QGAMCG), &
                         qaux(:,:,:,QGAMC), qaux(:,:,:,QCSML), qaux(:,:,:,QC), qa_lo, qa_hi, &
                         shk, shk_lo, shk_hi, &
                         2,lo(1),hi(1),lo(2),hi(2)+1,km,k3d-1,k3d-1,domlo,domhi)

             do j=lo(2)-1,hi(2)+2
                do i=lo(1)-1,hi(1)+1
                   q2(i,j,k3d-1,:) = qgdnvyf(i,j,km,:)
                end do
             end do

             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        HALF*(qgdnvxf(i+1,j,km,GDPRES) + qgdnvxf(i,j,km,GDPRES)) * &
                             (qgdnvxf(i+1,j,km,GDU) - qgdnvxf(i,j,km,GDU))/dx(1) + &
                        HALF*(qgdnvyf(i,j+1,km,GDPRES) + qgdnvyf(i,j,km,GDPRES)) * &
                             (qgdnvyf(i,j+1,km,GDV) - qgdnvyf(i,j,km,GDV))/dx(2)
                end do
             end do

          end if
       end if
    enddo

    ! Deallocate arrays
    call bl_deallocate(qgdnvx)
    call bl_deallocate(qgdnvxf)
    call bl_deallocate(qgdnvtmpx)

    call bl_deallocate(qgdnvy)
    call bl_deallocate(qgdnvyf)
    call bl_deallocate(qgdnvtmpy)

    call bl_deallocate(qgdnvz)
    call bl_deallocate(qgdnvzf)
    call bl_deallocate(qgdnvtmpz1)
    call bl_deallocate(qgdnvtmpz2)

    call bl_deallocate(qxm)
    call bl_deallocate(qxp)

    call bl_deallocate(qmxy)
    call bl_deallocate(qpxy)

    call bl_deallocate(qmxz)
    call bl_deallocate(qpxz)
    
    call bl_deallocate(qym)
    call bl_deallocate(qyp)

    call bl_deallocate(qmyx)
    call bl_deallocate(qpyx)

    call bl_deallocate(qmyz)
    call bl_deallocate(qpyz)

    call bl_deallocate(qzm)
    call bl_deallocate(qzp)

    call bl_deallocate(qxl)
    call bl_deallocate(qxr)
    call bl_deallocate(qyl)
    call bl_deallocate(qyr)
    call bl_deallocate(qzl)
    call bl_deallocate(qzr)
    
    call bl_deallocate(qmzx)
    call bl_deallocate(qpzx)

    call bl_deallocate(qmzy)
    call bl_deallocate(qpzy)

    call bl_deallocate(fx)
    call bl_deallocate(fy)
    call bl_deallocate(fz)

    call bl_deallocate(rfx)
    call bl_deallocate(rfy)
    call bl_deallocate(rfz)

    call bl_deallocate(fxy)
    call bl_deallocate(fxz)

    call bl_deallocate(rfxy)
    call bl_deallocate(rfxz)

    call bl_deallocate(fyx)
    call bl_deallocate(fyz)

    call bl_deallocate(rfyx)
    call bl_deallocate(rfyz)

    call bl_deallocate(fzx)
    call bl_deallocate(fzy)

    call bl_deallocate(rfzx)
    call bl_deallocate(rfzy)

    call bl_deallocate(Ip)
    call bl_deallocate(Im)

    call bl_deallocate(shk)

  end subroutine umeth3d_rad

  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine consup_rad(uin, uin_lo, uin_hi, &
                        q,q_lo,q_hi, &
                        uout, uout_lo, uout_hi, &
                        update,updt_lo,updt_hi, &
                        Erin, Erin_lo, Erin_hi, &
                        Erout, Erout_lo, Erout_hi, &
                        flux1, flux1_lo, flux1_hi, &
                        flux2, flux2_lo, flux2_hi, &
                        flux3, flux3_lo, flux3_hi, &
                        radflux1, radflux1_lo, radflux1_hi, &
                        radflux2, radflux2_lo, radflux2_hi, &
                        radflux3, radflux3_lo, radflux3_hi, &
                        q1, q1_lo, q1_hi, &
                        q2, q2_lo, q2_hi, &
                        q3, q3_lo, q3_hi, &
                        area1, area1_lo, area1_hi, &
                        area2, area2_lo, area2_hi, &
                        area3, area3_lo, area3_hi, &
                        vol, vol_lo, vol_hi, &
                        div, pdivu, &
                        lo,hi, dx, dt, nstep_fsp,verbose)

    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, NGDNV, QVAR, NQ, &
                                   fspace_type, comoving, &
                                   GDPRES, GDU, GDV, GDW, GDLAMS, GDERADS, &
                                   track_grid_losses, limit_fluxes_on_small_dens
    use bl_constants_module, only : ZERO, FOURTH, ONE
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
    use advection_util_3d_module, only : normalize_species_fluxes
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens
    use castro_util_module, only : position, linear_to_angular_momentum
    use prob_params_module, only : domlo_level, domhi_level, center
    use amrinfo_module, only : amr_level
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    implicit none

    integer nstep_fsp
    integer :: lo(3), hi(3)
    integer :: uin_lo(3), uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer :: uout_lo(3), uout_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer :: Erout_lo(3), Erout_hi(3)
    integer :: Erin_lo(3), Erin_hi(3)
    integer :: flux1_lo(3), flux1_hi(3)
    integer :: flux2_lo(3), flux2_hi(3)
    integer :: flux3_lo(3), flux3_hi(3)
    integer :: radflux1_lo(3), radflux1_hi(3)
    integer :: radflux2_lo(3), radflux2_hi(3)
    integer :: radflux3_lo(3), radflux3_hi(3)
    integer :: q1_lo(3), q1_hi(3)
    integer :: q2_lo(3), q2_hi(3)
    integer :: q3_lo(3), q3_hi(3)
    integer :: area1_lo(3), area1_hi(3)
    integer :: area2_lo(3), area2_hi(3)
    integer :: area3_lo(3), area3_hi(3)
    integer :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: verbose

    double precision uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    double precision, intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    double precision uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    double precision, intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)
    double precision  Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    double precision Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),Erout_lo(3):Erout_hi(3),0:ngroups-1)
    double precision flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    double precision radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),radflux1_lo(3):radflux1_hi(3),0:ngroups-1)
    double precision radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),radflux2_lo(3):radflux2_hi(3),0:ngroups-1)
    double precision radflux3(radflux3_lo(1):radflux3_hi(1),radflux3_lo(2):radflux3_hi(2),radflux3_lo(3):radflux3_hi(3),0:ngroups-1)
    double precision ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    double precision ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    double precision ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    double precision area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    double precision area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    double precision area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    double precision vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))

    double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    double precision dx(3), dt

    ! Local variables

    double precision :: div1
    integer          :: i, j, k, n, g

    double precision, dimension(0:ngroups-1) :: Erscale
    double precision, dimension(0:ngroups-1) :: ustar, af
    double precision :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp
    double precision :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    double precision :: Gf1E(3)
    double precision :: ux, uy, uz, divu, lamc, Egdc
    double precision :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    double precision :: dprdx, dprdy, dprdz, ek1, ek2, dek

    double precision :: mass_added_flux, E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux
    double precision :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    double precision :: eden_lost, xang_lost, yang_lost, zang_lost

    integer          :: domlo(3), domhi(3)
    double precision :: loc(3), ang_mom(3)

    double precision :: urho_new
    double precision :: umx_new1, umy_new1, umz_new1
    double precision :: umx_new2, umy_new2, umz_new2

    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if

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
                   div1 = difmag*min(ZERO, div1)

                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx(1)*div1*(uin(i,j,k,n)-uin(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                   div1 = difmag*min(ZERO, div1)

                   flux2(i,j,k,n) = flux2(i,j,k,n) + dx(2)*div1*(uin(i,j,k,n)-uin(i,j-1,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag*min(ZERO, div1)

                   flux3(i,j,k,n) = flux3(i,j,k,n) + dx(3)*div1*(uin(i,j,k,n)-uin(i,j,k-1,n))
                enddo
             enddo
          enddo

       endif

    enddo

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
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                       flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                       flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)

                ! Add the p div(u) source term to (rho e)
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
                lamc = (q1(i,j,k,GDLAMS+g) + q1(i+1,j,k,GDLAMS+g) + &
                        q2(i,j,k,GDLAMS+g) + q2(i,j+1,k,GDLAMS+g) + &
                        q3(i,j,k,GDLAMS+g) + q3(i,j,k+1,GDLAMS+g) ) / 6.d0
                dprdx = dprdx + lamc*(q1(i+1,j,k,GDERADS+g) - q1(i,j,k,GDERADS+g))/dx(1)
                dprdy = dprdy + lamc*(q2(i,j+1,k,GDERADS+g) - q2(i,j,k,GDERADS+g))/dx(2)
                dprdz = dprdz + lamc*(q3(i,j,k+1,GDERADS+g) - q3(i,j,k,GDERADS+g))/dx(3)
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

                ux = HALF*(q1(i,j,k,GDU) + q1(i+1,j,k,GDU))
                uy = HALF*(q2(i,j,k,GDV) + q2(i,j+1,k,GDV))
                uz = HALF*(q3(i,j,k,GDW) + q3(i,j,k+1,GDW))

                dudx(1) = (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1)
                dudx(2) = (q1(i+1,j,k,GDV) - q1(i,j,k,GDV))/dx(1)
                dudx(3) = (q1(i+1,j,k,GDW) - q1(i,j,k,GDW))/dx(1)

                dudy(1) = (q2(i,j+1,k,GDU) - q2(i,j,k,GDU))/dx(2)
                dudy(2) = (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2)
                dudy(3) = (q2(i,j+1,k,GDW) - q2(i,j,k,GDW))/dx(2)

                dudz(1) = (q3(i,j,k+1,GDU) - q3(i,j,k,GDU))/dx(3)
                dudz(2) = (q3(i,j,k+1,GDV) - q3(i,j,k,GDV))/dx(3)
                dudz(3) = (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)

                divu = dudx(1) + dudy(2) + dudz(3)

                ! Note that for single group, fspace_type is always 1
                do g=0, ngroups-1

                   nhat(1) = (q1(i+1,j,k,GDERADS+g) - q1(i,j,k,GDERADS+g))/dx(1)
                   nhat(2) = (q2(i,j+1,k,GDERADS+g) - q2(i,j,k,GDERADS+g))/dx(2)
                   nhat(3) = (q3(i,j,k+1,GDERADS+g) - q3(i,j,k,GDERADS+g))/dx(3)

                   GnDotu(1) = dot_product(nhat, dudx)
                   GnDotu(2) = dot_product(nhat, dudy)
                   GnDotu(3) = dot_product(nhat, dudz)

                   nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.d-50)

                   lamc = (q1(i,j,k,GDLAMS+g) + q1(i+1,j,k,GDLAMS+g) + &
                           q2(i,j,k,GDLAMS+g) + q2(i,j+1,k,GDLAMS+g) + &
                           q3(i,j,k,GDLAMS+g) + q3(i,j,k+1,GDLAMS+g) ) / 6.d0
                   Eddf = Edd_factor(lamc)
                   f1 = (ONE-Eddf)*HALF
                   f2 = (3.d0*Eddf-ONE)*HALF
                   af(g) = -(f1*divu + f2*nnColonDotGu)

                   if (fspace_type .eq. 1) then
                      Eddfxp = Edd_factor(q1(i+1,j  ,k  ,GDLAMS+g))
                      Eddfxm = Edd_factor(q1(i  ,j  ,k  ,GDLAMS+g))
                      Eddfyp = Edd_factor(q2(i  ,j+1,k  ,GDLAMS+g))
                      Eddfym = Edd_factor(q2(i  ,j  ,k  ,GDLAMS+g))
                      Eddfzp = Edd_factor(q3(i  ,j  ,k+1,GDLAMS+g))
                      Eddfzm = Edd_factor(q3(i  ,j  ,k  ,GDLAMS+g))

                      f1xp = HALF*(ONE-Eddfxp)
                      f1xm = HALF*(ONE-Eddfxm)
                      f1yp = HALF*(ONE-Eddfyp)
                      f1ym = HALF*(ONE-Eddfym)
                      f1zp = HALF*(ONE-Eddfzp)
                      f1zm = HALF*(ONE-Eddfzm)

                      Gf1E(1) = (f1xp*q1(i+1,j,k,GDERADS+g) - f1xm*q1(i,j,k,GDERADS+g)) / dx(1)
                      Gf1E(2) = (f1yp*q2(i,j+1,k,GDERADS+g) - f1ym*q2(i,j,k,GDERADS+g)) / dx(2)
                      Gf1E(3) = (f1zp*q3(i,j,k+1,GDERADS+g) - f1zm*q3(i,j,k,GDERADS+g)) / dx(3)

                      Egdc = (q1(i,j,k,GDERADS+g) + q1(i+1,j,k,GDERADS+g) &
                           +  q2(i,j,k,GDERADS+g) + q2(i,j+1,k,GDERADS+g) &
                           +  q3(i,j,k,GDERADS+g) + q3(i,j,k+1,GDERADS+g) ) / 6.d0

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

    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (verbose .eq. 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mass_added_flux = mass_added_flux + ( flux1(i,j,k,URHO) - flux1(i+1,j,k,URHO) + &
                                                      flux2(i,j,k,URHO) - flux2(i,j+1,k,URHO) + &
                                                      flux3(i,j,k,URHO) - flux3(i,j,k+1,URHO) )
                xmom_added_flux = xmom_added_flux + ( flux1(i,j,k,UMX) - flux1(i+1,j,k,UMX) + &
                                                      flux2(i,j,k,UMX) - flux2(i,j+1,k,UMX) + &
                                                      flux3(i,j,k,UMX) - flux3(i,j,k+1,UMX) )
                ymom_added_flux = ymom_added_flux + ( flux1(i,j,k,UMY) - flux1(i+1,j,k,UMY) + &
                                                      flux2(i,j,k,UMY) - flux2(i,j+1,k,UMY) + &
                                                      flux3(i,j,k,UMY) - flux3(i,j,k+1,UMY) )
                zmom_added_flux = zmom_added_flux + ( flux1(i,j,k,UMZ) - flux1(i+1,j,k,UMZ) + &
                                                      flux2(i,j,k,UMZ) - flux2(i,j+1,k,UMZ) + &
                                                      flux3(i,j,k,UMZ) - flux3(i,j,k+1,UMZ) )
                E_added_flux    = E_added_flux    + ( flux1(i,j,k,UEDEN) - flux1(i+1,j,k,UEDEN) + &
                                                      flux2(i,j,k,UEDEN) - flux2(i,j+1,k,UEDEN) + &
                                                      flux3(i,j,k,UEDEN) - flux3(i,j,k+1,UEDEN) )
             enddo
          enddo
       enddo

    endif

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


  end subroutine consup_rad

end module rad_advection_module
