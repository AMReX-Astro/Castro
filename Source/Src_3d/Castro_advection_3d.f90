module advection_module

  implicit none

  private

  public umeth3d, ctoprim, consup
  
contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH3D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  cound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: nz          => (const)  number of cells in Z direction
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dz          => (const)  grid spacing in Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
! ::: ----------------------------------------------------------------

  subroutine umeth3d(q, c, gamc, csml, flatn, qd_lo, qd_hi, &
                     srcQ, src_lo, src_hi, &
                     lo, hi, dx, dt, &
                     flux1, fd1_lo, fd1_hi, &
                     flux2, fd2_lo, fd2_hi, &
                     flux3, fd3_lo, fd3_hi, &
                     q1, q1_lo, q1_hi, &
                     q2, q2_lo, q2_hi, &
                     q3, q3_lo, q3_hi, &
                     pdivu, domlo, domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QV, QW, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   ppm_type, &
                                   use_pslope, ppm_trace_sources, ppm_temp_fix, &
                                   hybrid_riemann
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
    use trace_module, only : tracexy, tracez
    use transverse_module
    use ppm_module, only : ppm
    use slope_module, only : uslope, pslope
    use network
    use eos_module
    use riemann_module, only: cmpflx, shock
    use bl_constants_module

    implicit none

    integer :: qd_lo(3), qd_hi(3)
    integer :: src_lo(3), src_hi(3)
    integer :: lo(3), hi(3)
    integer :: fd1_lo(3), fd1_hi(3)
    integer :: fd2_lo(3), fd2_hi(3)
    integer :: fd3_lo(3), fd3_hi(3)
    integer :: q1_lo(3), q1_hi(3)
    integer :: q2_lo(3), q2_hi(3)
    integer :: q3_lo(3), q3_hi(3)
    integer :: domlo(3), domhi(3)

    double precision ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    double precision ::     c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision ::  gamc(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision ::  csml(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    double precision :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    double precision :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)
    double precision ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    double precision ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    double precision ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    double precision :: pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: dx(3), dt
    double precision :: dxinv, dyinv, dzinv
    double precision :: dtdx, dtdy, dtdz, hdt
    double precision :: cdtdx, cdtdy, cdtdz
    double precision :: hdtdx, hdtdy, hdtdz

    integer :: km, kc, kt, k3d, n
    integer :: i, j, iwave, idim
        
    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer :: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
    double precision, pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
    double precision, pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)
    
    double precision, pointer :: qmxy(:,:,:,:), qpxy(:,:,:,:)
    double precision, pointer :: qmxz(:,:,:,:), qpxz(:,:,:,:)
    
    double precision, pointer :: qmyx(:,:,:,:), qpyx(:,:,:,:)
    double precision, pointer :: qmyz(:,:,:,:), qpyz(:,:,:,:)
    
    double precision, pointer :: qmzx(:,:,:,:), qpzx(:,:,:,:)
    double precision, pointer :: qmzy(:,:,:,:), qpzy(:,:,:,:)
    
    double precision, pointer :: qxl(:,:,:,:), qxr(:,:,:,:)
    double precision, pointer :: qyl(:,:,:,:), qyr(:,:,:,:)
    double precision, pointer :: qzl(:,:,:,:), qzr(:,:,:,:)
    
    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    double precision, pointer ::  fx(:,:,:,:), fy(:,:,:,:), fz(:,:,:,:)
    
    double precision, pointer :: fxy(:,:,:,:), fxz(:,:,:,:)
    double precision, pointer :: fyx(:,:,:,:), fyz(:,:,:,:)
    double precision, pointer :: fzx(:,:,:,:), fzy(:,:,:,:)
    
    double precision, pointer :: qgdnvx(:,:,:,:), qgdnvxf(:,:,:,:), qgdnvtmpx(:,:,:,:)    
    double precision, pointer :: qgdnvy(:,:,:,:), qgdnvyf(:,:,:,:), qgdnvtmpy(:,:,:,:)    
    double precision, pointer :: qgdnvz(:,:,:,:), qgdnvzf(:,:,:,:), qgdnvtmpz1(:,:,:,:), qgdnvtmpz2(:,:,:,:)
    
    double precision, pointer :: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, pointer :: Ip_src(:,:,:,:,:,:), Im_src(:,:,:,:,:,:)
    double precision, pointer :: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    double precision, pointer :: shk(:,:,:)
    
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
    
    call bl_allocate ( qxm, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qxp, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qmxy, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qpxy, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qmxz, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qpxz, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qym, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qyp, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qmyx, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qpyx, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qmyz, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qpyz, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qzm, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qzp, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qxl, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qxr, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qyl, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qyr, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qzl, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qzr, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qmzx, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qpzx, qt_lo, qt_hi, QVAR)

    call bl_allocate ( qmzy, qt_lo, qt_hi, QVAR)
    call bl_allocate ( qpzy, qt_lo, qt_hi, QVAR)

    call bl_allocate ( fx, fx_lo, fx_hi, NVAR)
    call bl_allocate ( fy, fy_lo, fy_hi, NVAR)
    call bl_allocate ( fz, fz_lo, fz_hi, NVAR)

    call bl_allocate ( fxy, fx_lo, fx_hi, NVAR)
    call bl_allocate ( fxz, fx_lo, fx_hi, NVAR)

    call bl_allocate ( fyx, fy_lo, fy_hi, NVAR)
    call bl_allocate ( fyz, fy_lo, fy_hi, NVAR)

    call bl_allocate ( fzx, fz_lo, fz_hi, NVAR)
    call bl_allocate ( fzy, fz_lo, fz_hi, NVAR)

    if (ppm_type .gt. 0) then
       ! x-index, y-index, z-index, dim, characteristics, variables
       call bl_allocate ( Ip, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,QVAR)
       call bl_allocate ( Im, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,QVAR)
       
       ! for source terms
       call bl_allocate ( Ip_src, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,QVAR)
       call bl_allocate ( Im_src, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,QVAR)
       
       ! for gamc -- needed for the reference state in eigenvectors
       call bl_allocate ( Ip_gc, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,1)
       call bl_allocate ( Im_gc, It_lo(1),It_hi(1),It_lo(2),It_hi(2),It_lo(3),It_hi(3),1,3,1,3,1,1)
    else
       call bl_allocate ( dqx, qt_lo, qt_hi, QVAR)
       call bl_allocate ( dqy, qt_lo, qt_hi, QVAR)
       call bl_allocate ( dqz, qt_lo, qt_hi, QVAR)
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

    ! Initialize pdivu to zero
    pdivu(:,:,:) = ZERO


    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q,qd_lo,qd_hi,shk,shk_lo,shk_hi,lo,hi,dx)
    else
       shk(:,:,:) = ZERO
    endif
    

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

       if (ppm_type .gt. 0) then

          do n=1,QVAR
             call ppm(q(:,:,:,n  ),  qd_lo,qd_hi, &
                      q(:,:,:,QU:QW),c,qd_lo,qd_hi, &
                      flatn,qd_lo,qd_hi, &
                      Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n),It_lo,It_hi, &
                      lo(1),lo(2),hi(1),hi(2),dx,dt,k3d,kc)
          end do

          if (ppm_trace_sources .eq. 1) then
             do n=1,QVAR
                call ppm(srcQ(:,:,:,n),src_lo,src_hi, &
                         q(:,:,:,QU:QW),c,qd_lo,qd_hi, &
                         flatn,qd_lo,qd_hi, &
                         Ip_src(:,:,:,:,:,n),Im_src(:,:,:,:,:,n),It_lo,It_hi, &
                         lo(1),lo(2),hi(1),hi(2),dx,dt,k3d,kc)
             enddo
          endif

          if (ppm_temp_fix /= 1) then
             call ppm(gamc(:,:,:),qd_lo,qd_hi, &
                      q(:,:,:,QU:QW),c,qd_lo,qd_hi, &
                      flatn,qd_lo,qd_hi, &
                      Ip_gc(:,:,:,:,:,1),Im_gc(:,:,:,:,:,1),It_lo,It_hi, &
                      lo(1),lo(2),hi(1),hi(2),dx,dt,k3d,kc)
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
          call tracexy_ppm(q,c,flatn,qd_lo,qd_hi, &
                           Ip,Im,Ip_src,Im_src,Ip_gc,Im_gc,It_lo,It_hi, &
                           qxm,qxp,qym,qyp,qt_lo,qt_hi, &
                           gamc,qd_lo,qd_hi, &
                           lo(1),lo(2),hi(1),hi(2),dt,kc,k3d)

       else

          ! Compute all slopes at kc (k3d)
          call uslope(q,flatn,qd_lo,qd_hi, &
                      dqx,dqy,dqz,qt_lo,qt_hi, &
                      lo(1),lo(2),hi(1),hi(2),kc,k3d,QVAR)
          
          if (use_pslope .eq. 1) &
               call pslope(q(:,:,:,QPRES),q(:,:,:,QRHO), &
                           flatn,qd_lo,qd_hi, &
                           dqx(:,:,:,QPRES),dqy(:,:,:,QPRES),dqz(:,:,:,QPRES), &
                           qt_lo,qt_hi, &
                           srcQ,src_lo,src_hi, &
                           lo(1),lo(2),hi(1),hi(2),kc,k3d,dx)

          ! Compute U_x and U_y at kc (k3d)
          call tracexy(q,c,qd_lo,qd_hi, &
                       dqx,dqy,qt_lo,qt_hi, &
                       qxm,qxp,qym,qyp,qt_lo,qt_hi, &
                       lo(1),lo(2),hi(1),hi(2),dx,dt,kc,k3d)
          
       end if

       ! Compute \tilde{F}^x at kc (k3d)
       call cmpflx(qxm,qxp,qt_lo,qt_hi, &
                   fx,fx_lo,fx_hi, &
                   qgdnvx,qt_lo,qt_hi, &
                   gamc,csml,c,qd_lo,qd_hi, &
                   shk,shk_lo,shk_hi, &
                   1,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,kc,kc,k3d,domlo,domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym,qyp,qt_lo,qt_hi, &
                   fy,fy_lo,fy_hi, &
                   qgdnvy,qt_lo,qt_hi, &
                   gamc,csml,c,qd_lo,qd_hi, &
                   shk,shk_lo,shk_hi, &
                   2,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,kc,kc,k3d,domlo,domhi)
       
       ! Compute U'^y_x at kc (k3d)
       call transy1(qxm,qmxy,qxp,qpxy,qt_lo,qt_hi, &
                    fy,fy_lo,fy_hi, &
                    qgdnvy,qt_lo,qt_hi, &
                    gamc,qd_lo,qd_hi, &
                    cdtdy,lo(1)-1,hi(1)+1,lo(2),hi(2),kc,k3d)

       ! Compute U'^x_y at kc (k3d)
       call transx1(qym,qmyx,qyp,qpyx,qt_lo,qt_hi, &
                    fx,fx_lo,fx_hi, &
                    qgdnvx,qt_lo,qt_hi, &
                    gamc,qd_lo,qd_hi, &
                    cdtdx,lo(1),hi(1),lo(2)-1,hi(2)+1,kc,k3d)

       ! Compute F^{x|y} at kc (k3d)
       call cmpflx(qmxy,qpxy,qt_lo,qt_hi, &
                   fxy,fx_lo,fx_hi, &
                   qgdnvtmpx,qt_lo,qt_hi, &
                   gamc,csml,c,qd_lo,qd_hi, &
                   shk,shk_lo,shk_hi, &
                   1,lo(1),hi(1)+1,lo(2),hi(2),kc,kc,k3d,domlo,domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx,qpyx,qt_lo,qt_hi, &
                   fyx,fy_lo,fy_hi, &
                   qgdnvtmpy,qt_lo,qt_hi, &
                   gamc,csml,c,qd_lo,qd_hi, &
                   shk,shk_lo,shk_hi, &
                   2,lo(1),hi(1),lo(2),hi(2)+1,kc,kc,k3d,domlo,domhi)

       if (k3d.ge.lo(3)) then
          
          ! Compute U_z at kc (k3d)
          if (ppm_type .gt. 0) then
             call tracez_ppm(q,c,flatn,qd_lo,qd_hi, &
                             Ip,Im,Ip_src,Im_src,Ip_gc,Im_gc,It_lo,It_hi, &
                             qzm,qzp,qt_lo,qt_hi, &
                             gamc,qd_lo,qd_hi, &
                             lo(1),lo(2),hi(1),hi(2),dt,km,kc,k3d)
          else
             call tracez(q,c,qd_lo,qd_hi, &
                         dqz,qt_lo,qt_hi, &
                         qzm,qzp,qt_lo,qt_hi, &
                         lo(1),lo(2),hi(1),hi(2),dx,dt,km,kc,k3d)
          end if

          ! Compute \tilde{F}^z at kc (k3d)
          call cmpflx(qzm,qzp,qt_lo,qt_hi, &
                      fz,fz_lo,fz_hi, &
                      qgdnvz,qt_lo,qt_hi, &
                      gamc,csml,c,qd_lo,qd_hi, &
                      shk,shk_lo,shk_hi, &
                      3,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,kc,kc,k3d,domlo,domhi)

          ! Compute U'^y_z at kc (k3d)
          call transy2(qzm,qmzy,qzp,qpzy,qt_lo,qt_hi, &
                       fy,fy_lo,fy_hi, &
                       qgdnvy,qt_lo,qt_hi, &
                       gamc,qd_lo,qd_hi, &
                       cdtdy,lo(1)-1,hi(1)+1,lo(2),hi(2),kc,km,k3d)

          ! Compute U'^x_z at kc (k3d)
          call transx2(qzm,qmzx,qzp,qpzx,qt_lo,qt_hi, &
                       fx,fx_lo,fx_hi, &
                       qgdnvx,qt_lo,qt_hi, &
                       gamc,qd_lo,qd_hi, &
                       cdtdx,lo(1),hi(1),lo(2)-1,hi(2)+1,kc,km,k3d)

          ! Compute F^{z|x} at kc (k3d)
          call cmpflx(qmzx,qpzx,qt_lo,qt_hi, &
                      fzx,fz_lo,fz_hi, &
                      qgdnvtmpz1,qt_lo,qt_hi, &
                      gamc,csml,c,qd_lo,qd_hi, &
                      shk,shk_lo,shk_hi, &
                      3,lo(1),hi(1),lo(2)-1,hi(2)+1,kc,kc,k3d,domlo,domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy,qpzy,qt_lo,qt_hi, &
                      fzy,fz_lo,fz_hi, &
                      qgdnvtmpz2,qt_lo,qt_hi, &
                      gamc,csml,c,qd_lo,qd_hi, &
                      shk,shk_lo,shk_hi, &                       
                      3,lo(1)-1,hi(1)+1,lo(2),hi(2),kc,kc,k3d,domlo,domhi)
          
          ! Compute U''_z at kc (k3d)
          call transxy(qzm,qzl,qzp,qzr,qt_lo,qt_hi, &
                       fxy,fx_lo,fx_hi, &
                       fyx,fy_lo,fy_hi, &
                       qgdnvtmpx,qt_lo,qt_hi, &
                       qgdnvtmpy,qt_lo,qt_hi, &
                       gamc,qd_lo,qd_hi, &
                       srcQ,src_lo,src_hi,&
                       hdt,hdtdx,hdtdy,lo(1),hi(1),lo(2),hi(2),kc,km,k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl,qzr,qt_lo,qt_hi, &
                      flux3,fd3_lo,fd3_hi, &
                      qgdnvzf,qt_lo,qt_hi, &
                      gamc,csml,c,qd_lo,qd_hi, &
                      shk,shk_lo,shk_hi, &
                      3,lo(1),hi(1),lo(2),hi(2),kc,k3d,k3d,domlo,domhi)

          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                q3(i,j,k3d,:) = qgdnvzf(i,j,kc,:)
             end do
          end do

          if (k3d .ge. lo(3)+1 .and. k3d .le. hi(3)+1) then
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        HALF*(qgdnvzf(i,j,kc,GDPRES) + qgdnvzf(i,j,km,GDPRES)) * &
                             (qgdnvzf(i,j,kc,GDW) - qgdnvzf(i,j,km,GDW))*dzinv
                end do
             end do
          end if
          
          if (k3d.gt.lo(3)) then

             ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
             call transz(qxm,qmxz,qxp,qpxz,qym,qmyz,qyp,qpyz,qt_lo,qt_hi, &
                         fz,fz_lo,fz_hi, &
                         qgdnvz,qt_lo,qt_hi, &
                         gamc,qd_lo,qd_hi, &
                         cdtdz,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,km,kc,k3d)
         
             ! Compute F^{x|z} at km (k3d-1)
             call cmpflx(qmxz,qpxz,qt_lo,qt_hi, &
                         fxz,fx_lo,fx_hi, &
                         qgdnvx,qt_lo,qt_hi, &
                         gamc,csml,c,qd_lo,qd_hi, &
                         shk,shk_lo,shk_hi, &
                         1,lo(1),hi(1)+1,lo(2)-1,hi(2)+1,km,km,k3d-1,domlo,domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz,qpyz,qt_lo,qt_hi, &
                         fyz,fy_lo,fy_hi, &
                         qgdnvy,qt_lo,qt_hi, &
                         gamc,csml,c,qd_lo,qd_hi, &
                         shk,shk_lo,shk_hi, &
                         2,lo(1)-1,hi(1)+1,lo(2),hi(2)+1,km,km,k3d-1,domlo,domhi)

             ! Compute U''_x at km (k3d-1)
             call transyz(qxm,qxl,qxp,qxr,qt_lo,qt_hi, &
                          fyz,fy_lo,fy_hi, &
                          fzy,fz_lo,fz_hi, &
                          qgdnvy,qt_lo,qt_hi, &
                          qgdnvtmpz2,qt_lo,qt_hi, &
                          gamc,qd_lo,qd_hi, &
                          srcQ,src_lo,src_hi, &
                          hdt,hdtdy,hdtdz,lo(1)-1,hi(1)+1,lo(2),hi(2),km,kc,k3d-1)

             ! Compute U''_y at km (k3d-1)
             call transxz(qym,qyl,qyp,qyr,qt_lo,qt_hi, &
                          fxz,fx_lo,fx_hi, &
                          fzx,fz_lo,fz_hi, &
                          qgdnvx,qt_lo,qt_hi, &
                          qgdnvtmpz1,qt_lo,qt_hi, &
                          gamc,qd_lo,qd_hi, &
                          srcQ,src_lo,src_hi, &
                          hdt,hdtdx,hdtdz,lo(1),hi(1),lo(2)-1,hi(2)+1,km,kc,k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl,qxr,qt_lo,qt_hi, &
                         flux1,fd1_lo,fd1_hi, &
                         qgdnvxf,qt_lo,qt_hi, &
                         gamc,csml,c,qd_lo,qd_hi, &
                         shk,shk_lo,shk_hi, &
                         1,lo(1),hi(1)+1,lo(2),hi(2),km,k3d-1,k3d-1,domlo,domhi)
             
             do j=lo(2)-1,hi(2)+1
                do i=lo(1)-1,hi(1)+2
                   q1(i,j,k3d-1,:) = qgdnvxf(i,j,km,:)
                end do
             end do
             
             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl,qyr,qt_lo,qt_hi, &
                         flux2,fd2_lo,fd2_hi, &
                         qgdnvyf,qt_lo,qt_hi, &
                         gamc,csml,c,qd_lo,qd_hi, &
                         shk,shk_lo,shk_hi, &
                         2,lo(1),hi(1),lo(2),hi(2)+1,km,k3d-1,k3d-1,domlo,domhi)

             do j=lo(2)-1,hi(2)+2
                do i=lo(1)-1,hi(1)+1
                   q2(i,j,k3d-1,:) = qgdnvyf(i,j,km,:)
                end do
             end do

             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        HALF*(qgdnvxf(i+1,j,km,GDPRES) + qgdnvxf(i,j,km,GDPRES)) *  &
                             (qgdnvxf(i+1,j,km,GDU) - qgdnvxf(i,j,km,GDU))*dxinv + &
                        HALF*(qgdnvyf(i,j+1,km,GDPRES) + qgdnvyf(i,j,km,GDPRES)) *  &
                             (qgdnvyf(i,j+1,km,GDV) - qgdnvyf(i,j,km,GDV))*dyinv
                end do
             end do
               
          end if
       end if
    enddo

    ! Deallocate arrays
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
      
  end subroutine umeth3d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi,          uin,uin_lo,uin_hi, &
                     q,c,gamc,csml,flatn,  q_lo,  q_hi, &
                     src,                src_lo,src_hi, &
                     srcQ,               srQ_lo,srQ_hi, &
                     courno,dx,dt,ngp,ngf)
    !
    !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
    !     if use_flattening=1.  Declared dimensions of q,c,gamc,csml,flatn are given
    !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
    !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
    !     routine that computes flatn).  
    !
    use mempool_module, only : bl_allocate, bl_deallocate
    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UESGS, UTEMP, &
                                   QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QESGS, QPRES, QTEMP, QGAME, QFS, QFX, &
                                   use_flattening, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1, &
                                   allow_negative_energy
    use flatten_module
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer          :: lo(3), hi(3)
    integer          :: uin_lo(3), uin_hi(3)
    integer          :: q_lo(3), q_hi(3)
    integer          :: src_lo(3), src_hi(3)
    integer          :: srQ_lo(3), srQ_hi(3)
    
    double precision :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    double precision :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    double precision :: c(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: gamc(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: csml(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3))
    double precision ::  src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    double precision :: srcQ(srQ_lo(1):srQ_hi(1),srQ_lo(2):srQ_hi(2),srQ_lo(3):srQ_hi(3),QVAR)
    double precision :: dx(3), dt, courno

    double precision, pointer:: dpdrho(:,:,:)
    double precision, pointer:: dpde(:,:,:)
!    double precision, pointer:: dpdX_er(:,:,:,:)

    integer          :: i, j, k
    integer          :: ngp, ngf, loq(3), hiq(3)
    integer          :: n, nq, ipassive
    double precision :: courx, coury, courz, courmx, courmy, courmz
    double precision :: kineng, rhoinv
    double precision :: dtdx, dtdy, dtdz
    
    type (eos_t) :: eos_state

    dtdx = dt/dx(1)
    dtdy = dt/dx(2)
    dtdz = dt/dx(3)

    do i=1,3
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo    
    
    call bl_allocate( dpdrho, q_lo(1),q_hi(1),q_lo(2),q_hi(2),q_lo(3),q_hi(3))
    call bl_allocate(   dpde, q_lo(1),q_hi(1),q_lo(2),q_hi(2),q_lo(3),q_hi(3))
!    call bl_allocate(dpdX_er, q_lo(1),q_hi(1),q_lo(2),q_hi(2),q_lo(3),q_hi(3),1,nspec)

    !
    ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
    ! The temperature is used as an initial guess for the eos call and will be overwritten.
    !
    do k = loq(3),hiq(3)
       do j = loq(2),hiq(2)
          do i = loq(1),hiq(1)             
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                print *,'>>> ... negative density ',uin(i,j,k,URHO)
                call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
             end if
          end do

          do i = loq(1),hiq(1)

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)
             
             q(i,j,k,QU) = uin(i,j,k,UMX) * rhoinv
             q(i,j,k,QV) = uin(i,j,k,UMY) * rhoinv
             q(i,j,k,QW) = uin(i,j,k,UMZ) * rhoinv

             ! Get the internal energy, which we'll use for determining the pressure.
             ! We use a dual energy formalism. If (E - K) < eta1 and eta1 is suitably small, 
             ! then we risk serious numerical truncation error in the internal energy.
             ! Therefore we'll use the result of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .gt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv
             else
                q(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv
             endif

             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)
             
             ! convert "rho K" to "K"
             if (QESGS .gt. -1) &
                  q(i,j,k,QESGS) = uin(i,j,k,UESGS)*rhoinv

          enddo
       enddo
    enddo

    ! Load passively advected quatities into q
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do k = loq(3),hiq(3)
          do j = loq(2),hiq(2)
             do i = loq(1),hiq(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo

    if (allow_negative_energy .eq. 0) eos_state % reset = .true.    
    
    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)
             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e
             q(i,j,k,QPRES)  = eos_state % p

             dpdrho(i,j,k)   = eos_state % dpdr_e
             dpde(i,j,k)     = eos_state % dpde
             c(i,j,k)        = eos_state % cs
             gamc(i,j,k)     = eos_state % gam1

             csml(i,j,k)     = max(small, small * c(i,j,k))

             q(i,j,k,QREINT) = q(i,j,k,QREINT) * q(i,j,k,QRHO)
             
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE
             
          enddo
       enddo
    enddo

    ! compute srcQ terms
    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)
             rhoinv = ONE/q(i,j,k,QRHO)
             srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
             srcQ(i,j,k,QU    ) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv
             srcQ(i,j,k,QV    ) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv
             srcQ(i,j,k,QW    ) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv
             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) &
                                                   - q(i,j,k,QV)*src(i,j,k,UMY) &
                                                   - q(i,j,k,QW)*src(i,j,k,UMZ) &
                                    + HALF * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2) * srcQ(i,j,k,QRHO)

             srcQ(i,j,k,QPRES ) = dpde(i,j,k)*(srcQ(i,j,k,QREINT) - &
                  q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)*rhoinv) * rhoinv + &
                  dpdrho(i,j,k)*srcQ(i,j,k,QRHO)! + &
!                                    sum(dpdX_er(i,j,k,:)*(src(i,j,k,UFS:UFS+nspec-1) - &
!                                                          q(i,j,k,QFS:QFS+nspec-1)*srcQ(i,j,k,QRHO))) &
!                                    /q(i,j,k,QRHO)

             if (QESGS .gt. -1) &
                  srcQ(i,j,k,QESGS) = src(i,j,k,UESGS)*rhoinv - q(i,j,k,QESGS) * srcQ(i,j,k,QRHO)

          enddo
       enddo
    enddo

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do k = loq(3), hiq(3)
          do j = loq(2), hiq(2)
             do i = loq(1), hiq(1)
                srcQ(i,j,k,nq) = ( src(i,j,k,n) - q(i,j,k,nq) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             enddo
          enddo
       enddo

    enddo

    ! Compute running max of Courant number over grids
    courmx = courno
    courmy = courno
    courmz = courno
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dtdx
             coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dtdy
             courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dtdz
             
             courmx = max( courmx, courx )
             courmy = max( courmy, coury )
             courmz = max( courmz, courz )
             
             if (courx .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (u+c) * dt / dx > 1 ', courx
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... u, c                ',q(i,j,k,QU), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if
             
             if (coury .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (v+c) * dt / dx > 1 ', coury
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... v, c                ',q(i,j,k,QV), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if
             
             if (courz .gt. ONE) then
                print *,'   '
                call bl_warning("Warning:: Castro_advection_3d.f90 :: CFL violation in ctoprim")
                print *,'>>> ... (w+c) * dt / dx > 1 ', courz
                print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                print *,'>>> ... w, c                ',q(i,j,k,QW), c(i,j,k)
                print *,'>>> ... density             ',q(i,j,k,QRHO)
             end if
             
          enddo
       enddo
    enddo

    courno = max( courmx, courmy, courmz )

    ! Compute flattening coef for slope calculations
    if (use_flattening == 1) then
       do n=1,3
          loq(n)=lo(n)-ngf
          hiq(n)=hi(n)+ngf
       enddo
       call uflaten(loq,hiq, &
                    q(q_lo(1),q_lo(2),q_lo(3),QPRES), &
                    q(q_lo(1),q_lo(2),q_lo(3),QU), &
                    q(q_lo(1),q_lo(2),q_lo(3),QV), &
                    q(q_lo(1),q_lo(2),q_lo(3),QW), &
                    flatn,q_lo,q_hi)
    else
       flatn = ONE
    endif

    call bl_deallocate( dpdrho)
    call bl_deallocate(   dpde)
!    call bl_deallocate(dpdX_er)
    
  end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine consup(uin,uin_lo,uin_hi, &
                    uout,uout_lo,uout_hi, &
                    src,src_lo,src_hi, &
                    flux1,flux1_lo,flux1_hi, &
                    flux2,flux2_lo,flux2_hi, &
                    flux3,flux3_lo,flux3_hi, &
                    qx,qx_lo,qx_hi, &
                    qy,qy_lo,qy_hi, &
                    qz,qz_lo,qz_hi, &
                    area1,area1_lo,area1_hi, &
                    area2,area2_lo,area2_hi, &
                    area3,area3_lo,area3_hi, &
                    vol,vol_lo,vol_hi, &
                    div,pdivu,lo,hi,dx,dt,E_added_flux, &
                    xmom_added_flux,ymom_added_flux,zmom_added_flux, &
                    verbose)

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UTEMP, normalize_species, hybrid_hydro, QVAR, NGDNV
    use bl_constants_module
    use hybrid_advection_module, only : hybrid_update
    use advection_util_module, only : normalize_species_fluxes

    integer          ::       lo(3),       hi(3)
    integer          ::   uin_lo(3),   uin_hi(3)
    integer          ::  uout_lo(3),  uout_hi(3)
    integer          ::   src_lo(3),   src_hi(3) 
    integer          :: flux1_lo(3), flux1_hi(3)
    integer          :: flux2_lo(3), flux2_hi(3)
    integer          :: flux3_lo(3), flux3_hi(3)
    integer          :: area1_lo(3), area1_hi(3)
    integer          :: area2_lo(3), area2_hi(3)
    integer          :: area3_lo(3), area3_hi(3)
    integer          ::    qx_lo(3),    qx_hi(3)
    integer          ::    qy_lo(3),    qy_hi(3)
    integer          ::    qz_lo(3),    qz_hi(3)
    integer          ::   vol_lo(3),   vol_hi(3)

    integer          :: verbose

    double precision :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    double precision :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    double precision ::   src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    double precision :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    double precision :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    double precision :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    double precision ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    double precision ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    double precision ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    double precision :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    double precision :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    double precision :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    double precision :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    double precision :: div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    double precision :: pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: dx(3), dt
    double precision :: E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux

    double precision :: div1, volinv
    integer          :: i, j, k, n
    
    do n = 1, NVAR
         
       if ( n.eq.UTEMP ) then
          
          flux1(:,:,:,n) = ZERO
          flux2(:,:,:,n) = ZERO
          flux3(:,:,:,n) = ZERO
          
       else

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1                  
                   div1 = FOURTH*(div(i,j,k) + div(i,j+1,k) + div(i,j,k+1) + div(i,j+1,k+1))
                   div1 = difmag*min(ZERO,div1)

                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx(1) * div1 * (uin(i,j,k,n)-uin(i-1,j,k,n))
                   flux1(i,j,k,n) = flux1(i,j,k,n) * area1(i,j,k) * dt
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                   div1 = difmag*min(ZERO,div1)

                   flux2(i,j,k,n) = flux2(i,j,k,n) + dx(2) * div1 * (uin(i,j,k,n)-uin(i,j-1,k,n))
                   flux2(i,j,k,n) = flux2(i,j,k,n) * area2(i,j,k) * dt
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag*min(ZERO,div1)

                   flux3(i,j,k,n) = flux3(i,j,k,n) + dx(3) * div1 * (uin(i,j,k,n)-uin(i,j,k-1,n))
                   flux3(i,j,k,n) = flux3(i,j,k,n) * area3(i,j,k) * dt
                enddo
             enddo
          enddo
          
       endif

    enddo

    if (normalize_species .eq. 1) &
         call normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                       flux2,flux2_lo,flux2_hi, &
                                       flux3,flux3_lo,flux3_hi, &
                                       lo,hi)

    do n = 1, NVAR

       ! pass temperature through
       if (n .eq. UTEMP) then
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   uout(i,j,k,n) = uin(i,j,k,n)
                enddo
             enddo
          enddo
       else 
          ! update everything else with fluxes
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)

                   volinv = ONE/vol(i,j,k)

                   uout(i,j,k,n) = uin(i,j,k,n) &
                          + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) * volinv

                   !
                   ! Add the p div(u) source term to (rho e)
                   !
                   if (n .eq. UEINT) then
                      uout(i,j,k,n) = uout(i,j,k,n) - dt * pdivu(i,j,k)
                   endif

                enddo
             enddo
          enddo
       endif
         
    enddo

    ! Add up some diagnostic quantities. Note that these are volumetric sums
    ! so we are not dividing by the cell volume.

    if (verbose .eq. 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                xmom_added_flux = xmom_added_flux + ( flux1(i,j,k,UMX) - flux1(i+1,j,k,UMX) &
                                                  +   flux2(i,j,k,UMX) - flux2(i,j+1,k,UMX) &
                                                  +   flux3(i,j,k,UMX) - flux3(i,j,k+1,UMX) )
                ymom_added_flux = ymom_added_flux + ( flux1(i,j,k,UMY) - flux1(i+1,j,k,UMY) &
                                                  +   flux2(i,j,k,UMY) - flux2(i,j+1,k,UMY) &
                                                  +   flux3(i,j,k,UMY) - flux3(i,j,k+1,UMY) )
                zmom_added_flux = zmom_added_flux + ( flux1(i,j,k,UMZ) - flux1(i+1,j,k,UMZ) &
                                                  +   flux2(i,j,k,UMZ) - flux2(i,j+1,k,UMZ) &
                                                  +   flux3(i,j,k,UMZ) - flux3(i,j,k+1,UMZ) )
                E_added_flux = E_added_flux + ( flux1(i,j,k,UEDEN) - flux1(i+1,j,k,UEDEN) &
                                            +   flux2(i,j,k,UEDEN) - flux2(i,j+1,k,UEDEN) &
                                            +   flux3(i,j,k,UEDEN) - flux3(i,j,k+1,UEDEN))
             enddo
          enddo
       enddo

    endif

    ! Now update the hybrid momenta, and overwrite the linear momenta accordingly.
    
    if (hybrid_hydro .eq. 1) then
       
       call hybrid_update(lo, hi, dx, dt, &
                          uin, uin_lo, uin_hi, &
                          uout, uout_lo, uout_hi, &
                          qx, qx_lo, qx_hi, &
                          qy, qy_lo, qy_hi, &
                          qz, qz_lo, qz_hi)
       
    endif
    
    
  end subroutine consup

end module advection_module
