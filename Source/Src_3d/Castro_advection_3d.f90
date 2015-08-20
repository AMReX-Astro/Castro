module advection_module

  implicit none

  private

  public umeth3d, ctoprim, divu, consup, enforce_minimum_density, normalize_new_species, &
       normalize_species_fluxes
  
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

  subroutine umeth3d(q, c, gamc, csml, flatn, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                     srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                     grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                     rot,  rt_l1, rt_l2, rt_l3, rt_h1, rt_h2, rt_h3, &
                     ilo1, ilo2, ilo3, ihi1, ihi2, ihi3, dx, dy, dz, dt, &
                     flux1, fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3, &
                     flux2, fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                     flux3, fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
                     ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3, &
                     ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                     ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3, &
                     ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                     ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3, &
                     ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                     pdivu, domlo, domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QFS, QFX, QTEMP, QREINT, ppm_type, &
                                   use_pslope, ppm_trace_grav, ppm_trace_rot, ppm_temp_fix, &
                                   do_grav, do_rotation, hybrid_riemann
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
    use trace_module, only : tracexy, tracez
    use transverse_module
    use ppm_module, only : ppm
    use slope_module, only : uslope, pslope
    use network
    use eos_module
    use eos_type_module
    use riemann_module, only: cmpflx, shock
    use bl_constants_module

    implicit none

    integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
    integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    integer gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3
    integer rt_l1, rt_l2, rt_l3, rt_h1, rt_h2, rt_h3
    integer ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
    integer fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3
    integer fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3
    integer fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3
    integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
    integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
    integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
    integer domlo(3),domhi(3)
    integer km,kc,kt,k3d,n
    integer i,j,iwave,idim
    
    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision  grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
    double precision   rot(rt_l1:rt_h1,rt_l2:rt_h2,rt_l3:rt_h3,3)
    double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,NVAR)
    double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,NVAR)
    double precision flux3(fd3_l1:fd3_h1,fd3_l2:fd3_h2,fd3_l3:fd3_h3,NVAR)
    double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
    double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
    double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
    double precision pdivu(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3)
    double precision dx, dy, dz, dt
    double precision dxinv, dyinv, dzinv
    double precision dtdx, dtdy, dtdz, hdt
    double precision cdtdx, cdtdy, cdtdz
    double precision hdtdx, hdtdy, hdtdz
    
    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer:: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
    double precision, pointer::qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
    double precision, pointer::qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)
    
    double precision, pointer::qmxy(:,:,:,:),qpxy(:,:,:,:)
    double precision, pointer::qmxz(:,:,:,:),qpxz(:,:,:,:)
    
    double precision, pointer::qmyx(:,:,:,:),qpyx(:,:,:,:)
    double precision, pointer::qmyz(:,:,:,:),qpyz(:,:,:,:)
    
    double precision, pointer::qmzx(:,:,:,:),qpzx(:,:,:,:)
    double precision, pointer::qmzy(:,:,:,:),qpzy(:,:,:,:)
    
    double precision, pointer::qxl(:,:,:,:),qxr(:,:,:,:)
    double precision, pointer::qyl(:,:,:,:),qyr(:,:,:,:)
    double precision, pointer::qzl(:,:,:,:),qzr(:,:,:,:)
    
    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    double precision, pointer::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)
    
    double precision, pointer::fxy(:,:,:,:),fxz(:,:,:,:)
    double precision, pointer::fyx(:,:,:,:),fyz(:,:,:,:)
    double precision, pointer::fzx(:,:,:,:),fzy(:,:,:,:)
    
    double precision, pointer:: pgdnvx(:,:,:), ugdnvx(:,:,:), gegdnvx(:,:,:)
    double precision, pointer:: pgdnvxf(:,:,:), ugdnvxf(:,:,:), gegdnvxf(:,:,:)
    double precision, pointer:: pgdnvtmpx(:,:,:), ugdnvtmpx(:,:,:), gegdnvtmpx(:,:,:)
    
    double precision, pointer:: pgdnvy(:,:,:), ugdnvy(:,:,:), gegdnvy(:,:,:)
    double precision, pointer:: pgdnvyf(:,:,:), ugdnvyf(:,:,:), gegdnvyf(:,:,:)
    double precision, pointer:: pgdnvtmpy(:,:,:), ugdnvtmpy(:,:,:), gegdnvtmpy(:,:,:)
    
    double precision, pointer:: pgdnvz(:,:,:), ugdnvz(:,:,:), gegdnvz(:,:,:)
    double precision, pointer:: pgdnvzf(:,:,:), ugdnvzf(:,:,:), gegdnvzf(:,:,:)
    double precision, pointer:: pgdnvtmpz1(:,:,:), ugdnvtmpz1(:,:,:), gegdnvtmpz1(:,:,:)
    double precision, pointer:: pgdnvtmpz2(:,:,:), ugdnvtmpz2(:,:,:), gegdnvtmpz2(:,:,:)
    
    double precision, pointer:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, pointer:: Ip_g(:,:,:,:,:,:), Im_g(:,:,:,:,:,:)
    double precision, pointer:: Ip_r(:,:,:,:,:,:), Im_r(:,:,:,:,:,:)
    double precision, pointer:: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    double precision, pointer :: shk(:,:,:)
    
    type (eos_t_3D) :: eos_state
    double precision :: rhoInv

    call bl_allocate ( pgdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    
    call bl_allocate ( pgdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( gegdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( pgdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    
    call bl_allocate ( qxm, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qxp, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qmxy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qpxy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qmxz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qpxz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qym, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qyp, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qmyx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qpyx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qmyz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qpyz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qzm, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qzp, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qxl, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qxr, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qyl, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qyr, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qzl, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qzr, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qmzx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qpzx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( qmzy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    call bl_allocate ( qpzy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)

    call bl_allocate ( fx, ilo1,ihi1+1,ilo2-1,ihi2+1,1,2,1,NVAR)
    call bl_allocate ( fy, ilo1-1,ihi1+1,ilo2,ihi2+1,1,2,1,NVAR)
    call bl_allocate ( fz, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,NVAR)

    call bl_allocate ( fxy, ilo1,ihi1+1,ilo2-1,ihi2+1,1,2,1,NVAR)
    call bl_allocate ( fxz, ilo1,ihi1+1,ilo2-1,ihi2+1,1,2,1,NVAR)

    call bl_allocate ( fyx, ilo1-1,ihi1+1,ilo2,ihi2+1,1,2,1,NVAR)
    call bl_allocate ( fyz, ilo1-1,ihi1+1,ilo2,ihi2+1,1,2,1,NVAR)

    call bl_allocate ( fzx, ilo1,ihi1,ilo2-1,ihi2+1,1,2,1,NVAR)
    call bl_allocate ( fzy, ilo1-1,ihi1+1,ilo2,ihi2,1,2,1,NVAR)

    if (ppm_type .gt. 0) then
       ! x-index, y-index, z-index, dim, characteristics, variables
       call bl_allocate ( Ip, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,QVAR)
       call bl_allocate ( Im, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,QVAR)
       
       ! for gravity (last index is x,y,z component)
       call bl_allocate ( Ip_g, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,3)
       call bl_allocate ( Im_g, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,3)
       
       ! for rotation (last index is x,y,z component)
       call bl_allocate ( Ip_r, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,3)
       call bl_allocate ( Im_r, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,3)
       
       ! for gamc -- needed for the reference state in eigenvectors
       call bl_allocate ( Ip_gc, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,1)
       call bl_allocate ( Im_gc, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,1)

       eos_state = eos_t_3D( (/ ilo1-1, ilo2-1, 1 /), (/ ihi1+1, ihi2+1, 1 /) )

    else
       call bl_allocate ( dqx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
       call bl_allocate ( dqy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
       call bl_allocate ( dqz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2,1,QVAR)
    end if

    ! for the hybrid Riemann solver
    call bl_allocate(shk, ilo1-1,ihi1+1,ilo2-1,ihi2+1,ilo3-1,ihi3+1)
    
    ! Local constants
    dxinv = ONE/dx
    dyinv = ONE/dy
    dzinv = ONE/dz
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
       call shock(q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
            shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
            ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz)
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

    do k3d = ilo3-1, ihi3+1

       ! Swap pointers to levels
       kt = km
       km = kc
       kc = kt

       if (ppm_type .gt. 0) then

          do n=1,QVAR
             call ppm(q(:,:,:,n  ),  qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
          end do

          if (do_grav .eq. 1 .and. ppm_trace_grav .eq. 1) then
             do n=1,3
                call ppm(grav(:,:,:,n),gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                         q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip_g(:,:,:,:,:,n),Im_g(:,:,:,:,:,n), &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
             enddo
          endif

          if (do_rotation .eq. 1 .and. ppm_trace_rot .eq. 1) then
             do n=1,3
                call ppm(rot(:,:,:,n),rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3, &
                         q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip_r(:,:,:,:,:,n),Im_r(:,:,:,:,:,n), &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
             enddo
          endif


          if (ppm_temp_fix /= 1) then
             call ppm(gamc(:,:,:),qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      Ip_gc(:,:,:,:,:,1),Im_gc(:,:,:,:,:,1), &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
          else          

             do iwave = 1, 3
                do idim = 1, 3

                   do j = ilo2-1, ihi2+1
                      do i = ilo1-1, ihi1+1
                         rhoInv = ONE / Ip(i,j,kc,idim,iwave,QRHO)
                      
                         eos_state % rho(i,j,1)   = Ip(i,j,kc,idim,iwave,QRHO)
                         eos_state % T(i,j,1)     = Ip(i,j,kc,idim,iwave,QTEMP)
                         
                         eos_state % xn(i,j,1,:)  = Ip(i,j,kc,idim,iwave,QFS:QFS+nspec-1) * rhoInv
                         eos_state % aux(i,j,1,:) = Ip(i,j,kc,idim,iwave,QFX:QFX+naux-1) * rhoInv
                      enddo
                   enddo

                   call eos(eos_input_rt, eos_state)

                   do j = ilo2-1, ihi2+1
                      do i = ilo1-1, ihi1+1
                         Ip(i,j,kc,idim,iwave,QPRES)  = eos_state % p(i,j,1)
                         Ip(i,j,kc,idim,iwave,QREINT) = eos_state % e(i,j,1) * Ip(i,j,kc,idim,iwave,QRHO)
                         Ip_gc(i,j,kc,idim,iwave,1)   = eos_state % gam1(i,j,1)
                      enddo
                   enddo

                   do j = ilo2-1, ihi2+1
                      do i = ilo1-1, ihi1+1
                         rhoInv = ONE / Im(i,j,kc,idim,iwave,QRHO)
                         
                         eos_state % rho(i,j,1)   = Im(i,j,kc,idim,iwave,QRHO)
                         eos_state % T(i,j,1)     = Im(i,j,kc,idim,iwave,QTEMP)

                         eos_state % xn(i,j,1,:)  = Im(i,j,kc,idim,iwave,QFS:QFS+nspec-1) * rhoInv
                         eos_state % aux(i,j,1,:) = Im(i,j,kc,idim,iwave,QFX:QFX+naux-1) * rhoInv
                      enddo
                   enddo

                   call eos(eos_input_rt, eos_state)

                   do j = ilo2-1, ihi2+1
                      do i = ilo1-1, ihi1+1
                         Im(i,j,kc,idim,iwave,QPRES)  = eos_state % p(i,j,1)
                         Im(i,j,kc,idim,iwave,QREINT) = eos_state % e(i,j,1) * Im(i,j,kc,idim,iwave,QRHO)
                         Im_gc(i,j,kc,idim,iwave,1)   = eos_state % gam1(i,j,1)
                      enddo
                   enddo

                enddo
             enddo

          endif

          ! Compute U_x and U_y at kc (k3d)
          call tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip,Im,Ip_g,Im_g,Ip_r,Im_r,Ip_gc,Im_gc, &
                           qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           ilo1,ilo2,ihi1,ihi2,dt,kc,k3d)

       else

          ! Compute all slopes at kc (k3d)
          call uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      dqx,dqy,dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      ilo1,ilo2,ihi1,ihi2,kc,k3d,QVAR)
          
          if (use_pslope .eq. 1) &
               call pslope(q(:,:,:,QPRES),q(:,:,:,QRHO), &
                           flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           dqx(:,:,:,QPRES),dqy(:,:,:,QPRES),dqz(:,:,:,QPRES), &
                           ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                           ilo1,ilo2,ihi1,ihi2,kc,k3d,dx,dy,dz)

          ! Compute U_x and U_y at kc (k3d)
          call tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       dqx,dqy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)
          
       end if

       ! Compute \tilde{F}^x at kc (k3d)
       call cmpflx(qxm,qxp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                   ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   1,ilo1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                   ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   2,ilo1-1,ihi1+1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)
       
       ! Compute U'^y_x at kc (k3d)
       call transy1(qxm,qmxy,qxp,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                    fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                    ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                    gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                    cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,k3d)

       ! Compute U'^x_y at kc (k3d)
       call transx1(qym,qmyx,qyp,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                    fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                    ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                    gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                    cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,k3d)

       ! Compute F^{x|y} at kc (k3d)
       call cmpflx(qmxy,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                   ugdnvtmpx,pgdnvtmpx,gegdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                   ugdnvtmpy,pgdnvtmpy,gegdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)

       if (k3d.ge.ilo3) then
          
          ! Compute U_z at kc (k3d)
          if (ppm_type .gt. 0) then
             call tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             Ip,Im,Ip_g,Im_g,Ip_r,Im_r,Ip_gc,Im_gc, &
                             qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             ilo1,ilo2,ihi1,ihi2,dt,km,kc,k3d)
          else
             call tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)
          end if

          ! Compute \tilde{F}^z at kc (k3d)
          call cmpflx(qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                      ugdnvz,pgdnvz,gegdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                      3,ilo1-1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

          ! Compute U'^y_z at kc (k3d)
          call transy2(qzm,qmzy,qzp,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                       ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,km,k3d)

          ! Compute U'^x_z at kc (k3d)
          call transx2(qzm,qmzx,qzp,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                       ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,km,k3d)

          ! Compute F^{z|x} at kc (k3d)
          call cmpflx(qmzx,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                      ugdnvtmpz1,pgdnvtmpz1,gegdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                      3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                      ugdnvtmpz2,pgdnvtmpz2,gegdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &                       
                      3,ilo1-1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)
          
          ! Compute U''_z at kc (k3d)
          call transxy(qzm,qzl,qzp,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                       fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                       ugdnvtmpx,pgdnvtmpx,gegdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       ugdnvtmpy,pgdnvtmpy,gegdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                       gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                       grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
                       rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3,&
                       hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      flux3,fd3_l1,fd3_l2,fd3_l3,fd3_h1,fd3_h2,fd3_h3, &
                      ugdnvzf,pgdnvzf,gegdnvzf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                      3,ilo1,ihi1,ilo2,ihi2,kc,k3d,k3d,domlo,domhi)

          do j=ilo2-1,ihi2+1
             do i=ilo1-1,ihi1+1
                ugdnvz_out(i,j,k3d) = ugdnvzf(i,j,kc)
             end do
          end do

          if (k3d .ge. ilo3+1 .and. k3d .le. ihi3+1) then
             do j = ilo2,ihi2
                do i = ilo1,ihi1
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        HALF*(pgdnvzf(i,j,kc)+pgdnvzf(i,j,km)) * &
                              (ugdnvzf(i,j,kc)-ugdnvzf(i,j,km))*dzinv
                end do
             end do
          end if
          
          if (k3d.gt.ilo3) then

             ! Compute U'^z_x and U'^z_y at km (k3d-1) -- note flux3 has physical index
             call transz(qxm,qmxz,qxp,qpxz, &
                         qym,qmyz,qyp,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                         ugdnvz,pgdnvz,gegdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         cdtdz,ilo1-1,ihi1+1,ilo2-1,ihi2+1,km,kc,k3d)
         
             ! Compute F^{x|z} at km (k3d-1)
             call cmpflx(qmxz,qpxz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                         ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1,domlo,domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                         ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         2,ilo1-1,ihi1+1,ilo2,ihi2+1,km,km,k3d-1,domlo,domhi)

             ! Compute U''_x at km (k3d-1)
             call transyz(qxm,qxl,qxp,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                          fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                          ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          ugdnvtmpz2,pgdnvtmpz2,gegdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                          grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
                          rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3,&
                          hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1)

             ! Compute U''_y at km (k3d-1)
             call transxz(qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                          fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                          ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          ugdnvtmpz1,pgdnvtmpz1,gegdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                          grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3,&
                          rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3,&
                          hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         flux1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
                         ugdnvxf,pgdnvxf,gegdnvxf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         1,ilo1,ihi1+1,ilo2,ihi2,km,k3d-1,k3d-1,domlo,domhi)
             
             do j=ilo2-1,ihi2+1
                do i=ilo1-1,ihi1+2
                   ugdnvx_out(i,j,k3d-1) = ugdnvxf(i,j,km)
                end do
             end do
             
             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         flux2,fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3, &
                         ugdnvyf,pgdnvyf,gegdnvyf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         2,ilo1,ihi1,ilo2,ihi2+1,km,k3d-1,k3d-1,domlo,domhi)

             do j=ilo2-1,ihi2+2
                do i=ilo1-1,ihi1+1
                   ugdnvy_out(i,j,k3d-1) = ugdnvyf(i,j,km)
                end do
             end do

             do j = ilo2,ihi2
                do i = ilo1,ihi1
                   pdivu(i,j,k3d-1) = pdivu(i,j,k3d-1) +  &
                        HALF*(pgdnvxf(i+1,j,km) + pgdnvxf(i,j,km)) *  &
                              (ugdnvxf(i+1,j,km)-ugdnvxf(i,j,km))*dxinv + &
                        HALF*(pgdnvyf(i,j+1,km) + pgdnvyf(i,j,km)) *  &
                              (ugdnvyf(i,j+1,km)-ugdnvyf(i,j,km))*dyinv
                end do
             end do
               
          end if
       end if
    enddo

    ! Deallocate arrays
    call bl_deallocate ( pgdnvx)
    call bl_deallocate ( ugdnvx)
    call bl_deallocate (gegdnvx)

    call bl_deallocate ( pgdnvxf)
    call bl_deallocate ( ugdnvxf)
    call bl_deallocate (gegdnvxf)

    call bl_deallocate ( pgdnvtmpx)
    call bl_deallocate ( ugdnvtmpx)
    call bl_deallocate (gegdnvtmpx)
    
    call bl_deallocate ( pgdnvy)
    call bl_deallocate ( ugdnvy)
    call bl_deallocate ( gegdnvy)

    call bl_deallocate ( pgdnvyf)
    call bl_deallocate ( ugdnvyf)
    call bl_deallocate (gegdnvyf)

    call bl_deallocate ( pgdnvtmpy)
    call bl_deallocate ( ugdnvtmpy)
    call bl_deallocate (gegdnvtmpy)

    call bl_deallocate ( pgdnvz)
    call bl_deallocate ( ugdnvz)
    call bl_deallocate (gegdnvz)

    call bl_deallocate ( pgdnvzf)
    call bl_deallocate ( ugdnvzf)
    call bl_deallocate (gegdnvzf)

    call bl_deallocate ( pgdnvtmpz1)
    call bl_deallocate ( ugdnvtmpz1)
    call bl_deallocate (gegdnvtmpz1)

    call bl_deallocate ( pgdnvtmpz2)
    call bl_deallocate ( ugdnvtmpz2)
    call bl_deallocate (gegdnvtmpz2)
    
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
       
       call bl_deallocate ( Ip_g)
       call bl_deallocate ( Im_g)
       
       call bl_deallocate ( Ip_r)
       call bl_deallocate ( Im_r)
       
       call bl_deallocate ( Ip_gc)
       call bl_deallocate ( Im_gc)

       call eos_deallocate(eos_state)
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

  subroutine ctoprim(lo,hi,          uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                     q,c,gamc,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                     src,                src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                     srcQ,               srQ_l1,srQ_l2,srQ_l3,srQ_h1,srQ_h2,srQ_h3, &
                     courno,dx,dy,dz,dt,ngp,ngf)
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
    use eos_type_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UESGS, UTEMP, UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QESGS, QPRES, QTEMP, QGAME, QFA, QFS, QFX, &
                                   nadv, allow_negative_energy, small_temp, use_flattening, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1
    
    use flatten_module
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer lo(3), hi(3)
    integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer srQ_l1,srQ_l2,srQ_l3,srQ_h1,srQ_h2,srQ_h3
    
    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    double precision :: c(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: gamc(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision ::  src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,srQ_l3:srQ_h3,QVAR)
    double precision :: dx, dy, dz, dt, courno

    double precision, pointer:: dpdrho(:,:,:)
    double precision, pointer:: dpde(:,:,:)
    double precision, pointer:: dpdX_er(:,:,:,:)

    integer          :: i, j, k
    integer          :: ngp, ngf, loq(3), hiq(3)
    integer          :: n, nq
    integer          :: iadv, ispec, iaux
    double precision :: courx, coury, courz, courmx, courmy, courmz
    double precision :: kineng, rhoinv
    double precision :: dtdx, dtdy, dtdz

    integer :: ipassive

    type (eos_t_3D) :: eos_state

    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    call bl_allocate( dpdrho, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
    call bl_allocate(   dpde, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
!    call bl_allocate(dpdX_er, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3,1,nspec)

    do i=1,3
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo

    eos_state = eos_t_3D(loq, hiq)

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
             q(i,j,k,QU) = uin(i,j,k,UMX)*rhoinv
             q(i,j,k,QV) = uin(i,j,k,UMY)*rhoinv
             q(i,j,k,QW) = uin(i,j,k,UMZ)*rhoinv

             ! Get the internal energy, which we'll use for determining the pressure.
             ! We use a dual energy formalism. If (E - K) < eta1 and eta1 is suitably small, 
             ! then we risk serious numerical truncation error in the internal energy.
             ! Therefore we'll use the result of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .lt. dual_energy_eta1) then
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

    ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
    do iadv = 1, nadv
       n = UFA + iadv - 1
       nq = QFA + iadv - 1
       do k = loq(3),hiq(3)
          do j = loq(2),hiq(2)
             do i = loq(1),hiq(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo
      
    ! Load chemical species, c, into q, assuming they arrived in uin as rho.c
    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       do k = loq(3),hiq(3)
          do j = loq(2),hiq(2)
             do i = loq(1),hiq(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
                eos_state % xn(i,j,k,ispec) = q(i,j,k,nq)
             enddo
          enddo
       enddo
    enddo
      
    ! Load auxiliary variables which are needed in the EOS
    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       do k = loq(3),hiq(3)
          do j = loq(2),hiq(2)
             do i = loq(1),hiq(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
                eos_state % aux(i,j,k,iaux) = q(i,j,k,nq)
             enddo
          enddo
       enddo
    enddo

    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)
             eos_state % T(i,j,k)   = q(i,j,k,QTEMP )
             eos_state % rho(i,j,k) = q(i,j,k,QRHO  )
             eos_state % e(i,j,k)   = q(i,j,k,QREINT)
          enddo
       enddo
    enddo

    call eos(eos_input_re, eos_state)

    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)
             q(i,j,k,QTEMP)  = eos_state % T(i,j,k)
             q(i,j,k,QREINT) = eos_state % e(i,j,k)
             q(i,j,k,QPRES)  = eos_state % p(i,j,k)

             dpdrho(i,j,k)   = eos_state % dpdr_e(i,j,k)
             dpde(i,j,k)     = eos_state % dpde(i,j,k)
             c(i,j,k)        = eos_state % cs(i,j,k)
             gamc(i,j,k)     = eos_state % gam1(i,j,k)

             csml(i,j,k)     = max(small, small * c(i,j,k))

             q(i,j,k,QREINT) = q(i,j,k,QREINT) * q(i,j,k,QRHO)
             
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE
             
          enddo
       enddo
    enddo

    ! compute srcQ terms
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
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

       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
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
                    q(q_l1,q_l2,q_l3,QPRES), &
                    q(q_l1,q_l2,q_l3,QU), &
                    q(q_l1,q_l2,q_l3,QV), &
                    q(q_l1,q_l2,q_l3,QW), &
                    flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)
    else
       flatn = ONE
    endif

    call bl_deallocate( dpdrho)
    call bl_deallocate(   dpde)
!    call bl_deallocate(dpdX_er)
    call eos_deallocate(eos_state)
    
  end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                    uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                    src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                    flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                    flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                    flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                    area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                    area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                    area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                    vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                    div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux,&
                    xmom_added_flux,ymom_added_flux,zmom_added_flux)

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UTEMP, normalize_species
    use bl_constants_module

    implicit none

    integer lo(3), hi(3)
    integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    integer   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3 
    integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    integer area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3
    integer area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3
    integer area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3
    integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3

    double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision   src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
    double precision area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
    double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision dx, dy, dz, dt
    double precision E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux

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
                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx*div1*(uin(i,j,k,n)-uin(i-1,j,k,n))
                   flux1(i,j,k,n) = flux1(i,j,k,n) * area1(i,j,k) * dt
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j,k+1) + div(i+1,j,k+1))
                   div1 = difmag*min(ZERO,div1)
                   flux2(i,j,k,n) = flux2(i,j,k,n) + dy*div1*(uin(i,j,k,n)-uin(i,j-1,k,n))
                   flux2(i,j,k,n) = flux2(i,j,k,n) * area2(i,j,k) * dt
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag*min(ZERO,div1)
                   flux3(i,j,k,n) = flux3(i,j,k,n) + dz*div1*(uin(i,j,k,n)-uin(i,j,k-1,n))
                   flux3(i,j,k,n) = flux3(i,j,k,n) * area3(i,j,k) * dt
                enddo
             enddo
          enddo
          
       endif

    enddo

    if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
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
          ! update everything else with fluxes and source terms
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   volinv = ONE/vol(i,j,k)

                   uout(i,j,k,n) = uin(i,j,k,n) &
                          + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) * volinv &
                          +   dt * src(i,j,k,n)
                   !
                   ! Add the source term to (rho e)
                   !
                   if (n .eq. UEINT) then
                      uout(i,j,k,n) = uout(i,j,k,n) - dt * pdivu(i,j,k)
                   else if (n .eq. UMX) then
                      xmom_added_flux = xmom_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) * volinv
                   else if (n .eq. UMY) then
                      ymom_added_flux = ymom_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) * volinv
                   else if (n .eq. UMZ) then
                      zmom_added_flux = zmom_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) * volinv
                   else if (n .eq. UEDEN) then
                      E_added_flux = E_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) * volinv
                   endif
                enddo
             enddo
          enddo
       endif
         
    enddo

  end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,dx,dy,dz, &
                  div,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3)
    
    use meth_params_module, only : QU, QV, QW
    use bl_constants_module
    
    implicit none

    integer          :: lo(3),hi(3)
    integer          :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    integer          :: div_l1,div_l2,div_l3,div_h1,div_h2,div_h3
    double precision :: dx, dy, dz
    double precision :: div(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3)
    double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,*)

    integer          :: i, j, k
    double precision :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE/dx
    dyinv = ONE/dy
    dzinv = ONE/dz

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             
             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo
    
  end subroutine divu

! ::
! :: ----------------------------------------------------------
! ::

  subroutine normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_l3, &
                                      flux1_h1,flux1_h2,flux1_h3, &
                                      flux2,flux2_l1,flux2_l2,flux2_l3, &
                                      flux2_h1,flux2_h2,flux2_h3, &
                                      flux3,flux3_l1,flux3_l2,flux3_l3, &
                                      flux3_h1,flux3_h2,flux3_h3, &
                                      lo,hi)
    
    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer          :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer          :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: sum,fac
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux1(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux1(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux1(i,j,k,n) = flux1(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux2(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux2(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux2(i,j,k,n) = flux2(i,j,k,n) * fac
             end do
          end do
       end do
    end do

    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux3(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux3(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux3(i,j,k,n) = flux3(i,j,k,n) * fac
             end do
          end do
       end do
    end do

  end subroutine normalize_species_fluxes

! ::
! :: ----------------------------------------------------------
! ::

  subroutine enforce_minimum_density(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                                     uout,uout_l1,uout_l2,uout_l3, &
                                     uout_h1,uout_h2,uout_h3, &
                                     lo,hi,mass_added,eint_added,eden_added,verbose)
    
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, UFX, &
                                     UFA, small_dens, small_temp, nadv
    use bl_constants_module
    use eos_type_module
    use eos_module

    implicit none

    integer          :: lo(3), hi(3), verbose
    integer          ::  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
    integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    double precision ::  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision :: mass_added, eint_added, eden_added
    
    ! Local variables
    integer          :: i,ii,j,jj,k,kk,n
    integer          :: i_set, j_set, k_set
    double precision :: max_dens
    
    double precision :: initial_mass, final_mass
    double precision :: initial_eint, final_eint
    double precision :: initial_eden, final_eden

    type (eos_t) :: eos_state
    
    initial_mass = ZERO
      final_mass = ZERO

    initial_eint = ZERO
      final_eint = ZERO

    initial_eden = ZERO
      final_eden = ZERO

    max_dens = ZERO

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             initial_mass = initial_mass + uout(i,j,k,URHO )
             initial_eint = initial_eint + uout(i,j,k,UEINT)
             initial_eden = initial_eden + uout(i,j,k,UEDEN)
             
             if (uout(i,j,k,URHO) .eq. ZERO) then
                
                print *,'DENSITY EXACTLY ZERO AT CELL ',i,j,k
                print *,'  in grid ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
                call bl_error("Error:: Castro_3d.f90 :: enforce_minimum_density")
                
             else if (uout(i,j,k,URHO) < small_dens) then
                
                max_dens = uout(i,j,k,URHO)
                i_set = i
                j_set = j
                k_set = k
                do kk = -1,1
                   do jj = -1,1
                      do ii = -1,1
                         if (i+ii.ge.lo(1) .and. j+jj.ge.lo(2) .and. k+kk.ge.lo(3) .and. &
                             i+ii.le.hi(1) .and. j+jj.le.hi(2) .and. k+kk.le.hi(3)) then
                              if (uout(i+ii,j+jj,k+kk,URHO) .gt. max_dens) then
                                  i_set = i+ii
                                  j_set = j+jj
                                  k_set = k+kk
                                  max_dens = uout(i_set,j_set,k_set,URHO)
                              endif
                         endif
                      end do
                   end do
                end do

                ! If no neighboring zones are above small_dens, our only recourse 
                ! is to set the density equal to small_dens, and the temperature 
                ! equal to small_temp. We set the velocities to zero, 
                ! though any choice here would be arbitrary.

                if (max_dens < small_dens) then

                   do n = UFS, UFS+nspec-1
                      uout(i,j,k,n) = uout(i_set,j_set,k_set,n) * (small_dens / uout(i,j,k,URHO))
                   end do
                   do n = UFX, UFX+naux-1
                      uout(i,j,k,n) = uout(i_set,j_set,k_set,n) * (small_dens / uout(i,j,k,URHO))
                   end do
                   do n = UFA, UFA+nadv-1
                      uout(i,j,k,n) = uout(i_set,j_set,k_set,n) * (small_dens / uout(i,j,k,URHO))
                   end do

                   eos_state % rho = small_dens
                   eos_state % T   = small_temp
                   eos_state % xn  = uout(i,j,k,UFS:UFS+nspec-1) / uout(i,j,k,URHO)

                   call eos(eos_input_rt, eos_state)

                   uout(i,j,k,URHO ) = eos_state % rho
                   uout(i,j,k,UTEMP) = eos_state % T

                   uout(i,j,k,UMX  ) = ZERO
                   uout(i,j,k,UMY  ) = ZERO
                   uout(i,j,k,UMZ  ) = ZERO

                   uout(i,j,k,UEINT) = eos_state % rho * eos_state % e
                   uout(i,j,k,UEDEN) = uout(i,j,k,UEINT)

                endif
                
                if (verbose .gt. 0) then
                   if (uout(i,j,k,URHO) < ZERO) then
                      print *,'   '
                      print *,'>>> RESETTING NEG.  DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',uout(i_set,j_set,k_set,URHO)
                      print *,'   '
                   else
                      print *,'   '
                      print *,'>>> RESETTING SMALL DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',uout(i_set,j_set,k_set,URHO)
                      print *,'   '
                   end if
                end if
                
                uout(i,j,k,URHO ) = uout(i_set,j_set,k_set,URHO )
                uout(i,j,k,UTEMP) = uout(i_set,j_set,k_set,UTEMP)
                uout(i,j,k,UEINT) = uout(i_set,j_set,k_set,UEINT)
                uout(i,j,k,UEDEN) = uout(i_set,j_set,k_set,UEDEN)
                uout(i,j,k,UMX  ) = uout(i_set,j_set,k_set,UMX  )
                uout(i,j,k,UMY  ) = uout(i_set,j_set,k_set,UMY  )
                uout(i,j,k,UMZ  ) = uout(i_set,j_set,k_set,UMZ  )
   
                do n = UFS, UFS+nspec-1
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                do n = UFX, UFX+naux-1
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                do n = UFA, UFA+nadv-1
                   uout(i,j,k,n) = uout(i_set,j_set,k_set,n)
                end do
                
             end if

             final_mass = final_mass + uout(i,j,k,URHO )
             final_eint = final_eint + uout(i,j,k,UEINT)
             final_eden = final_eden + uout(i,j,k,UEDEN)                
             
          enddo
       enddo
    enddo

    if ( max_dens /= ZERO ) then
       mass_added = mass_added + final_mass - initial_mass
       eint_added = eint_added + final_eint - initial_eint
       eden_added = eden_added + final_eden - initial_eden
    endif

    call eos_deallocate(eos_state)
    
  end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
    double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: fac,sum
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + u(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = u(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                u(i,j,k,n) = u(i,j,k,n) * fac
             end do
          end do
       end do
    end do
    
  end subroutine normalize_new_species

end module advection_module
