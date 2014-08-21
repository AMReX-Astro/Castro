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

    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QFS, QFX, QTEMP, QREINT, ppm_type, &
                                   use_pslope, ppm_trace_grav, ppm_temp_fix
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
    use trace_module, only : tracexy, tracez
    use transverse_module
    use ppm_module, only : ppm
    use slope_module, only : uslope, pslope
    use network
    use eos_module
    use eos_type_module
    use riemann_module, only: cmpflx
    use bl_constants_module

    implicit none

    integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
    integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    integer gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3
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
    double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,NVAR)
    double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,NVAR)
    double precision flux3(fd3_l1:fd3_h1,fd3_l2:fd3_h2,fd3_l3:fd3_h3,NVAR)
    double precision ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
    double precision ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
    double precision ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
    double precision pdivu(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3)
    double precision dx, dy, dz, dt
    double precision dtdx, dtdy, dtdz, hdt
    double precision cdtdx, cdtdy, cdtdz
    double precision hdtdx, hdtdy, hdtdz
    
    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
    double precision, allocatable::qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
    double precision, allocatable::qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)
    
    double precision, allocatable::qmxy(:,:,:,:),qpxy(:,:,:,:)
    double precision, allocatable::qmxz(:,:,:,:),qpxz(:,:,:,:)
    
    double precision, allocatable::qmyx(:,:,:,:),qpyx(:,:,:,:)
    double precision, allocatable::qmyz(:,:,:,:),qpyz(:,:,:,:)
    
    double precision, allocatable::qmzx(:,:,:,:),qpzx(:,:,:,:)
    double precision, allocatable::qmzy(:,:,:,:),qpzy(:,:,:,:)
    
    double precision, allocatable::qxl(:,:,:,:),qxr(:,:,:,:)
    double precision, allocatable::qyl(:,:,:,:),qyr(:,:,:,:)
    double precision, allocatable::qzl(:,:,:,:),qzr(:,:,:,:)
    
    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    double precision, allocatable::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)
    
    double precision, allocatable::fxy(:,:,:,:),fxz(:,:,:,:)
    double precision, allocatable::fyx(:,:,:,:),fyz(:,:,:,:)
    double precision, allocatable::fzx(:,:,:,:),fzy(:,:,:,:)
    
    double precision, allocatable:: pgdnvx(:,:,:), ugdnvx(:,:,:), gegdnvx(:,:,:)
    double precision, allocatable:: pgdnvxf(:,:,:), ugdnvxf(:,:,:), gegdnvxf(:,:,:)
    double precision, allocatable:: pgdnvtmpx(:,:,:), ugdnvtmpx(:,:,:), gegdnvtmpx(:,:,:)
    
    double precision, allocatable:: pgdnvy(:,:,:), ugdnvy(:,:,:), gegdnvy(:,:,:)
    double precision, allocatable:: pgdnvyf(:,:,:), ugdnvyf(:,:,:), gegdnvyf(:,:,:)
    double precision, allocatable:: pgdnvtmpy(:,:,:), ugdnvtmpy(:,:,:), gegdnvtmpy(:,:,:)
    
    double precision, allocatable:: pgdnvz(:,:,:), ugdnvz(:,:,:), gegdnvz(:,:,:)
    double precision, allocatable:: pgdnvzf(:,:,:), ugdnvzf(:,:,:), gegdnvzf(:,:,:)
    double precision, allocatable:: pgdnvtmpz1(:,:,:), ugdnvtmpz1(:,:,:), gegdnvtmpz1(:,:,:)
    double precision, allocatable:: pgdnvtmpz2(:,:,:), ugdnvtmpz2(:,:,:), gegdnvtmpz2(:,:,:)
    
    double precision, allocatable:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, allocatable:: Ip_g(:,:,:,:,:,:), Im_g(:,:,:,:,:,:)
    double precision, allocatable:: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)
    
    type (eos_t) :: eos_state

    allocate ( pgdnvx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvxf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvxf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvxf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvtmpx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvtmpx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvtmpx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    
    allocate ( pgdnvy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( gegdnvy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvyf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvyf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvyf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvtmpy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvtmpy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvtmpy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvzf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvzf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvzf(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvtmpz1(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvtmpz1(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvtmpz1(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))

    allocate ( pgdnvtmpz2(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate ( ugdnvtmpz2(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    allocate (gegdnvtmpz2(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2))
    
    allocate ( dqx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( dqy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( dqz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    
    allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qmxy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qpxy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qmxz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qpxz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qmyx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qpyx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qmyz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qpyz(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qzm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qzp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qxl(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qxr(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qyl(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qyr(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qzl(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qzr(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qmzx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qpzx(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( qmzy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))
    allocate ( qpzy(ilo1-1:ihi1+2,ilo2-1:ihi2+2,2,QVAR))

    allocate ( fx(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
    allocate ( fy(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
    allocate ( fz(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))

    allocate ( fxy(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))
    allocate ( fxz(ilo1:ihi1+1,ilo2-1:ihi2+1,2,NVAR))

    allocate ( fyx(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))
    allocate ( fyz(ilo1-1:ihi1+1,ilo2:ihi2+1,2,NVAR))

    allocate ( fzx(ilo1:ihi1,ilo2-1:ihi2+1,2,NVAR))
    allocate ( fzy(ilo1-1:ihi1+1,ilo2:ihi2,2,NVAR))

    ! x-index, y-index, z-index, dim, characteristics, variables
    allocate ( Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,QVAR))
    allocate ( Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,QVAR))
    
    ! for gravity (last index is x,y,z component)
    allocate ( Ip_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,3))
    allocate ( Im_g(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,3))

    ! for gamc -- needed for the reference state in eigenvectors
    allocate ( Ip_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,1))
    allocate ( Im_gc(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,3,1))

    ! Local constants
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz
    hdt = HALF*dt
    hdtdx = HALF*dtdx
    hdtdy = HALF*dtdy
    hdtdz = HALF*dtdz
    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD

    ! Initialize pdivu to zero
    pdivu(:,:,:) = ZERO

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

          if (ppm_trace_grav .eq. 1) then
             do n=1,3
                call ppm(grav(:,:,:,n),gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                         q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip_g(:,:,:,:,:,n),Im_g(:,:,:,:,:,n), &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
             enddo
          endif


          if (ppm_temp_fix /= 1) then
             call ppm(gamc(:,:,:),qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      Ip_gc(:,:,:,:,:,1),Im_gc(:,:,:,:,:,1), &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
          endif
          
          ! temperature-based PPM
          if (ppm_temp_fix == 1) then
             do j = ilo2-1, ihi2+1
                do i = ilo1-1, ihi1+1
                   do idim = 1, 3
                      do iwave = 1, 3
                         eos_state % rho = Ip(i,j,kc,idim,iwave,QRHO)
                         eos_state % T   = Ip(i,j,kc,idim,iwave,QTEMP)
                         eos_state % xn  = Ip(i,j,kc,idim,iwave,QFS:QFS-1+nspec)
                         eos_state % aux = Ip(i,j,kc,idim,iwave,QFX:QFX-1+naux)

                         call eos(eos_input_rt, eos_state, .false.)
                         
                         Ip(i,j,kc,idim,iwave,QPRES) = eos_state%p
                         Ip(i,j,kc,idim,iwave,QREINT) = Ip(i,j,kc,idim,iwave,QRHO)*eos_state%e
                         Ip_gc(i,j,kc,idim,iwave,1) = eos_state%gam1

                         
                         eos_state % rho = Im(i,j,kc,idim,iwave,QRHO)
                         eos_state % T   = Im(i,j,kc,idim,iwave,QTEMP)
                         eos_state % xn  = Im(i,j,kc,idim,iwave,QFS:QFS-1+nspec)
                         eos_state % aux = Im(i,j,kc,idim,iwave,QFX:QFX-1+naux)
                         
                         call eos(eos_input_rt, eos_state, .false.)
                         
                         Im(i,j,kc,idim,iwave,QPRES) = eos_state%p
                         Im(i,j,kc,idim,iwave,QREINT) = Im(i,j,kc,idim,iwave,QRHO)*eos_state%e
                         Im_gc(i,j,kc,idim,iwave,1) = eos_state%gam1
                         
                      enddo
                   enddo
                enddo
             enddo
             
          endif

          ! Compute U_x and U_y at kc (k3d)
          call tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip,Im,Ip_g,Im_g,Ip_gc,Im_gc, &
                           qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                           gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

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
                   1,ilo1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                   ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
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
                   1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                   ugdnvtmpy,pgdnvtmpy,gegdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)

       if (k3d.ge.ilo3) then
          
          ! Compute U_z at kc (k3d)
          if (ppm_type .gt. 0) then
             call tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             Ip,Im,Ip_g,Im_g,Ip_gc,Im_gc, &
                             qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                             gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)
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
                      3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                      ugdnvtmpz2,pgdnvtmpz2,gegdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
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
                       hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      flux3,fd3_l1,fd3_l2,fd3_l3,fd3_h1,fd3_h2,fd3_h3, &
                      ugdnvzf,pgdnvzf,gegdnvzf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
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
                              (ugdnvzf(i,j,kc)-ugdnvzf(i,j,km))/dz
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
                         1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1,domlo,domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                         ugdnvy,pgdnvy,gegdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
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
                          hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         flux1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
                         ugdnvxf,pgdnvxf,gegdnvxf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
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
                              (ugdnvxf(i+1,j,km)-ugdnvxf(i,j,km))/dx + &
                        HALF*(pgdnvyf(i,j+1,km) + pgdnvyf(i,j,km)) *  &
                              (ugdnvyf(i,j+1,km)-ugdnvyf(i,j,km))/dy
                end do
             end do
               
          end if
       end if
    enddo

    ! Deallocate arrays
    deallocate(pgdnvx,ugdnvx,gegdnvx)
    deallocate(pgdnvxf,ugdnvxf,gegdnvxf)
    deallocate(pgdnvtmpx,ugdnvtmpx,gegdnvtmpx)
    deallocate(pgdnvy,ugdnvy,gegdnvy)
    deallocate(pgdnvyf,ugdnvyf,gegdnvyf)
    deallocate(pgdnvtmpy,ugdnvtmpy,gegdnvtmpy)
    deallocate(pgdnvz,ugdnvz,gegdnvz)
    deallocate(pgdnvtmpz1,ugdnvtmpz1,gegdnvtmpz1)
    deallocate(pgdnvtmpz2,ugdnvtmpz2,gegdnvtmpz2)
    deallocate(pgdnvzf,ugdnvzf,gegdnvzf)
    deallocate(dqx,dqy,dqz)
    deallocate(qxm,qxp)
    deallocate(qmxy,qpxy)
    deallocate(qmxz,qpxz)
    deallocate(qym,qyp)
    deallocate(qmyx,qpyx)
    deallocate(qmyz,qpyz)
    deallocate(qzm,qzp)
    deallocate(qxl,qxr,qyl,qyr,qzl,qzr)
    deallocate(qmzx,qpzx)
    deallocate(qmzy,qpzy)
    deallocate(fx,fy,fz)
    deallocate(fxy,fxz)
    deallocate(fyx,fyz)
    deallocate(fzx,fzy)
    deallocate(Ip,Im)
    deallocate(Ip_g,Im_g)
    deallocate(Ip_gc,Im_gc)
      
  end subroutine umeth3d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi,          uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                     q,c,gamc,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                     src,srcQ,           src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                     courno,dx,dy,dz,dt,ngp,ngf)
    !
    !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
    !     if use_flattening=1.  Declared dimensions of q,c,gamc,csml,flatn are given
    !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
    !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
    !     routine that computes flatn).  
    !
    use network, only : nspec, naux
    use eos_module
    use eos_type_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UESGS, UTEMP, UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QW, &
                                   QREINT, QESGS, QPRES, QTEMP, QGAME, QFA, QFS, QFX, &
                                   nadv, allow_negative_energy, small_temp, use_flattening
    use flatten_module
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer lo(3), hi(3)
    integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    
    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    double precision :: c(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: gamc(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
    double precision ::  src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision :: dx, dy, dz, dt, courno

    double precision, allocatable:: dpdrho(:,:,:)
    double precision, allocatable:: dpde(:,:,:)
    double precision, allocatable:: dpdX_er(:,:,:,:)

    integer          :: i, j, k
    integer          :: pt_index(3)
    integer          :: ngp, ngf, loq(3), hiq(3)
    integer          :: n, nq
    integer          :: iadv, ispec, iaux
    double precision :: courx, coury, courz, courmx, courmy, courmz

    type (eos_t) :: eos_state

    allocate( dpdrho(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
    allocate(   dpde(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3))
    allocate(dpdX_er(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,nspec))

    do i=1,3
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo
    !
    ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
    ! The temperature is used as an initial guess for the eos call and will be overwritten.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = loq(3),hiq(3)
       do j = loq(2),hiq(2)
          do i = loq(1),hiq(1)
             
             if (uin(i,j,k,URHO) .le. ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                print *,'>>> ... negative density ',uin(i,j,k,URHO)
                call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
             end if

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             q(i,j,k,QU) = uin(i,j,k,UMX)/uin(i,j,k,URHO)
             q(i,j,k,QV) = uin(i,j,k,UMY)/uin(i,j,k,URHO)
             q(i,j,k,QW) = uin(i,j,k,UMZ)/uin(i,j,k,URHO)
             ! convert "rho e" to "e"
             q(i,j,k,QREINT ) = uin(i,j,k,UEINT)/q(i,j,k,QRHO)
             q(i,j,k,QTEMP  ) = uin(i,j,k,UTEMP)
             
             ! convert "rho K" to "K"
             if (QESGS .gt. -1) &
                  q(i,j,k,QESGS) = uin(i,j,k,UESGS)/q(i,j,k,QRHO)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
    !$OMP PARALLEL DO PRIVATE(iadv,n,nq,i,j,k) IF(nadv.gt.1)
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
    !$OMP END PARALLEL DO
      
    ! Load chemical species, c, into q, assuming they arrived in uin as rho.c
    !$OMP PARALLEL DO PRIVATE(ispec,n,nq,i,j,k) IF(nspec.gt.1)
    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       do k = loq(3),hiq(3)
          do j = loq(2),hiq(2)
             do i = loq(1),hiq(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
      
    ! Load auxiliary variables which are needed in the EOS
    !$OMP PARALLEL DO PRIVATE(iaux,n,nq,i,j,k) IF(naux.gt.1)
    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       do k = loq(3),hiq(3)
          do j = loq(2),hiq(2)
             do i = loq(1),hiq(1)
                q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! Get gamc, p, T, c, csml using q state
    !$OMP PARALLEL DO PRIVATE(i,j,k,eos_state,pt_index)
    do k = loq(3), hiq(3)
       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)
             
             pt_index(1) = i
             pt_index(2) = j
             pt_index(3) = k

             eos_state % T   = q(i,j,k,QTEMP)
             eos_state % rho = q(i,j,k,QRHO)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)
             eos_state % aux = q(i,j,k,QFX:QFX+naux-1)             

             ! If necessary, reset the energy using small_temp
             if ((allow_negative_energy .eq. 0) .and. (q(i,j,k,QREINT) .lt. ZERO)) then
                q(i,j,k,QTEMP) = small_temp
                eos_state % T =  q(i,j,k,QTEMP)

                call eos(eos_input_rt, eos_state, .false., pt_index = pt_index)
                q(i,j,k,QREINT) = eos_state % e

                if (q(i,j,k,QREINT) .lt. ZERO) then
                   print *,'   '
                   print *,'>>> Error: Castro_advection_3d::ctoprim ',i,j,k
                   print *,'>>> ... new e from eos (input_rt) call is negative ' &
                        ,q(i,j,k,QREINT)
                   print *,'    '
                   call bl_error("Error:: Castro_advection_3d.f90 :: ctoprim")
                end if
             end if

             eos_state % e = q(i,j,k,QREINT)

             call eos(eos_input_re, eos_state, .false., pt_index = pt_index)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e
             q(i,j,k,QPRES)  = eos_state % p

             dpdrho(i,j,k) = eos_state % dpdr_e
             dpde(i,j,k)   = eos_state % dpde
             c(i,j,k)      = eos_state % cs
             gamc(i,j,k)   = eos_state % gam1

             csml(i,j,k) = max(small, small * c(i,j,k))

             ! convert "e" back to "rho e"
             q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

             q(i,j,k,QGAME) = q(i,j,k,QPRES)/q(i,j,k,QREINT) + ONE

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! compute srcQ terms
    !$OMP PARALLEL DO PRIVATE(i,j,k,ispec,iaux,iadv)
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             
             srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
             srcQ(i,j,k,QU    ) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
             srcQ(i,j,k,QV    ) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
             srcQ(i,j,k,QW    ) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) / q(i,j,k,QRHO)
             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) &
                                                   - q(i,j,k,QV)*src(i,j,k,UMY) &
                                                   - q(i,j,k,QW)*src(i,j,k,UMZ) &
                                    + HALF * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2) * srcQ(i,j,k,QRHO)

             srcQ(i,j,k,QPRES ) = dpde(i,j,k)*(srcQ(i,j,k,QREINT) - &
                  q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)/q(i,j,k,QRHO)) /q(i,j,k,QRHO) + &
                  dpdrho(i,j,k)*srcQ(i,j,k,QRHO)! + &
!                                    sum(dpdX_er(i,j,k,:)*(src(i,j,k,UFS:UFS+nspec-1) - &
!                                                          q(i,j,k,QFS:QFS+nspec-1)*srcQ(i,j,k,QRHO))) &
!                                    /q(i,j,k,QRHO)

             if (QESGS .gt. -1) &
                  srcQ(i,j,k,QESGS) = src(i,j,k,UESGS)/q(i,j,k,QRHO) - q(i,j,k,QESGS) * srcQ(i,j,k,QRHO)

             do ispec = 1,nspec
                srcQ(i,j,k,QFS+ispec-1) = ( src(i,j,k,UFS+ispec-1) - q(i,j,k,QFS+ispec-1) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             enddo

             do iaux = 1,naux
                srcQ(i,j,k,QFX+iaux-1) = ( src(i,j,k,UFX+iaux-1) - q(i,j,k,QFX+iaux-1) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             enddo
             
             do iadv = 1,nadv
                srcQ(i,j,k,QFA+iadv-1) = ( src(i,j,k,UFA+iadv-1) - q(i,j,k,QFA+iadv-1) * srcQ(i,j,k,QRHO) ) / &
                     q(i,j,k,QRHO)
             enddo
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! Compute running max of Courant number over grids
    courmx = courno
    courmy = courno
    courmz = courno
    !$OMP PARALLEL DO PRIVATE(i,j,k,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dt/dx
             coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dt/dy
             courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dt/dz
             
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
    !$OMP END PARALLEL DO

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

    deallocate(dpdrho,dpde)
    
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
                    div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux)

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
    double precision dx, dy, dz, dt, E_added_flux

    double precision :: div1
    integer          :: i, j, k, n

    do n = 1, NVAR
         
       if ( n.eq.UTEMP ) then
          
          flux1(:,:,:,n) = ZERO
          flux2(:,:,:,n) = ZERO
          flux3(:,:,:,n) = ZERO
          
       else

          !$OMP PARALLEL PRIVATE(i,j,k,div1)
          !$OMP DO
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
          !$OMP END DO NOWAIT
          !$OMP DO
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
          !$OMP END DO NOWAIT
          !$OMP DO
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
          !$OMP END DO
          !$OMP END PARALLEL
          
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
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   uout(i,j,k,n) = uin(i,j,k,n) &
                          + ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) &
                          +   dt * src(i,j,k,n)
                   !
                   ! Add the source term to (rho e)
                   !
                   if (n .eq. UEINT) then
                      uout(i,j,k,n) = uout(i,j,k,n) - dt * pdivu(i,j,k)
                   else if (n .eq. UEDEN) then
                      E_added_flux = E_added_flux + &
                           ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                         +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                         +   flux3(i,j,k,n) - flux3(i,j,k+1,n)) / vol(i,j,k) 
                   endif
                enddo
             enddo
          enddo
          !$OMP END PARALLEL DO
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
    double precision :: ux, vy, wz

    !$OMP PARALLEL DO PRIVATE(i,j,k,ux,vy,wz)
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             
             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) )/ dx

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) )/ dy

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) )/ dz

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
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
    
    !$OMP PARALLEL PRIVATE(i,j,k,sum,n,fac)

    !$OMP DO
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
    !$OMP END DO NOWAIT
    
    !$OMP DO
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
    !$OMP END DO NOWAIT
    
    !$OMP DO
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
    !$OMP END DO
    
    !$OMP END PARALLEL

  end subroutine normalize_species_fluxes

! ::
! :: ----------------------------------------------------------
! ::

  subroutine enforce_minimum_density(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                                     uout,uout_l1,uout_l2,uout_l3, &
                                     uout_h1,uout_h2,uout_h3, &
                                     lo,hi,mass_added,eint_added,eden_added,verbose)
    
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, &
                                     UFA, small_dens, nadv
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3), verbose
    integer          ::  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
    integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
    double precision ::  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision :: mass_added, eint_added, eden_added
    
    ! Local variables
    integer          :: i,ii,j,jj,k,kk,n
    double precision :: min_dens
    double precision, allocatable :: fac(:,:,:)
    
    double precision :: initial_mass, final_mass
    double precision :: initial_eint, final_eint
    double precision :: initial_eden, final_eden
    
    allocate(fac(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    
    initial_mass = ZERO
      final_mass = ZERO

    initial_eint = ZERO
      final_eint = ZERO

    initial_eden = ZERO
      final_eden = ZERO

    min_dens = ZERO

    !$OMP PARALLEL DO PRIVATE(i,j,k,ii,jj,kk,min_dens) reduction(+:initial_mass,initial_eint,initial_eden)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             initial_mass = initial_mass + uout(i,j,k,URHO)
             initial_eint = initial_eint + uout(i,j,k,UEINT)
             initial_eden = initial_eden + uout(i,j,k,UEDEN)
             
             if (uout(i,j,k,URHO) .eq. ZERO) then
                
                print *,'DENSITY EXACTLY ZERO AT CELL ',i,j,k
                print *,'  in grid ',lo(1),lo(2),lo(3),hi(1),hi(2),hi(3)
                call bl_error("Error:: Castro_3d.f90 :: enforce_minimum_density")
                
             else if (uout(i,j,k,URHO) < small_dens) then
                
                min_dens = uin(i,j,k,URHO)
                do kk = -1,1
                   do jj = -1,1
                      do ii = -1,1
                         min_dens = min(min_dens,uin(i+ii,j+jj,k+kk,URHO))
                         if (i+ii.ge.uout_l1 .and. j+jj.ge.uout_l2 .and. k+kk.ge.uout_l3 .and. &
                             i+ii.le.uout_h1 .and. j+jj.le.uout_h2 .and. k+kk.le.uout_h3) then
                              if (uout(i+ii,j+jj,k+kk,URHO).gt.small_dens) &
                                min_dens = min(min_dens,uout(i+ii,j+jj,k+kk,URHO))
                         endif
                      end do
                   end do
                end do
                
                if (verbose .gt. 0) then
                   if (uout(i,j,k,URHO) < ZERO) then
                      print *,'   '
                      print *,'>>> RESETTING NEG.  DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',min_dens
                      print *,'   '
                   else
                      print *,'   '
                      print *,'>>> RESETTING SMALL DENSITY AT ',i,j,k
                      print *,'>>> FROM ',uout(i,j,k,URHO),' TO ',min_dens
                      print *,'   '
                   end if
                end if
                
                fac(i,j,k) = min_dens / uout(i,j,k,URHO)
                
             end if
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,n) reduction(+:final_mass,final_eint,final_eden)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             if (uout(i,j,k,URHO) < small_dens) then

                uout(i,j,k,URHO ) = uout(i,j,k,URHO ) * fac(i,j,k)
                uout(i,j,k,UEINT) = uout(i,j,k,UEINT) * fac(i,j,k)
                uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) * fac(i,j,k)
                uout(i,j,k,UMX  ) = uout(i,j,k,UMX  ) * fac(i,j,k)
                uout(i,j,k,UMY  ) = uout(i,j,k,UMY  ) * fac(i,j,k)
                uout(i,j,k,UMZ  ) = uout(i,j,k,UMZ  ) * fac(i,j,k)
   
                do n = UFS, UFS+nspec-1
                   uout(i,j,k,n) = uout(i,j,k,n) * fac(i,j,k)
                end do
                do n = UFX, UFX+naux-1
                   uout(i,j,k,n) = uout(i,j,k,n) * fac(i,j,k)
                end do
                do n = UFA, UFA+nadv-1
                   uout(i,j,k,n) = uout(i,j,k,n) * fac(i,j,k)
                end do
                
             end if
             
             final_mass = final_mass + uout(i,j,k,URHO)
             final_eint = final_eint + uout(i,j,k,UEINT)
             final_eden = final_eden + uout(i,j,k,UEDEN)
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! When enabled with OpenMP sometimes there is a small numerical error
    ! in (final_mass - initial_mass) even if no cells have been reset.
    ! Guard against this by only taking the difference if any cell has been reset.

    if ( min_dens /= ZERO ) then
       mass_added = mass_added + final_mass - initial_mass
       eint_added = eint_added + final_eint - initial_eint
       eden_added = eden_added + final_eden - initial_eden
    endif
    
    deallocate(fac)
    
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
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,sum,n,fac)
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
    !$OMP END PARALLEL DO
    
  end subroutine normalize_new_species

end module advection_module
