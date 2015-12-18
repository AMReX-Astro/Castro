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

  subroutine umeth3d(q, c, gamc, csml, flatn, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                     srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                     ilo1, ilo2, ilo3, ihi1, ihi2, ihi3, dx, dy, dz, dt, &
                     flux1, fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3, &
                     flux2, fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                     flux3, fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
                     rgdnvx_out,rgdx_l1,rgdx_l2,rgdx_l3,rgdx_h1,rgdx_h2,rgdx_h3, &
                     rgdnvy_out,rgdy_l1,rgdy_l2,rgdy_l3,rgdy_h1,rgdy_h2,rgdy_h3, &
                     rgdnvz_out,rgdz_l1,rgdz_l2,rgdz_l3,rgdz_h1,rgdz_h2,rgdz_h3, &
                     ugdnvx_out,ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3, &
                     ugdnvy_out,ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3, &
                     ugdnvz_out,ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3, &
                     vgdnvx_out,vgdx_l1,vgdx_l2,vgdx_l3,vgdx_h1,vgdx_h2,vgdx_h3, &
                     vgdnvy_out,vgdy_l1,vgdy_l2,vgdy_l3,vgdy_h1,vgdy_h2,vgdy_h3, &
                     vgdnvz_out,vgdz_l1,vgdz_l2,vgdz_l3,vgdz_h1,vgdz_h2,vgdz_h3, &
                     wgdnvx_out,wgdx_l1,wgdx_l2,wgdx_l3,wgdx_h1,wgdx_h2,wgdx_h3, &
                     wgdnvy_out,wgdy_l1,wgdy_l2,wgdy_l3,wgdy_h1,wgdy_h2,wgdy_h3, &
                     wgdnvz_out,wgdz_l1,wgdz_l2,wgdz_l3,wgdz_h1,wgdz_h2,wgdz_h3, &
                     pgdnvx_out,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                     pgdnvy_out,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                     pgdnvz_out,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                     pdivu, domlo, domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QW, QFS, QFX, QTEMP, QREINT, ppm_type, &
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

    integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
    integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    integer ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
    integer fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3
    integer fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3
    integer fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3
    integer rgdx_l1,rgdx_l2,rgdx_l3,rgdx_h1,rgdx_h2,rgdx_h3
    integer rgdy_l1,rgdy_l2,rgdy_l3,rgdy_h1,rgdy_h2,rgdy_h3
    integer rgdz_l1,rgdz_l2,rgdz_l3,rgdz_h1,rgdz_h2,rgdz_h3
    integer ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3
    integer ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3
    integer ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3
    integer vgdx_l1,vgdx_l2,vgdx_l3,vgdx_h1,vgdx_h2,vgdx_h3
    integer vgdy_l1,vgdy_l2,vgdy_l3,vgdy_h1,vgdy_h2,vgdy_h3
    integer vgdz_l1,vgdz_l2,vgdz_l3,vgdz_h1,vgdz_h2,vgdz_h3
    integer wgdx_l1,wgdx_l2,wgdx_l3,wgdx_h1,wgdx_h2,wgdx_h3
    integer wgdy_l1,wgdy_l2,wgdy_l3,wgdy_h1,wgdy_h2,wgdy_h3
    integer wgdz_l1,wgdz_l2,wgdz_l3,wgdz_h1,wgdz_h2,wgdz_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3    
    integer domlo(3),domhi(3)
    integer km,kc,kt,k3d,n
    integer i,j,iwave,idim
    
    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,NVAR)
    double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,NVAR)
    double precision flux3(fd3_l1:fd3_h1,fd3_l2:fd3_h2,fd3_l3:fd3_h3,NVAR)
    double precision rgdnvx_out(rgdx_l1:rgdx_h1,rgdx_l2:rgdx_h2,rgdx_l3:rgdx_h3)
    double precision rgdnvy_out(rgdy_l1:rgdy_h1,rgdy_l2:rgdy_h2,rgdy_l3:rgdy_h3)
    double precision rgdnvz_out(rgdz_l1:rgdz_h1,rgdz_l2:rgdz_h2,rgdz_l3:rgdz_h3)
    double precision ugdnvx_out(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ugdx_l3:ugdx_h3)
    double precision ugdnvy_out(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ugdy_l3:ugdy_h3)
    double precision ugdnvz_out(ugdz_l1:ugdz_h1,ugdz_l2:ugdz_h2,ugdz_l3:ugdz_h3)
    double precision vgdnvx_out(vgdx_l1:vgdx_h1,vgdx_l2:vgdx_h2,vgdx_l3:vgdx_h3)
    double precision vgdnvy_out(vgdy_l1:vgdy_h1,vgdy_l2:vgdy_h2,vgdy_l3:vgdy_h3)
    double precision vgdnvz_out(vgdz_l1:vgdz_h1,vgdz_l2:vgdz_h2,vgdz_l3:vgdz_h3)
    double precision wgdnvx_out(wgdx_l1:wgdx_h1,wgdx_l2:wgdx_h2,wgdx_l3:wgdx_h3)
    double precision wgdnvy_out(wgdy_l1:wgdy_h1,wgdy_l2:wgdy_h2,wgdy_l3:wgdy_h3)
    double precision wgdnvz_out(wgdz_l1:wgdz_h1,wgdz_l2:wgdz_h2,wgdz_l3:wgdz_h3)
    double precision pgdnvx_out(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision pgdnvy_out(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision pgdnvz_out(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)    
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
    
    double precision, pointer:: rgdnvx(:,:,:), pgdnvx(:,:,:), ugdnvx(:,:,:), &
                                vgdnvx(:,:,:), wgdnvx(:,:,:), gegdnvx(:,:,:)
    double precision, pointer:: rgdnvxf(:,:,:), pgdnvxf(:,:,:), ugdnvxf(:,:,:), &
                                vgdnvxf(:,:,:), wgdnvxf(:,:,:), gegdnvxf(:,:,:)
    double precision, pointer:: rgdnvtmpx(:,:,:), pgdnvtmpx(:,:,:), ugdnvtmpx(:,:,:), &
                                vgdnvtmpx(:,:,:), wgdnvtmpx(:,:,:), gegdnvtmpx(:,:,:)
    
    double precision, pointer:: rgdnvy(:,:,:), pgdnvy(:,:,:), ugdnvy(:,:,:), &
                                vgdnvy(:,:,:), wgdnvy(:,:,:), gegdnvy(:,:,:)
    double precision, pointer:: rgdnvyf(:,:,:), pgdnvyf(:,:,:), ugdnvyf(:,:,:), &
                                vgdnvyf(:,:,:), wgdnvyf(:,:,:), gegdnvyf(:,:,:)
    double precision, pointer:: rgdnvtmpy(:,:,:), pgdnvtmpy(:,:,:), ugdnvtmpy(:,:,:), &
                                vgdnvtmpy(:,:,:), wgdnvtmpy(:,:,:), gegdnvtmpy(:,:,:)
    
    double precision, pointer:: rgdnvz(:,:,:), pgdnvz(:,:,:), ugdnvz(:,:,:), &
                                vgdnvz(:,:,:), wgdnvz(:,:,:), gegdnvz(:,:,:)
    double precision, pointer:: rgdnvzf(:,:,:), pgdnvzf(:,:,:), ugdnvzf(:,:,:), &
                                vgdnvzf(:,:,:), wgdnvzf(:,:,:), gegdnvzf(:,:,:)
    double precision, pointer:: rgdnvtmpz1(:,:,:), pgdnvtmpz1(:,:,:), ugdnvtmpz1(:,:,:), &
                                vgdnvtmpz1(:,:,:), wgdnvtmpz1(:,:,:), gegdnvtmpz1(:,:,:)
    double precision, pointer:: rgdnvtmpz2(:,:,:), pgdnvtmpz2(:,:,:), ugdnvtmpz2(:,:,:), &
                                vgdnvtmpz2(:,:,:), wgdnvtmpz2(:,:,:), gegdnvtmpz2(:,:,:)
    
    double precision, pointer:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    double precision, pointer:: Ip_src(:,:,:,:,:,:), Im_src(:,:,:,:,:,:)
    double precision, pointer:: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    double precision, pointer :: shk(:,:,:)
    
    type (eos_t) :: eos_state

    call bl_allocate ( rgdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvxf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate (gegdnvtmpx, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( gegdnvy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvyf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    
    call bl_allocate ( rgdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( pgdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvtmpy, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvz, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvzf, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate (gegdnvtmpz1, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)

    call bl_allocate ( rgdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
    call bl_allocate ( pgdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( ugdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( vgdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)
    call bl_allocate ( wgdnvtmpz2, ilo1-1,ihi1+2,ilo2-1,ihi2+2,1,2)    
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
       
       ! for source terms
       call bl_allocate ( Ip_src, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,QVAR)
       call bl_allocate ( Im_src, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,QVAR)
       
       ! for gamc -- needed for the reference state in eigenvectors
       call bl_allocate ( Ip_gc, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,1)
       call bl_allocate ( Im_gc, ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,1)
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
                      q(:,:,:,QU:QW),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
          end do

          if (ppm_trace_sources .eq. 1) then
             do n=1,QVAR
                call ppm(srcQ(:,:,:,n),src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         q(:,:,:,QU:QW),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip_src(:,:,:,:,:,n),Im_src(:,:,:,:,:,n), &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
             enddo
          endif

          if (ppm_temp_fix /= 1) then
             call ppm(gamc(:,:,:),qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      q(:,:,:,QU:QW),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      Ip_gc(:,:,:,:,:,1),Im_gc(:,:,:,:,:,1), &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc)
          else          

             do iwave = 1, 3
                do idim = 1, 3

                   do j = ilo2-1, ihi2+1
                      do i = ilo1-1, ihi1+1
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

                   do j = ilo2-1, ihi2+1
                      do i = ilo1-1, ihi1+1
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
          call tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip,Im,Ip_src,Im_src,Ip_gc,Im_gc, &
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
                           srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
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
                   rgdnvx,ugdnvx,vgdnvx,wgdnvx,pgdnvx,gegdnvx, &
                   ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   1,ilo1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

       ! Compute \tilde{F}^y at kc (k3d)
       call cmpflx(qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                   rgdnvy,ugdnvy,vgdnvy,wgdnvy,pgdnvy,gegdnvy, &
                   ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
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
                   rgdnvtmpx,ugdnvtmpx,vgdnvtmpx,wgdnvtmpx,pgdnvtmpx,gegdnvtmpx, &
                   ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d,domlo,domhi)

       ! Compute F^{y|x} at kc (k3d)
       call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                   rgdnvtmpy,ugdnvtmpy,vgdnvtmpy,wgdnvtmpy,pgdnvtmpy,gegdnvtmpy, &
                   ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                   gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                   2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d,domlo,domhi)

       if (k3d.ge.ilo3) then
          
          ! Compute U_z at kc (k3d)
          if (ppm_type .gt. 0) then
             call tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             Ip,Im,Ip_src,Im_src,Ip_gc,Im_gc, &
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
                      rgdnvz, ugdnvz,vgdnvz,wgdnvz,pgdnvz,gegdnvz, &
                      ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
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
                      rgdnvtmpz1,ugdnvtmpz1,vgdnvtmpz1,wgdnvtmpz1,pgdnvtmpz1,gegdnvtmpz1, &
                      ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                      3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d,domlo,domhi)

          ! Compute F^{z|y} at kc (k3d)
          call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                      rgdnvtmpz2,ugdnvtmpz2,vgdnvtmpz2,wgdnvtmpz2,pgdnvtmpz2,gegdnvtmpz2, &
                      ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
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
                       srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                       hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d)

          ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
          call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      flux3,fd3_l1,fd3_l2,fd3_l3,fd3_h1,fd3_h2,fd3_h3, &
                      rgdnvzf,ugdnvzf,vgdnvzf,wgdnvzf,pgdnvzf,gegdnvzf, &
                      ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                      gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                      shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                      3,ilo1,ihi1,ilo2,ihi2,kc,k3d,k3d,domlo,domhi)

          do j=ilo2-1,ihi2+1
             do i=ilo1-1,ihi1+1
                rgdnvz_out(i,j,k3d) = rgdnvzf(i,j,kc)
                ugdnvz_out(i,j,k3d) = ugdnvzf(i,j,kc)
                vgdnvz_out(i,j,k3d) = vgdnvzf(i,j,kc)
                wgdnvz_out(i,j,k3d) = wgdnvzf(i,j,kc)
                pgdnvz_out(i,j,k3d) = pgdnvzf(i,j,kc)
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
                         rgdnvx,ugdnvx,vgdnvx,wgdnvx,pgdnvx,gegdnvx, &
                         ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1,domlo,domhi)

             ! Compute F^{y|z} at km (k3d-1)
             call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                         rgdnvy,ugdnvy,vgdnvy,wgdnvy,pgdnvy,gegdnvy, &
                         ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
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
                          srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                          hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1)

             ! Compute U''_y at km (k3d-1)
             call transxz(qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                          fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                          ugdnvx,pgdnvx,gegdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          ugdnvtmpz1,pgdnvtmpz1,gegdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          gamc,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                          srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                          hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1)

             ! Compute F^x at km (k3d-1)
             call cmpflx(qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         flux1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
                         rgdnvxf,ugdnvxf,vgdnvxf,wgdnvxf,pgdnvxf,gegdnvxf, &
                         ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         1,ilo1,ihi1+1,ilo2,ihi2,km,k3d-1,k3d-1,domlo,domhi)
             
             do j=ilo2-1,ihi2+1
                do i=ilo1-1,ihi1+2
                   rgdnvx_out(i,j,k3d-1) = rgdnvxf(i,j,km)
                   ugdnvx_out(i,j,k3d-1) = ugdnvxf(i,j,km)
                   vgdnvx_out(i,j,k3d-1) = vgdnvxf(i,j,km)
                   wgdnvx_out(i,j,k3d-1) = wgdnvxf(i,j,km)
                   pgdnvx_out(i,j,k3d-1) = pgdnvxf(i,j,km)
                end do
             end do
             
             ! Compute F^y at km (k3d-1)
             call cmpflx(qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         flux2,fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3, &
                         rgdnvyf,ugdnvyf,vgdnvyf,wgdnvyf,pgdnvyf,gegdnvyf, &
                         ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         shk,ilo1-1,ilo2-1,ilo3-1,ihi1+1,ihi2+1,ihi3+1, &
                         2,ilo1,ihi1,ilo2,ihi2+1,km,k3d-1,k3d-1,domlo,domhi)

             do j=ilo2-1,ihi2+2
                do i=ilo1-1,ihi1+1
                   rgdnvy_out(i,j,k3d-1) = rgdnvyf(i,j,km)
                   ugdnvy_out(i,j,k3d-1) = ugdnvyf(i,j,km)
                   vgdnvy_out(i,j,k3d-1) = vgdnvyf(i,j,km)
                   wgdnvy_out(i,j,k3d-1) = wgdnvyf(i,j,km)
                   pgdnvy_out(i,j,k3d-1) = pgdnvyf(i,j,km)
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
    call bl_deallocate ( rgdnvx)
    call bl_deallocate ( pgdnvx)
    call bl_deallocate ( ugdnvx)
    call bl_deallocate ( vgdnvx)
    call bl_deallocate ( wgdnvx)    
    call bl_deallocate (gegdnvx)

    call bl_deallocate ( rgdnvxf)
    call bl_deallocate ( pgdnvxf)
    call bl_deallocate ( ugdnvxf)
    call bl_deallocate ( vgdnvxf)
    call bl_deallocate ( wgdnvxf)
    call bl_deallocate (gegdnvxf)

    call bl_deallocate ( rgdnvtmpx)
    call bl_deallocate ( pgdnvtmpx)
    call bl_deallocate ( ugdnvtmpx)
    call bl_deallocate ( vgdnvtmpx)
    call bl_deallocate ( wgdnvtmpx)
    call bl_deallocate (gegdnvtmpx)
    
    call bl_deallocate ( rgdnvy)
    call bl_deallocate ( pgdnvy)
    call bl_deallocate ( ugdnvy)
    call bl_deallocate ( vgdnvy)
    call bl_deallocate ( wgdnvy)
    call bl_deallocate ( gegdnvy)

    call bl_deallocate ( rgdnvyf)
    call bl_deallocate ( pgdnvyf)
    call bl_deallocate ( ugdnvyf)
    call bl_deallocate ( vgdnvyf)
    call bl_deallocate ( wgdnvyf)
    call bl_deallocate (gegdnvyf)

    call bl_deallocate ( rgdnvtmpy)
    call bl_deallocate ( pgdnvtmpy)
    call bl_deallocate ( ugdnvtmpy)
    call bl_deallocate ( vgdnvtmpy)
    call bl_deallocate ( wgdnvtmpy)
    call bl_deallocate (gegdnvtmpy)

    call bl_deallocate ( rgdnvz)
    call bl_deallocate ( pgdnvz)
    call bl_deallocate ( ugdnvz)
    call bl_deallocate ( vgdnvz)
    call bl_deallocate ( wgdnvz)
    call bl_deallocate (gegdnvz)

    call bl_deallocate ( rgdnvzf)
    call bl_deallocate ( pgdnvzf)
    call bl_deallocate ( ugdnvzf)
    call bl_deallocate ( vgdnvzf)
    call bl_deallocate ( wgdnvzf)
    call bl_deallocate (gegdnvzf)

    call bl_deallocate ( rgdnvtmpz1)
    call bl_deallocate ( pgdnvtmpz1)
    call bl_deallocate ( ugdnvtmpz1)
    call bl_deallocate ( vgdnvtmpz1)
    call bl_deallocate ( wgdnvtmpz1)
    call bl_deallocate (gegdnvtmpz1)

    call bl_deallocate ( rgdnvtmpz2)
    call bl_deallocate ( pgdnvtmpz2)
    call bl_deallocate ( ugdnvtmpz2)
    call bl_deallocate ( vgdnvtmpz2)
    call bl_deallocate ( wgdnvtmpz2)
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
!    double precision, pointer:: dpdX_er(:,:,:,:)

    integer          :: i, j, k
    integer          :: ngp, ngf, loq(3), hiq(3)
    integer          :: n, nq, ipassive
    double precision :: courx, coury, courz, courmx, courmy, courmz
    double precision :: kineng, rhoinv
    double precision :: dtdx, dtdy, dtdz
    
    type (eos_t) :: eos_state

    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    do i=1,3
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo    
    
    call bl_allocate( dpdrho, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
    call bl_allocate(   dpde, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3)
!    call bl_allocate(dpdX_er, q_l1,q_h1,q_l2,q_h2,q_l3,q_h3,1,nspec)

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
                    q(q_l1,q_l2,q_l3,QPRES), &
                    q(q_l1,q_l2,q_l3,QU), &
                    q(q_l1,q_l2,q_l3,QV), &
                    q(q_l1,q_l2,q_l3,QW), &
                    flatn,(/ q_l1,q_l2,q_l3 /), (/ q_h1,q_h2,q_h3 /))
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

  subroutine consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                    uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                    src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                    flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                    flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                    flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                    rgdnvx,rgdx_l1,rgdx_l2,rgdx_l3,rgdx_h1,rgdx_h2,rgdx_h3, &
                    rgdnvy,rgdy_l1,rgdy_l2,rgdy_l3,rgdy_h1,rgdy_h2,rgdy_h3, &
                    rgdnvz,rgdz_l1,rgdz_l2,rgdz_l3,rgdz_h1,rgdz_h2,rgdz_h3, &
                    ugdnvx,ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3, &
                    ugdnvy,ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3, &
                    ugdnvz,ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3, &
                    vgdnvx,vgdx_l1,vgdx_l2,vgdx_l3,vgdx_h1,vgdx_h2,vgdx_h3, &
                    vgdnvy,vgdy_l1,vgdy_l2,vgdy_l3,vgdy_h1,vgdy_h2,vgdy_h3, &
                    vgdnvz,vgdz_l1,vgdz_l2,vgdz_l3,vgdz_h1,vgdz_h2,vgdz_h3, &
                    wgdnvx,wgdx_l1,wgdx_l2,wgdx_l3,wgdx_h1,wgdx_h2,wgdx_h3, &
                    wgdnvy,wgdy_l1,wgdy_l2,wgdy_l3,wgdy_h1,wgdy_h2,wgdy_h3, &
                    wgdnvz,wgdz_l1,wgdz_l2,wgdz_l3,wgdz_h1,wgdz_h2,wgdz_h3, &
                    pgdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                    pgdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                    pgdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                    area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
                    area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
                    area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
                    vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                    div,pdivu,lo,hi,dx,dy,dz,dt,E_added_flux,&
                    xmom_added_flux,ymom_added_flux,zmom_added_flux)

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
         UEDEN, UEINT, UTEMP, normalize_species, hybrid_hydro
    use bl_constants_module
    use castro_util_module, only : position
    use advection_util_module, only : normalize_species_fluxes

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
    integer rgdx_l1,rgdx_l2,rgdx_l3,rgdx_h1,rgdx_h2,rgdx_h3
    integer rgdy_l1,rgdy_l2,rgdy_l3,rgdy_h1,rgdy_h2,rgdy_h3
    integer rgdz_l1,rgdz_l2,rgdz_l3,rgdz_h1,rgdz_h2,rgdz_h3
    integer ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3
    integer ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3
    integer ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3
    integer vgdx_l1,vgdx_l2,vgdx_l3,vgdx_h1,vgdx_h2,vgdx_h3
    integer vgdy_l1,vgdy_l2,vgdy_l3,vgdy_h1,vgdy_h2,vgdy_h3
    integer vgdz_l1,vgdz_l2,vgdz_l3,vgdz_h1,vgdz_h2,vgdz_h3
    integer wgdx_l1,wgdx_l2,wgdx_l3,wgdx_h1,wgdx_h2,wgdx_h3
    integer wgdy_l1,wgdy_l2,wgdy_l3,wgdy_h1,wgdy_h2,wgdy_h3
    integer wgdz_l1,wgdz_l2,wgdz_l3,wgdz_h1,wgdz_h2,wgdz_h3    
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
    
    double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
    double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision   src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
    double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    double precision rgdnvx(rgdx_l1:rgdx_h1,rgdx_l2:rgdx_h2,rgdx_l3:rgdx_h3)
    double precision rgdnvy(rgdy_l1:rgdy_h1,rgdy_l2:rgdy_h2,rgdy_l3:rgdy_h3)
    double precision rgdnvz(rgdz_l1:rgdz_h1,rgdz_l2:rgdz_h2,rgdz_l3:rgdz_h3)
    double precision ugdnvx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ugdx_l3:ugdx_h3)
    double precision ugdnvy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ugdy_l3:ugdy_h3)
    double precision ugdnvz(ugdz_l1:ugdz_h1,ugdz_l2:ugdz_h2,ugdz_l3:ugdz_h3)
    double precision vgdnvx(vgdx_l1:vgdx_h1,vgdx_l2:vgdx_h2,vgdx_l3:vgdx_h3)
    double precision vgdnvy(vgdy_l1:vgdy_h1,vgdy_l2:vgdy_h2,vgdy_l3:vgdy_h3)
    double precision vgdnvz(vgdz_l1:vgdz_h1,vgdz_l2:vgdz_h2,vgdz_l3:vgdz_h3)
    double precision wgdnvx(wgdx_l1:wgdx_h1,wgdx_l2:wgdx_h2,wgdx_l3:wgdx_h3)
    double precision wgdnvy(wgdy_l1:wgdy_h1,wgdy_l2:wgdy_h2,wgdy_l3:wgdy_h3)
    double precision wgdnvz(wgdz_l1:wgdz_h1,wgdz_l2:wgdz_h2,wgdz_l3:wgdz_h3)    
    double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)    
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
    double precision area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
    double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision dx, dy, dz, dt
    double precision E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux

    double precision ang_mom_flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3)
    double precision ang_mom_flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3)
    double precision ang_mom_flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3)

    double precision rad_mom_flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3)
    double precision rad_mom_flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3)
    double precision rad_mom_flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3)    

    double precision ang_mom, rad_mom
    
    double precision :: div1, volinv
    double precision :: rho, u, v, w, p
    integer          :: i, j, k, n
    double precision :: loc(3), R
    
    if (hybrid_hydro .eq. 1) then
       
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1

                loc = position(i,j,k,ccx=.false.)

                R = sqrt( loc(1)**2 + loc(2)**2 )
                
                rho = rgdnvx(i,j,k)
                u   = ugdnvx(i,j,k)
                v   = vgdnvx(i,j,k)
                w   = wgdnvx(i,j,k)
                p   = pgdnvx(i,j,k)

                rad_mom = (rho / R) * (loc(1) * u + loc(2) * v)
                ang_mom = rho * (loc(1) * v - loc(2) * u)
                
                rad_mom_flux1(i,j,k) = (rad_mom * u             ) * area1(i,j,k) * dt
                ang_mom_flux1(i,j,k) = (ang_mom * u + loc(2) * p) * area1(i,j,k) * dt

             enddo
          enddo
       enddo

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                
                loc = position(i,j,k,ccy=.false.)

                R = sqrt( loc(1)**2 + loc(2)**2 )
                
                rho = rgdnvy(i,j,k)
                u   = vgdnvy(i,j,k)
                v   = ugdnvy(i,j,k)
                w   = wgdnvy(i,j,k)
                p   = pgdnvy(i,j,k)

                rad_mom = (rho / R) * (loc(1) * u + loc(2) * v)
                ang_mom = rho * (loc(1) * v - loc(2) * u)
                
                rad_mom_flux2(i,j,k) = (rad_mom * v             ) * area2(i,j,k) * dt
                ang_mom_flux2(i,j,k) = (ang_mom * v - loc(1) * p) * area2(i,j,k) * dt
                
             enddo
          enddo
       enddo
       
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                
                loc = position(i,j,k,ccz=.false.)

                R = sqrt( loc(1)**2 + loc(2)**2 )                
                
                rho = rgdnvz(i,j,k)
                u   = vgdnvz(i,j,k)
                v   = wgdnvz(i,j,k)
                w   = ugdnvz(i,j,k)
                p   = pgdnvz(i,j,k)

                rad_mom = (rho / R) * (loc(1) * u + loc(2) * v)
                ang_mom = rho * (loc(1) * v - loc(2) * u)                
                
                rad_mom_flux3(i,j,k) = (rad_mom * w) * area3(i,j,k) * dt
                ang_mom_flux3(i,j,k) = (ang_mom * w) * area3(i,j,k) * dt

             enddo
          enddo
       enddo

    endif
       
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

                   ! Add up some diagnostic quantities.
                      
                   if (n .eq. UMX) then
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

    ! Now update the radial and angular momenta, and overwrite the
    ! x- and y- momenta accordingly.
    
    if (hybrid_hydro .eq. 1) then
       
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                
                loc = position(i,j,k)
                
                R = sqrt( loc(1)**2 + loc(2)**2 )
                
                rho = uin(i,j,k,URHO)

                volinv = ONE / vol(i,j,k)
                
                rad_mom = uin(i,j,k,UMX) * (loc(1) / R) + uin(i,j,k,UMY) * (loc(2) / R)
                ang_mom = uin(i,j,k,UMY) * loc(1)       - uin(i,j,k,UMX) * loc(2)

                rad_mom = rad_mom &
                        + ( rad_mom_flux1(i,j,k) - rad_mom_flux1(i+1,j,k) &
                        +   rad_mom_flux2(i,j,k) - rad_mom_flux2(i,j+1,k) &
                        +   rad_mom_flux3(i,j,k) - rad_mom_flux3(i,j,k+1) ) * volinv &
                        + dt * ( - (loc(1) / R) * (pgdnvx(i+1,j,k) - pgdnvx(i,j,k)) / dx &
                                 - (loc(2) / R) * (pgdnvy(i,j+1,k) - pgdnvy(i,j,k)) / dy &
                                 + ang_mom**2 / (rho * R**3) )                
                
                ang_mom = ang_mom &
                        + ( ang_mom_flux1(i,j,k) - ang_mom_flux1(i+1,j,k) &
                        +   ang_mom_flux2(i,j,k) - ang_mom_flux2(i,j+1,k) &
                        +   ang_mom_flux3(i,j,k) - ang_mom_flux3(i,j,k+1) ) * volinv

                uout(i,j,k,UMX) = rad_mom * (loc(1) / R)    - ang_mom * (loc(2) / R**2)
                uout(i,j,k,UMY) = ang_mom * (loc(1) / R**2) + rad_mom * (loc(2) / R)
                
             enddo
          enddo
       enddo
       
    endif
    
    
  end subroutine consup

end module advection_module
