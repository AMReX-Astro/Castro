module advection_module

  implicit none

  private

  public umeth2d, ctoprim, divu, consup, enforce_minimum_density, &
       normalize_new_species, &
       normalize_species_fluxes, uflaten
  
contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH2D     Compute hyperbolic fluxes using unsplit second
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
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: ----------------------------------------------------------------

  subroutine umeth2d(q, c, gamc, csml, flatn, qd_l1, qd_l2, qd_h1, qd_h2,&
                     srcQ, src_l1, src_l2, src_h1, src_h2, &
                     grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                     ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
                     flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                     flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                     pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                     pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                     ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                     ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                     area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                     area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                     pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                     dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                     domlo, domhi)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, ppm_type
    use trace_module, only : trace
    use trace_ppm_module, only : trace_ppm
    use riemann_module, only: cmpflx

    implicit none

    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer dloga_l1, dloga_l2, dloga_h1, dloga_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer gv_l1, gv_l2, gv_h1, gv_h2
    integer fd1_l1, fd1_l2, fd1_h1, fd1_h2
    integer fd2_l1, fd2_l2, fd2_h1, fd2_h2
    integer pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2
    integer pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2
    integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
    integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
    integer area1_l1, area1_l2, area1_h1, area1_h2
    integer area2_l1, area2_l2, area2_h1, area2_h2
    integer vol_l1, vol_l2, vol_h1, vol_h2
    integer ilo1, ilo2, ihi1, ihi2
    integer domlo(2), domhi(2)

    double precision dx, dy, dt
    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision  srcQ(src_l1:src_h1,src_l2:src_h2)
    double precision  grav( gv_l1: gv_h1, gv_l2: gv_h2)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    double precision pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
    double precision pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
    double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
    double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
    double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
    double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
    double precision pdivu(ilo1:ihi1,ilo2:ihi2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dq(:,:,:),  qm(:,:,:),   qp(:,:,:)
    double precision, allocatable::qxm(:,:,:),qym(:,:,:)
    double precision, allocatable::qxp(:,:,:),qyp(:,:,:)
    
    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    double precision, allocatable::   fx(:,:,:),  fy(:,:,:)
    double precision, allocatable::   pgdxtmp(:,:) ,  ugdxtmp(:,:)

    ! Local scalar variables
    double precision :: dtdx
    double precision :: hdtdx, hdt, hdtdy
    integer          :: i,j

    allocate ( pgdxtmp(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2))
    allocate ( ugdxtmp(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2))
    allocate ( dq(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( fx(ilo1:ihi1+1,ilo2-1:ihi2+1,NVAR))
    allocate ( fy(ilo1-1:ihi1+1,ilo2:ihi2+1,NVAR))

    ! Local constants
    dtdx = dt/dx
    hdtdx = 0.5d0*dtdx
    hdtdy = 0.5d0*dt/dy
    hdt = 0.5d0*dt
    
    ! NOTE: Geometry terms need to be punched through

    ! Trace to edges w/o transverse flux correction terms
    if (ppm_type .eq. 0) then
       call trace(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                  dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                  dq,qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                  grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                  ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    else
       call trace_ppm(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                      qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                      gamc,qd_l1,qd_l2,qd_h1,qd_h2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    end if

    call cmpflx(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                pgdxtmp, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                ugdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                1, ilo1, ihi1, ilo2-1, ihi2+1, domlo, domhi)

    call cmpflx(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                2, ilo1-1, ihi1+1, ilo2, ihi2, domlo, domhi)

    call transy(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ, src_l1, src_l2, src_h1, src_h2, &
                grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                hdt, hdtdy, &
                ilo1-1, ihi1+1, ilo2, ihi2)
    
    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                1, ilo1, ihi1, ilo2, ihi2, domlo, domhi)
      
    call transx(qym, qm,qyp,qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                pgdxtmp, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                ugdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
                grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                hdt, hdtdx, &
                area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                ilo1, ihi1, ilo2-1, ihi2+1)

    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                2, ilo1, ihi1, ilo2, ihi2, domlo, domhi)
      

    do j = ilo2,ihi2
       do i = ilo1,ihi1
          pdivu(i,j) = 0.5d0 * &
               ((pgdx(i+1,j)+pgdx(i,j))*(ugdx(i+1,j)*area1(i+1,j)-ugdx(i,j)*area1(i,j)) &
               +(pgdy(i,j+1)+pgdy(i,j))*(ugdy(i,j+1)*area2(i,j+1)-ugdy(i,j)*area2(i,j)) ) / vol(i,j)
       end do
    end do

    deallocate(dq,qm,qp,qxm,qxp,qym,qyp)
    deallocate(fx,fy)
    deallocate(pgdxtmp,ugdxtmp)
    
  end subroutine umeth2d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi, &
                     uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                     q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
                     src,srcQ,src_l1,src_l2,src_h1,src_h2, &
                     courno,dx,dy,dt,ngp,ngf,iflaten)
    
    ! Will give primitive variables on lo-ngp:hi+ngp, and flatn on
    ! lo-ngf:hi+ngf if iflaten=1.  Declared dimensions of
    ! q,c,gamc,csml,flatn are given by DIMS(q).  This declared region
    ! is assumed to encompass lo-ngp:hi+ngp.  Also, uflaten call
    ! assumes ngp>=ngf+3 (ie, primitve data is used by the routine
    ! that computes flatn).

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP,&
                                   UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QREINT, QPRES, QTEMP, &
                                   QFA, QFS, QFX, &
                                   nadv, allow_negative_energy, small_temp

    implicit none
    
    double precision, parameter:: small = 1.d-8

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer q_l1,q_l2,q_h1,q_h2
    integer src_l1,src_l2,src_h1,src_h2
    integer iflaten
    
    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    double precision :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
    double precision :: c(q_l1:q_h1,q_l2:q_h2)
    double precision :: gamc(q_l1:q_h1,q_l2:q_h2)
    double precision :: csml(q_l1:q_h1,q_l2:q_h2)
    double precision :: flatn(q_l1:q_h1,q_l2:q_h2)
    double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
    double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision :: dx, dy, dt, courno
    
    double precision, allocatable :: dpdrho(:,:)
    double precision, allocatable :: dpde(:,:)
    double precision, allocatable :: dpdX_er(:,:,:)

    integer          :: i, j
    integer          :: pt_index(2)
    integer          :: ngp, ngf, loq(2), hiq(2)
    integer          :: iadv, ispec, iaux, n, nq
    double precision :: courx, coury, courmx, courmy
    
    allocate(     dpdrho(q_l1:q_h1,q_l2:q_h2))
    allocate(     dpde(q_l1:q_h1,q_l2:q_h2))
    allocate(  dpdX_er(q_l1:q_h1,q_l2:q_h2,nspec))
    
    do i=1,2
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo
    
    ! Make q (all but p), except put e in slot for rho.e, fix after
    ! eos call The temperature is used as an initial guess for the eos
    ! call and will be overwritten
    do j = loq(2),hiq(2)
       do i = loq(1),hiq(1)

          if (uin(i,j,URHO) .le. 0.d0) then
             print *,'   '
             print *,'>>> Error: Castro_2d::ctoprim ',i,j
             print *,'>>> ... negative density ',uin(i,j,URHO)
             print *,'    '
             call bl_error("Error:: Castro_2d.f90 :: ctoprim")
          end if
          
          q(i,j,QRHO) = uin(i,j,URHO)
          q(i,j,QU) = uin(i,j,UMX)/uin(i,j,URHO)
          q(i,j,QV) = uin(i,j,UMY)/uin(i,j,URHO)
          q(i,j,QREINT ) = uin(i,j,UEINT)/q(i,j,QRHO)
          q(i,j,QTEMP  ) = uin(i,j,UTEMP)
       enddo
    enddo
    
    ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
    do iadv = 1, nadv
       n  = UFA + iadv - 1
       nq = QFA + iadv - 1
       do j = loq(2),hiq(2)
          do i = loq(1),hiq(1)
             q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
          enddo
       enddo
    enddo
    
    ! Load chemical species, c, into q, assuming they arrived in uin as rho.c
    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       do j = loq(2),hiq(2)
          do i = loq(1),hiq(1)
             q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
          enddo
       enddo
    enddo
    
    ! Load auxiliary variables which are needed in the EOS
    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       do j = loq(2),hiq(2)
          do i = loq(1),hiq(1)
             q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
          enddo
       enddo
    enddo

    ! Get gamc, p, T, c, csml using q state 
    do j = loq(2), hiq(2)
       do i = loq(1), hiq(1)

          ! If necessary, reset the energy using small_temp
          if ((allow_negative_energy .eq. 0) .and. (q(i,j,QREINT) .lt. 0)) then
             q(i,j,QTEMP) = small_temp
             call eos_given_RTX(q(i,j,QREINT),q(i,j,QPRES),q(i,j,QRHO), &
                                q(i,j,QTEMP),q(i,j,QFS:))
             if (q(i,j,QREINT) .lt. 0.d0) then
                print *,'   '
                print *,'>>> Error: Castro_2d::ctoprim ',i,j
                print *,'>>> ... new e from eos_given_RTX call is negative ',q(i,j,QREINT)
                print *,'    '
                call bl_error("Error:: Castro_2d.f90 :: ctoprim")
             end if
          end if
          
          pt_index(1) = i
          pt_index(2) = j
          call eos_given_ReX(gamc(i,j), q(i,j,QPRES), c(i,j), q(i,j,QTEMP), &
                             dpdrho(i,j), dpde(i,j), &
                             q(i,j,QRHO), q(i,j,QREINT), q(i,j,QFS:), &
                             pt_index=pt_index)!, &
          !                              dpdX_er=dpdX_er(i,j,:))
          csml(i,j) = max(small, small * c(i,j))
       end do
    end do

    ! Make this "rho e" instead of "e"
    do j = loq(2),hiq(2)
       do i = loq(1),hiq(1)
          q(i,j,QREINT) = q(i,j,QREINT)*q(i,j,QRHO)
       enddo
    enddo
    
    ! Compute sources in terms of Q
    do j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
          
          srcQ(i,j,QRHO  ) = src(i,j,URHO)
          srcQ(i,j,QU    ) = (src(i,j,UMX) - q(i,j,QU) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
          srcQ(i,j,QV    ) = (src(i,j,UMY) - q(i,j,QV) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
          ! S_rhoe = S_rhoE - u . (S_rhoU - 0.5 u S_rho)
          srcQ(i,j,QREINT) = src(i,j,UEDEN) - q(i,j,QU) *src(i,j,UMX)   &
                                            - q(i,j,QV) *src(i,j,UMY) + &
               0.5d0 * (q(i,j,QU)**2 + q(i,j,QV)**2) * srcQ(i,j,QRHO)
          srcQ(i,j,QPRES ) = dpde(i,j) * &
               (srcQ(i,j,QREINT) - q(i,j,QREINT)*srcQ(i,j,QRHO)/q(i,j,QRHO))/q(i,j,QRHO) + &
               dpdrho(i,j) * srcQ(i,j,QRHO)! + &
!                sum(dpdX_er(i,j,:)*(src(i,j,UFS:UFS+nspec-1) - &
!                    q(i,j,QFS:QFS+nspec-1)*srcQ(i,j,QRHO))) / q(i,j,QRHO)

          do ispec = 1,nspec
             srcQ(i,j,QFS+ispec-1) = ( src(i,j,UFS+ispec-1) - q(i,j,QFS+ispec-1) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
          enddo

          do iaux = 1,naux
             srcQ(i,j,QFX+iaux-1) = ( src(i,j,UFX+iaux-1) - q(i,j,QFX+iaux-1) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
          enddo
          
          do iadv = 1,nadv
             srcQ(i,j,QFA+iadv-1) = ( src(i,j,UFA+iadv-1) - q(i,j,QFA+iadv-1) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
          enddo
          
       end do
    end do

    ! Compute running max of Courant number over grids
    courmx = courno
    courmy = courno
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          courx =  ( c(i,j)+abs(q(i,j,QU)) ) * dt/dx
          coury =  ( c(i,j)+abs(q(i,j,QV)) ) * dt/dy
          courmx = max( courmx, courx )
          courmy = max( courmy, coury )
          
          if (courx .gt. 1.d0) then
             print *,'   '
             call bl_warning("Warning:: Castro_2d.f90 :: CFL violation in ctoprim")
             print *,'>>> ... (u+c) * dt / dx > 1 ', courx
             print *,'>>> ... at cell (i,j)     : ',i,j
             print *,'>>> ... u, c                ',q(i,j,QU), c(i,j)
             print *,'>>> ... density             ',q(i,j,QRHO)
          end if
          
          if (coury .gt. 1.d0) then
             print *,'   '
             call bl_warning("Warning:: Castro_2d.f90 :: CFL violation in ctoprim")
             print *,'>>> ... (v+c) * dt / dx > 1 ', coury
             print *,'>>> ... at cell (i,j)     : ',i,j
             print *,'>>> ... v, c                ',q(i,j,QV), c(i,j)
             print *,'>>> ... density             ',q(i,j,QRHO)
          end if
          
       enddo
    enddo
    courno = max( courmx, courmy )
      
    ! Compute flattening coef for slope calculations
    if(iflaten.eq.1)then
       do n=1,2
          loq(n)=lo(n)-ngf
          hiq(n)=hi(n)+ngf
       enddo
       call uflaten(loq,hiq, &
            q(q_l1,q_l2,QPRES), &
            q(q_l1,q_l2,QU), &
            q(q_l1,q_l2,QV), &
            flatn,q_l1,q_l2,q_h1,q_h2)
    else
       flatn = 1.d0
    endif
    
    deallocate(dpdrho,dpde)
    
  end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine consup( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                     uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                     pgdx,pgdx_l1,pgdx_l2,pgdx_h1,pgdx_h2, &
                     pgdy,pgdy_l1,pgdy_l2,pgdy_h1,pgdy_h2, &
                     src , src_l1, src_l2, src_h1, src_h2, &
                     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                     area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                     area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                     vol,vol_l1,vol_l2,vol_h1,vol_h2, &
                     div,pdivu,lo,hi,dx,dy,dt,E_added_flux)

    use eos_module
    use network, only : nspec, naux
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, &
                                   UEDEN, UEINT, UTEMP, &
                                   normalize_species

    implicit none

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer uout_l1,uout_l2,uout_h1,uout_h2
    integer pgdx_l1,pgdx_l2,pgdx_h1,pgdx_h2
    integer pgdy_l1,pgdy_l2,pgdy_h1,pgdy_h2
    integer   src_l1,  src_l2,  src_h1,  src_h2
    integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
    integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
    integer area1_l1,area1_l2,area1_h1,area1_h2
    integer area2_l1,area2_l2,area2_h1,area2_h2
    integer vol_l1,vol_l2,vol_h1,vol_h2
    
    double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
    double precision pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
    double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,NVAR)
    double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
    double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
    double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2))
    double precision dx, dy, dt, E_added_flux
    
    integer i, j, n

    double precision div1
    !double precision SrU, SrV
    !double precision rho, Up, Vp, SrE

    ! Normalize the species fluxes
    if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                lo,hi)

    ! correct the fluxes to include the effects of the artificial viscosity
    do n = 1, NVAR
       if ( n.eq.UTEMP) then
          flux1(:,:,n) = 0.d0
          flux2(:,:,n) = 0.d0
       else 
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                div1 = .5d0*(div(i,j) + div(i,j+1))
                div1 = difmag*min(0.d0,div1)
                flux1(i,j,n) = flux1(i,j,n) &
                     + dx*div1*(uin(i,j,n) - uin(i-1,j,n))
                flux1(i,j,n) = area1(i,j)*flux1(i,j,n)
             enddo
          enddo
          
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                div1 = .5d0*(div(i,j) + div(i+1,j))
                div1 = difmag*min(0.d0,div1)
                flux2(i,j,n) = flux2(i,j,n) &
                     + dy*div1*(uin(i,j,n) - uin(i,j-1,n))
                flux2(i,j,n) = area2(i,j)*flux2(i,j,n)
             enddo
          enddo
       endif
    enddo
    
    ! do the conservative updates
    do n = 1, NVAR
       if (n .eq. UTEMP) then
          uout(lo(1):hi(1),lo(2):hi(2),n) = uin(lo(1):hi(1),lo(2):hi(2),n)
       else 
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                uout(i,j,n) = uin(i,j,n) + dt * &
                     ( flux1(i,j,n) - flux1(i+1,j,n) &
                     +   flux2(i,j,n) - flux2(i,j+1,n) ) / vol(i,j) &
                     +   dt * src(i,j,n)
                
                if (n .eq. UEINT) then
                   ! Add source term to (rho e)
                   uout(i,j,UEINT) = uout(i,j,UEINT)  - dt * pdivu(i,j)
                else if (n .eq. UEDEN) then
                   E_added_flux = E_added_flux + dt * & 
                        ( flux1(i,j,n) - flux1(i+1,j,n) &
                        +   flux2(i,j,n) - flux2(i,j+1,n) ) / vol(i,j) 
                   
                end if
             enddo
          enddo
       end if
    enddo

    ! Add gradp term to momentum equation
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
!         uout(i,j,UMX) = uout(i,j,UMX)+ 0.5d0*(area1(i,j)+area1(i+1,j))* &
!            dt * ( pgdx(i,j)-pgdx(i+1,j) )/vol(i,j)
!         uout(i,j,UMY) = uout(i,j,UMY)+ 0.5d0*(area2(i,j)+area2(i,j+1))* &
!            dt * ( pgdy(i,j)-pgdy(i,j+1) )/vol(i,j)

          uout(i,j,UMX) = uout(i,j,UMX) - dt * (pgdx(i+1,j)-pgdx(i,j))/ dx
          uout(i,j,UMY) = uout(i,j,UMY) - dt * (pgdy(i,j+1)-pgdy(i,j))/ dy
       enddo
    enddo

    ! scale the fluxes (and correct the momentum flux with the grad p part)
    ! so we can use them in the flux correction at coarse-fine interfaces
    ! later.
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR)
          flux1(i,j,   UMX) = flux1(i,j,UMX) + dt*area1(i,j)*pgdx(i,j)
       enddo
    enddo
    
    do j = lo(2),hi(2)+1 
       do i = lo(1),hi(1)
          flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR)
          flux2(i,j,UMY) = flux2(i,j,UMY) + dt*area2(i,j)*pgdy(i,j)
       enddo
    enddo
    
  end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine uflaten(lo,hi,p,u,v,flatn, &
                     q_l1,q_l2,q_h1,q_h2)

    use meth_params_module, only : iorder, small_pres
    
    implicit none

    integer lo(2),hi(2)
    integer q_l1,q_l2,q_h1,q_h2
    double precision p(q_l1:q_h1,q_l2:q_h2)
    double precision u(q_l1:q_h1,q_l2:q_h2)
    double precision v(q_l1:q_h1,q_l2:q_h2)
    double precision flatn(q_l1:q_h1,q_l2:q_h2)
    
    ! Local arrays
    double precision, allocatable :: dp(:), z(:), chi(:)
    
    integer i, j, idx, ishft
    double precision shktst, zcut1, zcut2, dzcut
    double precision denom, zeta, tst, tmp, ftmp
    integer nx,ny,nmax
    
    ! Knobs for detection of strong shock
    data shktst /0.33d0/
    data zcut1 /0.75d0/
    data zcut2 /0.85d0/
    
    nx = hi(1)-lo(1)+3
    ny = hi(2)-lo(2)+3
    nmax = max(nx,ny)
    allocate(dp(0:nmax-1),z(0:nmax-1),chi(0:nmax-1))
    
    dzcut = 1.d0/(zcut2-zcut1)
    
    if (iorder .eq. 3) then
       do j = lo(2),hi(2) 
          do i = lo(1),hi(1) 
             flatn(i,j) = 1.d0
          enddo
       enddo
       return
    endif
    
    ! x-direction flattening coef
    do j = lo(2),hi(2) 
       do i = lo(1)-1,hi(1)+1
          idx = i-lo(1)+1
          dp(idx) = p(i+1,j) - p(i-1,j)
          denom = max(small_pres,abs(p(i+2,j)-p(i-2,j)))
          zeta = abs(dp(idx))/denom
          z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
          if (u(i-1,j)-u(i+1,j) .ge. 0.d0) then
             tst = 1.d0
          else
             tst = 0.d0
          endif
          tmp = min(p(i+1,j),p(i-1,j))
          if ((abs(dp(idx))/tmp).gt.shktst) then
             chi(idx) = tst
          else
             chi(idx) = 0.d0
          endif
       enddo
       do i = lo(1),hi(1)
          idx = i-lo(1)+1
          if(dp(idx).gt.0.d0)then
             ishft = 1
          else
             ishft = -1
          endif
          flatn(i,j) = 1.d0 - &
               max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
       enddo
    enddo
    
    ! y-direction flattening coef
    do i = lo(1),hi(1)
       do j = lo(2)-1,hi(2)+1
          idx = j-lo(2)+1
          dp(idx) = p(i,j+1) - p(i,j-1)
          denom = max(small_pres,abs(p(i,j+2)-p(i,j-2)))
          zeta = abs(dp(idx))/denom
          z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
          if (v(i,j-1)-v(i,j+1) .ge. 0.d0) then
             tst = 1.d0
          else
             tst = 0.d0
          endif
          tmp = min(p(i,j+1),p(i,j-1))
          if ((abs(dp(idx))/tmp).gt.shktst) then
             chi(idx) = tst
          else
             chi(idx) = 0.d0
          endif
       enddo
       do j = lo(2),hi(2)
          idx = j-lo(2)+1
          if(dp(idx).gt.0.d0)then
             ishft = 1
          else
             ishft = -1
          endif
          ftmp = 1.d0 - &
               max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
          flatn(i,j) = min( flatn(i,j), ftmp )
       enddo
    enddo
    
    deallocate(dp,z,chi)
    
  end subroutine uflaten

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,dx, &
                  div,div_l1,div_l2,div_h1,div_h2)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV
    
    implicit none
    
    integer          :: lo(2),hi(2)
    integer          :: q_l1,q_l2,q_h1,q_h2
    integer          :: div_l1,div_l2,div_h1,div_h2
    double precision :: q(q_l1:q_h1,q_l2:q_h2,*)
    double precision :: div(div_l1:div_h1,div_l2:div_h2)
    double precision :: dx(2)
    
    integer          :: i, j
    double precision :: rl, rr, rc, ul, ur
    double precision :: vb, vt
    double precision :: ux,vy
    
    if (coord_type .eq. 0) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             ux = 0.5d0*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/dx(1)
             vy = 0.5d0*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/dx(2)
             div(i,j) = ux + vy
          enddo
       enddo
    else
       do i=lo(1),hi(1)+1
          
          if (i.eq.0) then
             
             div(i,lo(2):hi(2)+1) = 0.d0
             
          else 

             rl = (dble(i)-0.5d0) * dx(1)
             rr = (dble(i)+0.5d0) * dx(1)
             rc = (dble(i)      ) * dx(1)
             
             do j=lo(2),hi(2)+1
                ! These are transverse averages in the y-direction
                ul = 0.5d0 * (q(i-1,j,QU)+q(i-1,j-1,QU))
                ur = 0.5d0 * (q(i  ,j,QU)+q(i  ,j-1,QU))
                
                ! Take 1/r d/dr(r*u)
                div(i,j) = (rr*ur - rl*ul) / dx(1) / rc
                
                ! These are transverse averages in the x-direction
                vb = 0.5d0 * (q(i,j-1,QV)+q(i-1,j-1,QV))
                vt = 0.5d0 * (q(i,j  ,QV)+q(i-1,j  ,QV))
                
                div(i,j) = div(i,j) + (vt - vb) / dx(2)
             enddo
             
          end if
       enddo
    end if

  end subroutine divu

! ::
! :: ----------------------------------------------------------
! ::

  subroutine normalize_species_fluxes(  &
                    flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                    flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                    lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    
    implicit none

    integer          :: lo(2),hi(2)
    integer          :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
    integer          :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    double precision :: sum,fac
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          sum = 0.d0
          do n = UFS, UFS+nspec-1
             sum = sum + flux1(i,j,n)
          end do
          if (sum .ne. 0.d0) then
             fac = flux1(i,j,URHO) / sum
          else
             fac = 1.d0
          end if
          do n = UFS, UFS+nspec-1
             flux1(i,j,n) = flux1(i,j,n) * fac
          end do
       end do
    end do
    do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          sum = 0.d0
          do n = UFS, UFS+nspec-1
             sum = sum + flux2(i,j,n)
          end do
          if (sum .ne. 0.d0) then
             fac = flux2(i,j,URHO) / sum
          else
             fac = 1.d0
          end if
          do n = UFS, UFS+nspec-1
             flux2(i,j,n) = flux2(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_species_fluxes

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine enforce_minimum_density( uin,  uin_l1, uin_l2, uin_h1, uin_h2, &
                                      uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                                      lo, hi, mass_added, eint_added, eden_added, verbose)
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UEDEN, &
                                   UFS, UFX, UFA, small_dens, nadv

    implicit none

    integer          :: lo(2), hi(2), verbose
    integer          :: uin_l1,uin_l2,uin_h1,uin_h2
    integer          :: uout_l1,uout_l2,uout_h1,uout_h2
    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision :: mass_added, eint_added, eden_added
    
    ! Local variables
    integer                       :: i,ii,j,jj,n
    double precision              :: min_dens
    double precision, allocatable :: fac(:,:)
    double precision              :: initial_mass, final_mass
    double precision              :: initial_eint, final_eint
    double precision              :: initial_eden, final_eden
    
    allocate(fac(lo(1):hi(1),lo(2):hi(2)))
    
    initial_mass = 0.d0
    final_mass = 0.d0
    initial_eint = 0.d0
    final_eint = 0.d0
    initial_eden = 0.d0
    final_eden = 0.d0
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          initial_mass = initial_mass + uout(i,j,URHO)
          initial_eint = initial_eint + uout(i,j,UEINT)
          initial_eden = initial_eden + uout(i,j,UEDEN)
          
          if (uout(i,j,URHO) .eq. 0.d0) then
             
             print *,'   '
             print *,'>>> Error: Castro_2d::enforce_minimum_density ',i,j
             print *,'>>> ... density exactly zero in grid ',lo(1),hi(1),lo(2),hi(2)
             print *,'    '
             call bl_error("Error:: Castro_2d.f90 :: enforce_minimum_density")
             
          else if (uout(i,j,URHO) < small_dens) then
             
             min_dens = uin(i,j,URHO)
             do jj = -1,1
                do ii = -1,1
                   min_dens = min(min_dens,uin(i+ii,j+jj,URHO))
                   if (i+ii.ge.uout_l1 .and. j+jj.ge.uout_l2 .and. &
                       i+ii.le.uout_h1 .and. j+jj.le.uout_h2) then
                      if (uout(i+ii,j+jj,URHO) > small_dens) &
                           min_dens = min(min_dens,uout(i+ii,j+jj,URHO))
                   endif
                end do
             end do
             
             if (verbose .gt. 0) then
                if (uout(i,j,URHO) < 0.d0) then
                   print *,'   '
                   print *,'>>> Warning: Castro_2d::enforce_minimum_density ',i,j
                   print *,'>>> ... resetting negative density '
                   print *,'>>> ... from ',uout(i,j,URHO),' to ',min_dens
                   print *,'    '
                else
                   print *,'   '
                   print *,'>>> Warning: Castro_2d::enforce_minimum_density ',i,j
                   print *,'>>> ... resetting small density '
                   print *,'>>> ... from ',uout(i,j,URHO),' to ',min_dens
                   print *,'    '
                end if
             end if
             
             fac(i,j) = min_dens / uout(i,j,URHO)
             
          end if
          
       enddo
    enddo
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          
          if (uout(i,j,URHO) < small_dens) then
             
             uout(i,j,URHO ) = uout(i,j,URHO ) * fac(i,j)
             uout(i,j,UEINT) = uout(i,j,UEINT) * fac(i,j)
             uout(i,j,UEDEN) = uout(i,j,UEDEN) * fac(i,j)
             uout(i,j,UMX  ) = uout(i,j,UMX  ) * fac(i,j)
             uout(i,j,UMY  ) = uout(i,j,UMY  ) * fac(i,j)
             
             do n = UFS, UFS+nspec-1
                uout(i,j,n) = uout(i,j,n) * fac(i,j)
             end do
             do n = UFX, UFX+naux-1
                uout(i,j,n) = uout(i,j,n) * fac(i,j)
             end do
             do n = UFA, UFA+nadv-1
                uout(i,j,n) = uout(i,j,n) * fac(i,j)
             end do
             
          end if
          
          final_mass = final_mass + uout(i,j,URHO)
          final_eint = final_eint + uout(i,j,UEINT)
          final_eden = final_eden + uout(i,j,UEDEN)
          
       enddo
    enddo

    mass_added = mass_added + (final_mass - initial_mass)
    eint_added = eint_added + (final_eint - initial_eint)
    eden_added = eden_added + (final_eden - initial_eden)
    
  end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_h1,u_h2,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    
    implicit none

    integer          :: lo(2), hi(2)
    integer          :: u_l1,u_l2,u_h1,u_h2
    double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    double precision :: fac,sum
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          sum = 0.d0
          do n = UFS, UFS+nspec-1
             sum = sum + u(i,j,n)
          end do
          if (sum .ne. 0.d0) then
             fac = u(i,j,URHO) / sum
          else
             fac = 1.d0
          end if
          do n = UFS, UFS+nspec-1
             u(i,j,n) = u(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_new_species
  
end module advection_module
