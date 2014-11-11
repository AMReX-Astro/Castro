module advection_module

  implicit none

  private

  public umeth2d, ctoprim, divu, consup, enforce_minimum_density, &
       normalize_new_species, &
       normalize_species_fluxes
  
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
    use transverse_module
    use riemann_module, only: cmpflx
    use bl_constants_module

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
    
    ! Work arrays to hold riemann state and conservative fluxes
    double precision, allocatable::   fx(:,:,:),  fy(:,:,:)
    double precision, allocatable::   pgdxtmp(:,:) ,  ugdxtmp(:,:)
    double precision, allocatable :: gegdxtmp(:,:), gegdx(:,:), gegdy(:,:)

    ! Local scalar variables
    double precision :: dtdx
    double precision :: hdtdx, hdt, hdtdy
    integer          :: i,j

    allocate ( pgdxtmp(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2))
    allocate ( ugdxtmp(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2))
    allocate ( gegdxtmp(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2))
    allocate ( gegdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2))
    allocate ( gegdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2))

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
    hdtdx = HALF*dtdx
    hdtdy = HALF*dt/dy
    hdt = HALF*dt
    
    ! NOTE: Geometry terms need to be punched through

    ! Trace to edges w/o transverse flux correction terms.  Here,
    !      qxm and qxp will be the states on either side of the x interfaces
    ! and  qym and qyp will be the states on either side of the y interfaces
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

    ! Solve the Riemann problem in the x-direction using these first
    ! guesses for the x-interface states.  This produces the flux fx
    call cmpflx(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                pgdxtmp, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                ugdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gegdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &                
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                1, ilo1, ihi1, ilo2-1, ihi2+1, domlo, domhi)

    ! Solve the Riemann problem in the y-direction using these first
    ! guesses for the y-interface states.  This produces the flux fy
    call cmpflx(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gegdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                2, ilo1-1, ihi1+1, ilo2, ihi2, domlo, domhi)

    ! Correct the x-interface states (qxm, qxp) by adding the
    ! transverse flux difference in the y-direction to the x-interface
    ! states.  This results in the new x-interface states qm and qp
    call transy(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gegdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ, src_l1, src_l2, src_h1, src_h2, &
                grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                hdt, hdtdy, &
                ilo1-1, ihi1+1, ilo2, ihi2)
    
    ! Solve the final Riemann problem across the x-interfaces with the
    ! full unsplit states.  The resulting flux through the x-interfaces
    ! is flux1
    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gegdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                1, ilo1, ihi1, ilo2, ihi2, domlo, domhi)
      
    ! Correct the y-interface states (qym, qyp) by adding the
    ! transverse flux difference in the x-direction to the y-interface
    ! states.  This results in the new y-interface states qm and qp
    call transx(qym, qm, qyp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                pgdxtmp, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                ugdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gegdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
                grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                hdt, hdtdx, &
                area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                ilo1, ihi1, ilo2-1, ihi2+1)

    ! Solve the final Riemann problem across the y-interfaces with the
    ! full unsplit states.  The resulting flux through the y-interfaces
    ! is flux2
    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gegdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                2, ilo1, ihi1, ilo2, ihi2, domlo, domhi)
      

    ! Construct p div{U} -- this will be used as a source to the internal
    ! energy update.  Note we construct this using the interface states
    ! returned from the Riemann solver.
    do j = ilo2,ihi2
       do i = ilo1,ihi1
          pdivu(i,j) = HALF * &
               ((pgdx(i+1,j)+pgdx(i,j))*(ugdx(i+1,j)*area1(i+1,j)-ugdx(i,j)*area1(i,j)) &
               +(pgdy(i,j+1)+pgdy(i,j))*(ugdy(i,j+1)*area2(i,j+1)-ugdy(i,j)*area2(i,j)) ) / vol(i,j)
       end do
    end do

    deallocate(dq,qm,qp,qxm,qxp,qym,qyp)
    deallocate(fx,fy)
    deallocate(pgdxtmp,ugdxtmp,gegdxtmp,gegdx,gegdy)
    
  end subroutine umeth2d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi, &
                     uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                     q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
                     src,srcQ,src_l1,src_l2,src_h1,src_h2, &
                     courno,dx,dy,dt,ngp,ngf)
    
    ! Will give primitive variables on lo-ngp:hi+ngp, and flatn on
    ! lo-ngf:hi+ngf if use_flattening=1.  Declared dimensions of
    ! q,c,gamc,csml,flatn are given by DIMS(q).  This declared region
    ! is assumed to encompass lo-ngp:hi+ngp.  Also, uflaten call
    ! assumes ngp>=ngf+3 (ie, primitve data is used by the routine
    ! that computes flatn).

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP,&
                                   UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QREINT, QPRES, QTEMP, QGAME, &
                                   QFA, QFS, QFX, &
                                   nadv, allow_negative_energy, small_temp, use_flattening
    use flatten_module
    use bl_constants_module

    implicit none
    
    double precision, parameter:: small = 1.d-8

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer q_l1,q_l2,q_h1,q_h2
    integer src_l1,src_l2,src_h1,src_h2
    
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

    type (eos_t) :: eos_state
    
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

          if (uin(i,j,URHO) .le. ZERO) then
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

          pt_index(1) = i
          pt_index(2) = j

          eos_state % T   = q(i,j,QTEMP)
          eos_state % rho = q(i,j,QRHO)
          eos_state % xn  = q(i,j,QFS:QFS+nspec-1)
          eos_state % aux = q(i,j,QFX:QFX+naux-1)

          ! If necessary, reset the energy using small_temp
          if ((allow_negative_energy .eq. 0) .and. (q(i,j,QREINT) .lt. ZERO)) then
             q(i,j,QTEMP) = small_temp
             eos_state % T = q(i,j,QTEMP)

             call eos(eos_input_rt, eos_state, pt_index = pt_index)
             q(i,j,QREINT) = eos_state % e

             if (q(i,j,QREINT) .lt. ZERO) then
                print *,'   '
                print *,'>>> Error: Castro_2d::ctoprim ',i,j
                print *,'>>> ... new e from eos (input_rt) call is negative ',q(i,j,QREINT)
                print *,'    '
                call bl_error("Error:: Castro_2d.f90 :: ctoprim")
             end if
          end if

          eos_state % e = q(i,j,QREINT)

          call eos(eos_input_re, eos_state, pt_index = pt_index)

          q(i,j,QTEMP)  = eos_state % T
          q(i,j,QREINT) = eos_state % e
          q(i,j,QPRES)  = eos_state % p

          dpdrho(i,j) = eos_state % dpdr_e
          dpde(i,j)   = eos_state % dpde
          c(i,j)      = eos_state % cs
          gamc(i,j)   = eos_state % gam1

          csml(i,j) = max(small, small * c(i,j))
       end do
    end do

    ! Make this "rho e" instead of "e"
    do j = loq(2),hiq(2)
       do i = loq(1),hiq(1)
          q(i,j,QREINT) = q(i,j,QREINT)*q(i,j,QRHO)

          q(i,j,QGAME) = q(i,j,QPRES)/q(i,j,QREINT) + ONE

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
               HALF * (q(i,j,QU)**2 + q(i,j,QV)**2) * srcQ(i,j,QRHO)
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
          
          if (courx .gt. ONE) then
             print *,'   '
             call bl_warning("Warning:: Castro_2d.f90 :: CFL violation in ctoprim")
             print *,'>>> ... (u+c) * dt / dx > 1 ', courx
             print *,'>>> ... at cell (i,j)     : ',i,j
             print *,'>>> ... u, c                ',q(i,j,QU), c(i,j)
             print *,'>>> ... density             ',q(i,j,QRHO)
          end if
          
          if (coury .gt. ONE) then
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
    if (use_flattening == 1) then
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
       flatn = ONE
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
    use bl_constants_module

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
          flux1(:,:,n) = ZERO
          flux2(:,:,n) = ZERO
       else 
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                div1 = HALF*(div(i,j) + div(i,j+1))
                div1 = difmag*min(ZERO,div1)
                flux1(i,j,n) = flux1(i,j,n) &
                     + dx*div1*(uin(i,j,n) - uin(i-1,j,n))
                flux1(i,j,n) = area1(i,j)*flux1(i,j,n)
             enddo
          enddo
          
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                div1 = HALF*(div(i,j) + div(i+1,j))
                div1 = difmag*min(ZERO,div1)
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
!         uout(i,j,UMX) = uout(i,j,UMX)+ HALF*(area1(i,j)+area1(i+1,j))* &
!            dt * ( pgdx(i,j)-pgdx(i+1,j) )/vol(i,j)
!         uout(i,j,UMY) = uout(i,j,UMY)+ HALF*(area2(i,j)+area2(i,j+1))* &
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

  subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,dx, &
                  div,div_l1,div_l2,div_h1,div_h2)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV
    use bl_constants_module
    
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
             ux = HALF*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/dx(1)
             vy = HALF*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/dx(2)
             div(i,j) = ux + vy
          enddo
       enddo
    else
       do i=lo(1),hi(1)+1
          
          if (i.eq.0) then
             
             div(i,lo(2):hi(2)+1) = ZERO
             
          else 

             rl = (dble(i)-HALF) * dx(1)
             rr = (dble(i)+HALF) * dx(1)
             rc = (dble(i)     ) * dx(1)
             
             do j=lo(2),hi(2)+1
                ! These are transverse averages in the y-direction
                ul = HALF * (q(i-1,j,QU)+q(i-1,j-1,QU))
                ur = HALF * (q(i  ,j,QU)+q(i  ,j-1,QU))
                
                ! Take 1/r d/dr(r*u)
                div(i,j) = (rr*ur - rl*ul) / dx(1) / rc
                
                ! These are transverse averages in the x-direction
                vb = HALF * (q(i,j-1,QV)+q(i-1,j-1,QV))
                vt = HALF * (q(i,j  ,QV)+q(i-1,j  ,QV))
                
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
    use bl_constants_module
    
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
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux1(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux1(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             flux1(i,j,n) = flux1(i,j,n) * fac
          end do
       end do
    end do
    do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + flux2(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = flux2(i,j,URHO) / sum
          else
             fac = ONE
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
    use bl_constants_module

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
    
    initial_mass = ZERO
    final_mass = ZERO
    initial_eint = ZERO
    final_eint = ZERO
    initial_eden = ZERO
    final_eden = ZERO
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          initial_mass = initial_mass + uout(i,j,URHO)
          initial_eint = initial_eint + uout(i,j,UEINT)
          initial_eden = initial_eden + uout(i,j,UEDEN)
          
          if (uout(i,j,URHO) .eq. ZERO) then
             
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
                if (uout(i,j,URHO) < ZERO) then
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
    use bl_constants_module
    
    implicit none

    integer          :: lo(2), hi(2)
    integer          :: u_l1,u_l2,u_h1,u_h2
    double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
    
    ! Local variables
    integer          :: i,j,n
    double precision :: fac,sum
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          sum = ZERO
          do n = UFS, UFS+nspec-1
             sum = sum + u(i,j,n)
          end do
          if (sum .ne. ZERO) then
             fac = u(i,j,URHO) / sum
          else
             fac = ONE
          end if
          do n = UFS, UFS+nspec-1
             u(i,j,n) = u(i,j,n) * fac
          end do
       end do
    end do
    
  end subroutine normalize_new_species
  
end module advection_module
