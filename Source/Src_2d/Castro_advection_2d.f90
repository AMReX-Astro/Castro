module advection_module

  implicit none

  private

  public umeth2d, ctoprim, consup
  
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
                     ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
                     flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                     flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                     q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                     q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                     area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                     area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                     pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                     dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                     domlo, domhi)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, ppm_type, hybrid_riemann, &
                                   GDU, GDV, GDPRES, ngdnv
    use trace_module, only : trace
    use trace_ppm_module, only : trace_ppm
    use transverse_module
    use riemann_module, only: cmpflx, shock
    use bl_constants_module

    implicit none

    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer dloga_l1, dloga_l2, dloga_h1, dloga_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer fd1_l1, fd1_l2, fd1_h1, fd1_h2
    integer fd2_l1, fd2_l2, fd2_h1, fd2_h2
    integer q1_l1, q1_l2, q1_h1, q1_h2
    integer q2_l1, q2_l2, q2_h1, q2_h2
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
    double precision  srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    double precision q1(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv)
    double precision q2(q2_l1:q2_h1,q2_l2:q2_h2,ngdnv)
    double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
    double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
    double precision pdivu(ilo1:ihi1,ilo2:ihi2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable::  qm(:,:,:),   qp(:,:,:)
    double precision, allocatable:: qxm(:,:,:), qym(:,:,:)
    double precision, allocatable:: qxp(:,:,:), qyp(:,:,:)
    
    ! Work arrays to hold riemann state and conservative fluxes
    double precision, allocatable ::  fx(:,:,:),  fy(:,:,:)
    double precision, allocatable ::  qgdxtmp(:,:,:)
    double precision, allocatable :: shk(:,:)

    ! Local scalar variables
    double precision :: dtdx
    double precision :: hdtdx, hdt, hdtdy
    integer          :: i,j

    allocate ( qgdxtmp(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv))

    allocate (  qm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate (  qp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
    allocate (  fx(ilo1  :ihi1+1,ilo2-1:ihi2+1,NVAR))
    allocate (  fy(ilo1-1:ihi1+1,ilo2  :ihi2+1,NVAR))

    allocate (shk(ilo1-1:ihi1+1,ilo2-1:ihi2+1))


    ! Local constants
    dtdx = dt/dx
    hdtdx = HALF*dtdx
    hdtdy = HALF*dt/dy
    hdt = HALF*dt


    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                  shk,ilo1-1,ilo2-1,ihi1+1,ihi2+1, &
                  ilo1,ilo2,ihi1,ihi2,dx,dy)
    else
       shk(:,:) = ZERO
    endif

    
    ! NOTE: Geometry terms need to be punched through

    ! Trace to edges w/o transverse flux correction terms.  Here,
    !      qxm and qxp will be the states on either side of the x interfaces
    ! and  qym and qyp will be the states on either side of the y interfaces
    if (ppm_type .eq. 0) then
       call trace(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                  dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                  qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                  srcQ,src_l1,src_l2,src_h1,src_h2, &
                  ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    else
       call trace_ppm(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                      qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      srcQ,src_l1,src_l2,src_h1,src_h2, &
                      gamc,qd_l1,qd_l2,qd_h1,qd_h2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    end if

    ! Solve the Riemann problem in the x-direction using these first
    ! guesses for the x-interface states.  This produces the flux fx
    call cmpflx(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                qgdxtmp, q1_l1, q1_l2, q1_h1, q1_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                1, ilo1, ihi1, ilo2-1, ihi2+1, domlo, domhi)

    ! Solve the Riemann problem in the y-direction using these first
    ! guesses for the y-interface states.  This produces the flux fy
    call cmpflx(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                2, ilo1-1, ihi1+1, ilo2, ihi2, domlo, domhi)

    ! Correct the x-interface states (qxm, qxp) by adding the
    ! transverse flux difference in the y-direction to the x-interface
    ! states.  This results in the new x-interface states qm and qp
    call transy(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ, src_l1, src_l2, src_h1, src_h2, &
                hdt, hdtdy, &
                ilo1-1, ihi1+1, ilo2, ihi2)
    
    ! Solve the final Riemann problem across the x-interfaces with the
    ! full unsplit states.  The resulting flux through the x-interfaces
    ! is flux1
    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                1, ilo1, ihi1, ilo2, ihi2, domlo, domhi)
      
    ! Correct the y-interface states (qym, qyp) by adding the
    ! transverse flux difference in the x-direction to the y-interface
    ! states.  This results in the new y-interface states qm and qp
    call transx(qym, qm, qyp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                qgdxtmp, q1_l1, q1_l2, q1_h1, q1_h2, &
                gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
                hdt, hdtdx, &
                area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                ilo1, ihi1, ilo2-1, ihi2+1)

    ! Solve the final Riemann problem across the y-interfaces with the
    ! full unsplit states.  The resulting flux through the y-interfaces
    ! is flux2
    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                2, ilo1, ihi1, ilo2, ihi2, domlo, domhi)
      
    ! Construct p div{U} -- this will be used as a source to the internal
    ! energy update.  Note we construct this using the interface states
    ! returned from the Riemann solver.
    do j = ilo2,ihi2
       do i = ilo1,ihi1
          pdivu(i,j) = HALF*( &
               (q1(i+1,j,GDPRES) + q1(i,j,GDPRES)) * &
               (q1(i+1,j,GDU)*area1(i+1,j) - q1(i,j,GDU)*area1(i,j)) + &
               (q2(i,j+1,GDPRES) + q2(i,j,GDPRES)) * &
               (q2(i,j+1,GDV)*area2(i,j+1) - q2(i,j,GDV)*area2(i,j)) ) / vol(i,j)
       end do
    end do

    deallocate(qm,qp,qxm,qxp,qym,qyp)
    deallocate(fx,fy)
    deallocate(shk)
    deallocate(qgdxtmp)
    
  end subroutine umeth2d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi, &
                     uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                     q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
                     src,src_l1,src_l2,src_h1,src_h2, &
                     srcQ,srQ_l1,srQ_l2,srQ_h1,srQ_h2, &
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
                                   QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QTEMP, QGAME, &
                                   QFS, QFX, &
                                   allow_negative_energy, small_temp, use_flattening, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1
    use flatten_module
    use bl_constants_module

    implicit none
    
    double precision, parameter:: small = 1.d-8

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer q_l1,q_l2,q_h1,q_h2
    integer src_l1,src_l2,src_h1,src_h2
    integer srQ_l1,srQ_l2,srQ_h1,srQ_h2
    
    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    double precision :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
    double precision :: c(q_l1:q_h1,q_l2:q_h2)
    double precision :: gamc(q_l1:q_h1,q_l2:q_h2)
    double precision :: csml(q_l1:q_h1,q_l2:q_h2)
    double precision :: flatn(q_l1:q_h1,q_l2:q_h2)
    double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
    double precision :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,QVAR)
    double precision :: dx, dy, dt, courno
    
    double precision, allocatable :: dpdrho(:,:)
    double precision, allocatable :: dpde(:,:)
    double precision, allocatable :: dpdX_er(:,:,:)

    integer          :: i, j
    integer          :: ngp, ngf, loq(2), hiq(2)
    integer          :: n, nq, ipassive
    double precision :: courx, coury, courmx, courmy
    double precision :: kineng

    type (eos_t) :: eos_state

    allocate( dpdrho(q_l1:q_h1,q_l2:q_h2))
    allocate(   dpde(q_l1:q_h1,q_l2:q_h2))
    allocate(dpdX_er(q_l1:q_h1,q_l2:q_h2,nspec))
    
    do i=1,2
       loq(i) = lo(i)-ngp
       hiq(i) = hi(i)+ngp
    enddo

    ! Make q (all but p), except put e in slot for rho.e, fix after
    ! eos call The temperature is used as an initial guess for the eos
    ! call and will be overwritten
    do j = loq(2),hiq(2)
       do i = loq(1),hiq(1)

          q(i,j,QRHO) = uin(i,j,URHO)
          
          ! Load passively-advected quatities, c, into q, assuming they
          ! arrived in uin as rho.c. Note that for DIM < 3, this includes
          ! the transverse velocities that are not explicitly evolved.
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nq = qpass_map(ipassive)

             q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
          enddo          
          
          if (uin(i,j,URHO) .le. ZERO) then
             print *,'   '
             print *,'>>> Error: Castro_2d::ctoprim ',i,j
             print *,'>>> ... negative density ',uin(i,j,URHO)
             print *,'    '
             call bl_error("Error:: Castro_2d.f90 :: ctoprim")
          end if
          
          q(i,j,QU:QV) = uin(i,j,UMX:UMY)/uin(i,j,URHO)

          ! Get the internal energy, which we'll use for determining the pressure.
          ! We use a dual energy formalism. If (E - K) < eta1 and eta1 is suitably small, 
          ! then we risk serious numerical truncation error in the internal energy.
          ! Therefore we'll use the result of the separately updated internal energy equation.
          ! Otherwise, we'll set e = E - K.

          kineng = HALF * q(i,j,QRHO) * sum(q(i,j,QU:QW)**2)
          
          if ( (uin(i,j,UEDEN) - kineng) / uin(i,j,UEDEN) .gt. dual_energy_eta1) then
             q(i,j,QREINT) = (uin(i,j,UEDEN) - kineng) / q(i,j,QRHO)
          else
             q(i,j,QREINT) = uin(i,j,UEINT) / q(i,j,QRHO)
          endif

          q(i,j,QTEMP  ) = uin(i,j,UTEMP)

          ! Get gamc, p, T, c, csml using q state 
          eos_state % T   = q(i,j,QTEMP)
          eos_state % rho = q(i,j,QRHO)
          eos_state % xn  = q(i,j,QFS:QFS+nspec-1)
          eos_state % aux = q(i,j,QFX:QFX+naux-1)

          ! if necessary, reset the energy using small_temp
          if ((allow_negative_energy .eq. 0) .and. (q(i,j,QREINT) .lt. ZERO)) then
             q(i,j,QTEMP) = small_temp
             eos_state % T = q(i,j,QTEMP)

             call eos(eos_input_rt, eos_state)
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

          call eos(eos_input_re, eos_state)

          q(i,j,QTEMP)  = eos_state % T
          q(i,j,QREINT) = eos_state % e
          q(i,j,QPRES)  = eos_state % p

          dpdrho(i,j) = eos_state % dpdr_e
          dpde(i,j)   = eos_state % dpde
          c(i,j)      = eos_state % cs
          gamc(i,j)   = eos_state % gam1

          csml(i,j) = max(small, small * c(i,j))

          ! Make this "rho e" instead of "e"
          q(i,j,QREINT) = q(i,j,QREINT)*q(i,j,QRHO)
          q(i,j,QGAME) = q(i,j,QPRES)/q(i,j,QREINT) + ONE
       enddo
    enddo
    
    ! Compute sources in terms of Q
    do j = loq(2), hiq(2)
       do i = loq(1), hiq(1)
          
          srcQ(i,j,QRHO  ) = src(i,j,URHO)
          srcQ(i,j,QU:QV ) = (src(i,j,UMX:UMY) - q(i,j,QU:QV) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
          ! S_rhoe = S_rhoE - u . (S_rhoU - 0.5 u S_rho)
          srcQ(i,j,QREINT) = src(i,j,UEDEN)                               &
                           - dot_product(q(i,j,QU:QV),src(i,j,UMX:UMY))   &
                           + HALF * sum(q(i,j,QU:QV)**2) * srcQ(i,j,QRHO)
          srcQ(i,j,QPRES ) = dpde(i,j) * &
               (srcQ(i,j,QREINT) - q(i,j,QREINT)*srcQ(i,j,QRHO)/q(i,j,QRHO))/q(i,j,QRHO) + &
               dpdrho(i,j) * srcQ(i,j,QRHO)! + &
!                sum(dpdX_er(i,j,:)*(src(i,j,UFS:UFS+nspec-1) - &
!                    q(i,j,QFS:QFS+nspec-1)*srcQ(i,j,QRHO))) / q(i,j,QRHO)

       enddo
    enddo

    ! and the passive advective quantities sources
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do j = loq(2), hiq(2)
          do i = loq(1), hiq(1)
             srcQ(i,j,nq) = ( src(i,j,n) - q(i,j,nq) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
          enddo
       enddo

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
       call uflaten((/ loq(1), loq(2), 0 /), (/ hiq(1), hiq(2), 0 /), &
            q(q_l1,q_l2,QPRES), &
            q(q_l1,q_l2,QU), &
            q(q_l1,q_l2,QV), &
            q(q_l1,q_l2,QW), &
            flatn,(/ q_l1, q_l2, 0 /), (/ q_h1, q_h2, 0 /))
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
                     q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                     q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                     src , src_l1, src_l2, src_h1, src_h2, &
                     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                     area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                     area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                     vol,vol_l1,vol_l2,vol_h1,vol_h2, &
                     div,pdivu,lo,hi,dx,dy,dt,E_added_flux, &
                     xmom_added_flux,ymom_added_flux,zmom_added_flux, &
                     verbose)

    use eos_module
    use network, only : nspec, naux
    use meth_params_module, only : difmag, NVAR, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, ngdnv, GDPRES, &
                                   normalize_species
    use prob_params_module, only : coord_type
    use bl_constants_module
    use advection_util_module, only : normalize_species_fluxes

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer uout_l1,uout_l2,uout_h1,uout_h2
    integer q1_l1, q1_l2, q1_h1, q1_h2
    integer q2_l1, q2_l2, q2_h1, q2_h2
    integer   src_l1,  src_l2,  src_h1,  src_h2
    integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
    integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
    integer area1_l1,area1_l2,area1_h1,area1_h2
    integer area2_l1,area2_l2,area2_h1,area2_h2
    integer vol_l1,vol_l2,vol_h1,vol_h2

    integer verbose

    double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision q1(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv)
    double precision q2(q2_l1:q2_h1,q2_l2:q2_h2,ngdnv)
    double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,NVAR)
    double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
    double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
    double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
    double precision pdivu(lo(1):hi(1),lo(2):hi(2))
    double precision dx, dy, dt, E_added_flux
    double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
    
    integer i, j, n

    double precision div1
    !double precision rho, Up, Vp, SrE

    ! Normalize the species fluxes
    if (normalize_species == 1) &
         call normalize_species_fluxes( &
                flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                lo,hi)

    ! correct the fluxes to include the effects of the artificial viscosity
    do n = 1, NVAR
       if (n == UTEMP) then
          flux1(:,:,n) = ZERO
          flux2(:,:,n) = ZERO
       else 
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1
                div1 = HALF*(div(i,j) + div(i,j+1))
                div1 = difmag*min(ZERO,div1)

                flux1(i,j,n) = flux1(i,j,n) + &
                     dx*div1*(uin(i,j,n) - uin(i-1,j,n))

                flux1(i,j,n) = area1(i,j)*flux1(i,j,n)
             enddo
          enddo
          
          do j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                div1 = HALF*(div(i,j) + div(i+1,j))
                div1 = difmag*min(ZERO,div1)

                flux2(i,j,n) = flux2(i,j,n) + &
                     dy*div1*(uin(i,j,n) - uin(i,j-1,n))

                flux2(i,j,n) = area2(i,j)*flux2(i,j,n)
             enddo
          enddo

       endif
    enddo
    
    ! do the conservative updates
    do n = 1, NVAR
       if (n == UTEMP) then
          uout(lo(1):hi(1),lo(2):hi(2),n) = uin(lo(1):hi(1),lo(2):hi(2),n)
       else 
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                uout(i,j,n) = uin(i,j,n) + dt * &
                     ( flux1(i,j,n) - flux1(i+1,j,n) + &
                       flux2(i,j,n) - flux2(i,j+1,n) ) / vol(i,j)
                
                if (n == UEINT) then
                   ! Add p div(u) source term to (rho e)
                   uout(i,j,UEINT) = uout(i,j,UEINT)  - dt * pdivu(i,j)
                endif
                   
             enddo
          enddo
       end if
    enddo

    ! Add up some diagnostic quantities. Note that these are volumetric sums
    ! so we are not dividing by the cell volume.
                   
    if (verbose .eq. 1) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             xmom_added_flux = xmom_added_flux + dt * ( flux1(i,j,UMX) - flux1(i+1,j,UMX) + &
                                                        flux2(i,j,UMX) - flux2(i,j+1,UMX) )

             ymom_added_flux = ymom_added_flux + dt * ( flux1(i,j,UMY) - flux1(i+1,j,UMY) + &
                                                        flux2(i,j,UMY) - flux2(i,j+1,UMY) )

             zmom_added_flux = zmom_added_flux + dt * ( flux1(i,j,UMZ) - flux1(i+1,j,UMZ) + &
                                                        flux2(i,j,UMZ) - flux2(i,j+1,UMZ) )

             E_added_flux = E_added_flux + dt * ( flux1(i,j,UEDEN) - flux1(i+1,j,UEDEN) + &
                                                  flux2(i,j,UEDEN) - flux2(i,j+1,UEDEN) )

          enddo
       enddo

    endif


    ! Add gradp term to momentum equation -- only for axisymmetric
    ! coords (and only for the radial flux)
    if (coord_type == 1) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,UMX) = uout(i,j,UMX) - dt * (q1(i+1,j,GDPRES) - q1(i,j,GDPRES))/ dx
             !uout(i,j,UMY) = uout(i,j,UMY) - dt * (pgdy(i,j+1)-pgdy(i,j))/ dy
          enddo
       enddo
    endif

    ! scale the fluxes (and correct the momentum flux with the grad p part)
    ! so we can use them in the flux correction at coarse-fine interfaces
    ! later.
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR)
          if (coord_type == 1) then
             flux1(i,j,UMX) = flux1(i,j,UMX) + dt*area1(i,j)*q1(i,j,GDPRES)
          endif
       enddo
    enddo
    
    do j = lo(2), hi(2)+1 
       do i = lo(1), hi(1)
          flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR)
          !flux2(i,j,UMY) = flux2(i,j,UMY) + dt*area2(i,j)*pgdy(i,j)
       enddo
    enddo
    
  end subroutine consup

end module advection_module
