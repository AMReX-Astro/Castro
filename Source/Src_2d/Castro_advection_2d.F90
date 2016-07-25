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

  subroutine umeth2d(q, flatn, qd_l1, qd_l2, qd_h1, qd_h2,&
                     srcQ, src_l1, src_l2, src_h1, src_h2, &
                     ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
                     uout, uout_l1, uout_l2, uout_h1, uout_h2, &
                     flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                     flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                     q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                     q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                     area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                     area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                     pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                     dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                     domlo, domhi)

    use meth_params_module, only : QVAR, NVAR, ppm_type, hybrid_riemann, &
                                   GDU, GDV, GDPRES, ngdnv, &
                                   QC, QCSML, QGAMC
    use trace_module, only : trace
    use trace_ppm_module, only : trace_ppm
    use transverse_module, only : transx, transy
    use riemann_module, only: cmpflx, shock
    use bl_constants_module, only : ZERO, HALF
#ifdef SHOCK_VAR
    use meth_params_,odule, only : USHK
#endif

    implicit none

    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer dloga_l1, dloga_l2, dloga_h1, dloga_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer uout_l1, uout_l2, uout_h1, uout_h2
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
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision  srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
    double precision q1(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv)
    double precision q2(q2_l1:q2_h1,q2_l2:q2_h2,ngdnv)
    double precision  uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
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

#ifdef SHOCK_VAR
    uout(:,:,USHK) = ZERO

    call shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
               shk,ilo1-1,ilo2-1,ihi1+1,ihi2+1, &
               ilo1,ilo2,ihi1,ihi2,dx,dy)

    ! Store the shock data for future use in the burning step.

    do j = ilo2, ihi2
       do i = ilo1, ihi1
          uout(i,j,USHK) = shk(i,j)
       enddo
    enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                  shk,ilo1-1,ilo2-1,ihi1+1,ihi2+1, &
                  ilo1,ilo2,ihi1,ihi2,dx,dy)
    else
       shk(:,:) = ZERO
    endif
#endif

    ! NOTE: Geometry terms need to be punched through

    ! Trace to edges w/o transverse flux correction terms.  Here,
    !      qxm and qxp will be the states on either side of the x interfaces
    ! and  qym and qyp will be the states on either side of the y interfaces
    if (ppm_type .eq. 0) then
       call trace(q,q(:,:,QC),flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                  dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                  qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                  srcQ,src_l1,src_l2,src_h1,src_h2, &
                  ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    else
       call trace_ppm(q,q(:,:,QC),flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                      qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      srcQ,src_l1,src_l2,src_h1,src_h2, &
                      q(:,:,QGAMC),qd_l1,qd_l2,qd_h1,qd_h2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
    end if

    ! Solve the Riemann problem in the x-direction using these first
    ! guesses for the x-interface states.  This produces the flux fx
    call cmpflx(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                qgdxtmp, q1_l1, q1_l2, q1_h1, q1_h2, &
                q(:,:,QGAMC), q(:,:,QCSML), q(:,:,QC), qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                1, ilo1, ihi1, ilo2-1, ihi2+1, domlo, domhi)

    ! Solve the Riemann problem in the y-direction using these first
    ! guesses for the y-interface states.  This produces the flux fy
    call cmpflx(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                q(:,:,QGAMC), q(:,:,QCSML), q(:,:,QC), qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                2, ilo1-1, ihi1+1, ilo2, ihi2, domlo, domhi)

    ! Correct the x-interface states (qxm, qxp) by adding the
    ! transverse flux difference in the y-direction to the x-interface
    ! states.  This results in the new x-interface states qm and qp
    call transy(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                q(:,:,QGAMC), qd_l1, qd_l2, qd_h1, qd_h2, &
                srcQ, src_l1, src_l2, src_h1, src_h2, &
                hdt, hdtdy, &
                ilo1-1, ihi1+1, ilo2, ihi2)

    ! Solve the final Riemann problem across the x-interfaces with the
    ! full unsplit states.  The resulting flux through the x-interfaces
    ! is flux1
    call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                q(:,:,QGAMC), q(:,:,QCSML), q(:,:,QC), qd_l1, qd_l2, qd_h1, qd_h2, &
                shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
                1, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    ! Correct the y-interface states (qym, qyp) by adding the
    ! transverse flux difference in the x-direction to the y-interface
    ! states.  This results in the new y-interface states qm and qp
    call transx(qym, qm, qyp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                qgdxtmp, q1_l1, q1_l2, q1_h1, q1_h2, &
                q(:,:,QGAMC), qd_l1, qd_l2, qd_h1, qd_h2, &
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
                q(:,:,QGAMC), q(:,:,QCSML), q(:,:,QC), qd_l1, qd_l2, qd_h1, qd_h2, &
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
                     q,q_l1,q_l2,q_h1,q_h2, &
                     src,src_l1,src_l2,src_h1,src_h2, &
                     srcQ,srQ_l1,srQ_l2,srQ_h1,srQ_h2, &
                     dx,dy,dt)

    use actual_network, only : nspec, naux
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP,&
                                   QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QTEMP, QGAME, &
                                   QFS, QFX, &
                                   QC, QCSML, QGAMC, QDPDR, QDPDE, &
                                   allow_negative_energy, &
                                   npassive, upass_map, qpass_map, dual_energy_eta1
    use bl_constants_module, only : ZERO, HALF, ONE

    implicit none

    double precision, parameter:: small = 1.d-8

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer q_l1,q_l2,q_h1,q_h2
    integer src_l1,src_l2,src_h1,src_h2
    integer srQ_l1,srQ_l2,srQ_h1,srQ_h2

    double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
    double precision :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
    double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
    double precision :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,QVAR)
    double precision :: dx, dy, dt

    integer          :: i, j
    integer          :: n, nq, ipassive
    double precision :: kineng

    type (eos_t) :: eos_state

    if (allow_negative_energy .eq. 0) eos_state % reset = .true.

    ! Make q (all but p), except put e in slot for rho.e, fix after
    ! eos call The temperature is used as an initial guess for the eos
    ! call and will be overwritten
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

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
          eos_state % e   = q(i,j,QREINT)
          eos_state % xn  = q(i,j,QFS:QFS+nspec-1)
          eos_state % aux = q(i,j,QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          q(i,j,QTEMP)  = eos_state % T
          q(i,j,QREINT) = eos_state % e
          q(i,j,QPRES)  = eos_state % p
          q(i,j,QDPDR)  = eos_state % dpdr_e
          q(i,j,QDPDE)  = eos_state % dpde
          q(i,j,QC)     = eos_state % cs
          q(i,j,QGAMC)  = eos_state % gam1
          q(i,j,QCSML)  = max(small, small * q(i,j,QC))

          ! Make this "rho e" instead of "e"
          q(i,j,QREINT) = q(i,j,QREINT)*q(i,j,QRHO)
          q(i,j,QGAME) = q(i,j,QPRES)/q(i,j,QREINT) + ONE
       enddo
    enddo

    srcQ = ZERO

    ! Compute sources in terms of Q
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          srcQ(i,j,QRHO  ) = src(i,j,URHO)
          srcQ(i,j,QU:QV ) = (src(i,j,UMX:UMY) - q(i,j,QU:QV) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
          ! S_rhoe = S_rhoE - u . (S_rhoU - 0.5 u S_rho)
          srcQ(i,j,QREINT) = src(i,j,UEDEN)                               &
                           - dot_product(q(i,j,QU:QV),src(i,j,UMX:UMY))   &
                           + HALF * sum(q(i,j,QU:QV)**2) * srcQ(i,j,QRHO)
          srcQ(i,j,QPRES ) = q(i,j,QDPDE) * &
               (srcQ(i,j,QREINT) - q(i,j,QREINT)*srcQ(i,j,QRHO)/q(i,j,QRHO))/q(i,j,QRHO) + &
               q(i,j,QDPDR) * srcQ(i,j,QRHO)

       enddo
    enddo

    ! and the passive advective quantities sources
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             srcQ(i,j,nq) = ( src(i,j,n) - q(i,j,nq) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
          enddo
       enddo

    end do

  end subroutine ctoprim

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine consup( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                     uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                     update,updt_l1,updt_l2,updt_h1,updt_h2, &
                     q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                     q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                     src , src_l1, src_l2, src_h1, src_h2, &
                     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                     area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                     area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                     vol,vol_l1,vol_l2,vol_h1,vol_h2, &
                     div,pdivu,lo,hi,dx,dy,dt,mass_added_flux,E_added_flux, &
                     xmom_added_flux,ymom_added_flux,zmom_added_flux, &
                     mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                     eden_lost,xang_lost,yang_lost,zang_lost, &
                     verbose)

    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, ngdnv, GDPRES, track_grid_losses
    use prob_params_module, only : coord_type, domlo_level, domhi_level, center
    use bl_constants_module, only : ZERO, HALF
    use advection_util_2d_module, only : normalize_species_fluxes
    use castro_util_module, only : position, linear_to_angular_momentum
    use amrinfo_module, only : amr_level
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    integer lo(2), hi(2)
    integer uin_l1,uin_l2,uin_h1,uin_h2
    integer uout_l1,uout_l2,uout_h1,uout_h2
    integer updt_l1,updt_l2,updt_h1,updt_h2
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
    double precision update(updt_l1:updt_h1,updt_l2:updt_h2,NVAR)
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
    double precision dx, dy, dt, E_added_flux, mass_added_flux
    double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
    double precision mass_lost, xmom_lost, ymom_lost, zmom_lost
    double precision eden_lost, xang_lost, yang_lost, zang_lost

    integer i, j, k, n

    double precision div1
    !double precision rho, Up, Vp, SrE
    integer domlo(3), domhi(3)
    double precision loc(3), ang_mom(3)

    ! Normalize the species fluxes.

    call normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                                  flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                                  lo,hi)

    ! Correct the fluxes to include the effects of the artificial viscosity.

    do n = 1, NVAR
       if (n == UTEMP) then
          flux1(:,:,n) = ZERO
          flux2(:,:,n) = ZERO
#ifdef SHOCK_VAR
       else if (n == USHK) then
          flux1(:,:,n) = ZERO
          flux2(:,:,n) = ZERO
#endif
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

    ! Fill the update array.

    do n = 1, NVAR
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             update(i,j,n) = update(i,j,n) + ( flux1(i,j,n) - flux1(i+1,j,n) + &
                                               flux2(i,j,n) - flux2(i,j+1,n) ) / vol(i,j)

             if (n == UEINT) then

                ! Add p div(u) source term to (rho e)

                update(i,j,n) = update(i,j,n) - pdivu(i,j)

             endif

          enddo
       enddo
    enddo

    ! Add gradp term to momentum equation -- only for axisymmetric
    ! coords (and only for the radial flux).

    if (coord_type == 1) then
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             update(i,j,UMX) = update(i,j,UMX) - (q1(i+1,j,GDPRES) - q1(i,j,GDPRES)) / dx
             !update(i,j,UMY) = update(i,j,UMY) - (pgdy(i,j+1)-pgdy(i,j)) / dy
          enddo
       enddo
    endif

    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (verbose .eq. 1) then

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             mass_added_flux = mass_added_flux + dt * ( flux1(i,j,URHO) - flux1(i+1,j,URHO) + &
                                                        flux2(i,j,URHO) - flux2(i,j+1,URHO) )

             xmom_added_flux = xmom_added_flux + dt * ( flux1(i,j,UMX) - flux1(i+1,j,UMX) + &
                                                        flux2(i,j,UMX) - flux2(i,j+1,UMX) )

             ymom_added_flux = ymom_added_flux + dt * ( flux1(i,j,UMY) - flux1(i+1,j,UMY) + &
                                                        flux2(i,j,UMY) - flux2(i,j+1,UMY) )

             zmom_added_flux = zmom_added_flux + dt * ( flux1(i,j,UMZ) - flux1(i+1,j,UMZ) + &
                                                        flux2(i,j,UMZ) - flux2(i,j+1,UMZ) )

             E_added_flux    = E_added_flux    + dt * ( flux1(i,j,UEDEN) - flux1(i+1,j,UEDEN) + &
                                                        flux2(i,j,UEDEN) - flux2(i,j+1,UEDEN) )

          enddo
       enddo

    endif

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)

       k = 0

       if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

          j = domlo(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.)

             mass_lost = mass_lost - dt * flux2(i,j,URHO)
             xmom_lost = xmom_lost - dt * flux2(i,j,UMX)
             ymom_lost = ymom_lost - dt * flux2(i,j,UMY)
             zmom_lost = zmom_lost - dt * flux2(i,j,UMZ)
             eden_lost = eden_lost - dt * flux2(i,j,UEDEN)

             ang_mom   = linear_to_angular_momentum(loc - center, dt * flux2(i,j,UMX:UMZ))
             xang_lost = xang_lost - ang_mom(1)
             yang_lost = yang_lost - ang_mom(2)
             zang_lost = zang_lost - ang_mom(3)

          enddo

       endif

       if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

          j = domhi(2) + 1
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.)

             mass_lost = mass_lost + dt * flux2(i,j,URHO)
             xmom_lost = xmom_lost + dt * flux2(i,j,UMX)
             ymom_lost = ymom_lost + dt * flux2(i,j,UMY)
             zmom_lost = zmom_lost + dt * flux2(i,j,UMZ)
             eden_lost = eden_lost + dt * flux2(i,j,UEDEN)

             ang_mom   = linear_to_angular_momentum(loc - center, dt * flux2(i,j,UMX:UMZ))
             xang_lost = xang_lost + ang_mom(1)
             yang_lost = yang_lost + ang_mom(2)
             zang_lost = zang_lost + ang_mom(3)

          enddo

       endif

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)
          do j = lo(2), hi(2)

             loc = position(i,j,k,ccx=.false.)

             mass_lost = mass_lost - dt * flux1(i,j,URHO)
             xmom_lost = xmom_lost - dt * flux1(i,j,UMX)
             ymom_lost = ymom_lost - dt * flux1(i,j,UMY)
             zmom_lost = zmom_lost - dt * flux1(i,j,UMZ)
             eden_lost = eden_lost - dt * flux1(i,j,UEDEN)

             ang_mom   = linear_to_angular_momentum(loc - center, dt * flux1(i,j,UMX:UMZ))
             xang_lost = xang_lost - ang_mom(1)
             yang_lost = yang_lost - ang_mom(2)
             zang_lost = zang_lost - ang_mom(3)

          enddo

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1
          do j = lo(2), hi(2)

             loc = position(i,j,k,ccx=.false.)

             mass_lost = mass_lost + dt * flux1(i,j,URHO)
             xmom_lost = xmom_lost + dt * flux1(i,j,UMX)
             ymom_lost = ymom_lost + dt * flux1(i,j,UMY)
             zmom_lost = zmom_lost + dt * flux1(i,j,UMZ)
             eden_lost = eden_lost + dt * flux1(i,j,UEDEN)

             ang_mom   = linear_to_angular_momentum(loc - center, dt * flux1(i,j,UMX:UMZ))
             xang_lost = xang_lost + ang_mom(1)
             yang_lost = yang_lost + ang_mom(2)
             zang_lost = zang_lost + ang_mom(3)

          enddo

       endif

    endif

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

             flux1(i,j,n) = dt * flux1(i,j,n)

             ! Correct the momentum flux with the grad p part.

             if (coord_type == 1) then
                if (n == UMX) then
                   flux1(i,j,n) = flux1(i,j,n) + dt * area1(i,j) * q1(i,j,GDPRES)
                endif
             endif

          enddo
       enddo
    enddo

    do n = 1, NVAR
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)

             flux2(i,j,n) = dt * flux2(i,j,n)

             !if (n == UMY) then
                !flux2(i,j,n) = flux2(i,j,n) + dt * area2(i,j) * pgdy(i,j)
             !endif

          enddo
       enddo
    enddo

  end subroutine consup

end module advection_module
