module advection_module

  implicit none

  private

  public umeth1d, consup

contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH1D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  sound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux      <=  (modify) flux in X direction on X edges
! ::: ----------------------------------------------------------------

  subroutine umeth1d(lo,hi,domlo,domhi, &
                     q,flatn,qd_l1,qd_h1, &
                     srcQ,src_l1,src_h1, &
                     ilo,ihi,dx,dt, &
                     flux ,   fd_l1,   fd_h1, &
                     pgdnv,pgdnv_l1,pgdnv_h1, &
                     ugdnv,ugdnv_l1,ugdnv_h1, &
                     dloga,dloga_l1,dloga_h1)

    use meth_params_module, only : QVAR, NVAR, ppm_type, QC, QCSML, QGAMC
    use riemann_module, only : cmpflx
    use trace_module, only : trace
    use trace_ppm_module, only : trace_ppm

    implicit none

    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer dloga_l1,dloga_h1
    integer qd_l1,qd_h1
    integer src_l1,src_h1
    integer fd_l1,fd_h1
    integer pgdnv_l1,pgdnv_h1
    integer ugdnv_l1,ugdnv_h1
    integer ilo,ihi
    double precision dx, dt
    double precision     q(   qd_l1:qd_h1,QVAR)
    double precision flatn(   qd_l1:qd_h1)
    double precision  flux(fd_l1   :fd_h1,NVAR)
    double precision  srcQ(src_l1  :src_h1,QVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision ugdnv(ugdnv_l1:ugdnv_h1)
    double precision dloga(dloga_l1:dloga_h1)
    
    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    allocate ( qm(ilo-1:ihi+1,QVAR))
    allocate ( qp(ilo-1:ihi+1,QVAR))

    ! Trace to edges w/o transverse flux correction terms
    if (ppm_type .gt. 0) then
       call trace_ppm(q,q(:,QC),flatn,q(:,QGAMC),qd_l1,qd_h1, &
                      dloga,dloga_l1,dloga_h1, &
                      srcQ,src_l1,src_h1, &
                      qm,qp,ilo-1,ihi+1, &
                      ilo,ihi,domlo,domhi,dx,dt)
    else
       allocate ( dq(ilo-1:ihi+1,QVAR))

       call trace(q,dq,q(:,QC),flatn,qd_l1,qd_h1, &
                  dloga,dloga_l1,dloga_h1, &
                  srcQ,src_l1,src_h1, &
                  qm,qp,ilo-1,ihi+1, &
                  ilo,ihi,domlo,domhi,dx,dt)

       deallocate(dq)
    end if

    ! Solve Riemann problem, compute xflux from improved predicted states 
    call cmpflx(lo, hi, domlo, domhi, &
                qm, qp, ilo-1,ihi+1, &
                flux ,  fd_l1, fd_h1, &
                pgdnv,pgdnv_l1,pgdnv_h1, &
                ugdnv,ugdnv_l1,ugdnv_h1, &
                q(:,QGAMC),q(:,QCSML),q(:,QC),qd_l1,qd_h1,ilo,ihi)

    deallocate (qm,qp)

  end subroutine umeth1d

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine consup( &
                    uin,  uin_l1,  uin_h1, &
                    uout, uout_l1 ,uout_h1, &
                    update,updt_l1,updt_h1, &
                    pgdnv,pgdnv_l1,pgdnv_h1, &
                    q, q_l1, q_h1, &
                    flux, flux_l1, flux_h1, &
                    area,area_l1,area_h1, &
                    vol,vol_l1,vol_h1, &
                    div,pdivu,lo,hi,dx,dt,mass_added_flux,E_added_flux, &
                    xmom_added_flux,ymom_added_flux,zmom_added_flux, &
                    mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                    eden_lost,xang_lost,yang_lost,zang_lost, &
                    verbose)

    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, track_grid_losses, &
                                   small_dens, limit_fluxes_on_small_dens, cfl, QVAR
    use bl_constants_module
    use advection_util_1d_module, only: normalize_species_fluxes
    use advection_util_module, only : dflux
    use prob_params_module, only : domlo_level, domhi_level, center, coord_type
    use castro_util_module, only : position, linear_to_angular_momentum
    use amrinfo_module, only : amr_level

    integer lo(1), hi(1)
    integer   uin_l1,  uin_h1
    integer  uout_l1, uout_h1
    integer  updt_l1, updt_h1
    integer pgdnv_l1,pgdnv_h1
    integer     q_l1,    q_h1
    integer  flux_l1, flux_h1
    integer  area_l1, area_h1
    integer   vol_l1,  vol_h1
    integer verbose
    double precision   uin(uin_l1:uin_h1,NVAR)
    double precision  uout(uout_l1:uout_h1,NVAR)
    double precision update(updt_l1:updt_h1,NVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision     q(q_l1:q_h1,QVAR)
    double precision  flux( flux_l1: flux_h1,NVAR)
    double precision  area( area_l1: area_h1)
    double precision    vol(vol_l1:vol_h1)
    double precision    div(lo(1):hi(1)+1)
    double precision  pdivu(lo(1):hi(1)  )
    double precision thetap(lo(1):hi(1)+1)
    double precision thetam(lo(1)-1:hi(1))
    double precision dx, dt
    double precision E_added_flux, mass_added_flux
    double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
    double precision mass_lost, xmom_lost, ymom_lost, zmom_lost
    double precision eden_lost, xang_lost, yang_lost, zang_lost

    integer          :: i, j, k, n
    double precision :: div1
    integer          :: domlo(3), domhi(3)
    double precision :: loc(3), ang_mom(3)

    double precision :: rho, fluxLF(NVAR), fluxL(NVAR), fluxR(NVAR), rhoLF, drhoLF, dtdx, theta
    integer          :: dir
    logical          :: include_pressure

    do n = 1, NVAR
       if ( n == UTEMP ) then
          flux(:,n) = ZERO
       else if ( n == UMY ) then
          flux(:,n) = ZERO
       else if ( n == UMZ ) then
          flux(:,n) = ZERO
       else
          do i = lo(1),hi(1)+1
             div1 = difmag*min(ZERO,div(i))
             flux(i,n) = flux(i,n) + dx*div1*(uin(i,n) - uin(i-1,n))
          enddo
       endif
    enddo

    ! Normalize the species fluxes.

    call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

    ! Limit the fluxes to avoid negative/small densities.

    if (limit_fluxes_on_small_dens .eq. 1) then

       thetap(:) = ONE
       thetam(:) = ONE

       dir = 1

       dtdx = dt / dx

       include_pressure = .false.

       ! The following algorithm comes from Hu, Adams, and Shu (2013), JCP, 242, 169,
       ! "Positivity-preserving method for high-order conservative schemes solving
       ! compressible Euler equations." It has been modified to enforce not only positivity
       ! but also the stronger requirement that rho > small_dens.

       do i = lo(1) - 1, hi(1) + 1

          ! Note that this loop includes one ghost zone on either side of the current
          ! bounds, but we only need a one-sided limiter for lo(1)-1 and hi(1)+1.

          ! First we'll do the plus state, which is on the left edge of the zone.

          if (i .ge. lo(1)) then

             ! Obtain the one-sided update to the density, based on Hu et al., Eq. 11.
             ! Note that the sign convention for the notation is opposite to our convention
             ! for the edge states for the flux limiter, that is, the "plus" limiter is on
             ! the left edge of the zone and so is the "minus" rho. The flux limiter convention
             ! is analogous to the convention for the hydro reconstruction edge states.

             rho = uin(i,URHO) + TWO * dt * (area(i) / vol(i)) * flux(i,URHO)

             if (rho < small_dens) then

                ! Construct the Lax-Friedrichs flux on the interface (Equation 12).
                ! Note that we are using the information from Equation 9 to obtain the
                ! effective maximum wave speed, (|u| + c)_max = CFL / lambda where lambda = dt/dx.

                fluxL = dflux(uin(i-1,:), q(i-1,:), dir, [i-1, 0, 0], include_pressure)
                fluxR = dflux(uin(i  ,:), q(i  ,:), dir, [i  , 0, 0], include_pressure)
                fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx) * (uin(i-1,:) - uin(i,:)))

                ! Limit the Lax-Friedrichs flux so that it doesn't cause a density < small_dens.
                ! To do this, first, construct the density change corresponding to the LF density flux.
                ! Then, if this update would create a density that is less than small_dens, scale all
                ! fluxes linearly such that the density flux gives small_dens when applied.

                drhoLF = TWO * dt * (area(i) / vol(i)) * fluxLF(URHO)

                if (uin(i,URHO) + drhoLF < small_dens) then
                   fluxLF = fluxLF * abs((small_dens - uin(i,URHO)) / drhoLF)
                endif

                ! Obtain the final density corresponding to the LF flux.

                rhoLF = uin(i,URHO) + TWO * dt * (area(i) / vol(i)) * fluxLF(URHO)

                ! Solve for theta from (1 - theta) * rhoLF + theta * rho = small_dens.

                thetap(i) = (small_dens - rhoLF) / (rho - rhoLF)

             endif

          endif

          ! Now do the minus state, which is on the right edge of the zone.
          ! This uses the same logic as the above.

          if (i .le. hi(1)) then

             rho = uin(i,URHO) - TWO * dt * (area(i+1) / vol(i)) * flux(i+1,URHO)

             if (rho < small_dens) then

                fluxL = dflux(uin(i  ,:), q(i  ,:), dir, [i  , 0, 0], include_pressure)
                fluxR = dflux(uin(i+1,:), q(i+1,:), dir, [i+1, 0, 0], include_pressure)
                fluxLF = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx) * (uin(i,:) - uin(i+1,:)))

                drhoLF = -TWO * dt * (area(i+1) / vol(i)) * fluxLF(URHO)

                if (uin(i,URHO) + drhoLF < small_dens) then
                   fluxLF = fluxLF * abs((small_dens - uin(i,URHO)) / drhoLF)
                endif

                rhoLF = uin(i,URHO) - TWO * dt * (area(i+1) / vol(i)) * fluxLF(URHO)

                thetam(i) = (small_dens - rhoLF) / (rho - rhoLF)

             endif

          endif

       enddo

       ! Now figure out the limiting values of theta. Each zone center has a thetap and thetam,
       ! but we want a nodal value of theta that is the strongest of the two limiters in each case.
       ! Then, limit the flux accordingly.

       do i = lo(1), hi(1) + 1

          theta = min(thetam(i-1), thetap(i))

          fluxL = dflux(uin(i-1,:), q(i-1,:), dir, [i-1, 0, 0], include_pressure)
          fluxR = dflux(uin(i  ,:), q(i  ,:), dir, [i  , 0, 0], include_pressure)
          fluxLF(:) = HALF * (fluxL(:) + fluxR(:) + (cfl / dtdx) * (uin(i-1,:) - uin(i,:)))

          drhoLF = TWO * dt * (area(i) / vol(i)) * fluxLF(URHO)

          if (uin(i,URHO) + drhoLF < small_dens) then
             fluxLF = fluxLF * abs((small_dens - uin(i,URHO)) / drhoLF)
          else if (uin(i-1,URHO) - drhoLF < small_dens) then
             fluxLF = fluxLF * abs((small_dens - uin(i-1,URHO)) / drhoLF)
          endif

          flux(i,:) = (ONE - theta) * fluxLF(:) + theta * flux(i,:)

       enddo

    endif

    ! Fill the update array.

    do n = 1, NVAR
       do i = lo(1), hi(1)

          update(i,n) = update(i,n) + ( flux(i,n) * area(i) - flux(i+1,n) * area(i+1) ) / vol(i)

          ! Add p div(u) source term to (rho e).

          if (n == UEINT) then

             update(i,n) = update(i,n) - pdivu(i)

          endif

       enddo
    enddo

    ! Add gradp term to momentum equation.

    do i = lo(1),hi(1)

       update(i,UMX) = update(i,UMX) - ( pgdnv(i+1) - pgdnv(i) ) / dx

    enddo

    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (verbose .eq. 1) then

       do i = lo(1), hi(1)

          mass_added_flux = mass_added_flux + dt * ( flux(i,URHO ) - flux(i+1,URHO ) )
          xmom_added_flux = xmom_added_flux + dt * ( flux(i,UMX  ) - flux(i+1,UMX  ) )
          ymom_added_flux = ymom_added_flux + dt * ( flux(i,UMY  ) - flux(i+1,UMY  ) )
          zmom_added_flux = zmom_added_flux + dt * ( flux(i,UMZ  ) - flux(i+1,UMZ  ) )
          E_added_flux    = E_added_flux    + dt * ( flux(i,UEDEN) - flux(i+1,UEDEN) )

       enddo

    endif

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)

       j = 0
       k = 0

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)

          loc = position(i,j,k,ccx=.false.)

          mass_lost = mass_lost - dt * flux(i,URHO)
          xmom_lost = xmom_lost - dt * flux(i,UMX)
          ymom_lost = ymom_lost - dt * flux(i,UMY)
          zmom_lost = zmom_lost - dt * flux(i,UMZ)
          eden_lost = eden_lost - dt * flux(i,UEDEN)

          ang_mom   = linear_to_angular_momentum(loc - center, dt * flux(i,UMX:UMZ))
          xang_lost = xang_lost - ang_mom(1)
          yang_lost = yang_lost - ang_mom(2)
          zang_lost = zang_lost - ang_mom(3)

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1

          loc = position(i,j,k,ccx=.false.)

          mass_lost = mass_lost + dt * flux(i,URHO)
          xmom_lost = xmom_lost + dt * flux(i,UMX)
          ymom_lost = ymom_lost + dt * flux(i,UMY)
          zmom_lost = zmom_lost + dt * flux(i,UMZ)
          eden_lost = eden_lost + dt * flux(i,UEDEN)

          ang_mom   = linear_to_angular_momentum(loc - center, dt * flux(i,UMX:UMZ))
          xang_lost = xang_lost + ang_mom(1)
          yang_lost = yang_lost + ang_mom(2)
          zang_lost = zang_lost + ang_mom(3)

       endif

    endif

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do i = lo(1), hi(1)+1

          flux(i,n) = dt * area(i) * flux(i,n)

          ! Correct the momentum flux with the grad p part.

          if (coord_type .eq. 0 .and. n == UMX) then
             flux(i,n) = flux(i,n) + dt * area(i) * pgdnv(i)
          endif

       enddo
    enddo

  end subroutine consup

end module advection_module
