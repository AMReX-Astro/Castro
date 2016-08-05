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
                                   UEDEN, UEINT, UTEMP, track_grid_losses
    use bl_constants_module
    use advection_util_1d_module, only: normalize_species_fluxes
    use prob_params_module, only : domlo_level, domhi_level, center
    use castro_util_module, only : position, linear_to_angular_momentum
    use amrinfo_module, only : amr_level

    integer lo(1), hi(1)
    integer   uin_l1,  uin_h1
    integer  uout_l1, uout_h1
    integer  updt_l1, updt_h1
    integer pgdnv_l1,pgdnv_h1
    integer  flux_l1, flux_h1
    integer  area_l1, area_h1
    integer   vol_l1,  vol_h1
    integer verbose
    double precision   uin(uin_l1:uin_h1,NVAR)
    double precision  uout(uout_l1:uout_h1,NVAR)
    double precision update(updt_l1:updt_h1,NVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision  flux( flux_l1: flux_h1,NVAR)
    double precision  area( area_l1: area_h1)
    double precision    vol(vol_l1:vol_h1)
    double precision    div(lo(1):hi(1)+1)
    double precision  pdivu(lo(1):hi(1)  )
    double precision dx, dt
    double precision E_added_flux, mass_added_flux
    double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
    double precision mass_lost, xmom_lost, ymom_lost, zmom_lost
    double precision eden_lost, xang_lost, yang_lost, zang_lost
    
    integer          :: i, j, k, n
    double precision :: div1
    integer          :: domlo(3), domhi(3)
    double precision :: loc(3), ang_mom(3)
    
    ! Normalize the species fluxes.

    call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)
    
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
             flux(i,n) = flux(i,n) &
                  + dx*div1*(uin(i,n) - uin(i-1,n))
             flux(i,n) = area(i) * flux(i,n)
          enddo
       endif
    enddo

    ! Fill the update array.

    do n = 1, NVAR
       do i = lo(1), hi(1)

          update(i,n) = update(i,n) + ( flux(i,n) - flux(i+1,n) ) / vol(i)

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

          flux(i,n) = dt * flux(i,n)

          ! Correct the momentum flux with the grad p part.

          if (n == UMX) then
             flux(i,n) = flux(i,n) + dt * area(i) * pgdnv(i)
          endif

       enddo
    enddo

  end subroutine consup

end module advection_module
