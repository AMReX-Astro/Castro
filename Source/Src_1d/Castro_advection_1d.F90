module advection_module

  use bl_constants_module, only : ZERO, HALF, ONE, FOURTH

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
! ::: :: qaux        => (const)  auxillary hydro info
! ::: :: flatn       => (const)  flattening parameter
! ::: :: srcQ        => (const)  source for primitive variables
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux      <=  (modify) flux in X direction on X edges
! ::: ----------------------------------------------------------------

  subroutine umeth1d(lo, hi, domlo, domhi, &
                     q, qd_l1, qd_h1, &
                     flatn, &
                     qaux, qa_l1, qa_h1, &
                     srcQ, src_l1,src_h1, &
                     ilo, ihi, dx, dt, &
                     uout, uout_l1, uout_h1, &
                     flux, fd_l1, fd_h1, &
#ifdef RADIATION
                     lam, lam_l1, lam_h1, &
                     rflux, rfd_l1, rfd_h1, &
#endif
                     q1, q1_l1, q1_h1, &
                     dloga, dloga_l1, dloga_h1)

    use meth_params_module, only : QVAR, NQ, NVAR, &
                                   QC, QCSML, QGAMC, NQAUX, NGDNV, &
#ifdef RADIATION
                                   QGAMCG, QCG, &
#endif
                                   ppm_type, hybrid_riemann
    use riemann_module, only : cmpflx, shock
    use trace_module, only : trace
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use trace_ppm_module, only : trace_ppm
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    implicit none

    integer, intent(in) :: lo(1),hi(1)
    integer, intent(in) :: domlo(1), domhi(1)
    integer, intent(in) :: dloga_l1, dloga_h1
    integer, intent(in) :: qd_l1, qd_h1
    integer, intent(in) :: qa_l1, qa_h1
    integer, intent(in) :: src_l1, src_h1
    integer, intent(in) :: uout_l1, uout_h1
    integer, intent(in) :: fd_l1, fd_h1
#ifdef RADIATION
    integer, intent(in) :: rfd_l1, rfd_h1
    integer, intent(in) :: lam_l1, lam_h1
#endif
    integer, intent(in) :: q1_l1, q1_h1
    integer, intent(in) :: ilo, ihi

    double precision, intent(in) :: dx, dt
    double precision, intent(in) :: q(   qd_l1:qd_h1,NQ)
    double precision, intent(in) :: qaux(   qa_l1:qa_h1,NQAUX)
    double precision, intent(in) :: flatn(   qd_l1:qd_h1)
    double precision, intent(inout) :: uout(uout_l1:uout_h1,NVAR)
    double precision, intent(inout) :: flux(fd_l1:fd_h1,NVAR)
#ifdef RADIATION
    double precision, intent(inout) :: rflux(rfd_l1:rfd_h1,0:ngroups-1)
    double precision, intent(in) :: lam(lam_l1:lam_h1,0:ngroups-1)
#endif
    double precision, intent(in) :: srcQ(src_l1  :src_h1,QVAR)
    double precision, intent(inout) :: q1(q1_l1:q1_h1,NGDNV)
    double precision, intent(in) :: dloga(dloga_l1:dloga_h1)

    double precision, allocatable :: shk(:)

    integer :: i

    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dq(:,:), qm(:,:), qp(:,:)


    allocate (shk(ilo-1:ihi+1))

#ifdef SHOCK_VAR
    uout(ilo:ihi,USHK) = ZERO

    call shock(q,qd_l1,qd_h1,shk,ilo-1,ihi+1,ilo,ihi,dx)

    ! Store the shock data for future use in the burning step.
    do i = ilo, ihi
       uout(i,USHK) = shk(i)
    enddo

    ! Discard it locally if we don't need it in the hydro update.
    if (hybrid_riemann /= 1) then
       shk(:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q,qd_l1,qd_h1,shk,ilo-1,ihi+1,ilo,ihi,dx)
    else
       shk(:) = ZERO
    endif
#endif



    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    allocate ( qm(ilo-1:ihi+1,NQ))
    allocate ( qp(ilo-1:ihi+1,NQ))

    ! Trace to edges w/o transverse flux correction terms
    if (ppm_type .gt. 0) then
#ifdef RADIATION
       call trace_ppm_rad(lam, lam_l1, lam_h1, &
                          q,qaux(:,QC),qaux(:,QCG),qaux(:,QGAMC),qaux(:,QGAMCG), &
                          flatn,qd_l1,qd_h1, &
                          dloga,dloga_l1,dloga_h1, &
                          srcQ,src_l1,src_h1, &
                          qm,qp,ilo-1,ihi+1, &
                          ilo,ihi,domlo,domhi,dx,dt)
#else
       call trace_ppm(q,qaux(:,QC),flatn,qaux(:,QGAMC),qd_l1,qd_h1, &
                      dloga,dloga_l1,dloga_h1, &
                      srcQ,src_l1,src_h1, &
                      qm,qp,ilo-1,ihi+1, &
                      ilo,ihi,domlo,domhi,dx,dt)
#endif
    else
#ifdef RADIATION
       call bl_error("ppm_type <=0 is not supported in umeth1d_rad")
#else
       allocate ( dq(ilo-1:ihi+1,NQ))

       call trace(q,dq,qaux(:,QC),flatn,qd_l1,qd_h1, &
                  dloga,dloga_l1,dloga_h1, &
                  srcQ,src_l1,src_h1, &
                  qm,qp,ilo-1,ihi+1, &
                  ilo,ihi,domlo,domhi,dx,dt)

       deallocate(dq)
#endif
    end if

    ! Solve Riemann problem, compute xflux from improved predicted states
    call cmpflx(lo, hi, domlo, domhi, &
                qm, qp, ilo-1, ihi+1, &
                flux, fd_l1, fd_h1, &
                q1, q1_l1, q1_h1, &
#ifdef RADIATION
                lam, lam_l1, lam_h1, &
                rflux, rfd_l1,rfd_h1, &
                qaux(:,QGAMCG), &
#endif
                qaux(:,QGAMC),qaux(:,QCSML),qaux(:,QC), &
                qd_l1,qd_h1,ilo,ihi)

    deallocate (qm,qp)

  end subroutine umeth1d

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine consup(uin, uin_l1, uin_h1, &
                    uout, uout_l1, uout_h1, &
                    update, updt_l1, updt_h1, &
                    q, q_l1, q_h1, &
                    flux, flux_l1, flux_h1, &
                    q1, q1_l1, q1_h1, &
#ifdef RADIATION
                    Erin, Erin_l1, Erin_h1, &
                    Erout, Erout_l1, Erout_h1, &
                    rflux, rflux_l1, rflux_h1, &
                    nstep_fsp, &
#endif
                    area, area_l1, area_h1, &
                    vol, vol_l1, vol_h1, &
                    div, pdivu, lo, hi, dx, dt, &
                    mass_added_flux, E_added_flux, &
                    xmom_added_flux, ymom_added_flux, zmom_added_flux, &
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost, &
                    verbose)

    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, track_grid_losses, &
#ifdef RADIATION
                                   GDLAMS, GDERADS, &
                                   comoving, fspace_type, &
                                   GDU, GDPRES, &
#endif
                                   limit_fluxes_on_small_dens, QVAR, NGDNV, GDPRES, NQ
    use advection_util_1d_module, only: normalize_species_fluxes
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens
    use prob_params_module, only : domlo_level, domhi_level, center, coord_type
    use castro_util_module, only : position, linear_to_angular_momentum
    use amrinfo_module, only : amr_level
#ifdef RADIATION
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    integer, intent(in) :: lo(1), hi(1)
    integer, intent(in) :: uin_l1, uin_h1
    integer, intent(in) :: uout_l1, uout_h1
    integer, intent(in) :: updt_l1, updt_h1
    integer, intent(in) :: q1_l1, q1_h1
    integer, intent(in) :: q_l1,  q_h1
    integer, intent(in) :: flux_l1, flux_h1
#ifdef RADIATION
    integer, intent(in) :: Erin_l1, Erin_h1
    integer, intent(in) :: Erout_l1, Erout_h1
    integer, intent(in) :: rflux_l1, rflux_h1
    integer, intent(inout) :: nstep_fsp
#endif
    integer, intent(in) :: area_l1, area_h1
    integer, intent(in) :: vol_l1,  vol_h1
    integer, intent(in) :: verbose

    double precision, intent(in) :: uin(uin_l1:uin_h1,NVAR)
    double precision, intent(in) :: uout(uout_l1:uout_h1,NVAR)
    double precision, intent(inout) :: update(updt_l1:updt_h1,NVAR)
    double precision, intent(in) :: q1(q1_l1:q1_h1,NGDNV)
    double precision, intent(in) :: q(q_l1:q_h1,NQ)
    double precision, intent(inout) :: flux(flux_l1:flux_h1,NVAR)
#ifdef RADIATION
    double precision, intent(in) :: Erin(Erin_l1:Erin_h1,0:ngroups-1)
    double precision, intent(inout) :: Erout(Erout_l1:Erout_h1,0:ngroups-1)
    double precision, intent(inout) :: rflux(rflux_l1:rflux_h1,0:ngroups-1)
#endif
    double precision, intent(in) :: area( area_l1: area_h1)
    double precision, intent(in) :: vol(vol_l1:vol_h1)
    double precision, intent(in) :: div(lo(1):hi(1)+1)
    double precision, intent(in) :: pdivu(lo(1):hi(1)  )
    double precision, intent(in) :: dx, dt
    double precision, intent(inout) :: E_added_flux, mass_added_flux
    double precision, intent(inout) :: xmom_added_flux, ymom_added_flux, zmom_added_flux
    double precision, intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    double precision, intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    integer          :: i, j, k, n
    double precision :: div1
    integer          :: domlo(3), domhi(3)
    double precision :: loc(3), ang_mom(3)
#ifdef RADIATION
    integer :: g
    double precision :: SrU, Up, SrE

    double precision, dimension(0:ngroups-1) :: Erscale
    double precision, dimension(0:ngroups-1) :: ustar, af
    double precision :: Eddf, Eddflft, Eddfrgt, f1, f2, f1lft, f1rgt
    double precision :: ux, divu, dudx, Egdc, lamc
    double precision :: dpdx, dprdx, ek1, ek2, dek

    double precision :: urho_new, umx_new1, umx_new2
#endif

#ifdef RADIATION
    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if
#endif

    do n = 1, NVAR
       if ( n == UTEMP ) then
          flux(lo(1):hi(1)+1,n) = ZERO
#ifdef SHOCK_VAR
       else if ( n == USHK) then
          flux(lo(1):hi(1)+1,n) = ZERO
#endif
       else if ( n == UMY ) then
          flux(lo(1):hi(1)+1,n) = ZERO
       else if ( n == UMZ ) then
          flux(lo(1):hi(1)+1,n) = ZERO
       else
          ! add the artifical viscosity
          do i = lo(1),hi(1)+1
             div1 = difmag*min(ZERO,div(i))
             flux(i,n) = flux(i,n) + dx*div1*(uin(i,n) - uin(i-1,n))
          enddo
       endif
    enddo

    ! Limit the fluxes to avoid negative/small densities.

    if (limit_fluxes_on_small_dens .eq. 1) then
       call limit_hydro_fluxes_on_small_dens(uin, [uin_l1, 0, 0], [uin_h1, 0, 0], &
                                             q, [q_l1, 0, 0], [q_h1, 0, 0], &
                                             vol, [vol_l1, 0, 0], [vol_h1, 0, 0], &
                                             flux, [flux_l1, 0, 0], [flux_h1, 0, 0], &
                                             area, [area_l1, 0, 0], [area_h1, 0, 0], &
                                             [lo(1), 0, 0], [hi(1), 0, 0], dt, [dx, ZERO, ZERO])
    endif

    ! Normalize the species fluxes.
    call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

#ifdef RADIATION
    ! now do the radiation energy groups parts
    do g=0, ngroups-1
       do i = lo(1),hi(1)+1
          div1 = difmag*min(ZERO, div(i))
          rflux(i,g) = rflux(i,g) + dx*div1*(Erin(i,g) - Erin(i-1,g))
       enddo
    enddo
#endif


    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update

    do n = 1, NVAR
       do i = lo(1), hi(1)
          update(i,n) = update(i,n) + ( flux(i,n) * area(i) - &
                                        flux(i+1,n) * area(i+1) ) / vol(i)

          ! Add p div(u) source term to (rho e)
          if (n == UEINT) then
             update(i,n) = update(i,n) - pdivu(i)
          endif

       enddo
    enddo

#ifndef RADIATION
    ! Add gradp term to momentum equation -- in 1-d this is not included in the
    ! flux
    do i = lo(1),hi(1)
       update(i,UMX) = update(i,UMX) - ( q1(i+1,GDPRES) - q1(i,GDPRES) ) / dx
    enddo

#else

    ! radiation energy update.  For the moment, we actually update things
    ! fully here, instead of creating a source term for the update
    do g=0, ngroups-1
       do i = lo(1),hi(1)
          Erout(i,g) = Erin(i,g) + dt * (rflux(i,g) * area(i) - &
                                         rflux(i+1,g) * area(i+1) ) / vol(i)
       enddo
    end do


    ! Add gradp term to momentum equation -- this includes hydro and radiation terms
    ! Note: in 1-d, the pressure term is not included in the momentum flux
    do i = lo(1),hi(1)

       ! hydrodynamics contribution
       dpdx  = (  q1(i+1,GDPRES)- q1(i,GDPRES) ) / dx
       update(i,UMX) = update(i,UMX) - dpdx

       ! radiation contribution -- this is sum{lambda E_r}
       dprdx = ZERO
       do g=0,ngroups-1
          lamc = HALF*(q1(i,GDLAMS+g) + q1(i+1,GDLAMS+g))
          dprdx = dprdx + lamc*(q1(i+1,GDERADS+g) - q1(i,GDERADS+g))/dx
       end do

       ! we now want to compute the change in kinetic energy -- we
       ! base this off of uout, since that has the source terms in it.
       ! But note that this does not have the fluxes (since we are
       ! using update)
       urho_new = uout(i,URHO) + dt * update(i,URHO)

       ! this update includes the hydro fluxes and dpdr from hydro
       umx_new1 = uout(i,UMX) + dt * update(i,UMX)
       ek1 = umx_new1**2/(2.d0*urho_new)

       ! now add the radiation pressure gradient
       update(i,UMX) = update(i,UMX) - dprdx
       umx_new2 = umx_new1 - dt * dprdx
       ek2 = umx_new2**2/(2.d0*urho_new)

       dek = ek2 - ek1

       ! the update is a source / dt, so scale accordingly
       update(i,UEDEN) = update(i,UEDEN) + dek/dt

       if (.not. comoving) then ! mixed-frame (single group only)
          Erout(i,0) = Erout(i,0) - dek
       end if
    enddo

    ! Add radiation source term to rho*u, rhoE, and Er
    if (comoving) then
       do i = lo(1),hi(1)

          ux = HALF*(q1(i,GDU) + q1(i+1,GDU))

          divu = (q1(i+1,GDU)*area(i+1) - q1(i,GDU)*area(i))/vol(i)
          dudx = (q1(i+1,GDU) - q1(i,GDU))/dx

          ! Note that for single group, fspace_type is always 1
          do g=0, ngroups-1

             lamc = HALF * (q1(i,GDLAMS+g) + q1(i+1,GDLAMS+g))
             Eddf = Edd_factor(lamc)
             f1 = (ONE - Eddf)*HALF
             f2 = (3.d0*Eddf - ONE)*HALF
             af(g) = -(f1*divu + f2*dudx)

             if (fspace_type .eq. 1) then
                Eddflft = Edd_factor(q1(i,GDLAMS+g))
                f1lft = HALF * (ONE - Eddflft)
                Eddfrgt = Edd_factor(q1(i+1,GDLAMS+g))
                f1rgt = HALF * (ONE - Eddfrgt)

                Egdc = HALF*(q1(i,GDERADS+g) + q1(i+1,GDERADS+g))
                Erout(i,g) = Erout(i,g) + &
                     dt*ux*(f1rgt*q1(i+1,GDERADS+g) - f1lft*q1(i,GDERADS+g))/dx &
                     - dt*f2*Egdc*dudx
             end if

          end do

          if (ngroups.gt.1) then
             ustar = Erout(i,:) / Erscale
             call advect_in_fspace(ustar, af, dt, nstep_fsp)
             Erout(i,:) = ustar * Erscale
          end if
       end do
    end if
#endif


    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do i = lo(1), hi(1)+1

          flux(i,n) = dt * area(i) * flux(i,n)

          ! Correct the momentum flux with the grad p part.
          if (coord_type .eq. 0 .and. n == UMX) then
             flux(i,n) = flux(i,n) + dt * area(i) * q1(i,GDPRES)
          endif

       enddo
    enddo

#ifdef RADIATION
    do g = 0, ngroups-1
       do i = lo(1), hi(1)+1
          rflux(i,g) = dt * area(i) * rflux(i,g)
       enddo
    enddo
#endif

    ! Add up some diagnostic quantities. Note that we are not dividing
    ! by the cell volume.

    if (verbose .eq. 1) then

       do i = lo(1), hi(1)

          mass_added_flux = mass_added_flux + ( flux(i,URHO ) - flux(i+1,URHO ) )
          xmom_added_flux = xmom_added_flux + ( flux(i,UMX  ) - flux(i+1,UMX  ) )
          ymom_added_flux = ymom_added_flux + ( flux(i,UMY  ) - flux(i+1,UMY  ) )
          zmom_added_flux = zmom_added_flux + ( flux(i,UMZ  ) - flux(i+1,UMZ  ) )
          E_added_flux    = E_added_flux    + ( flux(i,UEDEN) - flux(i+1,UEDEN) )

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

          mass_lost = mass_lost - flux(i,URHO)
          xmom_lost = xmom_lost - flux(i,UMX)
          ymom_lost = ymom_lost - flux(i,UMY)
          zmom_lost = zmom_lost - flux(i,UMZ)
          eden_lost = eden_lost - flux(i,UEDEN)

          ang_mom   = linear_to_angular_momentum(loc - center, flux(i,UMX:UMZ))
          xang_lost = xang_lost - ang_mom(1)
          yang_lost = yang_lost - ang_mom(2)
          zang_lost = zang_lost - ang_mom(3)

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1

          loc = position(i,j,k,ccx=.false.)

          mass_lost = mass_lost + flux(i,URHO)
          xmom_lost = xmom_lost + flux(i,UMX)
          ymom_lost = ymom_lost + flux(i,UMY)
          zmom_lost = zmom_lost + flux(i,UMZ)
          eden_lost = eden_lost + flux(i,UEDEN)

          ang_mom   = linear_to_angular_momentum(loc - center, flux(i,UMX:UMZ))
          xang_lost = xang_lost + ang_mom(1)
          yang_lost = yang_lost + ang_mom(2)
          zang_lost = zang_lost + ang_mom(3)

       endif

    endif

  end subroutine consup

end module advection_module
