module rad_advection_module

  use bl_constants_module, only : ZERO, HALF, ONE

  implicit none

  private

  public umeth1d_rad, consup_rad

contains


! ::: ---------------------------------------------------------------
! ::: :: UMETH1D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: qaux        => (const)  auxillary hydro info
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux      <=  (modify) flux in X direction on X edges
! ::: ----------------------------------------------------------------

  subroutine umeth1d_rad(lo,hi,domlo,domhi, &
                         lam, lam_l1, lam_h1, &
                         q, qd_l1, qd_h1, &
                         qaux, qa_l1, qa_h1, &
                         flatn, &
                         srcQ,src_l1,src_h1, &
                         ilo,ihi,dx,dt, &
                         flux ,   fd_l1,   fd_h1, &
                         rflux,  rfd_l1,  rfd_h1, &
                         q1, q1_l1, q1_h1, &
                         dloga,dloga_l1,dloga_h1)

    use meth_params_module, only : NVAR, ppm_type, QC, QCG, QCSML, QGAMC, QGAMCG, NQAUX, NGDNV
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
    use riemann_module, only : cmpflx
    use trace_ppm_rad_module, only : trace_ppm_rad

    implicit none

    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer dloga_l1,dloga_h1
    integer lam_l1,lam_h1
    integer qd_l1,qd_h1
    integer qa_l1,qa_h1
    integer src_l1,src_h1
    integer fd_l1,fd_h1
    integer rfd_l1,rfd_h1
    integer q1_l1, q1_h1
    integer ilo,ihi
    double precision dx, dt
    double precision lam(lam_l1:lam_h1, 0:ngroups-1)
    double precision     q(   qd_l1:qd_h1,QRADVAR)
    double precision  qaux(   qd_l1:qd_h1,NQAUX)
    double precision flatn(   qd_l1:qd_h1)
    double precision  flux(fd_l1   :fd_h1,NVAR)
    double precision rflux(rfd_l1:rfd_h1, 0:ngroups-1)
    double precision  srcQ(src_l1  :src_h1,NVAR)
    double precision    q1(   q1_l1:q1_h1, NGDNV)
    double precision dloga(dloga_l1:dloga_h1)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    allocate ( dq(ilo-1:ihi+1,QRADVAR))
    allocate ( qm(ilo-1:ihi+1,QRADVAR))
    allocate ( qp(ilo-1:ihi+1,QRADVAR))

    ! Trace to edges w/o transverse flux correction terms
    if (ppm_type .gt. 0) then
       call trace_ppm_rad(lam, lam_l1, lam_h1, &
            q,qaux(:,QC),qaux(:,QCG),qaux(:,QGAMC),qaux(:,QGAMCG),flatn,qd_l1,qd_h1, &
            dloga,dloga_l1,dloga_h1, &
            srcQ,src_l1,src_h1, &
            qm,qp,ilo-1,ihi+1, &
            ilo,ihi,domlo,domhi,dx,dt)
    else
       call bl_error("ppm_type <=0 is not supported in umeth1d_rad")
       ! call trace(q,dq,c,flatn,qd_l1,qd_h1, &
       !      dloga,dloga_l1,dloga_h1, &
       !      srcQ,src_l1,src_h1, &
       !      qm,qp,ilo-1,ihi+1, &
       !      ilo,ihi,domlo,domhi,dx,dt)
    end if

    ! Solve Riemann problem, compute xflux from improved predicted states
    call cmpflx(lo, hi, domlo, domhi, &
         qm, qp, ilo-1,ihi+1, &
         flux ,  fd_l1, fd_h1, &
         q1, q1_l1, q1_h1, &
         lam, lam_l1, lam_h1, &
         rflux, rfd_l1,rfd_h1, &
         qaux(:,QGAMCG),qaux(:,QGAMC),qaux(:,QCSML),qaux(:,QC),qd_l1,qd_h1,ilo,ihi)

    deallocate (dq,qm,qp)

  end subroutine umeth1d_rad


  ! :::
  ! ::: ------------------------------------------------------------------
  ! :::

  subroutine consup_rad(uin,  uin_l1,  uin_h1, &
                        uout, uout_l1 ,uout_h1, &
                        Erin,Erin_l1,Erin_h1, &
                        Erout,Erout_l1,Erout_h1, &
                        q1,q1_l1,q1_h1, &
                        flux, flux_l1, flux_h1, &
                        rflux,rflux_l1,rflux_h1, &
                        flat, flat_l1, flat_h1, &
                        area,area_l1,area_h1, &
                        vol,vol_l1,vol_h1, &
                        div,pdivu,lo,hi,dx,dt, &
                        nstep_fsp)

    use meth_params_module, only : difmag, NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
                                   NGDNV, GDU, GDPRES, GDLAMS, GDERADS
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_params_module, only : fspace_type, comoving
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
    use advection_util_1d_module, only : normalize_species_fluxes
    use prob_params_module, only : coord_type

    implicit none
    integer nstep_fsp
    integer lo(1), hi(1)
    integer   uin_l1,  uin_h1
    integer  uout_l1, uout_h1
    integer  Erin_l1, Erin_h1
    integer Erout_l1,Erout_h1
    integer    q1_l1,   q1_h1
    integer  flux_l1, flux_h1
    integer rflux_l1,rflux_h1
    integer  flat_l1, flat_h1
    integer  area_l1, area_h1
    integer   vol_l1,  vol_h1
    double precision   uin(uin_l1:uin_h1,NVAR)
    double precision  uout(uout_l1:uout_h1,NVAR)
    double precision  Erin( Erin_l1: Erin_h1, 0:ngroups-1)
    double precision Erout(Erout_l1:Erout_h1, 0:ngroups-1)
    double precision    q1(q1_l1:q1_h1, NGDNV)
    double precision  flux( flux_l1: flux_h1,NVAR)
    double precision rflux(rflux_l1:rflux_h1, 0:ngroups-1)
    double precision  flat( flat_l1: flat_h1)
    double precision  area( area_l1: area_h1)
    double precision    vol(vol_l1:vol_h1)
    double precision    div(lo(1):hi(1)+1)
    double precision  pdivu(lo(1):hi(1)  )
    double precision dx, dt

    integer          :: i, n, g
    double precision :: div1
    double precision :: SrU,Up,SrE

    double precision, dimension(0:ngroups-1) :: Erscale
    double precision, dimension(0:ngroups-1) :: ustar, af
    double precision :: Eddf, Eddflft, Eddfrgt, f1, f2, f1lft, f1rgt
    double precision :: ux, divu, dudx, Egdc, lamc
    double precision :: dpdx, dprdx, ek1, ek2, dek

    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if

    ! Normalize the species fluxes
    call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

    do n = 1, NVAR
       if ( n == UTEMP ) then
          flux(:,n) = ZERO
       else
          do i = lo(1),hi(1)+1
             div1 = difmag*min(ZERO,div(i))
             flux(i,n) = flux(i,n) &
                  + dx*div1*(uin(i,n) - uin(i-1,n))
             flux(i,n) = area(i) * flux(i,n) * dt
          enddo
       endif
    enddo

    do g=0, ngroups-1
       do i = lo(1),hi(1)+1
          div1 = difmag*min(ZERO,div(i))
          rflux(i,g) = rflux(i,g) &
               + dx*div1*(Erin(i,g) - Erin(i-1,g))
          rflux(i,g) = area(i) * rflux(i,g) * dt
       enddo
    end do

    do n = 1, NVAR
       do i = lo(1),hi(1)
          uout(i,n) = uout(i,n) + ( flux(i,n) - flux(i+1,n) ) / vol(i)
       enddo
    enddo

    do g=0, ngroups-1
       do i = lo(1),hi(1)
          Erout(i,g) = Erin(i,g) + (rflux(i,g) - rflux(i+1,g) ) / vol(i)
       enddo
    end do

    ! Add source term to (rho e)
    do i = lo(1),hi(1)
       uout(i,UEINT) = uout(i,UEINT) - dt * pdivu(i)
    enddo

    ! Add gradp term to momentum equation
    do i = lo(1),hi(1)
       dpdx  = (  q1(i+1,GDPRES)- q1(i,GDPRES) ) / dx

       dprdx = ZERO
       do g=0,ngroups-1
          lamc = HALF*(q1(i,GDLAMS+g) + q1(i+1,GDLAMS+g))
          dprdx = dprdx + lamc*(q1(i+1,GDERADS+g) - q1(i,GDERADS+g))/dx
       end do

       uout(i,UMX) = uout(i,UMX) - dt * dpdx
       ek1 = uout(i,UMX)**2/(2.d0*uout(i,URHO))

       uout(i,UMX) = uout(i,UMX) - dt * dprdx
       ek2 = uout(i,UMX)**2/(2.d0*uout(i,URHO))

       dek = ek2-ek1

       uout(i,UEDEN) = uout(i,UEDEN) +dek
       if (.not. comoving) then ! mixed-frame (single group only)
          Erout(i,0) = Erout(i,0) - dek
       end if
    enddo

    ! Add radiation source term to rho*u, rhoE, and Er
    if (comoving) then
       do i = lo(1),hi(1)

          ux = HALF*(q1(i,GDU) + q1(i+1,GDU))

          divu = (q1(i+1,GDU)*area(i+1)-q1(i,GDU)*area(i))/vol(i)
          dudx = (q1(i+1,GDU)-q1(i,GDU))/dx

          ! Note that for single group, fspace_type is always 1
          do g=0, ngroups-1

             lamc = HALF*(q1(i,GDLAMS+g) + q1(i+1,GDLAMS+g))
             Eddf = Edd_factor(lamc)
             f1 = (ONE-Eddf)*HALF
             f2 = (3.d0*Eddf-ONE)*HALF
             af(g) = -(f1*divu + f2*dudx)

             if (fspace_type .eq. 1) then
                Eddflft = Edd_factor(q1(i,GDLAMS+g))
                f1lft = HALF*(ONE-Eddflft)
                Eddfrgt = Edd_factor(q1(i+1,GDLAMS+g))
                f1rgt = HALF*(ONE-Eddfrgt)

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

    if (coord_type .eq. 0) then
       do i = lo(1),hi(1)+1
          flux(i,UMX) = flux(i,UMX) + dt*area(i)*q1(i,GDPRES)
       enddo
    end if

  end subroutine consup_rad

end module rad_advection_module
