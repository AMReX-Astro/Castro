! advection routines in support of the CTU unsplit advection scheme

module ctu_module

  implicit none

contains

  subroutine consup(uin, uin_lo, uin_hi, &
                    q, q_lo, q_hi, &
                    uout, uout_lo, uout_hi, &
                    update, updt_lo, updt_hi, &
                    flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                    flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                    flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                    Erin, Erin_lo, Erin_hi, &
                    Erout, Erout_lo, Erout_hi, &
                    radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                    radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                    radflux3, radflux3_lo, radflux3_hi, &
#endif
                    nstep_fsp, &
#endif
                    qx, qx_lo, qx_hi, &
#if AMREX_SPACEDIM >= 2
                    qy, qy_lo, qy_hi, &
#endif
#if AMREX_SPACEDIM == 3
                    qz, qz_lo, qz_hi, &
#endif
                    area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                    area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                    area3, area3_lo, area3_hi, &
#endif
                    vol,vol_lo,vol_hi, &
                    div, lo, hi, dx, dt, &
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost, &
                    verbose)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, NGDNV, NQ, &
                                   GDPRES, &
#ifdef RADIATION
                                   fspace_type, comoving, &
                                   GDU, GDV, GDW, GDLAMS, GDERADS, &
#endif
                                   track_grid_losses, limit_fluxes_on_small_dens
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens, normalize_species_fluxes, calc_pdivu
    use castro_util_module, only : position, linear_to_angular_momentum
    use prob_params_module, only : mom_flux_has_p, domlo_level, domhi_level, center, dg, coord_type
    use amrinfo_module, only : amr_level
#ifdef RADIATION
    use rad_params_module, only : ngroups, nugroup, dlognu
    use radhydro_nd_module, only : advect_in_fspace
    use fluxlimiter_module, only : Edd_factor
#endif
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use amrex_constants_module, only : ZERO, ONE, TWO, FOURTH, HALF

    use amrex_fort_module, only : rt => amrex_real

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) ::  uout_lo(3),  uout_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
#endif
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)

#ifdef RADIATION
    integer, intent(in) :: Erout_lo(3), Erout_hi(3)
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
#endif
    integer, intent(inout) :: nstep_fsp
#endif

    integer, intent(in) :: verbose


    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)

    real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt)        , intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt)        , intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
#endif

#if AMREX_SPACEDIM == 3
    real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt)        , intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
#endif

    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt)        , intent(in) :: div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(rt)        , intent(in) :: dx(3), dt

#ifdef RADIATION
    real(rt)          Erin(Erin_lo(1):Erin_hi(1),Erin_lo(2):Erin_hi(2),Erin_lo(3):Erin_hi(3),0:ngroups-1)
    real(rt)         Erout(Erout_lo(1):Erout_hi(1),Erout_lo(2):Erout_hi(2),Erout_lo(3):Erout_hi(3),0:ngroups-1)
    real(rt)         radflux1(radflux1_lo(1):radflux1_hi(1),radflux1_lo(2):radflux1_hi(2),radflux1_lo(3):radflux1_hi(3),0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt)         radflux2(radflux2_lo(1):radflux2_hi(1),radflux2_lo(2):radflux2_hi(2),radflux2_lo(3):radflux2_hi(3),0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt)         radflux3(radflux3_lo(1):radflux3_hi(1),radflux3_lo(2):radflux3_hi(2),radflux3_lo(3):radflux3_hi(3),0:ngroups-1)
#endif

#endif

    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    real(rt)         :: div1, volinv
    integer          :: i, j, g, k, n
    integer          :: domlo(3), domhi(3)
    real(rt)         :: loc(3), ang_mom(3)

#ifdef RADIATION
    real(rt)        , dimension(0:ngroups-1) :: Erscale
    real(rt)        , dimension(0:ngroups-1) :: ustar, af
    real(rt)         :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, Eddfzm, Eddfzp
    real(rt)         :: f1, f2, f1xm, f1xp, f1ym, f1yp, f1zm, f1zp
    real(rt)         :: Gf1E(3)
    real(rt)         :: ux, uy, uz, divu, lamc, Egdc
    real(rt)         :: dudx(3), dudy(3), dudz(3), nhat(3), GnDotu(3), nnColonDotGu
    real(rt)         :: dprdx, dprdy, dprdz, ek1, ek2, dek, dpdx
    real(rt)         :: urho_new
    real(rt)         :: umx_new1, umy_new1, umz_new1
    real(rt)         :: umx_new2, umy_new2, umz_new2
#endif
    real(rt)        , pointer:: pdivu(:,:,:)

#ifdef RADIATION
    if (ngroups .gt. 1) then
       if (fspace_type .eq. 1) then
          Erscale = dlognu
       else
          Erscale = nugroup*dlognu
       end if
    end if
#endif

    call bl_allocate(pdivu, lo, hi)

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                    qy, qy_lo, qy_hi, &
                    area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                    qz, qz_lo, qz_hi, &
                    area3, area3_lo, area3_hi, &
#endif
                    vol, vol_lo, vol_hi, &
                    dx, pdivu, lo, hi)

    do n = 1, NVAR

       if ( n == UTEMP ) then

          flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
          flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
          flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

#ifdef SHOCK_VAR
       else if ( n == USHK ) then

          flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
          flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
          flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif
#endif

       else

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   div1 = FOURTH*(div(i,j,k) + div(i,j+dg(2),k) + &
                                  div(i,j,k+dg(3)) + div(i,j+dg(2),k+dg(3)))
                   div1 = difmag*min(ZERO, div1)

                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx(1) * div1 * (uin(i,j,k,n)-uin(i-1,j,k,n))
                enddo
             enddo
          enddo

#if AMREX_SPACEDIM >= 2
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                  div(i,j,k+dg(3)) + div(i+1,j,k+dg(3)))
                   div1 = difmag*min(ZERO, div1)

                   flux2(i,j,k,n) = flux2(i,j,k,n) + dx(2) * div1 * (uin(i,j,k,n)-uin(i,j-1,k,n))
                enddo
             enddo
          enddo
#endif

#if AMREX_SPACEDIM == 3
          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                  div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag*min(ZERO, div1)

                   flux3(i,j,k,n) = flux3(i,j,k,n) + dx(3) * div1 * (uin(i,j,k,n)-uin(i,j,k-1,n))
                enddo
             enddo
          enddo
#endif

       endif

    enddo

#ifdef RADIATION
    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                div1 = FOURTH*(div(i,j,k) + div(i,j+dg(2),k) + &
                               div(i,j,k+dg(3)) + div(i,j+dg(2),k+dg(3)))
                div1 = difmag*min(ZERO, div1)

                radflux1(i,j,k,g) = radflux1(i,j,k,g) + dx(1)*div1*(Erin(i,j,k,g)-Erin(i-1,j,k,g))
             enddo
          enddo
       enddo
    enddo

#if AMREX_SPACEDIM >= 2
    do g=0,ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                               div(i,j,k+dg(3)) + div(i+1,j,k+dg(3)))
                div1 = difmag*min(ZERO, div1)

                radflux2(i,j,k,g) = radflux2(i,j,k,g) + dx(2)*div1*(Erin(i,j,k,g)-Erin(i,j-1,k,g))
             enddo
          enddo
       enddo
    enddo
#endif

#if AMREX_SPACEDIM == 3
    do g=0,ngroups-1
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                               div(i,j+1,k) + div(i+1,j+1,k))
                div1 = difmag*min(ZERO, div1)

                radflux3(i,j,k,g) = radflux3(i,j,k,g) + dx(3)*div1*(Erin(i,j,k,g)-Erin(i,j,k-1,g))
             enddo
          enddo
       enddo
    enddo
#endif
#endif

    if (limit_fluxes_on_small_dens == 1) then
       call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
                                             q,q_lo,q_hi, &
                                             vol,vol_lo,vol_hi, &
                                             flux1,flux1_lo,flux1_hi, &
                                             area1,area1_lo,area1_hi, &
#if AMREX_SPACEDIM >= 2
                                             flux2,flux2_lo,flux2_hi, &
                                             area2,area2_lo,area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                             flux3,flux3_lo,flux3_hi, &
                                             area3,area3_lo,area3_hi, &
#endif
                                             lo,hi,dt,dx)

    endif

    call normalize_species_fluxes(flux1_lo, flux1_hi, flux1, flux1_lo,flux1_hi)
#if AMREX_SPACEDIM >= 2
    call normalize_species_fluxes(flux2_lo, flux2_hi, flux2, flux2_lo,flux2_hi)
#endif
#if AMREX_SPACEDIM == 3
    call normalize_species_fluxes(flux3_lo, flux3_hi, flux3, flux3_lo, flux3_hi)
#endif

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
#if AMREX_SPACEDIM >= 2
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
#endif
#if AMREX_SPACEDIM == 3
                     + flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) &
#endif
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

#ifdef HYBRID_MOMENTUM
    call add_hybrid_advection_source(lo, hi, dt, &
                                     update, uout_lo, uout_hi, &
                                     qx, qx_lo, qx_hi, &
                                     qy, qy_lo, qy_hi, &
                                     qz, qz_lo, qz_hi)
#endif


#ifndef RADIATION
    ! Add gradp term to momentum equation -- only for axisymmetric
    ! coords (and only for the radial flux).

    if (.not. mom_flux_has_p(1)%comp(UMX)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                update(i,j,k,UMX) = update(i,j,k,UMX) - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx(1)
                !update(i,j,UMY) = update(i,j,UMY) - (pgdy(i,j+1)-pgdy(i,j)) / dy
             enddo
          enddo
       enddo
    endif

#else
    ! radiation energy update.  For the moment, we actually update things
    ! fully here, instead of creating a source term for the update
    do g = 0, ngroups-1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                Erout(i,j,k,g) = Erin(i,j,k,g) + dt * &
                     ( radflux1(i,j,k,g) * area1(i,j,k) - radflux1(i+1,j,k,g) * area1(i+1,j,k) &
#if AMREX_SPACEDIM >= 2
                     + radflux2(i,j,k,g) * area2(i,j,k) - radflux2(i,j+1,k,g) * area2(i,j+1,k) &
#endif
#if AMREX_SPACEDIM == 3
                     + radflux3(i,j,k,g) * area3(i,j,k) - radflux3(i,j,k+1,g) * area3(i,j,k+1) &
#endif
                ) / vol(i,j,k)
             enddo
          enddo
       enddo
    enddo

    ! Add gradp term to momentum equation -- only for axisymmetry coords
    ! (and only for the radial flux);  also add the radiation pressure gradient
    ! to the momentum for all directions

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! pgdnv from the Riemann solver is only the gas contribution,
             ! not the radiation contribution.  Note that we've already included
             ! the gas pressure in the momentum flux for all Cartesian coordinate
             ! directions
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                dpdx = ( qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES))/ dx(1)
                update(i,j,k,UMX) = update(i,j,k,UMX) - dpdx
             endif

             ! radiation contribution -- this is sum{lambda E_r}
             dprdx = ZERO
             dprdy = ZERO
             dprdz = ZERO

             do g = 0, ngroups-1
#if AMREX_SPACEDIM == 1
                lamc = HALF*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 2
                lamc = FOURTH*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                               qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + &
                        qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.e0_rt
#endif

                dprdx = dprdx + lamc*(qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx(1)
#if AMREX_SPACEDIM >= 2
                dprdy = dprdy + lamc*(qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx(2)
#endif
#if AMREX_SPACEDIM == 3
                dprdz = dprdz + lamc*(qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx(3)
#endif
             end do

             ! we now want to compute the change in the kinetic energy -- we
             ! base this off of uout, since that has the source terms in it.
             ! But note that this does not have the fluxes (since we are
             ! using update)

             ! note, we need to include the Z component here too (even
             ! for 1- and 2-d), since rotation might be in play

             urho_new = uout(i,j,k,URHO) + dt * update(i,j,k,URHO)

             ! this update includes the hydro fluxes and grad{p} from hydro
             umx_new1 = uout(i,j,k,UMX) + dt * update(i,j,k,UMX)
             umy_new1 = uout(i,j,k,UMY) + dt * update(i,j,k,UMY)
             umz_new1 = uout(i,j,k,UMZ) + dt * update(i,j,k,UMZ)

             ek1 = (umx_new1**2 + umy_new1**2 + umz_new1**2) / (TWO*urho_new)

             ! now add the radiation pressure gradient
             update(i,j,k,UMX) = update(i,j,k,UMX) - dprdx
             update(i,j,k,UMY) = update(i,j,k,UMY) - dprdy
             update(i,j,k,UMZ) = update(i,j,k,UMZ) - dprdz

             umx_new2 = umx_new1 - dt * dprdx
             umy_new2 = umy_new1 - dt * dprdy
             umz_new2 = umz_new1 - dt * dprdz

             ek2 = (umx_new2**2 + umy_new2**2 + umz_new2**2) / (TWO*urho_new)

             dek = ek2 - ek1

             ! the update is a source / dt, so scale accordingly
             update(i,j,k,UEDEN) = update(i,j,k,UEDEN) + dek/dt

             if (.not. comoving) then ! mixed-frame (single group only)
                Erout(i,j,k,0) = Erout(i,j,k,0) - dek
             end if

          end do
       end do
    end do

    ! Add radiation source terms to rho*u, rhoE, and Er
    if (comoving) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ux = HALF*(qx(i,j,k,GDU) + qx(i+1,j,k,GDU))
#if AMREX_SPACEDIM >= 2
                uy = HALF*(qy(i,j,k,GDV) + qy(i,j+dg(2),k,GDV))
#endif
#if AMREX_SPACEDIM == 3
                uz = HALF*(qz(i,j,k,GDW) + qz(i,j,k+dg(3),GDW))
#endif

                dudx(:) = ZERO
                dudy(:) = ZERO
                dudz(:) = ZERO

                dudx(1) = (qx(i+1,j,k,GDU) - qx(i,j,k,GDU))/dx(1)
                dudx(2) = (qx(i+1,j,k,GDV) - qx(i,j,k,GDV))/dx(1)
                dudx(3) = (qx(i+1,j,k,GDW) - qx(i,j,k,GDW))/dx(1)

#if AMREX_SPACEDIM >= 2
                dudy(1) = (qy(i,j+1,k,GDU) - qy(i,j,k,GDU))/dx(2)
                dudy(2) = (qy(i,j+1,k,GDV) - qy(i,j,k,GDV))/dx(2)
                dudy(3) = (qy(i,j+1,k,GDW) - qy(i,j,k,GDW))/dx(2)
#endif

#if AMREX_SPACEDIM == 3
                dudz(1) = (qz(i,j,k+1,GDU) - qz(i,j,k,GDU))/dx(3)
                dudz(2) = (qz(i,j,k+1,GDV) - qz(i,j,k,GDV))/dx(3)
                dudz(3) = (qz(i,j,k+1,GDW) - qz(i,j,k,GDW))/dx(3)
#endif

                divu = dudx(1) + dudy(2) + dudz(3)

                ! Note that for single group, fspace_type is always 1
                do g=0, ngroups-1

                   nhat(:) = ZERO

                   nhat(1) = (qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx(1)
#if AMREX_SPACEDIM >= 2
                   nhat(2) = (qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx(2)
#endif
#if AMREX_SPACEDIM == 3
                   nhat(3) = (qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx(3)
#endif

                   GnDotu(1) = dot_product(nhat, dudx)
                   GnDotu(2) = dot_product(nhat, dudy)
                   GnDotu(3) = dot_product(nhat, dudz)

                   nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.e-50_rt)

#if AMREX_SPACEDIM == 1
                   lamc = HALF*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 2
                   lamc = 0.25e0_rt*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                                     qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                   lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + &
                           qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + &
                           qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.e0_rt
#endif

                   Eddf = Edd_factor(lamc)
                   f1 = (ONE-Eddf)*HALF
                   f2 = (3.e0_rt*Eddf-ONE)*HALF
                   af(g) = -(f1*divu + f2*nnColonDotGu)

                   if (fspace_type .eq. 1) then
                      Eddfxp = Edd_factor(qx(i+1,j  ,k  ,GDLAMS+g))
                      Eddfxm = Edd_factor(qx(i  ,j  ,k  ,GDLAMS+g))
#if AMREX_SPACEDIM >= 2
                      Eddfyp = Edd_factor(qy(i  ,j+1,k  ,GDLAMS+g))
                      Eddfym = Edd_factor(qy(i  ,j  ,k  ,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                      Eddfzp = Edd_factor(qz(i  ,j  ,k+1,GDLAMS+g))
                      Eddfzm = Edd_factor(qz(i  ,j  ,k  ,GDLAMS+g))
#endif

                      f1xp = HALF*(ONE-Eddfxp)
                      f1xm = HALF*(ONE-Eddfxm)
#if AMREX_SPACEDIM >= 2
                      f1yp = HALF*(ONE-Eddfyp)
                      f1ym = HALF*(ONE-Eddfym)
#endif
#if AMREX_SPACEDIM == 3
                      f1zp = HALF*(ONE-Eddfzp)
                      f1zm = HALF*(ONE-Eddfzm)
#endif

                      Gf1E(1) = (f1xp*qx(i+1,j,k,GDERADS+g) - f1xm*qx(i,j,k,GDERADS+g)) / dx(1)
#if AMREX_SPACEDIM >= 2
                      Gf1E(2) = (f1yp*qy(i,j+1,k,GDERADS+g) - f1ym*qy(i,j,k,GDERADS+g)) / dx(2)
#endif
#if AMREX_SPACEDIM == 3
                      Gf1E(3) = (f1zp*qz(i,j,k+1,GDERADS+g) - f1zm*qz(i,j,k,GDERADS+g)) / dx(3)
#endif


#if AMREX_SPACEDIM == 1
                      Egdc = HALF*(qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g))
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*ux*Gf1E(1) &
                           - dt*f2*Egdc*nnColonDotGu
#endif
#if AMREX_SPACEDIM == 2
                      Egdc = 0.25e0_rt*(qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) + &
                                        qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g))
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)) &
                           - dt*f2*Egdc*nnColonDotGu
#endif
#if AMREX_SPACEDIM == 3
                      Egdc = (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) &
                           +  qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g) &
                           +  qz(i,j,k,GDERADS+g) + qz(i,j,k+1,GDERADS+g) ) / 6.e0_rt
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)+uz*Gf1E(3)) &
                           - dt*f2*Egdc*nnColonDotGu
#endif

                   end if

                end do

                if (ngroups.gt.1) then
                   ustar = Erout(i,j,k,:) / Erscale
                   call advect_in_fspace(ustar, af, dt, nstep_fsp)
                   Erout(i,j,k,:) = ustar * Erscale
                end if
             enddo
          enddo
       enddo
    endif
#endif

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
#if AMREX_SPACEDIM == 1
                ! Correct the momentum flux with the grad p part.
                if (coord_type .eq. 0 .and. n == UMX) then
                   flux1(i,j,k,n) = flux1(i,j,k,n) + dt * area1(i,j,k) * qx(i,j,k,GDPRES)
                endif
#endif
             enddo
          enddo
       enddo
    enddo

#if AMREX_SPACEDIM >= 2
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#if AMREX_SPACEDIM == 3
    do n = 1, NVAR
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#ifdef RADIATION
    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                radflux1(i,j,k,g) = dt * radflux1(i,j,k,g) * area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

#if AMREX_SPACEDIM >= 2
    do g = 0, ngroups-1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                radflux2(i,j,k,g) = dt * radflux2(i,j,k,g) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#if AMREX_SPACEDIM == 3
    do g = 0, ngroups-1
       do k = lo(3), hi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                radflux3(i,j,k,g) = dt * radflux3(i,j,k,g) * area3(i,j,k)
             enddo
          enddo
       enddo
    enddo
#endif

#endif


    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)

#if AMREX_SPACEDIM == 3
       if (lo(3) .le. domlo(3) .and. hi(3) .ge. domlo(3)) then

          k = domlo(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccz=.false.)

                mass_lost = mass_lost - flux3(i,j,k,URHO)
                xmom_lost = xmom_lost - flux3(i,j,k,UMX)
                ymom_lost = ymom_lost - flux3(i,j,k,UMY)
                zmom_lost = zmom_lost - flux3(i,j,k,UMZ)
                eden_lost = eden_lost - flux3(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux3(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(3) .le. domhi(3) .and. hi(3) .ge. domhi(3)) then

          k = domhi(3) + 1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccz=.false.)

                mass_lost = mass_lost + flux3(i,j,k,URHO)
                xmom_lost = xmom_lost + flux3(i,j,k,UMX)
                ymom_lost = ymom_lost + flux3(i,j,k,UMY)
                zmom_lost = zmom_lost + flux3(i,j,k,UMZ)
                eden_lost = eden_lost + flux3(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux3(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif
#endif

#if AMREX_SPACEDIM >= 2
       if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

          j = domlo(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost - flux2(i,j,k,URHO)
                xmom_lost = xmom_lost - flux2(i,j,k,UMX)
                ymom_lost = ymom_lost - flux2(i,j,k,UMY)
                zmom_lost = zmom_lost - flux2(i,j,k,UMZ)
                eden_lost = eden_lost - flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

          j = domhi(2) + 1
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost + flux2(i,j,k,URHO)
                xmom_lost = xmom_lost + flux2(i,j,k,UMX)
                ymom_lost = ymom_lost + flux2(i,j,k,UMY)
                zmom_lost = zmom_lost + flux2(i,j,k,UMZ)
                eden_lost = eden_lost + flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif
#endif

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost - flux1(i,j,k,URHO)
                xmom_lost = xmom_lost - flux1(i,j,k,UMX)
                ymom_lost = ymom_lost - flux1(i,j,k,UMY)
                zmom_lost = zmom_lost - flux1(i,j,k,UMZ)
                eden_lost = eden_lost - flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost + flux1(i,j,k,URHO)
                xmom_lost = xmom_lost + flux1(i,j,k,UMX)
                ymom_lost = ymom_lost + flux1(i,j,k,UMY)
                zmom_lost = zmom_lost + flux1(i,j,k,UMZ)
                eden_lost = eden_lost + flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

    endif

    call bl_deallocate(pdivu)

  end subroutine consup


  subroutine ca_ctu_update(lo, hi, is_finest_level, time, &
                           domlo, domhi, &
                           uin, uin_lo, uin_hi, &
                           uout, uout_lo, uout_hi, &
#ifdef RADIATION
                           Erin, Erin_lo, Erin_hi, &
                           Erout, Erout_lo, Erout_hi, &
#endif
                           q, q_lo, q_hi, &
                           qaux, qa_lo, qa_hi, &
                           srcQ, srQ_lo, srQ_hi, &
                           update, updt_lo, updt_hi, &
                           delta, dt, &
                           flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                           flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                           radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                           radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           radflux3, radflux3_lo, radflux3_hi, &
#endif
#endif
                           area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                           area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                           area3, area3_lo, area3_hi, &
#endif
#if AMREX_SPACEDIM <= 2
                           pradial, p_lo, p_hi, &
                           dloga, dloga_lo, dloga_hi, &
#endif
                           vol, vol_lo, vol_hi, &
                           verbose, &
#ifdef RADIATION
                           nstep_fsp, &
#endif
                           mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                           eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : NQ, QVAR, QPRES, NQAUX, NVAR, NHYP, NGDNV, UMX, GDPRES, &
#ifdef RADIATION
                                   QPTOT, &
#endif
                                   use_flattening, &
                                   first_order_hydro
    use advection_util_module, only : divu
    use amrex_constants_module, only : ZERO, ONE
    use flatten_module, only: uflatten
    use prob_params_module, only : mom_flux_has_p, dg, coord_type
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use flatten_module, only : rad_flatten
#endif
    use ctu_advection_module, only : umeth

    use amrex_fort_module, only : rt => amrex_real
    implicit none

#ifdef RADIATION
    integer, intent(inout) :: nstep_fsp
#endif
    integer, intent(in) :: is_finest_level
    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
#ifdef RADIATION
    integer, intent(in) :: Erin_lo(3), Erin_hi(3)
    integer, intent(in) :: Erout_lo(3), Erout_hi(3)
#endif
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)
    integer, intent(in) :: updt_lo(3), updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
#endif
#ifdef RADIATION
    integer, intent(in) :: radflux1_lo(3), radflux1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: radflux2_lo(3), radflux2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: radflux3_lo(3), radflux3_hi(3)
#endif
#endif
    integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#if AMREX_SPACEDIM <= 2
    integer, intent(in) :: p_lo(3), p_hi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif

    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
#ifdef RADIATION
    real(rt)        , intent(in) :: Erin(Erin_lo(1):Erin_hi(1), Erin_lo(2):Erin_hi(2), Erin_lo(3):Erin_hi(3), 0:ngroups-1)
    real(rt)        , intent(inout) :: Erout(Erout_lo(1):Erout_hi(1), Erout_lo(2):Erout_hi(2), Erout_lo(3):Erout_hi(3), 0:ngroups-1)
#endif
    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt)        , intent(in) :: srcQ(srQ_lo(1):srQ_hi(1), srQ_lo(2):srQ_hi(2), srQ_lo(3):srQ_hi(3), QVAR)
    real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
#endif
#ifdef RADIATION
    real(rt)        , intent(inout) :: radflux1(radflux1_lo(1):radflux1_hi(1), radflux1_lo(2):radflux1_hi(2), &
                                                radflux1_lo(3):radflux1_hi(3), 0:ngroups-1)
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(inout) :: radflux2(radflux2_lo(1):radflux2_hi(1), radflux2_lo(2):radflux2_hi(2), &
                                                radflux2_lo(3):radflux2_hi(3), 0:ngroups-1)
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(inout) :: radflux3(radflux3_lo(1):radflux3_hi(1), radflux3_lo(2):radflux3_hi(2), &
                                                radflux3_lo(3):radflux3_hi(3), 0:ngroups-1)
#endif
#endif
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))

#if AMREX_SPACEDIM < 3
    real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
    real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
#endif

    real(rt)        , intent(in) :: delta(3), dt, time

    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    ! Automatic arrays for workspace
    real(rt)        , pointer:: flatn(:,:,:)
    real(rt)        , pointer:: div(:,:,:)

    ! Edge-centered primitive variables (Riemann state)
    real(rt)        , pointer:: q1(:,:,:,:)
    real(rt)        , pointer:: q2(:,:,:,:)
    real(rt)        , pointer:: q3(:,:,:,:)

    integer :: ngf
    integer :: q1_lo(3), q1_hi(3), q2_lo(3), q2_hi(3), q3_lo(3), q3_hi(3)

    ngf = 1

    call bl_allocate(   div, lo, hi+dg)

    q1_lo = flux1_lo - dg
    q1_hi = flux1_hi + dg
#if AMREX_SPACEDIM >= 2
    q2_lo = flux2_lo - dg
    q2_hi = flux2_hi + dg
#endif
#if AMREX_SPACEDIM == 3
    q3_lo = flux3_lo - dg
    q3_hi = flux3_hi + dg
#endif

    call bl_allocate(q1, q1_lo, q1_hi, NGDNV)
#if AMREX_SPACEDIM >= 2
    call bl_allocate(q2, q2_lo, q2_hi, NGDNV)
#endif
#if AMREX_SPACEDIM == 3
    call bl_allocate(q3, q3_lo, q3_hi, NGDNV)
#endif

    ! Compute flattening coefficient for slope calculations.
    call bl_allocate( flatn, q_lo, q_hi)

    if (first_order_hydro == 1) then
       flatn = ZERO
    elseif (use_flattening == 1) then
#ifdef RADIATION
       call rad_flatten(lo-dg*ngf, hi+dg*ngf, &
                        q, flatn, q_lo, q_hi)
#else
       call uflatten(lo-dg*ngf, hi+dg*ngf, &
                     q, flatn, q_lo, q_hi, QPRES)
#endif
    else
       flatn = ONE
    endif

    ! Compute hyperbolic fluxes using unsplit Godunov
    call umeth(q, q_lo, q_hi, &
               flatn, &
               qaux, qa_lo, qa_hi, &
               srcQ, srQ_lo, srQ_hi, &
               lo, hi, delta, dt, &
               uout, uout_lo, uout_hi, &
               flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
               flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
               flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
               radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
               radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
               radflux3, radflux3_lo, radflux3_hi, &
#endif
#endif
               q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
               q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
               q3, q3_lo, q3_hi, &
#endif
#if AMREX_SPACEDIM < 3
                area1, area1_lo, area1_hi, &
#endif
#if AMREX_SPACEDIM == 2
                area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM < 3
                vol, vol_lo, vol_hi, &
                dloga, dloga_lo, dloga_hi, &
#endif
                domlo, domhi)


    call bl_deallocate( flatn)

    ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
    call divu(lo, hi+dg, q, q_lo, q_hi, delta, div, lo, hi+dg)

    ! Conservative update
    call consup(uin, uin_lo, uin_hi, &
                q, q_lo, q_hi, &
                uout, uout_lo, uout_hi, &
                update, updt_lo, updt_hi, &
                flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
                flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                flux3, flux3_lo, flux3_hi, &
#endif
#ifdef RADIATION
                Erin, Erin_lo, Erin_hi, &
                Erout, Erout_lo, Erout_hi, &
                radflux1, radflux1_lo, radflux1_hi, &
#if AMREX_SPACEDIM >= 2
                radflux2, radflux2_lo, radflux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                radflux3, radflux3_lo, radflux3_hi, &
#endif
                nstep_fsp, &
#endif
                q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
                q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                q3, q3_lo, q3_hi, &
#endif
                area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
                area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                area3, area3_lo, area3_hi, &
#endif
                vol, vol_lo, vol_hi, &
                div, lo, hi, delta, dt, &
                mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                eden_lost,xang_lost,yang_lost,zang_lost, &
                verbose)


#if AMREX_SPACEDIM == 1
    if (coord_type > 0) then
       pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
    end if
#endif
#if AMREX_SPACEDIM == 2
    if (.not. mom_flux_has_p(1)%comp(UMX)) then
       pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
    end if
#endif

    call bl_deallocate(   div)

    call bl_deallocate(    q1)
#if AMREX_SPACEDIM >= 2
    call bl_deallocate(    q2)
#endif
#if AMREX_SPACEDIM == 3
    call bl_deallocate(    q3)
#endif

  end subroutine ca_ctu_update

end module ctu_module
