#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::consup_hydro(const Box& bx,
                     Array4<Real const> const shk,
                     Array4<Real> const update,
                     Array4<Real> const flux0,
                     Array4<Real const> const qx,
                     Array4<Real const> const area0,
#if AMREX_SPACEDIM >= 2
                     Array4<Real> const flux1,
                     Array4<Real const> const qy,
                     Array4<Real const> const area1,
#endif
#if AMREX_SPACEDIM == 3
                     Array4<Real> const flux2,
                     Array4<Real const> const qz,
                     Array4<Real const> const area2,
#endif
                     Array4<Real const> const vol,
                     const Real* dx, const Real dt)
{



  // For hydro, we will create an update source term that is
  // essentially the flux divergence.  This can be added with dt to
  // get the update

  AMREX_PARALLEL_FOR_4D(bx, NUM_STATE, i, j, k, n,
  {

    Real volinv = 1.0 / vol(i,j,k);

    update(i,j,k,n) = update(i,j,k,n) +
      ( flux0(i,j,k,n) * area0(i,j,k) - flux0(i+1,j,k,n) * area0(i+1,j,k) +
#if AMREX_SPACEDIM >= 2
        flux1(i,j,k,n) * area1(i,j,k) - flux1(i,j+1,k,n) * area1(i,j+1,k) +
#endif
#if AMREX_SPACEDIM == 3
        flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j,k+1,n) * area2(i,j,k+1)
#endif
        ) * volinv;

   // Add the p div(u) source term to (rho e).
    if (n == Eint) {

      Real pdu = (qx(i+1,j,k,GDPRES) + qx(i,j,k,GDPRES)) *
                 (qx(i+1,j,k,GDU) * area0(i+1,j,k) - qx(i,j,k,GDU) * area0(i,j,k));

#if AMREX_SPACEDIM >= 2
      pdu = pdu +
        (qy(i,j+1,k,GDPRES) + qy(i,j,k,GDPRES)) *
        (qy(i,j+1,k,GDV) * area1(i,j+1,k) - qy(i,j,k,GDV) * area1(i,j,k));
#endif

#if AMREX_SPACEDIM == 3
      pdu = pdu +
        (qz(i,j,k+1,GDPRES) + qz(i,j,k,GDPRES)) *
        (qz(i,j,k+1,GDW) * area2(i,j,k+1) - qz(i,j,k,GDW) * area2(i,j,k));
#endif

      pdu = 0.5 * pdu * volinv;

      update(i,j,k,n) = update(i,j,k,n) - pdu;

#ifdef SHOCK_VAR
    } else if (n == Shock) {
      update(i,j,k,USHK) = shk(i,j,k,1) / dt;
#endif

#ifndef RADIATION
    } else if (n == Xmom) {
      // Add gradp term to momentum equation -- only for axisymmetric
      // coords (and only for the radial flux).

      if (! mom_flux_has_p[1][Xmom]) {
        update(i,j,k,Xmom) = update(i,j,k,Xmom) - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx[1];
      }
#endif

#ifdef HYBRID_MOMENTUM
    } else if (n == Rmom) {
      // update the radial momentum with the hybrid advection source
      add_hybrid_advection_source(i, j, k,
                                  dt, &dx,
                                  update,
                                  qx, qy, qz);
#endif

    }


  });

}



#ifdef RADIATION
void
Castro::consup_rad(const Box& bx,
                   Array4<Real> const shk,
                   Array4<Real> const update,
                   Array4<Real> const flux0,
#if AMREX_SPACEDIM >= 2
                   Array4<Real> const flux1,
#endif
#if AMREX_SPACEDIM == 3
                   Array4<Real> const flux2,
#endif
#ifdef RADIATION
                   Array4<Real> const Erin,
                   Array4<Real> const uout,
                   Array4<Real> const Erout,
                   Array4<Real> const radflux0,
#if AMREX_SPACEDIM >= 2
                   Array4<Real> const radflux1,
#endif
#if AMREX_SPACEDIM == 3
                   Array4<Real> const radflux2,
#endif
                   int nstep_fsp, &
#endif
                   Array4<Real> const qx,
#if AMREX_SPACEDIM >= 2
                   Array4<Real> const qy,
#endif
#if AMREX_SPACEDIM == 3
                   Array4<Real> const qz,
#endif
                   Array4<Real> const area0,
#if AMREX_SPACEDIM >= 2
                   Array4<Real> const area1,
#endif
#if AMREX_SPACEDIM == 3
                   Array4<Real> const area2,
#endif
                   Array4<Real> const vol,
                   Real dx, Real dt)
{


  if (Radiation::nGroups > 1) {
    if (fspace_type == 1) {
      Erscale = dlognu;
    } else {
      Erscale = nugroup*dlognu;
    }
  }

  // radiation energy update.  For the moment, we actually update things
  // fully here, instead of creating a source term for the update
  AMREX_PARALLEL_FOR_4D(bx, Radiation::nGroups, i, j, k, g,
  {
   Erout(i,j,k,g) = Erin(i,j,k,g) + dt *
       ( radflux0(i,j,k,g) * area0(i,j,k) - radflux0(i+1,j,k,g) * area0(i+1,j,k)
#if AMREX_SPACEDIM >= 2
       + radflux1(i,j,k,g) * area1(i,j,k) - radflux1(i,j+1,k,g) * area1(i,j+1,k)
#endif
#if AMREX_SPACEDIM == 3
       + radflux2(i,j,k,g) * area2(i,j,k) - radflux2(i,j,k+1,g) * area2(i,j,k+1)
#endif
       ) / vol(i,j,k); });

  // Add gradp term to momentum equation -- only for axisymmetry coords
  // (and only for the radial flux);  also add the radiation pressure gradient
  // to the momentum for all directions

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

              // ! pgdnv from the Riemann solver is only the gas contribution,
             // not the radiation contribution.  Note that we've already included
             // the gas pressure in the momentum flux for all Cartesian coordinate
             // directions
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                dpdx = ( qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES))/ dx(1)
                update(i,j,k,UMX) = update(i,j,k,UMX) - dpdx
             endif

             // radiation contribution -- this is sum{lambda E_r}
             dprdx = ZERO
             dprdy = ZERO
             dprdz = ZERO

             do g = 0, ngroups-1
#if AMREX_SPACEDIM == 1
                lamc = HALF*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 2
                lamc = FOURTH*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
                     qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) + 
                     qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) + 
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

             // we now want to compute the change in the kinetic energy -- we
             // base this off of uout, since that has the source terms in it.
             // But note that this does not have the fluxes (since we are
             // using update)

             // note, we need to include the Z component here too (even
             // for 1- and 2-d), since rotation might be in play

             urho_new = uout(i,j,k,URHO) + dt * update(i,j,k,URHO)

             // this update includes the hydro fluxes and grad{p} from hydro
             umx_new1 = uout(i,j,k,UMX) + dt * update(i,j,k,UMX)
             umy_new1 = uout(i,j,k,UMY) + dt * update(i,j,k,UMY)
             umz_new1 = uout(i,j,k,UMZ) + dt * update(i,j,k,UMZ)

             ek1 = (umx_new1**2 + umy_new1**2 + umz_new1**2) / (TWO*urho_new)

             // now add the radiation pressure gradient
             update(i,j,k,UMX) = update(i,j,k,UMX) - dprdx
             update(i,j,k,UMY) = update(i,j,k,UMY) - dprdy
             update(i,j,k,UMZ) = update(i,j,k,UMZ) - dprdz

             umx_new2 = umx_new1 - dt * dprdx
             umy_new2 = umy_new1 - dt * dprdy
             umz_new2 = umz_new1 - dt * dprdz

             ek2 = (umx_new2**2 + umy_new2**2 + umz_new2**2) / (TWO*urho_new)

             dek = ek2 - ek1

             // the update is a source / dt, so scale accordingly
             update(i,j,k,UEDEN) = update(i,j,k,UEDEN) + dek/dt

             if (.not. comoving) then ! mixed-frame (single group only)
                Erout(i,j,k,0) = Erout(i,j,k,0) - dek
             end if

          end do
       end do
    end do

    // Add radiation source terms to rho*u, rhoE, and Er
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

                // Note that for single group, fspace_type is always 1
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
                   lamc = 0.25e0_rt*(qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g))
#endif
#if AMREX_SPACEDIM == 3
                   lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) +
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
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*ux*Gf1E(1)
                           - dt*f2*Egdc*nnColonDotGu
#endif
#if AMREX_SPACEDIM == 2
                      Egdc = 0.25e0_rt*(qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) +
                           qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g))
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2))
                           - dt*f2*Egdc*nnColonDotGu
#endif
#if AMREX_SPACEDIM == 3
                      Egdc = (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g)
                           +  qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g)
                           +  qz(i,j,k,GDERADS+g) + qz(i,j,k+1,GDERADS+g) ) / 6.e0_rt
                      Erout(i,j,k,g) = Erout(i,j,k,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)+uz*Gf1E(3))
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

#ifdef HYBRID_MOMENTUM
void
Castro::add_hybrid_advection_source(const int i, const int j, const int k,
                                    const Real dt, const Real* dx,
                                    Array4<Real> const update,
                                    Array4<Real const> const qx,
                                    Array4<Real const> const qy,
                                    Array4<Real const> const qz)
  {

   Real center[3];
   ca_get_center(center);

   auto loc = position_c(i,j,k, false, false, false);

   for (int idir = 0; idir < 3; idir++) {
     loc[idir] -= center[idir];
   }

   Real R = std::sqrt( loc[1]*loc[1] + loc[2]*loc[2] );

   update(i,j,k,Rmom) += - ( (loc[1] / R) * (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx[1] +
                             (loc[2] / R) * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / dx[2] );

  }
#endif
