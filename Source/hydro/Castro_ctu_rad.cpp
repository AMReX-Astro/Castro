#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"
#include "Castro_hydro.H"

#include "Radiation.H"
#include "RadHydro.H"
#include "RAD_F.H"
#include "fluxlimiter.H"

using namespace amrex;

void
Castro::ctu_rad_consup(const Box& bx,
                       Array4<Real> const& update,
                       Array4<Real const> const& Erin,
                       Array4<Real> const& uout,
                       Array4<Real> const& Erout,
                       Array4<Real const> const& radflux1,
                       Array4<Real const> const& qx,
                       Array4<Real const> const& area1,
#if AMREX_SPACEDIM >= 2
                       Array4<Real const> const& radflux2,
                       Array4<Real const> const& qy,
                       Array4<Real const> const& area2,
#endif
#if AMREX_SPACEDIM == 3
                       Array4<Real const> const& radflux3,
                       Array4<Real const> const& qz,
                       Array4<Real const> const& area3,
#endif
                       int& nstep_fsp,
                       Array4<Real const> const& vol,
                       const Real dt)
{


  auto dx = geom.CellSizeArray();
  const int coord_type = geom.Coord();

  GpuArray<Real, NGROUPS> Erscale = {0.0};

  int fspace_type = Radiation::fspace_advection_type;

  GpuArray<Real, NGROUPS> dlognu = {0.0};

  GpuArray<Real, NGROUPS> nugroup = {0.0};



  if (NGROUPS > 1) {
    ca_get_nugroup(nugroup.begin());
    ca_get_dlognu(dlognu.begin());

    if (fspace_type == 1) {
      for (int g = 0; g < NGROUPS; g++) {
        Erscale[g] = dlognu[g];
      }
    } else {
      for (int g = 0; g < NGROUPS; g++) {
        Erscale[g] = nugroup[g] * dlognu[g];
      }
    }
  }


  int comov = Radiation::comoving;
  int limiter = Radiation::limiter;
  int closure = Radiation::closure;

  // radiation energy update.  For the moment, we actually update things
  // fully here, instead of creating a source term for the update
  amrex::ParallelFor(bx, NGROUPS,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int g)
  {

    Erout(i,j,k,g) = Erin(i,j,k,g) + dt *
      (  radflux1(i,j,k,g) * area1(i,j,k) - radflux1(i+1,j,k,g) * area1(i+1,j,k)
#if AMREX_SPACEDIM >= 2
       + radflux2(i,j,k,g) * area2(i,j,k) - radflux2(i,j+1,k,g) * area2(i,j+1,k)
#endif
#if AMREX_SPACEDIM == 3
       + radflux3(i,j,k,g) * area3(i,j,k) - radflux3(i,j,k+1,g) * area3(i,j,k+1)
#endif
         ) / vol(i,j,k);
  });

  // add the radiation pressure gradient to the momentum for all
  // directions

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    // radiation contribution -- this is sum{lambda E_r}
    Real dprdx = 0.0;
    Real dprdy = 0.0;
    Real dprdz = 0.0;

    Real lamc = 0.0;
    for (int g = 0; g < NGROUPS; g++) {

#if AMREX_SPACEDIM == 1
      lamc = 0.5_rt * (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g));
#endif
#if AMREX_SPACEDIM == 2
      lamc = 0.25_rt * (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
                        qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g));
#endif
#if AMREX_SPACEDIM == 3
      lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
              qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) +
              qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g)) / 6.0_rt;
#endif

      dprdx = dprdx + lamc * (qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx[0];
#if AMREX_SPACEDIM >= 2
      dprdy = dprdy + lamc * (qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx[1];
#endif
#if AMREX_SPACEDIM == 3
      dprdz = dprdz + lamc * (qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx[2];
#endif
    }

    // we now want to compute the change in the kinetic energy -- we
    // base this off of uout, since that has the source terms in it.
    // But note that this does not have the fluxes (since we are
    // using update)

    // note, we need to include the Z component here too (even
    // for 1- and 2-d), since rotation might be in play

    Real urho_new = uout(i,j,k,URHO) + dt * update(i,j,k,URHO);

    // this update includes the hydro fluxes and grad{p} from hydro
    Real umx_new1 = uout(i,j,k,UMX) + dt * update(i,j,k,UMX);
    Real umy_new1 = uout(i,j,k,UMY) + dt * update(i,j,k,UMY);
    Real umz_new1 = uout(i,j,k,UMZ) + dt * update(i,j,k,UMZ);

    Real ek1 = (umx_new1 * umx_new1 +
                umy_new1 * umy_new1 +
                umz_new1 * umz_new1) / (2.0_rt * urho_new);

    // now add the radiation pressure gradient
    update(i,j,k,UMX) = update(i,j,k,UMX) - dprdx;
    update(i,j,k,UMY) = update(i,j,k,UMY) - dprdy;
    update(i,j,k,UMZ) = update(i,j,k,UMZ) - dprdz;

    Real umx_new2 = umx_new1 - dt * dprdx;
    Real umy_new2 = umy_new1 - dt * dprdy;
    Real umz_new2 = umz_new1 - dt * dprdz;

    Real ek2 = (umx_new2 * umx_new2 +
                umy_new2 * umy_new2 +
                umz_new2 * umz_new2) / (2.0_rt * urho_new);

    Real dek = ek2 - ek1;

    // the update is a source / dt, so scale accordingly
    update(i,j,k,UEDEN) = update(i,j,k,UEDEN) + dek/dt;

    if (! comov) {
      // ! mixed-frame (single group only)
      Erout(i,j,k,0) = Erout(i,j,k,0) - dek;
    }

  });


  // Add radiation source terms to rho*u, rhoE, and Er

  if (comov) {

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    reduce_op.eval(bx, reduce_data,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
    {

      Real ux = 0.5_rt * (qx(i,j,k,GDU) + qx(i+1,j,k,GDU));
#if AMREX_SPACEDIM >= 2
      Real uy = 0.5_rt * (qy(i,j,k,GDV) + qy(i,j+1,k,GDV));
#endif
#if AMREX_SPACEDIM == 3
      Real uz = 0.5_rt * (qz(i,j,k,GDW) + qz(i,j,k+1,GDW));
#endif

      Real dudx[3] = {0, 0, 0};
      Real dudy[3] = {0, 0, 0};
      Real dudz[3] = {0, 0, 0};

      dudx[0] = (qx(i+1,j,k,GDU) - qx(i,j,k,GDU))/dx[0];
      dudx[1] = (qx(i+1,j,k,GDV) - qx(i,j,k,GDV))/dx[0];
      dudx[2] = (qx(i+1,j,k,GDW) - qx(i,j,k,GDW))/dx[0];

#if AMREX_SPACEDIM >= 2
      dudy[0] = (qy(i,j+1,k,GDU) - qy(i,j,k,GDU))/dx[1];
      dudy[1] = (qy(i,j+1,k,GDV) - qy(i,j,k,GDV))/dx[1];
      dudy[2] = (qy(i,j+1,k,GDW) - qy(i,j,k,GDW))/dx[1];
#endif

#if AMREX_SPACEDIM == 3
      dudz[0] = (qz(i,j,k+1,GDU) - qz(i,j,k,GDU))/dx[2];
      dudz[1] = (qz(i,j,k+1,GDV) - qz(i,j,k,GDV))/dx[2];
      dudz[2] = (qz(i,j,k+1,GDW) - qz(i,j,k,GDW))/dx[2];
#endif

      Real divu = dudx[0] + dudy[1] + dudz[2];

      // Note that for single group, fspace_type is always 1
      Real af[NGROUPS];

      for (int g = 0; g < NGROUPS; g++) {

        Real nhat[3] = {0., 0., 0.};

        nhat[0] = (qx(i+1,j,k,GDERADS+g) - qx(i,j,k,GDERADS+g))/dx[0];
#if AMREX_SPACEDIM >= 2
        nhat[1] = (qy(i,j+1,k,GDERADS+g) - qy(i,j,k,GDERADS+g))/dx[1];
#endif
#if AMREX_SPACEDIM == 3
        nhat[2] = (qz(i,j,k+1,GDERADS+g) - qz(i,j,k,GDERADS+g))/dx[2];
#endif

        Real GnDotu[3];

        GnDotu[0] = nhat[0] * dudx[0] + nhat[1] * dudx[1] + nhat[2] * dudx[2];
        GnDotu[1] = nhat[0] * dudy[0] + nhat[1] * dudy[1] + nhat[2] * dudy[2];
        GnDotu[2] = nhat[0] * dudz[0] + nhat[1] * dudz[1] + nhat[2] * dudz[2];

        Real nnColonDotGu = (nhat[0] * GnDotu[0] + nhat[1] * GnDotu[1] + nhat[2] * GnDotu[2]) /
          (nhat[0] * nhat[0] + nhat[1] * nhat[1] + nhat[2] * nhat[2] + 1.e-50_rt);

        Real lamc;
#if AMREX_SPACEDIM == 1
        lamc = 0.5_rt * (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g));
#endif
#if AMREX_SPACEDIM == 2
        lamc = 0.25_rt * (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
                          qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g));
#endif
#if AMREX_SPACEDIM == 3
        lamc = (qx(i,j,k,GDLAMS+g) + qx(i+1,j,k,GDLAMS+g) +
                qy(i,j,k,GDLAMS+g) + qy(i,j+1,k,GDLAMS+g) +
                qz(i,j,k,GDLAMS+g) + qz(i,j,k+1,GDLAMS+g) ) / 6.0_rt;
#endif

        Real Eddf = Edd_factor(lamc, limiter, closure);
        Real f1 = (1.0_rt - Eddf) * 0.5_rt;
        Real f2 = (3.0_rt * Eddf - 1.0_rt) * 0.5_rt;
        af[g] = -(f1*divu + f2*nnColonDotGu);

        if (fspace_type == 1) {
          Real Eddfxp = Edd_factor(qx(i+1,j,k,GDLAMS+g), limiter, closure);
          Real Eddfxm = Edd_factor(qx(i,j,k,GDLAMS+g), limiter, closure);
#if AMREX_SPACEDIM >= 2
          Real Eddfyp = Edd_factor(qy(i,j+1,k,GDLAMS+g), limiter, closure);
          Real Eddfym = Edd_factor(qy(i,j,k,GDLAMS+g), limiter, closure);
#endif
#if AMREX_SPACEDIM == 3
          Real Eddfzp = Edd_factor(qz(i,j,k+1,GDLAMS+g), limiter, closure);
          Real Eddfzm = Edd_factor(qz(i,j,k,GDLAMS+g), limiter, closure);
#endif

          Real f1xp = 0.5_rt * (1.0_rt - Eddfxp);
          Real f1xm = 0.5_rt * (1.0_rt - Eddfxm);
#if AMREX_SPACEDIM >= 2
          Real f1yp = 0.5_rt * (1.0_rt - Eddfyp);
          Real f1ym = 0.5_rt * (1.0_rt - Eddfym);
#endif
#if AMREX_SPACEDIM == 3
          Real f1zp = 0.5_rt * (1.0_rt - Eddfzp);
          Real f1zm = 0.5_rt * (1.0_rt - Eddfzm);
#endif

          Real Gf1E[3];

          Gf1E[0] = (f1xp*qx(i+1,j,k,GDERADS+g) - f1xm*qx(i,j,k,GDERADS+g)) / dx[0];
#if AMREX_SPACEDIM >= 2
          Gf1E[1] = (f1yp*qy(i,j+1,k,GDERADS+g) - f1ym*qy(i,j,k,GDERADS+g)) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
          Gf1E[2] = (f1zp*qz(i,j,k+1,GDERADS+g) - f1zm*qz(i,j,k,GDERADS+g)) / dx[2];
#endif


#if AMREX_SPACEDIM == 1
          Real Egdc = 0.5_rt * (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g));
          Erout(i,j,k,g) = Erout(i,j,k,g) + dt * ux * Gf1E[0]
            - dt * f2 * Egdc * nnColonDotGu;
#endif
#if AMREX_SPACEDIM == 2
          Real Egdc = 0.25_rt * (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) +
                                 qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g));
          Erout(i,j,k,g) = Erout(i,j,k,g) + dt * (ux * Gf1E[0] + uy * Gf1E[1])
            - dt * f2 * Egdc * nnColonDotGu;
#endif
#if AMREX_SPACEDIM == 3
          Real Egdc = (qx(i,j,k,GDERADS+g) + qx(i+1,j,k,GDERADS+g) +
                       qy(i,j,k,GDERADS+g) + qy(i,j+1,k,GDERADS+g) +
                       qz(i,j,k,GDERADS+g) + qz(i,j,k+1,GDERADS+g) ) / 6.0_rt;
          Erout(i,j,k,g) = Erout(i,j,k,g) + dt * (ux * Gf1E[0] + uy * Gf1E[1] + uz * Gf1E[2])
            - dt * f2 * Egdc * nnColonDotGu;
#endif

        }
      }

      int nstep_fsp_tmp = 1;

      if (NGROUPS > 1) {
        Real ustar[NGROUPS];
        for (int g = 0; g < NGROUPS; g++) { 
          ustar[g] = Erout(i,j,k,g) / Erscale[g];
        }

        update_one_species(NGROUPS, ustar, af, dlognu.begin(), dt, nstep_fsp_tmp);

        for (int g = 0; g < NGROUPS; g++) {
          Erout(i,j,k,g) = ustar[g] * Erscale[g];
        }
      }

      return nstep_fsp_tmp;

    });

    ReduceTuple hv = reduce_data.value();
    nstep_fsp = amrex::get<0>(hv);

  }
}

