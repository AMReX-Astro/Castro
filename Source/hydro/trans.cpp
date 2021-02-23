#include <Castro.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#include <fluxlimiter.H>
#endif

using namespace amrex;

// add the transverse flux difference in direction idir_t to the
// interface states in direction idir_n

void
Castro::trans_single(const Box& bx,
                     int idir_t, int idir_n,
                     Array4<Real const> const& qm,
                     Array4<Real> const& qmo,
                     Array4<Real const> const& qp,
                     Array4<Real> const& qpo,
                     Array4<Real const> const& qaux_arr,
                     Array4<Real const> const& flux_t,
#ifdef RADIATION
                     Array4<Real const> const& rflux_t,
#endif
                     Array4<Real const> const& q_t,
#if AMREX_SPACEDIM == 2
                     Array4<Real const> const& area_t,
                     Array4<Real const> const& vol,
#endif
                     Real hdt, Real cdtdx)
{
    // Evaluate the transverse terms for both
    // the minus and plus states.

    actual_trans_single(bx, idir_t, idir_n, -1,
                        qm, qmo,
                        qaux_arr,
                        flux_t,
#ifdef RADIATION
                        rflux_t,
#endif
                        q_t,
#if AMREX_SPACEDIM == 2
                        area_t,
                        vol,
#endif
                        hdt, cdtdx);

    actual_trans_single(bx, idir_t, idir_n, 0,
                        qp, qpo,
                        qaux_arr,
                        flux_t,
#ifdef RADIATION
                        rflux_t,
#endif
                        q_t,
#if AMREX_SPACEDIM == 2
                        area_t,
                        vol,
#endif
                        hdt, cdtdx);
}


void
Castro::actual_trans_single(const Box& bx,
                            int idir_t, int idir_n, int d,
                            Array4<Real const> const& q_arr,
                            Array4<Real> const& qo_arr,
                            Array4<Real const> const& qaux_arr,
                            Array4<Real const> const& flux_t,
#ifdef RADIATION
                            Array4<Real const> const& rflux_t,
#endif
                            Array4<Real const> const& q_t,
#if AMREX_SPACEDIM == 2
                            Array4<Real const> const& area_t,
                            Array4<Real const> const& vol,
#endif
                            Real hdt, Real cdtdx)
{

    //       qm|qp
    //         |
    // --------+--------
    //   i-1       i
    //        i-1/2
    //
    // the qp state will see the transverse flux in zone i
    // the qm state will see the transverse flux in zone i-1

    // we account for this with the 'd' variable:
    // d = 0 will do qp and d = -1 will do qm

    // idir_t is the transverse direction and we set il,jl,kl
    // and ir,jr,kr to be the face-centered indices needed for
    // the transverse flux difference

    int coord = geom.Coord();

    bool reset_density = transverse_reset_density;
    bool reset_rhoe = transverse_reset_rhoe;
    Real small_p = small_pres;

#ifdef RADIATION
    int fspace_t = Radiation::fspace_advection_type;
    int comov = Radiation::comoving;
    int limiter = Radiation::limiter;
    int closure = Radiation::closure;
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // We are handling the states at the interface of
        // (i, i+1) in the x-direction, and similarly for
        // the y- and z- directions.

        int il = i;
        int jl = j;
        int kl = k;

        int ir = i;
        int jr = j;
        int kr = k;

        // set the face indices in the transverse direction

        if (idir_t == 0) {
          ir = i+1;
          jr = j;
          kr = k;

        } else if (idir_t == 1) {
          ir = i;
          jr = j+1;
          kr = k;

        } else {
          ir = i;
          jr = j;
          kr = k+1;
        }

        // We're handling both the plus and minus states;
        // for the minus state we're shifting one zone to
        // the left in our chosen direction.

        if (idir_n == 0) {
          il += d;
          ir += d;

        } else if (idir_n == 1) {
          jl += d;
          jr += d;

        } else {
          kl += d;
          kr += d;
        }

        // Update all of the passively-advected quantities with the
        // transverse term and convert back to the primitive quantity.

#if AMREX_SPACEDIM == 2
        const Real volinv = 1.0_rt / vol(il,jl,kl);
#endif
        for (int ipassive = 0; ipassive < npassive; ipassive++) {
            int n = upassmap(ipassive);
            int nqp = qpassmap(ipassive);

#if AMREX_SPACEDIM == 2
            Real rrnew = q_arr(i,j,k,QRHO) - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,URHO) -
                                                area_t(il,jl,kl) * flux_t(il,jl,kl,URHO)) * volinv;
            Real compu = q_arr(i,j,k,QRHO) * q_arr(i,j,k,nqp) - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,n) -
                                                               area_t(il,jl,kl) * flux_t(il,jl,kl,n)) * volinv;
            qo_arr(i,j,k,nqp) = compu / rrnew;
#else
            Real rrnew = q_arr(i,j,k,QRHO) - cdtdx * (flux_t(ir,jr,kr,URHO) - flux_t(il,jl,kl,URHO));
            Real compu = q_arr(i,j,k,QRHO) * q_arr(i,j,k,nqp) - cdtdx * (flux_t(ir,jr,kr,n) - flux_t(il,jl,kl,n));
            qo_arr(i,j,k,nqp) = compu / rrnew;
#endif
        }

        Real pgp  = q_t(ir,jr,kr,GDPRES);
        Real pgm  = q_t(il,jl,kl,GDPRES);
        Real ugp  = q_t(ir,jr,kr,GDU+idir_t);
        Real ugm  = q_t(il,jl,kl,GDU+idir_t);

#ifdef RADIATION
        Real lambda[NGROUPS];
        Real ergp[NGROUPS];
        Real ergm[NGROUPS];

        for (int g = 0; g < NGROUPS; ++g) {
            lambda[g] = qaux_arr(il,jl,kl,QLAMS+g);
            ergp[g] = q_t(ir,jr,kr,GDERADS+g);
            ergm[g] = q_t(il,jl,kl,GDERADS+g);
        }
#endif

        // we need to augment our conserved system with either a p
        // equation or gammae (if we have ppm_predict_gammae = 1) to
        // be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
        Real dup = area_t(ir,jr,kr) * pgp * ugp - area_t(il,jl,kl) * pgm * ugm;
        Real du = area_t(ir,jr,kr) * ugp - area_t(il,jl,kl) * ugm;
#else
        Real dup = pgp * ugp - pgm * ugm;
        Real du = ugp - ugm;
#endif
        Real pav = 0.5_rt * (pgp + pgm);
#ifdef RADIATION
        Real uav = 0.5_rt * (ugp + ugm);
#endif

        // this is the gas gamma_1
#ifdef RADIATION
        Real gamc = qaux_arr(il,jl,kl,QGAMCG);
#else
        Real gamc = qaux_arr(il,jl,kl,QGAMC);
#endif

#ifdef RADIATION
        Real lamge[NGROUPS];
        Real luge[NGROUPS];
        Real der[NGROUPS];

        Real dmom = 0.0_rt;
        Real dre = 0.0_rt;

        for (int g = 0; g < NGROUPS; ++g) {
            lamge[g] = lambda[g] * (ergp[g] - ergm[g]);
            dmom += -cdtdx * lamge[g];
            luge[g] = uav * lamge[g];
            dre += -cdtdx * luge[g];
        }

        if (fspace_t == 1 && comov) {
            for (int g = 0; g < NGROUPS; ++g) {
                Real eddf = Edd_factor(lambda[g], limiter, closure);
                Real f1 = 0.5_rt * (1.0_rt - eddf);
                der[g] = cdtdx * uav * f1 * (ergp[g] - ergm[g]);
            }
        }
        else if (fspace_t == 2) {
#if AMREX_SPACEDIM == 2
            Real divu = (area_t(ir,jr,kr) * ugp - area_t(il,jl,kl) * ugm) * volinv;
            for (int g = 0; g < NGROUPS; g++) {
                Real eddf = Edd_factor(lambda[g], limiter, closure);
                Real f1 = 0.5_rt * (1.0_rt - eddf);
                der[g] = -hdt * f1 * 0.5_rt * (ergp[g] + ergm[g]) * divu;
            }
#else
            for (int g = 0; g < NGROUPS; g++) {
                Real eddf = Edd_factor(lambda[g], limiter, closure);
                Real f1 = 0.5_rt * (1.0_rt - eddf);
                der[g] = cdtdx * f1 * 0.5_rt * (ergp[g] + ergm[g]) * (ugm - ugp);
            }
#endif
        }
        else { // mixed frame
            for (int g = 0; g < NGROUPS; g++) {
                der[g] = cdtdx * luge[g];
            }
        }
#endif

        // Convert to conservation form
        Real rrn = q_arr(i,j,k,QRHO);
        Real run = rrn * q_arr(i,j,k,QU);
        Real rvn = rrn * q_arr(i,j,k,QV);
        Real rwn = rrn * q_arr(i,j,k,QW);
        Real ekenn = 0.5_rt * rrn * (q_arr(i,j,k,QU) * q_arr(i,j,k,QU) + q_arr(i,j,k,QV) * q_arr(i,j,k,QV) + q_arr(i,j,k,QW) * q_arr(i,j,k,QW));
        Real ren = q_arr(i,j,k,QREINT) + ekenn;
#ifdef RADIATION
        Real ern[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ern[g] = q_arr(i,j,k,QRAD+g);
        }
#endif

#if AMREX_SPACEDIM == 2
        // Add transverse predictor
        Real rrnewn = rrn - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,URHO) -
                                   area_t(il,jl,kl) * flux_t(il,jl,kl,URHO)) * volinv;

        // Note that pressure may be treated specially here, depending on
        // the geometry.  Our y-interface equation for (rho u) is:
        //
        //  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
        //
        // in cylindrical coords -- note that the p term is not
        // in a divergence for UMX in the x-direction, so there
        // are no area factors.  For this geometry, we do not
        // include p in our definition of the flux in the
        // x-direction, for we need to fix this now.
        Real runewn = run - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,UMX) -
                                   area_t(il,jl,kl) * flux_t(il,jl,kl,UMX)) * volinv;
        if (idir_t == 0 && !mom_flux_has_p(0, idir_t, coord)) {
            runewn = runewn - cdtdx * (pgp - pgm);
        }
        Real rvnewn = rvn - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,UMY) -
                                   area_t(il,jl,kl) * flux_t(il,jl,kl,UMY)) * volinv;
        Real rwnewn = rwn - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,UMZ) -
                                   area_t(il,jl,kl) * flux_t(il,jl,kl,UMZ)) * volinv;
        Real renewn = ren - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,UEDEN) -
                                   area_t(il,jl,kl) * flux_t(il,jl,kl,UEDEN)) * volinv;

#ifdef RADIATION
        if (idir_t == 0) {
            Real lamge_sum = 0.0_rt;
            for (int g = 0; g < NGROUPS; ++g)
                lamge_sum = lamge_sum + lamge[g];

            runewn = runewn - 0.5_rt * hdt * (area_t(ir,jr,kr) + area_t(il,jl,kl)) * lamge_sum * volinv;
        }
        else {
            rvnewn = rvnewn + dmom;
        }

        renewn = renewn + dre;

        Real ernewn[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ernewn[g] = ern[g] - hdt * (area_t(ir,jr,kr) * rflux_t(ir,jr,kr,g) -
                                        area_t(il,jl,kl) * rflux_t(il,jl,kl,g)) * volinv + der[g];
        }
#endif

#else
        // Add transverse predictor
        Real rrnewn = rrn - cdtdx * (flux_t(ir,jr,kr,URHO)  - flux_t(il,jl,kl,URHO));
        Real runewn = run - cdtdx * (flux_t(ir,jr,kr,UMX)   - flux_t(il,jl,kl,UMX));
        Real rvnewn = rvn - cdtdx * (flux_t(ir,jr,kr,UMY)   - flux_t(il,jl,kl,UMY));
        Real rwnewn = rwn - cdtdx * (flux_t(ir,jr,kr,UMZ)   - flux_t(il,jl,kl,UMZ));
        Real renewn = ren - cdtdx * (flux_t(ir,jr,kr,UEDEN) - flux_t(il,jl,kl,UEDEN));
#ifdef RADIATION
        runewn = runewn + dmom;
        renewn = renewn + dre;
        Real ernewn[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ernewn[g] = ern[g] - cdtdx * (rflux_t(ir,jr,kr,g) - rflux_t(il,jl,kl,g)) + der[g];
        }
#endif
#endif

        // Reset to original value if adding transverse terms made density negative
        bool reset_state = false;
        if (reset_density == 1 && rrnewn < 0.0_rt) {
            rrnewn = rrn;
            runewn = run;
            rvnewn = rvn;
            rwnewn = rwn;
            renewn = ren;
#ifdef RADIATION
            for (int g = 0; g < NGROUPS; ++g) {
                ernewn[g] = ern[g];
            }
#endif
            reset_state = true;
        }

        // Convert back to primitive form
        qo_arr(i,j,k,QRHO) = rrnewn;
        Real rhoinv = 1.0_rt / rrnewn;
        qo_arr(i,j,k,QU) = runewn * rhoinv;
        qo_arr(i,j,k,QV) = rvnewn * rhoinv;
        qo_arr(i,j,k,QW) = rwnewn * rhoinv;

        // note: we run the risk of (rho e) being negative here
        Real rhoekenn = 0.5_rt * (runewn * runewn + rvnewn * rvnewn + rwnewn * rwnewn) * rhoinv;
        qo_arr(i,j,k,QREINT) = renewn - rhoekenn;

        if (!reset_state) {
            // do the transverse terms for p, gamma, and rhoe, as necessary

            if (reset_rhoe == 1 && qo_arr(i,j,k,QREINT) <= 0.0_rt) {
                // If it is negative, reset the internal energy by
                // using the discretized expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                qo_arr(i,j,k,QREINT) = q_arr(i,j,k,QREINT) - hdt * (area_t(ir,jr,kr) * flux_t(ir,jr,kr,UEINT) -
                                                            area_t(il,jl,kl) * flux_t(il,jl,kl,UEINT) +
                                                            pav * du) * volinv;
#else
                qo_arr(i,j,k,QREINT) = q_arr(i,j,k,QREINT) - cdtdx * (flux_t(ir,jr,kr,UEINT) - flux_t(il,jl,kl,UEINT) + pav * du);
#endif
            }

            // Pretend QREINT has been fixed and transverse_use_eos != 1.
            // If we are wrong, we will fix it later.

            // Add the transverse term to the p evolution eq here.
#if AMREX_SPACEDIM == 2
            // the divergences here, dup and du, already have area factors
            Real pnewn = q_arr(i,j,k,QPRES) - hdt * (dup + pav * du * (gamc - 1.0_rt)) * volinv;
#else
            Real pnewn = q_arr(i,j,k,QPRES) - cdtdx * (dup + pav * du * (gamc - 1.0_rt));
#endif
            qo_arr(i,j,k,QPRES) = amrex::max(pnewn, small_p);

        }
        else {
            qo_arr(i,j,k,QPRES) = q_arr(i,j,k,QPRES);
        }

#ifdef RADIATION
        for (int g = 0; g < NGROUPS; ++g) {
            qo_arr(i,j,k,QRAD + g) = ernewn[g];
        }

        qo_arr(i,j,k,QPTOT) = qo_arr(i,j,k,QPRES);
        for (int g = 0; g < NGROUPS; ++g) {
            qo_arr(i,j,k,QPTOT) += lambda[g] * ernewn[g];
        }

        qo_arr(i,j,k,QREITOT) = qo_arr(i,j,k,QREINT);
        for (int g = 0; g < NGROUPS; ++g) {
            qo_arr(i,j,k,QREITOT) += qo_arr(i,j,k,QRAD + g);
        }
#endif

    });

}



void
Castro::trans_final(const Box& bx,
                    int idir_n, int idir_t1, int idir_t2,
                    Array4<Real const> const& qm,
                    Array4<Real> const& qmo,
                    Array4<Real const> const& qp,
                    Array4<Real> const& qpo,
                    Array4<Real const> const& qaux_arr,
                    Array4<Real const> const& flux_t1,
#ifdef RADIATION
                    Array4<Real const> const& rflux_t1,
#endif
                    Array4<Real const> const& flux_t2,
#ifdef RADIATION
                    Array4<Real const> const& rflux_t2,
#endif
                    Array4<Real const> const& q_t1,
                    Array4<Real const> const& q_t2,
                    Real cdtdx_t1, Real cdtdx_t2)
{
    // Evaluate the transverse terms for both
    // the minus and plus states.

    actual_trans_final(bx, idir_n, idir_t1, idir_t2, -1,
                       qm, qmo,
                       qaux_arr,
                       flux_t1,
#ifdef RADIATION
                       rflux_t1,
#endif
                       flux_t2,
#ifdef RADIATION
                       rflux_t2,
#endif
                       q_t1,
                       q_t2,
                       cdtdx_t1, cdtdx_t2);

    actual_trans_final(bx, idir_n, idir_t1, idir_t2, 0,
                       qp, qpo,
                       qaux_arr,
                       flux_t1,
#ifdef RADIATION
                       rflux_t1,
#endif
                       flux_t2,
#ifdef RADIATION
                       rflux_t2,
#endif
                       q_t1,
                       q_t2,
                       cdtdx_t1, cdtdx_t2);

}



void
Castro::actual_trans_final(const Box& bx,
                           int idir_n, int idir_t1, int idir_t2, int d,
                           Array4<Real const> const& q_arr,
                           Array4<Real> const& qo_arr,
                           Array4<Real const> const& qaux_arr,
                           Array4<Real const> const& flux_t1,
#ifdef RADIATION
                           Array4<Real const> const& rflux_t1,
#endif
                           Array4<Real const> const& flux_t2,
#ifdef RADIATION
                           Array4<Real const> const& rflux_t2,
#endif
                           Array4<Real const> const& q_t1,
                           Array4<Real const> const& q_t2,
                           Real cdtdx_t1, Real cdtdx_t2)
{

    bool reset_density = transverse_reset_density;
    bool reset_rhoe = transverse_reset_rhoe;
    Real small_p = small_pres;

#ifdef RADIATION
    int fspace_t = Radiation::fspace_advection_type;
    int comov = Radiation::comoving;
    int limiter = Radiation::limiter;
    int closure = Radiation::closure;
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // the normal state
        int iln = i;
        int jln = j;
        int kln = k;

        // the first transverse state
        int il_t1 = i;
        int jl_t1 = j;
        int kl_t1 = k;

        int ir_t1 = i;
        int jr_t1 = j;
        int kr_t1 = k;

        // the second transverse state
        int il_t2 = i;
        int jl_t2 = j;
        int kl_t2 = k;

        int ir_t2 = i;
        int jr_t2 = j;
        int kr_t2 = k;

        if (idir_n == 0) {

            // x is the normal direction

            // y is the first transverse state
            ir_t1 += d;
            jr_t1 += 1;

            // z is the second transverse state
            ir_t2 += d;
            kr_t2 += 1;

            // offset for the plus/minus state
            iln += d;
            il_t1 += d;
            il_t2 += d;

        }
        else if (idir_n == 1) {

            // y is the normal direction

            // x is the first transverse state
            ir_t1 += 1;
            jr_t1 += d;

            // z is the second transverse state
            jr_t2 += d;
            kr_t2 += 1;

            // offset for the plus/minus state
            jln += d;
            jl_t1 += d;
            jl_t2 += d;

        }
        else {

            // z is the normal direction

            // x is the first transverse state
            ir_t1 += 1;
            kr_t1 += d;

            // y is the second transverse state
            jr_t2 += 1;
            kr_t2 += d;

            // offset for the plus/minus state
            kln += d;
            kl_t1 += d;
            kl_t2 += d;

        }

        // Update all of the passively-advected quantities with the
        // transverse terms and convert back to the primitive quantity.

        for (int ipassive = 0; ipassive < npassive; ++ipassive) {
            int n = upassmap(ipassive);
            int nqp = qpassmap(ipassive);

            Real rrn = q_arr(i,j,k,QRHO);
            Real compn = rrn * q_arr(i,j,k,nqp);
            Real rrnewn = rrn - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,URHO) -
                                            flux_t1(il_t1,jl_t1,kl_t1,URHO))
                              - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,URHO) -
                                            flux_t2(il_t2,jl_t2,kl_t2,URHO));
            Real compnn = compn - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,n) -
                                              flux_t1(il_t1,jl_t1,kl_t1,n))
                                - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,n) -
                                              flux_t2(il_t2,jl_t2,kl_t2,n));

            qo_arr(i,j,k,nqp) = compnn / rrnewn;
        }

        // Add the transverse differences to the normal states for the
        // fluid variables.

        Real pgt1p  = q_t1(ir_t1,jr_t1,kr_t1,GDPRES);
        Real pgt1m  = q_t1(il_t1,jl_t1,kl_t1,GDPRES);
        Real ugt1p  = q_t1(ir_t1,jr_t1,kr_t1,GDU+idir_t1);
        Real ugt1m  = q_t1(il_t1,jl_t1,kl_t1,GDU+idir_t1);
#ifdef RADIATION
        Real ergt1p[NGROUPS];
        Real ergt1m[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ergt1p[g] = q_t1(ir_t1,jr_t1,kr_t1,GDERADS+g);
            ergt1m[g] = q_t1(il_t1,jl_t1,kl_t1,GDERADS+g);
        }
#endif

        Real pgt2p  = q_t2(ir_t2,jr_t2,kr_t2,GDPRES);
        Real pgt2m  = q_t2(il_t2,jl_t2,kl_t2,GDPRES);
        Real ugt2p  = q_t2(ir_t2,jr_t2,kr_t2,GDU+idir_t2);
        Real ugt2m  = q_t2(il_t2,jl_t2,kl_t2,GDU+idir_t2);
#ifdef RADIATION
        Real ergt2p[NGROUPS];
        Real ergt2m[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ergt2p[g] = q_t2(ir_t2,jr_t2,kr_t2,GDERADS+g);
            ergt2m[g] = q_t2(il_t2,jl_t2,kl_t2,GDERADS+g);
        }
#endif

        Real dupt1 = pgt1p * ugt1p - pgt1m * ugt1m;
        Real pt1av = 0.5_rt * (pgt1p + pgt1m);
        Real dut1 = ugt1p - ugt1m;
#ifdef RADIATION
        Real pt1new = cdtdx_t1 * (dupt1 + pt1av * dut1 * (qaux_arr(iln,jln,kln,QGAMCG) - 1.0_rt));
#else
        Real pt1new = cdtdx_t1 * (dupt1 + pt1av * dut1 * (qaux_arr(iln,jln,kln,QGAMC) - 1.0_rt));
#endif

        Real dupt2 = pgt2p * ugt2p - pgt2m * ugt2m;
        Real pt2av = 0.5_rt * (pgt2p + pgt2m);
        Real dut2 = ugt2p - ugt2m;
#ifdef RADIATION
        Real pt2new = cdtdx_t2 * (dupt2 + pt2av * dut2 * (qaux_arr(iln,jln,kln,QGAMCG) - 1.0_rt));
#else
        Real pt2new = cdtdx_t2 * (dupt2 + pt2av * dut2 * (qaux_arr(iln,jln,kln,QGAMC) - 1.0_rt));
#endif

#ifdef RADIATION
        Real lambda[NGROUPS];
        Real lget1[NGROUPS];
        Real lget2[NGROUPS];
        Real dmt1 = 0.0_rt;
        Real dmt2 = 0.0_rt;
        Real luget1[NGROUPS];
        Real luget2[NGROUPS];
        Real dre = 0.0_rt;

        for (int g = 0; g < NGROUPS; ++g) {
            lambda[g] = qaux_arr(iln,jln,kln,QLAMS+g);
            lget1[g] = lambda[g] * (ergt1p[g] - ergt1m[g]);
            lget2[g] = lambda[g] * (ergt2p[g] - ergt2m[g]);
            dmt1 -= cdtdx_t1 * lget1[g];
            dmt2 -= cdtdx_t2 * lget2[g];
            luget1[g] = 0.5_rt * (ugt1p + ugt1m) * lget1[g];
            luget2[g] = 0.5_rt * (ugt2p + ugt2m) * lget2[g];
            dre -= cdtdx_t1 * luget1[g] + cdtdx_t2 * luget2[g];
        }

        Real der[NGROUPS];

        if (fspace_t == 1 && comov) {
            for (int g = 0; g < NGROUPS; ++g) {
                Real eddf = Edd_factor(lambda[g], limiter, closure);
                Real f1 = 0.5_rt * (1.0_rt - eddf);
                der[g] = f1 * (cdtdx_t1 * 0.5_rt * (ugt1p + ugt1m) * (ergt1p[g] - ergt1m[g]) +
                               cdtdx_t2 * 0.5_rt * (ugt2p + ugt2m) * (ergt2p[g] - ergt2m[g]));
            }
        }
        else if (fspace_t == 2) {
            for (int g = 0; g < NGROUPS; ++g) {
                Real eddf = Edd_factor(lambda[g], limiter, closure);
                Real f1 = 0.5_rt * (1.0_rt - eddf);
                der[g] = f1 * (cdtdx_t1 * 0.5_rt * (ergt1p[g] + ergt1m[g]) * (ugt1m - ugt1p) +
                               cdtdx_t2 * 0.5_rt * (ergt2p[g] + ergt2m[g]) * (ugt2m - ugt2p));
            }
        }
        else { // mixed frame
            for (int g = 0; g < NGROUPS; ++g) {
                der[g] = cdtdx_t1 * luget1[g] + cdtdx_t2 * luget2[g];
            }
        }
#endif

        // Convert to conservation form
        Real rrn = q_arr(i,j,k,QRHO);
        Real run = rrn * q_arr(i,j,k,QU);
        Real rvn = rrn * q_arr(i,j,k,QV);
        Real rwn = rrn * q_arr(i,j,k,QW);
        Real ekenn = 0.5_rt * rrn * (q_arr(i,j,k,QU) * q_arr(i,j,k,QU) + q_arr(i,j,k,QV) * q_arr(i,j,k,QV) + q_arr(i,j,k,QW) * q_arr(i,j,k,QW));
        Real ren = q_arr(i,j,k,QREINT) + ekenn;
#ifdef RADIATION
        Real ern[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ern[g] = q_arr(i,j,k,QRAD+g);
        }
#endif

        // Add transverse predictor
        Real rrnewn = rrn - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,URHO) -
                                        flux_t1(il_t1,jl_t1,kl_t1,URHO))
                          - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,URHO) -
                                        flux_t2(il_t2,jl_t2,kl_t2,URHO));
        Real runewn = run - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,UMX) -
                                        flux_t1(il_t1,jl_t1,kl_t1,UMX))
                          - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,UMX) -
                                        flux_t2(il_t2,jl_t2,kl_t2,UMX));
        Real rvnewn = rvn - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,UMY) -
                                        flux_t1(il_t1,jl_t1,kl_t1,UMY))
                          - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,UMY) -
                                        flux_t2(il_t2,jl_t2,kl_t2,UMY));
        Real rwnewn = rwn - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,UMZ) -
                                        flux_t1(il_t1,jl_t1,kl_t1,UMZ))
                          - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,UMZ) -
                                        flux_t2(il_t2,jl_t2,kl_t2,UMZ));
        Real renewn = ren - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,UEDEN) -
                                        flux_t1(il_t1,jl_t1,kl_t1,UEDEN))
                          - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,UEDEN) -
                                        flux_t2(il_t2,jl_t2,kl_t2,UEDEN));
#ifdef RADIATION
        if (idir_n == 0) {
            rvnewn = rvnewn + dmt1;
            rwnewn = rwnewn + dmt2;
        }
        else if (idir_n == 1) {
            runewn = runewn + dmt1;
            rwnewn = rwnewn + dmt2;
        }
        else {
            runewn = runewn + dmt1;
            rvnewn = rvnewn + dmt2;
        }
        renewn = renewn + dre;

        Real ernewn[NGROUPS];
        for (int g = 0; g < NGROUPS; ++g) {
            ernewn[g] = ern[g] - cdtdx_t1 * (rflux_t1(ir_t1,jr_t1,kr_t1,g) -
                                             rflux_t1(il_t1,jl_t1,kl_t1,g))
                               - cdtdx_t2 * (rflux_t2(ir_t2,jr_t2,kr_t2,g) -
                                             rflux_t2(il_t2,jl_t2,kl_t2,g))
                               + der[g];
        }
#endif

        // Reset to original value if adding transverse terms
        // made density negative
        bool reset_state = false;
        if (reset_density == 1 && rrnewn < 0.0_rt) {
            rrnewn = rrn;
            runewn = run;
            rvnewn = rvn;
            rwnewn = rwn;
            renewn = ren;
#ifdef RADIATION
            for (int g = 0; g < NGROUPS; ++g) {
                ernewn[g] = ern[g];
            }
#endif
            reset_state = true;
        }

        qo_arr(i,j,k,QRHO) = rrnewn;
        qo_arr(i,j,k,QU) = runewn / rrnewn;
        qo_arr(i,j,k,QV) = rvnewn / rrnewn;
        qo_arr(i,j,k,QW) = rwnewn / rrnewn;

        // note: we run the risk of (rho e) being negative here
        Real rhoekenn = 0.5_rt * (runewn * runewn + rvnewn * rvnewn + rwnewn * rwnewn) / rrnewn;
        qo_arr(i,j,k,QREINT) = renewn - rhoekenn;

        if (!reset_state) {
            if (reset_rhoe == 1 && qo_arr(i,j,k,QREINT) <= 0.0_rt) {
                // If it is negative, reset the internal energy by
                // using the discretized expression for updating
                // (rho e).
                qo_arr(i,j,k,QREINT) = q_arr(i,j,k,QREINT)
                                 - cdtdx_t1 * (flux_t1(ir_t1,jr_t1,kr_t1,UEINT) -
                                               flux_t1(il_t1,jl_t1,kl_t1,UEINT) + pt1av * dut1)
                                 - cdtdx_t2 * (flux_t2(ir_t2,jr_t2,kr_t2,UEINT) -
                                               flux_t2(il_t2,jl_t2,kl_t2,UEINT) + pt2av * dut2);
            }

            // Pretend QREINT has been fixed and transverse_use_eos != 1.
            // If we are wrong, we will fix it later.

            // add the transverse term to the p evolution eq here
            Real pnewn = q_arr(i,j,k,QPRES) - pt1new - pt2new;
            qo_arr(i,j,k,QPRES) = pnewn;
        }
        else {
            qo_arr(i,j,k,QPRES) = q_arr(i,j,k,QPRES);
        }

        qo_arr(i,j,k,QPRES) = amrex::max(qo_arr(i,j,k,QPRES), small_p);

#ifdef RADIATION
        for (int g = 0; g < NGROUPS; ++g) {
            qo_arr(i,j,k,QRAD+g) = ernewn[g];
        }

        qo_arr(i,j,k,QPTOT) = qo_arr(i,j,k,QPRES);
        for (int g = 0; g < NGROUPS; ++g) {
            qo_arr(i,j,k,QPTOT) += lambda[g] * ernewn[g];
        }

        qo_arr(i,j,k,QREITOT) = qo_arr(i,j,k,QREINT);
        for (int g = 0; g < NGROUPS; ++g) {
            qo_arr(i,j,k,QREITOT) += qo_arr(i,j,k,QRAD+g);
        }
#endif

    });

}
