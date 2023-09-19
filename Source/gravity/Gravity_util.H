#ifndef GRAVITY_UTIL_H
#define GRAVITY_UTIL_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real factorial (int n)
{
    Real fact = 1.0_rt;

    for (int i = 2; i <= n; ++i) {
        fact *= static_cast<Real>(i);
    }

    return fact;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void calcAssocLegPolyLM (int l, int m, Real& assocLegPolyLM, Real& assocLegPolyLM1, Real& assocLegPolyLM2, Real x)
{
    // Calculate the associated Legendre polynomials. There are a number of
    // recurrence relations, but many are unstable. We'll use one that is known
    // to be stable for the reasonably low values of l we care about in a simulation:
    // (l-m)P_l^m(x) = x(2l-1)P_{l-1}^m(x) - (l+m-1)P_{l-2}^m(x).
    // This uses the following two expressions as initial conditions:
    // P_m^m(x) = (-1)^m (2m-1)! (1-x^2)^(m/2)
    // P_{m+1}^m(x) = x (2m+1) P_m^m (x)

    if (l == m) {

        // P_m^m

        assocLegPolyLM = std::pow(-1, m) * std::pow((1.0_rt - x) * (1.0_rt + x), m * 0.5_rt);

        for (int n = (2*m-1); n >= 3; n = n - 2) {

            assocLegPolyLM *= n;

        }

    }
    else if (l == m + 1) {

        // P_{m+1}^m

        assocLegPolyLM1 = assocLegPolyLM;
        assocLegPolyLM  = x * (2*m + 1) * assocLegPolyLM1;

    }
    else {

        assocLegPolyLM2 = assocLegPolyLM1;
        assocLegPolyLM1 = assocLegPolyLM;
        assocLegPolyLM  = (x * (2*l - 1) * assocLegPolyLM1 - (l + m - 1) * assocLegPolyLM2) / (l-m);

    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void calcLegPolyL (int l, Real& legPolyL, Real& legPolyL1, Real& legPolyL2, Real x)
{
    // Calculate the Legendre polynomials. We use a stable recurrence relation:
    // (l+1) P_{l+1}(x) = (2l+1) x P_l(x) - l P_{l-1}(x).
    // This uses initial conditions:
    // P_0(x) = 1
    // P_1(x) = x

    if (l == 0) {

        legPolyL = 1.0_rt;

    }
    else if (l == 1) {

        legPolyL1 = legPolyL;
        legPolyL  = x;

    }
    else {

        legPolyL2 = legPolyL1;
        legPolyL1 = legPolyL;
        legPolyL  = ((2*l - 1) * x * legPolyL1 - (l-1) * legPolyL2) / l;

    }
}

AMREX_GPU_DEVICE AMREX_INLINE
void multipole_add(Real cosTheta, Real phiAngle, Real r, Real rho, Real vol,
                   Array4<Real> const& qL0,
                   Array4<Real> const& qLC,
                   Array4<Real> const& qLS,
                   Array4<Real> const& qU0,
                   Array4<Real> const& qUC,
                   Array4<Real> const& qUS,
                   int npts, int nlo, int index,
                   amrex::Gpu::Handler const& handler,
                   bool parity = false)
{
    Real legPolyL, legPolyL1, legPolyL2;
    Real assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2;

    Real rho_r_L, rho_r_U;

    for (int n = nlo; n <= npts-1; ++n) {

        for (int l = 0; l <= gravity::lnum; ++l) {

            calcLegPolyL(l, legPolyL, legPolyL1, legPolyL2, cosTheta);

            Real dQL0 = 0.0_rt;
            Real dQU0 = 0.0_rt;

            if (index <= n) {

                rho_r_L = rho * (std::pow(r, l));

                dQL0 = legPolyL * rho_r_L * vol * multipole::volumeFactor;
                if (parity) {
                    dQL0 = dQL0 * multipole::parity_q0(l);
                }

            }
            else {

                rho_r_U = rho * (std::pow(r, -l-1));

                dQU0 = legPolyL * rho_r_U * vol * multipole::volumeFactor;
                if (parity) {
                    dQU0 = dQU0 * multipole::parity_q0(l);
                }

            }

            amrex::Gpu::deviceReduceSum(&qL0(l,0,n), dQL0, handler);
            amrex::Gpu::deviceReduceSum(&qU0(l,0,n), dQU0, handler);

        }

        // For the associated Legendre polynomial loop, we loop over m and then l.
        // It means that we have to recompute rho_r_L or rho_r_U again, but the
        // recursion relation we use for the polynomials depends on l being the
        // innermost loop index.

        for (int m = 1; m <= gravity::lnum; ++m) {
            for (int l = 1; l <= gravity::lnum; ++l) {

                if (m > l) {
                    continue;
                }

                calcAssocLegPolyLM(l, m, assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2, cosTheta);

                Real dQLC = 0.0_rt;
                Real dQLS = 0.0_rt;

                Real dQUC = 0.0_rt;
                Real dQUS = 0.0_rt;

                if (index <= n) {

                    rho_r_L = rho * (std::pow(r, l));

                    dQLC = assocLegPolyLM * std::cos(m * phiAngle) * rho_r_L * vol * multipole::factArray(l,m);
                    if (parity) {
                        dQLC = dQLC * multipole::parity_qC_qS(l,m);
                    }

                    dQLS = assocLegPolyLM * std::sin(m * phiAngle) * rho_r_L * vol * multipole::factArray(l,m);
                    if (parity) {
                        dQLS = dQLS * multipole::parity_qC_qS(l,m);
                    }

                }
                else {

                    rho_r_U = rho * (std::pow(r, -l-1));

                    dQUC = assocLegPolyLM * std::cos(m * phiAngle) * rho_r_U * vol * multipole::factArray(l,m);
                    if (parity) {
                        dQUC = dQUC * multipole::parity_qC_qS(l,m);
                    }

                    dQUS = assocLegPolyLM * std::sin(m * phiAngle) * rho_r_U * vol * multipole::factArray(l,m);
                    if (parity) {
                        dQUS = dQUS * multipole::parity_qC_qS(l,m);
                    }

                }

                amrex::Gpu::deviceReduceSum(&qLC(l,m,n), dQLC, handler);
                amrex::Gpu::deviceReduceSum(&qLS(l,m,n), dQLS, handler);
                amrex::Gpu::deviceReduceSum(&qUC(l,m,n), dQUC, handler);
                amrex::Gpu::deviceReduceSum(&qUS(l,m,n), dQUS, handler);

            }
        }

    }
}

AMREX_GPU_DEVICE AMREX_INLINE
void multipole_symmetric_add(Real x, Real y, Real z,
                             const GpuArray<Real, AMREX_SPACEDIM>& problo,
                             const GpuArray<Real, AMREX_SPACEDIM>& probhi,
                             Real rho, Real vol,
                             Array4<Real> const& qL0,
                             Array4<Real> const& qLC,
                             Array4<Real> const& qLS,
                             Array4<Real> const& qU0,
                             Array4<Real> const& qUC,
                             Array4<Real> const& qUS,
                             int npts, int nlo, int index,
                             amrex::Gpu::Handler const& handler)
{

    amrex::ignore_unused(probhi);

    Real xLo = (2.0_rt * (problo[0] - problem::center[0])) / multipole::rmax - x;

#if AMREX_SPACEDIM >= 2
    Real yLo = (2.0_rt * (problo[1] - problem::center[1])) / multipole::rmax - y;
#else
    Real yLo = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    Real zLo = (2.0_rt * (problo[2] - problem::center[2]) ) / multipole::rmax - z;
#else
    Real zLo = 0.0_rt;
#endif

    Real r, phiAngle, cosTheta;

    if (multipole::doSymmetricAddLo(0)) {

        r        = std::sqrt(xLo * xLo + y * y + z * z);
        phiAngle = std::atan2(y, xLo);
        cosTheta = z / r;

        multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

        if (multipole::doSymmetricAddLo(1)) {

            r        = std::sqrt(xLo * xLo + yLo * yLo + z * z);
            phiAngle = std::atan2(yLo, xLo);
            cosTheta = z / r;

            multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

        }

        if (multipole::doSymmetricAddLo(2)) {

            r        = std::sqrt(xLo * xLo + y * y + zLo * zLo);
            phiAngle = std::atan2(y, xLo);
            cosTheta = zLo / r;

            multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

        }

        if (multipole::doSymmetricAddLo(1) && multipole::doSymmetricAddLo(2)) {

            r        = std::sqrt(xLo * xLo + yLo * yLo + zLo * zLo);
            phiAngle = std::atan2(yLo, xLo);
            cosTheta = zLo / r;

            multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

        }

    }

    if (multipole::doSymmetricAddLo(1)) {

        r        = std::sqrt(x * x + yLo * yLo + z * z);
        phiAngle = std::atan2(yLo, x);
        cosTheta = z / r;

        multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

        if (multipole::doSymmetricAddLo(2)) {

            r        = std::sqrt(x * x + yLo * yLo + zLo * zLo);
            phiAngle = std::atan2(yLo, x);
            cosTheta = zLo / r;

            multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

        }

    }

    if (multipole::doSymmetricAddLo(2)) {

        r        = std::sqrt(x * x + y * y + zLo * zLo);
        phiAngle = std::atan2(y, x);
        cosTheta = zLo / r;

        multipole_add(cosTheta, phiAngle, r, rho, vol, qL0, qLC, qLS, qU0, qUC, qUS, npts, nlo, index, handler);

    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real direct_sum_symmetric_add(const GpuArray<Real, 3>& loc, const GpuArray<Real, 3>& locb,
                              const GpuArray<Real, 3>& problo, const GpuArray<Real, 3>& probhi,
                              Real rho, Real dV,
                              const GpuArray<bool, 3>& doSymmetricAddLo, const GpuArray<bool, 3>& doSymmetricAddHi)
{
    // Add contributions from any symmetric boundaries.

    Real bcTerm = 0.0_rt;

    Real x, y, z, r;

    if (doSymmetricAddLo[0]) {

        x = 2.0_rt * problo[0] - loc[0];
        y = loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddLo[1]) {

            x = 2.0_rt * problo[0] - loc[0];
            y = 2.0_rt * problo[1] - loc[1];
            z = loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

        if (doSymmetricAddLo[2]) {

            x = 2.0_rt * problo[0] - loc[0];
            y = loc[1];
            z = 2.0_rt * problo[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

        if (doSymmetricAddLo[1] && doSymmetricAddLo[2]) {

            x = 2.0_rt * problo[0] - loc[0];
            y = 2.0_rt * problo[1] - loc[1];
            z = 2.0_rt * problo[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

    }

    if (doSymmetricAddLo[1]) {

        x = loc[0];
        y = 2.0_rt * problo[1] - loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddLo[2]) {

            x = loc[0];
            y = 2.0_rt * problo[1] - loc[1];
            z = 2.0_rt * problo[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

    }

    if (doSymmetricAddLo[2]) {

        x = loc[0];
        y = loc[1];
        z = 2.0_rt * problo[2] - loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

    }



    if (doSymmetricAddHi[0]) {

        x = 2.0_rt * probhi[0] - loc[0];
        y = loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddHi[1]) {

            x = 2.0_rt * probhi[0] - loc[0];
            y = 2.0_rt * probhi[1] - loc[1];
            z = loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

       if (doSymmetricAddHi[2]) {

           x = 2.0_rt * probhi[0] - loc[0];
           y = loc[1];
           z = 2.0_rt * probhi[2] - loc[2];

           r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

           bcTerm -= C::Gconst * rho * dV / r;

       }

       if (doSymmetricAddHi[1] && doSymmetricAddHi[2]) {

           x = 2.0_rt * probhi[0] - loc[0];
           y = 2.0_rt * probhi[1] - loc[1];
           z = 2.0_rt * probhi[2] - loc[2];

           r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

           bcTerm -= C::Gconst * rho * dV / r;

       }

    }

    if (doSymmetricAddHi[1]) {

        x = loc[0];
        y = 2.0_rt * probhi[1] - loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddHi[2]) {

            x = loc[0];
            y = 2.0_rt * probhi[1] - loc[1];
            z = 2.0_rt * probhi[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

    }

    if (doSymmetricAddHi[2]) {

        x = loc[0];
        y = loc[1];
        z = 2.0_rt * probhi[2] - loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

    }

    return bcTerm;

}

#endif
