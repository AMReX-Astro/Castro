#ifndef blackbody_H
#define blackbody_H

#include <fundamental_constants.H>

namespace blackbody {
    const Real pi = 3.1415926535897932384626e0_rt;
    const Real bk_const = C::a_rad * 15.e0_rt / (pi * pi * pi * pi);
    const Real magic_const = (pi * pi * pi * pi) / 15.e0_rt;
    const Real tol = 1.0e-10_rt;

    const Real xmagic = 2.061981e0_rt;
    const Real xsmall = 1.e-5_rt;
    const Real xlarge = 100.e0_rt;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real Li(int n, Real z)
{
    const int kmax = 18;

    // Clark, Equation 19 / ALEGRA, Equation 2.1.64
    //
    // z here is really exp(-h nu / kT), and nu and T are positive,
    // so it satisfies the range of convergence.

    Real L = z;

    for (int k = 2; k <= kmax; ++k) {

        Real t = std::pow(z, k) / std::pow(k, n);

        L = L + t;

        // Terminate the series when the additional terms
        // become small enough that they approach roundoff error.

        if (t / L < blackbody::tol) {
            break;
        }

    }

    return L;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real integlarge(Real x)
{
    // ALEGRA, Equation 2.1.63
    //
    // Note that since we define z == exp(-x), then
    // x = -ln(z), so the signs below are consistent
    // with that relative to the paper (i.e. the terms
    // with odd powers of x have reversed sign).
    //
    // We are only evaluating this at a specific frequency,
    // and relying on the recursion relation to give the
    // correct absolute result.

    Real z = std::exp(-x);

    Real I = blackbody::magic_const -
             (x * x * x * Li(1,z) + 3.e0_rt * x * x * Li(2,z) +
              6.e0_rt * x * Li(3,z) + 6.e0_rt * Li(4,z));

    return I;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real integsmall(Real x)
{
    // ALEGRA, Equation 2.1.65

    Real x2 = x * x;
    Real x3 = x * x * x;

    Real I = x3 / 3.0_rt - x2 * x2 / 8.0_rt;

    Real xfoo = x3;

    const int Ncoefs = 8;
    const Real coefs[Ncoefs+1] = {1.e0_rt / 60.e0_rt,
                                  -1.e0_rt / 5040.e0_rt,
                                  1.e0_rt / 272160.e0_rt,
                                  -1.e0_rt / 13305600.e0_rt,
                                  1.e0_rt / 622702080.e0_rt,
                                  -691.e0_rt / 19615115520000.e0_rt,
                                  1.e0_rt / 1270312243200.e0_rt,
                                  -3617.e0_rt / 202741834014720000.e0_rt,
                                  43867.e0_rt / 107290978560589824000.e0_rt};

    for (int n = 0; n <= Ncoefs; ++n) {

        xfoo = xfoo * x2;
        Real t = coefs[n] * xfoo;
        I = I + t;

        // Terminate the series expansion early if the terms
        // approach numerical roundoff error.

        if (std::abs(t / I) < blackbody::tol) {
            break;
        }

    }

    return I;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void BdBdTIndefInteg (Real T, Real nu, Real& B, Real& dBdT)
{
    // This routine evaluates the incomplete integral of the Planck function:
    // integ(x) = int_{0}^{x} b(x') dx'
    //
    // We follow B. A. Clark, "Computing Multigroup Radiation Integrals using
    // Polylogarithm-Based Methods," JCP (1987), 70, 311.
    //
    // This incomplete integral can be used to recursively obtain the blackbody
    // function for all groups using Equation 35:
    // b_g = integ(x_{g+1}) - integ(x_{g})
    //
    // The incomplete integral is defined by Equation 25 in terms of polylogarithm
    // functions evaluated at exp(-x) where x = h * nu / k_B T. As Clark notes,
    // we evaluate those polylogarithm functions directly when x is sufficiently
    // large, approximating them with a truncated infinite series; but, when x is
    // sufficiently small, we instead evaluate a series expansion of the original
    // integral. The threshold between these two approaches is set so that there
    // is no discontinuity of the Planck function at the threshold.

    // Specifically, we follow Sandia Report SAND2005-6988:
    // Advances in Radiation Modeling in ALEGRA: A Final Report for LDRD-67120,
    // Efficient Implicit Multigroup Radiation Calculations
    // T. A. Brunner, T. Mehlhorn, R. McClarren, & C. J. Kurecka (2005)
    //
    // They set the threshold at the value 2.061981, above which we evaluate the
    // polylogarithm functions and below which we use the series expansion.
    //
    // We can also use this method to evaluate the derivative of the Planck function
    // with respect to temperature.

    Real x = C::hplanck * nu / (C::k_B * T);

    if (x > blackbody::xlarge) {

        // If x is sufficiently large, this is effectively an integral over
        // all frequencies, and we know analytically that this is Stefan-Boltzmann
        // radiation, so we can save time (and avoid difficulties with evaluating
        // the exponential numerically for large argument) by evaluating it directly.

        B = C::a_rad * std::pow(T, 4);
        dBdT = 4.0_rt * C::a_rad * std::pow(T, 3);

    }
    else if (x < blackbody::xsmall) {

        // When x is very small, the integral is effectively over zero frequency
        // space and we can approximate it to zero.

        B = 0.0_rt;
        dBdT = 0.0_rt;

    }
    else {

        Real integ;
        if (x > blackbody::xmagic) {
            integ = integlarge(x);
        }
        else {
            integ = integsmall(x);
        }

        // Clark, Equation 3

        B = blackbody::bk_const * std::pow(T, 4) * integ;

        // ALEGRA, Equation 2.1.66
        // Note that the second term in square brackets is the one
        // we care about for this group, and the third term is from
        // the recursive relationship between groups. Since we are
        // integrating upwards from zero frequency, the sign is reversed.

        Real part = std::pow(x, 4) / (std::exp(x) - 1.0_rt);
        dBdT = blackbody::bk_const * std::pow(T, 3) * (4.0_rt * integ - part);

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real BIndefInteg(Real T, Real nu)
{
    Real B;

    // See comments above; this function evaluates B but does not
    // evaluate dB/dT.

    Real x = C::hplanck * nu / (C::k_B * T);

    if (x > blackbody::xlarge) {

        B = C::a_rad * std::pow(T, 4);

    }
    else if (x < blackbody::xsmall) {

        B = 0.e0_rt;

    }
    else {

        Real integ;
        if (x > blackbody::xmagic) {
            integ = integlarge(x);
        }
        else {
            integ = integsmall(x);
        }

        B = blackbody::bk_const * std::pow(T, 4) * integ;

    }

    return B;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real BGroup(Real T, Real nu0, Real nu1)
{
    // Clark, Equation 35

    Real B = BIndefInteg(T, nu1) - BIndefInteg(T, nu0);

    return B;
}

#endif
