#include <wdmerger_util.H>

extern "C" {

// Given total mass of a binary system and the initial separation of
// two point particles, obtain the velocity at this separation
// assuming the point masses fell in from infinity. This will
// be the velocity in the frame where the center of mass is stationary.

void freefall_velocity (Real mass, Real distance, Real& vel)
{
    vel = std::sqrt(2.0_rt * C::Gconst * mass / distance);
}

// Accepts the masses of two stars (in solar masses)
// and the orbital period of a system,
// and returns the semimajor axis of the orbit (in cm),
// as well as the distances a_1 and a_2 from the center of mass.

void kepler_third_law (Real radius_1, Real mass_1, Real radius_2, Real mass_2,
                       Real& period, Real eccentricity, Real phi, Real& a,
                       Real& r_1, Real& r_2, Real& v_1r, Real& v_2r, Real& v_1p, Real& v_2p)
{
    Real length;

    Real M  = mass_1 + mass_2;     // Total mass
    Real mu = mass_1 * mass_2 / M; // Reduced mass

    // First, solve for the orbit in the reduced one-body problem, where
    // an object of mass mu orbits an object with mass M located at r = 0.
    // For this we follow Carroll and Ostlie, Chapter 2, but many texts discuss this.
    // Note that we use the convention that phi measures angle from aphelion,
    // which is opposite to the convention they use.

    if (period > 0.0_rt && a < 0.0_rt)
    {
        a = std::pow(C::Gconst * M * period * period / (4.0_rt * M_PI * M_PI), 1.0_rt / 3.0_rt); // C + O, Equation 2.37
    }
    else if (period < 0.0_rt && a > 0.0_rt)
    {
        period = std::sqrt(a * a * a * 4.0_rt * M_PI * M_PI / (C::Gconst * M));
    }
    else {
        amrex::Error("Error: overspecified Kepler's third law calculation.");
    }

    Real r = a * (1.0_rt - eccentricity * eccentricity) / (1.0_rt - eccentricity * std::cos(phi)); // C + O, Equation 2.3

    // To get the radial and azimuthal velocity, we take the appropriate derivatives of the above.
    // v_r = dr / dt = dr / d(phi) * d(phi) / dt, with d(phi) / dt being derived from
    // C + O, Equation 2.30 for the angular momentum, and the fact that L = mu * r**2 * d(phi) / dt.

    Real v_r   = -2.0_rt * M_PI * a * eccentricity * std::sin(phi) /
                  (period * std::sqrt(1.0_rt - eccentricity * eccentricity));
    Real v_phi =  2.0_rt * M_PI * a * (1.0_rt - eccentricity * std::cos(phi)) /
                  (period * std::sqrt(1.0_rt - eccentricity * eccentricity));

    // Now convert everything back to the binary frame, using C+O, Equation 2.23 and 2.24. This applies
    // to the velocities as well as the positions because the factor in front of r_1 and r_2 is constant.

    r_1  = -(mu / mass_1) * r;
    r_2  =  (mu / mass_2) * r;

    v_1r = -(mu / mass_1) * v_r;
    v_2r =  (mu / mass_2) * v_r;

    v_1p = -(mu / mass_1) * v_phi;
    v_2p =  (mu / mass_2) * v_phi;
}

}
