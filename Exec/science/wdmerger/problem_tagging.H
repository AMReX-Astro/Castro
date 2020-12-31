#ifndef problem_tagging_H
#define problem_tagging_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();

    Real loc[3] = {0.0};

    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
#if AMREX_SPACEDIM >= 2
    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#endif
#if AMREX_SPACEDIM == 3
    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#endif

    if (level < problem::max_stellar_tagging_level) {

        if (problem::problem == 0 || problem::problem == 2) {

            // For the collision, free-fall, and TDE problems, we just want to tag every
            // zone that meets the density criterion; we don't want to bother with
            // the Roche lobe radius as that doesn't mean much in these cases.

            if (state(i,j,k,URHO) > problem::stellar_density_threshold) {

                tag(i,j,k) = TagBox::SET;

            }

        }
        else {

            if (level == 0) {

                // On the coarse grid, tag all regions within the Roche radii of each star.
                // We'll add a buffer around each star to double the Roche
                // radius to ensure there aren't any sharp gradients in regions of
                // greater than ambient density.

                Real r_P = std::sqrt((loc[0] - problem::com_P[0]) * (loc[0] - problem::com_P[0]) +
                                     (loc[1] - problem::com_P[1]) * (loc[1] - problem::com_P[1]) +
                                     (loc[2] - problem::com_P[2]) * (loc[2] - problem::com_P[2]));

                Real r_S = std::sqrt((loc[0] - problem::com_S[0]) * (loc[0] - problem::com_S[0]) +
                                     (loc[1] - problem::com_S[1]) * (loc[1] - problem::com_S[1]) +
                                     (loc[2] - problem::com_S[2]) * (loc[2] - problem::com_S[2]));

                if (r_P <= problem::roche_tagging_factor * problem::roche_rad_P) {
                    tag(i,j,k) = TagBox::SET;
                }

                if (r_S <= problem::roche_tagging_factor * problem::roche_rad_S) {
                    tag(i,j,k) = TagBox::SET;
                }

            }
            else if (level >= 1) {

                // On more refined levels, tag all regions within the stars themselves (defined as
                // areas where the density is greater than some threshold).

                if (state(i,j,k,URHO) > problem::stellar_density_threshold) {

                    tag(i,j,k) = TagBox::SET;

                }

            }

        }

        // Tag all zones at all levels that are hotter than a specified temperature threshold.

        if (level < problem::max_temperature_tagging_level) {

            if (state(i,j,k,UTEMP) > problem::temperature_tagging_threshold) {

                tag(i,j,k) = TagBox::SET;

            }

        }

        // Clear all tagging that occurs outside the radius set by max_tagging_radius.

        Real r = std::sqrt((loc[0] - problem::center[0]) * (loc[0] - problem::center[0]) +
                           (loc[1] - problem::center[1]) * (loc[1] - problem::center[1]) +
                           (loc[2] - problem::center[2]) * (loc[2] - problem::center[2]));

        Real max_dist_lo = 0.0;
        Real max_dist_hi = 0.0;

        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
            max_dist_lo = amrex::max(max_dist_lo, std::abs(problo[dim] - problem::center[dim]));
            max_dist_hi = amrex::max(max_dist_hi, std::abs(probhi[dim] - problem::center[dim]));
        }

        if (r > problem::max_tagging_radius * amrex::max(max_dist_lo, max_dist_hi)) {

            tag(i,j,k) = TagBox::CLEAR;

        }

    }
}

#endif
