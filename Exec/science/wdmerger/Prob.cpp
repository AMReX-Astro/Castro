#include <Castro.H>
#include <Castro_util.H>

#include <Gravity.H>

#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>

#include <wdmerger_util.H>
#include <wdmerger_data.H>
#include <binary.H>

#include <fstream>

using namespace amrex;

void
Castro::problem_post_timestep()
{
    BL_PROFILE("Castro::problem_post_timestep()");

    using namespace wdmerger;
    using namespace problem;

    if (level != 0) return;

    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite

    // Update white dwarf masses, positions, velocities, and auxiliary data.

    wd_update(time, dt);

    // If we are doing the merger problem with an initial relaxation step,
    // perform any post-timestep updates to assist with the relaxation, then
    // determine whether the criterion for terminating the relaxation
    // has been satisfied.

    update_relaxation(time, dt);

    // Compute any integral quantities.

    problem_sums();
}



//
// This function updates the WD data, including the masses, center-of-mass locations,
// and velocities of the center of masses of the primary and secondary white dwarfs.
//

void
Castro::wd_update (Real time, Real dt)
{
    BL_PROFILE("Castro::wd_update()");

    using namespace wdmerger;
    using namespace problem;

    // Ensure we are either on the coarse level, or on the finest level
    // when we are not doing subcycling. The data should be synchronized
    // in both of these cases.

    BL_ASSERT(level == 0 || (!parent->subCycle() && level == parent->finestLevel()));

    for (int i = 0; i < 3; i++) {
      com_P[i] += vel_P[i] * dt;
      com_S[i] += vel_S[i] * dt;
    }

    // Now re-calculate a more accurate result using this as a starting point.

    // Save relevant current data.

    Real old_mass_P = mass_P;
    Real old_mass_S = mass_S;

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,
              ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,
              ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Real, Real, Real, Real,
               Real, Real, Real, Real, Real, Real,
               Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (int lev = 0; lev <= parent->finestLevel(); lev++) {

        Castro& c_lev = getLevel(lev);
        MultiFab& S_new = c_lev.get_new_data(State_Type);

        GeometryData geomdata = c_lev.geom.data();

        const auto problo = c_lev.geom.ProbLoArray();
        const auto probhi = c_lev.geom.ProbHiArray();
        int coord_type = c_lev.geom.Coord();

        GpuArray<bool, 3> symm_bound_lo{false};
        GpuArray<bool, 3> symm_bound_hi{false};

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            if (phys_bc.lo()[n] == amrex::PhysBCType::symmetry) {
                symm_bound_lo[n] = true;
            }
            if (phys_bc.hi()[n] == amrex::PhysBCType::symmetry) {
                symm_bound_hi[n] = true;
            }
        }

        bool mask_available = true;
        if (lev == parent->finestLevel()) {
            mask_available = false;
        }

        MultiFab tmp_mf;
        const MultiFab& mask_mf = mask_available ? getLevel(lev+1).build_fine_mask() : tmp_mf;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto rho   = S_new[mfi].array(URHO);
            auto xmom  = S_new[mfi].array(UMX);
            auto ymom  = S_new[mfi].array(UMY);
            auto zmom  = S_new[mfi].array(UMZ);
            auto vol   = c_lev.volume[mfi].array();
            auto level_mask = mask_available ? mask_mf[mfi].array() : Array4<Real>{};

            const Box& box  = mfi.tilebox();

            // Return the mass-weighted center of mass and velocity
            // for the primary and secondary, for a given FAB.
            // Note that ultimately what we are doing here is to use
            // an old guess at the effective potential of the primary
            // and secondary to generate a new estimate.

            reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                // Add to the COM locations and velocities of the primary and secondary
                // depending on which potential dominates, ignoring unbound material.
                // Note that in this routine we actually are summing mass-weighted
                // quantities for the COM and the velocity; we will account for this at
                // the end of the calculation in post_timestep() by dividing by the mass.

                // Our convention is that the COM locations for the WDs are
                // absolute positions on the grid, not relative to the center.

                GpuArray<Real, 3> r;
                position(i, j, k, geomdata, r);

                // We account for symmetric boundaries in this sum as usual,
                // by adding to the position the locations that would exist
                // on the opposite side of the symmetric boundary. Note that
                // in axisymmetric coordinates, some of this work is already
                // done for us in the definition of the zone volume.

                GpuArray<Real, 3> rSymmetric{r[0], r[1], r[2]};
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    if (symm_bound_lo[n]) {
                        rSymmetric[n] = rSymmetric[n] + (problo[n] - rSymmetric[n]);
                    }
                    if (symm_bound_hi[n]) {
                        rSymmetric[n] = rSymmetric[n] + (rSymmetric[n] - probhi[n]);
                    }
                }

                Real maskFactor = 1.0;

                if (mask_available) {
                    maskFactor = level_mask(i,j,k);
                }

                Real dm = rho(i,j,k) * vol(i,j,k) * maskFactor;

                Real dmSymmetric = dm;
                GpuArray<Real, 3> momSymmetric{xmom(i,j,k) * maskFactor,
                                               ymom(i,j,k) * maskFactor,
                                               zmom(i,j,k) * maskFactor};

                if (coord_type == 0) {

                    if (symm_bound_lo[0]) {
                        dmSymmetric *= 2.0_rt;
                        for (int n = 0; n < 3; ++n) {
                            momSymmetric[n] *= 2.0_rt;
                        }
                    }

                    if (symm_bound_lo[1]) {
                        dmSymmetric *= 2.0_rt;
                        for (int n = 0; n < 3; ++n) {
                            momSymmetric[n] *= 2.0_rt;
                        }
                    }

                    if (symm_bound_lo[2]) {
                        dmSymmetric *= 2.0_rt;
                        for (int n = 0; n < 3; ++n) {
                            momSymmetric[n] *= 2.0_rt;
                        }
                    }

                }

                Real primary_factor = 0.0_rt;
                Real secondary_factor = 0.0_rt;

                if (stellar_mask(i, j, k, geomdata, rho(i,j,k), true) * maskFactor > 0.0_rt) {

                    primary_factor = 1.0_rt;

                }
                else if (stellar_mask(i, j, k, geomdata, rho(i,j,k), false) * maskFactor > 0.0_rt) {

                    secondary_factor = 1.0_rt;

                }

                Real com_P_x = dmSymmetric * rSymmetric[0] * primary_factor;
                Real com_P_y = dmSymmetric * rSymmetric[1] * primary_factor;
                Real com_P_z = dmSymmetric * rSymmetric[2] * primary_factor;

                Real com_S_x = dmSymmetric * rSymmetric[0] * secondary_factor;
                Real com_S_y = dmSymmetric * rSymmetric[1] * secondary_factor;
                Real com_S_z = dmSymmetric * rSymmetric[2] * secondary_factor;

                Real vel_P_x = momSymmetric[0] * vol(i,j,k) * primary_factor;
                Real vel_P_y = momSymmetric[1] * vol(i,j,k) * primary_factor;
                Real vel_P_z = momSymmetric[2] * vol(i,j,k) * primary_factor;

                Real vel_S_x = momSymmetric[0] * vol(i,j,k) * secondary_factor;
                Real vel_S_y = momSymmetric[1] * vol(i,j,k) * secondary_factor;
                Real vel_S_z = momSymmetric[2] * vol(i,j,k) * secondary_factor;

                Real m_P = dmSymmetric * primary_factor;
                Real m_S = dmSymmetric * secondary_factor;

                return {com_P_x, com_P_y, com_P_z, com_S_x, com_S_y, com_S_z,
                        vel_P_x, vel_P_y, vel_P_z, vel_S_x, vel_S_y, vel_S_z,
                        m_P, m_S};
            });

        }

    }

    // Compute effective radii of stars at various density cutoffs

    bool local_flag = true;

    for (int i = 0; i <= 6; ++i) {
        Castro::volInBoundary(time, vol_P[i], vol_S[i], std::pow(10.0,i), local_flag);
    }

    // Do all of the reductions.

    ReduceTuple hv = reduce_data.value();

    com_P[0] = amrex::get<0>(hv);
    com_P[1] = amrex::get<1>(hv);
    com_P[2] = amrex::get<2>(hv);
    com_S[0] = amrex::get<3>(hv);
    com_S[1] = amrex::get<4>(hv);
    com_S[2] = amrex::get<5>(hv);
    vel_P[0] = amrex::get<6>(hv);
    vel_P[1] = amrex::get<7>(hv);
    vel_P[2] = amrex::get<8>(hv);
    vel_S[0] = amrex::get<9>(hv);
    vel_S[1] = amrex::get<10>(hv);
    vel_S[2] = amrex::get<11>(hv);
    mass_P   = amrex::get<12>(hv);
    mass_S   = amrex::get<13>(hv);

    const int nfoo_sum = 28;
    Real foo_sum[nfoo_sum] = { 0.0 };

    for (int i = 0; i <= 6; ++i) {
        foo_sum[i  ] = vol_P[i];
        foo_sum[i+7] = vol_S[i];
    }

    foo_sum[14] = mass_P;
    foo_sum[15] = mass_S;

    for (int i = 0; i <= 2; ++i) {
        foo_sum[i+16] = com_P[i];
        foo_sum[i+19] = com_S[i];
        foo_sum[i+22] = vel_P[i];
        foo_sum[i+25] = vel_S[i];
    }

    amrex::ParallelDescriptor::ReduceRealSum(foo_sum, nfoo_sum);

    for (int i = 0; i <= 6; ++i) {
        vol_P[i] = foo_sum[i  ];
        vol_S[i] = foo_sum[i+7];
    }

    mass_P = foo_sum[14];
    mass_S = foo_sum[15];

    for (int i = 0; i <= 2; ++i) {
        com_P[i] = foo_sum[i+16];
        com_S[i] = foo_sum[i+19];
        vel_P[i] = foo_sum[i+22];
        vel_S[i] = foo_sum[i+25];
    }

    // Compute effective WD radii

    for (int i = 0; i <= 6; ++i) {
        rad_P[i] = std::pow(vol_P[i] * 3.0 / 4.0 / M_PI, 1.0/3.0);
        rad_S[i] = std::pow(vol_S[i] * 3.0 / 4.0 / M_PI, 1.0/3.0);
    }

     // Complete calculations for center of mass quantities

    for (int i = 0; i < 3; i++) {
        if (mass_P > 0.0) {
            com_P[i] = com_P[i] / mass_P;
            vel_P[i] = vel_P[i] / mass_P;
        }

        if (mass_S > 0.0) {
            com_S[i] = com_S[i] / mass_S;
            vel_S[i] = vel_S[i] / mass_S;
        }
    }

    // For 1D we force the masses to remain constant

#if (AMREX_SPACEDIM == 1)
    mass_P = old_mass_P;
    mass_S = old_mass_S;
#endif

    if (mass_P > 0.0 && dt > 0.0) {
        mdot_P = (mass_P - old_mass_P) / dt;
    }
    else {
        mdot_P = 0.0;
    }

    if (mass_S > 0.0 && dt > 0.0) {
        mdot_S = (mass_S - old_mass_S) / dt;
    }
    else {
        mdot_S = 0.0;
    }

    // Free-fall timescale ~ 1 / sqrt(G * rho_avg}

    if (mass_P > 0.0 && vol_P[2] > 0.0) {
        rho_avg_P = mass_P / vol_P[2];
        t_ff_P = std::sqrt(3.0 * M_PI / (32.0 * C::Gconst * rho_avg_P));
    }

    if (mass_S > 0.0 && vol_S[2] > 0.0) {
        rho_avg_S = mass_S / vol_S[2];
        t_ff_S = std::sqrt(3.0 * M_PI / (32.0 * C::Gconst * rho_avg_S));
    }

    // Compute updated roche Radii

    update_roche_radii();
}



// This function uses the known center of mass of the two white dwarfs,
// and given a density cutoff, computes the total volume of all zones
// whose density is greater or equal to that density cutoff.
// We also impose a distance requirement so that we only look
// at zones that are within twice the original radius of the white dwarf.

void Castro::volInBoundary (Real time, Real& vol_P, Real& vol_S, Real rho_cutoff, bool local)
{

    amrex::ignore_unused(time);

    BL_PROFILE("Castro::volInBoundary()");

    using namespace wdmerger;
    using namespace problem;

    BL_ASSERT(level == 0);

    vol_P = 0.0;
    vol_S = 0.0;

    for (int lev = 0; lev <= parent->finestLevel(); lev++) {

      Castro& c_lev = getLevel(lev);

      auto geomdata = c_lev.geom.data();
      MultiFab& S_new = c_lev.get_new_data(State_Type);

      bool mask_available = true;
      if (lev == parent->finestLevel()) {
          mask_available = false;
      }

      MultiFab tmp_mf;
      const MultiFab& mask_mf = mask_available ? getLevel(lev+1).build_fine_mask() : tmp_mf;

      ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
      ReduceData<Real, Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          auto rho   = S_new[mfi].array(URHO);
          auto vol   = c_lev.volume[mfi].array();
          auto level_mask = mask_available ? mask_mf[mfi].array() : Array4<Real>{};

          const Box& box  = mfi.tilebox();

          // This function uses the known center of mass of the two white dwarfs,
          // and given a density cutoff, computes the total volume of all zones
          // whose density is greater or equal to that density cutoff.
          // We also impose a distance requirement so that we only look
          // at zones within the Roche lobe of the white dwarf.

          reduce_op.eval(box, reduce_data,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
          {
              Real primary_factor = 0.0_rt;
              Real secondary_factor = 0.0_rt;

              Real maskFactor = 1.0;
              if (mask_available) {
                  maskFactor = level_mask(i,j,k);
              }

              if (rho(i,j,k) * maskFactor > rho_cutoff) {

                  if (stellar_mask(i, j, k, geomdata, rho(i,j,k), true) * maskFactor > 0.0_rt) {

                      primary_factor = 1.0_rt;

                  }
                  else if (stellar_mask(i, j, k, geomdata, rho(i,j,k), false) * maskFactor > 0.0_rt) {

                      secondary_factor = 1.0_rt;

                  }

              }

              return {vol(i,j,k) * maskFactor * primary_factor,
                      vol(i,j,k) * maskFactor * secondary_factor};
          });

      }

      ReduceTuple hv = reduce_data.value();

      vol_P += amrex::get<0>(hv);
      vol_S += amrex::get<1>(hv);

    }

    if (!local)
      amrex::ParallelDescriptor::ReduceRealSum({vol_P, vol_S});

}



// Computes standard dot-product of two three-vectors.

Real Castro::dot_product(const Real a[], const Real b[]) {

  Real c = 0.0;

  c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  return c;

}




// Computes norm of a three-vector.

Real Castro::norm(const Real a[]) {

  Real n = 0.0;

  n = std::sqrt(dot_product(a, a));

  return n;

}



void Castro::problem_post_init() {

  using namespace wdmerger;
  using namespace problem;

  // Execute the post timestep diagnostics here,
  // so that the results at t = 0 and later are smooth.
  // This should generally be the last operation
  // in this function.

  problem_post_timestep();

}



void Castro::problem_post_restart() {

  using namespace wdmerger;
  using namespace problem;

  // If we're restarting from a checkpoint at t = 0 but don't yet
  // have diagnostics, we want to generate the headers and the t = 0
  // data at this time so that future timestep diagnostics make sense.

  Real time = state[State_Type].curTime();

  if (time == 0.0) {

      if (parent->NumDataLogs() > 0) {

          bool do_sums = false;

          const std::string datalogname = parent->DataLogName(0);

          std::ifstream log;
          log.open(datalogname, std::ios::ate);

          if (log.tellg() == 0) {
              do_sums = true;
          }
          log.close();

          if (do_sums && (sum_interval > 0 || sum_per > 0))
              sum_integrated_quantities();

      }

  }

  // Update Roche radii to ensure consistency of initial conditions

  update_roche_radii();

}



void Castro::writeGitHashes(std::ostream& log) {

  const char* castro_hash       = buildInfoGetGitHash(1);
  const char* amrex_hash        = buildInfoGetGitHash(2);
  const char* microphysics_hash = buildInfoGetGitHash(3);

  log << "# Castro       git hash: " << castro_hash       << std::endl;
  log << "# AMReX        git hash: " << amrex_hash        << std::endl;
  log << "# Microphysics git hash: " << microphysics_hash << std::endl;

}



void
Castro::update_relaxation(Real time, Real dt) {

    using namespace wdmerger;
    using namespace problem;

    // Check to make sure whether we should be doing the relaxation here.
    // Update the relaxation conditions if we are not stopping.

    if (problem::problem != 1 || relaxation_is_done || mass_P <= 0.0 || mass_S <= 0.0 || dt <= 0.0) return;

    // Construct the update to the rotation frequency. We calculate
    // the gravitational force at the end of the timestep, set the
    // rotational force to be equal to it, and infer the required
    // rotation frequency. Then we apply it so that the next timestep
    // has an updated rotation force that is a better balance against
    // the gravitational force.

    int coarse_level = 0;
    int finest_level = parent->finestLevel();
    int n_levs = finest_level + 1;

    Vector< std::unique_ptr<MultiFab> > force(n_levs);

    for (int lev = coarse_level; lev <= finest_level; ++lev) {

        const Real old_time = getLevel(lev).state[State_Type].prevTime();
        const Real new_time = getLevel(lev).state[State_Type].curTime();

        const Real ldt = new_time - old_time;

        force[lev] = std::make_unique<MultiFab>(getLevel(lev).grids, getLevel(lev).dmap, NUM_STATE, 0);
        force[lev]->setVal(0.0);

        MultiFab& S_new = getLevel(lev).get_new_data(State_Type);

        // Store the non-rotation forces. The inspiration for this approach is
        // the method of Rosswog, Speith & Wynn (2004) (which was in the context
        // of neutron star mergers; it was extended to white dwarf mergers by Dan et al.
        // (2011)). In that paper, the rotation force is calculated by exactly balancing
        // against the gravitational and hydrodynamic forces. We will just use the
        // gravitational force, since the desired equilibrium state is for the hydrodynamic
        // forces to be zero.

        // We'll use the "old" gravity source constructor, which is really just a first-order
        // predictor for rho * g, and apply it at the new time.

        getLevel(lev).construct_old_gravity_source(*force[lev], S_new, new_time, ldt);

        // Mask out regions covered by fine grids.

        if (lev < parent->finestLevel()) {
            const MultiFab& mask = getLevel(lev+1).build_fine_mask();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(*force[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();

                auto F = (*force[lev])[mfi].array();
                auto m = mask[mfi].array();

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    for (int n = 0; n < NUM_STATE; ++n) {
                        F(i,j,k,n) *= m(i,j,k);
                    }
                });
            }
        }

    }

    // Construct omega from this.

    Real force_P[3] = { 0.0 };
    Real force_S[3] = { 0.0 };

    for (int lev = coarse_level; lev <= finest_level; ++lev) {

        MultiFab& S_new = getLevel(lev).get_new_data(State_Type);

        MultiFab& vol = getLevel(lev).Volume();

        auto geomdata = getLevel(lev).geom.data();

        const int coord_type = geom.Coord();

        const int* lo_bc = phys_bc.lo();

        const bool symm_lo_x = (lo_bc[0] == amrex::PhysBCType::symmetry);
        const bool symm_lo_y = (lo_bc[1] == amrex::PhysBCType::symmetry);
        const bool symm_lo_z = (lo_bc[2] == amrex::PhysBCType::symmetry);

        ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<Real, Real, Real, Real, Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            // Compute the sum of the hydrodynamic and gravitational forces acting on the WDs.

            Array4<Real const> const rho_arr   = S_new[mfi].array(URHO);
            Array4<Real const> const force_arr = (*force[lev]).array(mfi);
            Array4<Real const> const vol_arr   = vol.array(mfi);

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                GpuArray<Real, 3> dF;
                for (int n = 0; n < 3; ++n) {
                    dF[n] = vol_arr(i,j,k) * force_arr(i,j,k,UMX+n);
                }

                // In the following we'll account for symmetry boundaries
                // by assuming they're on the lower boundary if they do exist.
                // In that case, assume the star's center is on the lower boundary.
                // Then we need to double the value of the force in the other dimensions,
                // and cancel out the force in that dimension.
                // We assume here that only one dimension at most has a symmetry boundary.

                if (coord_type == 0) {

                    if (symm_lo_x) {
                        dF[0] = 0.0_rt;
                        dF[1] = 2.0_rt * dF[1];
                        dF[2] = 2.0_rt * dF[2];
                    }

                    if (symm_lo_y) {
                        dF[0] = 2.0_rt * dF[0];
                        dF[1] = 0.0_rt;
                        dF[2] = 2.0_rt * dF[2];
                    }

                    if (symm_lo_z) {
                        dF[0] = 2.0_rt * dF[0];
                        dF[1] = 2.0_rt * dF[1];
                        dF[2] = 0.0_rt;
                    }

                } else if (coord_type == 1) {

                    dF[0] = 0.0_rt;

                }

                Real primary_factor = 0.0_rt;
                Real secondary_factor = 0.0_rt;

                if (stellar_mask(i, j, k, geomdata, rho_arr(i,j,k), true) > 0.0_rt) {

                    primary_factor = 1.0_rt;

                } else if (stellar_mask(i, j, k, geomdata, rho_arr(i,j,k), false) > 0.0_rt) {

                    secondary_factor = 1.0_rt;

                }

                return {dF[0] * primary_factor,
                        dF[1] * primary_factor,
                        dF[2] * primary_factor,
                        dF[0] * secondary_factor,
                        dF[1] * secondary_factor,
                        dF[2] * secondary_factor};
            });

        }

        ReduceTuple hv = reduce_data.value();

        force_P[0] += amrex::get<0>(hv);
        force_P[1] += amrex::get<1>(hv);
        force_P[2] += amrex::get<2>(hv);

        force_S[0] += amrex::get<3>(hv);
        force_S[1] += amrex::get<4>(hv);
        force_S[2] += amrex::get<5>(hv);
    }

    // Do the reduction over processors.

    amrex::ParallelDescriptor::ReduceRealSum({force_P[0], force_P[1], force_P[2], force_S[0], force_S[1], force_S[2]});

    // Divide by the mass of the stars to obtain the acceleration, and then get the new rotation frequency.

    Real fp = std::sqrt(std::pow(force_P[0], 2) + std::pow(force_P[1], 2) + std::pow(force_P[2], 2));
    Real fs = std::sqrt(std::pow(force_S[0], 2) + std::pow(force_S[1], 2) + std::pow(force_S[2], 2));

    Real ap = std::sqrt(std::pow(com_P[0], 2) + std::pow(com_P[1], 2) + std::pow(com_P[2], 2));
    Real as = std::sqrt(std::pow(com_S[0], 2) + std::pow(com_S[1], 2) + std::pow(com_S[2], 2));

    Real omega = 0.5 * (std::sqrt((fp / mass_P) / ap) + std::sqrt((fs / mass_S) / as));

    Real period = 2.0 * M_PI / omega;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cout << "\n";
        std::cout << "  Updating the rotational period from " << rotational_period << " s to " << period << " s." << "\n";
        std::cout << "\n";
    }

    rotational_period = period;

    // Check to see whether the relaxation should be turned off.
    // Note that at present the following check is only done on the
    // coarse grid but if we wanted more accuracy we could do a loop
    // over levels as above.

    // For the merger problem, we're going to turn the relaxation off
    // when we've reached the L1 Lagrange point.

    // First, calculate the location of the L1 Lagrange point.

    GpuArray<Real, 3> L1, L2, L3;
    get_lagrange_points(mass_P, mass_S, com_P, com_S, L1, L2, L3);

    const auto dx = geom.CellSizeArray();
    GeometryData geomdata = geom.data();

    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    Real potential;

    {
        // Then, figure out the effective potential corresponding to that
        // Lagrange point.

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(phi_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            Array4<Real const> const phi = phi_new[mfi].array();

            // Determine the critical Roche potential at the Lagrange point L1.
            // We will use a tri-linear interpolation that gets a contribution
            // from all the zone centers that bracket the Lagrange point.

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                GpuArray<Real, 3> r;
                position(i, j, k, geomdata, r);

                // Compute the effective potential.

                auto omega = get_omega_vec(j);
                Real phiEff = phi(i,j,k) + rotational_potential(r, omega);

                for (int n = 0; n < 3; ++n) {
                    r[n] -= L1[n];
                }

                // Scale r by dx (in dimensions we're actually simulating).

                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    r[n] /= dx[n];
                }
                for (int n = AMREX_SPACEDIM; n < 3; ++n) {
                    r[n] = 0.0_rt;
                }

                // We want a contribution from this zone if it is
                // less than one zone width away from the Lagrange point.

                Real dP = 0.0_rt;

                if ((r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) < 1.0_rt) {
                    dP = (1.0_rt - std::abs(r[0])) * (1.0_rt - std::abs(r[1])) * (1.0_rt - std::abs(r[2])) * phiEff;
                }

                return dP;

            });

        }

        ReduceTuple hv = reduce_data.value();
        potential = amrex::get<0>(hv);

        amrex::ParallelDescriptor::ReduceRealSum(potential);
    }

    {
        // Now cycle through the grids and determine if any zones
        // have crossed the density threshold outside the critical surface.

        MultiFab& S_new = get_new_data(State_Type);

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            Array4<Real const> const u = S_new.array(mfi);
            Array4<Real const> const phi = phi_new.array(mfi);

            // Check whether we should stop the initial relaxation.
            // The criterion is that we're outside the critical Roche surface
            // and the density is greater than a specified threshold.
            // If so, set do_initial_relaxation to false, which will effectively
            // turn off the external source terms.

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                GpuArray<Real, 3> r;
                position(i, j, k, geomdata, r);

                auto omega = get_omega_vec(j);
                Real phiEff = phi(i,j,k) + rotational_potential(r, omega);

                Real done = 0.0_rt;

                if (phiEff > potential && u(i,j,k,URHO) > relaxation_density_cutoff) {
                    done = 1.0_rt;
                }

                return done;
            });

        }

        ReduceTuple hv = reduce_data.value();
        Real is_done = amrex::get<0>(hv);

        amrex::ParallelDescriptor::ReduceRealSum(is_done);

        if (is_done > 0.0) {
            relaxation_is_done = 1;
            amrex::Print() << "Disabling relaxation at time " << time
                           << "s because the critical density threshold has been passed."
                           << std::endl;
        }

    }

    // We can also turn off the relaxation if we've passed
    // a certain number of dynamical timescales.

    if (relaxation_cutoff_time > 0.0 && time > relaxation_cutoff_time * std::max(t_ff_P, t_ff_S)) {
        relaxation_is_done = 1;
        amrex::Print() << "Disabling relaxation at time " << time
                       << "s because the maximum number of dynamical timescales has passed."
                       << std::endl;
    }

    if (relaxation_is_done == 1 && radial_damping_velocity_factor > 0.0_rt && radial_damping_is_done != 1) {
        amrex::Print() << "Radial damping force will now take effect." << std::endl;
    }
}



void
Castro::problem_sums ()
{
    using namespace wdmerger;
    using namespace problem;

    if (level > 0) return;

    Real time = state[State_Type].curTime();

    int timestep = parent->levelSteps(0);

    // Mass transfer rate

    Real mdot = 0.5 * (std::abs(mdot_P) + std::abs(mdot_S));

    // Distance between the WDs.

    Real wd_dist[3] = { 0.0 };
    Real wd_dist_init[3] = { 0.0 };

    Real separation = 0.0;
    Real angle = 0.0;

    // Stellar centers of mass and velocities.

    Real com_P_mag = 0.0;
    Real com_S_mag = 0.0;

    Real vel_P_mag = 0.0;
    Real vel_S_mag = 0.0;

    Real vel_P_rad = 0.0;
    Real vel_S_rad = 0.0;

    Real vel_P_phi = 0.0;
    Real vel_S_phi = 0.0;

    int dataprecision = 16; // Number of digits after the decimal point, for float data

    int datwidth      = 25; // Floating point data in scientific notation
    int fixwidth      = 25; // Floating point data not in scientific notation
    int intwidth      = 12; // Integer data

    wd_dist_init[problem::axis_1 - 1] = 1.0;

    com_P_mag += std::pow(std::pow(com_P[0],2) + std::pow(com_P[1],2) + std::pow(com_P[2],2), 0.5);
    com_S_mag += std::pow(std::pow(com_S[0],2) + std::pow(com_S[1],2) + std::pow(com_S[2],2), 0.5);
    vel_P_mag += std::pow(std::pow(vel_P[0],2) + std::pow(vel_P[1],2) + std::pow(vel_P[2],2), 0.5);
    vel_S_mag += std::pow(std::pow(vel_S[0],2) + std::pow(vel_S[1],2) + std::pow(vel_S[2],2), 0.5);

#if (AMREX_SPACEDIM == 3)
    if (mass_P > 0.0) {
      vel_P_rad = (com_P[problem::axis_1 - 1] / com_P_mag) * vel_P[problem::axis_1 - 1] +
                  (com_P[problem::axis_2 - 1] / com_P_mag) * vel_P[problem::axis_2 - 1];
      vel_P_phi = (com_P[problem::axis_1 - 1] / com_P_mag) * vel_P[problem::axis_2 - 1] -
                  (com_P[problem::axis_2 - 1] / com_P_mag) * vel_P[problem::axis_1 - 1];
    }

    if (mass_S > 0.0) {
      vel_S_rad = (com_S[problem::axis_1 - 1] / com_S_mag) * vel_S[problem::axis_1 - 1] +
                  (com_S[problem::axis_2 - 1] / com_S_mag) * vel_S[problem::axis_2 - 1];
      vel_S_phi = (com_S[problem::axis_1 - 1] / com_S_mag) * vel_S[problem::axis_2 - 1] -
                  (com_S[problem::axis_2 - 1] / com_S_mag) * vel_S[problem::axis_1 - 1];
    }
#else
    if (mass_P > 0.0) {
      vel_P_rad = vel_P[problem::axis_1 - 1];
      vel_P_phi = vel_P[problem::axis_3 - 1];
    }

    if (mass_S > 0.0) {
      vel_S_rad = vel_S[problem::axis_1 - 1];
      vel_S_phi = vel_S[problem::axis_3 - 1];
    }
#endif

    if (mass_P > 0.0 && mass_S > 0.0) {

      // Calculate the distance between the primary and secondary.

      for (int i = 0; i < 3; i++) {
          wd_dist[i] = com_S[i] - com_P[i];
      }

      separation = norm(wd_dist);

      // Calculate the angle between the initial stellar axis and
      // the line currently joining the two stars. Note that this
      // neglects any motion in the plane perpendicular to the initial orbit.

      angle = std::atan2(wd_dist[problem::axis_2 - 1] - wd_dist_init[problem::axis_2 - 1],
                         wd_dist[problem::axis_1 - 1] - wd_dist_init[problem::axis_1 - 1]) * 180.0 / M_PI;

      // Now let's transform from [-180, 180] to [0, 360].

      if (angle < 0.0) angle += 360.0;

    }

    // Write data out to the log.

    if (amrex::ParallelDescriptor::IOProcessor())
    {

      // The data logs are only defined on the IO processor
      // for parallel runs, so the stream should only be opened inside.

      if (parent->NumDataLogs() > 0) {

          std::ostream& log = parent->DataLog(0);

          if (log.good()) {

              if (time == 0.0) {

                  // Output the git commit hashes used to build the executable.

                  writeGitHashes(log);

                  int n = 0;

                  std::ostringstream header;

                  header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                  header << std::setw(fixwidth) << "                     TIME"; ++n;

                  header << std::setw(datwidth) << "              WD DISTANCE"; ++n;
                  header << std::setw(fixwidth) << "                 WD ANGLE"; ++n;
                  header << std::setw(datwidth) << "                     ADOT"; ++n;
                  header << std::setw(datwidth) << "                     MDOT"; ++n;

                  header << std::endl;

                  log << std::setw(intwidth) << "#   COLUMN 1";
                  log << std::setw(fixwidth) << "                        2";

                  for (int i = 3; i <= n; ++i) {
                      log << std::setw(datwidth) << i;
                  }

                  log << std::endl;

                  log << header.str();

              }

              log << std::fixed;

              log << std::setw(intwidth)                                     << timestep;
              log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

              log << std::scientific;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << separation;

              log << std::fixed;
              log << std::setw(fixwidth) << std::setprecision(dataprecision) << angle;

              log << std::scientific;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_rad + vel_S_rad;

              log << std::scientific;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot;

              log << std::endl;

          }
      }

      // Primary star

      if (parent->NumDataLogs() > 1) {

          std::ostream& log = parent->DataLog(1);

          if (log.good()) {

              if (time == 0.0) {

                  // Output the git commit hashes used to build the executable.

                  writeGitHashes(log);

                  int n = 0;

                  std::ostringstream header;

                  header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                  header << std::setw(fixwidth) << "                     TIME"; ++n;
                  header << std::setw(datwidth) << "             PRIMARY MASS"; ++n;
                  header << std::setw(datwidth) << "             PRIMARY MDOT"; ++n;
                  header << std::setw(datwidth) << "          PRIMARY MAG COM"; ++n;
#if (AMREX_SPACEDIM == 3)
                  header << std::setw(datwidth) << "            PRIMARY X COM"; ++n;
                  header << std::setw(datwidth) << "            PRIMARY Y COM"; ++n;
                  header << std::setw(datwidth) << "            PRIMARY Z COM"; ++n;
#else
                  header << std::setw(datwidth) << "            PRIMARY R COM"; ++n;
                  header << std::setw(datwidth) << "            PRIMARY Z COM"; ++n;
#endif
                  header << std::setw(datwidth) << "          PRIMARY MAG VEL"; ++n;
                  header << std::setw(datwidth) << "          PRIMARY RAD VEL"; ++n;
                  header << std::setw(datwidth) << "          PRIMARY ANG VEL"; ++n;
#if (AMREX_SPACEDIM == 3)
                  header << std::setw(datwidth) << "            PRIMARY X VEL"; ++n;
                  header << std::setw(datwidth) << "            PRIMARY Y VEL"; ++n;
                  header << std::setw(datwidth) << "            PRIMARY Z VEL"; ++n;
#else
                  header << std::setw(datwidth) << "            PRIMARY R VEL"; ++n;
                  header << std::setw(datwidth) << "            PRIMARY Z VEL"; ++n;
#endif
                  header << std::setw(datwidth) << "       PRIMARY T_FREEFALL"; ++n;
                  header << std::setw(datwidth) << "     PRIMARY ROCHE RADIUS"; ++n;
                  for (int i = 0; i <= 6; ++i) {
                      header << "       PRIMARY 1E" << i << " RADIUS";          ++n;
                  }

                  header << std::endl;

                  log << std::setw(intwidth) << "#   COLUMN 1";
                  log << std::setw(fixwidth) << "                        2";

                  for (int i = 3; i <= n; ++i) {
                      log << std::setw(datwidth) << i;
                  }

                  log << std::endl;

                  log << header.str();

              }

              log << std::fixed;

              log << std::setw(intwidth)                                     << timestep;
              log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

              log << std::scientific;

              log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_P;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot_P;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P_mag;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P[0];
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P[1];
#if (AMREX_SPACEDIM == 3)
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P[2];
#endif
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_mag;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_rad;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_phi;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P[0];
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P[1];
#if (AMREX_SPACEDIM == 3)
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P[2];
#endif
              log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_P;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << roche_rad_P;
              for (int i = 0; i <= 6; ++i) {
                  log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_P[i];
              }

              log << std::endl;

          }

      }

      // Secondary star

      if (parent->NumDataLogs() > 2) {

          std::ostream& log = parent->DataLog(2);

          if (log.good()) {

              if (time == 0.0) {

                  // Output the git commit hashes used to build the executable.

                  writeGitHashes(log);

                  int n = 0;

                  std::ostringstream header;

                  header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                  header << std::setw(fixwidth) << "                     TIME"; ++n;
                  header << std::setw(datwidth) << "           SECONDARY MASS"; ++n;
                  header << std::setw(datwidth) << "           SECONDARY MDOT"; ++n;
                  header << std::setw(datwidth) << "        SECONDARY MAG COM"; ++n;
#if (AMREX_SPACEDIM == 3)
                  header << std::setw(datwidth) << "          SECONDARY X COM"; ++n;
                  header << std::setw(datwidth) << "          SECONDARY Y COM"; ++n;
                  header << std::setw(datwidth) << "          SECONDARY Z COM"; ++n;
#else
                  header << std::setw(datwidth) << "          SECONDARY R COM"; ++n;
                  header << std::setw(datwidth) << "          SECONDARY Z COM"; ++n;
#endif
                  header << std::setw(datwidth) << "        SECONDARY MAG VEL"; ++n;
                  header << std::setw(datwidth) << "        SECONDARY RAD VEL"; ++n;
                  header << std::setw(datwidth) << "        SECONDARY ANG VEL"; ++n;
#if (AMREX_SPACEDIM == 3)
                  header << std::setw(datwidth) << "          SECONDARY X VEL"; ++n;
                  header << std::setw(datwidth) << "          SECONDARY Y VEL"; ++n;
                  header << std::setw(datwidth) << "          SECONDARY Z VEL"; ++n;
#else
                  header << std::setw(datwidth) << "          SECONDARY R VEL"; ++n;
                  header << std::setw(datwidth) << "          SECONDARY Z VEL"; ++n;
#endif
                  header << std::setw(datwidth) << "     SECONDARY T_FREEFALL"; ++n;
                  header << std::setw(datwidth) << "   SECONDARY ROCHE RADIUS"; ++n;
                  for (int i = 0; i <= 6; ++i) {
                      header << "     SECONDARY 1E" << i << " RADIUS";          ++n;
                  }

                  header << std::endl;

                  log << std::setw(intwidth) << "#   COLUMN 1";
                  log << std::setw(fixwidth) << "                        2";

                  for (int i = 3; i <= n; ++i) {
                      log << std::setw(datwidth) << i;
                  }

                  log << std::endl;

                  log << header.str();

              }

              log << std::fixed;

              log << std::setw(intwidth)                                     << timestep;
              log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

              log << std::scientific;

              log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_S;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot_S;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S_mag;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S[0];
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S[1];
#if (AMREX_SPACEDIM == 3)
              log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S[2];
#endif
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S_mag;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S_rad;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S_phi;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S[0];
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S[1];
#if (AMREX_SPACEDIM == 3)
              log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S[2];
#endif
              log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_S;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << roche_rad_S;
              for (int i = 0; i <= 6; ++i) {
                  log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_S[i];
              }

              log << std::endl;

          }

      }

      // Rotation period over time

      if (parent->NumDataLogs() > 3) {

          std::ostream& log = parent->DataLog(3);

          if (log.good()) {

              if (time == 0.0) {

                  // Output the git commit hashes used to build the executable.

                  writeGitHashes(log);

                  int n = 0;

                  std::ostringstream header;

                  header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                  header << std::setw(fixwidth) << "                     TIME"; ++n;
                  header << std::setw(datwidth) << "          ROTATION PERIOD"; ++n;
                  header << std::setw(datwidth) << "       ROTATION FREQUENCY"; ++n;

                  header << std::endl;

                  log << std::setw(intwidth) << "#   COLUMN 1";
                  log << std::setw(fixwidth) << "                        2";

                  for (int i = 3; i <= n; ++i) {
                      log << std::setw(datwidth) << i;
                  }

                  log << std::endl;

                  log << header.str();

              }

              log << std::fixed;

              log << std::setw(intwidth)                                     << timestep;
              log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

              log << std::scientific;

              log << std::setw(datwidth) << std::setprecision(dataprecision) << rotational_period;
              log << std::setw(datwidth) << std::setprecision(dataprecision) << (2.0 * M_PI / rotational_period);

              log << std::endl;

          }

      }

    }
}
