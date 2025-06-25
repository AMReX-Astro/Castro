#include <iomanip>

#include <Castro.H>
#include <Castro_util.H>

#ifdef GRAVITY
#include <Gravity.H>
#endif

#ifdef ROTATION
#include <Rotation.H>
#endif

using namespace amrex;

Real
Castro::volWgtSum (const std::string& name, Real time, bool local, bool finemask)
{
    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    return volWgtSum(*mf, 0, local, finemask);
}

Real
Castro::volWgtSum (const MultiFab& mf, int comp, bool local, bool finemask)
{
    BL_PROFILE("Castro::volWgtSum()");

    bool mask_available = level < parent->finestLevel() && finemask;

    MultiFab tmp_mf;
    const MultiFab& mask_mf = mask_available ? getLevel(level+1).build_fine_mask() : tmp_mf;

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = mf[mfi].array(comp);
        auto const& vol = volume.array(mfi);
        auto const& mask = mask_available ? mask_mf.array(mfi) : Array4<Real>{};

        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real maskFactor = mask_available ? mask(i,j,k) : 1.0_rt;

            return {fab(i,j,k) * vol(i,j,k) * maskFactor};
        });

    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) {
        ParallelDescriptor::ReduceRealSum(sum);
    }

    return sum;
}

Real
Castro::locWgtSum (const std::string& name, Real time, int idir, bool local)
{
    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    return locWgtSum(*mf, 0, idir, local);
}

Real
Castro::locWgtSum (const MultiFab& mf, int comp, int idir, bool local)
{
    BL_PROFILE("Castro::locWgtSum()");

    bool mask_available = level < parent->finestLevel();

    MultiFab tmp_mf;
    const MultiFab& mask_mf = mask_available ? getLevel(level+1).build_fine_mask() : tmp_mf;

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    auto dx     = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = mf[mfi].array(comp);
        auto const& vol = volume.array(mfi);
        auto const& mask = mask_available ? mask_mf.array(mfi) : Array4<Real>{};

        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real maskFactor = mask_available ? mask(i,j,k) : 1.0_rt;

            Real loc[3];

            loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
            loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
            loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
            loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
            loc[2] = 0.0_rt;
#endif

            Real ds;

            if (idir == 0) { // sum(mass * x)
                ds = fab(i,j,k) * vol(i,j,k) * maskFactor * loc[0];
            }
            else if (idir == 1) { // sum(mass * y)
                ds = fab(i,j,k) * vol(i,j,k) * maskFactor * loc[1];
            }
            else { // sum(mass * z)
                ds = fab(i,j,k) * vol(i,j,k) * maskFactor * loc[2];
            }

            return {ds};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) {
        ParallelDescriptor::ReduceRealSum(sum);
    }

    return sum;
}

Real
Castro::volProductSum (const std::string& name1,
                       const std::string& name2,
                       Real time, bool local)
{
    auto mf1 = derive(name1, time, 0);
    auto mf2 = derive(name2, time, 0);

    BL_ASSERT(mf1);
    BL_ASSERT(mf2);

    return volProductSum(*mf1, *mf2, 0, 0, local);
}

Real
Castro::volProductSum (const MultiFab& mf1,
                       const MultiFab& mf2,
                       int comp1, int comp2,
                       bool local)
{
    BL_PROFILE("Castro::volProductSum()");

    bool mask_available = level < parent->finestLevel();

    MultiFab tmp_mf;
    const MultiFab& mask_mf = mask_available ? getLevel(level+1).build_fine_mask() : tmp_mf;

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf1, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab1 = mf1[mfi].array(comp1);
        auto const& fab2 = mf2[mfi].array(comp2);
        auto const& vol  = volume.array(mfi);
        auto const& mask = mask_available ? mask_mf.array(mfi) : Array4<Real>{};

        const Box& box = mfi.tilebox();

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real maskFactor = mask_available ? mask(i,j,k) : 1.0_rt;

            return {fab1(i,j,k) * fab2(i,j,k) * vol(i,j,k) * maskFactor};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) {
        ParallelDescriptor::ReduceRealSum(sum);
    }

    return sum;
}

Real
Castro::locSquaredSum (const std::string& name,
                       Real               time,
                       int                idir,
                       bool               local)
{
    BL_PROFILE("Castro::locSquaredSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    bool mask_available = level < parent->finestLevel();

    MultiFab tmp_mf;
    const MultiFab& mask_mf = mask_available ? getLevel(level+1).build_fine_mask() : tmp_mf;

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    auto dx     = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = (*mf).array(mfi);
        auto const& mask = mask_available ? mask_mf.array(mfi) : Array4<Real>{};

        const Box& box = mfi.tilebox();

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real maskFactor = mask_available ? mask(i,j,k) : 1.0_rt;

            Real loc[3];

            loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
            loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
            loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
            loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
            loc[2] = 0.0_rt;
#endif

            Real ds;

            if (idir == 0) { // sum(mass * x^2)
                ds = fab(i,j,k) * maskFactor * loc[0] * loc[0];
            }
            else if (idir == 1) { // sum(mass * y^2)
                ds = fab(i,j,k) * maskFactor * loc[1] * loc[1];
            }
            else if (idir == 2) { // sum(mass * z^2)
                ds = fab(i,j,k) * maskFactor * loc[2] * loc[2];
            }
            else { // sum(mass * r^2)
                ds = fab(i,j,k) * maskFactor * (loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
            }

            return {ds};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) {
        ParallelDescriptor::ReduceRealSum(sum);
    }

    return sum;
}



#ifdef GRAVITY
void
Castro::gwstrain (Real time,
                  Real& h_plus_1, Real& h_cross_1,
                  Real& h_plus_2, Real& h_cross_2,
                  Real& h_plus_3, Real& h_cross_3,
                  bool local) {

    amrex::ignore_unused(time);

    BL_PROFILE("Castro::gwstrain()");

    // We have nothing to do if the user did not request the gravitational wave
    // strain (inferred from whether the observation distance is positive).
    if (castro::gw_dist <= 0.0_rt) {
        return;
    }

    GeometryData geomdata = geom.data();

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);

    bool mask_available = level < parent->finestLevel();

    MultiFab tmp_mf;
    const MultiFab& mask_mf = mask_available ? getLevel(level+1).build_fine_mask() : tmp_mf;

    // Qtt stores the second time derivative of the quadrupole moment.
    // We calculate it directly rather than computing the quadrupole moment
    // and differentiating it in time, because the latter method is less accurate
    // and requires the state at other timesteps. See, e.g., Equation 5 of
    // Loren-Aguilar et al. 2005.

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum,
              ReduceOpSum, ReduceOpSum, ReduceOpSum,
              ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Real,
               Real, Real, Real,
               Real, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& box = mfi.tilebox();

        auto rho = S_new[mfi].array(URHO);
        auto vol = volume.array(mfi);
        auto xmom = S_new[mfi].array(UMX);
        auto ymom = S_new[mfi].array(UMY);
        auto zmom = S_new[mfi].array(UMZ);
        auto gravx = grav_new[mfi].array(0);
        auto gravy = grav_new[mfi].array(1);
        auto gravz = grav_new[mfi].array(2);
        auto const& mask = mask_available ? mask_mf.array(mfi) : Array4<Real>{};

        // Calculate the second time derivative of the quadrupole moment tensor,
        // according to the formula in Equation 6.5 of Blanchet, Damour and Schafer 1990.
        // It involves integrating the mass distribution and then taking the symmetric
        // trace-free part of the tensor. We can do the latter operation here since the
        // integral is a linear operator and each part of the domain contributes independently.

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real maskFactor = mask_available ? mask(i,j,k) : 1.0_rt;

            GpuArray<Real, 3> r;
            position(i, j, k, geomdata, r);

            for (int n = 0; n < 3; ++n) {
                r[n] -= problem::center[n];
            }

            Real rhoInv;
            if (rho(i,j,k) * maskFactor > 0.0_rt) {
                rhoInv = 1.0_rt / rho(i,j,k);
            } else {
                rhoInv = 0.0_rt;
            }

            // Account for rotation, if there is any. These will leave
            // r and vel and changed, if not.

            GpuArray<Real, 3> pos{r};
#ifdef ROTATION
            auto omega = get_omega_vec(geomdata, j);
            pos = inertial_rotation(r, omega, time);
#endif

            // For constructing the velocity in the inertial frame, we need to
            // account for the fact that we have rotated the system already, so that
            // the r in omega x r is actually the position in the inertial frame, and
            // not the usual position in the rotating frame. It has to be on physical
            // grounds, because for binary orbits where the stars aren't moving, that
            // r never changes, and so the contribution from rotation would never change.
            // But it must, since the motion vector of the stars changes in the inertial
            // frame depending on where we are in the orbit.

            GpuArray<Real, 3> vel;
            vel[0] = xmom(i,j,k) * maskFactor * rhoInv;
            vel[1] = ymom(i,j,k) * maskFactor * rhoInv;
            vel[2] = zmom(i,j,k) * maskFactor * rhoInv;

            GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
            rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif

            GpuArray<Real, 3> g;
            g[0] = gravx(i,j,k) * maskFactor;
            g[1] = gravy(i,j,k) * maskFactor;
            g[2] = gravz(i,j,k) * maskFactor;

            // We need to rotate the gravitational field to be consistent with the rotated position.

            GpuArray<Real, 3> inertial_g{g};
#ifdef ROTATION
            inertial_g = inertial_rotation(g, omega, time);
#endif

            // Absorb the factor of 2 outside the integral into the zone mass, for efficiency.

            Real dM = 2.0_rt * rho(i,j,k) * vol(i,j,k) * maskFactor;

            Array2D<Real, 0, 2, 0, 2> dQtt{0.0};

            if (AMREX_SPACEDIM == 3) {

                for (int m = 0; m < 3; ++m) {
                    for (int l = 0; l < 3; ++l) {
                        dQtt(l,m) += dM * (inertial_vel[l] * inertial_vel[m] + pos[l] * inertial_g[m]);
                    }
                }

            } else {

                // For axisymmetric coordinates we need to be careful here.
                // We want to calculate the quadrupole tensor in terms of
                // Cartesian coordinates but our coordinates are cylindrical (R, z).
                // What we can do is to first express the Cartesian coordinates
                // as (x, y, z) = (R cos(phi), R sin(phi), z). Then we can integrate
                // out the phi coordinate for each component. The off-diagonal components
                // all then vanish automatically. The on-diagonal components xx and yy
                // pick up a factor of cos**2(phi) which when integrated from (0, 2*pi)
                // yields pi. Note that we're going to choose that the cylindrical z axis
                // coincides with the Cartesian x-axis, which is our default choice.

                // We also need to then divide by the volume by 2*pi since
                // it has already been integrated out.

                dM /= (2.0_rt * M_PI);

                dQtt(0,0) += dM * (2.0_rt * M_PI) * (inertial_vel[1] * inertial_vel[1] + pos[1] * inertial_g[1]);
                dQtt(1,1) += dM * M_PI * (inertial_vel[0] * inertial_vel[0] + pos[0] * g[0]);
                dQtt(2,2) += dM * M_PI * (inertial_vel[0] * inertial_vel[0] + pos[0] * g[0]);

            }

            // Now take the symmetric trace-free part of the quadrupole moment.
            // The operator is defined in Equation 6.7 of Blanchet et al. (1990):
            // STF(A^{ij}) = 1/2 A^{ij} + 1/2 A^{ji} - 1/3 delta^{ij} sum_{k} A^{kk}.

            Array2D<Real, 0, 2, 0, 2> dQ;

            for (int l = 0; l < 3; ++l) {
                for (int m = 0; m < 3; ++m) {
                    dQ(l,m) = 0.5_rt * dQtt(l,m) + 0.5_rt * dQtt(m,l);
                    if (l == m) {
                        dQ(l,m) -= (1.0_rt / 3.0_rt) * dQtt(m,m);
                    }
                }
            }

            return {dQ(0,0), dQ(1,0), dQ(2,0),
                    dQ(0,1), dQ(1,1), dQ(2,1),
                    dQ(0,2), dQ(1,2), dQ(2,2)};
        });
    }

    Array2D<Real, 0, 2, 0, 2> Qtt;

    ReduceTuple hv = reduce_data.value();
    Qtt(0,0) = amrex::get<0>(hv);
    Qtt(1,0) = amrex::get<1>(hv);
    Qtt(2,0) = amrex::get<2>(hv);
    Qtt(0,1) = amrex::get<3>(hv);
    Qtt(1,1) = amrex::get<4>(hv);
    Qtt(2,1) = amrex::get<5>(hv);
    Qtt(0,2) = amrex::get<6>(hv);
    Qtt(1,2) = amrex::get<7>(hv);
    Qtt(2,2) = amrex::get<8>(hv);

    // Now, do a global reduce over all processes.

    if (!local) {
        amrex::ParallelDescriptor::ReduceRealSum(Qtt.begin(), 9);
    }

    // Now that we have the second time derivative of the quadrupole
    // tensor, we can calculate the transverse-trace gauge strain tensor.

    // Standard Kronecker delta.

    Real delta[3][3] = {{0.0}};

    for (int i = 0; i < 3; ++i) {
        delta[i][i] = 1.0;
    }

    // Unit vector for the wave is simply the distance
    // vector to the observer normalized by the total distance.
    // We are going to repeat this process by looking along
    // all three coordinate axes.

    for (int dir = 0; dir < 3; ++dir) {

        Real dist[3] = {0.0};
        dist[dir] = castro::gw_dist;

        Real r = std::sqrt(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);

        Real n[3] = {dist[0] / r, dist[1] / r, dist[2] / r};

        // Projection operator onto the unit vector n.

        Real proj[3][3][3][3] = {{0.0}};

        for (int l = 0; l < 3; ++l) {
            for (int k = 0; k < 3; ++k) {
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        proj[l][k][j][i] = (delta[k][i] - n[i] * n[k]) * (delta[l][j] - n[j] * n[l]) -
                                            0.5_rt * (delta[j][i] - n[i] * n[j]) * (delta[l][k] - n[k] * n[l]);
                    }
                }
            }
        }

        // Now we can calculate the strain tensor.

        Real h[3][3] = {{0.0}};

        for (int l = 0; l < 3; ++l) {
            for (int k = 0; k < 3; ++k) {
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        h[j][i] += proj[l][k][j][i] * Qtt(k, l);
                    }
                }
            }
        }
        // Finally multiply by the coefficients.

        r *= C::parsec * 1.e3_rt; // Convert from kpc to cm

        for (int j = 0; j < 3; ++j) {   // NOLINT(modernize-loop-convert)
            for (int i = 0; i < 3; ++i) {
                h[j][i] *= 2.0_rt * C::Gconst / (std::pow(C::c_light, 4) * r);
            }
        }

        // We are adding here so that this calculation makes sense on multiple levels.

        if (dir == 0) {

            h_plus_1  += h[1][1];
            h_cross_1 += h[2][1];

        }
        else if (dir == 1) {

            h_plus_2  += h[2][2];
            h_cross_2 += h[0][2];

        }
        else if (dir == 2) {

            h_plus_3  += h[0][0];
            h_cross_3 += h[1][0];

        }

    }
}
#endif
