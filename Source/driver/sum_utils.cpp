#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>

#ifdef GRAVITY
#include <Gravity.H>
#include <Gravity_F.H>
#endif

#ifdef ROTATION
#include <Rotation.H>
#endif

using namespace amrex;

Real
Castro::volWgtSum (const std::string& name,
                   Real               time,
                   bool               local,
                   bool               finemask)
{
    BL_PROFILE("Castro::volWgtSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel() && finemask)
    {
        const MultiFab& mask = getLevel(level+1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = (*mf).array(mfi);
        auto const& vol = volume.array(mfi);

        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return {fab(i,j,k) * vol(i,j,k)};
        });

    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSquaredSum (const std::string& name,
                          Real               time,
                          bool               local)
{
    BL_PROFILE("Castro::volWgtSquaredSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
        const MultiFab& mask = getLevel(level+1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = (*mf).array(mfi);
        auto const& vol = volume.array(mfi);
    
        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return {fab(i,j,k) * fab(i,j,k) * vol(i,j,k)};
        });

    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSum (const std::string& name,
                   Real               time,
                   int                idir,
                   bool               local)
{
    BL_PROFILE("Castro::locWgtSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
        const MultiFab& mask = getLevel(level+1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    auto dx     = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = (*mf).array(mfi);
    
        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
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
                ds = fab(i,j,k) * loc[0];
            }
            else if (idir == 1) { // sum(mass * y)
                ds = fab(i,j,k) * loc[1];
            }
            else { // sum(mass * z)
                ds = fab(i,j,k) * loc[2];
            }

            return {ds};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volProductSum (const std::string& name1, 
                       const std::string& name2,
                       Real time, bool local)
{
    BL_PROFILE("Castro::volProductSum()");

    auto mf1 = derive(name1, time, 0);
    auto mf2 = derive(name2, time, 0);

    BL_ASSERT(mf1);
    BL_ASSERT(mf2);

    if (level < parent->finestLevel())
    {
        const MultiFab& mask = getLevel(level+1).build_fine_mask();
        MultiFab::Multiply(*mf1, mask, 0, 0, 1, 0);
        MultiFab::Multiply(*mf2, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(*mf1, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab1 = (*mf1).array(mfi);
        auto const& fab2 = (*mf2).array(mfi);
        auto const& vol  = volume.array(mfi);
    
        const Box& box = mfi.tilebox();

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return {fab1(i,j,k) * fab2(i,j,k) * vol(i,j,k)};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

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

    if (level < parent->finestLevel())
    {
        const MultiFab& mask = getLevel(level+1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    auto dx     = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto const& fab = (*mf).array(mfi);
    
        const Box& box = mfi.tilebox();

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
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
                ds = fab(i,j,k) * loc[0] * loc[0];
            }
            else if (idir == 1) { // sum(mass * y^2)
                ds = fab(i,j,k) * loc[1] * loc[1];
            }
            else if (idir == 2) { // sum(mass * z^2)
                ds = fab(i,j,k) * loc[2] * loc[2];
            }
            else { // sum(mass * r^2)
                ds = fab(i,j,k) * (loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
            }

            return {ds};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}



#ifdef GRAVITY
void
Castro::gwstrain (Real time,
		  Real& h_plus_1, Real& h_cross_1,
		  Real& h_plus_2, Real& h_cross_2,
		  Real& h_plus_3, Real& h_cross_3,
		  bool local) {

    BL_PROFILE("Castro::gwstrain()");

    // We have nothing to do if the user did not request the gravitational wave
    // strain (inferred from whether the observation distance is positive).
    if (castro::gw_dist <= 0.0_rt) {
        return;
    }

    GeometryData geomdata = geom.data();

    auto mfrho   = derive("density",time,0);
    auto mfxmom  = derive("xmom",time,0);
    auto mfymom  = derive("ymom",time,0);
    auto mfzmom  = derive("zmom",time,0);
    auto mfgravx = derive("grav_x",time,0);
    auto mfgravy = derive("grav_y",time,0);
    auto mfgravz = derive("grav_z",time,0);

    BL_ASSERT(mfrho   != nullptr);
    BL_ASSERT(mfxmom  != nullptr);
    BL_ASSERT(mfymom  != nullptr);
    BL_ASSERT(mfzmom  != nullptr);
    BL_ASSERT(mfgravx != nullptr);
    BL_ASSERT(mfgravy != nullptr);
    BL_ASSERT(mfgravz != nullptr);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();

	MultiFab::Multiply(*mfrho,   mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfxmom,  mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfymom,  mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfzmom,  mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravx, mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravy, mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravz, mask, 0, 0, 1, 0);
    }

    // Qtt stores the second time derivative of the quadrupole moment.
    // We calculate it directly rather than computing the quadrupole moment
    // and differentiating it in time, because the latter method is less accurate
    // and requires the state at other timesteps. See, e.g., Equation 5 of
    // Loren-Aguilar et al. 2005.

    // It is a 3x3 rank-2 tensor, but AMReX expects IntVect() to use BL_SPACEDIM
    // dimensions, so we add a redundant third index in 3D.

    Box bx( IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(2, 2, 0)) );

    FArrayBox Qtt(bx);

    Qtt.setVal<RunOn::Device>(0.0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    Vector< std::unique_ptr<FArrayBox> > priv_Qtt(nthreads);
    for (int i=0; i<nthreads; i++) {
	priv_Qtt[i].reset(new FArrayBox(bx));
    }
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
	priv_Qtt[tid]->setVal<RunOn::Device>(0.0);
#endif
	for (MFIter mfi(*mfrho, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	    const Box& bx = mfi.tilebox();

            auto rho = (*mfrho).array(mfi);
            auto vol = volume.array(mfi);
            auto xmom = (*mfxmom).array(mfi);
            auto ymom = (*mfymom).array(mfi);
            auto zmom = (*mfzmom).array(mfi);
            auto gravx = (*mfgravx).array(mfi);
            auto gravy = (*mfgravy).array(mfi);
            auto gravz = (*mfgravz).array(mfi);

            // Calculate the second time derivative of the quadrupole moment tensor,
            // according to the formula in Equation 6.5 of Blanchet, Damour and Schafer 1990.
            // It involves integrating the mass distribution and then taking the symmetric
            // trace-free part of the tensor. We can do the latter operation here since the
            // integral is a linear operator and each part of the domain contributes independently.

#ifdef _OPENMP
            auto Qtt_arr = priv_Qtt[tid]->array();
#else
            auto Qtt_arr = Qtt.array();
#endif

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                Array2D<Real, 0, 2, 0, 2> dQtt{};

                GpuArray<Real, 3> r;
                position(i, j, k, geomdata, r);

                for (int n = 0; n < 3; ++n) {
                    r[n] -= problem::center[n];
                }

                Real rhoInv;
                if (rho(i,j,k) > 0.0_rt) {
                    rhoInv = 1.0_rt / rho(i,j,k);
                } else {
                    rhoInv = 0.0_rt;
                }

                // Account for rotation, if there is any. These will leave
                // r and vel and changed, if not.

                GpuArray<Real, 3> pos{r};
#ifdef ROTATION
                pos = inertial_rotation(r, time);
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
                vel[0] = xmom(i,j,k) * rhoInv;
                vel[1] = ymom(i,j,k) * rhoInv;
                vel[2] = zmom(i,j,k) * rhoInv;

                GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
                rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif

                GpuArray<Real, 3> g;
                g[0] = gravx(i,j,k);
                g[1] = gravy(i,j,k);
                g[2] = gravz(i,j,k);

                // We need to rotate the gravitational field to be consistent with the rotated position.

                GpuArray<Real, 3> inertial_g{g};
#ifdef ROTATION
                inertial_g = inertial_rotation(g, time);
#endif

                // Absorb the factor of 2 outside the integral into the zone mass, for efficiency.

                Real dM = 2.0_rt * rho(i,j,k) * vol(i,j,k);

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

                for (int l = 0; l < 3; ++l) {
                    for (int m = 0; m < 3; ++m) {

                        Real dQ = 0.5_rt * dQtt(l,m) + 0.5_rt * dQtt(m,l);
                        if (l == m) {
                            dQ -= (1.0_rt / 3.0_rt) * dQtt(m,m);
                        }

                        Gpu::Atomic::Add(&Qtt_arr(l,m,0), dQ);

                    }
                }
            });
        }
    }

    // Do an OpenMP reduction on the tensor.

#ifdef _OPENMP
    int n = bx.numPts();
    Real* p = Qtt.dataPtr();
#pragma omp barrier
#pragma omp for nowait
    for (int i=0; i<n; ++i)
    {
        for (int it=0; it<nthreads; it++) {
            const Real* pq = priv_Qtt[it]->dataPtr();
            p[i] += pq[i];
        }
    }
#endif

    // Now, do a global reduce over all processes.

    if (!local) {
        amrex::ParallelDescriptor::ReduceRealSum(Qtt.dataPtr(),bx.numPts());
    }

    // Now that we have the second time derivative of the quadrupole
    // tensor, we can calculate the transverse-trace gauge strain tensor.

    // Standard Kronecker delta.

    Real delta[3][3] = {0.0};

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

        Real proj[3][3][3][3] = {0.0};

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

        Real h[3][3] = {0.0};

        for (int l = 0; l < 3; ++l) {
            for (int k = 0; k < 3; ++k) {
                for (int j = 0; j < 3; ++j) {
                    for (int i = 0; i < 3; ++i) {
                        h[j][i] += proj[l][k][j][i] * Qtt.array()(k, l, 0);
                    }
                }
            }
        }
        // Finally multiply by the coefficients.

        r *= C::parsec * 1.e3_rt; // Convert from kpc to cm

        for (int j = 0; j < 3; ++j) {
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
