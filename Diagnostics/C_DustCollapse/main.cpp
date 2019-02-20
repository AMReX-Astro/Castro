
#include <iostream>
#include "AMReX_DataServices.H"
#include <DustCollapse_F.H>

using namespace amrex;

std::string inputs_name = "";

template <typename T>
Vector<size_t> sort_indexes(const Vector<T> &v);

int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv, false);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	// Input arguments
	if (argc < 2)
		Abort("ERROR: Missing plotfiles");

	for (auto i = 1; i < argc; i++) {

		string pltfile = argv[i];

		// Start dataservices (no clue why we need to do this)
		DataServices::SetBatchMode();

		// Define the type of file
		Amrvis::FileType fileType(Amrvis::NEWPLT);
		DataServices dataServices (pltfile, fileType);

		if (!dataServices.AmrDataOk())
			DataServices::Dispatch(DataServices::ExitRequest, NULL);

		// get data from plot file
		AmrData& data = dataServices.AmrDataRef();

		int finestLevel = data.FinestLevel();

		// get variable names
		const Vector<string>& varNames = data.PlotVarNames();

		// get the index bounds and dx.
		Box domain = data.ProbDomain()[finestLevel];

		auto dx = data.CellSize(finestLevel);

		const Vector<Real>& problo = data.ProbLo();
		const Vector<Real>& probhi = data.ProbHi();
		Vector<int> rr = data.RefRatio();

        auto rmin = problo[0];

		double dx_fine = *(std::min_element(dx.begin(), dx.end()));

		int nbins = domain.length(0);

		Vector<Real> r(nbins);

		for (auto i = 0; i < nbins; i++)
			r[i] = (i + 0.5) * dx_fine;

		// find variable indices
		auto dens_comp = data.StateNumber("density");

		if (dens_comp < 0 )
			Abort("ERROR: density variable not found");

		// allocate storage for data
		Vector<Real> dens(nbins, 0.);

		auto r1 = 1.0;

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++)
			fill_comps[i] = i;

		// imask will be set to false if we've already output the data.
		// Note, imask is defined in terms of the finest level.  As we loop
		// over levels, we will compare to the finest level index space to
		// determine if we've already output here
		int mask_size = nbins;
		for (auto i = 0; i < finestLevel - 1; i++)
			mask_size *= rr[i];

		Vector<int> imask(pow(mask_size, AMREX_SPACEDIM), 1);

        int cnt = 0;

		// extract the 1d data
		for (int l = finestLevel; l >= 0; l--) {

			Vector<Real> level_dx = data.DxLevel()[l];

			const BoxArray& ba = data.boxArray(l);
			const DistributionMapping& dm = data.DistributionMap(l);

			MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
			data.FillVar(lev_data_mf, l, varNames, fill_comps);

#ifdef _OPENMP
#pragma omp parallel
#endif
			for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
				const Box& bx = mfi.tilebox();

#if (AMREX_SPACEDIM == 1)
				fdustcollapse1d(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				           BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				           nbins, dens.dataPtr(),
				           imask.dataPtr(), mask_size, r1, dens_comp, &cnt);
#elif (AMREX_SPACEDIM == 2)

				fdustcollapse2d(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				               nbins, dens_bin.dataPtr(),
				               vel_bin.dataPtr(), pres_bin.dataPtr(),
				               e_bin.dataPtr(), ncount.dataPtr(),
				               imask.dataPtr(), mask_size, r1,
				               dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp,
				               dx_fine, level_dx.dataPtr(),
				               yctr);
#else
				fdustcollapse3d(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				               nbins, dens_bin.dataPtr(),
				               vel_bin.dataPtr(), pres_bin.dataPtr(),
				               ncount.dataPtr(),
				               imask.dataPtr(), mask_size, r1,
				               dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp,
				               dx_fine, level_dx.dataPtr(), xctr, yctr);
#endif
			}

			// adjust r1 for the next lowest level
			if (l != 0) r1 *= rr[l-1];
		}

        // sort the data based on the coordinates
		auto isv = sort_indexes(r);

        // These are calculated analytically given initial density 1.e9 and the
        // analytic expression for the radius as a function of time t = 0.00
        Real max_dens = 1.e9;

        if (fabs(data.Time()) <= 1.e-8)
            max_dens = 1.e9;
        else if (fabs(data.Time() - 0.01) <= 1.e-8)
            max_dens = 1.043345e9;
        else if (fabs(data.Time() - 0.02) <= 1.e-8)
            max_dens = 1.192524e9;
        else if (fabs(data.Time() - 0.03) <= 1.e-8)
            max_dens = 1.527201e9;
        else if (fabs(data.Time() - 0.04) <= 1.e-8)
            max_dens = 2.312884e9;
        else if (fabs(data.Time() - 0.05) <= 1.e-8)
            max_dens = 4.779133e9;
        else if (fabs(data.Time() - 0.06) <= 1.e-8)
            max_dens = 24.472425e9;
        else if (fabs(data.Time() - 0.065) <= 1.e-8)
            max_dens = 423.447291e9;
        else {
            Print() << "Dont know the maximum density at this time: " << data.Time() <<std::endl;
            Abort();
        }

        // loop over the solution, from r = 0 outward, and find the first
        // place where the density drops below the threshold density
        auto index = -1;
        for (auto j = 0; j < cnt; j++) {
            if (dens[isv[j]] < 0.5 * max_dens) {
                index = i;
                break;
            }
        }

        Real r_interface = 0.0;

        if (index < 0)
            Abort("ERROR: density never fell below threshold");
        else if (index < 2)
            r_interface = r[index];
        else {
            auto rho_lo = dens[index];
            auto rho_hi = dens[index-1];

            auto x = ( (0.5 * max_dens) - rho_lo) / (rho_hi - rho_lo);

            r_interface = x * r[index-1] + (1.0-x) * r[index];
        }

		// output
        Print() << "time = " << data.Time() << ", r_interface = " << r_interface << std::endl;
	}


	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}

template <typename T>
Vector<size_t> sort_indexes(const Vector<T> &v) {

  // initialize original index locations
  Vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
