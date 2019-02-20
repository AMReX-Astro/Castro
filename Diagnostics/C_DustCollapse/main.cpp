
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

	Print() << "\nUsing a density threshhold of half the maximum analytic density" << std::endl;

	auto farg = 1;
	Real xctr, yctr, zctr = 0.0;

#if (AMREX_SPACEDIM >= 2)
	auto j = 1;
	while (j < argc) {
		if ( !strcmp(argv[j], "--xctr")  )
		{
			xctr = atof(argv[++j]);
			farg += 2;
		}
		else if ( !strcmp(argv[j], "--yctr")  )
		{
			yctr = atof(argv[++j]);
			farg += 2;
		}
		else if ( !strcmp(argv[j], "--zctr")  )
		{
			zctr = atof(argv[++j]);
			farg += 2;
		}
		// Go to the next parameter name
		++j;
	}

#endif

	for (auto f = farg; f < argc; f++) {

		string pltfile = argv[f];

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

		const auto dx_fine = *(std::min_element(dx.begin(), dx.end()));

#if (AMREX_SPACEDIM == 1)
		int nbins = domain.length(0);
#elif (AMREX_SPACEDIM == 2)
		// compute the size of the radially-binned array -- we'll do it to
		// the furtherest corner of the domain
		auto x_maxdist = fmax(fabs(problo[0]), fabs(probhi[0]));
		auto y_maxdist = fmax(fabs(problo[1] - yctr), fabs(probhi[1] - yctr));

		auto max_dist = sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist);

		int nbins = int(max_dist / dx_fine);
#else
		// compute the size of the radially-binned array -- we'll do it to
		// the furtherest corner of the domain
		auto x_maxdist = fmax(fabs(problo[0] - xctr), fabs(probhi[0] - xctr));
		auto y_maxdist = fmax(fabs(problo[1] - yctr), fabs(probhi[1] - yctr));
		auto z_maxdist = fmax(fabs(problo[2] - zctr), fabs(probhi[2] - zctr));

		auto max_dist = sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist +
		                     z_maxdist*z_maxdist);

		int nbins = int(max_dist / dx_fine);
#endif

		Vector<Real> r(nbins);

		for (auto i = 0; i < nbins; i++)
			r[i] = (i + 0.5) * dx_fine;

		// find variable indices
		auto dens_comp = data.StateNumber("density");

		if (dens_comp < 0 )
			Abort("ERROR: density variable not found");

		// allocate storage for data
		Vector<Real> dens(nbins, 0.);
		Vector<Real> volcount(nbins, 0);

		auto r1 = 1.0;

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++)
			fill_comps[i] = i;

		// imask will be set to false if we've already output the data.
		// Note, imask is defined in terms of the finest level.  As we loop
		// over levels, we will compare to the finest level index space to
		// determine if we've already output here
		int mask_size = 1;
		for (auto i = 0; i < AMREX_SPACEDIM; i++)
			mask_size = max(mask_size, domain.length(i));

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
				                nbins, dens.dataPtr(), volcount.dataPtr(),
				                imask.dataPtr(), mask_size, r1,
				                ZFILL(level_dx), dx_fine, yctr, dens_comp);
#else
				fdustcollapse3d(bx.loVect(), bx.hiVect(),
				                BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				                nbins, dens.dataPtr(), volcount.dataPtr(),
				                imask.dataPtr(), mask_size, r1,
				                level_dx.dataPtr(), dx_fine, xctr, yctr, zctr, dens_comp);
#endif
			}

			// adjust r1 for the next lowest level
			if (l != 0) r1 *= rr[l-1];
		}

#if (AMREX_SPACEDIM == 1)
		// sort the data based on the coordinates
		auto isv = sort_indexes(r);
#else
		//normalize
		for (int i = 0; i < nbins; i++)
			if (volcount[i] != 0.)
				dens[i] /= volcount[i];
#endif

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
#if (AMREX_SPACEDIM == 1)
		for (auto i = 0; i < cnt; i++) {
			if (dens[isv[i]] < 0.5 * max_dens) {
				index = i;
				break;
			}
#else
		for (auto i = 0; i < nbins; i++) {
			if (dens[i] < 0.5 * max_dens) {
				index = i;
				break;
			}
#endif
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
		Print() << "\ntime = " << data.Time() << ", r_interface = "
		        << std::setprecision(16) << r_interface << std::endl << std::endl;
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
	     [&v](size_t i1, size_t i2) {
		return v[i1] < v[i2];
	});

	return idx;
}
