
#include <iostream>
#include "AMReX_DataServices.H"
#include <Radiation_F.H>
#include <Radiation.H>

using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv, false);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	// wallclock time
	const Real strt_total = amrex::second();

	{
		// Input arguments
		string pltfile;
		string slcfile;
		double xctr = 0.0;
		double yctr = 0.0;
		int dir;

		GetInputArgs (argc, argv, pltfile, slcfile, xctr, yctr, dir);

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

		Vector<Real> dx(AMREX_SPACEDIM);
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			dx[i] = data.ProbSize()[i] / domain.length(i);

		const Vector<Real>& problo = data.ProbLo();
		const Vector<Real>& probhi = data.ProbHi();
		Vector<int> rr = data.RefRatio();

		// compute the size of the radially-binned array -- we'll do it to
		// the furtherest corner of the domain
		double x_maxdist = max(fabs(probhi[0] - xctr), fabs(problo[0] - xctr));
		double y_maxdist = max(fabs(probhi[1] - yctr), fabs(problo[1] - yctr));
		double maxdist = sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist);

		double dx_fine = *(std::min_element(dx.begin(), dx.end()));

		int nbins = int(maxdist / dx_fine);

		Vector<Real> r(nbins);

		for (auto i = 0; i < nbins; i++)
			r[i] = (i + 0.5) * dx_fine;

		// find variable indices
		Vector<std::string> slcvarNames = {"density"};
		const auto nvars = slcvarNames.size();

		auto varComps = GetComponents(data, slcvarNames);

		auto rad_comp = varComps[0];

		// allocate storage for data
		Vector<Real> rad_bin(nbins, 0.);
		Vector<int> ncount(nbins, 0);

		auto r1 = 1.0;

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++) {
			fill_comps[i] = i;
		}

		// imask will be set to false if we've already output the data.
		// Note, imask is defined in terms of the finest level.  As we loop
		// over levels, we will compare to the finest level index space to
		// determine if we've already output here
		int mask_size = domain.length().max();

#if (AMREX_SPACEDIM == 1)
             Vector<int> imask(mask_size);
#elif (AMREX_SPACEDIM >=2)
             Vector<int> imask(pow(mask_size, AMREX_SPACEDIM));
#endif

		for (auto it=imask.begin(); it!=imask.end(); ++it)
			*it = 1;

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

				fgaussian_pulse(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				                BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				                nbins, rad_bin.dataPtr(), ncount.dataPtr(),
				                imask.dataPtr(), mask_size, r1,
				                rad_comp, ZFILL(dx), dx_fine, xctr, yctr);

			}

			// adjust r1 for the next lowest level
			if (l != 0) r1 *= rr[l-1];
		}

		//normalize
		for (int i = 0; i < nbins; i++) {
			if (ncount[i] != 0)
				rad_bin[i] /= ncount[i];
		}


		Vector<Vector<Real> > vars(nvars);
		vars[0] = rad_bin;

		WriteSlicefile(nbins, r, slcvarNames, vars, slcfile);

	}


	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
