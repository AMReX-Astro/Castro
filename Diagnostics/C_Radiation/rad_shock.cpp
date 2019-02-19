
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
		int idir = 0;

		GetInputArgs (argc, argv, pltfile, slcfile, xctr, yctr, idir);

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
		Vector <int> varComps;

#if (AMREX_SPACEDIM == 1)
		Vector<std::string> compVarNames = {"density", "eint_E", "Temp",
			                            "pressure", "rad", "x_velocity"};
#elif (AMREX_SPACEDIM == 2)
		Vector<std::string> compVarNames = {"density",
			                            "eint_E", "Temp", "pressure", "rad",
			                            "x_velocity", "y_velocity"};
#else
		Vector<std::string> compVarNames = {"density", "eint_E", "Temp", "pressure",
			                            "rad", "x_velocity", "y_velocity","z_velocity"};
#endif
		GetComponents(data, compVarNames, varComps);

		auto cnt = 0;
		auto r1 = 1.0;

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++) {
			fill_comps[i] = i;
		}

		const BoxArray& ba_fine = data.boxArray(finestLevel);
		const DistributionMapping& dm_fine = data.DistributionMap(finestLevel);

		MultiFab data_mf(ba_fine, dm_fine, data.NComp(), data.NGrow());
		data.FillVar(data_mf, finestLevel, varNames, fill_comps);

		// allocate storage for data
		Vector<Real> vars_bin(nbins * (data.NComp()+1), 0.);

		// ! imask will be set to false if we've already output the data.
		// ! Note, imask is defined in terms of the finest level.  As we loop
		// ! over levels, we will compare to the finest level index space to
		// ! determine if we've already output here
		int mask_size = nbins;
		for (auto i = 0; i < finestLevel - 1; i++)
			mask_size *= rr[i];

#if (AMREX_SPACEDIM == 1)
		Vector<int> imask(mask_size);
#elif (AMREX_SPACEDIM >=2)
		Vector<int> imask(pow(mask_size, AMREX_SPACEDIM));
#endif

		for (auto it=imask.begin(); it!=imask.end(); ++it)
			*it = 1;

		cnt = 0;

		// extract the 1d data
		for (int l = finestLevel; l >= 0; l--) {

			int refratio = 1;
			for (auto lev = 0; lev < l; lev++) refratio *= rr[lev];

			Vector<Real> level_dx = data.DxLevel()[l];

			const BoxArray& ba = data.boxArray(l);
			const DistributionMapping& dm = data.DistributionMap(l);

			MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
			data.FillVar(lev_data_mf, l, varNames, fill_comps);

			for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
				const Box& bx = mfi.tilebox();

				fradshock(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				          ZFILL(problo), ZFILL(probhi),
				          BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				          nbins, vars_bin.dataPtr(),
				          imask.dataPtr(), mask_size, r1, refratio,
				          ZFILL(dx), idir, &cnt);

			}

			// adjust r1 for the next lowest level
			if (l != 0) r1 *= rr[l-1];
		}

#if (AMREX_SPACEDIM == 1)
		Vector<std::string> slcvarNames = {"density", "x-velocity",
			                           "int. energy", "temperature", "pressure",
			                           "rad energy", "rad temp"};
#elif (AMREX_SPACEDIM == 2)
		Vector<std::string> slcvarNames = {"density", "x-velocity", "y-velocity",
			                           "int. energy", "temperature", "pressure",
			                           "rad energy", "rad temp"};
#else
		Vector<std::string> slcvarNames = {"density", "x-velocity", "y-velocity",
			                           "z-velocity", "int. energy", "temperature", "pressure", "rad energy", "rad temp"};
#endif

		// reshape vars_bin into a vector of vectors
		// TODO: check that I've got the indexing correct here
		Vector<Vector<Real> > vars_array(data.NComp() + 1);
		for (auto i = 0; i < data.NComp() + 1; i++) {
			vars_array[i] = Vector<Real>(nbins);
			for (auto j = 0; j < nbins; j++) {
				vars_array[i][j] = vars_bin[i * nbins + j];
			}
		}

		WriteSlicefile(nbins, r, slcvarNames, vars_array, slcfile);

	}


	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
