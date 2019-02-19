//
// extract a 1-d slice of the data (all variables or a single variable)
// along the specified coordinate direction from a plotfile.  The
// plotfile can be 1-, 2-, or 3-d.
//
// This routine is a generalized version is based on fextract3d, but geared
// toward the CASTRO radiating shock problem
//
// We read in all the variables, but only output a subset
//
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
		int idir = 1;

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

		int nbins = probhi[AMREX_SPACEDIM-1] - problo[AMREX_SPACEDIM-1] + 1;

		double dx_fine = *(std::min_element(dx.begin(), dx.end()));

		Vector<Real> r(nbins);

		for (auto i = 0; i < nbins; i++)
			r[i] = (i + 0.5) * dx[idir-1] + problo[idir-1];

		// find variable indices
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
		auto varComps = GetComponents(data, compVarNames);

		auto cnt = 0;
		auto r1 = 1.0;

		const int nvars = data.NComp();

		// fill a multifab with the data
		Vector<int> fill_comps(nvars);
		for (auto i = 0; i < nvars; i++) {
			fill_comps[i] = i;
		}

		// allocate storage for data
		Vector<Real> vars_bin(nbins * (nvars + 1), 0.);

		// imask will be set to false if we've already output the data.
		// Note, imask is defined in terms of the finest level.  As we loop
		// over levels, we will compare to the finest level index space to
		// determine if we've already output here
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

			MultiFab lev_data_mf(ba, dm, nvars, data.NGrow());
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

		// sort the data based on the coordinates

		Vector<Real> coords(cnt);
		for (auto i = 0; i < cnt; i++)
			coords[i] = vars_bin[i];

		auto isv = sort_indexes(coords);

		// find variable indices
#if (AMREX_SPACEDIM == 1)
		compVarNames = {"density", "x_velocity", "pressure",
			        "eint_E", "Temp", "rad", "rad"};
#elif (AMREX_SPACEDIM == 2)
		compVarNames = {"density", "x_velocity", "y_velocity", "pressure",
			        "eint_E", "Temp", "rad", "rad"};
#else
		compVarNames = {"density", "x_velocity", "y_velocity", "z_velocity",
			        "pressure","eint_E", "Temp", "rad", "rad"};
#endif
		Vector<int> comps = GetComponents(data, compVarNames);

#if (AMREX_SPACEDIM == 1)
		Vector<std::string> slcvarNames = {"density", "x-velocity",
			                           "pressure", "int. energy", "temperature",
			                           "rad energy", "rad temp"};
#elif (AMREX_SPACEDIM == 2)
		Vector<std::string> slcvarNames = {"density", "x-velocity", "y-velocity",
			                           "pressure", "int. energy", "temperature",
			                           "rad energy", "rad temp"};
#else
		Vector<std::string> slcvarNames = {"density", "x-velocity", "y-velocity",
			                           "z-velocity", "pressure", "int. energy", "temperature", "rad energy", "rad temp"};
#endif

		// write to file
		std::ofstream slicefile;
		slicefile.open(slcfile);
		slicefile.setf(std::ios::scientific);
		slicefile.precision(12);
		const auto w = 24;

		// write the header
		slicefile << std::setw(w) << "x";
		for (auto it=slcvarNames.begin(); it!=slcvarNames.end(); ++it)
			slicefile << std::setw(w) << *it;

		slicefile << std::endl;

		const Real arad = 1.0;

		// write the data in columns
		for (auto i = 0; i < cnt; i++) {

			slicefile << std::setw(w) << vars_bin[isv[i]];

			for (auto it=comps.begin(); it!=(comps.end()-1); ++it)
				slicefile << std::setw(w) << vars_bin[isv[i]+((*it)+1)*nbins];

			auto it = comps.end()-1;
			slicefile << std::setw(w) << pow(vars_bin[isv[i]+((*it)+1)*nbins] / arad, 0.25);

			slicefile << std::endl;

		}

		slicefile.close();

	}


	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
