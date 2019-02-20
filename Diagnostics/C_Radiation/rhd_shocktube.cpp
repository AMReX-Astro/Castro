
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

	{
		// Input arguments
		string pltfile, groupfile, slcfile;
		int j = 1; // skip program name

		while ( j < argc)
		{

			if ( !strcmp(argv[j], "-p") || !strcmp(argv[j],"--pltfile") )
			{
				pltfile = argv[++j];
			}
			else if ( !strcmp(argv[j], "-g") || !strcmp(argv[j],"--groupfile") )
			{
				groupfile = argv[++j];
			}
			else if ( !strcmp(argv[j], "-s") || !strcmp(argv[j],"--slicefile") )
			{
				slcfile = argv[++j];
			}
			else
			{
				std::cout << "\n\nOption " << argv[j] << " not recognized" << std::endl;
				PrintHelp();
				exit ( EXIT_FAILURE );
			}

			// Go to the next parameter name
			++j;
		}

		if (pltfile.empty() || groupfile.empty())
		{
			PrintHelp();
			Abort("Missing input file");
		}

		Print() << "\nplotfile  = \"" << pltfile << "\"" << std::endl;
		Print() <<   "groupfile = \"" << groupfile << "\"" << std::endl;
		Print() <<   "slicefile = \"" << slcfile << "\"" << std::endl;
		Print() << std::endl;

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
		Box domain = data.ProbDomain()[0];
		Vector<Real> dx = data.CellSize(0);

		// Vector<Real> dx(AMREX_SPACEDIM);
		// for (int i = 0; i < AMREX_SPACEDIM; i++)
		//      dx[i] = data.ProbSize()[i] / domain.length(i);

		const Vector<Real>& problo = data.ProbLo();
		const Vector<Real>& probhi = data.ProbHi();
		Vector<int> rr = data.RefRatio();

		// compute the size of the radially-binned array -- we'll do it to
		// the furtherest corner of the domain
		// double maxdist = fabs(probhi[0] - problo[0]);

		double dx_coarse = *(std::min_element(dx.begin(), dx.end()));

		// int nbins = int(maxdist / dx_coarse);

		int nbins = domain.bigEnd()[0] - domain.smallEnd()[0] + 1;

		Vector<Real> x(nbins);

		for (auto i = 0; i < nbins; i++)
			x[i] = (i + 0.5) * dx_coarse;

		// find variable indices
		Vector<std::string> compVarNames = {"density", "x_velocity", "pressure", "rad0"};
		auto varComps = GetComponents(data, compVarNames);

		auto dens_comp = varComps[0];
		auto velx_comp = varComps[1];
		auto pres_comp = varComps[2];
		auto rad0_comp = varComps[3];

		// allocate storage for data
		Vector<Real> dens_bin(nbins, 0.);
		Vector<Real> vel_bin(nbins, 0.);
		Vector<Real> pres_bin(nbins, 0.);
		Vector<Real> rad_bin(nbins, 0.);

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++)
			fill_comps[i] = i;

		const BoxArray& ba_fine = data.boxArray(finestLevel);
		const DistributionMapping& dm_fine = data.DistributionMap(finestLevel);

		MultiFab data_mf(ba_fine, dm_fine, data.NComp(), data.NGrow());
		data.FillVar(data_mf, finestLevel, varNames, fill_comps);

		// extract the 1d data
		auto l = 0;

		Vector<Real> level_dx = data.DxLevel()[l];

		const BoxArray& ba = data.boxArray(l);
		const DistributionMapping& dm = data.DistributionMap(l);

		MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
		data.FillVar(lev_data_mf, l, varNames, fill_comps);

		for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
			const Box& bx = mfi.tilebox();

			frhdshocktube(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			              BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
			              nbins, dens_bin.dataPtr(), vel_bin.dataPtr(),
			              pres_bin.dataPtr(), rad_bin.dataPtr(),
			              dens_comp, velx_comp, pres_comp, rad0_comp);

		}

		Vector<Vector<Real> > vars = {dens_bin, vel_bin, pres_bin, rad_bin};
		Vector<std::string> slcvarNames = {"density", "velocity", "pressure", "rad"};
		WriteSlicefile(nbins, x, slcvarNames, vars, slcfile);

	}


	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
