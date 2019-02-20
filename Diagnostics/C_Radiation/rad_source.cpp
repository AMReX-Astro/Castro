
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

	// Input arguments
	if (argc < 2)
		Abort("ERROR: Missing plotfiles");

	// all the arguments are assumed to be plotfiles
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

		// get variable names
		const Vector<string>& varNames = data.PlotVarNames();

		// find variable indices
		Vector<std::string> compVarNames = {"rho_e", "rad"};

		auto varComps = GetComponents(data, compVarNames);
		auto rhoe_comp = varComps[0];
		auto rad_comp = varComps[1];

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto j = 0; j < data.NComp(); j++)
			fill_comps[j] = j;

		auto lev = 0;

		const BoxArray& ba = data.boxArray(lev);
		const DistributionMapping& dm = data.DistributionMap(lev);

		MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
		data.FillVar(lev_data_mf, lev, varNames, fill_comps);

		// we only care about a single zone, so take the first box
		MFIter mfi(lev_data_mf);
		const Box& bx = mfi.tilebox();

		Real rhoe, rad;

		fradsource(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		           BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
		           &rhoe, &rad, rhoe_comp, rad_comp);

		const auto w = 20;

		std::cout.setf(std::ios::scientific);
		std::cout.precision(12);

		if (i == 1)
			std::cout << std::setw(w) << "time" << std::setw(w) << "rho e"
			        << std::setw(w) << "rad" << std::endl;

		std::cout << std::setw(w) << data.Time() << std::setw(w) << rhoe << std::setw(w)
		        << rad << std::endl;

	}

	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
