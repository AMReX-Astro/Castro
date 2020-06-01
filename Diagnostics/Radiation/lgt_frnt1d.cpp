//
// Process a 1-d sedov problem to produce rho, u, and p as a
// function of r, for comparison to the analytic solution.
//
#include <iostream>
#include "AMReX_DataServices.H"
#include <Radiation_F.H>
#include <Radiation_utils.H>

using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv, false);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	// Input arguments
	string pltfile, slcfile;
	int dir;

	GetInputArgs(argc, argv, pltfile, slcfile, dir);

	// Start dataservices
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
	Vector<Real> dx = data.CellSize(finestLevel);
	const Vector<Real>& problo = data.ProbLo();
	const Vector<Real>& probhi = data.ProbHi();

	Vector<int> rr = data.RefRatio();

	// compute the size of the radially-binned array -- we'll do it to
	// the furtherest corner of the domain
	double maxdist = fabs(probhi[0] - problo[0]);
	double dx_fine = *(std::min_element(dx.begin(), dx.end()));
	int nbins = int(maxdist / dx_fine);

	// radial coordinate
	Vector<Real> r(nbins);

	for (auto i = 0; i < nbins; i++)
		r[i] = (i + 0.5) * dx_fine;

	// find variable indices
	Vector<std::string> compVarNames = {"density", "xmom", "pressure", "rad"};

	auto varComps = GetComponents(data, compVarNames);
	auto dens_comp = varComps[0];
	auto xmom_comp = varComps[1];
	auto pres_comp = varComps[2];
	auto rad_comp = varComps[3];

	// allocate storage for data
	Vector<Real> dens_bin(nbins, 0.);
	Vector<Real> vel_bin(nbins, 0.);
	Vector<Real> pres_bin(nbins, 0.);
	Vector<Real> rad_bin(nbins, 0.);
	Vector<int> ncount(nbins, 0);

	// r1 is the factor between the current level grid spacing and the
	// FINEST level
	auto r1 = 1.0;

	// fill a multifab with the data
	Vector<int> fill_comps(data.NComp());
	for (auto i = 0; i < data.NComp(); i++)
		fill_comps[i] = i;

	// imask will be set to false if we've already output the data.
	// Note, imask is defined in terms of the finest level.  As we loop
	// over levels, we will compare to the finest level index space to
	// determine if we've already output here
	int mask_size = domain.length().max();
	Vector<int> imask(pow(mask_size, AMREX_SPACEDIM), 1);

	// loop over the data, starting at the finest grid, and if we haven't
	// already stored data in that grid location (according to imask),
	// store it.
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

			flgt_frnt1d(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			            BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
			            nbins, dens_bin.dataPtr(), vel_bin.dataPtr(),
			            pres_bin.dataPtr(), rad_bin.dataPtr(),
			            imask.dataPtr(), mask_size, r1,
			            dens_comp, xmom_comp, pres_comp, rad_comp,
			            ZFILL(dx), dx_fine);
		}

		// adjust r1 for the next lowest level
		if (l != 0) r1 *= rr[l-1];
	}

	//normalize
	for (int i = 0; i < nbins; i++)
		if (ncount[i] != 0)
			rad_bin[i] /= ncount[i];

	// write data to slicefile
	Vector<Vector<Real> > vars = {dens_bin, vel_bin, pres_bin, rad_bin};
	Vector<std::string> slcvarNames = {"density", "velocity", "pressure", "rad"};

	WriteSlicefile(nbins, r, slcvarNames, vars, slcfile);

	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}
