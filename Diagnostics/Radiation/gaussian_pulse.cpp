//
// Process a 2-d gaussian radiation pulse
//
#include <iostream>
#include <AMReX_DataServices.H>
#include <Radiation_F.H>
#include <Radiation_utils.H>

using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv, false);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // abort if dim != 2
    if (AMREX_SPACEDIM != 2)
        Abort("ERROR: must compile gaussian pulse diagnostic with DIM=2");

    // Input arguments
    string pltfile, slcfile;
    int dir;

    GetInputArgs(argc, argv, pltfile, slcfile dir);

    auto center = GetCenter(pltfile);
    double xctr = center[0];
    double yctr = center[1];

    Print() << "xctr = " << xctr << std::endl;
    Print() << "yctr = " << yctr << std::endl;
    Print() << std::endl;

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
    double x_maxdist = std::max(std::abs(probhi[0] - xctr), std::abs(problo[0] - xctr));
    double y_maxdist = std::max(std::abs(probhi[1] - yctr), std::abs(problo[1] - yctr));
    double maxdist = std::sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist);

    double dx_fine = *(std::min_element(dx.begin(), dx.end()));

    int nbins = int(maxdist / dx_fine);

    // radial coordinate
    Vector<Real> r(nbins);
    for (auto i = 0; i < nbins; i++)
        r[i] = (i + 0.5) * dx_fine;

    // find variable indices
    Vector<std::string> slcvarNames = {"rad"};
    const auto nvars = slcvarNames.size();

    auto varComps = GetComponents(data, slcvarNames);
    auto rad_comp = varComps[0];

    // allocate storage for data
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
    Vector<int> imask(std::pow(mask_size, AMREX_SPACEDIM), 1);

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

            fgaussian_pulse(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
                            BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
                            nbins, rad_bin.dataPtr(), ncount.dataPtr(),
                            imask.dataPtr(), mask_size, r1,
                            rad_comp, AMREX_ZFILL(dx), dx_fine, xctr, yctr);

        }

        // adjust r1 for the next lowest level
        if (l != 0) r1 *= rr[l-1];
    }

    // normalize
    for (int i = 0; i < nbins; i++)
        if (ncount[i] != 0)
            rad_bin[i] /= ncount[i];

    // write data to slicefile
    Vector<Vector<Real> > vars(nvars);
    vars[0] = rad_bin;

    WriteSlicefile(nbins, r, slcvarNames, vars, slcfile);

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
