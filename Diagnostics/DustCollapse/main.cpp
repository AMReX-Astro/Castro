//
// Process a group of n-d plotfiles from the dustcollapse problem
// and output the position of the interface as a function of time.
//
// For 2d:
//      The initial dense sphere is assumed to be centered a r = 0 (x = 0).
//      We take as default that it is centered vertically at y = 0, but
//      this can be overridden with --yctr.
//
// For 3d:
//      The initial dense sphere is assumed to be centered a x = y = z = 0,
//      but this can be overridden with --{x,y,z}ctr.
//
#include <iostream>
#include <regex>
#include <AMReX_DataServices.H>
#include <DustCollapse_F.H>

using namespace amrex;

std::string inputs_name = "";

template <typename T>
Vector<size_t> sort_indexes(const Vector<T> &v);

string GetVarFromJobInfo (const string pltfile, const string varname);

Vector<Real> GetCenter (const string pltfile);

int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv, false);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // Input arguments
    if (argc < 2)
        Abort("ERROR: Missing plotfiles");

    Print() << "\nUsing a density threshold of half the maximum analytic density" << std::endl;

    auto farg = 1;
    bool profile = false;

#if (AMREX_SPACEDIM >= 2)
    auto j = 1;
    while (j < argc) {
        if ( !strcmp(argv[j], "--profile")  )
        {
            profile = true;;
            farg++;
        }
        // Go to the next parameter name
        ++j;
    }

#endif

    // loop over plotfiles
    for (auto f = farg; f < argc; f++) {

        string pltfile = argv[f];

        auto center = GetCenter(pltfile);
        Real xctr = center[0];
        Real yctr = center[1];
        Real zctr = center[2];

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
        auto dx = data.CellSize(finestLevel);
        const Vector<Real>& problo = data.ProbLo();
        const Vector<Real>& probhi = data.ProbHi();

        Vector<int> rr = data.RefRatio();

        auto rmin = problo[0];

        const auto dx_fine = *(std::min_element(dx.begin(), dx.end()));

        // compute the maximum number of zones, as if we were completely refined
#if (AMREX_SPACEDIM == 1)
        int nbins = domain.length(0);
#elif (AMREX_SPACEDIM == 2)
        auto x_maxdist = std::max(std::abs(problo[0]), std::abs(probhi[0]));
        auto y_maxdist = std::max(std::abs(problo[1] - yctr), std::abs(probhi[1] - yctr));

        auto max_dist = std::sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist);

        int nbins = int(max_dist / dx_fine);
#else
        auto x_maxdist = std::max(std::abs(problo[0] - xctr), std::abs(probhi[0] - xctr));
        auto y_maxdist = std::max(std::abs(problo[1] - yctr), std::abs(probhi[1] - yctr));
        auto z_maxdist = std::max(std::abs(problo[2] - zctr), std::abs(probhi[2] - zctr));

        auto max_dist = std::sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist +
                             z_maxdist*z_maxdist);

        int nbins = int(max_dist / dx_fine);
#endif

        // radial coordinate
        Vector<Real> r(nbins);
        for (auto i = 0; i < nbins; i++)
            r[i] = (i + 0.5) * dx_fine + rmin;

        // find variable indices
        auto dens_comp = data.StateNumber("density");

        if (dens_comp < 0 )
            Abort("ERROR: density variable not found");

        // allocate storage for data
        Vector<Real> dens(nbins, 0.);
        Vector<Real> volcount(nbins, 0);

        // r1 is the factor between the current level grid spacing and the
        // FINEST level
        auto r1 = 1.0;

        Vector<int> fill_comps(data.NComp());
        for (auto i = 0; i < data.NComp(); i++)
            fill_comps[i] = i;

        // imask will be set to false if we've already output the data.
        // Note, imask is defined in terms of the finest level.  As we loop
        // over levels, we will compare to the finest level index space to
        // determine if we've already output here
        int mask_size = domain.length().max();
        Vector<int> imask(std::pow(mask_size, AMREX_SPACEDIM), 1);

        // counter
        int cnt = 0;

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

#if (AMREX_SPACEDIM == 1)
                fdustcollapse1d(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
                                BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
                                nbins, dens.dataPtr(),
                                imask.dataPtr(), mask_size, r1, dens_comp, &cnt);
#elif (AMREX_SPACEDIM == 2)

                fdustcollapse2d(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
                                BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
                                nbins, dens.dataPtr(), volcount.dataPtr(),
                                imask.dataPtr(), mask_size, r1,
                                AMREX_ZFILL(level_dx), dx_fine, yctr, dens_comp);
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
        // normalize
        for (int i = 0; i < nbins; i++)
            if (volcount[i] != 0.)
                dens[i] /= volcount[i];
#endif

        // These are calculated analytically given initial density 1.e9 and the
        // analytic expression for the radius as a function of time t = 0.00
        Real max_dens = 1.e9;

        if (std::abs(data.Time()) <= 1.e-8)
            max_dens = 1.e9;
        else if (std::abs(data.Time() - 0.01) <= 1.e-8)
            max_dens = 1.043345e9;
        else if (std::abs(data.Time() - 0.02) <= 1.e-8)
            max_dens = 1.192524e9;
        else if (std::abs(data.Time() - 0.03) <= 1.e-8)
            max_dens = 1.527201e9;
        else if (std::abs(data.Time() - 0.04) <= 1.e-8)
            max_dens = 2.312884e9;
        else if (std::abs(data.Time() - 0.05) <= 1.e-8)
            max_dens = 4.779133e9;
        else if (std::abs(data.Time() - 0.06) <= 1.e-8)
            max_dens = 24.472425e9;
        else if (std::abs(data.Time() - 0.065) <= 1.e-8)
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

        // dump out profile file
        if (profile) {
            string outfile_name = pltfile;

            if (pltfile.back() == '/')
                outfile_name += "prof.profile";
            else
                outfile_name += "/prof.profile";

            std::ofstream outfile;
            outfile.open(outfile_name);
            outfile.setf(std::ios::scientific);
            outfile.precision(12);
            const auto w = 24;

            // write the header
            outfile << std::setw(w) << "r" << std::setw(w) << "density" << std::endl;

            for (auto i = 0; i < nbins; i++)
                outfile << std::setw(w) << r[i] << std::setw(w) << dens[i] << std::endl;

            outfile.close();

        }
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

///
/// Gets the variable ``varname`` from the ``job_info`` file and returns as a
/// string
///
string GetVarFromJobInfo (const string pltfile, const string varname) {
    string filename = pltfile + "/job_info";
    std::regex re("(?:[ \\t]*)" + varname + "\\s*:\\s*(.*)\\s*\\n");

    std::smatch m;

    std::ifstream jobfile(filename);
    if (jobfile.is_open()) {
        std::stringstream buf;
        buf << jobfile.rdbuf();
        string file_contents = buf.str();

        if (std::regex_search(file_contents, m, re)) {
            return m[1];
        } else {
            Print() << "Unable to find " << varname << " in job_info file!" << std::endl;
        }
    } else {
        Print() << "Could not open job_info file!" << std::endl;
    }

    return "";
}

// Get the center from the job info file and return as a Real Vector
Vector<Real> GetCenter (const string pltfile) {
    auto center_str = GetVarFromJobInfo(pltfile, "center");

    // split string
    std::istringstream iss {center_str};
    Vector<Real> center;

    std::string s;
    while (std::getline(iss, s, ','))
        center.push_back(stod(s));

    return center;
}
