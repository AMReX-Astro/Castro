#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>
#include <cmath>

// compute the integral of f dV, where f is one of the fields in the plotfile
// this understands axisymmetric geometry.

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile;
    std::string varname_arg;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-v" || name == "--variable") {
            varname_arg = amrex::get_command_argument(++farg);
        } else {
            break;
        }
        ++farg;
    }

    if (pltfile.empty() && farg <= narg) {
        pltfile = amrex::get_command_argument(farg);
    }

    if (pltfile.empty()) {
        amrex::Print()
            << "\n"
            << " Compute the integral of f dV, where f is the field specified via -v.\n"
            << " Works with 1-, 2-, or 3-d datasets, including 2-d axisymmetric.\n"
            << "\n"
            << " Usage:\n"
            << "    fvolumesum [-v variable] plotfile\n"
            << "\n"
            << " args [-v|--variable]  varname      : the field to integrate over\n"
            << "\n"
            << std::endl;
        return;
    }


    PlotFileData pf(pltfile);
    const Vector<std::string>& var_names_pf = pf.varNames();

    // find all of the variables that begin with burn_weights

    Vector<std::string> var_names;

    for (auto e : var_names_pf) {
        if (e.find("burn_weights") != std::string::npos) {
            var_names.push_back(e);
        }
    }

    int fine_level = pf.finestLevel();

    Vector<Real> lsum(var_names.size(), 0.0);
    Vector<Real> l2sum(var_names.size(), 0.0);
    Real nzones{0};

    int level = fine_level;

    for (int ivar = 0; ivar < var_names.size(); ++ivar) {
        const MultiFab& mf = pf.get(level, var_names[ivar]);
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            if (bx.ok()) {
                const auto& fab = mf.array(mfi);
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                for (int k = lo.z; k <= hi.z; ++k) {
                    for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            if (fab(i,j,k) > 0.0) {
                                if (ivar == 0) {
                                    nzones += 1;
                                }
                                lsum[ivar] += fab(i,j,k);
                                l2sum[ivar] += fab(i,j,k) * fab(i,j,k);                                
                            }
                        }
                    }
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(lsum.data(), lsum.size());
    ParallelDescriptor::ReduceRealSum(l2sum.data(), l2sum.size());
    ParallelDescriptor::ReduceRealSum(nzones);

    if (ParallelDescriptor::IOProcessor()) {

        // average

        for (int ivar = 0; ivar < var_names.size(); ++ivar) {
            lsum[ivar] /= nzones;
            l2sum[ivar] = std::sqrt(l2sum[ivar] / nzones);
        }

        for (int ivar = 0; ivar < var_names.size(); ++ivar) {
            std::cout << std::setw(25) << std::left << var_names[ivar] + "-avg" << " ";
            std::cout << std::setw(25) << std::left << var_names[ivar] + "-L2" << " ";
        }
        std::cout << std::endl;

        for (int ivar = 0; ivar < var_names.size(); ++ivar) {
            std::cout << std::setw(25) << std::right << lsum[ivar] << " ";
            std::cout << std::setw(25) << std::right << l2sum[ivar] << " ";
        }
        std::cout << std::endl;

    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
