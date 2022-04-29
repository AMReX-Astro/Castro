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

    // the executable name is the first arg

    int farg{1};

    if (narg == 1) {
        amrex::Print()
            << "\n"
            << " Compute the average and L2 norm of the burn weights on the finest grid.\n"
            << " Usage:\n"
            << "    fburn_weight_avg plotfile(s)\n"
            << "\n"
            << std::endl;
        return;
    }

    int ntime = narg;

    for (int f = 0; f < ntime; ++f) {
        const std::string& filename = amrex::get_command_argument(f+farg);

        PlotFileData pf(filename);

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
        Vector<Real> lmax(var_names.size(), -1.e30);
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
                                    lmax[ivar] = amrex::max(lmax[ivar], fab(i,j,k));
                                }
                            }
                        }
                    }
                }
            }
        }

        ParallelDescriptor::ReduceRealSum(lsum.data(), lsum.size());
        ParallelDescriptor::ReduceRealSum(l2sum.data(), l2sum.size());
        ParallelDescriptor::ReduceRealMax(lmax.data(), lmax.size());
        ParallelDescriptor::ReduceRealSum(nzones);

        if (ParallelDescriptor::IOProcessor()) {

            // average

            for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                if (nzones > 0) {
                    lsum[ivar] /= nzones;
                    l2sum[ivar] = std::sqrt(l2sum[ivar] / nzones);
                }
                if (lmax[ivar] < 0) {
                    lmax[ivar] = 0;
                }
            }

            if (f == 0) {
                // output the header
                std::cout << "# " << std::setw(23) << std::left << "time";
                for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                    std::cout << std::setw(25) << std::left << var_names[ivar] + "-avg" << " ";
                    std::cout << std::setw(25) << std::left << var_names[ivar] + "-L2" << " ";
                    std::cout << std::setw(25) << std::left << var_names[ivar] + "-max" << " ";
                }
                std::cout << std::endl;
            }

            std::cout << std::setw(25) << pf.time();
            for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                std::cout << std::setw(25) << std::right << lsum[ivar] << " ";
                std::cout << std::setw(25) << std::right << l2sum[ivar] << " ";
                std::cout << std::setw(25) << std::right << lmax[ivar] << " ";
            }
            std::cout << std::endl;

        }
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
