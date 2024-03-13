#include <bitset>

#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>

using namespace amrex;

void ca_derxash(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                const FArrayBox& datfab, const Geometry& /*geomdata*/,
                Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // determine which species should be considered ash

    std::bitset<NumSpec> is_ash{};

    for (int i = 0; i < NumSpec; ++i) {
        // include all elements beyond oxygen
        if (zion[i] > 8.0) {
            is_ash.set(i);
        }
    }

    // exclude all of the "ash" species from the input file; they're actually
    // used for the star composition and hiding them helps make the flame more
    // visible
    for (const std::string& ash_name : {problem::ash1_name,
                                        problem::ash2_name,
                                        problem::ash3_name}) {
        int i = network_spec_index(ash_name);
        if (i != -1) {
            is_ash.reset(i);
        }
    }

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real sum = 0.0_rt;
        for (int n = 0; n < NumSpec; ++n) {
            if (is_ash[n]) {
                sum += dat(i,j,k,1+n)/dat(i,j,k,0);
            }
        }
        der(i,j,k,0) = sum;
    });
}

