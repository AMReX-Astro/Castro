#include <iostream>

#include <model_parser.H>


int main(int argc, char *argv[]) {

    amrex::Initialize(argc, argv);

    // initialize the external runtime parameters in C++

    init_extern_parameters();

    // now initialize the C++ Microphysics
#ifdef REACTIONS
    network_init();
#endif

    Real small_temp = 1.e-200;
    Real small_dens = 1.e-200;
    eos_init(small_temp, small_dens);

    std::string model = "sub_chandra.M_WD-1.10.M_He-0.050.hse.CO.N14.N.10.00km";
    read_model_file(model);

    // test locate and ensure that the index we get is such that
    // profile.r(idx) < r < profile.r(idx+1)

    Real r{3.89e7};

    std::cout << "testing locate" << std::endl;

    auto idx = locate(r, 0);
    AMREX_ALWAYS_ASSERT(r >= model::profile(0).r(idx) &&
                        r <= model::profile(0).r(idx+1));

    std::cout << "testing interpolate" << std::endl;

    // density is monotonically decreasing
    // test to make sure that we see that

    auto dens = interpolate(r, model::idens);
    AMREX_ALWAYS_ASSERT(dens <= model::profile(0).state(idx, model::idens) &&
                        dens >= model::profile(0).state(idx+1, model::idens));

    // test the bounds.  Our model spans r = [5.e5, 4.0955e9]

    std::cout << "testing boundaries of the model" << std::endl;

    auto idx_lo_bnd = locate(-1.0, 0);
    AMREX_ALWAYS_ASSERT(idx_lo_bnd == 0);
    std::cout << "value a r < 0 = " << interpolate(-1.0, model::idens) << std::endl;

    auto idx_hi_bnd = locate(4.1e9, 0);
    AMREX_ALWAYS_ASSERT(idx_hi_bnd == model::npts-2);
    std::cout << "value a r > r_max = " << interpolate(4.1e9, model::idens) << std::endl;

    // test if we interpolate to exactly a point in the model that we
    // recover the data there to roundoff

    std::cout << "testing interpolating exactly on model point" << std::endl;

    int idx_test{100};
    Real r_test = model::profile(0).r(idx_test);
    auto dens_test = interpolate(r_test, model::idens);

    AMREX_ALWAYS_ASSERT(std::abs(dens_test - model::profile(0).state(idx_test, model::idens)) < 1.e-15_rt);

}
