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

    Real r{3.89e7};

    std::cout << "testing locate" << std::endl;

    auto idx = locate(r, 0);
    AMREX_ALWAYS_ASSERT(r >= model::profile(0).r(idx) &&
                        r <= model::profile(0).r(idx+1));

    std::cout << "testing interpolate" << std::endl;

    // density is monotonically decreasing
    auto dens = interpolate(r, model::idens);
    AMREX_ALWAYS_ASSERT(dens <= model::profile(0).state(idx, model::idens) &&
                        dens >= model::profile(0).state(idx+1, model::idens));

}
