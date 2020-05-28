/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

using namespace amrex;

void Castro::flame_width_properties(Real time, Real& T_max, Real& T_min,
                                    Real& grad_T_max) {
    BL_PROFILE("Castro::flame_width_properties()");

    const Real* dx = geom.CellSize();

    auto mf = derive("Temp", time, 1);

    BL_ASSERT(mf != nullptr);

#ifdef _OPENMP
#pragma omp parallel reduction(max : T_max, grad_T_max) reduction(min : T_min)
#endif
    for (MFIter mfi(*mf, true); mfi.isValid(); ++mfi) {
        FArrayBox& fab = (*mf)[mfi];

        const Box& box = mfi.tilebox();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();

        flame_width_temp(BL_TO_FORTRAN_ANYD(fab), ARLIM_3D(lo), ARLIM_3D(hi),
                         ZFILL(dx), &time, &T_max, &T_min, &grad_T_max);
    }
}

void Castro::flame_speed_properties(Real time, Real& rho_fuel_dot) {
    BL_PROFILE("Castro::flame_speed_properties()");

    const Real* dx = geom.CellSize();
    std::vector<std::string> spec_names;
    for (int i = 0; i < NumSpec; i++) {
        spec_names.push_back(short_spec_names_cxx[i]);
    }

    std::string name;

    for (auto nm : spec_names) {
        if (nm == "He4") {
            name = "omegadot_He4";
            break;
        }

        if (nm == "he4") {
            name = "omegadot_he4";
            break;
        }
    }

    auto mf = derive(name, time, 0);
    BL_ASSERT(mf != nullptr);
    Real rho_fuel_dot_temp = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+ : rho_fuel_dot_temp)
#endif
    for (MFIter mfi(*mf, true); mfi.isValid(); ++mfi) {
        FArrayBox& fab = (*mf)[mfi];

        const Box& box = mfi.tilebox();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();

        flame_speed_data(BL_TO_FORTRAN_ANYD(fab), ARLIM_3D(lo), ARLIM_3D(hi),
                         ZFILL(dx), &rho_fuel_dot_temp);
    }

    rho_fuel_dot += rho_fuel_dot_temp;
}
