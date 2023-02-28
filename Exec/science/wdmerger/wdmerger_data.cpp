#include <wdmerger_data.H>

using namespace amrex;

AMREX_GPU_MANAGED Real wdmerger::mdot_P = 0.0;
AMREX_GPU_MANAGED Real wdmerger::mdot_S = 0.0;

AMREX_GPU_MANAGED Real wdmerger::rad_P[7] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::rad_S[7] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::vol_P[7] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::vol_S[7] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::rho_avg_P = 0.0;
AMREX_GPU_MANAGED Real wdmerger::rho_avg_S = 0.0;
