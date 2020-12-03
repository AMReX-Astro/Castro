#include <wdmerger_data.H>

using namespace amrex;

AMREX_GPU_MANAGED int wdmerger::use_stopping_criterion = 1;
AMREX_GPU_MANAGED int wdmerger::use_energy_stopping_criterion = 0;
AMREX_GPU_MANAGED Real wdmerger::ts_te_stopping_criterion = 1.e200;
AMREX_GPU_MANAGED Real wdmerger::T_stopping_criterion = 1.e200;

AMREX_GPU_MANAGED Real wdmerger::mdot_P = 0.0;
AMREX_GPU_MANAGED Real wdmerger::mdot_S = 0.0;

AMREX_GPU_MANAGED Real wdmerger::rad_P[7] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::rad_S[7] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::vol_P[7] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::vol_S[7] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::rho_avg_P = 0.0;
AMREX_GPU_MANAGED Real wdmerger::rho_avg_S = 0.0;

AMREX_GPU_MANAGED Real wdmerger::T_curr_max = 0.0;
AMREX_GPU_MANAGED Real wdmerger::rho_curr_max = 0.0;
AMREX_GPU_MANAGED Real wdmerger::ts_te_curr_max = 0.0;
