#include <wdmerger_data.H>

using namespace amrex;

AMREX_GPU_MANAGED int wdmerger::relaxation_is_done = 0;
AMREX_GPU_MANAGED int wdmerger::problem = -1;
AMREX_GPU_MANAGED int wdmerger::use_stopping_criterion = 1;
AMREX_GPU_MANAGED int wdmerger::use_energy_stopping_criterion = 0;
AMREX_GPU_MANAGED Real wdmerger::ts_te_stopping_criterion = 1.e200;
AMREX_GPU_MANAGED Real wdmerger::T_stopping_criterion = 1.e200;

AMREX_GPU_MANAGED Real wdmerger::mass_p = 0.0;
AMREX_GPU_MANAGED Real wdmerger::mass_s = 0.0;

AMREX_GPU_MANAGED Real wdmerger::mdot_p = 0.0;
AMREX_GPU_MANAGED Real wdmerger::mdot_s = 0.0;

AMREX_GPU_MANAGED Real wdmerger::com_p[3] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::com_s[3] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::vel_p[3] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::vel_s[3] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::rad_p[7] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::rad_s[7] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::vol_p[7] = { 0.0 };
AMREX_GPU_MANAGED Real wdmerger::vol_s[7] = { 0.0 };

AMREX_GPU_MANAGED Real wdmerger::rho_avg_p = 0.0;
AMREX_GPU_MANAGED Real wdmerger::rho_avg_s = 0.0;

AMREX_GPU_MANAGED Real wdmerger::t_ff_p = 0.0;
AMREX_GPU_MANAGED Real wdmerger::t_ff_s = 0.0;

AMREX_GPU_MANAGED Real wdmerger::T_global_max = 0.0;
AMREX_GPU_MANAGED Real wdmerger::rho_global_max = 0.0;
AMREX_GPU_MANAGED Real wdmerger::ts_te_global_max = 0.0;

AMREX_GPU_MANAGED Real wdmerger::T_curr_max = 0.0;
AMREX_GPU_MANAGED Real wdmerger::rho_curr_max = 0.0;
AMREX_GPU_MANAGED Real wdmerger::ts_te_curr_max = 0.0;

AMREX_GPU_MANAGED Real wdmerger::total_ener_array[num_previous_ener_timesteps] = { 0.0 };
