&fortin

  rho_fuel = 2.d7
  T_fuel = 5.e7

  T_ash = 4.e9

  pert_frac = 0.15d0
  pert_delta = 0.02d0

/

&tagging

  denerr = 1.d-7
  dengrad = 0.01
  max_denerr_lev = 5
  max_dengrad_lev = 5

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = 5
  max_pressgrad_lev = 5

/

&sponge

  sponge_upper_density = 5.0d-8
  sponge_lower_density = 1.0d-8
  sponge_timescale     = 1.0d-6

/

&extern
  rtol_spec = 1.d-10
  atol_spec = 1.d-10

/
