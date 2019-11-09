&fortin

  rho_fuel = 2.e7
  T_fuel = 5.e7

  T_ash = 3.6e9

  ash1_name = "nickel-56"
  X_ash1 = 1.0

  smallx_init = 1.d-8

  pert_frac = 0.4d0
  pert_delta = 0.06d0

  v_inflow = 0.0d0

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

  small_x = 1.e-10
/
