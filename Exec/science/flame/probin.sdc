&fortin

  rho_fuel = 5.d6
  T_fuel = 5.e7

  T_ash = 2.5e9

  ash1_name = "silicon-28"
  X_ash1 = 0.25

  ash2_name = "sulfur-32"
  X_ash2 = 0.25

  ash3_name = "argon-36"
  X_ash3 = 0.25

  ash4_name = "calcium-40"
  X_ash4 = 0.25

  pert_frac = 0.3d0
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
