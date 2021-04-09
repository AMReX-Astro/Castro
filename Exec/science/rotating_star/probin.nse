&fortin

  model_name =  "15m_aprox19.6400"

/

&tagging

  max_dengrad_lev    = -1
  max_denerr_lev     =  3

  max_pressgrad_lev  = -1
  max_presserr_lev   = -1

  max_velgrad_lev    = -1
  max_velerr_lev     = -1

  max_temperr_lev    = -1
  max_tempgrad_lev   = -1

  dengrad            = 1.e18
  denerr             = 0.5e5

  velgrad            = 1.e18
  velerr             = 1.e18

  pressgrad          = 1.e27
  presserr           = 1.e27

  temperr            = 1.d20
  tempgrad           = 100.e0


/

&sponge

  sponge_upper_density = 1.d4
  sponge_lower_density = 1.d2
  sponge_timescale     = 1.d-3

/


&extern 

  rtol_spec = 1.d-6
  atol_spec = 1.d-6

  retry_burn = F
  abort_on_failure = F

  jacobian = 1

  rho_nse = 2.e6
  T_nse = 3.e9

/