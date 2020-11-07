&fortin
  rho0 = 5.45969027775d-13
  T0 = 100.0d0
  v0 = 235435.230236d0
  rho1 = 1.24812608999d-12
  T1 = 207.756999533d0
  v1 = 102986.665199d0
  idir=1
/


&extern

  eos_gamma = 1.6666666667

  const_kappa_p = 3.92663697758d-5
  kappa_p_exp_m = 0.0d0
  kappa_p_exp_n = 0.0d0

  const_kappa_r = 0.848902853095d0
  kappa_r_exp_m = 0.0d0
  kappa_r_exp_n = 0.0d0

/

&tagging

  max_tempgrad_rel_lev = 20
  tempgrad_rel = 0.1

/
