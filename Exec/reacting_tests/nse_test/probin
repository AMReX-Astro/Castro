&fortin
  T_l = 4.d9
  T_r = 5.d7

  dens = 2.d8
  cfrac = 0.d0
  
  smallx = 1.e-10

  idir = 1

  w_T = 5.d-4
  center_T = 3.d-1
/

&tagging
  denerr = 3
  dengrad = 0.01
  max_denerr_lev = 0
  max_dengrad_lev = 0

  presserr = 3
  pressgrad = 0.01
  max_presserr_lev = 0
  max_pressgrad_lev = 0

  temperr = 4.e9
  tempgrad = 1.e8
  max_temperr_lev = 5
  max_tempgrad_lev = 5
/

&extern
!  jacobian = 2
!  centered_diff_jac = T
!   rtol_enuc = 1.e-8
   call_eos_in_rhs = T
   do_constant_volume_burn = T
/
