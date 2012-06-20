program testeos

  use bl_types
  use network
  use eos_module

  implicit none

  real(kind=dp_t) :: dens, temp, pres, entr
  real(kind=dp_t), dimension(nspec) :: Xin
  
  integer :: ic12, io16, img24

  logical :: do_diag

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  dens = 2.6e9_dp_t
  temp = 1.e9_dp_t

  pres = 1.7e27_dp_t
  entr = 6.5e7_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  den_eos = dens
  temp_eos = temp
  xn_eos(:) = Xin(:)

  do_diag = .false.

  call eos(eos_input_rt, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, 'eos_input_rt:'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos,  ' ener: ', e_eos
  print *, 'h:    ', h_eos,  ' entr: ', s_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)

  p_eos = pres
  s_eos = entr

  call eos(eos_input_ps, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, 'eos_input_ps:'
  print *, 'dens: ', den_eos, ' temp: ', temp_eos
  print *, 'X: ', Xin
  print *, 'pres: ', pres,  ' entr: ', entr
  print *, 'p_eos: ', p_eos, ' s_eos: ', s_eos
  print *, 'h:    ', h_eos,  ' ener: ', e_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)
  

end program testeos
