program testeos

  use bl_types
  use network
  use eos_module

  implicit none

  real(kind=dp_t) :: dens, temp, pres, entr
  real(kind=dp_t), dimension(nspec) :: Xin
  
  integer :: ic12, io16, img24

  logical :: do_diag

  type (eos_t) :: state

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

  do_diag = .false.

  state % rho = dens
  state % T   = temp
  state % xn  = Xin

  call eos(eos_input_rt, state, do_diag)

  print *, 'eos_input_rt:'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', state % p,  ' ener: ', state % e
  print *, 'h:    ', state % h,  ' entr: ', state % s
  print *, 'c_v:  ', state % cv, ' c_p : ', state % cp
  print *, 'dpdT: ', state % dpdT, ' dpdr: ', state % dpdr
  print *, 'dedT: ', state % dedT, ' dedr: ', state % dedr
  print *, 'dpdX: ', state % dpdX
  print *, 'dhdX: ', state % dhdX

  call eos(eos_input_ps, state, do_diag)

  print *, 'eos_input_ps:'
  print *, 'dens: ', state % rho, ' temp: ', state % T
  print *, 'X: ', Xin
  print *, 'pres: ', state % p,  ' ener: ', state % e
  print *, 'h:    ', state % h,  ' entr: ', state % s
  print *, 'c_v:  ', state % cv, ' c_p : ', state % cp
  print *, 'dpdT: ', state % dpdT, ' dpdr: ', state % dpdr
  print *, 'dedT: ', state % dedT, ' dedr: ', state % dedr
  print *, 'dpdX: ', state % dpdX
  print *, 'dhdX: ', state % dhdX

end program testeos
