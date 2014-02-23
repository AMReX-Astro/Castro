program testeos

  use bl_types
  use network
  use eos_module

  implicit none

  real(kind=dp_t) :: dens_good, temp_good, pres_good, entr_good, enth_good, eint_good
  real(kind=dp_t), dimension(nspec) :: Xin
  
  integer :: ic12, io16, img24

  logical :: do_diag

  real(kind=dp_t) :: err1, err2

  type (eos_t) :: state

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  dens_good = 2.6e9_dp_t
  temp_good = 1.e9_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t

  do_diag = .false.

  !---------------------------------------------------------------------------
  ! get the initial state -- this will be considered the "right" answer --
  ! make sure we can get these values back when using the other input 
  ! methods
  !---------------------------------------------------------------------------

  state % rho = dens_good
  state % T   = temp_good
  state % xn  = Xin

  call eos(eos_input_rt, state, do_diag)

  pres_good = state%p
  entr_good = state%s
  enth_good = state%h
  eint_good = state%e

  print *, 'eos_input_rt:'
  print *, 'dens: ', state%rho, ' temp: ', state%T
  print *, 'X: ', Xin
  print *, 'pres: ', state % p,  ' ener: ', state % e
  print *, 'h:    ', state % h,  ' entr: ', state % s
  print *, 'c_v:  ', state % cv, ' c_p : ', state % cp
  print *, 'dpdT: ', state % dpdT, ' dpdr: ', state % dpdr
  print *, 'dedT: ', state % dedT, ' dedr: ', state % dedr
  print *, 'dpdX: ', state % dpdX
  print *, 'dhdX: ', state % dhdX

  
  !---------------------------------------------------------------------------
  ! now try the other input methods and make sure we get back what we 
  ! started with
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! density, enthalpy
  !---------------------------------------------------------------------------
  state%rho = dens_good
  state%h = enth_good

  state%T = 0.1*temp_good

  call eos(eos_input_rh, state, do_diag)

  ! compute the error in T
  err1 = abs(temp_good - state%T)/temp_good
  print *, " "
  print *, "eos_input_rh, err: ", err1


  !---------------------------------------------------------------------------
  ! temp, pres
  !---------------------------------------------------------------------------
  state%T = temp_good
  state%p = pres_good

  state%rho = 0.1*dens_good

  call eos(eos_input_tp, state, do_diag)

  ! compute the error in rho
  err1 = abs(dens_good - state%rho)/dens_good
  print *, " "
  print *, "eos_input_tp, err: ", err1
  

  !---------------------------------------------------------------------------
  ! dens, pres
  !---------------------------------------------------------------------------
  state%rho = dens_good
  state%p = pres_good

  state%T = 0.1*temp_good

  call eos(eos_input_rp, state, do_diag)

  ! compute the error in T
  err1 = abs(temp_good - state%T)/temp_good
  print *, " "
  print *, "eos_input_rp, err: ", err1


  !---------------------------------------------------------------------------
  ! dens, energy
  !---------------------------------------------------------------------------
  state%rho = dens_good
  state%e = eint_good

  state%T = 0.1*temp_good

  call eos(eos_input_re, state, do_diag)

  ! compute the error in T
  err1 = abs(temp_good - state%T)/temp_good
  print *, " "
  print *, "eos_input_re, err: ", err1


  !---------------------------------------------------------------------------
  ! temp, enthalpy
  !---------------------------------------------------------------------------
  state%rho = dens_good
  state%h = enth_good

  state%rho = 0.1*dens_good

  call eos(eos_input_th, state, do_diag)

  ! compute the error in rho
  err1 = abs(dens_good - state%rho)/dens_good
  print *, " "
  print *, "eos_input_th, err: ", err1


  !---------------------------------------------------------------------------
  ! pres, entropy
  !---------------------------------------------------------------------------
  state%p = pres_good
  state%s = entr_good

  state%T = 0.1*temp_good
  state%rho = 0.1*dens_good

  call eos(eos_input_ps, state, do_diag)

  ! compute the error in T and rho
  err1 = abs(temp_good - state%T)/temp_good
  err2 = abs(dens_good - state%rho)/dens_good
  print *, " "
  print *, "eos_input_ps, err: ", err1, err2


  !---------------------------------------------------------------------------
  ! pres, enthalpy
  !---------------------------------------------------------------------------
  state%p = pres_good
  state%h = enth_good

  state%T = 0.1*temp_good
  state%rho = 0.1*dens_good

  call eos(eos_input_ph, state, do_diag)

  ! compute the error in T and rho
  err1 = abs(temp_good - state%T)/temp_good
  err2 = abs(dens_good - state%rho)/dens_good
  print *, " "
  print *, "eos_input_ph, err: ", err1, err2



end program testeos
