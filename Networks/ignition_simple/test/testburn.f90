program testburn
  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module
  use castro_burner_module

  implicit none

  real(kind=dp_t) :: dens, temp, t, dt, ein, eout
  real(kind=dp_t), dimension(nspec) :: Xin, Xout
  
  integer :: ic12, io16, img24

  type (eos_t) :: eos_state

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  if (ic12 < 0 .or. io16 < 0 .or. img24 < 0) then
     call bl_error("ERROR: species index not defined")
  endif

  dens = 2.6e9_dp_t
  temp = 6.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t

  t = 0.0_dp_t
  dt = 0.06_dp_t

  eos_state%rho = dens
  eos_state%T = temp
  eos_state%xn(:) = Xin(:)

  call eos(eos_input_rt, eos_state)

  ein = eos_state%e

  print *, 'calling the burner...'


  call burner(dens, temp, Xin, ein, dt, t, Xout, eout)

  print *, 'done!'

  print *, 'Xin:  ', Xin
  print *, 'Xout: ', Xout
  print *, 'Hnuc (erg/g/s): ', eout-ein

end program testburn
