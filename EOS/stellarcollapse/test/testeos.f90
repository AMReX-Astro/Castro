program testeos

  use bl_types
  use network
  use eos_module
  use eos_type_module
  use runtime_init_module
  use eos_aux_data_module
  use fundamental_constants_module
  use bl_constants_module

  implicit none

  integer :: iye
  
  logical :: do_diag

  type (eos_t) :: state

  call runtime_init()
  call network_init()
  call eos_init()

  iye = network_species_index('ye')
  if (iye < 0) call bl_error('bad iye')

  ! initial conditions
  state%rho = 1.0e13_dp_t
  state%T = 23_dp_t * MeV2eV * ev2erg / k_B
  state%p = ONE
  state%e = ONE
  state%aux(iye) = 0.2660725_dp_t

  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'
  print *, 'eos_input_rt'
  call eos(eos_input_rt,state)
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'
  
  state%T = state%T * 1.5
  print *, '----------------------------------------------------------------------'
  print *, "changed state%T to", state%T
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'
  print *, 'eos_input_re'
  call eos(eos_input_re,state)
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'

  state%T = state%T * 0.95
  print *, '----------------------------------------------------------------------'
  print *, "changed state%T to", state%T
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'
  print *, 'eos_input_rp'
  call eos(eos_input_rp,state)
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'

  state%rho = state%rho * 2.0
  print *, '----------------------------------------------------------------------'
  print *, "changed state%rho to", state%rho
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'
  print *, 'eos_input_tp'
  call eos(eos_input_tp,state)
  print *, '----------------------------------------------------------------------'
  print *, state%rho, state%T, state%aux
  print *, state%p, state%e, state%s
  print *, state%cs, state%gam1, state%dedT, state%dpdr_e
  print *, '======================================================================'


end program testeos
