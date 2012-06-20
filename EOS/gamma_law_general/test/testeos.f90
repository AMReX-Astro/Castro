program testburn

  use bl_types
  use network
  use eos_module

  implicit none

  real(kind=dp_t) :: dens, temp, pres, entr
  real(kind=dp_t), dimension(nspec) :: Xin
  
  integer :: ic12, io16, img24

  logical :: do_diag

  call network_init()
  call eos_init(gamma_in=5.0d0/3.0d0)

  do_diag = .false.

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  dens = 1.e4_dp_t
  temp = 1.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  den_eos(1) = dens
  temp_eos(1) = temp
  xn_eos(1,:) = Xin(:)

  
  ! test eos_input_rt
  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, 'input: eos_input_rt'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos(1),  ' ener: ', e_eos(1)
  print *, 'h:    ', h_eos(1)
  print *, 'c_v:  ', cv_eos(1), ' c_p : ', cp_eos(1)
  print *, 'dpdT: ', dpdt_eos(1), ' dpdr: ', dpdr_eos(1)
  print *, 'dedT: ', dedt_eos(1), ' dedr: ', dedr_eos(1)
  print *, 'dpdX: ', dpdX_eos(1,:)
  print *, 'dhdX: ', dhdX_eos(1,:)

  print *, ' '
  print *, 'setting pres = ', p_eos(1)
  pres = p_eos(1)
  print *, 'setting entr = ', s_eos(1)
  entr = s_eos(1)

  
  ! test eos_input_rh
  temp_eos(1) = 0.d0
  call eos(eos_input_rh, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_rh'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos(1),  ' ener: ', e_eos(1)
  print *, 'entropy: ', s_eos(1)
  print *, 'h:    ', h_eos(1)
  print *, 'c_v:  ', cv_eos(1), ' c_p : ', cp_eos(1)
  print *, 'dpdT: ', dpdt_eos(1), ' dpdr: ', dpdr_eos(1)
  print *, 'dedT: ', dedt_eos(1), ' dedr: ', dedr_eos(1)
  print *, 'dpdX: ', dpdX_eos(1,:)
  print *, 'dhdX: ', dhdX_eos(1,:)

  ! test eos_input_tp
  den_eos(1) = 0.d0
  call eos(eos_input_tp, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_tp'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos(1),  ' ener: ', e_eos(1)
  print *, 'h:    ', h_eos(1)
  print *, 'c_v:  ', cv_eos(1), ' c_p : ', cp_eos(1)
  print *, 'dpdT: ', dpdt_eos(1), ' dpdr: ', dpdr_eos(1)
  print *, 'dedT: ', dedt_eos(1), ' dedr: ', dedr_eos(1)
  print *, 'dpdX: ', dpdX_eos(1,:)
  print *, 'dhdX: ', dhdX_eos(1,:)

  ! test eos_input_rp
  temp_eos(1) = 0.d0
  call eos(eos_input_rp, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_rp'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos(1),  ' ener: ', e_eos(1)
  print *, 'h:    ', h_eos(1)
  print *, 'c_v:  ', cv_eos(1), ' c_p : ', cp_eos(1)
  print *, 'dpdT: ', dpdt_eos(1), ' dpdr: ', dpdr_eos(1)
  print *, 'dedT: ', dedt_eos(1), ' dedr: ', dedr_eos(1)
  print *, 'dpdX: ', dpdX_eos(1,:)
  print *, 'dhdX: ', dhdX_eos(1,:)

  ! test eos_input_re
  temp_eos(1) = 0.d0
  call eos(eos_input_re, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_re'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos(1),  ' ener: ', e_eos(1)
  print *, 'h:    ', h_eos(1)
  print *, 'c_v:  ', cv_eos(1), ' c_p : ', cp_eos(1)
  print *, 'dpdT: ', dpdt_eos(1), ' dpdr: ', dpdr_eos(1)
  print *, 'dedT: ', dedt_eos(1), ' dedr: ', dedr_eos(1)
  print *, 'dpdX: ', dpdX_eos(1,:)
  print *, 'dhdX: ', dhdX_eos(1,:)


  ! test eos_input_ps
  p_eos(1) = pres
  s_eos(1) = entr
  call eos(eos_input_ps, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_ps'
  print *, 'pres: ', pres, 'entr: ', entr
  print *, 'dens: ', den_eos(1), ' temp: ', temp_eos(1)
  print *, 'X: ', Xin
  print *, 'pres_eos: ', p_eos(1),  'entr_eos: ', s_eos(1)
  print *, 'h:    ', h_eos(1), ' ener: ', e_eos(1)
  print *, 'c_v:  ', cv_eos(1), ' c_p : ', cp_eos(1)
  print *, 'dpdT: ', dpdt_eos(1), ' dpdr: ', dpdr_eos(1)
  print *, 'dedT: ', dedt_eos(1), ' dedr: ', dedr_eos(1)
  print *, 'dpdX: ', dpdX_eos(1,:)
  print *, 'dhdX: ', dhdX_eos(1,:)

end program testburn
