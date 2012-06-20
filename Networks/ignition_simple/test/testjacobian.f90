! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module
  use burner_module
  use burner_aux_module

  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(nspec_advance+1) :: y, ydot
  real(kind=dp_t), dimension(nspec_advance+1) :: yp, ym
  real(kind=dp_t), dimension(nspec_advance+1) :: ydotp, ydotm
  real(kind=dp_t), dimension(nspec_advance+1,nspec_advance+1) :: pd
  real(kind=dp_t) :: enucdot

  real(kind=dp_t) :: rpar
  integer :: ipar

  integer :: ic12, io16, img24
  integer :: i, j, n

  real(kind=dp_t), parameter :: delta = 0.001d0
  real(kind=dp_t) :: num_jac

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  if (ic12 < 0 .or. io16 < 0 .or. img24 < 0) then
     call bl_error("ERROR: species index not defined")
  endif
  
  dens = 2.6e9_dp_t
  temp = 7.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  den_eos(1) = dens
  temp_eos(1) = temp
  xn_eos(1,:) = Xin(:)
  
  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)


  print *, 'evaluating the RHS...'

  ! load the state
  y(1) = Xin(ic12)
  y(nspec_advance+1) = temp

  ! set the burner_aux variables
  dens_pass = dens
  c_p_pass = cp_eos(1)
  dhdx_pass(:) = dhdX_eos(1,:)
  X_O16_pass = Xin(io16)
  

  call f_rhs(nspec_advance+1, ZERO, y, ydot, rpar, ipar)

  call jac(nspec_advance+1, ZERO, y, 0, 0, pd, nspec_advance+1, rpar, ipar)

999 format(1x, "df(",i1,")/dy(",i1,")", g20.10, g20.10)

  do j = 1, nspec_advance+1

     yp(:) = y(:)
     ym(:) = y(:)

     yp(j) = (1.d0 + delta)*y(j)
     call f_rhs(nspec_advance+1, ZERO, yp, ydotp, rpar, ipar)
     
     ym(j) = (1.d0 - delta)*y(j)
     call f_rhs(nspec_advance+1, ZERO, ym, ydotm, rpar, ipar)        

     do i = 1, nspec_advance+1
        
        num_jac = (ydotp(i) - ydotm(i))/(yp(j) - ym(j))

        write (*,999) i, j, num_jac, pd(i,j)

     enddo
  enddo

end program testjacobian
