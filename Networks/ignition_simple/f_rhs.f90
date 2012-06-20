subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  
  ! we get the thermodynamic state through the burner_aux module -- we freeze
  ! these to the values are the top of the timestep to avoid costly
  ! EOS calls
  use burner_aux_module, only : dens_pass, c_p_pass, dhdx_pass, X_O16_pass

  implicit none

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  integer :: n
  real(kind=dp_t) :: y(n), ydot(n)

  integer :: k
  real(kind=dp_t) :: ymass(nspec)

  real(kind=dp_t) :: rpar
  integer :: ipar

  real(kind=dp_t) :: t

  real(kind=dp_t) :: dens, temp, T9, T9a, dT9dt, dT9adt

  real(kind=dp_t) :: rate, dratedt
  real(kind=dp_t) :: sc1212, dsc1212dt
  real(kind=dp_t) :: xc12tmp
  common /rate_info/ rate, dratedt, sc1212, dsc1212dt, xc12tmp

  real(kind=dp_t), PARAMETER :: &
                     one_twelvth = 1.0d0/12.0d0, &
                     five_sixths = 5.0d0/ 6.0d0, &
                       one_third = 1.0d0/ 3.0d0, &
                      two_thirds = 2.0d0/ 3.0d0

  real(kind=dp_t) :: scratch, dscratchdt
  
  real(kind=dp_t) :: a, b, dadt, dbdt

  integer, save :: ic12, io16, img24

  logical, save :: firstCall = .true.

  if (firstCall) then
     ic12 = network_species_index("carbon-12")
     io16 = network_species_index("oxygen-16")
     img24 = network_species_index("magnesium-24")

     firstCall = .false.
  end if



  ! compute the molar fractions -- needed for the screening
  ymass(ic12) = y(1)/aion(ic12)
  ymass(io16) = X_O16_pass/aion(io16)
  ymass(img24) = (ONE - y(1) - X_O16_pass)/aion(img24)

  ! for convinence, carry a temp variable and dens variable
  dens = dens_pass
  temp = y(nspec_advance+1)

  ! call the screening routine
  call screenz(temp,dens,6.0d0,6.0d0,12.0d0,12.0d0,ymass,aion,zion,nspec,     &
               sc1212, dsc1212dt)

  
  ! compute some often used temperature constants
  T9     = temp/1.d9
  dT9dt  = ONE/1.d9
  T9a    = T9/(1.0d0 + 0.0396d0*T9)
  dT9adt = (T9a / T9 - (T9a / (1.0d0 + 0.0396d0*T9)) * 0.0396d0) * dT9dt

  ! compute the CF88 rate
  scratch    = T9a**one_third
  dscratchdt = one_third * T9a**(-2.0d0 * one_third) * dT9adt

  a       = 4.27d26*T9a**five_sixths*T9**(-1.5d0)
  dadt    = five_sixths * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt

  b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
  dbdt    = (84.165d0 * dscratchdt/ scratch**2.0d0                            &
             - 3.0d0 * 2.12d-3 * T9 * T9 * dT9dt) * b

  rate    = a *  b
  dratedt = dadt * b + a * dbdt

  ! The change in number density of C12 is
  ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
  !
  ! where <sigma v> is the average of the relative velocity times the cross
  ! section for the reaction, and the factor accounting for the total number
  ! of particle pairs has a 1/2 because we are considering a reaction involving 
  ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
  ! each reaction, we lose 2 carbon nuclei.
  !
  ! The corresponding Mg24 change is
  ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
  !
  ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
  !
  ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
  ! Avagadro's number, and A is the mass number of the nucleon, we get
  !
  ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
  !
  ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
  !
  ! these are equal and opposite.
  !
  ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

  ! we will always refer to the species by integer indices that come from
  ! the network module -- this makes things robust to a shuffling of the 
  ! species ordering

  xc12tmp = max(y(ic12),0.d0)
  ydot(ic12) = -one_twelvth*dens*sc1212*rate*xc12tmp**2

  ! now compute the change in temperature, using the evolution equation
  ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
  ! 
  ! we make use of the fact that omega(Mg24) = - omega(C12), and that
  ! omega(O16) = 0 in our simplified burner
  ydot(nspec_advance+1) =  ( (dhdx_pass(img24) - dhdx_pass(ic12)) + &
                             (ebin(img24) - ebin(ic12)) )*ydot(ic12)/c_p_pass

  return

end subroutine f_rhs


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network

  ! we get the thermodynamic state through the burner_aux module -- we freeze
  ! these to the values are the top of the timestep to avoid costly
  ! EOS calls
  use burner_aux_module, only : dens_pass, c_p_pass, dhdx_pass, X_O16_pass

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(kind=dp_t), intent(IN   ) :: y(neq), rpar, t
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

  real(kind=dp_t) :: rate, dratedt, scorr, dscorrdt, xc12tmp
  common /rate_info/ rate, dratedt, scorr, dscorrdt, xc12tmp

  integer :: itemp

  integer, save :: ic12, io16, img24

  logical, save :: firstCall = .true.

  if (firstCall) then
     ic12 = network_species_index("carbon-12")
     io16 = network_species_index("oxygen-16")
     img24 = network_species_index("magnesium-24")

     firstCall = .false.
  end if


  ! initialize
  pd(:,:)  = ZERO

  itemp = nspec_advance + 1


  ! carbon jacobian elements
  pd(ic12, ic12) = -(1.0d0/6.0d0)*dens_pass*scorr*rate*xc12tmp



  ! add the temperature derivatives: df(y_i) / dT
  pd(ic12,itemp) = -(1.0d0/12.0d0)*(dens_pass*rate*xc12tmp**2*dscorrdt    &
                                    + dens_pass*scorr*xc12tmp**2*dratedt)  

  ! add the temperature jacobian elements df(T) / df(y)
  pd(itemp,ic12) =  ( (dhdx_pass(img24) - dhdx_pass(ic12)) + &
                      (ebin(img24) - ebin(ic12)) )*pd(ic12,ic12)/c_p_pass


  ! add df(T) / dT
  pd(itemp,itemp) = ( (dhdx_pass(img24) - dhdx_pass(ic12)) + &
                      (ebin(img24) - ebin(ic12)) )*pd(ic12,itemp)/c_p_pass


  return
end subroutine jac

