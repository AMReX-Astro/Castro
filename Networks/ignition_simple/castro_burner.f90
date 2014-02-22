module castro_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_data_module
  use eos_type_module
  use network

  private
  public :: burner

contains

  subroutine burner(dens, temp, Xin, ein, dt, time_in, Xout, eout)

    use burner_aux_module, only : dens_pass, c_p_pass, dhdx_pass, X_O16_pass

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), ein, dt, time_in
    real(kind=dp_t), intent(  out) :: Xout(nspec), eout
  
    real(kind=dp_t) :: enuc

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be temperature
    ! + the number of species
    integer, parameter :: NEQ = 1 + nspec_advance


    ! allocate storage for the input state
    real(kind=dp_t), dimension(NEQ) :: y


    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering
    integer, save :: ic12, io16, img24

    ! our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
    ! function is supplied, 22 means stiff, figure out my jacobian through 
    ! differencing
    integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22


    ! tolerance parameters:
    !
    !  itol specifies whether to use an single absolute tolerance for
    !  all variables (1), or to pass an array of absolute tolerances, one
    !  for each variable with a scalar relative tol (2), a scalar absolute
    !  and array of relative tolerances (3), or arrays for both (4)
    !  
    !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
    !  be > 0.  Since we have some compositions that may be 0 initially,
    !  we will specify both an absolute and a relative tolerance.
    !
    ! We will use arrays for both the absolute and relative tolerances, 
    ! since we want to be easier on the temperature than the species
    integer, parameter :: ITOL = 4
    real(kind=dp_t), dimension(NEQ) :: atol, rtol

    real(kind=dp_t) :: integration_time

    ! we want to do a normal computation, and get the output values of y(t)
    ! after stepping though dt
    integer, PARAMETER :: ITASK = 1

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate

    ! we will override the maximum number of steps, so turn on the 
    ! optional arguments flag
    integer, parameter :: IOPT = 1
    
    ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
    ! integer work array of since 30 + NEQ
    integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
    real(kind=dp_t), dimension(LRW) :: rwork
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    

    real(kind=dp_t) :: rpar
    integer :: ipar

    EXTERNAL jac, f_rhs
    
    logical, save :: firstCall = .true.

    type (eos_t) :: eos_state

    if (firstCall) then

       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif
     
       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")
       img24 = network_species_index("magnesium-24")
       
       if (ic12 < 0 .OR. io16 < 0 .OR. img24 < 0) then
          call bl_error("ERROR in burner: species undefined")
       endif
       
       firstCall = .false.
    endif

    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    atol(1:nspec_advance) = 1.d-12    ! mass fractions
    atol(nspec_advance+1) = 1.d-8     ! temperature
       
    rtol(1:nspec_advance) = 1.d-12    ! mass fractions
    rtol(nspec_advance+1) = 1.d-5     ! temperature
    

    ! we want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = ZERO
    iwork(:) = 0
    
    
    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000
    
    ! initialize the integration time
    integration_time = ZERO
    
    ! abundances are the first nspec_advance values and temperature is the last
    y(ic12) = Xin(ic12)
    y(nspec_advance+1) = temp

    ! density, specific heat at constant pressure, c_p, and dhdX are needed
    ! in the righthand side routine, so we will pass these in through the
    ! burner_aux module.
    !
    ! Since evaluating the EOS is expensive, we don't call it for every RHS
    ! call -- instead we freeze these values over the timestep.
    ! Since we are only integrating C12, we will need the O16 mass fraction
    ! in the RHS routine to compute the screening (and we know that the
    ! Mg24 abundance is constraint so things add to 1).
    dens_pass = dens
    
    ! we need the specific heat at constant pressure and dhdX |_p.  Take
    ! T, rho, Xin as input
    eos_state%rho = dens
    eos_state%T = temp
    eos_state%xn(:) = Xin(:)
       
    call eos(eos_input_rt, eos_state,.false.)
    
    c_p_pass = eos_state%cp
    dhdx_pass(:) = eos_state%dhdX(:)

    X_O16_pass = Xin(io16)

    ! call the integration routine
    call dvode(f_rhs, NEQ, y, integration_time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_ANALYTIC_JAC, &
               rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', integration_time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions -- note, we discard the temperature
    ! here and instead compute the energy release from the binding
    ! energy -- make sure that they are positive
    Xout(ic12)  = max(y(ic12), ZERO)
    Xout(io16)  = Xin(io16)
    Xout(img24) = ONE - Xout(ic12) - Xout(io16)

    ! compute the energy release and update the enthalpy.  Our convention
    ! is that the binding energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    !
    ! since this version of the network only evolves C12, we can
    ! compute the energy release easily
    enuc = (ebin(img24) - ebin(ic12))*(Xout(ic12) - Xin(ic12))
    
    eout = ein + enuc

    if (verbose) then
       
       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'dens: ', dens, ' temp: ', temp
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif

  end subroutine burner

end module castro_burner_module
