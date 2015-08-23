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

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), ein, dt, time_in
    real(kind=dp_t), intent(  out) :: Xout(nspec), eout

    real(kind=dp_t) :: enuc

    logical, parameter :: verbose = .false.

    ! set the number of independent variable
    integer, parameter :: NEQ = nspec_advance


    ! allocate storage for the input state
    real(kind=dp_t), dimension(NEQ) :: y


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
    ! We will use arrays for both the absolute and relative tolerances
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

       firstCall = .false.
    endif

    ! set the tolerances. 
    atol(1:nspec_advance) = 1.d-12    ! mass fractions

    rtol(1:nspec_advance) = 1.d-12    ! mass fractions


    ! we want VODE to re-initialize each time we call it
    istate = 1

    rwork(:) = ZERO
    iwork(:) = 0


    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000

    ! initialize the integration time
    integration_time = ZERO

    y(:) = Xin(:)

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


    ! store the new mass fractions -- make sure that they are positive
    Xout(_ifuel)= max(y(_ifuel), ZERO)
    Xout(_iash) = min(y(_iash), ONE)

    ! compute the energy release and update the enthalpy.  Our convention
    ! is that the binding energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    !
    ! since this version of the network only evolves C12, we can
    ! compute the energy release easily
    enuc = A_burning*(Xout(_ifuel) - Xin(_ifuel))

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
