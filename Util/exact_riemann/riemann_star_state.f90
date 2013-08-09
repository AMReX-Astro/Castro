module riemann_star_module

contains

subroutine riemann_star_state(rho_l, u_l, p_l, &
                              rho_r, u_r, p_r, &
                              xn_l, xn_r, &
                              ustar, pstar, W_l, W_r, &
                              verbose_in)

  use bl_types
  use bl_constants_module
  use eos_module
  use eos_type_module
  use network, only: nspec
  use riemann_support

  implicit none

  real (kind=dp_t), intent(in) :: rho_l, u_l, p_l
  real (kind=dp_t), intent(in) :: rho_r, u_r, p_r
  real (kind=dp_t), intent(in) :: xn_l(nspec), xn_r(nspec)
  real (kind=dp_t), intent(out) :: ustar, pstar, W_l, W_r
  logical, optional, intent(in) :: verbose_in



  real (kind=dp_t) :: cs_l, cs_r

  real (kind=dp_t) :: pstar_new, ustar_l, ustar_r

  real (kind=dp_t) :: rho_tmp, u_tmp, p_tmp
  real (kind=dp_t) :: rhostar

  real (kind=dp_t) :: csbar, rhobar
  real (kind=dp_t) :: cs_star

  real (kind=dp_t) :: gammaE_l, gammaC_l, gammaE_r, gammaC_r
  real (kind=dp_t) :: gammaE_bar, gammaC_bar

  real (kind=dp_t) :: Z_l, Z_r
  real (kind=dp_t) :: W_temp, Z_temp

  real (kind=dp_t) :: err1, err2
  real (kind=dp_t), parameter :: tol = 1.e-10_dp_t
  logical :: converged

  integer :: iter
  integer, parameter :: max_iter = 10

  type (eos_t) :: eos_state

  character (len=16) :: lwave, rwave

  real (kind=dp_t), parameter :: smallp = 1.d-8
  real (kind=dp_t), parameter :: SMALL = 1.d-13

  logical, parameter :: debug = .false.

  logical :: verbose

  if (present(verbose_in)) then
     verbose = verbose_in
  else
     verbose = .false.
  endif


  ! get the initial sound speeds
  eos_state%rho = rho_l
  eos_state%p = p_l
  eos_state%xn(:) = xn_l(:)
  eos_state%T = 100000.0_dp_t   ! initial guess

  call eos(eos_input_rp, eos_state, .false.)

  cs_l = sqrt(eos_state%gam1*p_l/rho_l)

  if (verbose) print *, 'T_l = ', eos_state%T

  gammaE_l = p_l/(rho_l*eos_state%e) + ONE
  gammaC_l = eos_state%gam1

  eos_state%rho = rho_r
  eos_state%p = p_r
  eos_state%xn(:) = xn_r(:)
  eos_state%T = 100000.0_dp_t   ! initial guess

  call eos(eos_input_rp, eos_state, .false.)

  cs_r = sqrt(eos_state%gam1*p_r/rho_r)

  if (verbose) print *, 'T_r = ', eos_state%T

  gammaE_r = p_r/(rho_r*eos_state%e) + ONE
  gammaC_r = eos_state%gam1

  gammaE_bar = HALF*(gammaE_l + gammaE_r)
  gammaC_bar = HALF*(gammaC_l + gammaC_r)

  !---------------------------------------------------------------------------
  ! create an initial guess for pstar using a primitive variable
  ! Riemann solver
  ! ---------------------------------------------------------------------------

  ! We follow the PVRS solver from Toro (Chapter 9)

  ! alternative: two shock solver (see Toro 9.42)
  W_l = rho_l*cs_l
  W_r = rho_r*cs_r

  ! prevent roundoff errors from giving us a pstar that is unphysical
  ! if our input states are the same
  if (W_l == W_r) then
     pstar = 0.5d0*(p_l + p_r + W_l*(u_l - u_r))
  else
     pstar = ((W_r*p_l + W_l*p_r) + W_l*W_r*(u_l - u_r))/(W_l + W_r)
  endif


  ! if (abs(p_l - pstar) < SMALL*pstar .and. abs(p_r - pstar) < SMALL*pstar) then
  !    rhobar = HALF*(rho_l + rho_r)
  !    csbar = HALF*(cs_l + cs_r)

  !    pstar = HALF*(p_l + p_r) + HALF*(u_l - u_r)*rhobar*csbar
  ! endif


  pstar = max(pstar, smallp)


  !---------------------------------------------------------------------------
  ! find the exact pstar and ustar
  !---------------------------------------------------------------------------

  !if (verbose) 
  print *, 'solving for star state: ', rho_l, u_l, p_l, rho_r, u_r, p_r
  if (verbose) print *, 'pstar: ', pstar, (pstar-p_l)/pstar, (pstar-p_r)/pstar

  ! this procedure follows directly from Colella & Glaz 1985, section 1

  converged = .false.
  iter = 1
  do while (.not. converged .and. iter < max_iter)

     ! compute Z_l and Z_r -- the form of these depend on whether the
     ! wave is a shock or a rarefaction
     
     ! left wave
     if (pstar - p_l > SMALL*p_l) then
        ! left shock
        call shock(pstar, rho_l, u_l, p_l, xn_l, gammaE_bar, gammaC_bar, Z_l, W_l)
        lwave = "shock"
     else
        ! left rarefaction
        call rarefaction(pstar, rho_l, u_l, p_l, xn_l, 1, Z_l, W_l)
        lwave = "rarefaction"
     endif

     ! right wave
     if (pstar - p_r > SMALL*p_r) then
        ! right shock
        call shock(pstar, rho_r, u_r, p_r, xn_r, gammaE_bar, gammaC_bar, Z_r, W_r)
        rwave = "shock"
     else
        ! right rarefaction
        call rarefaction(pstar, rho_r, u_r, p_r, xn_r, 3, Z_r, W_r)
        rwave = "rarefaction"
     endif

     if (verbose) print *, 'left wave: ', trim(lwave), ';   rightwave: ', trim(rwave)

     ustar_l = u_l - (pstar - p_l)/W_l
     ustar_r = u_r + (pstar - p_r)/W_r

     pstar_new = pstar - Z_l*Z_r*(ustar_r - ustar_l)/(Z_l + Z_r)

     if (verbose) print *, "done with iteration", iter
     if (verbose) print *, "ustar_l/r, pstar: ", &
          real(ustar_l), real(ustar_r), real(pstar_new)
     if (verbose) print * , " "

     ! estimate the error in the current star solution
     err1 = abs(ustar_r - ustar_l)
     err2 = pstar_new - pstar

     if (verbose) print *, "ERRORS: ", err1, err2
     if (verbose) print *, " "

     if (err1 < tol*max(abs(ustar_l),abs(ustar_r)) .and. err2 < tol*pstar) then
        converged = .true.
     endif

     ! get ready for the next iteration
     pstar = pstar_new

     iter = iter + 1
     
  enddo

  ustar = HALF*(ustar_l + ustar_r)

  if (verbose) print *, 'found pstar, ustar: ', pstar, ustar

  ! let's test if our integration across the rarefaction works, as expected
  if (lwave == "rarefaction" .and. debug) then
     call rarefaction(pstar, rho_l, u_l, p_l, xn_l, 1, Z_temp, W_temp, rhostar)

     ! get the soundspeed via the EOS (C&G suggest getting it from
     ! the jump conditions)
     eos_state%rho = rhostar
     eos_state%p = pstar
     eos_state%xn(:) = xn_l(:)
     eos_state%T = 100000.0_dp_t   ! initial guess
     
     call eos(eos_input_rp, eos_state, .false.)
        
     cs_star = sqrt(eos_state%gam1*pstar/rhostar)

     call rarefaction_to_u(rho_l, u_l, p_l, xn_l, 1, ustar-cs_star, rho_tmp, p_tmp, u_tmp)

     print *, "left rarefaction integration test"
     print *, "(this should agree with pstar/ustar if the 1-wave is a rarefaction): ", p_tmp, u_tmp

  endif

  if (rwave == "rarefaction" .and. debug) then
     call rarefaction(pstar, rho_r, u_r, p_r, xn_r, 3, Z_temp, W_temp, rhostar)

     ! get the soundspeed via the EOS (C&G suggest getting it from
     ! the jump conditions)
     eos_state%rho = rhostar
     eos_state%p = pstar
     eos_state%xn(:) = xn_r(:)
     eos_state%T = 100000.0_dp_t   ! initial guess
     
     call eos(eos_input_rp, eos_state, .false.)
        
     cs_star = sqrt(eos_state%gam1*pstar/rhostar)

     call rarefaction_to_u(rho_r, u_r, p_r, xn_r, 3, ustar+cs_star, rho_tmp, p_tmp, u_tmp)

     print *, "left rarefaction integration test"
     print *, "(this should agree with pstar/ustar if the 3-wave is a rarefaction): ", p_tmp, u_tmp

  endif

end subroutine riemann_star_state
end module riemann_star_module
