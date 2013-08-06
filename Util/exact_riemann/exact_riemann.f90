program riemann_exact

  use bl_types
  use bl_constants_module
  use eos_module
  use eos_type_module
  use network, only: nspec
  use probin_module, only: rho_l, u_l, p_l, T_l, rho_r, u_r, p_r, T_r, &
                           xmin, xmax, xjump, t, npts, use_Tinit
  use riemann_support
  use runtime_init_module

  implicit none

  real (kind=dp_t) :: xn(nspec)

  real (kind=dp_t) :: cs_l, cs_r

  real (kind=dp_t) :: pstar, ustar
  real (kind=dp_t) :: pstar_new, ustar_l, ustar_r

  real (kind=dp_t) :: csbar, rhobar
  
  real (kind=dp_t) :: gammaE_l, gammaC_l, gammaE_r, gammaC_r
  real (kind=dp_t) :: gammaE_bar, gammaC_bar

  real (kind=dp_t) :: W_l, W_r, Z_l, Z_r
  real (kind=dp_t) :: W_temp, Z_temp

  real (kind=dp_t) :: err1, err2
  real (kind=dp_t), parameter :: tol = 1.e-10_dp_t
  logical :: converged

  integer :: iter
  integer, parameter :: max_iter = 10

  type (eos_t) :: eos_state

  character (len=16) :: lwave, rwave

  real (kind=dp_t) :: x, dx, xi, chi
  real (kind=dp_t) :: rhostar
  real (kind=dp_t) :: xihat, uhat_s, uhat_star
  real (kind=dp_t) :: lambdahat_s, lambdahat_star
  real (kind=dp_t) :: cs_star
  real (kind=dp_t) :: rho_s, p_s, u_s, W_s, cs_s
  real (kind=dp_t) :: rho, u, p

  real (kind=dp_t) :: smallp = 1.d-8

  integer :: i, lun

  ! general Maestro initializations
  call runtime_init()

  ! microphysics
  call network_init()
  call eos_init(gamma_in=1.4d0)

  ! we need a composition to interface with our EOS, but we are not
  ! exploring composition jumps here.  We'll take a constant
  ! composition.
  xn(:) = ZERO
  xn(1) = ONE


  ! if we are using T as the independent variable (rather than p), then
  ! get p now
  if (use_Tinit) then
     eos_state%rho = rho_l
     eos_state%T = T_l
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state, .false.)

     p_l = eos_state%p

     print *, 'p_l = ', p_l

     eos_state%rho = rho_r
     eos_state%T = T_r
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state, .false.)

     p_r = eos_state%p

     print *, 'p_r = ', p_r
  endif


  ! get the initial sound speeds
  eos_state%rho = rho_l
  eos_state%p = p_l
  eos_state%xn(:) = xn(:)
  eos_state%T = 100000.0_dp_t   ! initial guess

  call eos(eos_input_rp, eos_state, .false.)

  cs_l = sqrt(eos_state%gam1*p_l/rho_l)

  print *, 'T_l = ', eos_state%T

  gammaE_l = p_l/(rho_l*eos_state%e) + ONE
  gammaC_l = eos_state%gam1

  eos_state%rho = rho_r
  eos_state%p = p_r
  eos_state%xn(:) = xn(:)
  eos_state%T = 100000.0_dp_t   ! initial guess

  call eos(eos_input_rp, eos_state, .false.)

  cs_r = sqrt(eos_state%gam1*p_r/rho_r)

  print *, 'T_r = ', eos_state%T

  gammaE_r = p_r/(rho_r*eos_state%e) + ONE
  gammaC_r = eos_state%gam1

  gammaE_bar = HALF*(gammaE_l + gammaE_r)
  gammaC_bar = HALF*(gammaC_l + gammaC_r)

  !---------------------------------------------------------------------------
  ! create an initial guess for pstar using a primitive variable
  ! Riemann solver
  ! ---------------------------------------------------------------------------

  ! We follow the PVRS solver from Toro (Chapter 9)
  !rhobar = HALF*(rho_l + rho_r)
  !csbar = HALF*(cs_l + cs_r)

  !pstar = HALF*(p_l + p_r) + HALF*(u_l - u_r)*rhobar*csbar

  ! alternative: two shock solver (see Toro 9.42)
  W_l = rho_l*cs_l
  W_r = rho_r*cs_r

  pstar = ((W_r*p_l + W_l*p_r) + W_l*W_r*(u_l - u_r))/(W_l + W_r)

  pstar = max(pstar, smallp)
  

  !---------------------------------------------------------------------------
  ! find the exact pstar and ustar
  !---------------------------------------------------------------------------

  ! this procedure follows directly from Colella & Glaz 1985, section 1

  converged = .false.
  iter = 1
  do while (.not. converged .and. iter < max_iter)

     ! compute Z_l and Z_r -- the form of these depend on whether the
     ! wave is a shock or a rarefaction
     
     ! left wave
     if (pstar > p_l) then
        ! left shock
        call shock(pstar, rho_l, u_l, p_l, xn, gammaE_bar, gammaC_bar, Z_l, W_l)
        lwave = "shock"
     else
        ! left rarefaction
        call rarefaction(pstar, rho_l, u_l, p_l, xn, 1, Z_l, W_l)
        lwave = "rarefaction"
     endif

     ! right wave
     if (pstar > p_r) then
        ! right shock
        call shock(pstar, rho_r, u_r, p_r, xn, gammaE_bar, gammaC_bar, Z_r, W_r)
        rwave = "shock"
     else
        ! right rarefaction
        call rarefaction(pstar, rho_r, u_r, p_r, xn, 3, Z_r, W_r)
        rwave = "rarefaction"
     endif

     print *, 'left wave: ', trim(lwave), ';   rightwave: ', trim(rwave)

     ustar_l = u_l - (pstar - p_l)/W_l
     ustar_r = u_r + (pstar - p_r)/W_r

     pstar_new = pstar - Z_l*Z_r*(ustar_r - ustar_l)/(Z_l + Z_r)

     print *, "done with iteration", iter
     print *, "ustar_l/r, pstar: ", &
          real(ustar_l), real(ustar_r), real(pstar_new)
     print * , " "

     ! estimate the error in the current star solution
     err1 = abs(ustar_r - ustar_l)
     err2 = pstar_new - pstar

     print *, "ERRORS: ", err1, err2
     print *, " "

     if (err1 < tol*max(abs(ustar_l),abs(ustar_r)) .and. err2 < tol*pstar) then
        converged = .true.
     endif

     ! get ready for the next iteration
     pstar = pstar_new

     iter = iter + 1
     
  enddo

  ustar = HALF*(ustar_l + ustar_r)

  print *, 'found pstar, ustar: ', pstar, ustar

  ! let's test if our integration across the rarefaction works, as expected
  if (lwave == "rarefaction") then
     call rarefaction(pstar, rho_l, u_l, p_l, xn, 1, Z_temp, W_temp, rhostar)

     print *, "here"

     ! get the soundspeed via the EOS (C&G suggest getting it from
     ! the jump conditions)
     eos_state%rho = rhostar
     eos_state%p = pstar
     eos_state%xn(:) = xn(:)
     eos_state%T = 100000.0_dp_t   ! initial guess
     
     call eos(eos_input_rp, eos_state, .false.)
        
     cs_star = sqrt(eos_state%gam1*pstar/rhostar)

     print *, "about to test"

     call rarefaction_to_u(rho_l, u_l, p_l, xn, 1, ustar-cs_star, rho, p, u)

     print *, "left rarefaction integration test"
     print *, "(this should agree with pstar/ustar if the 1-wave is a rarefaction): ", p, u

  endif

  if (rwave == "rarefaction") then
     call rarefaction(pstar, rho_r, u_r, p_r, xn, 3, Z_temp, W_temp, rhostar)

     ! get the soundspeed via the EOS (C&G suggest getting it from
     ! the jump conditions)
     eos_state%rho = rhostar
     eos_state%p = pstar
     eos_state%xn(:) = xn(:)
     eos_state%T = 100000.0_dp_t   ! initial guess
     
     call eos(eos_input_rp, eos_state, .false.)
        
     cs_star = sqrt(eos_state%gam1*pstar/rhostar)

     call rarefaction_to_u(rho_r, u_r, p_r, xn, 3, ustar+cs_star, rho, p, u)

     print *, "left rarefaction integration test"
     print *, "(this should agree with pstar/ustar if the 3-wave is a rarefaction): ", p, u

  endif



  !---------------------------------------------------------------------------
  ! find the solution as a function of xi = x/t
  !---------------------------------------------------------------------------

  open(newunit=lun, file='riemann.out', status="unknown", action="write")

  ! This follows from the discussion around C&G Eq. 15

  dx = (xmax - xmin)/npts

  ! loop over xi space
  do i = 1, npts

     ! compute xi = x/t -- this is the similarity variable for the
     ! solution
     x  = xmin + (dble(i) - HALF)*dx
     xi = (x - xjump)/t

     ! check which side of the contact we need to worry about
     chi = sign(ONE, xi - ustar)

     if (chi == -ONE) then
        rho_s = rho_l
        u_s = u_l
        p_s = p_l

        W_s = W_l

        cs_s = cs_l

        uhat_s = chi*u_s
        xihat = chi*xi
        uhat_star = chi*ustar

     else if (chi == ONE) then
        rho_s = rho_r
        u_s = u_r
        p_s = p_r

        W_s = W_r

        cs_s = cs_r

        uhat_s = chi*u_s
        xihat = chi*xi
        uhat_star = chi*ustar

     else
        ! we should average in this case
        call bl_error("Not implemented")

     endif

     ! are we a shock or rarefaction?
     if (pstar > p_s) then
        ! shock
        
        rhostar = ONE/(ONE/rho_s - (pstar - p_s)/W_s**2)

        !print *, 'shock: rhostar = ', rhostar
     else
        ! rarefaction.  Here we need to integrate the Riemann
        ! invariant curves to our pstar to get rhostar and cs_star
        if (chi == -ONE) then
           call rarefaction(pstar, rho_s, u_s, p_s, xn, 1, Z_temp, W_temp, rhostar)
        else
           call rarefaction(pstar, rho_s, u_s, p_s, xn, 3, Z_temp, W_temp, rhostar)
        endif

        !print *, 'rare: rhostar = ', rhostar
     endif

     ! get the soundspeed via the EOS (C&G suggest getting it from
     ! the jump conditions)
     eos_state%rho = rhostar
     eos_state%p = pstar
     eos_state%xn(:) = xn(:)
     eos_state%T = 100000.0_dp_t   ! initial guess

     call eos(eos_input_rp, eos_state, .false.)
        
     cs_star = sqrt(eos_state%gam1*pstar/rhostar)

     ! now deal with the cases where we are not spanning a rarefaction
     if (pstar <= p_s) then
        lambdahat_s = uhat_s + cs_s
        lambdahat_star = uhat_star + cs_star
     else
        lambdahat_s = uhat_s + W_s/rho_s
        lambdahat_star = lambdahat_s
     endif

     if (xihat <= lambdahat_star) then
        p = pstar
        rho = rhostar
        u = ustar

     else if (xihat > lambdahat_s) then
        p = p_s
        rho = rho_s
        u = u_s

     else
        ! we are inside the rarefaction.  To find the solution here,
        ! we need to integrate up to the point where uhat + c = xihat
        ! starting from U = U_s
        
        ! for the 1-rarefaction, chi = -1, so this is -u + c = -xi,
        ! or u - c = xi, meaning we integrate to u = xi + c
        
        ! for the 3-rarefaction, chi = 1, so this is u + c = xi,
        ! so we integrate to u = xi - c

        ! Note that c here is c(rho,p) -- we need to compute that
        ! self-consistently as we integrate
        if (chi == -ONE) then
           call rarefaction_to_u(rho_s, u_s, p_s, xn, 1, xi, rho, p, u)

        else if (chi == ONE) then
           call rarefaction_to_u(rho_s, u_s, p_s, xn, 3, xi, rho, p, u)

        endif
     endif

     ! get the thermodynamics for this state for output
     eos_state%rho = rho
     eos_state%p = p
     eos_state%xn(:) = xn(:)
     eos_state%T = 100000.0_dp_t   ! initial guess

     call eos(eos_input_rp, eos_state, .false.)
        
     if (i == 1) then
        write (unit=lun, fmt="(a1, a3, 8(1x, a20))") "#", "i", "x", "rho", "u", "p", "T", "e", "gamma_1"
     endif

     write (unit=lun, fmt="(1x, i3, 8(1x, g20.10))") i, x, rho, u, p, eos_state%T, eos_state%e, eos_state%gam1

  enddo

  close (unit=lun)

  call runtime_close()

end program riemann_exact

