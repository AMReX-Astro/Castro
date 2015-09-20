program riemann_exact

  use bl_types
  use bl_constants_module
  use eos_module
  use eos_type_module
  use network, only: nspec
  use probin_module, only: rho_l, u_l, p_l, T_l, rho_r, u_r, p_r, T_r, &
                           xmin, xmax, xjump, t, npts, use_Tinit, &
                           co_moving_frame
  use riemann_support
  use runtime_init_module
  use riemann_sample_module
  use riemann_star_module

  implicit none

  real (kind=dp_t) :: xn(nspec)

  real (kind=dp_t) :: pstar, ustar
  real (kind=dp_t) :: W_avg

  real (kind=dp_t) :: W_l, W_r

  type (eos_t) :: eos_state

  real (kind=dp_t) :: x, dx
  real (kind=dp_t) :: rho, u, p, xn_s(nspec)

  integer :: i, lun

  ! general Maestro initializations
  call runtime_init()

  ! microphysics
  call network_init()
  call eos_init()  !gamma_in=1.4d0)

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

  if (co_moving_frame) then
     W_avg = 0.5d0*(u_l + u_r)
     u_l = u_l - W_avg
     u_r = u_r - W_avg
  endif

  call riemann_star_state(rho_l, u_l, p_l, rho_r, u_r, p_r, xn, xn, &
                          ustar, pstar, W_l, W_r, .true.)

  !if (co_moving_frame) then
  !   u_l = u_l + W_avg
  !   u_r = u_r + W_avg
  !   ustar = ustar + W_avg
  !endif

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
     !if (co_moving_frame) then
     !   x = x + W_avg*t
     !endif

     call riemann_sample(rho_l, u_l, p_l, rho_r, u_r, p_r, xn, xn, &
                         ustar, pstar, W_l, W_r, &
                         x, xjump, t, &
                         rho, u, p, xn_s)

     if (co_moving_frame) then
        u = u + W_avg
        x = x + t*W_avg
     endif

     ! get the thermodynamics for this state for output
     eos_state%rho = rho
     eos_state%p = p
     eos_state%xn(:) = xn_s(:)
     eos_state%T = 100000.0_dp_t   ! initial guess

     call eos(eos_input_rp, eos_state, .false.)
        
     if (i == 1) then
        write (unit=lun, fmt="(a1, a3, 8(1x, a25))") &
             "#", "i", "x", "rho", "u", "p", "T", "e", "gamma_1"
     endif

     write (unit=lun, fmt="(1x, i3, 8(1x, g25.15))") &
          i, x, rho, u, p, eos_state%T, eos_state%e, eos_state%gam1

  enddo

  close (unit=lun)

  call runtime_close()

end program riemann_exact

