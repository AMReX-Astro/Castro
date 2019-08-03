module riemann_sample_module

contains

subroutine riemann_sample(rho_l, u_l, p_l, &
                          rho_r, u_r, p_r, &
                          xn_l, xn_r, &
                          ustar, pstar, &
                          W_l, W_r, &
                          x, xjump, time, &
                          rho, u, p, xn)

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use eos_module
  use eos_type_module
  use network, only: nspec
  use riemann_support

  implicit none

  real (rt), intent(in) :: rho_l, u_l, p_l
  real (rt), intent(in) :: rho_r, u_r, p_r
  real (rt), intent(in) :: xn_l(nspec), xn_r(nspec)
  real (rt), intent(in) :: ustar, pstar
  real (rt), intent(in) :: W_l, W_r
  real (rt), intent(in) :: x, xjump, time
  real (rt), intent(out) :: rho, u, p, xn(nspec)

  real (rt) :: cs_l, cs_r

  real (rt) :: W_temp, Z_temp

  type (eos_t) :: eos_state

  real (rt) :: xi, chi
  real (rt) :: rhostar
  real (rt) :: xihat, uhat_s, uhat_star
  real (rt) :: lambdahat_s, lambdahat_star
  real (rt) :: cs_star
  real (rt) :: rho_s, p_s, u_s, W_s, cs_s

  ! get the initial sound speeds
  eos_state%rho = rho_l
  eos_state%p = p_l
  eos_state%xn(:) = xn_l(:)
  eos_state%T = initial_temp_guess

  call eos(eos_input_rp, eos_state)

  cs_l = sqrt(eos_state%gam1*p_l/rho_l)

  eos_state%rho = rho_r
  eos_state%p = p_r
  eos_state%xn(:) = xn_r(:)
  eos_state%T = initial_temp_guess

  call eos(eos_input_rp, eos_state)

  cs_r = sqrt(eos_state%gam1*p_r/rho_r)


  !---------------------------------------------------------------------------
  ! find the solution as a function of xi = x/t
  !---------------------------------------------------------------------------

  ! This follows from the discussion around C&G Eq. 15

  ! compute xi = x/t -- this is the similarity variable for the
  ! solution
  xi = (x - xjump)/time

  ! check which side of the contact we need to worry about
  chi = sign(ONE, xi - ustar)

  if (chi == -ONE) then
     rho_s = rho_l
     u_s = u_l
     p_s = p_l
     
     W_s = W_l

     cs_s = cs_l

     xn = xn_l

     uhat_s = chi*u_s
     xihat = chi*xi
     uhat_star = chi*ustar

  else if (chi == ONE) then
     rho_s = rho_r
     u_s = u_r
     p_s = p_r

     W_s = W_r

     cs_s = cs_r

     xn = xn_r

     uhat_s = chi*u_s
     xihat = chi*xi
     uhat_star = chi*ustar

  else
     ! we should average in this case
     call castro_error("Not implemented")

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

  endif

  ! get the soundspeed via the EOS (C&G suggest getting it from
  ! the jump conditions)
  eos_state%rho = rhostar
  eos_state%p = pstar
  eos_state%xn(:) = xn(:)
  eos_state%T = initial_temp_guess

  call eos(eos_input_rp, eos_state)
        
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

end subroutine riemann_sample

end module riemann_sample_module
