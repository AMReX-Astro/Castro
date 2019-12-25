! Compute the analytic solution for the relaxation to thermal
! equilibrium radiation problem.  We follow the problem outline from
! Swesty and Myra (2009), section 4.4.1.
!
! This code returns the radiation energy as a function of time.


module constants_module

  use amrex_fort_module, only : rt => amrex_real

  ! fundamental constants
  real(rt), parameter :: k_B      = 1.3806504e-16_rt   ! erg / K
  real(rt), parameter :: c_light  = 2.99792458e10_rt   ! cm / s
  real(rt), parameter :: sigma_SB = 5.670400e-5_rt     ! erg/s/cm^2/K^4
  real(rt), parameter :: a_rad    = 4.0e0_rt*sigma_SB/c_light
  real(rt), parameter :: m_p      = 1.672621637e-24_rt ! g

  ! problem parameters

  ! matter density
  real(rt), parameter :: rho_ambient = 1.e-7   ! g / cm**3

  ! mean molecular weight
  real(rt), parameter :: mu = 0.6

  ! ratio of specific heats
  real(rt), parameter :: gamma = 5.0e0_rt/3.0e0_rt

  ! starting time
  real(rt), parameter :: t_0 = 1.e-16_rt   ! s

  ! starting radiation energy density
  real(rt), parameter :: E_rad = 1.e12_rt  ! erg / cm**3

  ! starting matter energy density
  !
  ! This is one of the main parameters that controls the behavior of
  ! the problem.  The matter energy sets the matter temperature (via
  ! the EOS).  If (a T**4) > E_rad, then we are doing a cooling
  ! problem.  If (a T**4) < E_rad, then we are doing a heating problem
  ! (i.e. the matter temperature will increase with time to the
  ! equilibrium temperature).

  !real(rt), parameter :: rhoe_i = 1.e9_rt  ! erg / cm**3
  real(rt), parameter :: rhoe_i = 1.e2_rt  ! erg / cm**3

  ! we will specify the opacity as a constant kappa.  In Castro lingo,
  ! this is the Planck mean opacity, which appears in the source term.
  real(rt), parameter :: kappa_p = 4.e-8_rt


end module constants_module

function f(x,beta,eta) result (fval)

  real(rt), intent(in)  :: x, beta, eta
  real(rt) :: fval

  ! There are different analytic solutions depending on whether we are
  ! heating (x < beta/eta) or cooling (x > beta/eta)

  if (x < beta/eta) then

     ! note: this form is different than what appears in Swesty and Myra,
     ! but instead, this is the form that Doug seems to have used in his
     ! code.
     fval = (1.e0_rt/(2.e0_rt*eta*beta**3)) * &
          (0.5e0_rt*log((beta + eta*x)/(beta - eta*x)) + &
          atan(eta*x/beta))
  else
     fval = (1.e0_rt/(2.e0_rt*eta*beta**3)) * &
          (0.5e0_rt*log((x + beta/eta)/(x - beta/eta)) + &
          atan(eta*x/beta))
  endif

  return
end function f

program analytic

  ! compute the analytic solution to the thermal relaxation problem,
  ! as given by Swesty and Myra (2009), Eq. 83, 84.

  ! The analytic solution provides t(rho*e) for the gas component.

  use constants_module

  implicit none

  real(rt) :: safe_print
  real(rt) :: f

  real(rt) :: beta, c_v, eta
  real(rt) :: rhoe_f, rhoe, dlog_rhoe
  real(rt) :: t

  integer :: n
  integer, parameter :: nsteps = 2000

  ! the beta constant as implied in Swesty and Myra, eq. 81
  beta = (c_light*kappa_p*E_rad)**0.25e0_rt

  ! specific heat at constant volume for an idea gas (Note: Swesty and
  ! Myra work with temperature in MeV, so Boltzmann's constant, k_B,
  ! does not appear in their derivation.
  c_v = k_B/(mu*m_p*(gamma - 1.e0_rt))

  ! the eta constant as implied in Swesty and Myra, eq. 81
  eta = ((c_light*kappa_p*a_rad)**0.25e0_rt)/(c_v*rho_ambient)

  ! final energy density (see Swesty and Myra, text just before
  ! Eq. 82)
  rhoe_f = beta/eta

  ! energy increment
  dlog_rhoe = (log(rhoe_f) - log(rhoe_i))/(nsteps-1)


  print *, '# t, rho*e (matter)'

  ! loop in energy to find corresponding time
  do n = 1, nsteps
     rhoe = exp(log(rhoe_i) + dble(n-1)*dlog_rhoe)

     t = t_0 + f(rhoe,beta,eta) - f(rhoe_i,beta,eta)

     print *, t, rhoe
  enddo

end program analytic
