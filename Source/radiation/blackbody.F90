module blackbody_module

  use fundamental_constants_module, only: a_rad, k_B, hplanck
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter :: pi = 3.1415926535897932384626e0_rt
  real(rt), parameter :: bk_const = a_rad * 15.e0_rt/pi**4
  real(rt), parameter :: magic_const = pi**4/15.e0_rt
  integer, parameter :: Ncoefs = 8
  real(rt), dimension(0:Ncoefs) :: coefs = &
                                           [1.e0_rt/60.e0_rt, &
                                           -1.e0_rt/5040.e0_rt, &
                                            1.e0_rt/272160.e0_rt, &
                                           -1.e0_rt/13305600.e0_rt, &
                                            1.e0_rt/622702080.e0_rt, &
                                           -691.e0_rt/19615115520000.e0_rt, &
                                            1.e0_rt/1270312243200.e0_rt, &
                                           -3617.e0_rt/202741834014720000.e0_rt, &
                                            43867.e0_rt/107290978560589824000.e0_rt]
  real(rt), parameter :: tol = 1.0e-10_rt

  real(rt), parameter :: xmagic = 2.061981e0_rt
  real(rt), parameter :: xsmall = 1.e-5_rt
  real(rt), parameter :: xlarge = 100.e0_rt

  private :: pi, bk_const, magic_const, Ncoefs, coefs, tol, xmagic, xsmall, xlarge

contains

  subroutine BdBdTIndefInteg(T, nu, B, dBdT)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in   ) :: T, nu
    real(rt), intent(  out) :: B, dBdT

    real(rt) :: x, integ, part

    ! This routine evaluates the incomplete integral of the Planck function:
    ! integ(x) = int_{0}^{x} b(x') dx'
    !
    ! We follow B. A. Clark, "Computing Multigroup Radiation Integrals using
    ! Polylogarithm-Based Methods," JCP (1987), 70, 311.
    !
    ! This incomplete integral can be used to recursively obtain the blackbody
    ! function for all groups using Equation 35:
    ! b_g = integ(x_{g+1}) - integ(x_{g})
    !
    ! The incomplete integral is defined by Equation 25 in terms of polylogarithm
    ! functions evaluated at exp(-x) where x = h * nu / k_B T. As Clark notes,
    ! we evaluate those polylogarithm functions directly when x is sufficiently
    ! large, approximating them with a truncated infinite series; but, when x is
    ! sufficiently small, we instead evaluate a series expansion of the original
    ! integral. The threshold between these two approaches is set so that there
    ! is no discontinuity of the Planck function at the threshold.

    ! Specifically, we follow Sandia Report SAND2005-6988:
    ! Advances in Radiation Modeling in ALEGRA: A Final Report for LDRD-67120,
    ! Efficient Implicit Multigroup Radiation Calculations
    ! T. A. Brunner, T. Mehlhorn, R. McClarren, & C. J. Kurecka (2005)
    !
    ! They set the threshold at the value 2.061981, above which we evaluate the
    ! polylogarithm functions and below which we use the series expansion.
    !
    ! We can also use this method to evaluate the derivative of the Planck function
    ! with respect to temperature.

    x = hplanck * nu / (k_B * T)

    if (x .gt. xlarge) then

       ! If x is sufficiently large, this is effectively an integral over
       ! all frequencies, and we know analytically that this is Stefan-Boltzmann
       ! radiation, so we can save time (and avoid difficulties with evaluating
       ! the exponential numerically for large argument) by evaluating it directly.

       B = a_rad * T**4
       dBdT = 4.0_rt * a_rad * T**3

    else if (x .lt. xsmall) then

       ! When x is very small, the integral is effectively over zero frequency
       ! space and we can approximate it to zero.

       B = 0.0_rt
       dBdT = 0.0_rt

    else

       if (x .gt. xmagic) then
          integ = integlarge(x)
       else
          integ = integsmall(x)
       end if

       ! Clark, Equation 3

       B = bk_const * T**4 * integ

       ! ALEGRA, Equation 2.1.66
       ! Note that the second term in square brackets is the one
       ! we care about for this group, and the third term is from
       ! the recursive relationship between groups. Since we are
       ! integrating upwards from zero frequency, the sign is reversed.

       part = x**4 / (exp(x) - 1.0_rt)
       dBdT = bk_const * T**3 * (4.0_rt * integ - part)

    end if

  end subroutine BdBdTIndefInteg



  function BIndefInteg(T, nu)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: T, nu

    real(rt) :: BIndefInteg
    real(rt) :: x, integ

    ! See comments above; this function evaluates B but does not
    ! evaluate dB/dT.

    x = hplanck * nu / (k_B * T)

    if (x .gt. xlarge) then

       BIndefInteg = a_rad * T**4

    else if (x .lt. xsmall) then

       BIndefInteg = 0.e0_rt

    else

       if (x .gt. xmagic) then
          integ = integlarge(x)
       else
          integ = integsmall(x)
       end if

       BIndefInteg = bk_const * T**4 * integ

    end if

  end function BIndefInteg



  function BGroup(T, nu0, nu1)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: T, nu0, nu1

    real(rt) :: BGroup

    ! Clark, Equation 35

    BGroup = BIndefInteg(T, nu1) - BIndefInteg(T, nu0)

  end function BGroup



  function Li(n, z)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in) :: n
    real(rt), intent(in) :: z

    real(rt) :: Li, t
    integer :: k
    integer, parameter :: kmax = 18

    ! Clark, Equation 19 / ALEGRA, Equation 2.1.64
    !
    ! z here is really exp(-h nu / kT), and nu and T are positive,
    ! so it satisfies the range of convergence.

    Li = z

    do k = 2, kmax

       t = z**k / k**n
       Li = Li + t

       ! Terminate the series when the additional terms
       ! become small enough that they approach roundoff error.

       if (t / Li .lt. tol) then
          exit
       end if

    end do

  end function Li



  function integlarge(x)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: x

    real(rt) :: integlarge, z

    ! ALEGRA, Equation 2.1.63
    !
    ! Note that since we define z == exp(-x), then
    ! x = -ln(z), so the signs below are consistent
    ! with that relative to the paper (i.e. the terms
    ! with odd powers of x have reversed sign).
    !
    ! We are only evaluating this at a specific frequency,
    ! and relying on the recursion relation to give the
    ! correct absolute result.

    z = exp(-x)

    integlarge = magic_const &
                 - (x**3 * Li(1,z) + 3.e0_rt * x**2 * Li(2,z) &
                 + 6.e0_rt* x * Li(3,z) + 6.e0_rt * Li(4,z))

  end function integlarge


  function integsmall(x)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: x
    real(rt) :: integsmall, t
    integer :: i
    real(rt) :: x2, x3, xfoo

    ! ALEGRA, Equation 2.1.65

    x2 = x**2
    x3 = x**3

    integsmall = x3 / 3.0_rt - x2**2 / 8.0_rt

    xfoo = x3

    do i = 0, Ncoefs

       xfoo = xfoo * x2
       t = coefs(i) * xfoo
       integsmall = integsmall + t

       ! Terminate the series expansion early if the terms
       ! approach numerical roundoff error.

       if (abs(t / integsmall) .lt. tol) then
          exit
       end if

    end do

  end function integsmall

end module blackbody_module
