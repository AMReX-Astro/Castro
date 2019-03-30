module binary_module

  ! This module contains a number of routines that are designed to
  ! calculate generic properties of binary orbits.
  
  implicit none

contains

  ! Given the mass ratio q of two stars (assumed to be q = M_1 / M_2), 
  ! compute the effective Roche radii of the stars, normalized to unity, 
  ! using the approximate formula of Eggleton (1983). Optionally we can
  ! pass in a distance scale.
  
  subroutine get_roche_radii(mass_ratio, r_1, r_2, a)

    use amrex_constants_module, only: ONE, TWO3RD, THIRD

    implicit none

    double precision, intent(in   ) :: mass_ratio
    double precision, intent(inout) :: r_1, r_2
    double precision, intent(in   ), optional :: a

    double precision :: q
    double precision :: c1, c2

    double precision :: scale

    if (present(a)) then
       scale = a
    else
       scale = ONE
    endif

    c1 = 0.49d0
    c2 = 0.60d0

    q = mass_ratio

    r_1 = scale * c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

    q = ONE / q

    r_2 = scale * c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

  end subroutine get_roche_radii


  
  ! Calculate Lagrange points. In each case we give the zone index
  ! closest to it (assuming we're on the coarse grid).

  subroutine get_lagrange_points(mass_1, mass_2, com_1, com_2, &
                                 L1, L2, L3) bind(C, name="get_lagrange_points")

    use amrex_constants_module

    implicit none

    double precision, intent(in   ), value :: mass_1, mass_2
    double precision, intent(in   ) :: com_1(3), com_2(3)
    double precision, intent(inout) :: L1(3), L2(3), L3(3)
    
    double precision :: r ! Distance from Lagrange point to primary
    double precision :: a ! Distance between secondary and primary

    double precision :: r1, r2

    ! Don't try to calculate the Lagrange points if the secondary
    ! is already gone.

    if (mass_2 < ZERO) return

    a = sqrt(sum((com_2 - com_1)**2))

    r1 = -sqrt(sum(com_1**2))
    r2 = sqrt(sum(com_2**2))

    ! Do a root-find over the quintic equation for L1. 

    r = ZERO

    call lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_min = r1, r_max = r2)
    
    ! Now turn this radial distance into a grid coordinate.

    L1 = r * (com_2 - com_1) / a

    ! Repeat for L2 Lagrange point.

    r = r2 + HALF * a

    call lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_min = r2)
    
    L2 = r * (com_2 - com_1) / a

    ! Repeat for L3 Lagrange point.

    r = r1 - HALF * a

    call lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_max = r1)
    
    L3 = r * (com_2 - com_1) / a
    
  end subroutine get_lagrange_points



  ! Iterate over the force balance equation to find a Lagrange point.
  ! r_min and r_max set the domain over which we want to find the answer,
  ! since in general the equation has multiple roots.
  ! We assume that r comes in with a valid starting guess value.
  
  subroutine lagrange_iterate(r, mass_1, mass_2, r1, r2, a, r_min, r_max)

    use amrex_constants_module, only: ZERO, HALF
    use amrex_error_module, only: amrex_error
    
    implicit none

    ! Input variables
    
    double precision :: r
    double precision :: mass_1, mass_2, r1, r2, a
    double precision, optional :: r_min, r_max
    
    ! Root-find parameters
    
    double precision :: tolerance = 1.0d-8
    integer          :: max_iters = 200

    ! Local variables
    
    double precision :: rm, rp, rc, width, fm, fp, fc
    integer :: i

    if (.not. (present(r_min) .or. present(r_max))) then
       call amrex_error("Lagrange point iteration must have at least one bound provided.")
    else if (present(r_min) .and. present(r_max)) then
       rm = r_min
       rp = r_max
    else if (present(r_min) .and. (.not. present(r_max))) then
       rm = r_min
       rp = abs(r_min) * 1000.0d0
    else if (present(r_max) .and. (.not. present(r_min))) then
       rm = -abs(r_max) * 1000.0d0
       rp = r_max
    endif

    width = (rp - rm)
    
    rm = rm + width / 1000.0d0
    rp = rp - width / 1000.0d0

    ! Use a bisection search to find the root of the force-balance equation.
    ! The reason we don't use something faster like Newton-Raphson is that
    ! the force terms have terms like 1/r^2, so evaluating the derivative and
    ! then dividing it can be numerically dangerous. Since we only call this routine
    ! once per timestep and the force evaluation is cheap, it shouldn't matter
    ! that the method itself converges slowly.
    
    do i = 1, max_iters

       fm = fL(mass_1, mass_2, r1, r2, rm, a)
       fp = fL(mass_1, mass_2, r1, r2, rp, a)       
       
       rc = HALF * (rm + rp)
       fc = fL(mass_1, mass_2, r1, r2, rc, a)              

       if (fm * fc < ZERO) then
          rp = rc
       else
          rm = rc
       endif

       ! Check to see if we've convernged.

       if (abs( (rp - rm) / a ) < tolerance) then
          exit
       endif

    enddo

    r = rc

    ! If we didn't converge, set the value to zero.
    
    if (i > max_iters) then
       r = ZERO
    endif    
    
  end subroutine lagrange_iterate
  


  function fL(M1, M2, r1, r2, r, a)

    use amrex_constants_module

    implicit none

    double precision :: M1, M2, r1, r2, r, a
    double precision :: fL

    double precision :: g1, g2, c

    g1 = gforce(M1, r - r1)
    g2 = gforce(M2, r - r2)

    c  = cforce(M1 + M2, a, r)
    
    fL = g1 + g2 + c

  end function fL



  function fdLdr(M1, M2, r1, r2, r, a)

    use amrex_constants_module

    implicit none

    double precision :: M1, M2, r1, r2, r, a
    double precision :: fdLdr

    double precision :: dg1, dg2, dc
    
    dg1 = dgforcedr(M1, r - r1)
    dg2 = dgforcedr(M2, r - r2)

    dc  = dcforcedr(M1 + M2, a, r)
    
    fdLdr = dg1 + dg2 + dc

  end function fdLdr



  function gforce(M, r)

    use amrex_constants_module, only: ONE    
    use fundamental_constants_module, only: Gconst

    double precision :: M, r
    
    double precision :: gforce

    gforce = -Gconst * M / r**2 * sign(ONE,r)

  end function gforce


  function dgforcedr(M, r)

    use amrex_constants_module, only: ONE, TWO
    use fundamental_constants_module, only: Gconst

    double precision :: M, r
    
    double precision :: dgforcedr

    dgforcedr = TWO * Gconst * M / r**3 * sign(ONE,r)

  end function dgforcedr


  
  function cforce(M, a, r)

    use amrex_constants_module, only: ONE    
    use fundamental_constants_module, only: Gconst

    double precision :: M, a, r
    
    double precision :: cforce

    cforce = Gconst * M / a**3 * r

  end function cforce



  function dcforcedr(M, a, r)

    use amrex_constants_module, only: ONE
    use fundamental_constants_module, only: Gconst

    double precision :: M, a, r    
    
    double precision :: dcforcedr

    dcforcedr = Gconst * M / a**3

  end function dcforcedr
  
end module binary_module
