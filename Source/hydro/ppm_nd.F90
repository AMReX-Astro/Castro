module ppm_module
  use amrex_constants_module, only: ZERO, HALF, ONE, TWO, SIXTH, &
                                    TWO3RD, THREE, SIX, SEVEN12TH, TWELFTH
  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ppm_reconstruct_f90(s, flatn, sm, sp)
    ! This routine does the reconstruction of the zone data into a parabola.

    implicit none

    real(rt), intent(in   ) :: s(-2:2), flatn
    real(rt), intent(inout) :: sm, sp

    ! local
    real(rt) :: dsl, dsr, dsc, dsvl_l, dsvl_r

    !$gpu

    ! Compute van Leer slopes

    dsl = TWO * (s(-1) - s(-2))
    dsr = TWO * (s(0) - s(-1))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(0) - s(-2))
       dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
    else
       dsvl_l = ZERO
    end if

    dsl = TWO * (s(0) - s(-1))
    dsr = TWO * (s(1) - s(0))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(1) - s(-1))
       dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
    else
       dsvl_r = ZERO
    end if

    ! Interpolate s to edges

    sm = HALF*(s(0) + s(-1)) - SIXTH*(dsvl_r - dsvl_l)

    ! Make sure sedge lies in between adjacent cell-centered values

    sm = max(sm, min(s(0), s(-1)))
    sm = min(sm, max(s(0), s(-1)))



    ! Compute van Leer slopes

    dsl = TWO  * (s(0) - s(-1))
    dsr = TWO  * (s(1) - s(0))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(1) - s(-1))
       dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
    else
       dsvl_l = ZERO
    end if

    dsl = TWO  * (s(1) - s(0))
    dsr = TWO  * (s(2) - s(1))
    if (dsl*dsr .gt. ZERO) then
       dsc = HALF * (s(2) - s(0))
       dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
    else
       dsvl_r = ZERO
    end if

    ! Interpolate s to edges

    sp = HALF*(s(1) + s(0)) - SIXTH*(dsvl_r - dsvl_l)

    ! Make sure sedge lies in between adjacent cell-centered values

    sp = max(sp, min(s(1), s(0)))
    sp = min(sp, max(s(1), s(0)))



    ! Flatten the parabola

    sm = flatn * sm + (ONE - flatn) * s(0)
    sp = flatn * sp + (ONE - flatn) * s(0)

    ! Modify using quadratic limiters -- note this version of the limiting comes
    ! from Colella and Sekora (2008), not the original PPM paper.

    if ((sp - s(0)) * (s(0) - sm) .le. ZERO) then

       sp = s(0)
       sm = s(0)

    else if (abs(sp - s(0)) .ge. TWO * abs(sm - s(0))) then

       sp = THREE * s(0) - TWO * sm

    else if (abs(sm - s(0)) .ge. TWO * abs(sp - s(0))) then

       sm = THREE * s(0) - TWO * sp

    end if

  end subroutine ppm_reconstruct_f90



  subroutine ppm_int_profile_f90(sm, sp, sc, u, c, dtdx, Ip, Im)
    ! Integrate the parabolic profile to the edge of the cell.

    implicit none

    real(rt), intent(in   ) :: sm, sp, sc, u, c, dtdx
    real(rt), intent(inout) :: Ip(1:3), Im(1:3)

    ! local
    real(rt) :: speed, sigma, s6

    !$gpu

    ! compute x-component of Ip and Im
    s6 = SIX * sc - THREE * (sm + sp)

    ! Ip/m is the integral under the parabola for the extent
    ! that a wave can travel over a timestep
    !
    ! Ip integrates to the right edge of a cell
    ! Im integrates to the left edge of a cell

    ! u-c wave
    speed = u - c
    sigma = abs(speed) * dtdx

    ! if speed == ZERO, then either branch is the same
    if (speed <= ZERO) then
       Ip(1) = sp
       Im(1) = sm + HALF*sigma * (sp - sm + (ONE - TWO3RD * sigma) * s6)
    else
       Ip(1) = sp - HALF*sigma * (sp - sm - (ONE - TWO3RD * sigma) * s6)
       Im(1) = sm
    endif

    ! u wave
    speed = u
    sigma = abs(speed) * dtdx

    if (speed <= ZERO) then
       Ip(2) = sp
       Im(2) = sm + HALF * sigma * (sp - sm + (ONE - TWO3RD * sigma) * s6)
    else
       Ip(2) = sp - HALF * sigma * (sp - sm - (ONE - TWO3RD * sigma) * s6)
       Im(2) = sm
    endif

    ! u+c wave
    speed = u + c
    sigma = abs(speed) * dtdx

    if (speed <= ZERO) then
       Ip(3) = sp
       Im(3) = sm + HALF * sigma * (sp - sm + (ONE - TWO3RD * sigma) * s6)
    else
       Ip(3) = sp - HALF * sigma * (sp - sm - (ONE - TWO3RD * sigma) * s6)
       Im(3) = sm
    endif

  end subroutine ppm_int_profile_f90

end module ppm_module
