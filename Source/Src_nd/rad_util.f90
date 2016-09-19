module rad_util_module

  use bl_types
  use bl_error_module
  use rad_params_module, only : ngroups
  use bl_constants_module

  implicit none

contains

  subroutine compute_ptot_ctot(lam, q, cg, ptot, ctot, gamc_tot)

    use meth_params_module, only : QPRES, QRHO
    use radhydro_params_module, only: comoving, QRAD, QPTOT, QRADVAR
    use fluxlimiter_module, only : Edd_factor

    real (kind=dp_t), intent(in) :: lam(0:ngroups-1)
    real (kind=dp_t), intent(in) :: q(QRADVAR)
    real (kind=dp_t), intent(in) :: cg
    real (kind=dp_t), intent(out) :: ptot
    real (kind=dp_t), intent(out) :: ctot
    real (kind=dp_t), intent(out) :: gamc_tot

    integer :: g

    real (kind=dp_t) :: csrad2, Eddf, gamr, prad

    csrad2 = ZERO
    prad = ZERO

    do g = 0, ngroups-1
       if (comoving) then
          Eddf = Edd_factor(lam(g))
          gamr = (THREE - Eddf)/TWO
       else
          gamr = lam(g) + ONE
       end if

       prad = prad + lam(g)*q(QRAD+g)
       csrad2 = csrad2 + gamr * (lam(g)*q(QRAD+g)) / q(QRHO)
    end do

    ptot = q(QPRES) + prad

    ctot = cg**2 + csrad2
    gamc_tot = ctot * q(QRHO) / ptot

    ctot = sqrt(ctot)

  end subroutine compute_ptot_ctot

  function FLDlambda(r, limiter) result (lambda)

    double precision :: r
    integer :: limiter

    double precision :: lambda

    if (limiter .eq. 0) then
       ! no limiter
       lambda = 1.d0/3.d0

    else if (limiter < 10) then
       ! approximate LP
       lambda = (2.d0 + r) / (6.d0 + r * (3.d0 + r))

    else if (limiter < 20) then
       ! Bruenn
       lambda = 1.d0 / (3.d0 + r)

    else if (limiter < 30) then
       ! Larsen's square root
       lambda = 1.d0 / sqrt(9.d0 + r**2)

    else if (limiter < 40) then 
       ! Minerbo
       if (r .lt. 1.5d0) then
          lambda = 2.d0/(3.d0 + sqrt(9.d0+12.d0*r**2))
       else 
          lambda = 1.d0/(1.d0+r+sqrt(1.d0+2.d0*r))
       end if

    else
       call bl_error("Unknown limiter ", limiter)
    endif
  end function FLDlambda

end module rad_util_module
