module rad_util_module

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

end module rad_util_module
