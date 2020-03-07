  subroutine compute_ptot_ctot(lam, q, cg, ptot, ctot, gamc_tot)

    use meth_params_module, only : QPRES, QRHO, comoving, QRAD, QPTOT, NQ
    use fluxlimiter_module, only : Edd_factor ! function
    use amrex_fort_module, only : rt => amrex_real

    real(rt), intent(in) :: lam(0:ngroups-1)
    real(rt), intent(in) :: q(NQ)
    real(rt), intent(in) :: cg
    real(rt), intent(out) :: ptot
    real(rt), intent(out) :: ctot
    real(rt), intent(out) :: gamc_tot

    integer :: g

    real(rt) :: csrad2, Eddf, gamr, prad

    !$gpu

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

