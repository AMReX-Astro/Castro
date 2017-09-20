module riemann_module

  use bl_types
  use bl_constants_module

  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QPRES, QREINT, &
                                 QFS, QFX, &
                                 NGDNV, GDU, GDPRES, QGAMC, QC, QCSML, &
#ifdef RADIATION
                                 GDERADS, GDLAMS, QGAMCG, QLAMS, &
#endif
                                 small_temp, small_dens, small_pres, &
                                 npassive, upass_map, qpass_map, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 fix_mass_flux, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(lo, hi, domlo, domhi, &
                    qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, qg_lo, qg_hi, &
#ifdef RADIATION
                    rflx, rflx_lo, rflx_hi, &
#endif
                    qaux, qa_lo, qa_hi, &
                    ilo, ihi)


    use eos_type_module, only: eos_input_re, eos_input_rt, eos_t
    use eos_module, only: eos
    use network, only: nspec, naux
#ifdef RADIATION
    use rad_params_module, only : ngroups
#endif
    use actual_riemann_module, only : riemanncg, riemannus

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer qpd_lo(3), qpd_hi(3)
    integer flx_lo(3), flx_hi(3)
    integer  qg_lo(3), qg_hi(3)
    integer  qa_lo(3), qa_hi(3)

    real(rt)            qm(qpd_lo(1):qpd_hi(1), NQ)
    real(rt)            qp(qpd_lo(1):qpd_hi(1), NQ)

    real(rt)           flx(flx_lo(1):flx_hi(1), NVAR)
    real(rt)          qint( qg_lo(1): qg_hi(1), NGDNV)
    real(rt)          qaux( qa_lo(1): qa_hi(1), NQAUX)

#ifdef RADIATION
    integer rflx_lo(3), rflx_hi(3)
    real(rt)         rflx(rflx_lo(1):rflx_hi(1), 0:ngroups-1)
#endif

    ! Local variables
    integer i

    type (eos_t) :: eos_state


    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent -- THIS PROBABLY DOESN"T WORK WITH RADIATION

       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are
       ! thermodynamically consistent.

       do i = ilo, ihi+1

          ! this is an initial guess for iterations, since we
          ! can't be certain that temp is on interfaces
          eos_state%T = 10000.0e0_rt

          ! minus state
          eos_state % rho = qm(i,QRHO)
          eos_state % p   = qm(i,QPRES)
          eos_state % e   = qm(i,QREINT)/qm(i,QRHO)
          eos_state % xn  = qm(i,QFS:QFS-1+nspec)
          eos_state % aux = qm(i,QFX:QFX-1+naux)

          ! Protect against negative energies

          if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
             eos_state % T = small_temp
             call eos(eos_input_rt, eos_state)
          else
             call eos(eos_input_re, eos_state)
          endif

          qm(i,QREINT) = qm(i,QRHO)*eos_state%e
          qm(i,QPRES) = eos_state%p
          !gamcm(i) = eos_state%gam1


          ! plus state
          eos_state % rho = qp(i,QRHO)
          eos_state % p   = qp(i,QPRES)
          eos_state % e   = qp(i,QREINT)/qp(i,QRHO)
          eos_state % xn  = qp(i,QFS:QFS-1+nspec)
          eos_state % aux = qp(i,QFX:QFX-1+naux)
          
          ! Protect against negative energies

          if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
             eos_state % T = small_temp
             call eos(eos_input_rt, eos_state)
          else
             call eos(eos_input_re, eos_state)
          endif
          
          qp(i,QREINT) = qp(i,QRHO)*eos_state%e
          qp(i,QPRES) = eos_state%p
          !gamcp(i) = eos_state%gam1
          
       enddo

    endif


    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson
       call riemannus(qm, qp, qpd_lo, qpd_hi, &
                      qaux, qa_lo, qa_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
#ifdef RADIATION
                      rflx, rflx_lo, rflx_hi, &
#endif
                      1, ilo, ihi+1, 0, 0, 0, 0, 0, &
                      [domlo(1), 0, 0], [domhi(1), 0, 0])

    elseif (riemann_solver == 1) then
       ! Colella & Glaz
       call riemanncg(qm, qp, qpd_lo, qpd_hi, &
                      qaux, qa_lo, qa_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, qg_lo, qg_hi, &
                      1, ilo, ihi+1, 0, 0, 0, 0, 0, &
                      [domlo(1), 0, 0], [domhi(1), 0, 0])
    else
       call bl_error("ERROR: HLLC not support in 1-d yet")
    endif

  end subroutine cmpflx


end module riemann_module
