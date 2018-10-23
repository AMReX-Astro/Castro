module ctu_advection_module

  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH

  use amrex_error_module, only : amrex_error
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public umeth

contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH     Compute hyperbolic fluxes using unsplit second
! ::: ::            order Godunov integrator.
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: qaux        => (const)  auxiliary hydro info
! ::: :: flatn       => (const)  flattening parameter
! ::: :: srcQ        => (const)  source for primitive variables
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux      <=  (modify) flux in X direction on X edges
! ::: ----------------------------------------------------------------

  subroutine umeth(q, q_lo, q_hi, &
                   flatn, &
                   qaux, qa_lo, qa_hi, &
                   srcQ, src_lo, src_hi, &
                   lo, hi, dx, dt, &
                   uout, uout_lo, uout_hi, &
                   flux, fd_lo, fd_hi, &
#ifdef RADIATION
                   rflux, rfd_lo, rfd_hi, &
#endif
                   q1, q1_lo, q1_hi, &
                   area1, area1_lo, area1_hi, &
                   vol, vol_lo, vol_hi, &
                   dloga, dloga_lo, dloga_hi, &
                   domlo, domhi)

    use meth_params_module, only : QVAR, NQ, NVAR, &
                                   QC, QFS, QFX, QGAMC, QU, QRHO, QTEMP, QPRES, QREINT, QGAME, &
                                   NQAUX, NGDNV, &
                                   ppm_type, hybrid_riemann, ppm_predict_gammae, &
                                   use_pslope, plm_iorder, ppm_temp_fix
    use riemann_module, only : cmpflx
    use trace_plm_module, only : trace_plm
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use trace_ppm_module, only : trace_ppm
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif

    use amrex_fort_module, only : rt => amrex_real
    use advection_util_module, only : shock
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux
    use ppm_module, only : ppm_reconstruct, ppm_int_profile
    use slope_module, only : uslope, pslope
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: fd_lo(3), fd_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rfd_lo(3), rfd_hi(3)
#endif
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)

    real(rt)        , intent(in) :: dx(3), dt
    real(rt)        , intent(in) :: q(   q_lo(1):q_hi(1),NQ)
    real(rt)        , intent(in) :: qaux(   qa_lo(1):qa_hi(1),NQAUX)
    real(rt)        , intent(in) :: flatn(   q_lo(1):q_hi(1))
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),NVAR)
    real(rt)        , intent(inout) :: flux(fd_lo(1):fd_hi(1),NVAR)
#ifdef RADIATION
    real(rt)        , intent(inout) :: rflux(rfd_lo(1):rfd_hi(1),0:ngroups-1)
#endif
    real(rt)        , intent(in) :: srcQ(src_lo(1)  :src_hi(1),QVAR)
    real(rt)        , intent(inout) :: q1(q1_lo(1):q1_hi(1),NGDNV)
    real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1))
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1))
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1))

    real(rt)        , allocatable :: shk(:)

    integer :: i, n, iwave

    ! Left and right state arrays (edge centered, cell centered)
    real(rt)        , allocatable:: dq(:,:), qm(:,:), qp(:,:)

    integer :: qp_lo(3), qp_hi(3)
    integer :: shk_lo(3), shk_hi(3)
    integer :: I_lo(3), I_hi(3)

    real(rt), allocatable :: Ip(:,:,:)
    real(rt), allocatable :: Im(:,:,:)

    real(rt), allocatable :: Ip_src(:,:,:)
    real(rt), allocatable :: Im_src(:,:,:)

    ! gamma_c/1 on the interfaces
    real(rt), allocatable :: Ip_gc(:,:,:)
    real(rt), allocatable :: Im_gc(:,:,:)

    ! temporary interface values of the parabola
    real(rt), allocatable :: sxm(:), sxp(:)

    type(eos_t) :: eos_state

    logical :: source_nonzero(QVAR)
    logical :: reconstruct_state(NQ)

    qp_lo = lo - dg
    qp_hi = hi + 2*dg

    shk_lo = lo - dg
    shk_hi = hi + dg

    I_lo = lo - dg
    I_hi = hi + dg


    if (ppm_type > 0) then
       ! indices: (x, y, dimension, wave, variable)
       allocate(Ip(I_lo(1):I_hi(1), 3, NQ))
       allocate(Im(I_lo(1):I_hi(1), 3, NQ))

       allocate(Ip_src(I_lo(1):I_hi(1), 3, QVAR))
       allocate(Im_src(I_lo(1):I_hi(1), 3, QVAR))

       allocate(Ip_gc(I_lo(1):I_hi(1), 3, 1))
       allocate(Im_gc(I_lo(1):I_hi(1), 3, 1))

       allocate(sxm(q_lo(1):q_hi(1)))
       allocate(sxp(q_lo(1):q_hi(1)))
    endif

    allocate (shk(shk_lo(1):shk_hi(1)))


#ifdef SHOCK_VAR
    uout(lo(1):hi(1),USHK) = ZERO

    call shock(q, q_lo, q_hi, &
               shk, shk_lo, shk_hi, &
               lo, hi, dx)

    ! Store the shock data for future use in the burning step.
    do i = lo(1), hi(1)
       uout(i,USHK) = shk(i)
    enddo

    ! Discard it locally if we don't need it in the hydro update.
    if (hybrid_riemann /= 1) then
       shk(:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, q_lo, q_hi, &
                  shk, shk_lo, shk_hi, &
                  lo, hi, dx)
    else
       shk(:) = ZERO
    endif
#endif

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    allocate ( qm(qp_lo(1):qp_hi(1),NQ))
    allocate ( qp(qp_lo(1):qp_hi(1),NQ))

    ! we don't need to reconstruct all of the NQ state variables,
    ! depending on how we are tracing
    reconstruct_state(:) = .true.
    if (ppm_predict_gammae /= 1) then
       reconstruct_state(QGAME) = .false.
    else
       reconstruct_state(QREINT) = .false.
    endif
    if (ppm_temp_fix < 3) then
       reconstruct_state(QTEMP) = .false.
    endif

    if (ppm_type == 0) then
       allocate(dq(q_lo(1):q_hi(1),NQ))

       if (plm_iorder == 1) then
          dq(lo(1)-1:hi(1)+1,1:NQ) = ZERO

       else
          ! piecewise linear slope
          do n = 1, NQ
             if (.not. reconstruct_state(n)) cycle
             call uslope(q, q_lo, q_hi, n, &
                         flatn, q_lo, q_hi, &
                         dq, q_lo, q_hi, &
                         lo, hi)
          end do

          if (use_pslope == 1) &
               call pslope(q, q_lo, q_hi, &
                           flatn, q_lo, q_hi, &
                           dq, q_lo, q_hi, &
                           srcQ, src_lo, src_hi, &
                           lo, hi, dx)

       endif

    else

       do n = 1, QVAR
          if (minval(srcQ(lo(1)-2:hi(1)+2,n)) == ZERO .and. &
              maxval(srcQ(lo(1)-2:hi(1)+2,n)) == ZERO) then
             source_nonzero(n) = .false.
          else
             source_nonzero(n) = .true.
          endif
       enddo

       ! Compute Ip and Im -- this does the parabolic reconstruction,
       ! limiting, and returns the integral of each profile under each
       ! wave to each interface
       do n = 1, NQ
          if (.not. reconstruct_state(n)) cycle
          call ppm_reconstruct(q, q_lo, q_hi, NQ, n, &
                               flatn, q_lo, q_hi, &
                               sxm, sxp, q_lo, q_hi, &
                               lo, hi, dx)

          call ppm_int_profile(q, q_lo, q_hi, NQ, n, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, q_lo, q_hi, &
                               Ip, Im, I_lo, I_hi, NQ, n, &
                               lo, hi, dx, dt)
       enddo

       ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
       ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
       if (ppm_temp_fix == 1) then
          do i = lo(1)-1, hi(1)+1
             do iwave = 1, 3
                eos_state%rho   = Ip(i,iwave,QRHO)
                eos_state%T     = Ip(i,iwave,QTEMP)
                eos_state%xn(:) = Ip(i,iwave,QFS:QFS-1+nspec)
                eos_state%aux   = Ip(i,iwave,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                Ip(i,iwave,QPRES) = eos_state%p
                Ip(i,iwave,QREINT) = Ip(i,iwave,QRHO)*eos_state%e
                Ip_gc(i,iwave,1) = eos_state%gam1

                eos_state%rho   = Im(i,iwave,QRHO)
                eos_state%T     = Im(i,iwave,QTEMP)
                eos_state%xn(:) = Im(i,iwave,QFS:QFS-1+nspec)
                eos_state%aux   = Im(i,iwave,QFX:QFX-1+naux)

                call eos(eos_input_rt, eos_state)

                Im(i,iwave,QPRES) = eos_state%p
                Im(i,iwave,QREINT) = Im(i,iwave,QRHO)*eos_state%e
                Im_gc(i,iwave,1) = eos_state%gam1
             enddo
          enddo

       endif


       ! get an edge-based gam1 here if we didn't get it from the EOS
       ! call above (for ppm_temp_fix = 1)
       if (ppm_temp_fix /= 1) then
          call ppm_reconstruct(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               flatn, q_lo, q_hi, &
                               sxm, sxp, q_lo, q_hi, &
                               lo, hi, dx)

          call ppm_int_profile(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, q_lo, q_hi, &
                               Ip_gc, Im_gc, I_lo, I_hi, 1, 1, &
                               lo, hi, dx, dt)
       endif

       do n = 1, QVAR
          if (source_nonzero(n)) then
             call ppm_reconstruct(srcQ, src_lo, src_hi, QVAR, n, &
                                  flatn, q_lo, q_hi, &
                                  sxm, sxp, q_lo, q_hi, &
                                  lo, hi, dx)

             call ppm_int_profile(srcQ, src_lo, src_hi, QVAR, n, &
                                  q, q_lo, q_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, q_lo, q_hi, &
                                  Ip_src, Im_src, I_lo, I_hi, QVAR, n, &
                                  lo, hi, dx, dt)
          else
             Ip_src(I_lo(1):I_hi(1),:,n) = ZERO
             Im_src(I_lo(1):I_hi(1),:,n) = ZERO
          endif
       enddo

       deallocate(sxm, sxp)
    endif

    ! Trace to edges
    if (ppm_type .gt. 0) then
#ifdef RADIATION
       call trace_ppm_rad(1, q, q_lo, q_hi, &
                          qaux, qa_lo, qa_hi, &
                          Ip, Im, Ip_src, Im_src, I_lo, I_hi, &
                          qm, qp, qp_lo, qp_hi, &
                          dloga, dloga_lo, dloga_hi, &
                          lo, hi, domlo, domhi, &
                          dx, dt)
#else
       call trace_ppm(1, q, q_lo, q_hi, &
                      qaux, qa_lo, qa_hi, &
                      Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                      qm, qp, qp_lo, qp_hi, &
                      dloga, dloga_lo, dloga_hi, &
                      lo, hi, domlo, domhi, &
                      dx, dt)
#endif
    else

       ! PLM

#ifdef RADIATION
#ifndef AMREX_USE_CUDA
       call amrex_error("ppm_type <=0 is not supported in umeth with radiation")
#endif
#else

       call trace_plm(1, q, q_lo, q_hi, &
                      qaux, qa_lo, qa_hi, &
                      dq, q_lo, q_hi, &
                      qm, qp, qp_lo, qp_hi, &
                      dloga, dloga_lo, dloga_hi, &
                      srcQ, src_lo, src_hi, &
                      lo, hi, domlo, domhi, &
                      dx, dt)

       deallocate(dq)
#endif
    end if

    ! Solve Riemann problem, compute xflux from improved predicted states
    call cmpflx(qm, qp, qp_lo, qp_hi, &
                flux, fd_lo, fd_hi, &
                q1, q1_lo, q1_hi, &
#ifdef RADIATION
                rflux, rfd_lo,rfd_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, lo, hi+dg(:), &
                domlo, domhi)

    deallocate (qm,qp)

  end subroutine umeth

end module ctu_advection_module
