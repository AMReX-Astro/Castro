module ctu_module
  !
  ! advection routines in support of the CTU unsplit advection scheme

  use amrex_constants_module
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine ctu_ppm_states(lo, hi, &
                            vlo, vhi, &
                            q, qd_lo, qd_hi, &
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            srcQ, src_lo, src_hi, &
                            qxm, qxm_lo, qxm_hi, &
                            qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM >= 2
                            qym, qym_lo, qym_hi, &
                            qyp, qyp_lo, qyp_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            qzm, qzm_lo, qzm_hi, &
                            qzp, qzp_lo, qzp_hi, &
#endif
                            dx, dt, &
#if AMREX_SPACEDIM < 3
                            dloga, dloga_lo, dloga_hi, &
#endif
                            domlo, domhi) bind(C, name="ctu_ppm_states")
    ! Compute the normal interface states by reconstructing
    ! the primitive variables using the piecewise parabolic method
    ! and doing characteristic tracing.  We do not apply the
    ! transverse terms here.

    use meth_params_module, only : NQSRC, NQ, NVAR, NQAUX


#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use trace_ppm_module, only : trace_ppm
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
#endif
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3)   ! grid spacing in X, Y, Z direction
    real(rt), intent(in), value :: dt    ! time stepsize
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)   ! input state, primitives
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)   ! auxiliary hydro data
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))   ! flattening parameter
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NQSRC)   ! primitive variable source

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), NQ)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), NQ)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), NQ)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n, idir

    !$gpu

    hdt = HALF*dt

    do idir = 1, AMREX_SPACEDIM


       ! compute the interface states

#ifdef RADIATION
       if (idir == 1) then
          call trace_ppm_rad(lo, hi, &
                             1, &
                             q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             srcQ, src_lo, src_hi, &
                             flatn, f_lo, f_hi, &
                             qxm, qxm_lo, qxm_hi, &
                             qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM <= 2
                             dloga, dloga_lo, dloga_hi, &
#endif
                             vlo, vhi, domlo, domhi, &
                             dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_ppm_rad(lo, hi, &
                             2, &
                             q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             srcQ, src_lo, src_hi, &
                             flatn, f_lo, f_hi, &
                             qym, qym_lo, qym_hi, &
                             qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM == 2
                             dloga, dloga_lo, dloga_hi, &
#endif
                             vlo, vhi, domlo, domhi, &
                             dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_ppm_rad(lo, hi, &
                             3, &
                             q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             srcQ, src_lo, src_hi, &
                             flatn, f_lo, f_hi, &
                             qzm, qzm_lo, qzm_hi, &
                             qzp, qzp_lo, qzp_hi, &
                             vlo, vhi, domlo, domhi, &
                             dx, dt)
#endif
       endif
#else
       ! hydro (no radiation)
       if (idir == 1) then
          call trace_ppm(lo, hi, &
                         1, &
                         q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         srcQ, src_lo, src_hi, &
                         flatn, f_lo, f_hi, &
                         qxm, qxm_lo, qxm_hi, &
                         qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM <= 2
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_ppm(lo, hi, &
                         2, &
                         q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         srcQ, src_lo, src_hi, &
                         flatn, f_lo, f_hi, &
                         qym, qym_lo, qym_hi, &
                         qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM == 2
                         dloga, dloga_lo, dloga_hi, &
#endif
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_ppm(lo, hi, &
                         3, &
                         q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         srcQ, src_lo, src_hi, &
                         flatn, f_lo, f_hi, &
                         qzm, qzm_lo, qzm_hi, &
                         qzp, qzp_lo, qzp_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif
       end if
#endif

    end do

  end subroutine ctu_ppm_states


  subroutine ctu_plm_states(lo, hi, &
                            vlo, vhi, &
                            q, qd_lo, qd_hi, &
                            flatn, f_lo, f_hi, &
                            qaux, qa_lo, qa_hi, &
                            srcQ, src_lo, src_hi, &
                            dq, dq_lo, dq_hi, &
                            qxm, qxm_lo, qxm_hi, &
                            qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM >= 2
                            qym, qym_lo, qym_hi, &
                            qyp, qyp_lo, qyp_hi, &
#endif
#if AMREX_SPACEDIM == 3
                            qzm, qzm_lo, qzm_hi, &
                            qzp, qzp_lo, qzp_hi, &
#endif
                            dx, dt, &
#if AMREX_SPACEDIM < 3
                            dloga, dloga_lo, dloga_hi, &
#endif
                            domlo, domhi) bind(C, name="ctu_plm_states")
    ! Compute the normal interface states by reconstructing
    ! the primitive variables using piecewise linear slopes and doing
    ! characteristic tracing.  We do not apply the transverse terms here.
    !
    ! .. todo::
    !    we can get rid of the the different temporary q Godunov
    !    state arrays
    !

    use meth_params_module, only : NQSRC, NQ, NVAR, &
                                   QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, QGAME, QREINT, &
                                   use_pslope
    use trace_plm_module, only : trace_plm
    use slope_module, only : uslope, pslope

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
#endif
#if AMREX_SPACEDIM < 3
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
    real(rt), intent(in) :: dx(3)   ! grid spacing in X, Y, Z direction
    real(rt), intent(in), value :: dt   ! time stepsize
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)   ! input state, primitives
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)   ! auxiliary hydro data
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))   ! flattening parameter
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NQSRC)   ! primitive variable source

    real(rt), intent(inout) :: dq(dq_lo(1):dq_hi(1), dq_lo(2):dq_hi(2), dq_lo(3):dq_hi(3), NQ)

    real(rt), intent(inout) :: qxm(qxm_lo(1):qxm_hi(1), qxm_lo(2):qxm_hi(2), qxm_lo(3):qxm_hi(3), NQ)
    real(rt), intent(inout) :: qxp(qxp_lo(1):qxp_hi(1), qxp_lo(2):qxp_hi(2), qxp_lo(3):qxp_hi(3), NQ)
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: qym(qym_lo(1):qym_hi(1), qym_lo(2):qym_hi(2), qym_lo(3):qym_hi(3), NQ)
    real(rt), intent(inout) :: qyp(qyp_lo(1):qyp_hi(1), qyp_lo(2):qyp_hi(2), qyp_lo(3):qyp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: qzm(qzm_lo(1):qzm_hi(1), qzm_lo(2):qzm_hi(2), qzm_lo(3):qzm_hi(3), NQ)
    real(rt), intent(inout) :: qzp(qzp_lo(1):qzp_hi(1), qzp_lo(2):qzp_hi(2), qzp_lo(3):qzp_hi(3), NQ)
#endif
#if AMREX_SPACEDIM < 3
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
#endif
    real(rt) :: hdt
    integer :: i, j, k, n, idir

    logical :: reconstruct_state(NQ)

    !$gpu

    hdt = HALF*dt

    ! we don't need to reconstruct all of the NQ state variables,
    ! depending on how we are tracing
    reconstruct_state(:) = .true.
    reconstruct_state(QGAME) = .false.
    reconstruct_state(QTEMP) = .false.


#ifdef RADIATION
#ifndef AMREX_USE_CUDA
    call castro_error("ppm_type <=0 is not supported in with radiation")
#endif
#endif

    ! Compute all slopes
    do idir = 1, AMREX_SPACEDIM

       do n = 1, NQ
          if (.not. reconstruct_state(n)) cycle
          call uslope(lo, hi, idir, &
                      q, qd_lo, qd_hi, n, &
                      flatn, f_lo, f_hi, &
                      dq, dq_lo, dq_hi, &
                      dx, domlo, domhi)
       end do

       if (use_pslope == 1) then
          call pslope(lo, hi, idir, &
                      q, qd_lo, qd_hi, &
                      flatn, f_lo, f_hi, &
                      dq, dq_lo, dq_hi, &
                      srcQ, src_lo, src_hi, &
                      dx)
       endif


       ! compute the interface states

       if (idir == 1) then
          call trace_plm(lo, hi, &
                         1, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq, dq_lo, dq_hi, &
                         qxm, qxm_lo, qxm_hi, &
                         qxp, qxp_lo, qxp_hi, &
#if AMREX_SPACEDIM < 3
                         dloga, dloga_lo, dloga_hi, &
#endif
                         SrcQ, src_lo, src_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then
          call trace_plm(lo, hi, &
                         2, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq, dq_lo, dq_hi, &
                         qym, qym_lo, qym_hi, &
                         qyp, qyp_lo, qyp_hi, &
#if AMREX_SPACEDIM < 3
                         dloga, dloga_lo, dloga_hi, &
#endif
                         SrcQ, src_lo, src_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif

#if AMREX_SPACEDIM == 3
       else
          call trace_plm(lo, hi, &
                         3, q, qd_lo, qd_hi, &
                         qaux, qa_lo, qa_hi, &
                         dq, dq_lo, dq_hi, &
                         qzm, qzm_lo, qzm_hi, &
                         qzp, qzp_lo, qzp_hi, &
                         SrcQ, src_lo, src_hi, &
                         vlo, vhi, domlo, domhi, &
                         dx, dt)
#endif
       end if

    end do

  end subroutine ctu_plm_states


  subroutine ca_track_grid_losses(lo, hi, &
       flux1, flux1_lo, flux1_hi, &
#if AMREX_SPACEDIM >= 2
       flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
       flux3, flux3_lo, flux3_hi, &
#endif
       mass_lost, xmom_lost, ymom_lost, zmom_lost, &
       eden_lost, xang_lost, yang_lost, zang_lost) &
       bind(C, name="ca_track_grid_losses")


    use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN, NVAR
    use amrinfo_module, only : amr_level
    use prob_params_module, only : domlo_level, domhi_level, center
    use castro_util_module, only: position ! function
    use castro_util_module, only: linear_to_angular_momentum ! function
    use reduction_module, only: reduce_add

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
#if AMREX_SPACEDIM >= 2
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
#endif
#if AMREX_SPACEDIM == 3
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
#endif

    real(rt), intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt), intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    real(rt) :: loc(3), ang_mom(3), flux(3)
    integer :: domlo(3), domhi(3)
    integer :: i, j, k

    !$gpu

    domlo = domlo_level(:,amr_level)
    domhi = domhi_level(:,amr_level)

#if AMREX_SPACEDIM == 3
    if (lo(3) .le. domlo(3) .and. hi(3) .ge. domlo(3)) then

       k = domlo(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccz=.false.) - center

             call reduce_add(mass_lost, -flux3(i,j,k,URHO))
             call reduce_add(xmom_lost, -flux3(i,j,k,UMX))
             call reduce_add(ymom_lost, -flux3(i,j,k,UMY))
             call reduce_add(zmom_lost, -flux3(i,j,k,UMZ))
             call reduce_add(eden_lost, -flux3(i,j,k,UEDEN))

             flux(:) = flux3(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call reduce_add(xang_lost, -ang_mom(1))
             call reduce_add(yang_lost, -ang_mom(2))
             call reduce_add(zang_lost, -ang_mom(3))

          enddo
       enddo

    endif

    if (lo(3) .le. domhi(3) .and. hi(3) .ge. domhi(3)) then

       k = domhi(3) + 1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccz=.false.) - center

             call reduce_add(mass_lost, flux3(i,j,k,URHO))
             call reduce_add(xmom_lost, flux3(i,j,k,UMX))
             call reduce_add(ymom_lost, flux3(i,j,k,UMY))
             call reduce_add(zmom_lost, flux3(i,j,k,UMZ))
             call reduce_add(eden_lost, flux3(i,j,k,UEDEN))

             flux(:) = flux3(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call reduce_add(xang_lost, ang_mom(1))
             call reduce_add(yang_lost, ang_mom(2))
             call reduce_add(zang_lost, ang_mom(3))

          enddo
       enddo

    endif
#endif

#if AMREX_SPACEDIM >= 2
    if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

       j = domlo(2)
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.) - center

             call reduce_add(mass_lost, -flux2(i,j,k,URHO))
             call reduce_add(xmom_lost, -flux2(i,j,k,UMX))
             call reduce_add(ymom_lost, -flux2(i,j,k,UMY))
             call reduce_add(zmom_lost, -flux2(i,j,k,UMZ))
             call reduce_add(eden_lost, -flux2(i,j,k,UEDEN))

             flux(:) = flux2(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call reduce_add(xang_lost, -ang_mom(1))
             call reduce_add(yang_lost, -ang_mom(2))
             call reduce_add(zang_lost, -ang_mom(3))

          enddo
       enddo

    endif

    if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

       j = domhi(2) + 1
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.) - center

             call reduce_add(mass_lost, flux2(i,j,k,URHO))
             call reduce_add(xmom_lost, flux2(i,j,k,UMX))
             call reduce_add(ymom_lost, flux2(i,j,k,UMY))
             call reduce_add(zmom_lost, flux2(i,j,k,UMZ))
             call reduce_add(eden_lost, flux2(i,j,k,UEDEN))

             flux(:) = flux2(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call reduce_add(xang_lost, ang_mom(1))
             call reduce_add(yang_lost, ang_mom(2))
             call reduce_add(zang_lost, ang_mom(3))

          enddo
       enddo

    endif
#endif

    if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

       i = domlo(1)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             loc = position(i,j,k,ccx=.false.) - center

             call reduce_add(mass_lost, -flux1(i,j,k,URHO))
             call reduce_add(xmom_lost, -flux1(i,j,k,UMX))
             call reduce_add(ymom_lost, -flux1(i,j,k,UMY))
             call reduce_add(zmom_lost, -flux1(i,j,k,UMZ))
             call reduce_add(eden_lost, -flux1(i,j,k,UEDEN))

             flux(:) = flux1(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call reduce_add(xang_lost, -ang_mom(1))
             call reduce_add(yang_lost, -ang_mom(2))
             call reduce_add(zang_lost, -ang_mom(3))

          enddo
       enddo

    endif

    if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

       i = domhi(1) + 1
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             loc = position(i,j,k,ccx=.false.) - center

             call reduce_add(mass_lost, flux1(i,j,k,URHO))
             call reduce_add(xmom_lost, flux1(i,j,k,UMX))
             call reduce_add(ymom_lost, flux1(i,j,k,UMY))
             call reduce_add(zmom_lost, flux1(i,j,k,UMZ))
             call reduce_add(eden_lost, flux1(i,j,k,UEDEN))

             flux(:) = flux1(i,j,k,UMX:UMZ)
             ang_mom = linear_to_angular_momentum(loc, flux)
             call reduce_add(xang_lost, ang_mom(1))
             call reduce_add(yang_lost, ang_mom(2))
             call reduce_add(zang_lost, ang_mom(3))

          enddo
       enddo

    endif

  end subroutine ca_track_grid_losses

end module ctu_module
