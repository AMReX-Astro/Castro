module ctu_advection_module

  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH, TWO
  use amrex_error_module
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
! ::: :: qaux        => (const)  auxiliary hydro data
! ::: :: flatn       => (const)  flattening parameter
! ::: :: srcQ        => (const)  primitive variable source
! ::: :: dx          => (const)  grid spacing in X, Y, Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: q1         <=  (modify) Godunov interface state in X
! ::: :: q2         <=  (modify) Godunov interface state in Y
! ::: ----------------------------------------------------------------

  subroutine umeth(q, q_lo, q_hi, &
                   flatn, &
                   qaux, qa_lo, qa_hi, &
                   srcQ, src_lo, src_hi, &
                   lo, hi, dx, dt, &
                   uout, uout_lo, uout_hi, &
                   flux1, f1_lo, f1_hi, &
                   flux2, f2_lo, f2_hi, &
#ifdef RADIATION
                   rflux1, rf1_lo, rf1_hi, &
                   rflux2, rf2_lo, rf2_hi, &
#endif
                   q1, q1_lo, q1_hi, &
                   q2, q2_lo, q2_hi, &
                   area1, area1_lo, area1_hi, &
                   area2, area2_lo, area2_hi, &
                   vol, vol_lo, vol_hi, &
                   dloga, dloga_lo, dloga_hi, &
                   domlo, domhi)

    use meth_params_module, only : NQ, QVAR, NVAR, ppm_type, hybrid_riemann, &
                                   QC, QFS, QFX, QGAMC, QU, QV, QRHO, QTEMP, QPRES, QREINT, &
                                   GDU, GDV, GDPRES, NGDNV, NQ, QGAME, &
                                   NQAUX, &
                                   ppm_type, ppm_predict_gammae, &
                                   use_pslope, plm_iorder, ppm_temp_fix
    use network, only : nspec, naux
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use trace_plm_module, only : trace_plm
    use ppm_module, only : ppm_reconstruct, ppm_int_profile
    use slope_module, only : uslope, pslope
    use multid_slope_module, only : multid_slope
    use riemann_module, only: cmpflx
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : trace_ppm_rad
#else
    use trace_ppm_module, only : trace_ppm, trace_ppm_gammae, trace_ppm_temp
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use advection_util_module, only : shock

    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: f1_lo(3), f1_hi(3)
    integer, intent(in) :: f2_lo(3), f2_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rf1_lo(3), rf1_hi(3)
    integer, intent(in) :: rf2_lo(3), rf2_hi(3)
#endif
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: dx(3), dt
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
    real(rt), intent(in) :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2))
    real(rt), intent(in) :: srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),QVAR)
    real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2))
    real(rt), intent(inout) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),NGDNV)
    real(rt), intent(inout) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),NGDNV)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),NVAR)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),NVAR)
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),NVAR)
#ifdef RADIATION
    real(rt), intent(inout) :: rflux1(rf1_lo(1):rf1_hi(1),rf1_lo(2):rf1_hi(2),0:ngroups-1)
    real(rt), intent(inout) :: rflux2(rf2_lo(1):rf2_hi(1),rf2_lo(2):rf2_hi(2),0:ngroups-1)
#endif
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))

    ! Left and right state arrays (edge centered, cell centered)
    real(rt), allocatable::  qm(:,:,:),  qp(:,:,:)
    real(rt), allocatable:: qxm(:,:,:), qym(:,:,:)
    real(rt), allocatable:: qxp(:,:,:), qyp(:,:,:)

    ! Work arrays to hold riemann state and conservative fluxes

    real(rt), allocatable ::  fx(:,:,:),  fy(:,:,:)
#ifdef RADIATION
    real(rt), allocatable ::  rfx(:,:,:),  rfy(:,:,:)
#endif
    real(rt), allocatable ::  qgdxtmp(:,:,:)
    real(rt), allocatable :: shk(:,:)

    ! Local scalar variables
    real(rt) :: dtdx
    real(rt) :: hdtdx, hdt, hdtdy

    integer :: i, j, iwave, idim, n

    integer :: tflx_lo(3), tflx_hi(3)
    integer :: tfly_lo(3), tfly_hi(3)
    integer :: shk_lo(3), shk_hi(3)
    integer :: qp_lo(3), qp_hi(3)

    real(rt), allocatable :: Ip(:,:,:,:,:)
    real(rt), allocatable :: Im(:,:,:,:,:)

    real(rt), allocatable :: Ip_src(:,:,:,:,:)
    real(rt), allocatable :: Im_src(:,:,:,:,:)

    ! gamma_c/1 on the interfaces
    real(rt), allocatable :: Ip_gc(:,:,:,:,:)
    real(rt), allocatable :: Im_gc(:,:,:,:,:)

    ! temporary interface values of the parabola
    real(rt), allocatable :: sxm(:,:), sxp(:,:), sym(:,:), syp(:,:)

    real(rt), allocatable :: dqx(:,:,:), dqy(:,:,:)

    integer :: I_lo(3), I_hi(3)

    type(eos_t) :: eos_state

    logical :: source_nonzero(QVAR)
    logical :: reconstruct_state(NQ)

    tflx_lo = [lo(1), lo(2)-1, 0]
    tflx_hi = [hi(1)+1, hi(2)+1, 0]

    tfly_lo = [lo(1)-1, lo(2), 0]
    tfly_hi = [hi(1)+1, hi(2)+1, 0]

    shk_lo = [lo(1)-1, lo(2)-1, 0]
    shk_hi = [hi(1)+1, hi(2)+1, 0]

    qp_lo = [lo(1)-1, lo(2)-1, 0]
    qp_hi = [hi(1)+2, hi(2)+2, 0]


    I_lo = [lo(1)-1, lo(2)-1, 0]
    I_hi = [hi(1)+1, hi(2)+1, 0]

    if (ppm_type > 0) then
       ! indices: (x, y, dimension, wave, variable)
       allocate(Ip(I_lo(1):I_hi(1), I_lo(2):I_hi(2), 2, 3, NQ))
       allocate(Im(I_lo(1):I_hi(1), I_lo(2):I_hi(2), 2, 3, NQ))

       allocate(Ip_src(I_lo(1):I_hi(1), I_lo(2):I_hi(2), 2, 3, QVAR))
       allocate(Im_src(I_lo(1):I_hi(1), I_lo(2):I_hi(2), 2, 3, QVAR))

       allocate(Ip_gc(I_lo(1):I_hi(1), I_lo(2):I_hi(2), 2, 3, 1))
       allocate(Im_gc(I_lo(1):I_hi(1), I_lo(2):I_hi(2), 2, 3, 1))

       allocate(sxm(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
       allocate(sxp(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
       allocate(sym(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
       allocate(syp(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
    endif

    allocate ( qgdxtmp(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),NGDNV))

    allocate (  qm(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),NQ) )
    allocate (  qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),NQ) )
    allocate ( qxm(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),NQ) )
    allocate ( qxp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),NQ) )
    allocate ( qym(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),NQ) )
    allocate ( qyp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),NQ) )

    allocate (  fx(tflx_lo(1):tflx_hi(1),tflx_lo(2):tflx_hi(2),NVAR) )
    allocate (  fy(tfly_lo(1):tfly_hi(1),tfly_lo(2):tfly_hi(2),NVAR) )
#ifdef RADIATION
    allocate ( rfx(tflx_lo(1):tflx_hi(1),tflx_lo(2):tflx_hi(2),0:ngroups-1) )
    allocate ( rfy(tfly_lo(1):tfly_hi(1),tfly_lo(2):tfly_hi(2),0:ngroups-1) )
#endif

    allocate (shk(shk_lo(1):shk_hi(1),shk_lo(2):shk_hi(2)))


    ! Local constants
    dtdx = dt/dx(1)
    hdtdx = HALF*dtdx
    hdtdy = HALF*dt/dx(2)
    hdt = HALF*dt

    ! multidimensional shock detection

#ifdef SHOCK_VAR
    uout(lo(1):hi(1),lo(2):hi(2),USHK) = ZERO

    call shock(q, q_lo, q_hi, &
               shk, shk_lo, shk_hi, &
               lo, hi, dx)

    ! Store the shock data for future use in the burning step.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          uout(i,j,USHK) = shk(i,j)
       enddo
    enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, q_lo, q_hi, &
                  shk, shk_lo, shk_hi, &
                  lo, hi, dx)
    else
       shk(:,:) = ZERO
    endif
#endif

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

    ! preprocess the sources -- we don't want to trace under a source that is empty
    if (ppm_type > 0) then
       do n = 1, QVAR
          if (minval(srcQ(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,n)) == ZERO .and. &
               maxval(srcQ(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,n)) == ZERO) then
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
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               lo, hi, dx)

          call ppm_int_profile(q, q_lo, q_hi, NQ, n, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               Ip, Im, I_lo, I_hi, NQ, n, &
                               lo, hi, dx, dt)
       end do

       if (ppm_temp_fix /= 1) then
          call ppm_reconstruct(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               flatn, q_lo, q_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               lo, hi, dx)

          call ppm_int_profile(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               Ip_gc, Im_gc, I_lo, I_hi, 1, 1, &
                               lo, hi, dx, dt)
       else
          ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
          ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
          ! get an edge-based gam1 here if we didn't get it from the EOS
          ! call above (for ppm_temp_fix = 1)
          do iwave = 1, 3
             do idim = 1, 2

                do j = lo(2)-1, hi(2)+1
                   do i = lo(1)-1, hi(1)+1

                      eos_state%rho   = Ip(i,j,idim,iwave,QRHO)
                      eos_state%T     = Ip(i,j,idim,iwave,QTEMP)

                      eos_state%xn(:) = Ip(i,j,idim,iwave,QFS:QFS-1+nspec)
                      eos_state%aux   = Ip(i,j,idim,iwave,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      Ip(i,j,idim,iwave,QPRES) = eos_state % p
                      Ip(i,j,idim,iwave,QREINT) = Ip(i,j,idim,iwave,QRHO) * eos_state % e
                      Ip_gc(i,j,idim,iwave,1) = eos_state%gam1

                   end do
                end do

                do j = lo(2)-1, hi(2)+1
                   do i = lo(1)-1, hi(1)+1

                      eos_state%rho   = Im(i,j,idim,iwave,QRHO)
                      eos_state%T     = Im(i,j,idim,iwave,QTEMP)

                      eos_state%xn(:) = Im(i,j,idim,iwave,QFS:QFS-1+nspec)
                      eos_state%aux   = Im(i,j,idim,iwave,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      Im(i,j,idim,iwave,QPRES) = eos_state % p
                      Im(i,j,idim,iwave,QREINT) = Im(i,j,idim,iwave,QRHO) * eos_state % e
                      Im_gc(i,j,idim,iwave,1) = eos_state%gam1
                   enddo
                enddo
             enddo
          enddo

       endif

       ! source terms
       do n = 1, QVAR
          if (source_nonzero(n)) then
             call ppm_reconstruct(srcQ, src_lo, src_hi, QVAR, n, &
                                  flatn, q_lo, q_hi, &
                                  sxm, sxp, sym, syp, q_lo, q_hi, &
                                  lo, hi, dx)

             call ppm_int_profile(srcQ, src_lo, src_hi, QVAR, n, &
                                  q, q_lo, q_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, q_lo, q_hi, &
                                  Ip_src, Im_src, I_lo, I_hi, QVAR, n, &
                                  lo, hi, dx, dt)
          else
             Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),:,:,n) = ZERO
             Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),:,:,n) = ZERO
          endif
       enddo

       deallocate(sxm, sxp, sym, syp)

#ifdef RADIATION
       call trace_ppm_rad(1, q, q_lo, q_hi, &
                          qaux, qa_lo, qa_hi, &
                          Ip, Im, Ip_src, Im_src, I_lo, I_hi, &
                          qxm, qxp, qp_lo, qp_hi, &
                          dloga, dloga_lo, dloga_hi, &
                          lo, hi, domlo, domhi, &
                          dx, dt)

       call trace_ppm_rad(2, q, q_lo, q_hi, &
                          qaux, qa_lo, qa_hi, &
                          Ip, Im, Ip_src, Im_src, I_lo, I_hi, &
                          qym, qyp, qp_lo, qp_hi, &
                          dloga, dloga_lo, dloga_hi, &
                          lo, hi, domlo, domhi, &
                          dx, dt)

#else
       if (ppm_temp_fix < 3) then
          if (ppm_predict_gammae == 0) then
             call trace_ppm(1, q, q_lo, q_hi, &
                            qaux, qa_lo, qa_hi, &
                            Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                            qxm, qxp, qp_lo, qp_hi, &
                            dloga, dloga_lo, dloga_hi, &
                            lo, hi, domlo, domhi, &
                            dx, dt)

             call trace_ppm(2, q, q_lo, q_hi, &
                            qaux, qa_lo, qa_hi, &
                            Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                            qym, qyp, qp_lo, qp_hi, &
                            dloga, dloga_lo, dloga_hi, &
                            lo, hi, domlo, domhi, &
                            dx, dt)
          else
             call trace_ppm_gammae(1, q, q_lo, q_hi, &
                                   qaux, qa_lo, qa_hi, &
                                   Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                                   qxm, qxp, qp_lo, qp_hi, &
                                   dloga, dloga_lo, dloga_hi, &
                                   lo, hi, domlo, domhi, &
                                   dx, dt)

             call trace_ppm_gammae(2, q, q_lo, q_hi, &
                                   qaux, qa_lo, qa_hi, &
                                   Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                                   qym, qyp, qp_lo, qp_hi, &
                                   dloga, dloga_lo, dloga_hi, &
                                   lo, hi, domlo, domhi, &
                                   dx, dt)
          endif
       else
          call trace_ppm_temp(1, q, q_lo, q_hi, &
                              qaux, qa_lo, qa_hi, &
                              Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                              qxm, qxp, qp_lo, qp_hi, &
                              dloga, dloga_lo, dloga_hi, &
                              lo, hi, domlo, domhi, &
                              dx, dt)

          call trace_ppm_temp(2, q, q_lo, q_hi, &
                              qaux, qa_lo, qa_hi, &
                              Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                              qym, qyp, qp_lo, qp_hi, &
                              dloga, dloga_lo, dloga_hi, &
                              lo, hi, domlo, domhi, &
                              dx, dt)
       endif
#endif

    else
#ifdef RADIATION

#ifndef AMREX_USE_CUDA
       call amrex_error("ppm_type <=0 is not supported in with radiation")
#endif
#endif

       ! Compute all slopes

       allocate(dqx(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ))
       allocate(dqy(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ))

       ! Compute slopes
       if (plm_iorder == 1) then
          dqx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:NQ) = ZERO
          dqy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:NQ) = ZERO

       else if (plm_iorder == 2) then
          ! these are piecewise linear slopes.  The limiter is a 4th order
          ! limiter, but the overall method will be second order.

          do n = 1, NQ
             if (.not. reconstruct_state(n)) cycle
             call uslope(q, q_lo, q_hi, n, &
                         flatn, q_lo, q_hi, &
                         dqx, dqy, q_lo, q_hi, &
                         lo, hi)
          end do

          if (use_pslope == 1) then
             call pslope(q, q_lo, q_hi, &
                         flatn, q_lo, q_hi, &
                         dqx, dqy, q_lo, q_hi, &
                         srcQ, src_lo, src_hi, &
                         lo, hi, dx)

          endif

       elseif (plm_iorder == -2) then
          ! these are also piecewise linear, but it uses a multidimensional
          ! reconstruction based on the BDS advection method to construct
          ! the x- and y-slopes together
          do n = 1, NQ
             if (.not. reconstruct_state(n)) cycle
             call multid_slope(q, q_lo, q_hi, NQ, n, &
                               flatn, &
                               dqx, dqy, q_lo, q_hi, &
                               dx(1), dx(2), &
                               lo(1), lo(2), hi(1), hi(2))
          enddo
#ifndef AMREX_USE_CUDA
       else
          call amrex_error("ERROR: invalid value of islope")
#endif
       endif

       call trace_plm(1, q, q_lo, q_hi, &
                      qaux, qa_lo, qa_hi, &
                      dqx, q_lo, q_hi, &
                      qxm, qxp, qp_lo, qp_hi, &
                      dloga, dloga_lo, dloga_hi, &
                      lo, hi, domlo, domhi, &
                      dx, dt)

       call trace_plm(2, q, q_lo, q_hi, &
                      qaux, qa_lo, qa_hi, &
                      dqy, q_lo, q_hi, &
                      qym, qyp, qp_lo, qp_hi, &
                      dloga, dloga_lo, dloga_hi, &
                      lo, hi, domlo, domhi, &
                      dx, dt)

    endif

    if (ppm_type > 0) then
       deallocate(Ip,Im)
       deallocate(Ip_gc,Im_gc)
       deallocate(Ip_src,Im_src)
    else
       deallocate(dqx)
       deallocate(dqy)
    end if

    !-------------------------------------------------------------------------!
    ! Some notes on the work index (i.e., lo and hi arguments near the end    !
    !                               of the argument list).                    !
    ! * For cmpflx, we use face index in the flux direction and cell-centered !
    !   index for others.                                                     !
    ! * For trans*, we use cell-centered index of the valid region.           !
    !-------------------------------------------------------------------------!

    ! Compute F^x
    ! Inputs: qxm, qxp                     : xface, +-1 at y
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fx, ugdnvx, pgdnvx, gegdnvx : xface, +-1 at y
    call cmpflx(qxm, qxp, qp_lo, qp_hi, &
                fx, tflx_lo, tflx_hi, &
                qgdxtmp, q1_lo, q1_hi, &
#ifdef RADIATION
                rfx, tflx_lo, tflx_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, [lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0], &
                domlo, domhi)


    ! Compute F^y
    ! Inputs: qym, qyp                     : yface, +-1 at x
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fy, ugdnvy, pgdnvy, gegdnvy : yface, +-1 at x
    call cmpflx(qym, qyp, qp_lo, qp_hi, &
                fy, tfly_lo, tfly_hi, &
                q2, q2_lo, q2_hi, &
#ifdef RADIATION
                rfy, tfly_lo, tfly_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0], &
                domlo, domhi)

    ! add the transverse flux difference in y to the x states
    ! Inputs: qxm, qxp                     : xface, +-1 at y
    !         fy, ugdnvy, pgdnvy, gegdnvy  : yface, +-1 at x
    !         gamc                         : +-4
    ! Outputs: qm, qp                      : xface, +-0 at y
    call transy(qxm, qm, qxp, qp, qp_lo, qp_hi, &
                qaux, qa_lo, qa_hi, &
                fy, tfly_lo, tfly_hi, &
#ifdef RADIATION
                rfy, tfly_lo, tfly_hi, &
#endif
                q2, q2_lo, q2_hi, &
                srcQ, src_lo, src_hi, &
                hdt, hdtdy, &
                lo(1)-1, hi(1)+1, lo(2), hi(2))

    ! Solve the final Riemann problem across the x-interfaces with the
    ! full unsplit states.  The resulting flux through the x-interfaces
    ! is flux1
    call cmpflx(qm, qp, qp_lo, qp_hi, &
                flux1, f1_lo, f1_hi, &
                q1, q1_lo, q1_hi, &
#ifdef RADIATION
                rflux1, rf1_lo, rf1_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0], &
                domlo, domhi)

    ! add the transverse flux difference in x to the y states
    ! Inputs: qym, qyp                     : yface, +-1 at x
    !         fx, ugdnvx, pgdnvx, gegdnvx  : xface, +-1 at y
    !         gamc                         : +-4
    ! Outputs: qm, qp                      : yface, +-0 at x
    call transx(qym, qm, qyp, qp, qp_lo, qp_hi, &
                qaux, qa_lo, qa_hi, &
                fx, tflx_lo, tflx_hi, &
#ifdef RADIATION
                rfx, tflx_lo, tflx_hi, &
#endif
                qgdxtmp, q1_lo, q1_hi, &
                srcQ, src_lo, src_hi, &
                hdt, hdtdx, &
                area1, area1_lo, area1_hi, &
                vol, vol_lo, vol_hi, &
                lo(1), hi(1), lo(2)-1, hi(2)+1)

    ! Solve the final Riemann problem across the y-interfaces with the
    ! full unsplit states.  The resulting flux through the y-interfaces
    ! is flux2
    call cmpflx(qm, qp, qp_lo, qp_hi, &
                flux2, f2_lo, f2_hi, &
                q2, q2_lo, q2_hi, &
#ifdef RADIATION
                rflux2, rf2_lo, rf2_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, [lo(1), lo(2), 0], [hi(1), hi(2)+1, 0], &
                domlo, domhi)

    deallocate(qm,qp,qxm,qxp,qym,qyp)
    deallocate(fx,fy)
#ifdef RADIATION
    deallocate(rfx,rfy)
#endif
    deallocate(shk)
    deallocate(qgdxtmp)

  end subroutine umeth

  subroutine transx(qm, qmo, qp, qpo, qd_lo, qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    fx, fx_lo, fx_hi, &
#ifdef RADIATION
                    rfx, rfx_lo, rfx_hi, &
#endif
                    qgdx, qgdx_lo, qgdx_hi, &
                    srcQ, src_lo, src_hi, &
                    hdt, cdtdx,  &
                    area1, area1_lo, area1_hi, &
                    vol, vol_lo, vol_hi, &
                    ilo, ihi, jlo, jhi)

    use amrex_constants_module
    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NQAUX, &
                                 NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QGAME, &
                                 URHO, UMX, UMY, UEDEN, UEINT, QFS, QFX, &
                                 GDU, GDV, GDPRES, GDGAME, &
                                 NGDNV, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 GDERADS, QGAMCG, QLAMS, &
                                 fspace_type, comoving, &
#endif
                                 small_pres, small_temp, &
                                 npassive, qpass_map, upass_map, &
                                 transverse_use_eos, ppm_type, &
                                 transverse_reset_density, transverse_reset_rhoe, &
                                 ppm_predict_gammae
    use prob_params_module, only : mom_flux_has_p
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer qd_lo(3), qd_hi(3)
    integer qa_lo(3), qa_hi(3)
    integer fx_lo(3), fx_hi(3)
    integer qgdx_lo(3), qgdx_hi(3)
    integer src_lo(3), src_hi(3)
    integer area1_lo(3), area1_hi(3)
    integer vol_lo(3), vol_hi(3)
    integer ilo, ihi, jlo, jhi

#ifdef RADIATION
    integer rfx_lo(3), rfx_hi(3)
    real(rt)         rfx(rfx_lo(1):rfx_hi(1),rfx_lo(2):rfx_hi(2),0:ngroups-1)
#endif

    real(rt)         qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)

    real(rt)         fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),NVAR)
    real(rt)         qgdx(qgdx_lo(1):qgdx_hi(1),qgdx_lo(2):qgdx_hi(2),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),QVAR)
    real(rt)         area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
    real(rt)         vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))

    real(rt)         hdt, cdtdx

    integer          :: i, j, g
    integer          :: n, nqp, ipassive

    real(rt)         :: rr, rrnew, compo, compn
    real(rt)         :: rrr, rur, rvr, rer, ekinr, rhoekinr
    real(rt)         :: rrnewr, runewr, rvnewr, renewr
    real(rt)         :: rrl, rul, rvl, rel, ekinl, rhoekinl
    real(rt)         :: rrnewl, runewl, rvnewl, renewl

    ! here, pggp/pggm is the Godunov gas pressure (not radiation contribution)
    real(rt)         :: pggp, pggm, ugp, ugm, dAup, pav, uav, dAu, pnewl,pnewr
    real(rt)         :: geav, dge, gegp, gegm, gamc
    real(rt)         :: rhotmp

#ifdef RADIATION
    real(rt)         dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    real(rt)         eddf, f1, ugc, divu
#endif

    type (eos_t) :: eos_state

    logical :: reset_state

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    ! here cdtdx = 0.5 dt / dx

    ! update all of the passively-advected quantities in a single loop
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do j = jlo, jhi
          do i = ilo, ihi

             compn = hdt*(area1(i+1,j)*fx(i+1,j,n) - &
                          area1(i  ,j)*fx(i  ,j,n))/vol(i,j)

             if (j >= jlo+1) then
                rr = qp(i,j,  QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                  area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)
                compo = rr*qp(i,j,nqp) - compn
                qpo(i,j,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if

             if (j <= jhi-1) then
                rr = qm(i,j+1,QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                  area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)
                compo = rr*qm(i,j+1,nqp) - compn
                qmo(i,j+1,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if
          enddo
       enddo
    enddo

    ! hydro variables
    do j = jlo, jhi
       do i = ilo, ihi

          pggp = qgdx(i+1,j,GDPRES)
          pggm = qgdx(i,  j,GDPRES)
          ugp = qgdx(i+1,j,GDU)
          ugm = qgdx(i,  j,GDU)
          gegp = qgdx(i+1,j,GDGAME)
          gegm = qgdx(i,  j,GDGAME)

#ifdef RADIATION
          lambda(:) = qaux(i,j,QLAMS:QLAMS+ngroups-1)
          ugc = 0.5e0_rt*(ugp+ugm)
          ergp(:) = qgdx(i+1,j,GDERADS:GDERADS-1+ngroups)
          ergm(:) = qgdx(i  ,j,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dAup = area1(i+1,j)*pggp*ugp - area1(i,j)*pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          dAu = area1(i+1,j)*ugp-area1(i,j)*ugm
          geav = HALF*(gegp+gegm)
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,QGAMCG)
#else
          gamc = qaux(i,j,QGAMC)
#endif

#ifdef RADIATION
          lamge(:) = lambda(:) * (ergp(:)-ergm(:))
          luge(:) = ugc * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             divu = (area1(i+1,j)*ugp-area1(i,j)*ugm)/vol(i,j)
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = -hdt * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * divu
             end do
          else ! mixed frame
             der(:) = cdtdx * luge(:)
          end if
#endif

          !-------------------------------------------------------------------
          ! qp state
          !-------------------------------------------------------------------

          ! "right" state on the j-1/2 interface
          if (j >= jlo+1) then

             ! Convert to conservation form
             rrr = qp(i,j,QRHO)
             rur = rrr*qp(i,j,QU)
             rvr = rrr*qp(i,j,QV)
             ekinr = HALF*rrr*sum(qp(i,j,QU:QW)**2)
             rer = qp(i,j,QREINT) + ekinr
#ifdef RADIATION
             err(:) = qp(i,j,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                 area1(i,j)*fx(i,j,URHO))/vol(i,j)

             ! Note that pressure may be treated specially here, depending on 
             ! the geometry.  Our y-interface equation for (rho u) is:
             !
             !  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
             !
             ! in cylindrical coords -- note that the p term is not in
             ! a divergence, so there are no area factors.  For this
             ! geometry, we do not include p in our definition of the
             ! flux in the x-direction, for we need to fix this now.
             runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                 area1(i,j)*fx(i,j,UMX))/vol(i,j)
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewr = runewr -cdtdx *(pggp-pggm)
             endif
             rvnewr = rvr - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                                 area1(i,j)*fx(i,j,UMY))/vol(i,j)
             renewr = rer - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                                 area1(i,j)*fx(i,j,UEDEN))/vol(i,j)

#ifdef RADIATION
             runewr = runewr - HALF*hdt*(area1(i+1,j)+area1(i,j))*sum(lamge)/vol(i,j)
             renewr = renewr + dre
             ernewr(:) = err(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
                  area1(i,j)*rfx(i,j,:))/vol(i,j) &
                  + der(:)
#endif

             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                renewr = rer
#ifdef RADIATION
                ernewr(:) = err(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             rhotmp = rrnewr
             qpo(i,j,QRHO) = rhotmp
             qpo(i,j,QU  ) = runewr/rhotmp
             qpo(i,j,QV  ) = rvnewr/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qpo(i,j,QRHO) = qpo(i,j,QRHO) + hdt*srcQ(i,j,QRHO)
                qpo(i,j,QU:QV) = qpo(i,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekinr = HALF*(runewr**2+rvnewr**2+(rhotmp*qpo(i,j,QW))**2)/rhotmp
             qpo(i,j,QREINT) = renewr - rhoekinr

             if (ppm_type == 0) then
                qpo(i,j,QREINT) = qpo(i,j,QREINT) + hdt*srcQ(i,j,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qpo(i,j,QREINT) = qp(i,j,QREINT) - &
                        hdt*(area1(i+1,j)*fx(i+1,j,UEINT)-  &
                             area1(i,j)*fx(i,j,UEINT) + pav*dAu)/vol(i,j)

                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,QREINT) < ZERO .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qpo(i,j,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qpo(i,j,QREINT) = qpo(i,j,QRHO)*eos_state % e
                      qpo(i,j,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % e   = qpo(i,j,QREINT) / qpo(i,j,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qpo(i,j,QFS:QFS+nspec-1)
                      eos_state % aux = qpo(i,j,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qpo(i,j,QPRES ) = pnewr
                      qpo(i,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      ! we are expressing the pressure evolution as:
                      !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
                      ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
                      ! but these are divergences, so we need area factors
                      pnewr = qp(i,j,QPRES) - hdt*(dAup + pav*dAu*(gamc - ONE))/vol(i,j)
                      qpo(i,j,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qpo(i,j,QGAME) = qp(i,j,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*dAu)/vol(i,j) - cdtdx*uav*dge

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,QPRES) = qpo(i,j,QREINT)*(qpo(i,j,QGAME)-ONE)

                endif
             else
                qpo(i,j,QPRES) = qp(i,j,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qpo(i,j,QGAME) = qp(i,j,QGAME)
             endif

#ifdef RADIATION
             qpo(i,j,qrad:qradhi) = ernewr(:)
             qpo(i,j,qptot)   = sum(lambda*ernewr) + qpo(i,j,QPRES)
             qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qm state
          !-------------------------------------------------------------------

          ! "left" state on the j+1/2 interface
          if (j <= jhi-1) then

             rrl = qm(i,j+1,QRHO)
             rul = rrl*qm(i,j+1,QU)
             rvl = rrl*qm(i,j+1,QV)
             ekinl = HALF*rrl*sum(qm(i,j+1,QU:QW)**2)
             rel = qm(i,j+1,QREINT) + ekinl
#ifdef RADIATION
             erl(:) = qm(i,j+1,qrad:qradhi)
#endif

             rrnewl = rrl - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                 area1(i,j)*fx(i,j,URHO))/vol(i,j)
             runewl = rul - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                 area1(i,j)*fx(i,j,UMX))/vol(i,j) 
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewl = runewl -cdtdx *(pggp-pggm)
             endif
             rvnewl = rvl - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                                 area1(i,j)*fx(i,j,UMY))/vol(i,j)
             renewl = rel - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                                 area1(i,j)*fx(i,j,UEDEN))/vol(i,j)


#ifdef RADIATION
             runewl = runewl - HALF*hdt*(area1(i+1,j)+area1(i,j))*sum(lamge)/vol(i,j)
             renewl = renewl + dre
             ernewl(:) = erl(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
                  area1(i,j)*rfx(i,j,:))/vol(i,j) &
                  + der(:)
#endif
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                renewl = rel
#ifdef RADIATION
                ernewl(:) = erl(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             rhotmp = rrnewl
             qmo(i,j+1,QRHO) = rhotmp
             qmo(i,j+1,QU  ) = runewl/rhotmp
             qmo(i,j+1,QV  ) = rvnewl/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qmo(i,j+1,QRHO)  = qmo(i,j+1,QRHO) + hdt*srcQ(i,j,QRHO)
                qmo(i,j+1,QU:QV) = qmo(i,j+1,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinl = HALF*(runewl**2+rvnewl**2+(rhotmp*qmo(i,j+1,QW))**2)/rhotmp
             qmo(i,j+1,QREINT) = renewl - rhoekinl
             if (ppm_type == 0) then
                qmo(i,j+1,QREINT) = qmo(i,j+1,QREINT) + hdt*srcQ(i,j,QREINT)
             endif


             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j+1,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using
                   ! the discretized expression for updating (rho e).

                   qmo(i,j+1,QREINT) = qm(i,j+1,QREINT) - &
                        hdt*(area1(i+1,j)*fx(i+1,j,UEINT)-  &
                             area1(i,j)*fx(i,j,UEINT) + pav*dAu)/vol(i,j)

                   ! if we are still negative, then we need to reset
                   if (qmo(i,j+1,QREINT) < ZERO .and. qmo(i,j+1,QRHO) > ZERO) then
                      eos_state % rho = qmo(i,j+1,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i,j+1,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qmo(i,j+1,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qmo(i,j+1,QREINT) = qmo(i,j+1,QRHO)*eos_state % e
                      qmo(i,j+1,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qmo(i,j+1,QRHO) > ZERO) then
                      eos_state % rho = qmo(i,j+1,QRHO)
                      eos_state % e   = qmo(i,j+1,QREINT) / qmo(i,j+1,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qmo(i,j+1,QFS:QFS+nspec-1)
                      eos_state % aux = qmo(i,j+1,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qmo(i,j+1,QPRES ) = pnewr
                      qmo(i,j+1,QREINT) = eos_state % e * eos_state % rho
                   else
                      ! we are expressing the pressure evolution as:
                      !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
                      ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
                      ! but these are divergences, so we need area factors
                      pnewl = qm(i,j+1,QPRES) - hdt*(dAup + pav*dAu*(gamc - ONE))/vol(i,j)
                      qmo(i,j+1,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i,j+1,QPRES) = qmo(i,j+1,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qmo(i,j+1,QPRES) = max(qmo(i,j+1,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qmo(i,j+1,QGAME) = qm(i,j+1,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*dAu)/vol(i,j) - cdtdx*uav*dge

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j+1,QPRES) = qmo(i,j+1,QREINT)*(qmo(i,j+1,QGAME)-ONE)

                endif
             else
                qmo(i,j+1,QPRES) = qm(i,j+1,QPRES)
                if (ppm_type == 0) then
                   qmo(i,j+1,QPRES) = qmo(i,j+1,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qmo(i,j+1,QGAME) = qm(i,j+1,QGAME)
             endif

#ifdef RADIATION
             qmo(i,j+1,qrad:qradhi) = ernewl(:)
             qmo(i,j+1,qptot)   = sum(lambda*ernewl) + qmo(i,j+1,QPRES)
             qmo(i,j+1,qreitot) = sum(qmo(i,j+1,qrad:qradhi)) + qmo(i,j+1,QREINT)
#endif

          end if

       enddo
    enddo

  end subroutine transx


  subroutine transy(qm, qmo, qp, qpo, qd_lo, qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    fy, fy_lo, fy_hi, &
#ifdef RADIATION
                    rfy, rfy_lo, rfy_hi, &
#endif
                    qgdy, qgdy_lo, qgdy_hi, &
                    srcQ, src_lo, src_hi, &
                    hdt, cdtdy, ilo, ihi, jlo, jhi)

    use amrex_constants_module
    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NQAUX, &
                                 NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QGAME, &
                                 URHO, UMX, UMY, UEDEN, UEINT, QFS, QFX, &
                                 GDU, GDV, GDPRES, GDGAME, &
                                 NGDNV, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 GDERADS, QGAMCG, QLAMS, &
                                 fspace_type, comoving, &
#endif
                                 small_pres, small_temp, &
                                 npassive, qpass_map, upass_map, &
                                 transverse_use_eos, ppm_type, &
                                 transverse_reset_density, transverse_reset_rhoe, &
                                 ppm_predict_gammae
    use prob_params_module, only : mom_flux_has_p
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer qd_lo(3), qd_hi(3)
    integer qa_lo(3), qa_hi(3)
    integer fy_lo(3), fy_hi(3)
    integer qgdy_lo(3), qgdy_hi(3)
    integer src_lo(3), src_hi(3)
    integer ilo, ihi, jlo, jhi

#ifdef RADIATION
    integer rfy_lo(3), rfy_hi(3)
    real(rt)         rfy(rfy_lo(1):rfy_hi(1),rfy_lo(2):rfy_hi(2),0:ngroups-1)
#endif

    real(rt)         qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)

    real(rt)         fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),NVAR)
    real(rt)         qgdy(qgdy_lo(1):qgdy_hi(1),qgdy_lo(2):qgdy_hi(2),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),QVAR)

    real(rt)         hdt, cdtdy

    integer          :: i, j, g
    integer          :: n, nqp, ipassive

    real(rt)         :: rr,rrnew
    real(rt)         :: pggp, pggm, ugp, ugm, dup, pav, uav, du, pnewr,pnewl
    real(rt)         :: gegp, gegm, geav, dge, gamc
    real(rt)         :: rrr, rur, rvr, rer, ekinr, rhoekinr
    real(rt)         :: rrnewr, runewr, rvnewr, renewr
    real(rt)         :: rrl, rul, rvl, rel, ekinl, rhoekinl
    real(rt)         :: rrnewl, runewl, rvnewl, renewl
    real(rt)         :: rhotmp
    real(rt)         :: compo, compn

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    real(rt)         :: eddf, f1, ugc
#endif

    type (eos_t) :: eos_state

    logical :: reset_state

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    ! update all of the passively-advected quantities in a single loop
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdy*(fy(i,j+1,n)-fy(i,j,n))

             if (i >= ilo+1) then
                rr = qp(i,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
                compo = rr*qp(i,j,nqp) - compn
                qpo(i,j,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if

             if (i <= ihi-1) then
                rr = qm(i+1,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
                compo = rr*qm(i+1,j,nqp) - compn
                qmo(i+1,j,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if
          enddo
       enddo
    enddo

    ! hydro variables
    do j = jlo, jhi
       do i = ilo, ihi

          pggp = qgdy(i,j+1,GDPRES)
          pggm = qgdy(i,j  ,GDPRES)
          ugp = qgdy(i,j+1,GDV)
          ugm = qgdy(i,j  ,GDV)
          gegp = qgdy(i,j+1,GDGAME)
          gegm = qgdy(i,j  ,GDGAME)

#ifdef RADIATION
          lambda(:) = qaux(i,j,QLAMS:QLAMS+ngroups-1)
          ugc = 0.5e0_rt*(ugp+ugm)
          ergp(:) = qgdy(i,j+1,GDERADS:GDERADS-1+ngroups)
          ergm(:) = qgdy(i,j  ,GDERADS:GDERADS-1+ngroups)
#endif
          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          du = ugp-ugm
          geav = HALF*(gegp+gegm)
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,QGAMCG)
#else
          gamc = qaux(i,j,QGAMC)
#endif


#ifdef RADIATION
          lamge(:) = lambda(:) * (ergp(:)-ergm(:))
          luge(:) = ugc * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = cdtdy * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qp state
          !-------------------------------------------------------------------

          ! right state on the i-1/2 interface
          if (i >= ilo+1) then

             ! Convert to conservation form
             rrr = qp(i,j,QRHO)
             rur = rrr*qp(i,j,QU)
             rvr = rrr*qp(i,j,QV)
             ekinr = HALF*rrr*sum(qp(i,j,QU:QW)**2)
             rer = qp(i,j,QREINT) + ekinr
#ifdef RADIATION
             err(:) = qp(i,j,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))

             runewr = rur - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
             ! note: we are always Cartesian in the y-direction, so the
             ! pressure term is already accounted for in the flux
             rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) 
             renewr = rer - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN))

#ifdef RADIATION
             rvnewr = rvnewr - cdtdy*sum(lamge)
             renewr = renewr + dre
             ernewr(:) = err(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
                  + der(:)
#endif

             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr <= ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                renewr = rer
#ifdef RADIATION
                ernewr(:) = err(:)
#endif
                reset_state = .true.
             end if

             ! convert back to non-conservation form
             rhotmp =  rrnewr
             qpo(i,j,QRHO  ) = rhotmp
             qpo(i,j,QU    ) = runewr/rhotmp
             qpo(i,j,QV    ) = rvnewr/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qpo(i,j,QRHO  ) = qpo(i,j,QRHO  ) + hdt*srcQ(i,j,QRHO)
                qpo(i,j,QU:QV) = qpo(i,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinr = HALF*(runewr**2+rvnewr**2+(rhotmp*qpo(i,j,QW))**2)/rhotmp
             qpo(i,j,QREINT) = renewr - rhoekinr
             if (ppm_type == 0) then
                qpo(i,j,QREINT) = qpo(i,j,QREINT) + hdt*srcQ(i,j,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,QREINT) <= ZERO) then
                   qpo(i,j,QREINT) = qp(i,j,QREINT) - &
                        cdtdy*(fy(i,j+1,UEINT)- fy(i,j,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,QREINT) < ZERO .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qpo(i,j,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qpo(i,j,QREINT) = qpo(i,j,QRHO) * eos_state % e
                      qpo(i,j,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % e   = qpo(i,j,QREINT) / qpo(i,j,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qpo(i,j,QFS:QFS+nspec-1)
                      eos_state % aux = qpo(i,j,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qpo(i,j,QPRES ) = pnewr
                      qpo(i,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      pnewr = qp(i  ,j,QPRES)-cdtdy*(dup + pav*du*(gamc - ONE))
                      qpo(i,j,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qpo(i,j,QGAME) = qp(i,j,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,QPRES) = qpo(i,j,QREINT)*(qpo(i,j,QGAME)-ONE)

                endif
             else
                qpo(i,j,QPRES) = qp(i,j,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qpo(i,j,QGAME) = qp(i,j,QGAME)
             endif

#ifdef RADIATION
             qpo(i,j,qrad:qradhi) = ernewr(:)
             qpo(i,j,qptot  ) = sum(lambda*ernewr) + qpo(i,j,QPRES)
             qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qm state
          !-------------------------------------------------------------------

          ! left state on the i+1/2 interface
          if (i <= ihi-1) then

             rrl = qm(i+1,j,QRHO)
             rul = rrl*qm(i+1,j,QU)
             rvl = rrl*qm(i+1,j,QV)
             ekinl = HALF*rrl*sum(qm(i+1,j,QU:QW)**2)
             rel = qm(i+1,j,QREINT) + ekinl
#ifdef RADIATION
             erl(:) = qm(i+1,j,qrad:qradhi)
#endif

             rrnewl = rrl - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))
             runewl = rul - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
             rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) 
             renewl = rel - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN))

#ifdef RADIATION
             rvnewl = rvnewl - cdtdy*sum(lamge)
             renewl = renewl + dre
             ernewl(:) = erl(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
                  + der(:)
#endif

             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl <= ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                renewl = rel
#ifdef RADIATION
                ernewl(:) = erl(:)
#endif
                reset_state = .true.
             endif

             rhotmp =  rrnewl
             qmo(i+1,j,QRHO  ) = rhotmp            
             qmo(i+1,j,QU    ) = runewl/rhotmp
             qmo(i+1,j,QV    ) = rvnewl/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qmo(i+1,j,QRHO  ) = qmo(i+1,j,QRHO  ) + hdt*srcQ(i,j,QRHO)
                qmo(i+1,j,QU:QV) = qmo(i+1,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinl = HALF*(runewl**2+rvnewl**2+(rhotmp*qmo(i+1,j,QW))**2)/rhotmp
             qmo(i+1,j,QREINT) = renewl - rhoekinl
             if (ppm_type == 0) then
                qmo(i+1,j,QREINT) = qmo(i+1,j,QREINT) + hdt*srcQ(i,j,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i+1,j,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).

                   qmo(i+1,j,QREINT) = qm(i+1,j,QREINT) - &
                        cdtdy*(fy(i,j+1,UEINT) - fy(i,j,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qmo(i+1,j,QREINT) < ZERO .and. qmo(i+1,j,QRHO) > ZERO) then
                      eos_state % rho = qmo(i+1,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i+1,j,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qmo(i+1,j,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qmo(i+1,j,QREINT) = qmo(i+1,j,QRHO)*eos_state % e
                      qmo(i+1,j,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qmo(i+1,j,QRHO) > ZERO) then
                      eos_state % rho = qmo(i+1,j,QRHO)
                      eos_state % e   = qmo(i+1,j,QREINT) / qmo(i+1,j,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qmo(i+1,j,QFS:QFS+nspec-1)
                      eos_state % aux = qmo(i+1,j,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qmo(i+1,j,QPRES ) = pnewr
                      qmo(i+1,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      pnewl = qm(i+1,j,QPRES)-cdtdy*(dup + pav*du*(gamc - ONE))
                      qmo(i+1,j,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i+1,j,QPRES) = qmo(i+1,j,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qmo(i+1,j,QPRES) = max(qmo(i+1,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qmo(i+1,j,QGAME) = qm(i+1,j,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qmo(i+1,j,QPRES) = qmo(i+1,j,QREINT)*(qmo(i+1,j,QGAME)-ONE)

                endif
             else
                qmo(i+1,j,QPRES) = qm(i+1,j,QPRES)
                if (ppm_type == 0) then
                   qmo(i+1,j,QPRES) = qmo(i+1,j,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qmo(i+1,j,QGAME) = qm(i+1,j,QGAME)
             endif

#ifdef RADIATION
             qmo(i+1,j,qrad:qradhi) = ernewl(:)
             qmo(i+1,j,qptot  ) = sum(lambda*ernewl) + qmo(i+1,j,QPRES)
             qmo(i+1,j,qreitot) = sum(qmo(i+1,j,qrad:qradhi)) + qmo(i+1,j,QREINT)
#endif

          end if

       enddo
    enddo

  end subroutine transy

end module ctu_advection_module
