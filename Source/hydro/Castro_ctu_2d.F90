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
! ::: :: qaux        => (const)  auxiliary hydro info
! ::: :: flatn       => (const)  flattening parameter
! ::: :: srcQ        => (const)  primitive variable source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: dx          => (const)  grid spacing in all 3 directions (array)
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: ----------------------------------------------------------------

  subroutine umeth(q, q_lo, q_hi, &
                   flatn, &
                   qaux, qa_lo, qa_hi, &
                   srcQ, src_lo, src_hi, &
                   lo, hi, dx, dt, &
                   uout, uout_lo, uout_hi, &
                   flux1, fd1_lo, fd1_hi, &
                   flux2, fd2_lo, fd2_hi, &
#ifdef RADIATION
                   rflux1, rfd1_lo, rfd1_hi, &
                   rflux2, rfd2_lo, rfd2_hi, &
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
                                   GDU, GDV, GDPRES, NGDNV, NQ, &
                                   NQAUX, &
                                   ppm_type, &
                                   use_pslope, plm_iorder, ppm_temp_fix
    use trace_module, only : tracexy
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : tracexy_ppm_rad
#else
    use trace_ppm_module, only : tracexy_ppm
#endif
    use transverse_module, only : transx, transy
    use riemann_module, only: cmpflx
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
    use multid_slope_module, only : multid_slope

    implicit none

    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: fd1_lo(3), fd1_hi(3)
    integer, intent(in) :: fd2_lo(3), fd2_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rfd1_lo(3), rfd1_hi(3)
    integer, intent(in) :: rfd2_lo(3), rfd2_hi(3)
#endif
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt)        , intent(in) :: dx(3), dt
    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ)
    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
    real(rt)        , intent(in) :: flatn(q_lo(1):q_hi(1),q_lo(2):q_hi(2))
    real(rt)        , intent(in) :: srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),QVAR)
    real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2))
    real(rt)        , intent(inout) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),NGDNV)
    real(rt)        , intent(inout) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),NGDNV)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),NVAR)
    real(rt)        , intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),NVAR)
    real(rt)        , intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),NVAR)
#ifdef RADIATION
    real(rt)        , intent(inout) :: rflux1(rfd1_lo(1):rfd1_hi(1),rfd1_lo(2):rfd1_hi(2),0:ngroups-1)
    real(rt)        , intent(inout) :: rflux2(rfd2_lo(1):rfd2_hi(1),rfd2_lo(2):rfd2_hi(2),0:ngroups-1)
#endif
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2))
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))

    ! Left and right state arrays (edge centered, cell centered)
    real(rt)        , allocatable::  qm(:,:,:),  qp(:,:,:)
    real(rt)        , allocatable:: qxm(:,:,:), qym(:,:,:)
    real(rt)        , allocatable:: qxp(:,:,:), qyp(:,:,:)

    ! Work arrays to hold riemann state and conservative fluxes

    real(rt)        , allocatable ::  fx(:,:,:),  fy(:,:,:)
#ifdef RADIATION
    real(rt)        , allocatable ::  rfx(:,:,:),  rfy(:,:,:)
#endif
    real(rt)        , allocatable ::  qgdxtmp(:,:,:)
    real(rt)        , allocatable :: shk(:,:)

    ! Local scalar variables
    real(rt)         :: dtdx
    real(rt)         :: hdtdx, hdt, hdtdy
    integer          :: i, j, idim, iwave, n

    integer :: tflx_lo(3), tflx_hi(3)
    integer :: tfly_lo(3), tfly_hi(3)
    integer :: shk_lo(3), shk_hi(3)
    integer :: qp_lo(3), qp_hi(3)

    real(rt)        , allocatable :: Ip(:,:,:,:,:)
    real(rt)        , allocatable :: Im(:,:,:,:,:)

    real(rt)        , allocatable :: Ip_src(:,:,:,:,:)
    real(rt)        , allocatable :: Im_src(:,:,:,:,:)

    ! gamma_c/1 on the interfaces
    real(rt)        , allocatable :: Ip_gc(:,:,:,:,:)
    real(rt)        , allocatable :: Im_gc(:,:,:,:,:)

    ! temporary interface values of the parabola
    real(rt), allocatable :: sxm(:,:), sxp(:,:), sym(:,:), syp(:,:)

    real(rt)        , allocatable :: dqx(:,:,:), dqy(:,:,:)

    integer :: I_lo(3), I_hi(3)

    type(eos_t) :: eos_state

    logical :: source_nonzero(QVAR)

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

    if (ppm_type == 0) then

       allocate(dqx(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ))
       allocate(dqy(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ))

       ! Compute slopes
       if (plm_iorder == 1) then
          dqx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:NQ) = ZERO
          dqy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:NQ) = ZERO

       elseif (plm_iorder == 2) then
          ! these are piecewise linear slopes.  The limiter is a 4th order
          ! limiter, but the overall method will be second order.
          call uslope(q, flatn, q_lo, q_hi, &
                      dqx, dqy, dqx, q_lo, q_hi, &  ! second dqx is dummy
                      lo(1), lo(2), hi(1), hi(2), 0, 0)

          if (use_pslope == 1) then
             call pslope(q, flatn, q_lo, q_hi, &
                         dqx, dqy, dqx, q_lo, q_hi, &  ! second dqx is dummy
                         srcQ, src_lo, src_hi, &
                         lo(1), lo(2), hi(1), hi(2), 0, 0, dx)

          endif

       elseif (plm_iorder == -2) then
          ! these are also piecewise linear, but it uses a multidimensional
          ! reconstruction based on the BDS advection method to construct
          ! the x- and y-slopes together
          do n = 1, NQ
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

    else

       ! preprocess the sources -- we don't want to trace under a source that is empty
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
          call ppm_reconstruct(q, q_lo, q_hi, NQ, n, &
                               flatn, q_lo, q_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               lo(1), lo(2), hi(1), hi(2), dx, 0, 0)

          call ppm_int_profile(q, q_lo, q_hi, NQ, n, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               Ip, Im, I_lo, I_hi, NQ, n, &
                               lo(1), lo(2), hi(1), hi(2), dx, dt, 0, 0)
       end do

       ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
       ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
       if (ppm_temp_fix == 1) then
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                do idim = 1, 2
                   do iwave = 1, 3
                      eos_state%rho   = Ip(i,j,idim,iwave,QRHO)
                      eos_state%T     = Ip(i,j,idim,iwave,QTEMP)
                      eos_state%xn(:) = Ip(i,j,idim,iwave,QFS:QFS-1+nspec)
                      eos_state%aux   = Ip(i,j,idim,iwave,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      Ip(i,j,idim,iwave,QPRES) = eos_state%p
                      Ip(i,j,idim,iwave,QREINT) = Ip(i,j,idim,iwave,QRHO)*eos_state%e
                      Ip_gc(i,j,idim,iwave,1) = eos_state%gam1

                      eos_state%rho   = Im(i,j,idim,iwave,QRHO)
                      eos_state%T     = Im(i,j,idim,iwave,QTEMP)
                      eos_state%xn(:) = Im(i,j,idim,iwave,QFS:QFS-1+nspec)
                      eos_state%aux   = Im(i,j,idim,iwave,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      Im(i,j,idim,iwave,QPRES) = eos_state%p
                      Im(i,j,idim,iwave,QREINT) = Im(i,j,idim,iwave,QRHO)*eos_state%e
                      Im_gc(i,j,idim,iwave,1) = eos_state%gam1
                   enddo
                enddo
             enddo
          enddo

       endif

       ! get an edge-based gam1 here if we didn't get it from the EOS
       ! call above (for ppm_temp_fix = 1)
       if (ppm_temp_fix /= 1) then
          call ppm_reconstruct(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               flatn, q_lo, q_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               lo(1), lo(2), hi(1), hi(2), dx, 0, 0)
          
          call ppm_int_profile(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, q_lo, q_hi, &
                               Ip_gc, Im_gc, I_lo, I_hi, 1, 1, &
                               lo(1), lo(2), hi(1), hi(2), dx, dt, 0, 0)
       endif

       do n = 1, QVAR
          if (source_nonzero(n)) then
             call ppm_reconstruct(srcQ, src_lo, src_hi, QVAR, n, &
                                  flatn, q_lo, q_hi, &
                                  sxm, sxp, sym, syp, q_lo, q_hi, &
                                  lo(1), lo(2), hi(1), hi(2), dx, 0, 0)

             call ppm_int_profile(srcQ, src_lo, src_hi, QVAR, n, &
                                  q, q_lo, q_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, q_lo, q_hi, &
                                  Ip_src, Im_src, I_lo, I_hi, QVAR, n, &
                                  lo(1), lo(2), hi(1), hi(2), dx, dt, 0, 0)
          else
             Ip_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),:,:,n) = ZERO
             Im_src(I_lo(1):I_hi(1),I_lo(2):I_hi(2),:,:,n) = ZERO
          endif
       enddo

       deallocate(sxm, sxp, sym, syp)
    endif

    ! Trace to edges w/o transverse flux correction terms.  Here,
    !      qxm and qxp will be the states on either side of the x interfaces
    ! and  qym and qyp will be the states on either side of the y interfaces
    if (ppm_type .eq. 0) then
#ifdef RADIATION
#ifndef AMREX_USE_CUDA
       call amrex_error("ppm_type <=0 is not supported in umeth for radiation")
#endif
#else
       call tracexy(q, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    dqx, dqy, q_lo, q_hi, &
                    qxm, qxp, qym, qyp, qp_lo, qp_hi, &
                    dloga, dloga_lo, dloga_hi, &
                    lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                    dx, dt, 0, 0)
#endif
    else
#ifdef RADIATION
       call tracexy_ppm_rad(q, q_lo, q_hi, &
                            qaux, qa_lo, qa_hi, &
                            Ip, Im, Ip_src, Im_src, I_lo, I_hi, &            
                            qxm, qxp, qym, qyp, qp_lo, qp_hi, &
                            dloga, dloga_lo, dloga_hi, &
                            lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                            dx, dt, 0, 0)
#else
       call tracexy_ppm(q, q_lo, q_hi, &
                        qaux, qa_lo, qa_hi, &
                        Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, I_lo, I_hi, &
                        qxm, qxp, qym, qyp, qp_lo, qp_hi, &
                        dloga, dloga_lo, dloga_hi, &
                        lo(1), lo(2), hi(1), hi(2), domlo, domhi, &
                        dx, dt, 0, 0)
#endif


       deallocate(Ip,Im)
       deallocate(Ip_gc,Im_gc)
       deallocate(Ip_src,Im_src)

    end if

    ! Solve the Riemann problem in the x-direction using these first
    ! guesses for the x-interface states.  This produces the flux fx
    call cmpflx(qxm, qxp, qp_lo, qp_hi, &
                fx, tflx_lo, tflx_hi, &
                qgdxtmp, q1_lo, q1_hi, &
#ifdef RADIATION
                rfx, tflx_lo, tflx_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, lo(1), hi(1)+1, lo(2)-1, hi(2)+1, 0, 0, 0, &
                domlo, domhi)

    ! Solve the Riemann problem in the y-direction using these first
    ! guesses for the y-interface states.  This produces the flux fy
    call cmpflx(qym, qyp, qp_lo, qp_hi, &
                fy, tfly_lo, tfly_hi, &
                q2, q2_lo, q2_hi, &
#ifdef RADIATION
                rfy, tfly_lo, tfly_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, lo(1)-1, hi(1)+1, lo(2), hi(2)+1, 0, 0, 0, &
                domlo, domhi)

    ! Correct the x-interface states (qxm, qxp) by adding the
    ! transverse flux difference in the y-direction to the x-interface
    ! states.  This results in the new x-interface states qm and qp
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
                flux1, fd1_lo, fd1_hi, &
                q1, q1_lo, q1_hi, &
#ifdef RADIATION
                rflux1, rfd1_lo, rfd1_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                1, lo(1), hi(1)+1, lo(2), hi(2), 0, 0, 0, &
                domlo, domhi)

    ! Correct the y-interface states (qym, qyp) by adding the
    ! transverse flux difference in the x-direction to the y-interface
    ! states.  This results in the new y-interface states qm and qp
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
                flux2, fd2_lo, fd2_hi, &
                q2, q2_lo, q2_hi, &
#ifdef RADIATION
                rflux2, rfd2_lo, rfd2_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, shk_lo, shk_hi, &
                2, lo(1), hi(1), lo(2), hi(2)+1, 0, 0, 0, &
                domlo, domhi)

    deallocate(qm,qp,qxm,qxp,qym,qyp)
    deallocate(fx,fy)
#ifdef RADIATION
    deallocate(rfx,rfy)
#endif
    deallocate(shk)
    deallocate(qgdxtmp)

  end subroutine umeth

! :::
! ::: ------------------------------------------------------------------
! :::


end module ctu_advection_module
