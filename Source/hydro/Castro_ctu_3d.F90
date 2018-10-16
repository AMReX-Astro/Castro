module ctu_advection_module

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public umeth

contains

! ::: ---------------------------------------------------------------
! ::: :: UMETH     Compute hyperbolic fluxes using unsplit second
! ::: ::           order Godunov integrator.
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: qaux        => (const)  auxiliary hydro data
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: nz          => (const)  number of cells in Z direction
! ::: :: dx          => (const)  grid spacing in X, Y, Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
! ::: ----------------------------------------------------------------

  !! TODO: we can get rid of the the different temporary q Godunov
  !! state arrays

  subroutine umeth(q, qd_lo, qd_hi, &
                   flatn, &
                   qaux, qa_lo, qa_hi, &
                   srcQ, src_lo, src_hi, &
                   lo, hi, dx, dt, &
                   uout, uout_lo, uout_hi, &
                   flux1, f1_lo, f1_hi, &
                   flux2, f2_lo, f2_hi, &
                   flux3, f3_lo, f3_hi, &
#ifdef RADIATION
                   rflux1, rf1_lo, rf1_hi, &
                   rflux2, rf2_lo, rf2_hi, &
                   rflux3, rf3_lo, rf3_hi, &
#endif
                   q1, q1_lo, q1_hi, &
                   q2, q2_lo, q2_hi, &
                   q3, q3_lo, q3_hi, &
                   domlo, domhi)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NQ, NVAR, QPRES, QRHO, QU, QW, &
                                   QFS, QFX, QTEMP, QREINT, &
                                   QC, QGAMC, NQAUX, &
                                   NGDNV, GDU, GDV, GDW, GDPRES, &
                                   ppm_type, &
                                   use_pslope, ppm_temp_fix, &
                                   hybrid_riemann
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
    use trace_module, only : tracexy, tracez
    use transverse_module
    use ppm_module, only : ppm_reconstruct, ppm_int_profile
    use slope_module, only : uslope, pslope
    use actual_network, only : nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use riemann_module, only: cmpflx
    use amrex_constants_module
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use trace_ppm_rad_module, only : tracexy_ppm_rad, tracez_ppm_rad
#else
    use trace_ppm_module, only : tracexy_ppm, tracez_ppm
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use advection_util_module, only : shock
    use prob_params_module, only : dg

    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: f1_lo(3), f1_hi(3)
    integer, intent(in) :: f2_lo(3), f2_hi(3)
    integer, intent(in) :: f3_lo(3), f3_hi(3)
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)
#ifdef RADIATION
    integer, intent(in) :: rf1_lo(3), rf1_hi(3)
    integer, intent(in) :: rf2_lo(3), rf2_hi(3)
    integer, intent(in) :: rf3_lo(3), rf3_hi(3)
#endif

    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    real(rt), intent(in) ::  srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)

    real(rt), intent(inout) ::  uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(inout) :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)
    real(rt), intent(inout) ::    q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    real(rt), intent(inout) ::    q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt), intent(inout) ::    q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    real(rt), intent(in) :: dx(3), dt

#ifdef RADIATION
    real(rt), intent(inout) :: rflux1(rf1_lo(1):rf1_hi(1),rf1_lo(2):rf1_hi(2),rf1_lo(3):rf1_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: rflux2(rf2_lo(1):rf2_hi(1),rf2_lo(2):rf2_hi(2),rf2_lo(3):rf2_hi(3),0:ngroups-1)
    real(rt), intent(inout) :: rflux3(rf3_lo(1):rf3_hi(1),rf3_lo(2):rf3_hi(2),rf3_lo(3):rf3_hi(3),0:ngroups-1)
#endif

    real(rt) :: dxinv, dyinv, dzinv
    real(rt) :: dtdx, dtdy, dtdz, hdt
    real(rt) :: cdtdx, cdtdy, cdtdz
    real(rt) :: hdtdx, hdtdy, hdtdz

    integer :: n
    integer :: i, j, k, iwave, idim, d

    real(rt), pointer :: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)

    real(rt), pointer :: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
    real(rt), pointer :: Ip_src(:,:,:,:,:,:), Im_src(:,:,:,:,:,:)
    real(rt), pointer :: Ip_gc(:,:,:,:,:,:), Im_gc(:,:,:,:,:,:)

    real(rt)        , pointer :: shk(:,:,:)

    real(rt)        , pointer :: sxm(:,:,:), sym(:,:,:), szm(:,:,:)
    real(rt)        , pointer :: sxp(:,:,:), syp(:,:,:), szp(:,:,:)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, dimension(:,:,:,:), pointer :: &
         qxm, qym, qzm, qxp, qyp, qzp, ql, qr, &
         qmxy, qpxy, qmxz, qpxz, qmyx, qpyx, &
         qmyz, qpyz, qmzx, qpzx, qmzy, qpzy, &
         qxl, qxr, qyl, qyr, qzl, qzr

    double precision, dimension(:,:,:,:), pointer:: &
         fx, fy, fz, fxy, fxz, fyx, fyz, fzx, fzy

#ifdef RADIATION
    double precision, dimension(:,:,:,:), pointer:: &
         rfx, rfy, rfz, rfxy, rfxz, rfyx, rfyz, rfzx, rfzy
#endif

    double precision, dimension(:,:,:,:), pointer:: &
         qgdnvx, qgdnvy, qgdnvz, &
         qgdnvxy, qgdnvxz, &
         qgdnvyx, qgdnvyz, &
         qgdnvzx, qgdnvzy

    ! these will be the temporary arrays we actually allocate space for
    double precision, dimension(:,:,:,:), pointer :: ftmp1, ftmp2, rftmp1, rftmp2
    double precision, dimension(:,:,:,:), pointer :: qgdnvtmp1, qgdnvtmp2

    type (eos_t) :: eos_state

    logical :: source_nonzero(QVAR)

    integer :: fglo(3), fghi(3), glo(3), ghi(3)

    ! Local constants
    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    dtdx = dt*dxinv
    dtdy = dt*dyinv
    dtdz = dt*dzinv

    hdt = HALF*dt

    hdtdx = HALF*dtdx
    hdtdy = HALF*dtdy
    hdtdz = HALF*dtdz

    cdtdx = dtdx*THIRD
    cdtdy = dtdy*THIRD
    cdtdz = dtdz*THIRD

    fglo = lo - dg(:)  ! face + one ghost
    fghi = hi + 2*dg(:)

    glo = lo - dg(:)  ! one ghost,  this can be used for face-based arrays too
    ghi = hi + dg(:)


    call bl_allocate ( qxm, fglo, fghi, NQ)
    call bl_allocate ( qxp, fglo, fghi, NQ)

    call bl_allocate ( qym, fglo, fghi, NQ)
    call bl_allocate ( qyp, fglo, fghi, NQ)

    call bl_allocate ( qzm, fglo, fghi, NQ)
    call bl_allocate ( qzp, fglo, fghi, NQ)


    if (ppm_type .gt. 0) then
       ! x-index, y-index, z-index, dim, characteristics, variables
       call bl_allocate ( Ip, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,3,1,3,1,NQ)
       call bl_allocate ( Im, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,3,1,3,1,NQ)

       ! for source terms
       call bl_allocate ( Ip_src, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,3,1,3,1,QVAR)
       call bl_allocate ( Im_src, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,3,1,3,1,QVAR)

       ! for gamc -- needed for the reference state in eigenvectors
       call bl_allocate ( Ip_gc, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,3,1,3,1,1)
       call bl_allocate ( Im_gc, glo(1),ghi(1),glo(2),ghi(2),glo(3),ghi(3),1,3,1,3,1,1)
    else
       call bl_allocate ( dqx, glo, ghi, NQ)
       call bl_allocate ( dqy, glo, ghi, NQ)
       call bl_allocate ( dqz, glo, ghi, NQ)
    end if

    ! for the hybrid Riemann solver
    call bl_allocate(shk, glo, ghi)


#ifdef SHOCK_VAR
    uout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),USHK) = ZERO

    call shock(q, qd_lo, qd_hi, &
               shk, glo, ghi, &
               lo, hi, dx)

    ! Store the shock data for future use in the burning step.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,k,USHK) = shk(i,j,k)
          enddo
       enddo
    enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:,:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, qd_lo, qd_hi, &
                  shk, glo, ghi, &
                  lo, hi, dx)
    else
       shk(:,:,:) = ZERO
    endif
#endif


    call bl_allocate(sxm, glo, ghi)
    call bl_allocate(sxp, glo, ghi)
    call bl_allocate(sym, glo, ghi)
    call bl_allocate(syp, glo, ghi)
    call bl_allocate(szm, glo, ghi)
    call bl_allocate(szp, glo, ghi)

    ! preprocess the sources -- we don't want to trace under a source that is empty
    if (ppm_type > 0) then
       do n = 1, QVAR
          if (minval(srcQ(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,n)) == ZERO .and. &
              maxval(srcQ(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,n)) == ZERO) then
             source_nonzero(n) = .false.
          else
             source_nonzero(n) = .true.
          endif
       enddo
    endif

    if (ppm_type > 0) then

       do n = 1, NQ
          call ppm_reconstruct(q, qd_lo, qd_hi, NQ, n, &
                               flatn, qd_lo, qd_hi, &
                               sxm, sxp, sym, syp, szm, szp, glo, ghi, &
                               lo, hi, dx)

          call ppm_int_profile(q, qd_lo, qd_hi, NQ, n, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, szm, szp, glo, ghi, &
                               Ip, Im, glo, ghi, NQ, n, &
                               lo, hi, dx, dt)
       end do

       ! source terms
       do n = 1, QVAR
          if (source_nonzero(n)) then
             call ppm_reconstruct(srcQ, src_lo, src_hi, QVAR, n, &
                                  flatn, qd_lo, qd_hi, &
                                  sxm, sxp, sym, syp, szm, szp, glo, ghi, &
                                  lo, hi, dx)

             call ppm_int_profile(srcQ, src_lo, src_hi, QVAR, n, &
                                  q, qd_lo, qd_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  sxm, sxp, sym, syp, szm, szp, glo, ghi, &
                                  Ip_src, Im_src, glo, ghi, QVAR, n, &
                                  lo, hi, dx, dt)
          else
             Ip_src(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),:,:,n) = ZERO
             Im_src(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),:,:,n) = ZERO
          endif

       enddo

       ! this probably doesn't support radiation
       if (ppm_temp_fix /= 1) then
          call ppm_reconstruct(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               flatn, qd_lo, qd_hi, &
                               sxm, sxp, sym, syp, szm, szp, glo, ghi, &
                               lo, hi, dx)

          call ppm_int_profile(qaux, qa_lo, qa_hi, NQAUX, QGAMC, &
                               q, qd_lo, qd_hi, &
                               qaux, qa_lo, qa_hi, &
                               sxm, sxp, sym, syp, szm, szp, glo, ghi, &
                               Ip_gc, Im_gc, glo, ghi, 1, 1, &
                               lo, hi, dx, dt)
       else

          do iwave = 1, 3
             do idim = 1, 3

                do k = lo(3)-1, hi(3)+1
                   do j = lo(2)-1, hi(2)+1
                      do i = lo(1)-1, hi(1)+1

                         eos_state % rho = Ip(i,j,k,idim,iwave,QRHO)
                         eos_state % T   = Ip(i,j,k,idim,iwave,QTEMP)

                         eos_state % xn  = Ip(i,j,k,idim,iwave,QFS:QFS+nspec-1)
                         eos_state % aux = Ip(i,j,k,idim,iwave,QFX:QFX+naux-1)

                         call eos(eos_input_rt, eos_state)

                         Ip(i,j,k,idim,iwave,QPRES)  = eos_state % p
                         Ip(i,j,k,idim,iwave,QREINT) = eos_state % e * Ip(i,j,k,idim,iwave,QRHO)
                         Ip_gc(i,j,k,idim,iwave,1)   = eos_state % gam1
                      end do
                   end do
                end do

                do k = lo(3)-1, hi(3)+1
                   do j = lo(2)-1, hi(2)+1
                      do i = lo(1)-1, hi(1)+1
                         eos_state % rho = Im(i,j,k,idim,iwave,QRHO)
                         eos_state % T   = Im(i,j,k,idim,iwave,QTEMP)

                         eos_state % xn  = Im(i,j,k,idim,iwave,QFS:QFS+nspec-1)
                         eos_state % aux = Im(i,j,k,idim,iwave,QFX:QFX+naux-1)

                         call eos(eos_input_rt, eos_state)

                         Im(i,j,k,idim,iwave,QPRES)  = eos_state % p
                         Im(i,j,k,idim,iwave,QREINT) = eos_state % e * Im(i,j,k,idim,iwave,QRHO)
                         Im_gc(i,j,k,idim,iwave,1)   = eos_state % gam1
                      end do
                   end do
                end do

             end do
          end do

       end if


       ! Compute U_x and U_y at kc (k3d)

#ifdef RADIATION
       call tracexy_ppm_rad(q, qd_lo, qd_hi, &
                            qaux, qa_lo, qa_hi, &
                            Ip, Im, Ip_src, Im_src, glo, ghi, &
                            qxm, qxp, qym, qyp, fglo, fghi, &
                            lo, hi, domlo, domhi, &
                            dx, dt)

       call tracez_ppm_rad(q, qd_lo, qd_hi, &
                           qaux, qa_lo, qa_hi, &
                           Ip, Im, Ip_src, Im_src, glo, ghi, &
                           qzm, qzp, fglo, fghi, &
                           lo, hi, domlo, domhi, &
                           dt)

#else
       call tracexy_ppm(q, qd_lo, qd_hi, &
                        qaux, qa_lo, qa_hi, &
                        Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, glo, ghi, &
                        qxm, qxp, qym, qyp, fglo, fghi, &
                        lo, hi, domlo, domhi, &
                        dx, dt)

       call tracez_ppm(q, qd_lo, qd_hi, &
                       qaux, qa_lo, qa_hi, &
                       Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc, glo, ghi, &
                       qzm, qzp, fglo, fghi, &
                       lo, hi, domlo, domhi, &
                       dt)
#endif

    else

#ifdef RADIATION
#ifndef AMREX_USE_CUDA
       call amrex_error("ppm_type <=0 is not supported in with radiation")
#endif
#endif

       ! Compute all slopes at kc (k3d)
       call uslope(q, flatn, qd_lo, qd_hi, &
                   dqx, dqy, dqz, glo, ghi, &
                   lo, hi)

       if (use_pslope .eq. 1) &
            call pslope(q, flatn, qd_lo, qd_hi, &
                        dqx, dqy, dqz, glo, ghi, &
                        srcQ, src_lo, src_hi, &
                        lo, hi, dx)

       ! Compute U_x and U_y at kc (k3d)
       call tracexy(q, qd_lo, qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    dqx, dqy, glo, ghi, &
                    qxm, qxp, qym, qyp, fglo, fghi, &
                    lo, hi, domlo, domhi, &
                    dx, dt)

       ! we should not land here with radiation
       call tracez(q, qd_lo, qd_hi, &
                   qaux, qa_lo, qa_hi, &
                   dqz, glo, ghi, &
                   qzm, qzp, fglo, fghi, &
                   lo, hi, domlo, domhi, &
                   dx, dt)

    end if  ! ppm test

    if (ppm_type .gt. 0) then
       call bl_deallocate ( Ip)
       call bl_deallocate ( Im)

       call bl_deallocate ( Ip_src)
       call bl_deallocate ( Im_src)

       call bl_deallocate ( Ip_gc)
       call bl_deallocate ( Im_gc)
    else
       call bl_deallocate ( dqx)
       call bl_deallocate ( dqy)
       call bl_deallocate ( dqz)
    end if

    ! Deallocate arrays
    call bl_deallocate(sxm)
    call bl_deallocate(sxp)
    call bl_deallocate(sym)
    call bl_deallocate(syp)
    call bl_deallocate(szm)
    call bl_deallocate(szp)


    call bl_allocate( qmxy, fglo, fghi, NQ)
    call bl_allocate( qpxy, fglo, fghi, NQ)

    call bl_allocate( qmxz, fglo, fghi, NQ)
    call bl_allocate( qpxz, fglo, fghi, NQ)

    call bl_allocate( qmyx, fglo, fghi, NQ)
    call bl_allocate( qpyx, fglo, fghi, NQ)

    call bl_allocate( qmyz, fglo, fghi, NQ)
    call bl_allocate( qpyz, fglo, fghi, NQ)

    call bl_allocate( qmzx, fglo, fghi, NQ)
    call bl_allocate( qpzx, fglo, fghi, NQ)

    call bl_allocate( qmzy, fglo, fghi, NQ)
    call bl_allocate( qpzy, fglo, fghi, NQ)

    call bl_allocate( ql, fglo, fghi, NQ)
    call bl_allocate( qr, fglo, fghi, NQ)

    call bl_allocate( ftmp1, glo, ghi, NVAR)
    call bl_allocate( ftmp2, glo, ghi, NVAR)
#ifdef RADIATION
    call bl_allocate ( rftmp1, glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3), 0, ngroups-1)
    call bl_allocate ( rftmp2, glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3), 0, ngroups-1)
#endif
    call bl_allocate ( qgdnvtmp1, fglo, fghi, NGDNV)
    call bl_allocate ( qgdnvtmp2, fglo, fghi, NGDNV)


    !-------------------------------------------------------------------------!
    ! Some notes on the work index (i.e., lo and hi arguments near the end    !
    !                               of the argument list).                    !
    ! * For cmpflx, we use face index in the flux direction and cell-centered !
    !   index for others.                                                     !
    ! * For trans*, we use cell-centered index of the valid region.           !
    !-------------------------------------------------------------------------!

    fx     =>     ftmp1
#ifdef RADIATION
    rfx    =>    rftmp1
#endif
    qgdnvx  =>  qgdnvtmp1

    ! Compute F^x
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fx, ugdnvx, pgdnvx, gegdnvx : xface, +-1 at y & z
    call cmpflx(qxm, qxp, fglo, fghi, &
                fx, glo, ghi, &
                qgdnvx, fglo, fghi, &
#ifdef RADIATION
                rfx, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                1, [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! add the transverse flux difference in x to the y and z states
    ! Inputs: qym, qyp                     : yface, +-1 at x & z
    !         qzm, qzp                     : zface, +-1 at x & y
    !         fx, ugdnvx, pgdnvx, gegdnvx  : xface, +-1 at y & z
    !         gamc                         : +-4
    ! Outputs: qmyx, qpyx                  : yface, +-0 at x, +-1 at z
    !          qmzx, qpzx                  : zface, +-0 at x, +-1 at y
    call transx(qym, qmyx, qyp, qpyx, &
                qzm, qmzx, qzp, qpzx, fglo, fghi, &
                qaux, qa_lo, qa_hi, &
                fx, &
#ifdef RADIATION
                rfx, &
#endif
                glo, ghi, &
                qgdnvx, fglo, fghi, &
                cdtdx, lo, hi)

    nullify(fx, qgdnvx)
#ifdef RADIATION
    nullify(rfx)
#endif

    fy     =>     ftmp1
#ifdef RADIATION
    rfy    =>    rftmp1
#endif
    qgdnvy => qgdnvtmp1

    ! Compute F^y
    ! Inputs: qym, qyp                     : yface, +-1 at x & z
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fy, ugdnvy, pgdnvy, gegdnvy : yface, +-1 at x & z
    call cmpflx(qym, qyp, fglo, fghi, &
                fy, glo, ghi, &
                qgdnvy, fglo, fghi, &
#ifdef RADIATION
                rfy, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                2, [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! add the transverse flux difference in y to the x and z states
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         qzm, qzp                     : zface, +-1 at x & y
    !         fy, ugdnvy, pgdnvy, gegdnvy  : yface, +-1 at x & z
    !         gamc                         : +-4
    ! Outputs: qmxy, qpxy                  : xface, +-0 at y, +-1 at z
    !          qmzy, qpzy                  : zface, +-0 at y, +-1 at x
    call transy(qxm, qmxy, qxp, qpxy, &
                qzm, qmzy, qzp, qpzy, fglo, fghi, &
                qaux, qa_lo, qa_hi, &
                fy, &
#ifdef RADIATION
                rfy, &
#endif
                glo, ghi, &
                qgdnvy, fglo, fghi, &
                cdtdy, lo, hi)

    nullify(fy, qgdnvy)
#ifdef RADIATION
    nullify(rfy)
#endif

    fz      =>     ftmp1
#ifdef RADIATION
    rfz     =>    rftmp1
#endif
    qgdnvz  =>  qgdnvtmp1

    ! Compute F^z
    ! Inputs: qzm, qzp                     : zface, +-1 at x & y
    !         gamc, csml, c                : +-4
    !         shk                          : +-1
    ! Outputs: fz, ugdnvz, pgdnvz, gegdnvz : zface, +-1 at x & y
    call cmpflx(qzm, qzp, fglo, fghi, &
                fz, glo, ghi, &
                qgdnvz, fglo, fghi, &
#ifdef RADIATION
                rfz, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                3, [lo(1)-1, lo(2)-1, lo(3)], [ hi(1)+1, hi(2)+1, hi(3)+1], &
                domlo, domhi)

    ! add the transverse flux difference in z to the x and y states
    ! Inputs: qxm, qxp                     : xface, +-1 at y & z
    !         qym, qyp                     : yface, +-1 at x & z
    !         fz, ugdnvz, pgdnvz, gegdnvz  : zface, +-1 at x & y
    !         gamc                         : +-4
    ! Outputs: qmxz, qpxz                  : xface, +-0 at z, +-1 at y
    !          qmyz, qpyz                  : yface, +-0 at z, +-1 at x
    call transz(qxm, qmxz, qxp, qpxz, &
                qym, qmyz, qyp, qpyz, fglo, fghi, &
                qaux, qa_lo, qa_hi, &
                fz, &
#ifdef RADIATION
                rfz, &
#endif
                glo, ghi, &
                qgdnvz, fglo, fghi, &
                cdtdz, lo, hi)

    nullify(fz, qgdnvz)
#ifdef RADIATION
    nullify(rfz)
#endif

    ! We now have qx?, qy?, qz?
    !         and q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

    !
    ! Use qx?, q?yz, q?zy to compute final x-flux
    !

    fyz      =>      ftmp1
#ifdef RADIATION
    rfyz     =>     rftmp1
#endif
    qgdnvyz  =>  qgdnvtmp1

    ! Compute F^{y|z}
    ! Inputs: qmyz, qpyz                       : yface, +-1 at x, +-0 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fyz, ugdnvyz, pgdnvyz, gegdnvyz : yface, +-1 at x, +-0 at z
    call cmpflx(qmyz, qpyz, fglo, fghi, &
                fyz, glo, ghi, &
                qgdnvyz, fglo, fghi, &
#ifdef RADIATION
                rfyz, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                2, [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)], &
                domlo, domhi)

    fzy      =>      ftmp2
#ifdef RADIATION
    rfzy     =>     rftmp2
#endif
    qgdnvzy  =>  qgdnvtmp2

    ! Compute F^{z|y}
    ! Inputs: qmzy, qpzy                       : zface, +-1 at x, +-0 at y
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fzy, ugdnvzy, pgdnvzy, gegdnvzy : zface, +-1 at x, +-0 at y
    call cmpflx(qmzy, qpzy, fglo, fghi, &
                fzy, glo, ghi, &
                qgdnvzy, fglo, fghi, &
#ifdef RADIATION
                rfzy, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                3, [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1], &
                domlo, domhi)

    qxl => ql
    qxr => qr

    ! Compute the corrected x interface states
    ! Inputs: qxm, qxp                        : xface, +-1 at y & z
    !         fyz, ugdnvyz, pgdnvyz, gegdnvyz : yface, +-1 at x, +-0 at z
    !         fzy, ugdnvzy, pgdnvzy, gegdnvzy : zface, +-1 at x, +-0 at y
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qxl, qxr                       : xface, +-0 at y & z
    call transyz(qxm, qxl, qxp, qxr, fglo, fghi, &
                 qaux, qa_lo, qa_hi, &
                 fyz, &
#ifdef RADIATION
                 rfyz, &
#endif
                 glo, ghi, &
                 fzy, &
#ifdef RADIATION
                 rfzy, &
#endif
                 glo, ghi, &
                 qgdnvyz, fglo, fghi, &
                 qgdnvzy, fglo, fghi, &
                 srcQ, src_lo, src_hi, &
                 hdt, hdtdy, hdtdz, lo, hi)

    nullify(fyz, qgdnvyz)
    nullify(fzy, qgdnvzy)
#ifdef RADIATION
    nullify(rfyz, rfzy)
#endif

    qgdnvx  =>  qgdnvtmp1

    ! now compute the final x fluxes, F^x
    call cmpflx(qxl, qxr, fglo, fghi, &
                flux1, f1_lo, f1_hi, &
                qgdnvx, fglo, fghi, &
#ifdef RADIATION
                rflux1, rf1_lo, rf1_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                1, lo, [hi(1)+1, hi(2), hi(3)], domlo, domhi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i=lo(1), hi(1)+1
             q1(i,j,k,:) = qgdnvx(i,j,k,:)
          end do
       end do
    end do

    nullify(qgdnvx)
    nullify(qxl, qxr)

    !
    ! Use qy?, q?zx, q?xz to compute final y-flux
    !

    fzx      =>      ftmp1
#ifdef RADIATION
    rfzx     =>     rftmp1
#endif
    qgdnvzx  =>  qgdnvtmp1

    ! Compute F^{z|x}
    ! Inputs: qmzx, qpzx                       : zface, +-0 at x, +-1 at y
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fzx, ugdnvzx, pgdnvzx, gegdnvzx : zface, +-0 at x, +-1 at y
    call cmpflx(qmzx, qpzx, fglo, fghi, &
                fzx, glo, ghi, &
                qgdnvzx, fglo, fghi, &
#ifdef RADIATION
                rfzx, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                3, [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1], &
                domlo, domhi)

    fxz      =>      ftmp2
#ifdef RADIATION
    rfxz     =>     rftmp2
#endif
    qgdnvxz  =>  qgdnvtmp2

    ! Compute F^{x|z} at km
    ! Inputs: qmxz, qpxz                       : xface, +-1 at y, +-0 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fxz, ugdnvxz, pgdnvxz, gegdnvxz : xface, +-1 at y, +-0 at z
    call cmpflx(qmxz, qpxz, fglo, fghi, &
                fxz, glo, ghi, &
                qgdnvxz, fglo, fghi, &
#ifdef RADIATION
                rfxz, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                1, [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)], &
                domlo, domhi)

    qyl => ql
    qyr => qr

    ! Compute the corrected y interface states
    ! Inputs: qym, qyp                        : yface, +-1 at x & z
    !         fxz, ugdnvxz, pgdnvxz, gegdnvxz : xface, +-1 at y, +-0 at z
    !         fzx, ugdnvzx, pgdnvzx, gegdnvzx : zface, +-0 at x, +-1 at y
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qyl, qyr                       : yface, +-0 at x & z
    call transxz(qym, qyl, qyp, qyr, fglo, fghi, &
                 qaux, qa_lo, qa_hi, &
                 fxz, &
#ifdef RADIATION
                 rfxz, &
#endif
                 glo, ghi, &
                 fzx, &
#ifdef RADIATION
                 rfzx, &
#endif
                 glo, ghi, &
                 qgdnvxz, fglo, fghi, &
                 qgdnvzx, fglo, fghi, &
                 srcQ, src_lo, src_hi, &
                 hdt, hdtdx, hdtdz, lo, hi)

    nullify(fzx, qgdnvzx)
    nullify(fxz, qgdnvxz)
#ifdef RADIATION
    nullify(rfzx, rfxz)
#endif

    qgdnvy  =>  qgdnvtmp1

    ! now compute the final y fluxes F^y
    ! Inputs: qyl, qyr                        : yface, +-0 at x & y
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux2, ugdnvy, pgdnvy, gegdnvy : yface, +-0 at x & y
    call cmpflx(qyl, qyr, fglo, fghi, &
                flux2, f2_lo, f2_hi, &
                qgdnvy, fglo, fghi, &
#ifdef RADIATION
                rflux2, rf2_lo, rf2_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                2, [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], domlo, domhi)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             q2(i,j,k,:) = qgdnvy(i,j,k,:)
          end do
       end do
    end do

    nullify(qgdnvy)
    nullify(qyl,qyr)

    !
    ! Use qz?, q?xy, q?yx to compute final z-flux
    !

    fxy      =>      ftmp1
#ifdef RADIATION
    rfxy     =>     rftmp1
#endif
    qgdnvxy  =>   qgdnvtmp1


    ! Compute F^{x|y}
    ! Inputs: qmxy, qpxy                       : xface, +-0 at y, +-1 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fxy, ugdnvxy, pgdnvxy, gegdnvxy : xface, +-0 at y, +-1 at z
    call cmpflx(qmxy, qpxy, fglo, fghi, &
                fxy, glo, ghi, &
                qgdnvxy, fglo, fghi, &
#ifdef RADIATION
                rfxy, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                1, [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1], &
                domlo, domhi)

    fyx      =>      ftmp2
#ifdef RADIATION
    rfyx     =>     rftmp2
#endif
    qgdnvyx  =>  qgdnvtmp2

    ! Compute F^{y|x}
    ! Inputs: qmyx, qpyx                       : yface, +-0 at x, +-1 at z
    !         gamc, csml, c                    : +-4
    !         shk                              : +-1
    ! Outputs: fyx, ugdnvyx, pgdnvyx, gegdnvyx : yface, +-0 at x, +-1 at z
    call cmpflx(qmyx, qpyx, fglo, fghi, &
                fyx, glo, ghi, &
                qgdnvyx, fglo, fghi, &
#ifdef RADIATION
                rfyx, glo, ghi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                2, [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1], &
                domlo, domhi)

    qzl => ql
    qzr => qr

    ! Compute the corrected z interface states
    ! Inputs: qzm, qzp                        : zface, +-1 at x & y
    !         fxy, ugdnvxy, pgdnvxy, gegdnvxy : xface, +-0 at y, +-1 at z
    !         fyx, ugdnvyx, pgdnvyx, gegdnvyx : yface, +-0 at x, +-1 at z
    !         gamc, grav, rot                 : +-4
    !         srcQ                            : +-1
    ! Outputs: qzl, qzr                       : zface, +-0 at x & y
    call transxy(qzm, qzl, qzp, qzr, fglo, fghi, &
                 qaux, qa_lo, qa_hi, &
                 fxy, &
#ifdef RADIATION
                 rfxy, &
#endif
                 glo, ghi, &
                 fyx,&
#ifdef RADIATION
                 rfyx, &
#endif
                 glo, ghi, &
                 qgdnvxy, fglo, fghi, &
                 qgdnvyx, fglo, fghi, &
                 srcQ, src_lo, src_hi,&
                 hdt, hdtdx, hdtdy, lo, hi)

    nullify(fxy, qgdnvxy)
    nullify(fyx, qgdnvyx)
#ifdef RADIATION
    nullify(rfxy, rfyz)
#endif

    qgdnvz  =>  qgdnvtmp1

    ! Compute the final z fluxes F^z
    ! Inputs: qzl, qzr                        : zface, +-0 at x & y
    !         gamc, csml, c                   : +-4
    !         shk                             : +-1
    ! Outputs: flux3, ugdnvz, pgdnvz, gegdnvz : zface, +-0 at x & y
    call cmpflx(qzl, qzr, fglo, fghi, &
                flux3, f3_lo, f3_hi, &
                qgdnvz, fglo, fghi, &
#ifdef RADIATION
                rflux3, rf3_lo, rf3_hi, &
#endif
                qaux, qa_lo, qa_hi, &
                shk, glo, ghi, &
                3, [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1], &
                domlo, domhi)

    do k = lo(3), hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             q3(i,j,k,:) = qgdnvz(i,j,k,:)
          end do
       end do
    end do

    nullify(qgdnvz)
    nullify(qzl,qzr)

    call bl_deallocate(shk)

    call bl_deallocate ( qxm)
    call bl_deallocate ( qxp)

    call bl_deallocate ( qym)
    call bl_deallocate ( qyp)

    call bl_deallocate ( qzm)
    call bl_deallocate ( qzp)

    call bl_deallocate ( ql)
    call bl_deallocate ( qr)

    call bl_deallocate ( qmxy)
    call bl_deallocate ( qpxy)

    call bl_deallocate ( qmxz)
    call bl_deallocate ( qpxz)

    call bl_deallocate ( qmzx)
    call bl_deallocate ( qpzx)

    call bl_deallocate ( qmzy)
    call bl_deallocate ( qpzy)

    call bl_deallocate ( qmyx)
    call bl_deallocate ( qpyx)

    call bl_deallocate ( qmyz)
    call bl_deallocate ( qpyz)

    call bl_deallocate ( ftmp1)
    call bl_deallocate ( ftmp2)

#ifdef RADIATION
    call bl_deallocate (rftmp1)
    call bl_deallocate (rftmp2)
#endif

    call bl_deallocate(qgdnvtmp1)
    call bl_deallocate(qgdnvtmp2)

  end subroutine umeth

end module ctu_advection_module
