module transverse_module

  use bl_constants_module

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, GDERADS, GDLAMS, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_trace_sources, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif

  use eos_module

  implicit none

  private
  public :: transx1, transx2, transy1, transy2, transz, transxy, transyz, transxz

contains

  subroutine reset_edge_state_thermo(qedge, qd_lo, qd_hi, ii, jj, kk)

    integer, intent(in) :: ii, jj, kk
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    double precision, intent(inout) :: qedge(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)

    logical :: reset
    type (eos_t) :: eos_state

    reset = .false.

    if (transverse_reset_rhoe == 1) then
       ! if we are still negative, then we need to reset
       if (qedge(ii,jj,kk,QREINT) < ZERO) then
          reset = .true.

          eos_state % rho = qedge(ii,jj,kk,QRHO)
          eos_state % T = small_temp
          eos_state % xn(:) = qedge(ii,jj,kk,QFS:QFS-1+nspec)
          eos_state % aux(:) = qedge(ii,jj,kk,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          qedge(ii,jj,kk,QREINT) = qedge(ii,jj,kk,QRHO)*eos_state % e
          qedge(ii,jj,kk,QPRES) = eos_state % p
       endif

    end if

    if (ppm_predict_gammae == 0 ) then

       if (transverse_use_eos == 1) then
          eos_state % rho = qedge(ii,jj,kk,QRHO)
          eos_state % e   = qedge(ii,jj,kk,QREINT) / qedge(ii,jj,kk,QRHO)
          eos_state % T   = small_temp
          eos_state % xn  = qedge(ii,jj,kk,QFS:QFS+nspec-1)
          eos_state % aux = qedge(ii,jj,kk,QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          qedge(ii,jj,kk,QREINT) = eos_state % e * eos_state % rho
          qedge(ii,jj,kk,QPRES) = max(eos_state % p, small_pres)
       end if

    else
       if (reset) then
          ! recompute the p edge state from this and (rho e), since we reset
          ! qreint  (actually, is this code even necessary?)
          qedge(ii,jj,kk,QPRES) = qedge(ii,jj,kk,QREINT)*(qedge(ii,jj,kk,QGAME)-ONE)
          qedge(ii,jj,kk,QPRES) = max(qedge(ii,jj,kk,QPRES), small_pres)
       end if
    end if

  end subroutine reset_edge_state_thermo


  !===========================================================================
  ! transx1
  !===========================================================================
  subroutine transx1( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qym,qymo,qyp,qypo,qd_lo,qd_hi, &
                     fx, &
#ifdef RADIATION
                     rfx, &
#endif
                     fx_lo,fx_hi, &
                     qx,qx_lo,qx_hi, &
                     gamc,gd_lo,gd_hi, &
                     cdtdx,ilo,ihi,jlo,jhi,kc,k3d)

    ! Note that what we call ilo here is ilo = lo(1)
    ! Note that what we call ihi here is ihi = hi(1)
    ! Note that what we call jlo here is jlo = lo(2) - 1
    ! Note that what we call jhi here is jhi = hi(2) + 1

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer :: gd_lo(3),gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
#endif

    double precision  qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qymo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qypo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdx

    integer i, j, n, nq, ipassive

    double precision rhoinv
    double precision rrnew, rr
    double precision rrry, rrly
    double precision rury, ruly
    double precision rvry, rvly
    double precision rwry, rwly
    double precision ekenry, ekenly
    double precision rery, rely
    double precision rrnewry, rrnewly
    double precision runewry, runewly
    double precision rvnewry, rvnewly
    double precision rwnewry, rwnewly
    double precision renewry, renewly
    double precision pnewry, pnewly
    double precision rhoekenry, rhoekenly
    double precision compn, compu
    double precision pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

#ifdef RADIATION
    double precision :: dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state


    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transverse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             if (j >= jlo+1) then
                rr = qyp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qyp(i,j,kc,nq) - compn
                qypo(i,j,kc,nq) = compu/rrnew
             end if

             if (j <= jhi-1) then
                rr = qym(i,j+1,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qym(i,j+1,kc,nq) - compn
                qymo(i,j+1,kc,nq) = compu/rrnew
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to y-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          pggp  = qx(i+1,j,kc,GDPRES)
          pggm  = qx(i  ,j,kc,GDPRES)
          ugp  = qx(i+1,j,kc,GDU   )
          ugm  = qx(i  ,j,kc,GDU   )
          gegp = qx(i+1,j,kc,GDGAME)
          gegm = qx(i  ,j,kc,GDGAME)

#ifdef RADIATION
          lambda = lam(i,j,k3d,:)
          ugc = HALF*(ugp+ugm)
          ergp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          luge = ugc * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          if (j >= jlo+1) then

             ! Convert to conservation form
             rrry = qyp(i,j,kc,QRHO)
             rury = rrry*qyp(i,j,kc,QU)
             rvry = rrry*qyp(i,j,kc,QV)
             rwry = rrry*qyp(i,j,kc,QW)
             ekenry = HALF*rrry* &
                  (qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
             rery = qyp(i,j,kc,QREINT) + ekenry
#ifdef RADIATION
             err  = qyp(i,j,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewry = rury - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
             rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
#ifdef RADIATION
             runewry = runewry + dmom
             renewry = renewry + dre
             ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewry < ZERO) then
                rrnewry = rrry
                runewry = rury
                rvnewry = rvry
                rwnewry = rwry
                renewry = rery
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             qypo(i,j,kc,QRHO) = rrnewry
             rhoinv = ONE/rrnewry
             qypo(i,j,kc,QU) = runewry*rhoinv
             qypo(i,j,kc,QV) = rvnewry*rhoinv
             qypo(i,j,kc,QW) = rwnewry*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
             qypo(i,j,kc,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe .eq. 1 .and. qypo(i,j,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qypo(i,j,kc,QREINT) = qyp(i,j,kc,QREINT) - &
                        cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-ONE))
                   qypo(i,j,kc,QPRES) = max(pnewry,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qypo(i,j,kc,QGAME) = qyp(i,j,kc,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qypo(i,j,kc,QPRES) = qypo(i,j,kc,QREINT)*(qypo(i,j,kc,QGAME)-ONE)
                   qypo(i,j,kc,QPRES) = max(qypo(i,j,kc,QPRES),small_pres)
                end if
             else
                qypo(i,j,kc,QPRES) = qyp(i,j,kc,QPRES)
                qypo(i,j,kc,QGAME) = qyp(i,j,kc,QGAME)
             endif

             call reset_edge_state_thermo(qypo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
             qypo(i,j,kc,qrad:qradhi) = ernewr(:)
             qypo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,kc,QPRES)
             qypo(i,j,kc,qreitot) = sum(qypo(i,j,kc,qrad:qradhi)) + qypo(i,j,kc,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------

          if (j <= jhi-1) then

             ! Convert to conservation form
             rrly = qym(i,j+1,kc,QRHO)
             ruly = rrly*qym(i,j+1,kc,QU)
             rvly = rrly*qym(i,j+1,kc,QV)
             rwly = rrly*qym(i,j+1,kc,QW)
             ekenly = HALF*rrly* &
                  (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
             rely = qym(i,j+1,kc,QREINT) + ekenly
#ifdef RADIATION
             erl  = qym(i,j+1,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
             rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
#ifdef RADIATION
             runewly = runewly + dmom
             renewly = renewly + dre
             ernewl  = erl(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewly < ZERO) then
                rrnewly = rrly
                runewly = ruly
                rvnewly = rvly
                rwnewly = rwly
                renewly = rely
#ifdef RADIATION
                ernewl  = erl(:)
#endif
                reset_state = .true.
             endif

             qymo(i,j+1,kc,QRHO) = rrnewly
             rhoinv = ONE/rrnewly
             qymo(i,j+1,kc,QU) = runewly*rhoinv
             qymo(i,j+1,kc,QV) = rvnewly*rhoinv
             qymo(i,j+1,kc,QW) = rwnewly*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
             qymo(i,j+1,kc,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary
                if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qymo(i,j+1,kc,QREINT) = qym(i,j+1,kc,QREINT) - &
                        cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-ONE))
                   qymo(i,j+1,kc,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qymo(i,j+1,kc,QGAME) = qym(i,j+1,kc,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qymo(i,j+1,kc,QPRES) = qymo(i,j+1,kc,QREINT)*(qymo(i,j+1,kc,QGAME)-ONE)
                   qymo(i,j+1,kc,QPRES) = max(qymo(i,j+1,kc,QPRES), small_pres)
                end if
             else
                qymo(i,j+1,kc,QPRES) = qym(i,j+1,kc,QPRES)
                qymo(i,j+1,kc,QGAME) = qym(i,j+1,kc,QGAME)
             endif

             call reset_edge_state_thermo(qymo, qd_lo, qd_hi, i, j+1, kc)

#ifdef RADIATION
             qymo(i,j+1,kc,qrad:qradhi) = ernewl(:)
             qymo(i,j+1,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j+1,kc,QPRES)
             qymo(i,j+1,kc,qreitot) = sum(qymo(i,j+1,kc,qrad:qradhi)) + qymo(i,j+1,kc,QREINT)
#endif
          endif
       enddo

    enddo

  end subroutine transx1


  !===========================================================================
  ! transx2
  !===========================================================================
  subroutine transx2( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qzm,qzmo,qzp,qzpo,qd_lo,qd_hi, &
                     fx, &
#ifdef RADIATION
                     rfx, &
#endif
                     fx_lo,fx_hi, &
                     qx,qx_lo,qx_hi, &
                     gamc,gd_lo,gd_hi, &
                     cdtdx,ilo,ihi,jlo,jhi,kc,km,k3d)

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer :: gd_lo(3),gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,km,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
#endif

    double precision  qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qzmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qzpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdx

    integer i, j, n, nq, ipassive

    double precision rhoinv
    double precision rrnew, rr
    double precision rrrz, rrlz
    double precision rurz, rulz
    double precision rvrz, rvlz
    double precision rwrz, rwlz
    double precision ekenrz, ekenlz
    double precision rerz, relz
    double precision rrnewrz, rrnewlz
    double precision runewrz, runewlz
    double precision rvnewrz, rvnewlz
    double precision rwnewrz, rwnewlz
    double precision renewrz, renewlz
    double precision pnewrz, pnewlz
    double precision rhoekenrz, rhoekenlz
    double precision compn, compu
    double precision pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

#ifdef RADIATION
    double precision :: dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nq) - compn
             qzpo(i,j,kc,nq) = compu/rrnew

             compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nq) - compn
             qzmo(i,j,kc,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to z-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------

          pggp  = qx(i+1,j,kc,GDPRES)
          pggm  = qx(i  ,j,kc,GDPRES)
          ugp  = qx(i+1,j,kc,GDU   )
          ugm  = qx(i  ,j,kc,GDU   )
          gegp = qx(i+1,j,kc,GDGAME)
          gegm = qx(i  ,j,kc,GDGAME)
#ifdef RADIATION
          lambda = lam(i,j,k3d,:)
          ugc = HALF*(ugp+ugm)
          ergp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if
#endif


          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
#ifdef RADIATION
          err  = qzp(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
#ifdef RADIATION
          runewrz = runewrz + dmom
          renewrz = renewrz + dre
          ernewr  = err(:) - cdtdx*(rfx(i+1,j,kc,:) - rfx(i,j,kc,:)) &
               + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewrz < ZERO) then
             rrnewrz = rrrz
             runewrz = rurz
             rvnewrz = rvrz
             rwnewrz = rwrz
             renewrz = rerz
#ifdef RADIATION
             ernewr  = err(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzpo(i,j,kc,QRHO) = rrnewrz
          rhoinv = ONE/rrnewrz
          qzpo(i,j,kc,QU) = runewrz*rhoinv
          qzpo(i,j,kc,QV) = rvnewrz*rhoinv
          qzpo(i,j,kc,QW) = rwnewrz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz

          if (.not. reset_state) then
             ! do the transverse terms for p, gamma, and rhoe, as necessary

             if (transverse_reset_rhoe == 1 .and. qzpo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzpo(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) - &
                     cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
             end if

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-ONE))
                qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT)*(qzpo(i,j,kc,QGAME)-ONE)
                qzpo(i,j,kc,QPRES) = max(qzpo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzpo(i,j,kc,QPRES) = qzp(i,j,kc,QPRES)
             qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)
#endif

          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------

          pggp  = qx(i+1,j,km,GDPRES)
          pggm  = qx(i  ,j,km,GDPRES)
          ugp  = qx(i+1,j,km,GDU   )
          ugm  = qx(i  ,j,km,GDU   )
          gegp = qx(i+1,j,km,GDGAME)
          gegm = qx(i  ,j,km,GDGAME)
#ifdef RADIATION
          lambda = lam(i,j,k3d-1,:)
          ugc = HALF*(ugp+ugm)
          ergp = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdx*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdx * luge
          end if
#endif

          ! Convert to conservation form
          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
#ifdef RADIATION
          erl  = qzm(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
          runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX) - fx(i,j,km,UMX))
          rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY) - fx(i,j,km,UMY))
          rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ) - fx(i,j,km,UMZ))
          renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN))
#ifdef RADIATION
          runewlz = runewlz + dmom
          renewlz = renewlz + dre
          ernewl  = erl(:) - cdtdx*(rfx(i+1,j,km,:) - rfx(i,j,km,:)) &
               +der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewlz < ZERO) then
             rrnewlz = rrlz
             runewlz = rulz
             rvnewlz = rvlz
             rwnewlz = rwlz
             renewlz = relz
#ifdef RADIATION
             ernewl  = erl(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzmo(i,j,kc,QRHO) = rrnewlz
          rhoinv = ONE/rrnewlz
          qzmo(i,j,kc,QU) = runewlz*rhoinv
          qzmo(i,j,kc,QV) = rvnewlz*rhoinv
          qzmo(i,j,kc,QW) = rwnewlz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz

          if (.not. reset_state) then
             ! do the transverse terms for p, gamma, and rhoe, as necessary

             if (transverse_reset_rhoe == 1 .and. qzmo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzmo(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) - &
                     cdtdx*(fx(i+1,j,km,UEINT) - fx(i,j,km,UEINT) + pav*du)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d-1)-ONE))
                qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT)*(qzmo(i,j,kc,QGAME)-ONE)
                qzmo(i,j,kc,QPRES) = max(qzmo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzmo(i,j,kc,QPRES) = qzm(i,j,kc,QPRES)
             qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzmo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)
#endif

       enddo
    enddo

  end subroutine transx2


  !===========================================================================
  ! transy1
  !===========================================================================
  subroutine transy1( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qxm,qxmo,qxp,qxpo,qd_lo,qd_hi, &
                     fy, &
#ifdef RADIATION
                     rfy, &
#endif
                     fy_lo,fy_hi, &
                     qy,qy_lo,qy_hi, &
                     gamc,gd_lo,gd_hi, &
                     cdtdy,ilo,ihi,jlo,jhi,kc,k3d)

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer :: gd_lo(3),gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    double precision  qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qxmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qxpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdy

    integer i, j, n, nq, ipassive

    double precision rhoinv
    double precision rrnew, rr
    double precision compn, compu
    double precision rrrx, rrlx
    double precision rurx, rulx
    double precision rvrx, rvlx
    double precision rwrx, rwlx
    double precision ekenrx, ekenlx
    double precision rerx, relx
    double precision rrnewrx, rrnewlx
    double precision runewrx, runewlx
    double precision rvnewrx, rvnewlx
    double precision rwnewrx, rwnewlx
    double precision renewrx, renewlx
    double precision pnewrx, pnewlx
    double precision rhoekenrx, rhoekenlx
    double precision pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

#ifdef RADIATION
    double precision :: dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             if (i >= ilo+1) then
                rr = qxp(i,j,kc,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                compu = rr*qxp(i,j,kc,nq) - compn
                qxpo(i,j,kc,nq) = compu/rrnew
             end if

             if (i <= ihi-1) then
                rr = qxm(i+1,j,kc,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                compu = rr*qxm(i+1,j,kc,nq) - compn
                qxmo(i+1,j,kc,nq) = compu/rrnew
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to x-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          pggp  = qy(i,j+1,kc,GDPRES)
          pggm  = qy(i,j  ,kc,GDPRES)
          ugp  = qy(i,j+1,kc,GDV   )
          ugm  = qy(i,j  ,kc,GDV   )
          gegp = qy(i,j+1,kc,GDGAME)
          gegm = qy(i,j  ,kc,GDGAME)
#ifdef RADIATION
          lambda = lam(i,j,k3d,:)
          ugc = HALF*(ugp+ugm)
          ergp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          luge = ugc * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------

          if (i >= ilo+1) then
             ! Convert to conservation form
             rrrx = qxp(i,j,kc,QRHO)
             rurx = rrrx*qxp(i,j,kc,QU)
             rvrx = rrrx*qxp(i,j,kc,QV)
             rwrx = rrrx*qxp(i,j,kc,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
                  + qxp(i,j,kc,QW)**2)
             rerx = qxp(i,j,kc,QREINT) + ekenrx
#ifdef RADIATION
             err  = qxp(i,j,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
             rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
             rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
             renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))
#ifdef RADIATION
             rvnewrx = rvnewrx + dmom
             renewrx = renewrx + dre
             ernewr = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewrx < ZERO) then
                rrnewrx = rrrx
                runewrx = rurx
                rvnewrx = rvrx
                rwnewrx = rwrx
                renewrx = rerx
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             qxpo(i,j,kc,QRHO) = rrnewrx
             rhoinv = ONE/rrnewrx
             qxpo(i,j,kc,QU) = runewrx*rhoinv
             qxpo(i,j,kc,QV) = rvnewrx*rhoinv
             qxpo(i,j,kc,QW) = rwnewrx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
             qxpo(i,j,kc,QREINT) = renewrx - rhoekenrx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxpo(i,j,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qxpo(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) - &
                        cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - ONE))
                   qxpo(i,j,kc,QPRES) = max(pnewrx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxpo(i,j,kc,QGAME) = qxp(i,j,kc,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxpo(i,j,kc,QPRES) = qxpo(i,j,kc,QREINT)*(qxpo(i,j,kc,QGAME)-ONE)
                   qxpo(i,j,kc,QPRES) = max(qxpo(i,j,kc,QPRES), small_pres)
                endif
             else
                qxpo(i,j,kc,QPRES) = qxp(i,j,kc,QPRES)
                qxpo(i,j,kc,QGAME) = qxp(i,j,kc,QGAME)
             endif

             call reset_edge_state_thermo(qxpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
             qxpo(i,j,kc,qrad:qradhi) = ernewr(:)
             qxpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,kc,QPRES)
             qxpo(i,j,kc,qreitot) = sum(qxpo(i,j,kc,qrad:qradhi)) + qxpo(i,j,kc,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i <= ihi-1) then
             ! Convert to conservation form
             rrlx = qxm(i+1,j,kc,QRHO)
             rulx = rrlx*qxm(i+1,j,kc,QU)
             rvlx = rrlx*qxm(i+1,j,kc,QV)
             rwlx = rrlx*qxm(i+1,j,kc,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
                  + qxm(i+1,j,kc,QW)**2)
             relx = qxm(i+1,j,kc,QREINT) + ekenlx
#ifdef RADIATION
             erl  = qxm(i+1,j,kc,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
             rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
             rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
             renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN)- fy(i,j,kc,UEDEN))
#ifdef RADIATION
             rvnewlx = rvnewlx + dmom
             renewlx = renewlx + dre
             ernewl  = erl(:) + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewlx < ZERO) then
                rrnewlx = rrlx
                runewlx = rulx
                rvnewlx = rvlx
                rwnewlx = rwlx
                renewlx = relx
#ifdef RADIATION
                ernewl  = erl(:)
#endif
                reset_state = .true.
             endif

             qxmo(i+1,j,kc,QRHO) = rrnewlx
             rhoinv = ONE/rrnewlx
             qxmo(i+1,j,kc,QU) = runewlx*rhoinv
             qxmo(i+1,j,kc,QV) = rvnewlx*rhoinv
             qxmo(i+1,j,kc,QW) = rwnewlx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
             qxmo(i+1,j,kc,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,kc,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) - &
                        cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - ONE))
                   qxmo(i+1,j,kc,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo(i+1,j,kc,QGAME) = qxm(i+1,j,kc,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo(i+1,j,kc,QPRES) = qxmo(i+1,j,kc,QREINT)*(qxmo(i+1,j,kc,QGAME)-ONE)
                   qxmo(i+1,j,kc,QPRES) = max(qxmo(i+1,j,kc,QPRES), small_pres)
                endif
             else
                qxmo(i+1,j,kc,QPRES) = qxm(i+1,j,kc,QPRES)
                qxmo(i+1,j,kc,QGAME) = qxm(i+1,j,kc,QGAME)
             endif

             call reset_edge_state_thermo(qxmo, qd_lo, qd_hi, i+1, j, kc)

#ifdef RADIATION
             qxmo(i+1,j,kc,qrad:qradhi) = ernewl(:)
             qxmo(i+1,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i+1,j,kc,QPRES)
             qxmo(i+1,j,kc,qreitot) = sum(qxmo(i+1,j,kc,qrad:qradhi)) + qxmo(i+1,j,kc,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transy1


  !===========================================================================
  ! transy2
  !===========================================================================
  subroutine transy2( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qzm,qzmo,qzp,qzpo,qd_lo,qd_hi, &
                     fy, &
#ifdef RADIATION
                     rfy, &
#endif
                     fy_lo,fy_hi, &
                     qy,qy_lo,qy_hi, &
                     gamc,gd_lo,gd_hi, &
                     cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer :: gd_lo(3),gd_hi(3)
    integer ilo,ihi,jlo,jhi,kc,km,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    double precision  qzm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qzp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qzmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qzpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdy

    integer i, j, n, nq, ipassive

    double precision rhoinv
    double precision rrnew, rr
    double precision compn, compu
    double precision rrrz, rrlz
    double precision rurz, rulz
    double precision rvrz, rvlz
    double precision rwrz, rwlz
    double precision ekenrz, ekenlz
    double precision rerz, relz
    double precision rrnewrz, rrnewlz
    double precision runewrz, runewlz
    double precision rvnewrz, rvnewlz
    double precision rwnewrz, rwnewlz
    double precision renewrz, renewlz
    double precision pnewrz, pnewlz
    double precision rhoekenrz, rhoekenlz
    double precision pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

#ifdef RADIATION
    double precision :: dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    double precision eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))

             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nq) - compn
             qzpo(i,j,kc,nq) = compu/rrnew

             compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))

             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nq) - compn
             qzmo(i,j,kc,nq) = compu/rrnew

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to z-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          !-------------------------------------------------------------------
          ! qzpo states
          !-------------------------------------------------------------------

          pggp  = qy(i,j+1,kc,GDPRES)
          pggm  = qy(i,j  ,kc,GDPRES)
          ugp  = qy(i,j+1,kc,GDV   )
          ugm  = qy(i,j  ,kc,GDV   )
          gegp = qy(i,j+1,kc,GDGAME)
          gegm = qy(i,j  ,kc,GDGAME)
#ifdef RADIATION
          lambda = lam(i,j,k3d,:)
          ugc = HALF*(ugp+ugm)
          ergp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
               + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
#ifdef RADIATION
          err  = qzp(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))
#ifdef RADIATION
          rvnewrz = rvnewrz + dmom
          renewrz = renewrz + dre
          ernewr  = err(:) - cdtdy*(rfy(i,j+1,kc,:) - rfy(i,j,kc,:)) &
               + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewrz < ZERO) then
             rrnewrz = rrrz
             runewrz = rurz
             rvnewrz = rvrz
             rwnewrz = rwrz
             renewrz = rerz
#ifdef RADIATION
             ernewr  = err(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzpo(i,j,kc,QRHO) = rrnewrz
          rhoinv = ONE/rrnewrz
          qzpo(i,j,kc,QU) = runewrz*rhoinv
          qzpo(i,j,kc,QV) = rvnewrz*rhoinv
          qzpo(i,j,kc,QW) = rwnewrz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qzpo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzpo(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) - &
                     cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - ONE))
                qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT)*(qzpo(i,j,kc,QGAME)-ONE)
                qzpo(i,j,kc,QPRES) = max(qzpo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzpo(i,j,kc,QPRES) = qzp(i,j,kc,QPRES)
             qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qzpo(i,j,kc,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,kc,QPRES)
          qzpo(i,j,kc,qreitot) = sum(qzpo(i,j,kc,qrad:qradhi)) + qzpo(i,j,kc,QREINT)
#endif

          !-------------------------------------------------------------------
          ! qzmo states
          !-------------------------------------------------------------------

          pggp  = qy(i,j+1,km,GDPRES)
          pggm  = qy(i,j  ,km,GDPRES)
          ugp  = qy(i,j+1,km,GDV   )
          ugm  = qy(i,j  ,km,GDV   )
          gegp = qy(i,j+1,km,GDGAME)
          gegm = qy(i,j  ,km,GDGAME)
#ifdef RADIATION
          lambda = lam(i,j,k3d-1,:)
          ugc = HALF*(ugp+ugm)
          ergp = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergm = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmom = - cdtdy*sum(lamge(:))
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdy * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          ! Convert to conservation form
          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
               + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
#ifdef RADIATION
          erl  = qzm(i,j,kc,qrad:qradhi)
#endif

          ! Add transverse predictor
          rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
          runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX) - fy(i,j,km,UMX))
          rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY) - fy(i,j,km,UMY))
          rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ) - fy(i,j,km,UMZ))
          renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN)- fy(i,j,km,UEDEN))
#ifdef RADIATION
          rvnewlz = rvnewlz + dmom
          renewlz = renewlz + dre
          ernewl  = erl(:) - cdtdy*(rfy(i,j+1,km,:)- rfy(i,j,km,:)) &
               + der
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewlz < ZERO) then
             rrnewlz = rrlz
             runewlz = rulz
             rvnewlz = rvlz
             rwnewlz = rwlz
             renewlz = relz
#ifdef RADIATION
             ernewl  = erl(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qzmo(i,j,kc,QRHO) = rrnewlz
          rhoinv = ONE/rrnewlz
          qzmo(i,j,kc,QU) = runewlz*rhoinv
          qzmo(i,j,kc,QV) = rvnewlz*rhoinv
          qzmo(i,j,kc,QW) = rwnewlz*rhoinv

          ! note: we run the risk of (rho e) being negative here
          rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qzmo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qzmo(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) - &
                     cdtdy*(fy(i,j+1,km,UEINT) - fy(i,j,km,UEINT) + pav*du)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)
             else
                ! Update gammae with its transverse terms
                qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT)*(qzmo(i,j,kc,QGAME)-ONE)
                qzmo(i,j,kc,QPRES) = max(qzmo(i,j,kc,QPRES), small_pres)
             endif
          else
             qzmo(i,j,kc,QPRES) = qzm(i,j,kc,QPRES)
             qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME)
          endif

          call reset_edge_state_thermo(qzmo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qzmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qzmo(i,j,kc,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,kc,QPRES)
          qzmo(i,j,kc,qreitot) = sum(qzmo(i,j,kc,qrad:qradhi)) + qzmo(i,j,kc,QREINT)
#endif

       enddo
    enddo

  end subroutine transy2


  !===========================================================================
  ! transz
  !===========================================================================
  subroutine transz( &
#ifdef RADIATION
                    lam, lam_lo, lam_hi, &
#endif
                    qxm,qxmo,qxp,qxpo,qym,qymo,qyp,qypo,qd_lo,qd_hi, &
                    fz, &
#ifdef RADIATION
                    rfz, &
#endif
                    fz_lo,fz_hi, &
                    qz,qz_lo,qz_hi, &
                    gamc,gd_lo,gd_hi, &
                    cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fz_lo(3),fz_hi(3)
    integer :: qz_lo(3),qz_hi(3)
    integer :: gd_lo(3),gd_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    double precision  qxm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qxp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qym(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qyp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qxmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qxpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qymo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qypo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    double precision qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision cdtdz

    integer n, nq, i, j, ipassive

    double precision rrnew, rr
    double precision compn, compu
    double precision rrrx, rrry, rrlx, rrly
    double precision rurx, rury, rulx, ruly
    double precision rvrx, rvry, rvlx, rvly
    double precision rwrx, rwry, rwlx, rwly
    double precision ekenrx, ekenry, ekenlx, ekenly
    double precision rerx, rery, relx, rely
    double precision rrnewrx, rrnewry, rrnewlx, rrnewly
    double precision runewrx, runewry, runewlx, runewly
    double precision rvnewrx, rvnewry, rvnewlx, rvnewly
    double precision rwnewrx, rwnewry, rwnewlx, rwnewly
    double precision renewrx, renewry, renewlx, renewly
    double precision pnewrx, pnewry, pnewlx, pnewly
    double precision rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
    double precision pggp, pggm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

#ifdef RADIATION
    double precision :: dmz, dre
    double precision, dimension(0:ngroups-1) :: der, lambda, luge, lamge, &
         ergp, errx, ernewrx, erry, ernewry, ergm, erlx, ernewlx, erly, ernewly
    double precision eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))

             if (i >= ilo+1) then
                rr = qxp(i,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qxp(i,j,km,nq) - compn
                qxpo(i,j,km,nq) = compu/rrnew
             end if

             if (j >= jlo+1) then
                rr = qyp(i,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qyp(i,j,km,nq) - compn
                qypo(i,j,km,nq) = compu/rrnew
             end if

             if (i <= ihi-1) then
                rr = qxm(i+1,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qxm(i+1,j,km,nq) - compn
                qxmo(i+1,j,km,nq) = compu/rrnew
             end if

             if (j <= jhi-1) then
                rr = qym(i,j+1,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qym(i,j+1,km,nq) - compn
                qymo(i,j+1,km,nq) = compu/rrnew
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add transverse flux difference in the z-direction to the x- and
          ! y-states for the fluid variables
          !-------------------------------------------------------------------

          pggp  = qz(i,j,kc,GDPRES)
          pggm  = qz(i,j,km,GDPRES)
          ugp  = qz(i,j,kc,GDW   )
          ugm  = qz(i,j,km,GDW   )
          gegp = qz(i,j,kc,GDGAME)
          gegm = qz(i,j,km,GDGAME)
#ifdef RADIATION
          lambda = lam(i,j,k3d-1,:)
          ergp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

#ifdef RADIATION
          lamge = lambda(:) * (ergp(:)-ergm(:))
          dmz = - cdtdz*sum(lamge)
          luge = HALF*(ugp+ugm) * lamge(:)
          dre = -cdtdz*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdz*HALF*(ugp+ugm)*(ergp(g)-ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = cdtdz*HALF*(ergp(g)+ergm(g))*(ugm-ugp)*f1
             end do
          else ! mixed frame
             der(:) = cdtdz * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------

          if (i >= ilo+1) then
             ! Convert to conservation form
             rrrx = qxp(i,j,km,QRHO)
             rurx = rrrx*qxp(i,j,km,QU)
             rvrx = rrrx*qxp(i,j,km,QV)
             rwrx = rrrx*qxp(i,j,km,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,km,QU)**2 + qxp(i,j,km,QV)**2 &
                  + qxp(i,j,km,QW)**2)
             rerx = qxp(i,j,km,QREINT) + ekenrx
#ifdef RADIATION
             errx = qxp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewrx = rrrx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewrx = rurx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewrx = rvrx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewrx = rwrx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewrx = rerx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewrx = rwnewrx + dmz
             renewrx = renewrx + dre
             ernewrx = errx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewrx < ZERO) then
                rrnewrx = rrrx
                runewrx = rurx
                rvnewrx = rvrx
                rwnewrx = rwrx
                renewrx = rerx
#ifdef RADIATION
                ernewrx = errx(:)
#endif
                reset_state = .true.
             end if

             qxpo(i,j,km,QRHO) = rrnewrx
             qxpo(i,j,km,QU) = runewrx/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QV) = rvnewrx/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QW) = rwnewrx/qxpo(i,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QREINT) = renewrx - rhoekenrx

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qxpo(i,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxpo(i,j,km,QREINT) = qxp(i,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                   qxpo(i,j,km,QPRES) = max(pnewrx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxpo(i,j,km,QGAME) = qxp(i,j,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxpo(i,j,km,QPRES) = qxpo(i,j,km,QREINT)*(qxpo(i,j,km,QGAME)-ONE)
                   qxpo(i,j,km,QPRES) = max(qxpo(i,j,km,QPRES), small_pres)
                endif
             else
                qxpo(i,j,km,QPRES) = qxp(i,j,km,QPRES)
                qxpo(i,j,km,QGAME) = qxp(i,j,km,QGAME)
             endif

             call reset_edge_state_thermo(qxpo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qxpo(i,j,km,qrad:qradhi) = ernewrx(:)
             qxpo(i,j,km,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,km,QPRES)
             qxpo(i,j,km,qreitot) = sum(qxpo(i,j,km,qrad:qradhi)) + qxpo(i,j,km,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          if (j >= jlo+1) then
             ! Convert to conservation form
             rrry = qyp(i,j,km,QRHO)
             rury = rrry*qyp(i,j,km,QU)
             rvry = rrry*qyp(i,j,km,QV)
             rwry = rrry*qyp(i,j,km,QW)
             ekenry = HALF*rrry*(qyp(i,j,km,QU)**2 + qyp(i,j,km,QV)**2 &
                  + qyp(i,j,km,QW)**2)
             rery = qyp(i,j,km,QREINT) + ekenry
#ifdef RADIATION
             erry = qyp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewry = rrry - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewry = rury - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewry = rvry - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewry = rwry - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewry = rery - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewry = rwnewry + dmz
             renewry = renewry + dre
             ernewry = erry(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewry < ZERO) then
                rrnewry = rrry
                runewry = rury
                rvnewry = rvry
                rwnewry = rwry
                renewry = rery
#ifdef RADIATION
                ernewry = erry(:)
#endif
                reset_state = .true.
             end if

             qypo(i,j,km,QRHO) = rrnewry
             qypo(i,j,km,QU) = runewry/qypo(i,j,km,QRHO)
             qypo(i,j,km,QV) = rvnewry/qypo(i,j,km,QRHO)
             qypo(i,j,km,QW) = rwnewry/qypo(i,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,km,QRHO)
             qypo(i,j,km,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qypo(i,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qypo(i,j,km,QREINT) = qyp(i,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                   qypo(i,j,km,QPRES) = max(pnewry,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qypo(i,j,km,QGAME) = qyp(i,j,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qypo(i,j,km,QPRES) = qypo(i,j,km,QREINT)*(qypo(i,j,km,QGAME)-ONE)
                   qypo(i,j,km,QPRES) = max(qypo(i,j,km,QPRES), small_pres)
                endif
             else
                qypo(i,j,km,QPRES) = qyp(i,j,km,QPRES)
                qypo(i,j,km,QGAME) = qyp(i,j,km,QGAME)
             endif

             call reset_edge_state_thermo(qypo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qypo(i,j,km,qrad:qradhi) = ernewry(:)
             qypo(i,j,km,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,km,QPRES)
             qypo(i,j,km,qreitot) = sum(qypo(i,j,km,qrad:qradhi)) + qypo(i,j,km,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i <= ihi-1) then

             ! Convert to conservation form
             rrlx = qxm(i+1,j,km,QRHO)
             rulx = rrlx*qxm(i+1,j,km,QU)
             rvlx = rrlx*qxm(i+1,j,km,QV)
             rwlx = rrlx*qxm(i+1,j,km,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,km,QU)**2 + qxm(i+1,j,km,QV)**2 &
                  + qxm(i+1,j,km,QW)**2)
             relx = qxm(i+1,j,km,QREINT) + ekenlx
#ifdef RADIATION
             erlx = qxm(i+1,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlx = rrlx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewlx = rulx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewlx = rvlx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewlx = rwlx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewlx = relx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewlx = rwnewlx + dmz
             renewlx = renewlx + dre
             ernewlx = erlx(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewlx < ZERO) then
                rrnewlx = rrlx
                runewlx = rulx
                rvnewlx = rvlx
                rwnewlx = rwlx
                renewlx = relx
#ifdef RADIATION
                ernewlx = erlx(:)
#endif
                reset_state = .true.
             end if

             qxmo(i+1,j,km,QRHO) = rrnewlx
             qxmo(i+1,j,km,QU) = runewlx/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QV) = rvnewlx/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QW) = rwnewlx/qxmo(i+1,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo(i+1,j,km,QREINT) = qxm(i+1,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                   qxmo(i+1,j,km,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo(i+1,j,km,QGAME) = qxm(i+1,j,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo(i+1,j,km,QPRES) = qxmo(i+1,j,km,QREINT)*(qxmo(i+1,j,km,QGAME)-ONE)
                   qxmo(i+1,j,km,QPRES) = max(qxmo(i+1,j,km,QPRES), small_pres)
                end if
             else
                qxmo(i+1,j,km,QPRES) = qxm(i+1,j,km,QPRES)
                qxmo(i+1,j,km,QGAME) = qxm(i+1,j,km,QGAME)
             endif

             call reset_edge_state_thermo(qxmo, qd_lo, qd_hi, i+1, j, km)

#ifdef RADIATION
             qxmo(i+1,j,km,qrad:qradhi) = ernewlx(:)
             qxmo(i+1,j,km,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i+1,j,km,QPRES)
             qxmo(i+1,j,km,qreitot) = sum(qxmo(i+1,j,km,qrad:qradhi)) + qxmo(i+1,j,km,QREINT)
#endif

          endif


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------

          if (j <= jhi-1) then
             ! Convert to conservation form
             rrly = qym(i,j+1,km,QRHO)
             ruly = rrly*qym(i,j+1,km,QU)
             rvly = rrly*qym(i,j+1,km,QV)
             rwly = rrly*qym(i,j+1,km,QW)
             ekenly = HALF*rrly*(qym(i,j+1,km,QU)**2 + qym(i,j+1,km,QV)**2 &
                  + qym(i,j+1,km,QW)**2)
             rely = qym(i,j+1,km,QREINT) + ekenly
#ifdef RADIATION
             erly = qym(i,j+1,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewly = rrly - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewly = ruly - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewly = rvly - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewly = rwly - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewly = rely - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
#ifdef RADIATION
             rwnewly = rwnewly + dmz
             renewly = renewly + dre
             ernewly = erly(:) - cdtdz*(rfz(i,j,kc,:) - rfz(i,j,km,:)) &
                  + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewly < ZERO) then
                rrnewly = rrly
                runewly = ruly
                rvnewly = rvly
                rwnewly = rwly
                renewly = rely
#ifdef RADIATION
                ernewly = erly(:)
#endif
                reset_state = .true.
             endif

             qymo(i,j+1,km,QRHO) = rrnewly
             qymo(i,j+1,km,QU) = runewly/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QV) = rvnewly/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QW) = rwnewly/qymo(i,j+1,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qymo(i,j+1,km,QREINT) = qym(i,j+1,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j+1,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                   qymo(i,j+1,km,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qymo(i,j+1,km,QGAME) = qym(i,j+1,km,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qymo(i,j+1,km,QPRES) = qymo(i,j+1,km,QREINT)*(qymo(i,j+1,km,QGAME)-ONE)
                   qymo(i,j+1,km,QPRES) = max(qymo(i,j+1,km,QPRES), small_pres)
                endif
             else
                qymo(i,j+1,km,QPRES) = qym(i,j+1,km,QPRES)
                qymo(i,j+1,km,QGAME) = qym(i,j+1,km,QGAME)
             endif

             call reset_edge_state_thermo(qymo, qd_lo, qd_hi, i, j+1, km)

#ifdef RADIATION
             qymo(i,j+1,km,qrad:qradhi) = ernewly(:)
             qymo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j+1,km,QPRES)
             qymo(i,j+1,km,qreitot) = sum(qymo(i,j+1,km,qrad:qradhi)) + qymo(i,j+1,km,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transz


  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qm,qmo,qp,qpo,qd_lo,qd_hi, &
                     fxy, &
#ifdef RADIATION
                     rfxy, &
#endif
                     fx_lo,fx_hi, &
                     fyx, &
#ifdef RADIATION
                     rfyx, &
#endif
                     fy_lo,fy_hi, &
                     qx,qx_lo,qx_hi, &
                     qy,qy_lo,qy_hi, &
                     gamc,gd_lo,gd_hi, &
                     srcQ,src_lo,src_hi, &
                     hdt,cdtdx,cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer :: gd_lo(3),gd_hi(3)
    integer :: src_lo(3),src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    double precision rfyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    double precision  qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision fyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision  qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    double precision  qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    double precision gamc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2),gd_lo(3):gd_hi(3))
    double precision srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision hdt,cdtdx,cdtdy

    integer i, j, n, nq, ipassive

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pggxp, pggxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    double precision pggyp, pggym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    double precision uxav, gexav, dgex, uyav, geyav, dgey
    double precision pggxpm, pggxmm, ugxpm, ugxmm, gegxpm, gegxmm, duxpm, pxavm, duxm, pxnewm, gexnewm
    double precision pggypm, pggymm, ugypm, ugymm, gegypm, gegymm, duypm, pyavm, duym, pynewm, geynewm
    double precision uxavm, gexavm, dgexm, uyavm, geyavm, dgeym
    double precision compr, compl, compnr, compnl

#ifdef RADIATION
    double precision :: dmx, dmy, dre
    double precision, dimension(0:ngroups-1) :: der, lamc, lamm, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym, ergxpm, ergypm, ergxmm, ergymm
    double precision eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             rrr = qp(i,j,kc,QRHO)
             rrl = qm(i,j,kc,QRHO)

             compr = rrr*qp(i,j,kc,nq)
             compl = rrl*qm(i,j,kc,nq)

             rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                          - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
             rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                          - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                            - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
             compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                            - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))

             qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d  ,nq)
             qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse xy and yx differences to the z-states for the
          ! fluid variables
          !-------------------------------------------------------------------

          pggxp  = qx(i+1,j,kc,GDPRES)
          pggxm  = qx(i  ,j,kc,GDPRES)
          ugxp  = qx(i+1,j,kc,GDU   )
          ugxm  = qx(i  ,j,kc,GDU   )
          gegxp = qx(i+1,j,kc,GDGAME)
          gegxm = qx(i  ,j,kc,GDGAME)

          pggyp  = qy(i,j+1,kc,GDPRES)
          pggym  = qy(i,j  ,kc,GDPRES)
          ugyp  = qy(i,j+1,kc,GDV   )
          ugym  = qy(i,j  ,kc,GDV   )
          gegyp = qy(i,j+1,kc,GDGAME)
          gegym = qy(i,j  ,kc,GDGAME)

#ifdef RADIATION
          lamc = lam(i,j,k3d,:)
          ergxp = qx(i+1,j,kc,GDERADS:GDERADS-1+ngroups)
          ergxm = qx(i  ,j,kc,GDERADS:GDERADS-1+ngroups)
          ergyp = qy(i,j+1,kc,GDERADS:GDERADS-1+ngroups)
          ergym = qy(i,j  ,kc,GDERADS:GDERADS-1+ngroups)
#endif

          pggxpm  = qx(i+1,j,km,GDPRES)
          pggxmm  = qx(i  ,j,km,GDPRES)
          ugxpm  = qx(i+1,j,km,GDU   )
          ugxmm  = qx(i  ,j,km,GDU   )
          gegxpm = qx(i+1,j,km,GDGAME)
          gegxmm = qx(i  ,j,km,GDGAME)

          pggypm  = qy(i,j+1,km,GDPRES)
          pggymm  = qy(i,j  ,km,GDPRES)
          ugypm  = qy(i,j+1,km,GDV   )
          ugymm  = qy(i,j  ,km,GDV   )
          gegypm = qy(i,j+1,km,GDGAME)
          gegymm = qy(i,j  ,km,GDGAME)

#ifdef RADIATION
          lamm = lam(i,j,k3d-1,:)
          ergxpm = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergxmm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)
          ergypm = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergymm = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)
#endif

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = HALF*(pggxp+pggxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav-gamc(i,j,k3d))*dux - uxav*dgex )

          duxpm = pggxpm*ugxpm - pggxmm*ugxmm
          pxavm = HALF*(pggxpm+pggxmm)
          uxavm = HALF*(ugxpm+ugxmm)
          gexavm = HALF*(gegxpm+gegxmm)
          duxm = ugxpm-ugxmm
          dgexm = gegxpm-gegxmm
          pxnewm = cdtdx*(duxpm + pxavm*duxm*(gamc(i,j,k3d-1)-ONE))
          gexnewm = cdtdx*( (gexavm-ONE)*(gexavm-gamc(i,j,k3d-1))*duxm - uxavm*dgexm )

          duyp = pggyp*ugyp - pggym*ugym
          pyav = HALF*(pggyp+pggym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav-gamc(i,j,k3d))*duy - uyav*dgey )

          duypm = pggypm*ugypm - pggymm*ugymm
          pyavm = HALF*(pggypm+pggymm)
          uyavm = HALF*(ugypm+ugymm)
          geyavm = HALF*(gegypm+gegymm)
          duym = ugypm-ugymm
          dgeym = gegypm-gegymm
          pynewm = cdtdy*(duypm + pyavm*duym*(gamc(i,j,k3d-1)-ONE))
          geynewm = cdtdy*( (geyavm-ONE)*(geyavm-gamc(i,j,k3d-1))*duym - uyavm*dgeym )


          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------

          ! Convert to conservation form
          rrr = qp(i,j,kc,QRHO)
          rur = rrr*qp(i,j,kc,QU)
          rvr = rrr*qp(i,j,kc,QV)
          rwr = rrr*qp(i,j,kc,QW)
          ekenr = HALF*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
               qp(i,j,kc,QW)**2)
          rer = qp(i,j,kc,QREINT) + ekenr
#ifdef RADIATION
          err = qp(i,j,kc,qrad:qradhi)

          lgex = lamc(:) * (ergxp(:)-ergxm(:))
          lgey = lamc(:) * (ergyp(:)-ergym(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          lugex = HALF*(ugxp+ugxm) * lgex(:)
          lugey = HALF*(ugyp+ugym) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamc(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if
#endif

          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                       - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
          runewr = rur - cdtdx*(fxy(i+1,j,kc,UMX) - fxy(i,j,kc,UMX)) &
                       - cdtdy*(fyx(i,j+1,kc,UMX) - fyx(i,j,kc,UMX))
          rvnewr = rvr - cdtdx*(fxy(i+1,j,kc,UMY) - fxy(i,j,kc,UMY)) &
                       - cdtdy*(fyx(i,j+1,kc,UMY) - fyx(i,j,kc,UMY))
          rwnewr = rwr - cdtdx*(fxy(i+1,j,kc,UMZ) - fxy(i,j,kc,UMZ)) &
                       - cdtdy*(fyx(i,j+1,kc,UMZ) - fyx(i,j,kc,UMZ))
          renewr = rer - cdtdx*(fxy(i+1,j,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
                       - cdtdy*(fyx(i,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))
#ifdef RADIATION
          runewr = runewr + dmx
          rvnewr = rvnewr + dmy
          renewr = renewr + dre
          ernewr = err(:) - cdtdx*(rfxy(i+1,j,kc,:) - rfxy(i,j,kc,:)) &
               &          - cdtdy*(rfyx(i,j+1,kc,:) - rfyx(i,j,kc,:))  &
               &          + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
             rrnewr = rrr
             runewr = rur
             rvnewr = rvr
             rwnewr = rwr
             renewr = rer
#ifdef RADIATION
             ernewr = err(:)
#endif
             reset_state = .true.
          endif

          ! Convert back to primitive form
          qpo(i,j,kc,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,kc,QU    ) = runewr/rrnewr
          qpo(i,j,kc,QV    ) = rvnewr/rrnewr
          qpo(i,j,kc,QW    ) = rwnewr/rrnewr

          ! if ppm_trace_sources == 1, then we already added the piecewise parabolic traced
          ! source terms to the normal edge states.
          if (ppm_trace_sources == 0 .or. ppm_type == 0) then
             qpo(i,j,kc,QU:QW) = qpo(i,j,kc,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
          endif

          ! note: we run the risk of (rho e) being negative here
          rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
          qpo(i,j,kc,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qpo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qpo(i,j,kc,QREINT) = qp(i,j,kc,QREINT) &
                     - cdtdx*(fxy(i+1,j,kc,UEINT) - fxy(i,j,kc,UEINT) + pxav*dux) &
                     - cdtdy*(fyx(i,j+1,kc,UEINT) - fyx(i,j,kc,UEINT) + pyav*duy) &
                     + hdt*srcQ(i,j,k3d,QREINT)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
                qpo(i,j,kc,QPRES) = pnewr + hdt*srcQ(i,j,k3d,QPRES)
             else
                ! Update gammae with its transverse terms
                qpo(i,j,kc,QGAME) = qp(i,j,kc,QGAME) + gexnew + geynew

                ! and compute the p edge state from this and (rho e)
                qpo(i,j,kc,QPRES) = qpo(i,j,kc,QREINT)*(qpo(i,j,kc,QGAME)-ONE)
             endif
          else
             qpo(i,j,kc,QPRES) = qp(i,j,kc,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
             qpo(i,j,kc,QGAME) = qp(i,j,kc,QGAME)
          endif

          qpo(i,j,kc,QPRES) = max(qpo(i,j,kc,QPRES), small_pres)

          call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qpo(i,j,kc,qrad:qradhi) = ernewr(:)
          qpo(i,j,kc,qptot  ) = sum(lamc(:)*ernewr(:)) + qpo(i,j,kc,QPRES)
          qpo(i,j,kc,qreitot) = sum(qpo(i,j,kc,qrad:qradhi)) + qpo(i,j,kc,QREINT)
#endif

          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------

          ! Convert to conservation form
          rrl = qm(i,j,kc,QRHO)
          rul = rrl*qm(i,j,kc,QU)
          rvl = rrl*qm(i,j,kc,QV)
          rwl = rrl*qm(i,j,kc,QW)
          ekenl = HALF*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
               qm(i,j,kc,QW)**2)
          rel = qm(i,j,kc,QREINT) + ekenl
#ifdef RADIATION
          erl = qm(i,j,kc,qrad:qradhi)

          lgex = lamm(:) * (ergxpm(:)-ergxmm(:))
          lgey = lamm(:) * (ergypm(:)-ergymm(:))
          dmx = - cdtdx*sum(lgex)
          dmy = - cdtdy*sum(lgey)
          lugex = HALF*(ugxpm+ugxmm) * lgex(:)
          lugey = HALF*(ugypm+ugymm) * lgey(:)
          dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ugxpm+ugxmm)*(ergxpm(g)-ergxmm(g)) &
                     +       cdtdy*HALF*(ugypm+ugymm)*(ergypm(g)-ergymm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lamm(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ergxpm(g)+ergxmm(g))*(ugxmm-ugxpm) &
                     +       cdtdy*HALF*(ergypm(g)+ergymm(g))*(ugymm-ugypm) )
             end do
          else ! mixed frame
             der(:) = cdtdx * lugex + cdtdy * lugey
          end if
#endif

          ! Add transverse predictor
          rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                       - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxy(i+1,j,km,UMX) - fxy(i,j,km,UMX)) &
                       - cdtdy*(fyx(i,j+1,km,UMX) - fyx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxy(i+1,j,km,UMY) - fxy(i,j,km,UMY)) &
                       - cdtdy*(fyx(i,j+1,km,UMY) - fyx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxy(i+1,j,km,UMZ) - fxy(i,j,km,UMZ)) &
                       - cdtdy*(fyx(i,j+1,km,UMZ) - fyx(i,j,km,UMZ))
          renewl = rel - cdtdx*(fxy(i+1,j,km,UEDEN) - fxy(i,j,km,UEDEN)) &
                       - cdtdy*(fyx(i,j+1,km,UEDEN) - fyx(i,j,km,UEDEN))
#ifdef RADIATION
          runewl = runewl + dmx
          rvnewl = rvnewl + dmy
          renewl = renewl + dre
          ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,km,:) - rfxy(i,j,km,:)) &
               &          - cdtdy*(rfyx(i  ,j+1,km,:) - rfyx(i,j,km,:)) &
               &          + der(:)
#endif

          ! Reset to original value if adding transverse terms made density negative
          reset_state = .false.
          if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
             rrnewl = rrl
             runewl = rul
             rvnewl = rvl
             rwnewl = rwl
             renewl = rel
#ifdef RADIATION
             ernewl = erl(:)
#endif
             reset_state = .true.
          endif

          qmo(i,j,kc,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k3d-1,QRHO)
          qmo(i,j,kc,QU    ) = runewl/rrnewl
          qmo(i,j,kc,QV    ) = rvnewl/rrnewl
          qmo(i,j,kc,QW    ) = rwnewl/rrnewl

          ! if ppm_trace_sources == 1, then we already added the piecewise parabolic traced
          ! source terms to the normal edge states.
          if (ppm_trace_sources == 0 .or. ppm_type == 0) then
             qmo(i,j,kc,QU:QW) = qmo(i,j,kc,QU:QW) + hdt * srcQ(i,j,k3d-1,QU:QW)
          endif

          ! note: we run the risk of (rho e) being negative here
          rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
          qmo(i,j,kc,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d-1,QREINT)

          if (.not. reset_state) then
             if (transverse_reset_rhoe == 1 .and. qmo(i,j,kc,QREINT) <= ZERO) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                qmo(i,j,kc,QREINT) = qm(i,j,kc,QREINT) &
                     - cdtdx*(fxy(i+1,j,km,UEINT) - fxy(i,j,km,UEINT) + pxavm*duxm) &
                     - cdtdy*(fyx(i,j+1,km,UEINT) - fyx(i,j,km,UEINT) + pyavm*duym) &
                     + hdt*srcQ(i,j,k3d-1,QREINT)
             endif

             ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
             ! If we are wrong, we will fix it later

             if (ppm_predict_gammae == 0) then
                ! add the transverse term to the p evolution eq here
                pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm
                qmo(i,j,kc,QPRES) = pnewl + hdt*srcQ(i,j,k3d-1,QPRES)
             else
                ! Update gammae with its transverse terms
                qmo(i,j,kc,QGAME) = qm(i,j,kc,QGAME) + gexnewm + geynewm

                ! and compute the p edge state from this and (rho e)
                qmo(i,j,kc,QPRES) = qmo(i,j,kc,QREINT)*(qmo(i,j,kc,QGAME)-ONE)
             endif
          else
             qmo(i,j,kc,QPRES) = qm(i,j,kc,QPRES) + hdt*srcQ(i,j,k3d-1,QPRES)
             qmo(i,j,kc,QGAME) = qm(i,j,kc,QGAME)
          endif

          qmo(i,j,kc,QPRES) = max(qmo(i,j,kc,QPRES), small_pres)

          call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i, j, kc)

#ifdef RADIATION
          qmo(i,j,kc,qrad:qradhi) = ernewl(:)
          qmo(i,j,kc,qptot  ) = sum(lamm(:)*ernewl(:)) + qmo(i,j,kc,QPRES)
          qmo(i,j,kc,qreitot) = sum(qmo(i,j,kc,qrad:qradhi)) + qmo(i,j,kc,QREINT)
#endif

       enddo
    enddo

  end subroutine transxy


  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qm,qmo,qp,qpo,qd_lo,qd_hi, &
                     fxz, &
#ifdef RADIATION
                     rfxz, &
#endif
                     fx_lo,fx_hi, &
                     fzx, &
#ifdef RADIATION
                     rfzx, &
#endif
                     fz_lo,fz_hi, &
                     qx,qx_lo,qx_hi, &
                     qz,qz_lo,qz_hi, &
                     gamc,gc_lo,gc_hi, &
                     srcQ,src_lo,src_hi, &
                     hdt,cdtdx,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

#ifdef RADIATION
    integer :: lam_lo(3), lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fx_lo(3),fx_hi(3)
    integer :: fz_lo(3),fz_hi(3)
    integer :: qx_lo(3),qx_hi(3)
    integer :: qz_lo(3),qz_hi(3)
    integer :: gc_lo(3),gc_hi(3)
    integer :: src_lo(3),src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    double precision rfzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    double precision  qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision  qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    double precision fzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    double precision  qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    double precision  qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    double precision gamc(gc_lo(1):gc_hi(1),gc_lo(2):gc_hi(2),gc_lo(3):gc_hi(3))
    double precision srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision hdt,cdtdx,cdtdz

    integer i, j, n, nq, ipassive

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pggxp, pggxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    double precision pggzp, pggzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    double precision uxav, gexav, dgex, uzav, gezav, dgez
    double precision compr, compl, compnr, compnl, drr, dcompn

#ifdef RADIATION
    double precision :: dmx, dmz, dre
    double precision, dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm
    double precision eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             drr    = - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                      - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             dcompn = - cdtdx*(fxz(i+1,j,km,n   ) - fxz(i,j,km,n)) &
                      - cdtdz*(fzx(i  ,j,kc,n   ) - fzx(i,j,km,n))

             if (j >= jlo+1) then
                rrr = qp(i,j,km,QRHO)
                compr = rrr*qp(i,j,km,nq)

                rrnewr = rrr + drr
                compnr = compr + dcompn

                qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             end if

             if (j <= jhi-1) then
                rrl = qm(i,j+1,km,QRHO)
                compl = rrl*qm(i,j+1,km,nq)

                rrnewl = rrl + drr
                compnl = compl + dcompn

                qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse xz and zx differences to the y-states for the
          ! fluid variables
          !-------------------------------------------------------------------

#ifdef RADIATION
          lambda = lam(i,j,k3d,:)
#endif

          pggxp  = qx(i+1,j,km,GDPRES)
          pggxm  = qx(i  ,j,km,GDPRES)
          ugxp  = qx(i+1,j,km,GDU   )
          ugxm  = qx(i  ,j,km,GDU   )
          gegxp = qx(i+1,j,km,GDGAME)
          gegxm = qx(i  ,j,km,GDGAME)
#ifdef RADIATION
          ergxp = qx(i+1,j,km,GDERADS:GDERADS-1+ngroups)
          ergxm = qx(i  ,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          pggzp  = qz(i,j,kc,GDPRES)
          pggzm  = qz(i,j,km,GDPRES)
          ugzp  = qz(i,j,kc,GDW   )
          ugzm  = qz(i,j,km,GDW   )
          gegzp = qz(i,j,kc,GDGAME)
          gegzm = qz(i,j,km,GDGAME)
#ifdef RADIATION
          ergzp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergzm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          duxp = pggxp*ugxp - pggxm*ugxm
          pxav = HALF*(pggxp+pggxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav-gamc(i,j,k3d))*dux - uxav*dgex )

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = HALF*(pggzp+pggzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav-gamc(i,j,k3d))*duz - uzav*dgez )

#ifdef RADIATION
          lgex = lambda(:) * (ergxp(:)-ergxm(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmx = - cdtdx*sum(lgex)
          dmz = - cdtdz*sum(lgez)
          lugex = HALF*(ugxp+ugxm) * lgex(:)
          lugez = HALF*(ugzp+ugzm) * lgez(:)
          dre = -cdtdx * sum(lugex) - cdtdz * sum(lugez)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                     +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                     +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdx*lugex + cdtdz*lugez
          end if
#endif

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          if (j >= jlo+1) then
             ! Convert to conservation form
             rrr = qp(i,j,km,QRHO)
             rur = rrr*qp(i,j,km,QU)
             rvr = rrr*qp(i,j,km,QV)
             rwr = rrr*qp(i,j,km,QW)
             ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
             rer = qp(i,j,km,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                          - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
             runewr = rur - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                          - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
             rvnewr = rvr - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                          - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
             rwnewr = rwr - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                          - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
             renewr = rer - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                          - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))
#ifdef RADIATION
             runewr = runewr + dmx
             rwnewr = rwnewr + dmz
             renewr = renewr + dre
             ernewr = err(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
                             - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             endif

             qpo(i,j,km,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
             qpo(i,j,km,QU    ) = runewr/rrnewr
             qpo(i,j,km,QV    ) = rvnewr/rrnewr
             qpo(i,j,km,QW    ) = rwnewr/rrnewr

             ! if ppm_trace_sources == 1, then we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qpo(i,j,km,QU:QW) = qpo(i,j,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,km,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,km,QREINT) <= ZERO) then
                   qpo(i,j,km,QREINT) = qp(i,j,km,QREINT) &
                        - cdtdx*(fxz(i+1,j,km,UEINT) - fxz(i,j,km,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i  ,j,kc,UEINT) - fzx(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,km,QPRES) - pxnew - pznew
                   qpo(i,j,km,QPRES) = pnewr + hdt*srcQ(i,j,k3d,QPRES)
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,km,QGAME) = qp(i,j,km,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT)*(qpo(i,j,km,QGAME)-ONE)
                endif
             else
                qpo(i,j,km,QPRES) = qp(i,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                qpo(i,j,km,QGAME) = qp(i,j,km,QGAME)
             endif

             qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qpo(i,j,km,qrad:qradhi) = ernewr(:)
             qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
             qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)
#endif
          end if


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------

          if (j <= jhi-1) then
             ! Convert to conservation form
             rrl = qm(i,j+1,km,QRHO)
             rul = rrl*qm(i,j+1,km,QU)
             rvl = rrl*qm(i,j+1,km,QV)
             rwl = rrl*qm(i,j+1,km,QW)
             ekenl = HALF*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
             rel = qm(i,j+1,km,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i,j+1,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                          - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
             runewl = rul - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                          - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
             rvnewl = rvl - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                          - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
             rwnewl = rwl - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                          - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
             renewl = rel - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                          - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxz(i+1,j,km,:) - rfxz(i,j,km,:)) &
                             - cdtdz*(rfzx(i  ,j,kc,:) - rfzx(i,j,km,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
#ifdef RADIATION
                ernewl = erl(:)
#endif
                reset_state = .true.
             endif

             qmo(i,j+1,km,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k3d,QRHO)
             qmo(i,j+1,km,QU    ) = runewl/rrnewl
             qmo(i,j+1,km,QV    ) = rvnewl/rrnewl
             qmo(i,j+1,km,QW    ) = rwnewl/rrnewl

             ! if ppm_trace_sources == 1, then we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qmo(i,j+1,km,QU:QW) = qmo(i,j+1,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j+1,km,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j+1,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j+1,km,QREINT) = qm(i,j+1,km,QREINT) &
                        - cdtdx*(fxz(i+1,j,km,UEINT) - fxz(i,j,km,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i,j,kc,UEINT) - fzx(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew
                   qmo(i,j+1,km,QPRES) = pnewl + hdt*srcQ(i,j,k3d,QPRES)
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j+1,km,QGAME) = qm(i,j+1,km,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j+1,km,QPRES) = qmo(i,j+1,km,QREINT)*(qmo(i,j+1,km,QGAME)-ONE)
                endif
             else
                qmo(i,j+1,km,QPRES) = qm(i,j+1,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                qmo(i,j+1,km,QGAME) = qm(i,j+1,km,QGAME)
             endif

             qmo(i,j+1,km,QPRES) = max(qmo(i,j+1,km,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i, j+1, km)

#ifdef RADIATION
             qmo(i,j+1,km,qrad:qradhi) = ernewl(:)
             qmo(i,j+1,km,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j+1,km,QPRES)
             qmo(i,j+1,km,qreitot) = sum(qmo(i,j+1,km,qrad:qradhi)) + qmo(i,j+1,km,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transxz


  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz( &
#ifdef RADIATION
                     lam, lam_lo, lam_hi, &
#endif
                     qm,qmo,qp,qpo,qd_lo,qd_hi, &
                     fyz, &
#ifdef RADIATION
                     rfyz, &
#endif
                     fy_lo, fy_hi, &
                     fzy, &
#ifdef RADIATION
                     rfzy, &
#endif
                     fz_lo,fz_hi, &
                     qy,qy_lo,qy_hi, &
                     qz,qz_lo,qz_hi, &
                     gamc,gc_lo,gc_hi, &
                     srcQ,src_lo,src_hi, &
                     hdt,cdtdy,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

#ifdef RADIATION
    integer :: lam_lo(3),lam_hi(3)
#endif
    integer :: qd_lo(3),qd_hi(3)
    integer :: fy_lo(3),fy_hi(3)
    integer :: fz_lo(3),fz_hi(3)
    integer :: qy_lo(3),qy_hi(3)
    integer :: qz_lo(3),qz_hi(3)
    integer :: gc_lo(3),gc_hi(3)
    integer :: src_lo(3),src_hi(3)
    integer ilo,ihi,jlo,jhi,km,kc,k3d

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rfyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    double precision rfzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    double precision qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)
    double precision qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    double precision fyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    double precision fzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    double precision  qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    double precision  qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    double precision gamc(gc_lo(1):gc_hi(1),gc_lo(2):gc_hi(2),gc_lo(3):gc_hi(3))
    double precision srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
    double precision hdt,cdtdy,cdtdz

    integer i, j, n, nq, ipassive

    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pggyp, pggym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    double precision pggzp, pggzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    double precision uyav, geyav, dgey, uzav, gezav, dgez
    double precision compr, compl, compnr, compnl
    double precision drr, dcompn

#ifdef RADIATION
    double precision :: dmy, dmz, dre
    double precision, dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
    double precision eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi

             drr    = - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                      - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             dcompn = - cdtdy*(fyz(i,j+1,km,n   ) - fyz(i,j,km,n)) &
                      - cdtdz*(fzy(i,j  ,kc,n   ) - fzy(i,j,km,n))

             if (i >= ilo+1) then
                rrr = qp(i,j,km,QRHO)
                compr = rrr*qp(i,j,km,nq)

                rrnewr = rrr +drr
                compnr = compr +dcompn

                qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             end if

             if (i <= ihi-1) then
                rrl = qm(i+1,j,km,QRHO)
                compl = rrl*qm(i+1,j,km,nq)

                rrnewl = rrl + drr
                compnl = compl +dcompn

                qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)
             end if
          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse yz and zy differences to the x-states for the
          ! fluid variables
          !-------------------------------------------------------------------

#ifdef RADIATION
          lambda = lam(i,j,k3d,:)
#endif

          pggyp  = qy(i,j+1,km,GDPRES)
          pggym  = qy(i,j  ,km,GDPRES)
          ugyp  = qy(i,j+1,km,GDV   )
          ugym  = qy(i,j  ,km,GDV   )
          gegyp = qy(i,j+1,km,GDGAME)
          gegym = qy(i,j  ,km,GDGAME)
#ifdef RADIATION
          ergyp = qy(i,j+1,km,GDERADS:GDERADS-1+ngroups)
          ergym = qy(i,j  ,km,GDERADS:GDERADS-1+ngroups)
#endif

          pggzp  = qz(i,j,kc,GDPRES)
          pggzm  = qz(i,j,km,GDPRES)
          ugzp  = qz(i,j,kc,GDW   )
          ugzm  = qz(i,j,km,GDW   )
          gegzp = qz(i,j,kc,GDGAME)
          gegzm = qz(i,j,km,GDGAME)
#ifdef RADIATION
          ergzp = qz(i,j,kc,GDERADS:GDERADS-1+ngroups)
          ergzm = qz(i,j,km,GDERADS:GDERADS-1+ngroups)
#endif

          duyp = pggyp*ugyp - pggym*ugym
          pyav = HALF*(pggyp+pggym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav-gamc(i,j,k3d))*duy - uyav*dgey )

          duzp = pggzp*ugzp - pggzm*ugzm
          pzav = HALF*(pggzp+pggzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav-gamc(i,j,k3d))*duz - uzav*dgez )

#ifdef RADIATION
          lgey = lambda(:) * (ergyp(:)-ergym(:))
          lgez = lambda(:) * (ergzp(:)-ergzm(:))
          dmy = - cdtdy*sum(lgey)
          dmz = - cdtdz*sum(lgez)
          lugey = HALF*(ugyp+ugym) * lgey(:)
          lugez = HALF*(ugzp+ugzm) * lgez(:)
          dre = -cdtdy*sum(lugey) - cdtdz*sum(lugez)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) &
                     +       cdtdz*HALF*(ugzp+ugzm)*(ergzp(g)-ergzm(g)) )
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = HALF*(ONE-eddf)
                der(g) = f1*(cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) &
                     +       cdtdz*HALF*(ergzp(g)+ergzm(g))*(ugzm-ugzp) )
             end do
          else ! mixed frame
             der(:) = cdtdy*lugey + cdtdz*lugez
          end if
#endif

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------

          if (i >= ilo+1) then
             ! Convert to conservation form
             rrr = qp(i,j,km,QRHO)
             rur = rrr*qp(i,j,km,QU)
             rvr = rrr*qp(i,j,km,QV)
             rwr = rrr*qp(i,j,km,QW)
             ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
                  qp(i,j,km,QW)**2)
             rer = qp(i,j,km,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                          - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
             runewr = rur - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                          - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
             rvnewr = rvr - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                          - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
             rwnewr = rwr - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                          - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
             renewr = rer - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                          - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))
#ifdef RADIATION
             rvnewr = rvnewr + dmy
             rwnewr = rwnewr + dmz
             renewr = renewr + dre
             ernewr = err(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
                             - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
                             + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                rwnewr = rwr
                renewr = rer
#ifdef RADIATION
                ernewr = err(:)
#endif
                reset_state = .true.
             end if

             qpo(i,j,km,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
             qpo(i,j,km,QU    ) = runewr/rrnewr
             qpo(i,j,km,QV    ) = rvnewr/rrnewr
             qpo(i,j,km,QW    ) = rwnewr/rrnewr

             ! if ppm_trace_sources == 1, then we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qpo(i,j,km,QU:QW) = qpo(i,j,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,km,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qpo(i,j,km,QREINT) = qp(i,j,km,QREINT) &
                        - cdtdy*(fyz(i,j+1,km,UEINT) - fyz(i,j,km,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,kc,UEINT) - fzy(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,km,QPRES) - pynew - pznew
                   qpo(i,j,km,QPRES) = pnewr + hdt*srcQ(i,j,k3d,QPRES)
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,km,QGAME) = qp(i,j,km,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT)*(qpo(i,j,km,QGAME)-ONE)
                end if
             else
                qpo(i,j,km,QPRES) = qp(i,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                qpo(i,j,km,QGAME) = qp(i,j,km,QGAME)
             endif

             qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qd_lo, qd_hi, i, j, km)

#ifdef RADIATION
             qpo(i,j,km,qrad:qradhi) = ernewr(:)
             qpo(i,j,km,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,km,QPRES)
             qpo(i,j,km,qreitot) = sum(qpo(i,j,km,qrad:qradhi)) + qpo(i,j,km,QREINT)
#endif

          endif

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i <=ihi-1) then
             ! Convert to conservation form
             rrl = qm(i+1,j,km,QRHO)
             rul = rrl*qm(i+1,j,km,QU)
             rvl = rrl*qm(i+1,j,km,QV)
             rwl = rrl*qm(i+1,j,km,QW)
             ekenl = HALF*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
                  qm(i+1,j,km,QW)**2)
             rel = qm(i+1,j,km,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i+1,j,km,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                          - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
             runewl = rul - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                          - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
             rvnewl = rvl - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                          - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
             rwnewl = rwl - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                          - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
             renewl = rel - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                          - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))
#ifdef RADIATION
             rvnewl = rvnewl + dmy
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdy*(rfyz(i,j+1,km,:) - rfyz(i,j,km,:)) &
                  &          - cdtdz*(rfzy(i,j  ,kc,:) - rfzy(i,j,km,:)) &
                  &          + der(:)
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                rwnewl = rwl
                renewl = rel
#ifdef RADIATION
                ernewl = erl(:)
#endif
                reset_state = .true.
             endif

             qmo(i+1,j,km,QRHO   ) = rrnewl        + hdt*srcQ(i,j,k3d,QRHO)
             qmo(i+1,j,km,QU     ) = runewl/rrnewl
             qmo(i+1,j,km,QV     ) = rvnewl/rrnewl
             qmo(i+1,j,km,QW     ) = rwnewl/rrnewl

             ! if ppm_trace_sources == 1, then we already added the piecewise parabolic traced
             ! source terms to the normal edge states.
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qmo(i+1,j,km,QU:QW) = qmo(i+1,j,km,QU:QW) + hdt * srcQ(i,j,k3d,QU:QW)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i+1,j,km,QREINT ) = renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i+1,j,km,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i+1,j,km,QREINT ) = qm(i+1,j,km,QREINT) &
                        - cdtdy*(fyz(i,j+1,km,UEINT) - fyz(i,j,km,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,kc,UEINT) - fzy(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i+1,j,km,QPRES) - pynew - pznew
                   qmo(i+1,j,km,QPRES  ) = pnewl + hdt*srcQ(i,j,k3d,QPRES)
                else
                   ! Update gammae with its transverse terms
                   qmo(i+1,j,km,QGAME) = qm(i+1,j,km,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i+1,j,km,QPRES) = qmo(i+1,j,km,QREINT)*(qmo(i+1,j,km,QGAME)-ONE)
                end if
             else
                qmo(i+1,j,km,QPRES  ) = qm(i+1,j,km,QPRES) + hdt*srcQ(i,j,k3d,QPRES)
                qmo(i+1,j,km,QGAME) = qm(i+1,j,km,QGAME)
             endif

             qmo(i+1,j,km,QPRES) = max(qmo(i+1,j,km,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qd_lo, qd_hi, i+1, j, km)

#ifdef RADIATION
             qmo(i+1,j,km,qrad:qradhi) = ernewl(:)
             qmo(i+1,j,km,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i+1,j,km,QPRES)
             qmo(i+1,j,km,qreitot) = sum(qmo(i+1,j,km,qrad:qradhi)) + qmo(i+1,j,km,QREINT)
#endif

          endif

       enddo
    enddo

  end subroutine transyz

end module transverse_module
