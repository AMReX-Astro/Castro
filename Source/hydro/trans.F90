module transverse_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  !===========================================================================
  ! transx
  !===========================================================================
  subroutine transx(qym, qymo, qyp, qypo, &
                    qzm, qzmo, qzp, qzpo, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    fx, &
#ifdef RADIATION
                    rfx, &
#endif
                    fx_lo, fx_hi, &
                    qx, qx_lo, qx_hi, &
                    hdt, cdtdx, lo, hi)

  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: lo(3), hi(3)

#ifdef RADIATION
    real(rt) :: rfx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: cdtdx

    real(rt), intent(in) :: qym(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qyp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qzm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qzp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

    real(rt), intent(out) :: qymo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qypo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qzmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qzpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)


    integer i, j, k, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         rrry, rrly, rrrz, rrlz
    real(rt)         rury, ruly, rurz, rulz
    real(rt)         rvry, rvly, rvrz, rvlz
    real(rt)         rwry, rwly, rwrz, rwlz
    real(rt)         ekenry, ekenly, ekenrz, ekenlz
    real(rt)         rery, rely, rerz, relz
    real(rt)         rrnewry, rrnewly, rrnewrz, rrnewlz
    real(rt)         runewry, runewly, runewrz, runewlz
    real(rt)         rvnewry, rvnewly, rvnewrz, rvnewlz
    real(rt)         rwnewry, rwnewly, rwnewrz, rwnewlz
    real(rt)         renewry, renewly, renewrz, renewlz
    real(rt)         pnewry, pnewly, pnewrz, pnewlz
    real(rt)         rhoekenry, rhoekenly, rhoekenrz, rhoekenlz
    real(rt)         pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         compu
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !==========================================================================
    ! work on qy*
    !==========================================================================

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transverse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3)-dg(3), hi(3)+dg(3)
          do j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                if (j >= lo(2)) then
#if AMREX_SPACEDIM == 2
                   rrnew = qyp(i,j,k,QRHO) - hdt*(area1(i+1,j,k)*fx(i+1,j,k,URHO) - &
                                                  area1(i,j,k)*fx(i,j,k,URHO))/vol(i,j,k)
                   compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nqp) - &
                        hdt*(area1(i+1,j,k)*fx(i+1,j,k,n) - &
                             area1(i,j,k)*fx(i,j,k,n))/vol(i,j,k)
                   qypo(i,j,k,nqp) = compu/rrnew + hdt*srcQ(i,j,k,nqp)
#else
                   rrnew = qyp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                   compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                   qypo(i,j,k,nqp) = compu/rrnew
#endif
                end if

                if (j <= hi(2)) then
#if AMREX_SPACEDIM == 2
                   rrnew = qym(i,j+1,k,QRHO) - hdt*(area1(i+1,j,k)*fx(i+1,j,k,URHO) - &
                                                    area1(i,j,k)*fx(i,j,k,URHO))/vol(i,j,k)
                   compu = qym(i,j+1,k,QRHO)*qym(i,j+1,k,nqp) - &
                        hdt*(area1(i+1,j,k)*fx(i+1,j,k,n) - &
                             area1(i,j,k)*fx(i,j,k,n))/vol(i,j,k)
                   qymo(i,j+1,k,nqp) = compu/rrnew + hdt*srcQ(i,j,nqp)
#else
                   rrnew = qym(i,j+1,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                   compu = qym(i,j+1,k,QRHO)*qym(i,j+1,k,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                   qymo(i,j+1,k,nqp) = compu/rrnew
#endif
                end if
             end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-1, hi(2)+1
          do i = lo(1), hi(1)

             pgp  = qx(i+1,j,k,GDPRES)
             pgm  = qx(i  ,j,k,GDPRES)
             ugp  = qx(i+1,j,k,GDU   )
             ugm  = qx(i  ,j,k,GDU   )
             gegp = qx(i+1,j,k,GDGAME)
             gegm = qx(i  ,j,k,GDGAME)

#ifdef RADIATION
             lambda = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

############
             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
#endif

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

             !----------------------------------------------------------------
             ! qypo state
             !----------------------------------------------------------------

             if (j >= lo(2)) then

                ! Convert to conservation form
                rrry = qyp(i,j,k,QRHO)
                rury = rrry*qyp(i,j,k,QU)
                rvry = rrry*qyp(i,j,k,QV)
                rwry = rrry*qyp(i,j,k,QW)
                ekenry = HALF*rrry*sum(qyp(i,j,k,QU:QW)**2)
                rery = qyp(i,j,k,QREINT) + ekenry
#ifdef RADIATION
                err  = qyp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewry = rrry - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                runewry = rury - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
                rvnewry = rvry - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
                rwnewry = rwry - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
                renewry = rery - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
                runewry = runewry + dmom
                renewry = renewry + dre
                ernewr  = err(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qypo(i,j,k,QRHO) = rrnewry
                rhoinv = ONE/rrnewry
                qypo(i,j,k,QU) = runewry*rhoinv
                qypo(i,j,k,QV) = rvnewry*rhoinv
                qypo(i,j,k,QW) = rwnewry*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
                qypo(i,j,k,QREINT) = renewry - rhoekenry

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qypo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating (rho e).
                      qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                           cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewry = qyp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                      qypo(i,j,k,QPRES) = max(pnewry, small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qypo(i,j,k,QPRES) = qypo(i,j,k,QREINT)*(qypo(i,j,k,QGAME)-ONE)
                      qypo(i,j,k,QPRES) = max(qypo(i,j,k,QPRES),small_pres)
                   end if
                else
                   qypo(i,j,k,QPRES) = qyp(i,j,k,QPRES)
                   qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qypo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qypo(i,j,k,qrad:qradhi) = ernewr(:)
                qypo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,k,QPRES)
                qypo(i,j,k,qreitot) = sum(qypo(i,j,k,qrad:qradhi)) + qypo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             if (j <= hi(2)) then

                ! Convert to conservation form
                rrly = qym(i,j+1,k,QRHO)
                ruly = rrly*qym(i,j+1,k,QU)
                rvly = rrly*qym(i,j+1,k,QV)
                rwly = rrly*qym(i,j+1,k,QW)
                ekenly = HALF*rrly*sum(qym(i,j+1,k,QU:QW)**2)
                rely = qym(i,j+1,k,QREINT) + ekenly
#ifdef RADIATION
                erl  = qym(i,j+1,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewly = rrly - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                runewly = ruly - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
                rvnewly = rvly - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
                rwnewly = rwly - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
                renewly = rely - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
                runewly = runewly + dmom
                renewly = renewly + dre
                ernewl  = erl(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qymo(i,j+1,k,QRHO) = rrnewly
                rhoinv = ONE/rrnewly
                qymo(i,j+1,k,QU) = runewly*rhoinv
                qymo(i,j+1,k,QV) = rvnewly*rhoinv
                qymo(i,j+1,k,QW) = rwnewly*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
                qymo(i,j+1,k,QREINT) = renewly - rhoekenly

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qymo(i,j+1,k,QREINT) = qym(i,j+1,k,QREINT) - &
                           cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewly = qym(i,j+1,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                      qymo(i,j+1,k,QPRES) = max(pnewly,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qymo(i,j+1,k,QPRES) = qymo(i,j+1,k,QREINT)*(qymo(i,j+1,k,QGAME)-ONE)
                      qymo(i,j+1,k,QPRES) = max(qymo(i,j+1,k,QPRES), small_pres)
                   end if
                else
                   qymo(i,j+1,k,QPRES) = qym(i,j+1,k,QPRES)
                   qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME)
                endif

                call reset_edge_state_thermo(qymo, q_lo, q_hi, i, j+1, k)

#ifdef RADIATION
                qymo(i,j+1,k,qrad:qradhi) = ernewl(:)
                qymo(i,j+1,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j+1,k,QPRES)
                qymo(i,j+1,k,qreitot) = sum(qymo(i,j+1,k,qrad:qradhi)) + qymo(i,j+1,k,QREINT)
#endif

             end if

          end do
       end do
    end do

    !==========================================================================
    ! work on qz*
    !==========================================================================

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3)-1, hi(3)+1
          do j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                if (k >= lo(3)) then
                   rrnew = qzp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                   compu = qzp(i,j,k,QRHO)*qzp(i,j,k,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                   qzpo(i,j,k,nqp) = compu/rrnew
                endif

                if (k <= hi(3)) then
                   rrnew = qzm(i,j,k+1,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                   compu = qzm(i,j,k+1,QRHO)*qzm(i,j,k+1,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                   qzmo(i,j,k+1,nqp) = compu/rrnew
                endif
             end do
          end do
       end do

    end do


    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to z-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1), hi(1)

             pgp  = qx(i+1,j,k,GDPRES)
             pgm  = qx(i  ,j,k,GDPRES)
             ugp  = qx(i+1,j,k,GDU   )
             ugm  = qx(i  ,j,k,GDU   )
             gegp = qx(i+1,j,k,GDGAME)
             gegm = qx(i  ,j,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
#endif

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
             ! qzpo state
             !-------------------------------------------------------------------

             if (k >= lo(3)) then

                ! Convert to conservation form
                rrrz = qzp(i,j,k,QRHO)
                rurz = rrrz*qzp(i,j,k,QU)
                rvrz = rrrz*qzp(i,j,k,QV)
                rwrz = rrrz*qzp(i,j,k,QW)
                ekenrz = HALF*rrrz*sum(qzp(i,j,k,QU:QW)**2)
                rerz = qzp(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
                err  = qzp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewrz = rrrz - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                runewrz = rurz - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
                rvnewrz = rvrz - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
                rwnewrz = rwrz - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
                renewrz = rerz - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
                runewrz = runewrz + dmom
                renewrz = renewrz + dre
                ernewr  = err(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) + der(:)
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
                qzpo(i,j,k,QRHO) = rrnewrz
                rhoinv = ONE/rrnewrz
                qzpo(i,j,k,QU) = runewrz*rhoinv
                qzpo(i,j,k,QV) = rvnewrz*rhoinv
                qzpo(i,j,k,QW) = rwnewrz*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
                qzpo(i,j,k,QREINT) = renewrz - rhoekenrz

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qzpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating
                      ! (rho e).
                      qzpo(i,j,k,QREINT) = qzp(i,j,k,QREINT) - &
                           cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   end if

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewrz = qzp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                      qzpo(i,j,k,QPRES) = max(pnewrz,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qzpo(i,j,k,QPRES) = qzpo(i,j,k,QREINT)*(qzpo(i,j,k,QGAME)-ONE)
                      qzpo(i,j,k,QPRES) = max(qzpo(i,j,k,QPRES), small_pres)
                   endif
                else
                   qzpo(i,j,k,QPRES) = qzp(i,j,k,QPRES)
                   qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qzpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qzpo(i,j,k,qrad:qradhi) = ernewr(:)
                qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
                qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qradhi)) + qzpo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             if (k <= hi(3)) then

                ! Convert to conservation form
                rrlz = qzm(i,j,k+1,QRHO)
                rulz = rrlz*qzm(i,j,k+1,QU)
                rvlz = rrlz*qzm(i,j,k+1,QV)
                rwlz = rrlz*qzm(i,j,k+1,QW)
                ekenlz = HALF*rrlz*sum(qzm(i,j,k+1,QU:QW)**2)
                relz = qzm(i,j,k+1,QREINT) + ekenlz
#ifdef RADIATION
                erl  = qzm(i,j,k+1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewlz = rrlz - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                runewlz = rulz - cdtdx*(fx(i+1,j,k,UMX) - fx(i,j,k,UMX))
                rvnewlz = rvlz - cdtdx*(fx(i+1,j,k,UMY) - fx(i,j,k,UMY))
                rwnewlz = rwlz - cdtdx*(fx(i+1,j,k,UMZ) - fx(i,j,k,UMZ))
                renewlz = relz - cdtdx*(fx(i+1,j,k,UEDEN) - fx(i,j,k,UEDEN))
#ifdef RADIATION
                runewlz = runewlz + dmom
                renewlz = renewlz + dre
                ernewl  = erl(:) - cdtdx*(rfx(i+1,j,k,:) - rfx(i,j,k,:)) + der(:)
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
                qzmo(i,j,k+1,QRHO) = rrnewlz
                rhoinv = ONE/rrnewlz
                qzmo(i,j,k+1,QU) = runewlz*rhoinv
                qzmo(i,j,k+1,QV) = rvnewlz*rhoinv
                qzmo(i,j,k+1,QW) = rwnewlz*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
                qzmo(i,j,k+1,QREINT) = renewlz - rhoekenlz

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qzmo(i,j,k+1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qzmo(i,j,k+1,QREINT) = qzm(i,j,k+1,QREINT) - &
                           cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewlz = qzm(i,j,k+1,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                      qzmo(i,j,k+1,QPRES) = max(pnewlz,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qzmo(i,j,k+1,QGAME) = qzm(i,j,k+1,QGAME) + &
                           cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qzmo(i,j,k+1,QPRES) = qzmo(i,j,k+1,QREINT)*(qzmo(i,j,k+1,QGAME)-ONE)
                      qzmo(i,j,k+1,QPRES) = max(qzmo(i,j,k+1,QPRES), small_pres)
                   endif
                else
                   qzmo(i,j,k+1,QPRES) = qzm(i,j,k+1,QPRES)
                   qzmo(i,j,k+1,QGAME) = qzm(i,j,k+1,QGAME)
                endif

                call reset_edge_state_thermo(qzmo, q_lo, q_hi, i, j, k+1)

#ifdef RADIATION
                qzmo(i,j,k+1,qrad:qradhi) = ernewl(:)
                qzmo(i,j,k+1,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k+1,QPRES)
                qzmo(i,j,k+1,qreitot) = sum(qzmo(i,j,k+1,qrad:qradhi)) + qzmo(i,j,k+1,QREINT)
#endif

             end if

          end do
       end do
    end do

  end subroutine transx

  !===========================================================================
  ! transy
  !===========================================================================
  subroutine transy(qxm, qxmo, qxp, qxpo, &
                    qzm, qzmo, qzp, qzpo, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    fy, &
#ifdef RADIATION
                    rfy, &
#endif
                    fy_lo, fy_hi, &
                    qy, qy_lo, qy_hi, &
                    cdtdy, lo, hi)


  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer, intent(in) :: q_lo(3),q_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fy_lo(3),fy_hi(3)
    integer, intent(in) :: qy_lo(3),qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

#ifdef RADIATION
    real(rt), intent(in) :: rfy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: cdtdy

    real(rt), intent(in) :: qxm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qxp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qzm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qzp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    real(rt), intent(out) :: qxmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qxpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qzmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qzpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         rrrx, rrlx, rrrz, rrlz
    real(rt)         rurx, rulx, rurz, rulz
    real(rt)         rvrx, rvlx, rvrz, rvlz
    real(rt)         rwrx, rwlx, rwrz, rwlz
    real(rt)         ekenrx, ekenlx, ekenrz, ekenlz
    real(rt)         rerx, relx, rerz, relz
    real(rt)         rrnewrx, rrnewlx, rrnewrz, rrnewlz
    real(rt)         runewrx, runewlx, runewrz, runewlz
    real(rt)         rvnewrx, rvnewlx, rvnewrz, rvnewlz
    real(rt)         rwnewrx, rwnewlx, rwnewrz, rwnewlz
    real(rt)         renewrx, renewlx, renewrz, renewlz
    real(rt)         pnewrx, pnewlx, pnewrz, pnewlz
    real(rt)         rhoekenrx, rhoekenlx, rhoekenrz, rhoekenlz
    real(rt)         pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         compu
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1, ugc
    integer :: g
#endif

    logical :: reset_state

    !=========================================================================
    ! work on qx*
    !=========================================================================

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3)-1, hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                if (i >= lo(1)) then
                   rrnew = qxp(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                   compu = qxp(i,j,k,QRHO)*qxp(i,j,k,nqp) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                   qxpo(i,j,k,nqp) = compu/rrnew
                endif

                if (i <= hi(1)) then
                   rrnew = qxm(i+1,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                   compu = qxm(i+1,j,k,QRHO)*qxm(i+1,j,k,nqp) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                   qxmo(i+1,j,k,nqp) = compu/rrnew
                end if
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2), hi(2)
          do i = lo(1)-1, hi(1)+1

             pgp  = qy(i,j+1,k,GDPRES)
             pgm  = qy(i,j  ,k,GDPRES)
             ugp  = qy(i,j+1,k,GDV   )
             ugm  = qy(i,j  ,k,GDV   )
             gegp = qy(i,j+1,k,GDGAME)
             gegm = qy(i,j  ,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
#endif

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

             !-------------------------------------------------------------
             ! qxpo state
             !-------------------------------------------------------------

             if (i >= lo(1)) then

                ! Convert to conservation form
                rrrx = qxp(i,j,k,QRHO)
                rurx = rrrx*qxp(i,j,k,QU)
                rvrx = rrrx*qxp(i,j,k,QV)
                rwrx = rrrx*qxp(i,j,k,QW)
                ekenrx = HALF*rrrx*sum(qxp(i,j,k,QU:QW)**2)
                rerx = qxp(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
                err  = qxp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewrx = rrrx - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                runewrx = rurx - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
                rvnewrx = rvrx - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
                rwnewrx = rwrx - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
                renewrx = rerx - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewrx = rvnewrx + dmom
                renewrx = renewrx + dre
                ernewr = err(:) - cdtdy*(rfy(i,j+1,k,:) - rfy(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qxpo(i,j,k,QRHO) = rrnewrx
                rhoinv = ONE/rrnewrx
                qxpo(i,j,k,QU) = runewrx*rhoinv
                qxpo(i,j,k,QV) = rvnewrx*rhoinv
                qxpo(i,j,k,QW) = rwnewrx*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
                qxpo(i,j,k,QREINT) = renewrx - rhoekenrx

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qxpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating (rho e).
                      qxpo(i,j,k,QREINT) = qxp(i,j,k,QREINT) - &
                           cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewrx = qxp(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                      qxpo(i,j,k,QPRES) = max(pnewrx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME) + &
                           cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxpo(i,j,k,QPRES) = qxpo(i,j,k,QREINT)*(qxpo(i,j,k,QGAME)-ONE)
                      qxpo(i,j,k,QPRES) = max(qxpo(i,j,k,QPRES), small_pres)
                   endif
                else
                   qxpo(i,j,k,QPRES) = qxp(i,j,k,QPRES)
                   qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qxpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qxpo(i,j,k,qrad:qradhi) = ernewr(:)
                qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,k,QPRES)
                qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qradhi)) + qxpo(i,j,k,QREINT)
#endif

             endif

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             if (i <= hi(1)) then

                ! Convert to conservation form
                rrlx = qxm(i+1,j,k,QRHO)
                rulx = rrlx*qxm(i+1,j,k,QU)
                rvlx = rrlx*qxm(i+1,j,k,QV)
                rwlx = rrlx*qxm(i+1,j,k,QW)
                ekenlx = HALF*rrlx*sum(qxm(i+1,j,k,QU:QW)**2)
                relx = qxm(i+1,j,k,QREINT) + ekenlx
#ifdef RADIATION
                erl  = qxm(i+1,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewlx = rrlx - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                runewlx = rulx - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
                rvnewlx = rvlx - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
                rwnewlx = rwlx - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
                renewlx = relx - cdtdy*(fy(i,j+1,k,UEDEN)- fy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewlx = rvnewlx + dmom
                renewlx = renewlx + dre
                ernewl  = erl(:) - cdtdy*(rfy(i,j+1,k,:) - rfy(i,j,k,:)) + der(:)
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

                qxmo(i+1,j,k,QRHO) = rrnewlx
                rhoinv = ONE/rrnewlx
                qxmo(i+1,j,k,QU) = runewlx*rhoinv
                qxmo(i+1,j,k,QV) = rvnewlx*rhoinv
                qxmo(i+1,j,k,QW) = rwnewlx*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
                qxmo(i+1,j,k,QREINT) = renewlx - rhoekenlx

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qxmo(i+1,j,k,QREINT) = qxm(i+1,j,k,QREINT) - &
                           cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewlx = qxm(i+1,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                      qxmo(i+1,j,k,QPRES) = max(pnewlx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME) + &
                           cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxmo(i+1,j,k,QPRES) = qxmo(i+1,j,k,QREINT)*(qxmo(i+1,j,k,QGAME)-ONE)
                      qxmo(i+1,j,k,QPRES) = max(qxmo(i+1,j,k,QPRES), small_pres)
                   endif
                else
                   qxmo(i+1,j,k,QPRES) = qxm(i+1,j,k,QPRES)
                   qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qxmo, q_lo, q_hi, i+1, j, k)

#ifdef RADIATION
                qxmo(i+1,j,k,qrad:qradhi) = ernewl(:)
                qxmo(i+1,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i+1,j,k,QPRES)
                qxmo(i+1,j,k,qreitot) = sum(qxmo(i+1,j,k,qrad:qradhi)) + qxmo(i+1,j,k,QREINT)
#endif

             end if

          end do
       end do
    end do


    !=========================================================================
    ! work on qz*
    !=========================================================================

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3)-1, hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1
                if (k >= lo(3)) then
                   rrnew = qzp(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                   compu = qzp(i,j,k,QRHO)*qzp(i,j,k,nqp) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                   qzpo(i,j,k,nqp) = compu/rrnew
                end if

                if (k <= hi(3)) then
                   rrnew = qzm(i,j,k+1,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                   compu = qzm(i,j,k+1,QRHO)*qzm(i,j,k+1,nqp) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                   qzmo(i,j,k+1,nqp) = compu/rrnew
                end if
             end do
          end do
       end do
    end do


    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to z-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3)-1, hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1)-1, hi(1)+1

             pgp  = qy(i,j+1,k,GDPRES)
             pgm  = qy(i,j  ,k,GDPRES)
             ugp  = qy(i,j+1,k,GDV   )
             ugm  = qy(i,j  ,k,GDV   )
             gegp = qy(i,j+1,k,GDGAME)
             gegm = qy(i,j  ,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
#endif


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

             !-------------------------------------------------------------------
             ! qzpo states
             !-------------------------------------------------------------------

             if (k >= lo(3)) then

                ! Convert to conservation form
                rrrz = qzp(i,j,k,QRHO)
                rurz = rrrz*qzp(i,j,k,QU)
                rvrz = rrrz*qzp(i,j,k,QV)
                rwrz = rrrz*qzp(i,j,k,QW)
                ekenrz = HALF*rrrz*sum(qzp(i,j,k,QU:QW)**2)
                rerz = qzp(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
                err  = qzp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewrz = rrrz - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                runewrz = rurz - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
                rvnewrz = rvrz - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
                rwnewrz = rwrz - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
                renewrz = rerz - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewrz = rvnewrz + dmom
                renewrz = renewrz + dre
                ernewr  = err(:) - cdtdy*(rfy(i,j+1,k,:) - rfy(i,j,k,:)) + der(:)
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
                qzpo(i,j,k,QRHO) = rrnewrz
                rhoinv = ONE/rrnewrz
                qzpo(i,j,k,QU) = runewrz*rhoinv
                qzpo(i,j,k,QV) = rvnewrz*rhoinv
                qzpo(i,j,k,QW) = rwnewrz*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
                qzpo(i,j,k,QREINT) = renewrz - rhoekenrz

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qzpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qzpo(i,j,k,QREINT) = qzp(i,j,k,QREINT) - &
                           cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewrz = qzp(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                      qzpo(i,j,k,QPRES) = max(pnewrz,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME) + &
                           cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qzpo(i,j,k,QPRES) = qzpo(i,j,k,QREINT)*(qzpo(i,j,k,QGAME)-ONE)
                      qzpo(i,j,k,QPRES) = max(qzpo(i,j,k,QPRES), small_pres)
                   endif
                else
                   qzpo(i,j,k,QPRES) = qzp(i,j,k,QPRES)
                   qzpo(i,j,k,QGAME) = qzp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qzpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qzpo(i,j,k,qrad:qradhi) = ernewr(:)
                qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
                qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qradhi)) + qzpo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             if (k <= hi(3)) then

                ! Convert to conservation form
                rrlz = qzm(i,j,k+1,QRHO)
                rulz = rrlz*qzm(i,j,k+1,QU)
                rvlz = rrlz*qzm(i,j,k+1,QV)
                rwlz = rrlz*qzm(i,j,k+1,QW)
                ekenlz = HALF*rrlz*sum(qzm(i,j,k+1,QU:QW)**2)
                relz = qzm(i,j,k+1,QREINT) + ekenlz
#ifdef RADIATION
                erl  = qzm(i,j,k+1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewlz = rrlz - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                runewlz = rulz - cdtdy*(fy(i,j+1,k,UMX) - fy(i,j,k,UMX))
                rvnewlz = rvlz - cdtdy*(fy(i,j+1,k,UMY) - fy(i,j,k,UMY))
                rwnewlz = rwlz - cdtdy*(fy(i,j+1,k,UMZ) - fy(i,j,k,UMZ))
                renewlz = relz - cdtdy*(fy(i,j+1,k,UEDEN) - fy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewlz = rvnewlz + dmom
                renewlz = renewlz + dre
                ernewl  = erl(:) - cdtdy*(rfy(i,j+1,k,:) - rfy(i,j,k,:)) + der
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
                qzmo(i,j,k+1,QRHO) = rrnewlz
                rhoinv = ONE/rrnewlz
                qzmo(i,j,k+1,QU) = runewlz*rhoinv
                qzmo(i,j,k+1,QV) = rvnewlz*rhoinv
                qzmo(i,j,k+1,QW) = rwnewlz*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
                qzmo(i,j,k+1,QREINT) = renewlz - rhoekenlz

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qzmo(i,j,k+1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qzmo(i,j,k+1,QREINT) = qzm(i,j,k+1,QREINT) - &
                           cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewlz = qzm(i,j,k+1,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                      qzmo(i,j,k+1,QPRES) = max(pnewlz,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qzmo(i,j,k+1,QGAME) = qzm(i,j,k+1,QGAME) + &
                           cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qzmo(i,j,k+1,QPRES) = qzmo(i,j,k+1,QREINT)*(qzmo(i,j,k+1,QGAME)-ONE)
                      qzmo(i,j,k+1,QPRES) = max(qzmo(i,j,k+1,QPRES), small_pres)
                   endif
                else
                   qzmo(i,j,k+1,QPRES) = qzm(i,j,k+1,QPRES)
                   qzmo(i,j,k+1,QGAME) = qzm(i,j,k+1,QGAME)
                endif

                call reset_edge_state_thermo(qzmo, q_lo, q_hi, i, j, k+1)

#ifdef RADIATION
                qzmo(i,j,k+1,qrad:qradhi) = ernewl(:)
                qzmo(i,j,k+1,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k+1,QPRES)
                qzmo(i,j,k+1,qreitot) = sum(qzmo(i,j,k+1,qrad:qradhi)) + qzmo(i,j,k+1,QREINT)
#endif

             end if

          end do
       end do
    end do

  end subroutine transy

  !===========================================================================
  ! transz
  !===========================================================================
  subroutine transz(qxm, qxmo, qxp, qxpo, &
                    qym, qymo, qyp, qypo, q_lo, q_hi, &
                    qaux, qa_lo, qa_hi, &
                    fz, &
#ifdef RADIATION
                    rfz, &
#endif
                    fz_lo, fz_hi, &
                    qz, qz_lo, qz_hi, &
                    cdtdz, lo, hi)

    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, ppm_type, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    integer, intent(in) :: qz_lo(3), qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

#ifdef RADIATION
    real(rt), intent(in) :: rfz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: cdtdz

    real(rt), intent(in) :: qxm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qxp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qym(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qyp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    real(rt), intent(out) :: qxmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qxpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qymo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qypo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew, rr
    real(rt)         rrrx, rrry, rrlx, rrly
    real(rt)         rurx, rury, rulx, ruly
    real(rt)         rvrx, rvry, rvlx, rvly
    real(rt)         rwrx, rwry, rwlx, rwly
    real(rt)         ekenrx, ekenry, ekenlx, ekenly
    real(rt)         rerx, rery, relx, rely
    real(rt)         rrnewrx, rrnewry, rrnewlx, rrnewly
    real(rt)         runewrx, runewry, runewlx, runewly
    real(rt)         rvnewrx, rvnewry, rvnewlx, rvnewly
    real(rt)         rwnewrx, rwnewry, rwnewlx, rwnewly
    real(rt)         renewrx, renewry, renewlx, renewly
    real(rt)         pnewrx, pnewry, pnewlx, pnewly
    real(rt)         rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
    real(rt)         pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         compu
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, luge, lamge, &
         ergp, errx, ernewrx, erry, ernewry, ergm, erlx, ernewlx, erly, ernewly
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !=========================================================================
    ! work on qx*
    !=========================================================================

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transverse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

        do k = lo(3), hi(3)
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                if (i >= lo(1)) then
                   rrnew = qxp(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                   compu = qxp(i,j,k,QRHO)*qxp(i,j,k,nqp) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                   qxpo(i,j,k,nqp) = compu/rrnew
                end if

                if (i <= hi(1)) then
                   rrnew = qxm(i+1,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                   compu = qxm(i+1,j,k,QRHO)*qxm(i+1,j,k,nqp) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                   qxmo(i+1,j,k,nqp) = compu/rrnew
                end if
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the z-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             pgp  = qz(i,j,k+1,GDPRES)
             pgm  = qz(i,j,k,GDPRES)
             ugp  = qz(i,j,k+1,GDW   )
             ugm  = qz(i,j,k,GDW   )
             gegp = qz(i,j,k+1,GDGAME)
             gegm = qz(i,j,k,GDGAME)

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ergp = qz(i,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergm = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
#endif

#ifdef RADIATION
             lamge = lambda(:) * (ergp(:)-ergm(:))
             dmz = - cdtdz*sum(lamge)
             luge = HALF*(ugp+ugm) * lamge(:)
             dre = -cdtdz*sum(luge)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdz*HALF*(ugp+ugm) * f1 * (ergp(g) - ergm(g))
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdz * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
                end do
             else ! mixed frame
                der(:) = cdtdz * luge
             end if
#endif

             !-------------------------------------------------------------------
             ! qxpo state
             !-------------------------------------------------------------------

             if (i >= lo(1)) then

                ! Convert to conservation form
                rrrx = qxp(i,j,k,QRHO)
                rurx = rrrx*qxp(i,j,k,QU)
                rvrx = rrrx*qxp(i,j,k,QV)
                rwrx = rrrx*qxp(i,j,k,QW)
                ekenrx = HALF*rrrx*sum(qxp(i,j,k,QU:QW)**2)
                rerx = qxp(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
                errx = qxp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewrx = rrrx - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                runewrx = rurx - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
                rvnewrx = rvrx - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
                rwnewrx = rwrx - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
                renewrx = rerx - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))
#ifdef RADIATION
                rwnewrx = rwnewrx + dmz
                renewrx = renewrx + dre
                ernewrx = errx(:) - cdtdz*(rfz(i,j,k+1,:) - rfz(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qxpo(i,j,k,QRHO) = rrnewrx
                rhoinv = ONE/rrnewrx
                qxpo(i,j,k,QU) = runewrx*rhoinv
                qxpo(i,j,k,QV) = rvnewrx*rhoinv
                qxpo(i,j,k,QW) = rwnewrx*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
                qxpo(i,j,k,QREINT) = renewrx - rhoekenrx

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qxpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating (rho e).
                      qxpo(i,j,k,QREINT) = qxp(i,j,k,QREINT) - &
                           cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewrx = qxp(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qxpo(i,j,k,QPRES) = max(pnewrx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxpo(i,j,k,QPRES) = qxpo(i,j,k,QREINT)*(qxpo(i,j,k,QGAME)-ONE)
                      qxpo(i,j,k,QPRES) = max(qxpo(i,j,k,QPRES), small_pres)
                   endif
                else
                   qxpo(i,j,k,QPRES) = qxp(i,j,k,QPRES)
                   qxpo(i,j,k,QGAME) = qxp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qxpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qxpo(i,j,k,qrad:qradhi) = ernewrx(:)
                qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,k,QPRES)
                qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qradhi)) + qxpo(i,j,k,QREINT)
#endif
             end if

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             if (i <= hi(1)) then

                ! Convert to conservation form
                rrlx = qxm(i+1,j,k,QRHO)
                rulx = rrlx*qxm(i+1,j,k,QU)
                rvlx = rrlx*qxm(i+1,j,k,QV)
                rwlx = rrlx*qxm(i+1,j,k,QW)
                ekenlx = HALF*rrlx*sum(qxm(i+1,j,k,QU:QW)**2)
                relx = qxm(i+1,j,k,QREINT) + ekenlx
#ifdef RADIATION
                erlx = qxm(i+1,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewlx = rrlx - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                runewlx = rulx - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
                rvnewlx = rvlx - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
                rwnewlx = rwlx - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
                renewlx = relx - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))
#ifdef RADIATION
                rwnewlx = rwnewlx + dmz
                renewlx = renewlx + dre
                ernewlx = erlx(:) - cdtdz*(rfz(i,j,k+1,:) - rfz(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qxmo(i+1,j,k,QRHO) = rrnewlx
                rhoinv = ONE/rrnewlx
                qxmo(i+1,j,k,QU) = runewlx*rhoinv
                qxmo(i+1,j,k,QV) = rvnewlx*rhoinv
                qxmo(i+1,j,k,QW) = rwnewlx*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
                qxmo(i+1,j,k,QREINT) = renewlx - rhoekenlx

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qxmo(i+1,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qxmo(i+1,j,k,QREINT) = qxm(i+1,j,k,QREINT) - &
                           cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewlx = qxm(i+1,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qxmo(i+1,j,k,QPRES) = max(pnewlx,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qxmo(i+1,j,k,QPRES) = qxmo(i+1,j,k,QREINT)*(qxmo(i+1,j,k,QGAME)-ONE)
                      qxmo(i+1,j,k,QPRES) = max(qxmo(i+1,j,k,QPRES), small_pres)
                   end if
                else
                   qxmo(i+1,j,k,QPRES) = qxm(i+1,j,k,QPRES)
                   qxmo(i+1,j,k,QGAME) = qxm(i+1,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qxmo, q_lo, q_hi, i+1, j, k)

#ifdef RADIATION
                qxmo(i+1,j,k,qrad:qradhi) = ernewlx(:)
                qxmo(i+1,j,k,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i+1,j,k,QPRES)
                qxmo(i+1,j,k,qreitot) = sum(qxmo(i+1,j,k,qrad:qradhi)) + qxmo(i+1,j,k,QREINT)
#endif

             end if

          end do
       end do
    end do

    !=========================================================================
    ! work on qy*
    !=========================================================================

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                if (j >= lo(2)) then
                   rrnew = qyp(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                   compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nqp) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                   qypo(i,j,k,nqp) = compu/rrnew
                end if

                if (j <= hi(2)) then
                   rrnew = qym(i,j+1,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                   compu = qym(i,j+1,k,QRHO)*qym(i,j+1,k,nqp) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                   qymo(i,j+1,k,nqp) = compu/rrnew
                end if
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the z-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             pgp  = qz(i,j,k+1,GDPRES)
             pgm  = qz(i,j,k,GDPRES)
             ugp  = qz(i,j,k+1,GDW   )
             ugm  = qz(i,j,k,GDW   )
             gegp = qz(i,j,k+1,GDGAME)
             gegm = qz(i,j,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ergp = qz(i,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergm = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

             dup = pgp*ugp - pgm*ugm
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             du = ugp-ugm
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j,k,QGAMCG)
#else
             gamc = qaux(i,j,k,QGAMC)
#endif

#ifdef RADIATION
             lamge = lambda(:) * (ergp(:)-ergm(:))
             dmz = - cdtdz*sum(lamge)
             luge = HALF*(ugp+ugm) * lamge(:)
             dre = -cdtdz*sum(luge)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdz*HALF*(ugp+ugm) * f1 * (ergp(g) - ergm(g))
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdz* f1 * HALF*(ergp(g)+ergm(g))*(ugm-ugp)
                end do
             else ! mixed frame
                der(:) = cdtdz * luge
             end if
#endif

             !-------------------------------------------------------------------
             ! qypo state
             !-------------------------------------------------------------------

             if (j >= lo(2)) then

                ! Convert to conservation form
                rrry = qyp(i,j,k,QRHO)
                rury = rrry*qyp(i,j,k,QU)
                rvry = rrry*qyp(i,j,k,QV)
                rwry = rrry*qyp(i,j,k,QW)
                ekenry = HALF*rrry*sum(qyp(i,j,k,QU:QW)**2)
                rery = qyp(i,j,k,QREINT) + ekenry
#ifdef RADIATION
                erry = qyp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewry = rrry - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                runewry = rury - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
                rvnewry = rvry - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
                rwnewry = rwry - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
                renewry = rery - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))
#ifdef RADIATION
                rwnewry = rwnewry + dmz
                renewry = renewry + dre
                ernewry = erry(:) - cdtdz*(rfz(i,j,k+1,:) - rfz(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qypo(i,j,k,QRHO) = rrnewry
                rhoinv = ONE/rrnewry
                qypo(i,j,k,QU) = runewry*rhoinv
                qypo(i,j,k,QV) = rvnewry*rhoinv
                qypo(i,j,k,QW) = rwnewry*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
                qypo(i,j,k,QREINT) = renewry - rhoekenry

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qypo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                           cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewry = qyp(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qypo(i,j,k,QPRES) = max(pnewry,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qypo(i,j,k,QPRES) = qypo(i,j,k,QREINT)*(qypo(i,j,k,QGAME)-ONE)
                      qypo(i,j,k,QPRES) = max(qypo(i,j,k,QPRES), small_pres)
                   endif
                else
                   qypo(i,j,k,QPRES) = qyp(i,j,k,QPRES)
                   qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME)
                endif

                call reset_edge_state_thermo(qypo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qypo(i,j,k,qrad:qradhi) = ernewry(:)
                qypo(i,j,k,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,k,QPRES)
                qypo(i,j,k,qreitot) = sum(qypo(i,j,k,qrad:qradhi)) + qypo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             if (j <= hi(2)) then

                ! Convert to conservation form
                rrly = qym(i,j+1,k,QRHO)
                ruly = rrly*qym(i,j+1,k,QU)
                rvly = rrly*qym(i,j+1,k,QV)
                rwly = rrly*qym(i,j+1,k,QW)
                ekenly = HALF*rrly*sum(qym(i,j+1,k,QU:QW)**2)
                rely = qym(i,j+1,k,QREINT) + ekenly
#ifdef RADIATION
                erly = qym(i,j+1,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewly = rrly - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                runewly = ruly - cdtdz*(fz(i,j,k+1,UMX) - fz(i,j,k,UMX))
                rvnewly = rvly - cdtdz*(fz(i,j,k+1,UMY) - fz(i,j,k,UMY))
                rwnewly = rwly - cdtdz*(fz(i,j,k+1,UMZ) - fz(i,j,k,UMZ))
                renewly = rely - cdtdz*(fz(i,j,k+1,UEDEN) - fz(i,j,k,UEDEN))
#ifdef RADIATION
                rwnewly = rwnewly + dmz
                renewly = renewly + dre
                ernewly = erly(:) - cdtdz*(rfz(i,j,k+1,:) - rfz(i,j,k,:)) + der(:)
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

                ! Convert back to primitive form
                qymo(i,j+1,k,QRHO) = rrnewly
                rhoinv = ONE/rrnewly
                qymo(i,j+1,k,QU) = runewly*rhoinv
                qymo(i,j+1,k,QV) = rvnewly*rhoinv
                qymo(i,j+1,k,QW) = rwnewly*rhoinv

                ! note: we run the risk of (rho e) being negative here
                rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
                qymo(i,j+1,k,QREINT) = renewly - rhoekenly

                if (.not. reset_state) then
                   ! do the transverse terms for p, gamma, and rhoe, as necessary

                   if (transverse_reset_rhoe == 1 .and. qymo(i,j+1,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qymo(i,j+1,k,QREINT) = qym(i,j+1,k,QREINT) - &
                           cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewly = qym(i,j+1,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                      qymo(i,j+1,k,QPRES) = max(pnewly,small_pres)
                   else
                      ! Update gammae with its transverse terms
                      qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME) + &
                           cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                      ! and compute the p edge state from this and (rho e)
                      qymo(i,j+1,k,QPRES) = qymo(i,j+1,k,QREINT)*(qymo(i,j+1,k,QGAME)-ONE)
                      qymo(i,j+1,k,QPRES) = max(qymo(i,j+1,k,QPRES), small_pres)
                   endif
                else
                   qymo(i,j+1,k,QPRES) = qym(i,j+1,k,QPRES)
                   qymo(i,j+1,k,QGAME) = qym(i,j+1,k,QGAME)
                endif

                call reset_edge_state_thermo(qymo, q_lo, q_hi, i, j+1, k)

#ifdef RADIATION
                qymo(i,j+1,k,qrad:qradhi) = ernewly(:)
                qymo(i,j+1,k,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j+1,k,QPRES)
                qymo(i,j+1,k,qreitot) = sum(qymo(i,j+1,k,qrad:qradhi)) + qymo(i,j+1,k,QREINT)
#endif
             end if

          end do
       end do
    end do

  end subroutine transz

  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz(qm, qmo, qp, qpo, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi, &
                     fyz, &
#ifdef RADIATION
                     rfyz, &
#endif
                     fy_lo, fy_hi, &
                     fzy, &
#ifdef RADIATION
                     rfzy, &
#endif
                     fz_lo, fz_hi, &
                     qy, qy_lo, qy_hi, &
                     qz, qz_lo, qz_hi, &
                     srcQ, src_lo, src_hi, &
                     hdt, cdtdy, cdtdz, lo, hi)


    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, ppm_type, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer, intent(in) :: q_lo(3),q_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fy_lo(3),fy_hi(3)
    integer, intent(in) :: fz_lo(3),fz_hi(3)
    integer, intent(in) :: qy_lo(3),qy_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: hdt, cdtdy, cdtdz

#ifdef RADIATION
    real(rt), intent(in) :: rfyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fyz(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: fzy(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt), intent(in) :: srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt)         rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt)         rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt)         rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt)         pnewr, pnewl
    real(rt)         pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt)         pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt)         uyav, geyav, dgey, uzav, gezav, dgez
    real(rt)         compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt)         :: dmy, dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1)-1, hi(1)+1

                if (i >= lo(1)) then
                   rrr = qp(i,j,k,QRHO)
                   compr = rrr*qp(i,j,k,nqp)
                   rrnewr = rrr - cdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                                - cdtdz*(fzy(i,j  ,k+1,URHO) - fzy(i,j,k,URHO))
                   compnr = compr - cdtdy*(fyz(i,j+1,k,n   ) - fyz(i,j,k,n)) &
                                 - cdtdz*(fzy(i,j  ,k+1,n   ) - fzy(i,j,k,n))

                   qpo(i  ,j,k,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k,nqp)
                end if

                if (i <= hi(1)) then
                   rrl = qm(i+1,j,k,QRHO)
                   compl = rrl*qm(i+1,j,k,nqp)
                   rrnewl = rrl - cdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                                - cdtdz*(fzy(i,j  ,k+1,URHO) - fzy(i,j,k,URHO))
                   compnl = compl - cdtdy*(fyz(i,j+1,k,n   ) - fyz(i,j,k,n)) &
                                - cdtdz*(fzy(i,j  ,k+1,n   ) - fzy(i,j,k,n))

                   qmo(i+1,j,k,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k,nqp)
                end if
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse yz and zy differences to the x-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1)-1, hi(1)+1

             pgyp  = qy(i,j+1,k,GDPRES)
             pgym  = qy(i,j,k,GDPRES)
             ugyp  = qy(i,j+1,k,GDV)
             ugym  = qy(i,j,k,GDV)
             gegyp = qy(i,j+1,k,GDGAME)
             gegym = qy(i,j,k,GDGAME)
#ifdef RADIATION
             ergyp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgzp  = qz(i,j,k+1,GDPRES)
             pgzm  = qz(i,j,k,GDPRES)
             ugzp  = qz(i,j,k+1,GDW)
             ugzm  = qz(i,j,k,GDW)
             gegzp = qz(i,j,k+1,GDGAME)
             gegzm = qz(i,j,k,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMC))*duy - uyav*dgey )
#endif

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)

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

             if (i >= lo(1)) then
                ! Convert to conservation form
                rrr = qp(i,j,k,QRHO)
                rur = rrr*qp(i,j,k,QU)
                rvr = rrr*qp(i,j,k,QV)
                rwr = rrr*qp(i,j,k,QW)
                ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
                rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
                err = qp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewr = rrr - cdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                             - cdtdz*(fzy(i,j,k+1,URHO) - fzy(i,j,k,URHO))
                runewr = rur - cdtdy*(fyz(i,j+1,k,UMX) - fyz(i,j,k,UMX)) &
                             - cdtdz*(fzy(i,j,k+1,UMX) - fzy(i,j,k,UMX))
                rvnewr = rvr - cdtdy*(fyz(i,j+1,k,UMY) - fyz(i,j,k,UMY)) &
                             - cdtdz*(fzy(i,j,k+1,UMY) - fzy(i,j,k,UMY))
                rwnewr = rwr - cdtdy*(fyz(i,j+1,k,UMZ) - fyz(i,j,k,UMZ)) &
                             - cdtdz*(fzy(i,j,k+1,UMZ) - fzy(i,j,k,UMZ))
                renewr = rer - cdtdy*(fyz(i,j+1,k,UEDEN) - fyz(i,j,k,UEDEN)) &
                             - cdtdz*(fzy(i,j,k+1,UEDEN) - fzy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewr = rvnewr + dmy
                rwnewr = rwnewr + dmz
                renewr = renewr + dre
                ernewr = err(:) - cdtdy*(rfyz(i,j+1,k,:) - rfyz(i,j,k,:)) &
                                - cdtdz*(rfzy(i,j  ,k+1,:) - rfzy(i,j,k,:)) &
                                + der(:)
#endif

                ! Reset to original value if adding transverse terms
                ! made density negative
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

                qpo(i,j,k,QRHO  ) = rrnewr
                qpo(i,j,k,QU    ) = runewr/rrnewr
                qpo(i,j,k,QV    ) = rvnewr/rrnewr
                qpo(i,j,k,QW    ) = rwnewr/rrnewr

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qpo(i,j,k,QRHO  ) = qpo(i,j,k,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qpo(i,j,k,QU:QW) = qpo(i,j,k,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                qpo(i,j,k,QREINT) = renewr - rhoekenr
                if (ppm_type == 0) then
                   qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating
                      ! (rho e).
                      qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                           - cdtdy*(fyz(i,j+1,k,UEINT) - fyz(i,j,k,UEINT) + pyav*duy) &
                           - cdtdz*(fzy(i,j  ,k+1,UEINT) - fzy(i,j,k,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewr = qp(i,j,k,QPRES) - pynew - pznew
                      qpo(i,j,k,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + geynew + geznew

                      ! and compute the p edge state from this and (rho e)
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                   end if
                else
                   qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                   if (ppm_type == 0) then
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
                endif

                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

                call reset_edge_state_thermo(qpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qpo(i,j,k,qrad:qradhi) = ernewr(:)
                qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
                qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             if (i <= hi(1)) then

                ! Convert to conservation form
                rrl = qm(i+1,j,k,QRHO)
                rul = rrl*qm(i+1,j,k,QU)
                rvl = rrl*qm(i+1,j,k,QV)
                rwl = rrl*qm(i+1,j,k,QW)
                ekenl = HALF*rrl*sum(qm(i+1,j,k,QU:QW)**2)
                rel = qm(i+1,j,k,QREINT) + ekenl
#ifdef RADIATION
                erl = qm(i+1,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewl = rrl - cdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                             - cdtdz*(fzy(i,j,k+1,URHO) - fzy(i,j,k,URHO))
                runewl = rul - cdtdy*(fyz(i,j+1,k,UMX) - fyz(i,j,k,UMX)) &
                             - cdtdz*(fzy(i,j,k+1,UMX) - fzy(i,j,k,UMX))
                rvnewl = rvl - cdtdy*(fyz(i,j+1,k,UMY) - fyz(i,j,k,UMY)) &
                             - cdtdz*(fzy(i,j,k+1,UMY) - fzy(i,j,k,UMY))
                rwnewl = rwl - cdtdy*(fyz(i,j+1,k,UMZ) - fyz(i,j,k,UMZ)) &
                             - cdtdz*(fzy(i,j,k+1,UMZ) - fzy(i,j,k,UMZ))
                renewl = rel - cdtdy*(fyz(i,j+1,k,UEDEN) - fyz(i,j,k,UEDEN)) &
                             - cdtdz*(fzy(i,j,k+1,UEDEN) - fzy(i,j,k,UEDEN))
#ifdef RADIATION
                rvnewl = rvnewl + dmy
                rwnewl = rwnewl + dmz
                renewl = renewl + dre
                ernewl = erl(:) - cdtdy*(rfyz(i,j+1,k,:) - rfyz(i,j,k,:)) &
                                - cdtdz*(rfzy(i,j  ,k+1,:) - rfzy(i,j,k,:)) &
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

                qmo(i+1,j,k,QRHO   ) = rrnewl
                qmo(i+1,j,k,QU     ) = runewl/rrnewl
                qmo(i+1,j,k,QV     ) = rvnewl/rrnewl
                qmo(i+1,j,k,QW     ) = rwnewl/rrnewl

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qmo(i+1,j,k,QRHO   ) = qmo(i+1,j,k,QRHO   ) + hdt*srcQ(i,j,k,QRHO)
                   qmo(i+1,j,k,QU:QW) = qmo(i+1,j,k,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
                qmo(i+1,j,k,QREINT ) = renewl - rhoekenl
                if (ppm_type == 0) then
                   qmo(i+1,j,k,QREINT ) = qmo(i+1,j,k,QREINT ) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qmo(i+1,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qmo(i+1,j,k,QREINT ) = qm(i+1,j,k,QREINT) &
                           - cdtdy*(fyz(i,j+1,k,UEINT) - fyz(i,j,k,UEINT) + pyav*duy) &
                           - cdtdz*(fzy(i,j  ,k+1,UEINT) - fzy(i,j,k,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qmo(i+1,j,k,QREINT ) = qmo(i+1,j,k,QREINT ) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewl = qm(i+1,j,k,QPRES) - pynew - pznew
                      qmo(i+1,j,k,QPRES  ) = pnewl
                      if (ppm_type == 0) then
                         qmo(i+1,j,k,QPRES  ) = qmo(i+1,j,k,QPRES  ) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qmo(i+1,j,k,QGAME) = qm(i+1,j,k,QGAME) + geynew + geznew

                      ! and compute the p edge state from this and (rho e)
                      qmo(i+1,j,k,QPRES) = qmo(i+1,j,k,QREINT)*(qmo(i+1,j,k,QGAME)-ONE)
                   end if
                else
                   qmo(i+1,j,k,QPRES  ) = qm(i+1,j,k,QPRES)
                   if (ppm_type == 0) then
                      qmo(i+1,j,k,QPRES  ) = qmo(i+1,j,k,QPRES  ) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qmo(i+1,j,k,QGAME) = qm(i+1,j,k,QGAME)
                endif

                qmo(i+1,j,k,QPRES) = max(qmo(i+1,j,k,QPRES), small_pres)

                call reset_edge_state_thermo(qmo, q_lo, q_hi, i+1, j, k)

#ifdef RADIATION
                qmo(i+1,j,k,qrad:qradhi) = ernewl(:)
                qmo(i+1,j,k,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i+1,j,k,QPRES)
                qmo(i+1,j,k,qreitot) = sum(qmo(i+1,j,k,qrad:qradhi)) + qmo(i+1,j,k,QREINT)
#endif

             end if

          end do
       end do
    end do

  end subroutine transyz

  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz(qm, qmo, qp, qpo, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxz, &
#ifdef RADIATION
                     rfxz, &
#endif
                     fx_lo, fx_hi, &
                     fzx, &
#ifdef RADIATION
                     rfzx, &
#endif
                     fz_lo, fz_hi, &
                     qx, qx_lo, qx_hi, &
                     qz, qz_lo, qz_hi, &
                     srcQ, src_lo, src_hi, &
                     hdt, cdtdx, cdtdz, lo, hi)


    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, ppm_type, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer, intent(in) :: q_lo(3),q_hi(3)
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fx_lo(3),fx_hi(3)
    integer, intent(in) :: fz_lo(3),fz_hi(3)
    integer, intent(in) :: qx_lo(3),qx_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: src_lo(3),src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: hdt, cdtdx, cdtdz

#ifdef RADIATION
    real(rt), intent(in) :: rfxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) :: qp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxz(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: fzx(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt), intent(in) :: srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt)         rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt)         rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt)         rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt)         pnewr, pnewl
    real(rt)         pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt)         pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    real(rt)         uxav, gexav, dgex, uzav, gezav, dgez
    real(rt)         compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt)         :: dmx, dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2)-1, hi(2)+1
             do i = lo(1), hi(1)
                if (j >= lo(2)) then
                   rrr = qp(i,j,k,QRHO)
                   compr = rrr*qp(i,j,k,nqp)
                   rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                                - cdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
                   compnr = compr - cdtdx*(fxz(i+1,j,k,n) - fxz(i,j,k,n)) &
                                  - cdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))

                   qpo(i,j  ,k,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k,nqp)
                end if

                if (j <= hi(2)) then
                   rrl = qm(i,j+1,k,QRHO)
                   compl = rrl*qm(i,j+1,k,nqp)
                   rrnewl = rrl - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                                - cdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
                   compnl = compl - cdtdx*(fxz(i+1,j,k,n) - fxz(i,j,k,n)) &
                                  - cdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))

                   qmo(i,j+1,k,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k,nqp)
                endif
             end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add the transverse xz and zx differences to the y-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2)-1, hi(2)+1
          do i = lo(1), hi(1)

             pgxp  = qx(i+1,j,k,GDPRES)
             pgxm  = qx(i,j,k,GDPRES)
             ugxp  = qx(i+1,j,k,GDU)
             ugxm  = qx(i,j,k,GDU)
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i,j,k,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgzp  = qz(i,j,k+1,GDPRES)
             pgzm  = qz(i,j,k,GDPRES)
             ugzp  = qz(i,j,k+1,GDW)
             ugzm  = qz(i,j,k,GDW)
             gegzp = qz(i,j,k+1,GDGAME)
             gegzm = qz(i,j,k,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )
#endif

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j,k,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)

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

             if (j >= lo(2)) then

                ! Convert to conservation form
                rrr = qp(i,j,k,QRHO)
                rur = rrr*qp(i,j,k,QU)
                rvr = rrr*qp(i,j,k,QV)
                rwr = rrr*qp(i,j,k,QW)
                ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
                rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
                err = qp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                             - cdtdz*(fzx(i,j,k+1,URHO) - fzx(i,j,k,URHO))
                runewr = rur - cdtdx*(fxz(i+1,j,k,UMX) - fxz(i,j,k,UMX)) &
                             - cdtdz*(fzx(i,j,k+1,UMX) - fzx(i,j,k,UMX))
                rvnewr = rvr - cdtdx*(fxz(i+1,j,k,UMY) - fxz(i,j,k,UMY)) &
                             - cdtdz*(fzx(i,j,k+1,UMY) - fzx(i,j,k,UMY))
                rwnewr = rwr - cdtdx*(fxz(i+1,j,k,UMZ) - fxz(i,j,k,UMZ)) &
                             - cdtdz*(fzx(i,j,k+1,UMZ) - fzx(i,j,k,UMZ))
                renewr = rer - cdtdx*(fxz(i+1,j,k,UEDEN) - fxz(i,j,k,UEDEN)) &
                             - cdtdz*(fzx(i,j,k+1,UEDEN) - fzx(i,j,k,UEDEN))
#ifdef RADIATION
                runewr = runewr + dmx
                rwnewr = rwnewr + dmz
                renewr = renewr + dre
                ernewr = err(:) - cdtdx*(rfxz(i+1,j,k,:) - rfxz(i,j,k,:)) &
                                - cdtdz*(rfzx(i  ,j,k+1,:) - rfzx(i,j,k,:)) &
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

                ! Convert back to primitive form
                qpo(i,j,k,QRHO  ) = rrnewr
                qpo(i,j,k,QU    ) = runewr/rrnewr
                qpo(i,j,k,QV    ) = rvnewr/rrnewr
                qpo(i,j,k,QW    ) = rwnewr/rrnewr

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qpo(i,j,k,QRHO  ) = qpo(i,j,k,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qpo(i,j,k,QU:QW) = qpo(i,j,k,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                qpo(i,j,k,QREINT) = renewr - rhoekenr
                if (ppm_type == 0) then
                   qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating
                      ! (rho e).
                      qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                           - cdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                           - cdtdz*(fzx(i  ,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewr = qp(i,j,k,QPRES) - pxnew - pznew
                      qpo(i,j,k,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geznew

                      ! and compute the p edge state from this and (rho e)
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                   endif
                else
                   qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                   if (ppm_type == 0) then
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
                endif

                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

                call reset_edge_state_thermo(qpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qpo(i,j,k,qrad:qradhi) = ernewr(:)
                qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
                qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
#endif
             endif

             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             if (j <= hi(2)) then

                ! Convert to conservation form
                rrl = qm(i,j+1,k,QRHO)
                rul = rrl*qm(i,j+1,k,QU)
                rvl = rrl*qm(i,j+1,k,QV)
                rwl = rrl*qm(i,j+1,k,QW)
                ekenl = HALF*rrl*sum(qm(i,j+1,k,QU:QW)**2)
                rel = qm(i,j+1,k,QREINT) + ekenl
#ifdef RADIATION
                erl = qm(i,j+1,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewl = rrl - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                             - cdtdz*(fzx(i,j,k+1,URHO) - fzx(i,j,k,URHO))
                runewl = rul - cdtdx*(fxz(i+1,j,k,UMX) - fxz(i,j,k,UMX)) &
                             - cdtdz*(fzx(i,j,k+1,UMX) - fzx(i,j,k,UMX))
                rvnewl = rvl - cdtdx*(fxz(i+1,j,k,UMY) - fxz(i,j,k,UMY)) &
                             - cdtdz*(fzx(i,j,k+1,UMY) - fzx(i,j,k,UMY))
                rwnewl = rwl - cdtdx*(fxz(i+1,j,k,UMZ) - fxz(i,j,k,UMZ)) &
                             - cdtdz*(fzx(i,j,k+1,UMZ) - fzx(i,j,k,UMZ))
                renewl = rel - cdtdx*(fxz(i+1,j,k,UEDEN) - fxz(i,j,k,UEDEN)) &
                             - cdtdz*(fzx(i,j,k+1,UEDEN) - fzx(i,j,k,UEDEN))
#ifdef RADIATION
                runewl = runewl + dmx
                rwnewl = rwnewl + dmz
                renewl = renewl + dre
                ernewl = erl(:) - cdtdx*(rfxz(i+1,j,k,:) - rfxz(i,j,k,:)) &
                                - cdtdz*(rfzx(i  ,j,k+1,:) - rfzx(i,j,k,:)) &
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

                ! Convert back to primitive form
                qmo(i,j+1,k,QRHO  ) = rrnewl
                qmo(i,j+1,k,QU    ) = runewl/rrnewl
                qmo(i,j+1,k,QV    ) = rvnewl/rrnewl
                qmo(i,j+1,k,QW    ) = rwnewl/rrnewl

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qmo(i,j+1,k,QRHO  ) = qmo(i,j+1,k,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qmo(i,j+1,k,QU:QW) = qmo(i,j+1,k,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
                qmo(i,j+1,k,QREINT) = renewl - rhoekenl
                if (ppm_type == 0) then
                   qmo(i,j+1,k,QREINT) = qmo(i,j+1,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qmo(i,j+1,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qmo(i,j+1,k,QREINT) = qm(i,j+1,k,QREINT) &
                           - cdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                           - cdtdz*(fzx(i,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz)
                      if (ppm_type == 0) then
                         qmo(i,j+1,k,QREINT) = qmo(i,j+1,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewl = qm(i,j+1,k,QPRES) - pxnew - pznew
                      qmo(i,j+1,k,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i,j+1,k,QPRES) = qmo(i,j+1,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qmo(i,j+1,k,QGAME) = qm(i,j+1,k,QGAME) + gexnew + geznew

                      ! and compute the p edge state from this and (rho e)
                      qmo(i,j+1,k,QPRES) = qmo(i,j+1,k,QREINT)*(qmo(i,j+1,k,QGAME)-ONE)
                   endif
                else
                   qmo(i,j+1,k,QPRES) = qm(i,j+1,k,QPRES)
                   if (ppm_type == 0) then
                      qmo(i,j+1,k,QPRES) = qmo(i,j+1,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qmo(i,j+1,k,QGAME) = qm(i,j+1,k,QGAME)
                endif

                qmo(i,j+1,k,QPRES) = max(qmo(i,j+1,k,QPRES), small_pres)

                call reset_edge_state_thermo(qmo, q_lo, q_hi, i, j+1, k)

#ifdef RADIATION
                qmo(i,j+1,k,qrad:qradhi) = ernewl(:)
                qmo(i,j+1,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j+1,k,QPRES)
                qmo(i,j+1,k,qreitot) = sum(qmo(i,j+1,k,qrad:qradhi)) + qmo(i,j+1,k,QREINT)
#endif

             end if

          end do
       end do
    end do

  end subroutine transxz

  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy(qm, qmo, qp, qpo, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxy, &
#ifdef RADIATION
                     rfxy, &
#endif
                     fx_lo, fx_hi, &
                     fyx, &
#ifdef RADIATION
                     rfyx, &
#endif
                     fy_lo, fy_hi, &
                     qx, qx_lo, qx_hi, &
                     qy, qy_lo, qy_hi, &
                     srcQ, src_lo, src_hi, &
                     hdt, cdtdx, cdtdy, lo, hi)


    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, ppm_type, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t


    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: src_lo(3), src_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) :: hdt, cdtdx, cdtdy

#ifdef RADIATION
    real(rt), intent(in) :: rfxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) ::   qm(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in) ::   qp(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qmo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(out) :: qpo(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxy(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: fyx(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)


    integer i, j, k, n, nqp, ipassive

    real(rt)         rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    real(rt)         rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    real(rt)         rrnewr, runewr, rvnewr, rwnewr, renewr
    real(rt)         rrnewl, runewl, rvnewl, rwnewl, renewl
    real(rt)         pnewr, pnewl
    real(rt)         pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    real(rt)         pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    real(rt)         uxav, gexav, dgex, uyav, geyav, dgey
    real(rt)         compr, compl, compnr, compnl

#ifdef RADIATION
    real(rt)         :: dmx, dmy, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3)-1, hi(3)+1
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (k >= lo(3)) then
                   rrr = qp(i,j,k,QRHO)
                   compr = rrr*qp(i,j,k,nqp)
                   rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                                - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                   compnr = compr - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                                  - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))

                   qpo(i,j,k,nqp) = compnr/rrnewr + hdt*srcQ(i,j,k  ,nqp)
                end if

                if (k <= hi(3)) then
                   rrl = qm(i,j,k+1,QRHO)
                   compl = rrl*qm(i,j,k+1,nqp)
                   rrnewl = rrl - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                                - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                   compnl = compl - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                                  - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))

                   qmo(i,j,k+1,nqp) = compnl/rrnewl + hdt*srcQ(i,j,k,nqp)
                end if
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse xy and yx differences to the z-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! qzpo state
    !-------------------------------------------------------------------

    do k = lo(3)-1, hi(3)+1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             pgxp = qx(i+1,j,k,GDPRES)
             pgxm = qx(i,j,k,GDPRES)
             ugxp = qx(i+1,j,k,GDU)
             ugxm = qx(i,j,k,GDU)
             gegxp = qx(i+1,j,k,GDGAME)
             gegxm = qx(i,j,k,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgyp = qy(i,j+1,k,GDPRES)
             pgym = qy(i,j,k,GDPRES)
             ugyp = qy(i,j+1,k,GDV)
             ugym = qy(i,j,k,GDV)
             gegyp = qy(i,j+1,k,GDGAME)
             gegym = qy(i,j,k,GDGAME)
#ifdef RADIATION
             ergyp = qy(i,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k,QGAMC))*dux - uxav*dgex )
#endif

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k,QGAMC))*duy - uyav*dgey )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)

             lgex = lambda(:) * (ergxp(:)-ergxm(:))
             lgey = lambda(:) * (ergyp(:)-ergym(:))
             dmx = - cdtdx*sum(lgex)
             dmy = - cdtdy*sum(lgey)
             lugex = HALF*(ugxp+ugxm) * lgex(:)
             lugey = HALF*(ugyp+ugym) * lgey(:)
             dre = -cdtdx*sum(lugex) - cdtdy*sum(lugey)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ugxp+ugxm)*(ergxp(g)-ergxm(g)) &
                        +       cdtdy*HALF*(ugyp+ugym)*(ergyp(g)-ergym(g)) )
                end do
             else if (fspace_type .eq. 2) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = f1*(cdtdx*HALF*(ergxp(g)+ergxm(g))*(ugxm-ugxp) &
                        +       cdtdy*HALF*(ergyp(g)+ergym(g))*(ugym-ugyp) )
                end do
             else ! mixed frame
                der(:) = cdtdx * lugex + cdtdy * lugey
             end if
#endif

             if (k >= lo(3)) then

                ! Convert to conservation form
                rrr = qp(i,j,k,QRHO)
                rur = rrr*qp(i,j,k,QU)
                rvr = rrr*qp(i,j,k,QV)
                rwr = rrr*qp(i,j,k,QW)
                ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
                rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
                err = qp(i,j,k,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                runewr = rur - cdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                             - cdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
                rvnewr = rvr - cdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                             - cdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
                rwnewr = rwr - cdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                             - cdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
                renewr = rer - cdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                             - cdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))
#ifdef RADIATION
                runewr = runewr + dmx
                rvnewr = rvnewr + dmy
                renewr = renewr + dre
                ernewr = err(:) - cdtdx*(rfxy(i+1,j,k,:) - rfxy(i,j,k,:)) &
                                - cdtdy*(rfyx(i,j+1,k,:) - rfyx(i,j,k,:))  &
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

                ! Convert back to primitive form
                qpo(i,j,k,QRHO  ) = rrnewr
                qpo(i,j,k,QU    ) = runewr/rrnewr
                qpo(i,j,k,QV    ) = rvnewr/rrnewr
                qpo(i,j,k,QW    ) = rwnewr/rrnewr

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qpo(i,j,k,QRHO ) = qpo(i,j,k,QRHO ) + hdt*srcQ(i,j,k,QRHO)
                   qpo(i,j,k,QU:QW) = qpo(i,j,k,QU:QW) + hdt*srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
                qpo(i,j,k,QREINT) = renewr - rhoekenr
                if (ppm_type == 0) then
                   qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by
                      ! using the discretized expression for updating
                      ! (rho e).
                      qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                           - cdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                           - cdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy)
                      if (ppm_type == 0) then
                         qpo(i,j,k,QREINT) = qpo(i,j,k,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif


                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewr = qp(i,j,k,QPRES) - pxnew - pynew
                      qpo(i,j,k,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geynew

                      ! and compute the p edge state from this and (rho e)
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                   endif
                else
                   qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                   if (ppm_type == 0) then
                      qpo(i,j,k,QPRES) = qpo(i,j,k,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
                endif

                qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

                call reset_edge_state_thermo(qpo, q_lo, q_hi, i, j, k)

#ifdef RADIATION
                qpo(i,j,k,qrad:qradhi) = ernewr(:)
                qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
                qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
#endif

             end if

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             if (k <= hi(3)) then

                ! Convert to conservation form
                rrl = qm(i,j,k+1,QRHO)
                rul = rrl*qm(i,j,k+1,QU)
                rvl = rrl*qm(i,j,k+1,QV)
                rwl = rrl*qm(i,j,k+1,QW)
                ekenl = HALF*rrl*sum(qm(i,j,k+1,QU:QW)**2)
                rel = qm(i,j,k+1,QREINT) + ekenl
#ifdef RADIATION
                erl = qm(i,j,k+1,qrad:qradhi)
#endif

                ! Add transverse predictor
                rrnewl = rrl - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                runewl = rul - cdtdx*(fxy(i+1,j,k,UMX) - fxy(i,j,k,UMX)) &
                             - cdtdy*(fyx(i,j+1,k,UMX) - fyx(i,j,k,UMX))
                rvnewl = rvl - cdtdx*(fxy(i+1,j,k,UMY) - fxy(i,j,k,UMY)) &
                             - cdtdy*(fyx(i,j+1,k,UMY) - fyx(i,j,k,UMY))
                rwnewl = rwl - cdtdx*(fxy(i+1,j,k,UMZ) - fxy(i,j,k,UMZ)) &
                             - cdtdy*(fyx(i,j+1,k,UMZ) - fyx(i,j,k,UMZ))
                renewl = rel - cdtdx*(fxy(i+1,j,k,UEDEN) - fxy(i,j,k,UEDEN)) &
                             - cdtdy*(fyx(i,j+1,k,UEDEN) - fyx(i,j,k,UEDEN))
#ifdef RADIATION
                runewl = runewl + dmx
                rvnewl = rvnewl + dmy
                renewl = renewl + dre
                ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,k,:) - rfxy(i,j,k,:)) &
                                - cdtdy*(rfyx(i  ,j+1,k,:) - rfyx(i,j,k,:)) &
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

                ! Convert back to primitive form
                qmo(i,j,k+1,QRHO  ) = rrnewl
                qmo(i,j,k+1,QU    ) = runewl/rrnewl
                qmo(i,j,k+1,QV    ) = rvnewl/rrnewl
                qmo(i,j,k+1,QW    ) = rwnewl/rrnewl

                ! for ppm_type > 0 we already added the piecewise parabolic traced
                ! source terms to the normal edge states.
                if (ppm_type == 0) then
                   qmo(i,j,k+1,QRHO  ) = qmo(i,j,k+1,QRHO  ) + hdt*srcQ(i,j,k,QRHO)
                   qmo(i,j,k+1,QU:QW) = qmo(i,j,k+1,QU:QW) + hdt * srcQ(i,j,k,QU:QW)
                endif

                ! note: we run the risk of (rho e) being negative here
                rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
                qmo(i,j,k+1,QREINT) = renewl - rhoekenl
                if (ppm_type == 0) then
                   qmo(i,j,k+1,QREINT) = qmo(i,j,k+1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                endif

                if (.not. reset_state) then
                   if (transverse_reset_rhoe == 1 .and. qmo(i,j,k+1,QREINT) <= ZERO) then
                      ! If it is negative, reset the internal energy by using the discretized
                      ! expression for updating (rho e).
                      qmo(i,j,k+1,QREINT) = qm(i,j,k+1,QREINT) &
                           - cdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                           - cdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy)
                      if (ppm_type == 0) then
                         qmo(i,j,k+1,QREINT) = qmo(i,j,k+1,QREINT) + hdt*srcQ(i,j,k,QREINT)
                      endif
                   endif

                   ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                   ! If we are wrong, we will fix it later

                   if (ppm_predict_gammae == 0) then
                      ! add the transverse term to the p evolution eq here
                      pnewl = qm(i,j,k+1,QPRES) - pxnew - pynew
                      qmo(i,j,k+1,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i,j,k+1,QPRES) = qmo(i,j,k+1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                      endif
                   else
                      ! Update gammae with its transverse terms
                      qmo(i,j,k+1,QGAME) = qm(i,j,k+1,QGAME) + gexnew + geynew

                      ! and compute the p edge state from this and (rho e)
                      qmo(i,j,k+1,QPRES) = qmo(i,j,k+1,QREINT)*(qmo(i,j,k+1,QGAME)-ONE)
                   endif
                else
                   qmo(i,j,k+1,QPRES) = qm(i,j,k+1,QPRES)
                   if (ppm_type == 0) then
                      qmo(i,j,k+1,QPRES) = qmo(i,j,k+1,QPRES) + hdt*srcQ(i,j,k,QPRES)
                   endif
                   qmo(i,j,k+1,QGAME) = qm(i,j,k+1,QGAME)
                endif

                qmo(i,j,k+1,QPRES) = max(qmo(i,j,k+1,QPRES), small_pres)

                call reset_edge_state_thermo(qmo, q_lo, q_hi, i, j, k+1)

#ifdef RADIATION
                qmo(i,j,k+1,qrad:qradhi) = ernewl(:)
                qmo(i,j,k+1,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k+1,QPRES)
                qmo(i,j,k+1,qreitot) = sum(qmo(i,j,k+1,qrad:qradhi)) + qmo(i,j,k+1,QREINT)
#endif

             end if
          end do
       end do
    end do

  end subroutine transxy

  subroutine reset_edge_state_thermo(qedge, qd_lo, qd_hi, ii, jj, kk)

  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, QVAR, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 fspace_type, comoving, &
                                 GDERADS, GDLAMS, &
                                 QCG, QGAMCG, QLAMS, &
#endif
                                 URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, &
                                 NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                 small_pres, small_temp, &
                                 npassive, upass_map, qpass_map, &
                                 ppm_predict_gammae, ppm_type, &
                                 transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer, intent(in) :: ii, jj, kk
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    real(rt)        , intent(inout) :: qedge(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)

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

end module transverse_module
