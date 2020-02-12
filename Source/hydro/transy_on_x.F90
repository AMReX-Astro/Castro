  subroutine transy_on_xstates(lo, hi, &
                               qxm, qxm_lo, qxm_hi, &
                               qxmo, qxmo_lo, qxmo_hi, &
                               qxp, qxp_lo, qxp_hi, &
                               qxpo, qxpo_lo, qxpo_hi, &
                               qaux, qa_lo, qa_hi, &
                               fy, fy_lo, fy_hi, &
#ifdef RADIATION
                               rfy, rfy_lo, rfy_hi, &
#endif
                               qy, qy_lo, qy_hi, &
                               cdtdy) bind(C, name="transy_on_xstates")

    ! here, lo and hi are the bounds of the edges we are looping over

    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qptot, qreitot, &
                                   fspace_type, comoving, &
                                   GDERADS, GDLAMS, &
                                   QCG, QGAMCG, QLAMS, &
#endif
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NGDNV, GDPRES, GDU, GDV, GDW, GDGAME, &
                                   small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_predict_gammae

#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor ! function
#endif

    integer, intent(in) :: qxm_lo(3), qxm_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
    integer, intent(in) :: qxmo_lo(3), qxmo_hi(3)
    integer, intent(in) :: qxpo_lo(3), qxpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rfy_lo(3), rfy_hi(3)
#endif
    integer, intent(in) :: qy_lo(3),qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

#ifdef RADIATION
    real(rt), intent(in) :: rfy(rfy_lo(1):rfy_hi(1),rfy_lo(2):rfy_hi(2),rfy_lo(3):rfy_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in), value :: cdtdy

    real(rt), intent(in) :: qxm(qxm_lo(1):qxm_hi(1),qxm_lo(2):qxm_hi(2),qxm_lo(3):qxm_hi(3),NQ)
    real(rt), intent(in) :: qxp(qxp_lo(1):qxp_hi(1),qxp_lo(2):qxp_hi(2),qxp_lo(3):qxp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    real(rt), intent(out) :: qxmo(qxmo_lo(1):qxmo_hi(1),qxmo_lo(2):qxmo_hi(2),qxmo_lo(3):qxmo_hi(3),NQ)
    real(rt), intent(out) :: qxpo(qxpo_lo(1):qxpo_hi(1),qxpo_lo(2):qxpo_hi(2),qxpo_lo(3):qxpo_hi(3),NQ)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew
    real(rt)         rrrx, rrlx
    real(rt)         rurx, rulx
    real(rt)         rvrx, rvlx
    real(rt)         rwrx, rwlx
    real(rt)         ekenrx, ekenlx
    real(rt)         rerx, relx
    real(rt)         rrnewrx, rrnewlx
    real(rt)         runewrx, runewlx
    real(rt)         rvnewrx, rvnewlx
    real(rt)         rwnewrx, rwnewlx
    real(rt)         renewrx, renewlx
    real(rt)         pnewrx, pnewlx
    real(rt)         rhoekenrx, rhoekenlx
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

    !$gpu


    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qxp(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = qxp(i,j,k,QRHO)*qxp(i,j,k,nqp) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                qxpo(i,j,k,nqp) = compu/rrnew

                rrnew = qxm(i,j,k,QRHO) - cdtdy*(fy(i-1,j+1,k,URHO) - fy(i-1,j,k,URHO))
                compu = qxm(i,j,k,QRHO)*qxm(i,j,k,nqp) - cdtdy*(fy(i-1,j+1,k,n) - fy(i-1,j,k,n))
                qxmo(i,j,k,nqp) = compu/rrnew
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------
             ! qxpo state
             !-------------------------------------------------------------

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

             ! Convert to conservation form
             rrrx = qxp(i,j,k,QRHO)
             rurx = rrrx*qxp(i,j,k,QU)
             rvrx = rrrx*qxp(i,j,k,QV)
             rwrx = rrrx*qxp(i,j,k,QW)
             ekenrx = HALF*rrrx*sum(qxp(i,j,k,QU:QW)**2)
             rerx = qxp(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
             err  = qxp(i,j,k,qrad:qrad-1+ngroups)
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

             call reset_edge_state_thermo(qxpo, qxpo_lo, qxpo_hi, i, j, k)

#ifdef RADIATION
             qxpo(i,j,k,qrad:qrad-1+ngroups) = ernewr(:)
             qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,k,QPRES)
             qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qrad-1+ngroups)) + qxpo(i,j,k,QREINT)
#endif


             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             pgp  = qy(i-1,j+1,k,GDPRES)
             pgm  = qy(i-1,j  ,k,GDPRES)
             ugp  = qy(i-1,j+1,k,GDV   )
             ugm  = qy(i-1,j  ,k,GDV   )
             gegp = qy(i-1,j+1,k,GDGAME)
             gegm = qy(i-1,j  ,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i-1,j,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i-1,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i-1,j  ,k,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i-1,j,k,QGAMCG)
#else
             gamc = qaux(i-1,j,k,QGAMC)
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

             ! Convert to conservation form
             rrlx = qxm(i,j,k,QRHO)
             rulx = rrlx*qxm(i,j,k,QU)
             rvlx = rrlx*qxm(i,j,k,QV)
             rwlx = rrlx*qxm(i,j,k,QW)
             ekenlx = HALF*rrlx*sum(qxm(i,j,k,QU:QW)**2)
             relx = qxm(i,j,k,QREINT) + ekenlx
#ifdef RADIATION
             erl  = qxm(i,j,k,qrad:qrad-1+ngroups)
#endif

             ! Add transverse predictor
             rrnewlx = rrlx - cdtdy*(fy(i-1,j+1,k,URHO) - fy(i-1,j,k,URHO))
             runewlx = rulx - cdtdy*(fy(i-1,j+1,k,UMX) - fy(i-1,j,k,UMX))
             rvnewlx = rvlx - cdtdy*(fy(i-1,j+1,k,UMY) - fy(i-1,j,k,UMY))
             rwnewlx = rwlx - cdtdy*(fy(i-1,j+1,k,UMZ) - fy(i-1,j,k,UMZ))
             renewlx = relx - cdtdy*(fy(i-1,j+1,k,UEDEN)- fy(i-1,j,k,UEDEN))
#ifdef RADIATION
             rvnewlx = rvnewlx + dmom
             renewlx = renewlx + dre
             ernewl  = erl(:) - cdtdy*(rfy(i-1,j+1,k,:) - rfy(i-1,j,k,:)) + der(:)
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

             qxmo(i,j,k,QRHO) = rrnewlx
             rhoinv = ONE/rrnewlx
             qxmo(i,j,k,QU) = runewlx*rhoinv
             qxmo(i,j,k,QV) = rvnewlx*rhoinv
             qxmo(i,j,k,QW) = rwnewlx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
             qxmo(i,j,k,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo(i,j,k,QREINT) = qxm(i,j,k,QREINT) - &
                        cdtdy*(fy(i-1,j+1,k,UEINT) - fy(i-1,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qxmo(i,j,k,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo(i,j,k,QGAME) = qxm(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo(i,j,k,QPRES) = qxmo(i,j,k,QREINT)*(qxmo(i,j,k,QGAME)-ONE)
                   qxmo(i,j,k,QPRES) = max(qxmo(i,j,k,QPRES), small_pres)
                endif
             else
                qxmo(i,j,k,QPRES) = qxm(i,j,k,QPRES)
                qxmo(i,j,k,QGAME) = qxm(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qxmo, qxmo_lo, qxmo_hi, i, j, k)

#ifdef RADIATION
             qxmo(i,j,k,qrad:qrad-1+ngroups) = ernewl(:)
             qxmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i,j,k,QPRES)
             qxmo(i,j,k,qreitot) = sum(qxmo(i,j,k,qrad:qrad-1+ngroups)) + qxmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transy_on_xstates
