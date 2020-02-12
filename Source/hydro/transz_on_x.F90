  subroutine transz_on_xstates(lo, hi, &
                               qxm, qxm_lo, qxm_hi, &
                               qxmo, qxmo_lo, qxmo_hi, &
                               qxp, qxp_lo, qxp_hi, &
                               qxpo, qxpo_lo, qxpo_hi, &
                               qaux, qa_lo, qa_hi, &
                               fz, fz_lo, fz_hi, &
#ifdef RADIATION
                               rfz, rfz_lo, rfz_hi, &
#endif
                               qz, qz_lo, qz_hi, &
                               cdtdz) bind(C, name="transz_on_xstates")

    ! here, lo and hi are the bounds of edges we are looping over

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
    integer, intent(in) :: qxmo_lo(3), qxmo_hi(3)
    integer, intent(in) :: qxp_lo(3), qxp_hi(3)
    integer, intent(in) :: qxpo_lo(3), qxpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rfz_lo(3), rfz_hi(3)
#endif
    integer, intent(in) :: qz_lo(3), qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

#ifdef RADIATION
    real(rt), intent(in) :: rfz(rfz_lo(1):rfz_hi(1),rfz_lo(2):rfz_hi(2),rfz_lo(3):rfz_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in), value :: cdtdz

    real(rt), intent(in) :: qxm(qxm_lo(1):qxm_hi(1),qxm_lo(2):qxm_hi(2),qxm_lo(3):qxm_hi(3),NQ)
    real(rt), intent(in) :: qxp(qxp_lo(1):qxp_hi(1),qxp_lo(2):qxp_hi(2),qxp_lo(3):qxp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

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
    real(rt)         :: dmz, dre
    real(rt)        , dimension(0:ngroups-1) :: der, lambda, luge, lamge, &
         ergp, errx, ernewrx, erry, ernewry, ergm, erlx, ernewlx, erly, ernewly
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    !$gpu

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transverse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

        do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qxp(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                compu = qxp(i,j,k,QRHO)*qxp(i,j,k,nqp) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                qxpo(i,j,k,nqp) = compu/rrnew

                rrnew = qxm(i,j,k,QRHO) - cdtdz*(fz(i-1,j,k+1,URHO) - fz(i-1,j,k,URHO))
                compu = qxm(i,j,k,QRHO)*qxm(i,j,k,nqp) - cdtdz*(fz(i-1,j,k+1,n) - fz(i-1,j,k,n))
                qxmo(i,j,k,nqp) = compu/rrnew
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the z-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qxpo state
             !-------------------------------------------------------------------

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

             ! Convert to conservation form
             rrrx = qxp(i,j,k,QRHO)
             rurx = rrrx*qxp(i,j,k,QU)
             rvrx = rrrx*qxp(i,j,k,QV)
             rwrx = rrrx*qxp(i,j,k,QW)
             ekenrx = HALF*rrrx*sum(qxp(i,j,k,QU:QW)**2)
             rerx = qxp(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
             errx = qxp(i,j,k,qrad:qrad-1+ngroups)
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

             call reset_edge_state_thermo(qxpo, qxpo_lo, qxpo_hi, i, j, k)

#ifdef RADIATION
             qxpo(i,j,k,qrad:qrad-1+ngroups) = ernewrx(:)
             qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,k,QPRES)
             qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qrad-1+ngroups)) + qxpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qxmo state
             !-------------------------------------------------------------------

             pgp  = qz(i-1,j,k+1,GDPRES)
             pgm  = qz(i-1,j,k,GDPRES)
             ugp  = qz(i-1,j,k+1,GDW   )
             ugm  = qz(i-1,j,k,GDW   )
             gegp = qz(i-1,j,k+1,GDGAME)
             gegm = qz(i-1,j,k,GDGAME)

#ifdef RADIATION
             lambda(:) = qaux(i-1,j,k,QLAMS:QLAMS+ngroups-1)
             ergp = qz(i-1,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergm = qz(i-1,j,k,GDERADS:GDERADS-1+ngroups)
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

             ! Convert to conservation form
             rrlx = qxm(i,j,k,QRHO)
             rulx = rrlx*qxm(i,j,k,QU)
             rvlx = rrlx*qxm(i,j,k,QV)
             rwlx = rrlx*qxm(i,j,k,QW)
             ekenlx = HALF*rrlx*sum(qxm(i,j,k,QU:QW)**2)
             relx = qxm(i,j,k,QREINT) + ekenlx
#ifdef RADIATION
             erlx = qxm(i,j,k,qrad:qrad-1+ngroups)
#endif

             ! Add transverse predictor
             rrnewlx = rrlx - cdtdz*(fz(i-1,j,k+1,URHO) - fz(i-1,j,k,URHO))
             runewlx = rulx - cdtdz*(fz(i-1,j,k+1,UMX) - fz(i-1,j,k,UMX))
             rvnewlx = rvlx - cdtdz*(fz(i-1,j,k+1,UMY) - fz(i-1,j,k,UMY))
             rwnewlx = rwlx - cdtdz*(fz(i-1,j,k+1,UMZ) - fz(i-1,j,k,UMZ))
             renewlx = relx - cdtdz*(fz(i-1,j,k+1,UEDEN) - fz(i-1,j,k,UEDEN))
#ifdef RADIATION
             rwnewlx = rwnewlx + dmz
             renewlx = renewlx + dre
             ernewlx = erlx(:) - cdtdz*(rfz(i-1,j,k+1,:) - rfz(i-1,j,k,:)) + der(:)
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
             qxmo(i,j,k,QRHO) = rrnewlx
             rhoinv = ONE/rrnewlx
             qxmo(i,j,k,QU) = runewlx*rhoinv
             qxmo(i,j,k,QV) = rvnewlx*rhoinv
             qxmo(i,j,k,QW) = rwnewlx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
             qxmo(i,j,k,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qxmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo(i,j,k,QREINT) = qxm(i,j,k,QREINT) - &
                        cdtdz*(fz(i-1,j,k+1,UEINT) - fz(i-1,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qxmo(i,j,k,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo(i,j,k,QGAME) = qxm(i,j,k,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo(i,j,k,QPRES) = qxmo(i,j,k,QREINT)*(qxmo(i,j,k,QGAME)-ONE)
                   qxmo(i,j,k,QPRES) = max(qxmo(i,j,k,QPRES), small_pres)
                end if
             else
                qxmo(i,j,k,QPRES) = qxm(i,j,k,QPRES)
                qxmo(i,j,k,QGAME) = qxm(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qxmo, qxmo_lo, qxmo_hi, i, j, k)

#ifdef RADIATION
             qxmo(i,j,k,qrad:qrad-1+ngroups) = ernewlx(:)
             qxmo(i,j,k,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i,j,k,QPRES)
             qxmo(i,j,k,qreitot) = sum(qxmo(i,j,k,qrad:qrad-1+ngroups)) + qxmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transz_on_xstates
