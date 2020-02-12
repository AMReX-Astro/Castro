  subroutine transz_on_ystates(lo, hi, &
                               qym, qym_lo, qym_hi, &
                               qymo, qymo_lo, qymo_hi, &
                               qyp, qyp_lo, qyp_hi, &
                               qypo, qypo_lo, qypo_hi, &
                               qaux, qa_lo, qa_hi, &
                               fz, fz_lo, fz_hi, &
#ifdef RADIATION
                               rfz, rfz_lo, rfz_hi, &
#endif
                               qz, qz_lo, qz_hi, &
                               cdtdz) bind(C, name="transz_on_ystates")

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

    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qymo_lo(3), qymo_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
    integer, intent(in) :: qypo_lo(3), qypo_hi(3)
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

    real(rt), intent(in) :: qym(qym_lo(1):qym_hi(1),qym_lo(2):qym_hi(2),qym_lo(3):qym_hi(3),NQ)
    real(rt), intent(in) :: qyp(qyp_lo(1):qyp_hi(1),qyp_lo(2):qyp_hi(2),qyp_lo(3):qyp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    real(rt), intent(out) :: qymo(qymo_lo(1):qymo_hi(1),qymo_lo(2):qymo_hi(2),qymo_lo(3):qymo_hi(3),NQ)
    real(rt), intent(out) :: qypo(qypo_lo(1):qypo_hi(1),qypo_lo(2):qypo_hi(2),qypo_lo(3):qypo_hi(3),NQ)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew
    real(rt)         rrry, rrly
    real(rt)         rury, ruly
    real(rt)         rvry, rvly
    real(rt)         rwry, rwly
    real(rt)         ekenry, ekenly
    real(rt)         rery, rely
    real(rt)         rrnewry, rrnewly
    real(rt)         runewry, runewly
    real(rt)         rvnewry, rvnewly
    real(rt)         rwnewry, rwnewly
    real(rt)         renewry, renewly
    real(rt)         pnewry, pnewly
    real(rt)         rhoekenry, rhoekenly
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
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1, npassive
       n = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qyp(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nqp) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                qypo(i,j,k,nqp) = compu/rrnew

                rrnew = qym(i,j,k,QRHO) - cdtdz*(fz(i,j-1,k+1,URHO) - fz(i,j-1,k,URHO))
                compu = qym(i,j,k,QRHO)*qym(i,j,k,nqp) - cdtdz*(fz(i,j-1,k+1,n) - fz(i,j-1,k,n))
                qymo(i,j,k,nqp) = compu/rrnew
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the z-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qypo state
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
                   der(g) = cdtdz* f1 * HALF*(ergp(g)+ergm(g))*(ugm-ugp)
                end do
             else ! mixed frame
                der(:) = cdtdz * luge
             end if
#endif

             ! Convert to conservation form
             rrry = qyp(i,j,k,QRHO)
             rury = rrry*qyp(i,j,k,QU)
             rvry = rrry*qyp(i,j,k,QV)
             rwry = rrry*qyp(i,j,k,QW)
             ekenry = HALF*rrry*sum(qyp(i,j,k,QU:QW)**2)
             rery = qyp(i,j,k,QREINT) + ekenry
#ifdef RADIATION
             erry = qyp(i,j,k,qrad:qrad-1+ngroups)
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

             call reset_edge_state_thermo(qypo, qypo_lo, qypo_hi, i, j, k)

#ifdef RADIATION
             qypo(i,j,k,qrad:qrad-1+ngroups) = ernewry(:)
             qypo(i,j,k,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,k,QPRES)
             qypo(i,j,k,qreitot) = sum(qypo(i,j,k,qrad:qrad-1+ngroups)) + qypo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------

             pgp  = qz(i,j-1,k+1,GDPRES)
             pgm  = qz(i,j-1,k,GDPRES)
             ugp  = qz(i,j-1,k+1,GDW   )
             ugm  = qz(i,j-1,k,GDW   )
             gegp = qz(i,j-1,k+1,GDGAME)
             gegm = qz(i,j-1,k,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j-1,k,QLAMS:QLAMS+ngroups-1)
             ergp = qz(i,j-1,k+1,GDERADS:GDERADS-1+ngroups)
             ergm = qz(i,j-1,k,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j-1,k,QGAMCG)
#else
             gamc = qaux(i,j-1,k,QGAMC)
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

             ! Convert to conservation form
             rrly = qym(i,j,k,QRHO)
             ruly = rrly*qym(i,j,k,QU)
             rvly = rrly*qym(i,j,k,QV)
             rwly = rrly*qym(i,j,k,QW)
             ekenly = HALF*rrly*sum(qym(i,j,k,QU:QW)**2)
             rely = qym(i,j,k,QREINT) + ekenly
#ifdef RADIATION
             erly = qym(i,j,k,qrad:qrad-1+ngroups)
#endif

             ! Add transverse predictor
             rrnewly = rrly - cdtdz*(fz(i,j-1,k+1,URHO) - fz(i,j-1,k,URHO))
             runewly = ruly - cdtdz*(fz(i,j-1,k+1,UMX) - fz(i,j-1,k,UMX))
             rvnewly = rvly - cdtdz*(fz(i,j-1,k+1,UMY) - fz(i,j-1,k,UMY))
             rwnewly = rwly - cdtdz*(fz(i,j-1,k+1,UMZ) - fz(i,j-1,k,UMZ))
             renewly = rely - cdtdz*(fz(i,j-1,k+1,UEDEN) - fz(i,j-1,k,UEDEN))
#ifdef RADIATION
             rwnewly = rwnewly + dmz
             renewly = renewly + dre
             ernewly = erly(:) - cdtdz*(rfz(i,j-1,k+1,:) - rfz(i,j-1,k,:)) + der(:)
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
             qymo(i,j,k,QRHO) = rrnewly
             rhoinv = ONE/rrnewly
             qymo(i,j,k,QU) = runewly*rhoinv
             qymo(i,j,k,QV) = rvnewly*rhoinv
             qymo(i,j,k,QW) = rwnewly*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
             qymo(i,j,k,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qymo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qymo(i,j,k,QREINT) = qym(i,j,k,QREINT) - &
                        cdtdz*(fz(i,j-1,k+1,UEINT) - fz(i,j-1,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qymo(i,j,k,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qymo(i,j,k,QGAME) = qym(i,j,k,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qymo(i,j,k,QPRES) = qymo(i,j,k,QREINT)*(qymo(i,j,k,QGAME)-ONE)
                   qymo(i,j,k,QPRES) = max(qymo(i,j,k,QPRES), small_pres)
                endif
             else
                qymo(i,j,k,QPRES) = qym(i,j,k,QPRES)
                qymo(i,j,k,QGAME) = qym(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qymo, qymo_lo, qymo_hi, i, j, k)

#ifdef RADIATION
             qymo(i,j,k,qrad:qrad-1+ngroups) = ernewly(:)
             qymo(i,j,k,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j,k,QPRES)
             qymo(i,j,k,qreitot) = sum(qymo(i,j,k,qrad:qrad-1+ngroups)) + qymo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transz_on_ystates
