  subroutine transxy(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxy, fxy_lo, fxy_hi, &
#ifdef RADIATION
                     rfxy, rfxy_lo, rfxy_hi, &
#endif
                     fyx, fyx_lo, fyx_hi, &
#ifdef RADIATION
                     rfyx, rfyx_lo, rfyx_hi, &
#endif
                     qx, qx_lo, qx_hi, &
                     qy, qy_lo, qy_hi, &
                     hdt, cdtdx, cdtdy) bind(C, name="transxy")

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

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fxy_lo(3), fxy_hi(3)
    integer, intent(in) :: fyx_lo(3), fyx_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdx, cdtdy

#ifdef RADIATION
    integer, intent(in) :: rfxy_lo(3), rfxy_hi(3)
    integer, intent(in) :: rfyx_lo(3), rfyx_hi(3)
    real(rt), intent(in) :: rfxy(rfxy_lo(1):rfxy_hi(1),rfxy_lo(2):rfxy_hi(2),rfxy_lo(3):rfxy_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfyx(rfyx_lo(1):rfyx_hi(1),rfyx_lo(2):rfyx_hi(2),rfyx_lo(3):rfyx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) ::   qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) ::   qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxy(fxy_lo(1):fxy_hi(1),fxy_lo(2):fxy_hi(2),fxy_lo(3):fxy_hi(3),NVAR)
    real(rt), intent(in) :: fyx(fyx_lo(1):fyx_hi(1),fyx_lo(2):fyx_hi(2),fyx_lo(3):fyx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

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
    real(rt) :: dmx, dmy, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, lugex, lugey, lgex, lgey, &
         err, ernewr, erl, ernewl, ergxp, ergyp, ergxm, ergym
    real(rt)         eddf, f1
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
                rrr = qp(i,j,k,QRHO)
                compr = rrr*qp(i,j,k,nqp)
                rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                               - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))

                qpo(i,j,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                             - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
                compnl = compl - cdtdx*(fxy(i+1,j,k-1,n) - fxy(i,j,k-1,n)) &
                               - cdtdy*(fyx(i,j+1,k-1,n) - fyx(i,j,k-1,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
             end do
          end do
       end do

    end do

    !-------------------------------------------------------------------
    ! add the transverse xy and yx differences to the z-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qzpo state
             !-------------------------------------------------------------------

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

             ! Convert to conservation form
             rrr = qp(i,j,k,QRHO)
             rur = rrr*qp(i,j,k,QU)
             rvr = rrr*qp(i,j,k,QV)
             rwr = rrr*qp(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp(i,j,k,QU:QW)**2)
             rer = qp(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp(i,j,k,qrad:qrad-1+ngroups)
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

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy)
                endif


                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pxnew - pynew
                   qpo(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geynew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                endif
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif

             qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo, qpo_lo, qpo_hi, i, j, k)

#ifdef RADIATION
             qpo(i,j,k,qrad:qrad-1+ngroups) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qrad-1+ngroups)) + qpo(i,j,k,QREINT)
#endif


             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             pgxp = qx(i+1,j,k-1,GDPRES)
             pgxm = qx(i,j,k-1,GDPRES)
             ugxp = qx(i+1,j,k-1,GDU)
             ugxm = qx(i,j,k-1,GDU)
             gegxp = qx(i+1,j,k-1,GDGAME)
             gegxm = qx(i,j,k-1,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             pgyp = qy(i,j+1,k-1,GDPRES)
             pgym = qy(i,j,k-1,GDPRES)
             ugyp = qy(i,j+1,k-1,GDV)
             ugym = qy(i,j,k-1,GDV)
             gegyp = qy(i,j+1,k-1,GDGAME)
             gegym = qy(i,j,k-1,GDGAME)
#ifdef RADIATION
             ergyp = qy(i,j+1,k-1,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i,j  ,k-1,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k-1,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k-1,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j,k-1,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j,k-1,QGAMC))*dux - uxav*dgex )
#endif

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k-1,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k-1,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i,j,k-1,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i,j,k-1,QGAMC))*duy - uyav*dgey )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)

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

             ! Convert to conservation form
             rrl = qm(i,j,k,QRHO)
             rul = rrl*qm(i,j,k,QU)
             rvl = rrl*qm(i,j,k,QV)
             rwl = rrl*qm(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm(i,j,k,QU:QW)**2)
             rel = qm(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm(i,j,k,qrad:qrad-1+ngroups)
#endif

             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                          - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
             runewl = rul - cdtdx*(fxy(i+1,j,k-1,UMX) - fxy(i,j,k-1,UMX)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMX) - fyx(i,j,k-1,UMX))
             rvnewl = rvl - cdtdx*(fxy(i+1,j,k-1,UMY) - fxy(i,j,k-1,UMY)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMY) - fyx(i,j,k-1,UMY))
             rwnewl = rwl - cdtdx*(fxy(i+1,j,k-1,UMZ) - fxy(i,j,k-1,UMZ)) &
                          - cdtdy*(fyx(i,j+1,k-1,UMZ) - fyx(i,j,k-1,UMZ))
             renewl = rel - cdtdx*(fxy(i+1,j,k-1,UEDEN) - fxy(i,j,k-1,UEDEN)) &
                          - cdtdy*(fyx(i,j+1,k-1,UEDEN) - fyx(i,j,k-1,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rvnewl = rvnewl + dmy
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxy(i+1,j  ,k-1,:) - rfxy(i,j,k-1,:)) &
                             - cdtdy*(rfyx(i  ,j+1,k-1,:) - rfyx(i,j,k-1,:)) &
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
             qmo(i,j,k,QRHO  ) = rrnewl
             qmo(i,j,k,QU    ) = runewl/rrnewl
             qmo(i,j,k,QV    ) = rvnewl/rrnewl
             qmo(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k-1,UEINT) - fxy(i,j,k-1,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k-1,UEINT) - fyx(i,j,k-1,UEINT) + pyav*duy)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j,k,QPRES) - pxnew - pynew
                   qmo(i,j,k,QPRES) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + gexnew + geynew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QREINT)*(qmo(i,j,k,QGAME)-ONE)
                endif
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif

             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qmo_lo, qmo_hi, i, j, k)

#ifdef RADIATION
             qmo(i,j,k,qrad:qrad-1+ngroups) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qrad-1+ngroups)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxy
