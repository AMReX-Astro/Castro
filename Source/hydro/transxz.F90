  subroutine transxz(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fxz, fxz_lo, fxz_hi, &
#ifdef RADIATION
                     rfxz, rfxz_lo, rfxz_hi, &
#endif
                     fzx, fzx_lo, fzx_hi, &
#ifdef RADIATION
                     rfzx, rfzx_lo, rfzx_hi, &
#endif
                     qx, qx_lo, qx_hi, &
                     qz, qz_lo, qz_hi, &
                     hdt, cdtdx, cdtdz) bind(C, name="transxz")

    ! here, lo and hi are the bounds of the y interfaces we are looping over

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
    integer, intent(in) :: qa_lo(3),qa_hi(3)
    integer, intent(in) :: fxz_lo(3), fxz_hi(3)
    integer, intent(in) :: fzx_lo(3), fzx_hi(3)
    integer, intent(in) :: qx_lo(3),qx_hi(3)
    integer, intent(in) :: qz_lo(3),qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdx, cdtdz

#ifdef RADIATION
    integer, intent(in) :: rfxz_lo(3), rfxz_hi(3)
    integer, intent(in) :: rfzx_lo(3), rfzx_hi(3)
    real(rt), intent(in) :: rfxz(rfxz_lo(1):rfxz_hi(1),rfxz_lo(2):rfxz_hi(2),rfxz_lo(3):rfxz_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfzx(rfzx_lo(1):rfzx_hi(1),rfzx_lo(2):rfzx_hi(2),rfzx_lo(3):rfzx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxz(fxz_lo(1):fxz_hi(1),fxz_lo(2):fxz_hi(2),fxz_lo(3):fxz_hi(3),NVAR)
    real(rt), intent(in) :: fzx(fzx_lo(1):fzx_hi(1),fzx_lo(2):fzx_hi(2),fzx_lo(3):fzx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

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
    real(rt) :: dmx, dmz, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, lugex, lugez, lgex, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergxp, ergzm,  ergxm
    real(rt) :: eddf, f1
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
                rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                             - cdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxz(i+1,j,k,n) - fxz(i,j,k,n)) &
                               - cdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))

                qpo(i,j  ,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                             - cdtdz*(fzx(i  ,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
                compnl = compl - cdtdx*(fxz(i+1,j-1,k,n) - fxz(i,j-1,k,n)) &
                               - cdtdz*(fzx(i  ,j-1,k+1,n) - fzx(i,j-1,k,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
             end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add the transverse xz and zx differences to the y-states for the
    ! fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qypo state
             !-------------------------------------------------------------------

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

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - cdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i  ,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pxnew - pznew
                   qpo(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + gexnew + geznew

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
             ! qymo state
             !-------------------------------------------------------------------

             pgxp  = qx(i+1,j-1,k,GDPRES)
             pgxm  = qx(i,j-1,k,GDPRES)
             ugxp  = qx(i+1,j-1,k,GDU)
             ugxm  = qx(i,j-1,k,GDU)
             gegxp = qx(i+1,j-1,k,GDGAME)
             gegxm = qx(i,j-1,k,GDGAME)
#ifdef RADIATION
             ergxp = qx(i+1,j-1,k,GDERADS:GDERADS-1+ngroups)
             ergxm = qx(i  ,j-1,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgzp  = qz(i,j-1,k+1,GDPRES)
             pgzm  = qz(i,j-1,k,GDPRES)
             ugzp  = qz(i,j-1,k+1,GDW)
             ugzm  = qz(i,j-1,k,GDW)
             gegzp = qz(i,j-1,k+1,GDGAME)
             gegzm = qz(i,j-1,k,GDGAME)
#ifdef RADIATION
             ergzp = qz(i,j-1,k+1,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i,j-1,k,GDERADS:GDERADS-1+ngroups)
#endif

             duxp = pgxp*ugxp - pgxm*ugxm
             pxav = HALF*(pgxp+pgxm)
             uxav = HALF*(ugxp+ugxm)
             gexav = HALF*(gegxp+gegxm)
             dux = ugxp-ugxm
             dgex = gegxp-gegxm
#ifdef RADIATION
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j-1,k,QGAMCG) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j-1,k,QGAMCG))*dux - uxav*dgex )
#else
             pxnew = cdtdx*(duxp + pxav*dux*(qaux(i,j-1,k,QGAMC) - ONE))
             gexnew = cdtdx*( (gexav-ONE)*(gexav - qaux(i,j-1,k,QGAMC))*dux - uxav*dgex )
#endif

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j-1,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j-1,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i,j-1,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i,j-1,k,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i,j-1,k,QLAMS:QLAMS+ngroups-1)

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
             rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                          - cdtdz*(fzx(i,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
             runewl = rul - cdtdx*(fxz(i+1,j-1,k,UMX) - fxz(i,j-1,k,UMX)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMX) - fzx(i,j-1,k,UMX))
             rvnewl = rvl - cdtdx*(fxz(i+1,j-1,k,UMY) - fxz(i,j-1,k,UMY)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMY) - fzx(i,j-1,k,UMY))
             rwnewl = rwl - cdtdx*(fxz(i+1,j-1,k,UMZ) - fxz(i,j-1,k,UMZ)) &
                          - cdtdz*(fzx(i,j-1,k+1,UMZ) - fzx(i,j-1,k,UMZ))
             renewl = rel - cdtdx*(fxz(i+1,j-1,k,UEDEN) - fxz(i,j-1,k,UEDEN)) &
                          - cdtdz*(fzx(i,j-1,k+1,UEDEN) - fzx(i,j-1,k,UEDEN))
#ifdef RADIATION
             runewl = runewl + dmx
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdx*(rfxz(i+1,j-1,k,:) - rfxz(i,j-1,k,:)) &
                             - cdtdz*(rfzx(i  ,j-1,k+1,:) - rfzx(i,j-1,k,:)) &
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
                        - cdtdx*(fxz(i+1,j-1,k,UEINT) - fxz(i,j-1,k,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i,j-1,k+1,UEINT) - fzx(i,j-1,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j,k,QPRES) - pxnew - pznew
                   qmo(i,j,k,QPRES) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + gexnew + geznew

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

  end subroutine transxz
