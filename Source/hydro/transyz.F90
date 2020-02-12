  subroutine transyz(lo, hi, &
                     qm, qm_lo, qm_hi, &
                     qmo, qmo_lo, qmo_hi, &
                     qp, qp_lo, qp_hi, &
                     qpo, qpo_lo, qpo_hi, &
                     qaux, qa_lo, qa_hi, &
                     fyz, fyz_lo, fyz_hi, &
#ifdef RADIATION
                     rfyz, rfyz_lo, rfyz_hi, &
#endif
                     fzy, fzy_lo, fzy_hi, &
#ifdef RADIATION
                     rfzy, rfzy_lo, rfzy_hi, &
#endif
                     qy, qy_lo, qy_hi, &
                     qz, qz_lo, qz_hi, &
                     hdt, cdtdy, cdtdz) bind(C, name="transyz")

    ! here, lo and hi are the bounds of the x interfaces we are looping over

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
    integer, intent(in) :: fyz_lo(3), fyz_hi(3)
    integer, intent(in) :: fzy_lo(3), fzy_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: qz_lo(3), qz_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in), value :: hdt, cdtdy, cdtdz

#ifdef RADIATION
    integer, intent(in) :: rfyz_lo(3), rfyz_hi(3)
    integer, intent(in) :: rfzy_lo(3), rfzy_hi(3)
    real(rt), intent(in) :: rfyz(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3),0:ngroups-1)
    real(rt), intent(in) :: rfzy(fzy_lo(1):fzy_hi(1),fzy_lo(2):fzy_hi(2),fzy_lo(3):fzy_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(in) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fyz(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3),NVAR)
    real(rt), intent(in) :: fzy(fzy_lo(1):fzy_hi(1),fzy_lo(2):fzy_hi(2),fzy_lo(3):fzy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

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
    real(rt) :: dmy, dmz, dre
    real(rt), dimension(0:ngroups-1) :: der, lambda, lugey, lugez, lgey, lgez, &
         err, ernewr, erl, ernewl, ergzp, ergyp, ergzm, ergym
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
                rrnewr = rrr - cdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                             - cdtdz*(fzy(i,j  ,k+1,URHO) - fzy(i,j,k,URHO))
                compnr = compr - cdtdy*(fyz(i,j+1,k,n   ) - fyz(i,j,k,n)) &
                               - cdtdz*(fzy(i,j  ,k+1,n   ) - fzy(i,j,k,n))

                qpo(i  ,j,k,nqp) = compnr/rrnewr

                rrl = qm(i,j,k,QRHO)
                compl = rrl*qm(i,j,k,nqp)
                rrnewl = rrl - cdtdy*(fyz(i-1,j+1,k,URHO) - fyz(i-1,j,k,URHO)) &
                             - cdtdz*(fzy(i-1,j  ,k+1,URHO) - fzy(i-1,j,k,URHO))
                compnl = compl - cdtdy*(fyz(i-1,j+1,k,n   ) - fyz(i-1,j,k,n)) &
                               - cdtdz*(fzy(i-1,j  ,k+1,n   ) - fzy(i-1,j,k,n))

                qmo(i,j,k,nqp) = compnl/rrnewl
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
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qxpo state
             !-------------------------------------------------------------------

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

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) &
                        - cdtdy*(fyz(i,j+1,k,UEINT) - fyz(i,j,k,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,k+1,UEINT) - fzy(i,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,k,QPRES) - pynew - pznew
                   qpo(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                end if
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
             ! qxmo state
             !-------------------------------------------------------------------

             pgyp  = qy(i-1,j+1,k,GDPRES)
             pgym  = qy(i-1,j,k,GDPRES)
             ugyp  = qy(i-1,j+1,k,GDV)
             ugym  = qy(i-1,j,k,GDV)
             gegyp = qy(i-1,j+1,k,GDGAME)
             gegym = qy(i-1,j,k,GDGAME)
#ifdef RADIATION
             ergyp = qy(i-1,j+1,k,GDERADS:GDERADS-1+ngroups)
             ergym = qy(i-1,j  ,k,GDERADS:GDERADS-1+ngroups)
#endif

             pgzp  = qz(i-1,j,k+1,GDPRES)
             pgzm  = qz(i-1,j,k,GDPRES)
             ugzp  = qz(i-1,j,k+1,GDW)
             ugzm  = qz(i-1,j,k,GDW)
             gegzp = qz(i-1,j,k+1,GDGAME)
             gegzm = qz(i-1,j,k,GDGAME)
#ifdef RADIATION
             ergzp = qz(i-1,j,k+1,GDERADS:GDERADS-1+ngroups)
             ergzm = qz(i-1,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             duyp = pgyp*ugyp - pgym*ugym
             pyav = HALF*(pgyp+pgym)
             uyav = HALF*(ugyp+ugym)
             geyav = HALF*(gegyp+gegym)
             duy = ugyp-ugym
             dgey = gegyp-gegym
#ifdef RADIATION
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i-1,j,k,QGAMCG) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i-1,j,k,QGAMCG))*duy - uyav*dgey )
#else
             pynew = cdtdy*(duyp + pyav*duy*(qaux(i-1,j,k,QGAMC) - ONE))
             geynew = cdtdy*( (geyav-ONE)*(geyav - qaux(i-1,j,k,QGAMC))*duy - uyav*dgey )
#endif

             duzp = pgzp*ugzp - pgzm*ugzm
             pzav = HALF*(pgzp+pgzm)
             uzav = HALF*(ugzp+ugzm)
             gezav = HALF*(gegzp+gegzm)
             duz = ugzp-ugzm
             dgez = gegzp-gegzm
#ifdef RADIATION
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i-1,j,k,QGAMCG) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i-1,j,k,QGAMCG))*duz - uzav*dgez )
#else
             pznew = cdtdz*(duzp + pzav*duz*(qaux(i-1,j,k,QGAMC) - ONE))
             geznew = cdtdz*( (gezav-ONE)*(gezav - qaux(i-1,j,k,QGAMC))*duz - uzav*dgez )
#endif

#ifdef RADIATION
             lambda(:) = qaux(i-1,j,k,QLAMS:QLAMS+ngroups-1)

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
             rrnewl = rrl - cdtdy*(fyz(i-1,j+1,k,URHO) - fyz(i-1,j,k,URHO)) &
                          - cdtdz*(fzy(i-1,j,k+1,URHO) - fzy(i-1,j,k,URHO))
             runewl = rul - cdtdy*(fyz(i-1,j+1,k,UMX) - fyz(i-1,j,k,UMX)) &
                          - cdtdz*(fzy(i-1,j,k+1,UMX) - fzy(i-1,j,k,UMX))
             rvnewl = rvl - cdtdy*(fyz(i-1,j+1,k,UMY) - fyz(i-1,j,k,UMY)) &
                          - cdtdz*(fzy(i-1,j,k+1,UMY) - fzy(i-1,j,k,UMY))
             rwnewl = rwl - cdtdy*(fyz(i-1,j+1,k,UMZ) - fyz(i-1,j,k,UMZ)) &
                          - cdtdz*(fzy(i-1,j,k+1,UMZ) - fzy(i-1,j,k,UMZ))
             renewl = rel - cdtdy*(fyz(i-1,j+1,k,UEDEN) - fyz(i-1,j,k,UEDEN)) &
                          - cdtdz*(fzy(i-1,j,k+1,UEDEN) - fzy(i-1,j,k,UEDEN))
#ifdef RADIATION
             rvnewl = rvnewl + dmy
             rwnewl = rwnewl + dmz
             renewl = renewl + dre
             ernewl = erl(:) - cdtdy*(rfyz(i-1,j+1,k,:) - rfyz(i-1,j,k,:)) &
                             - cdtdz*(rfzy(i-1,j  ,k+1,:) - rfzy(i-1,j,k,:)) &
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

             qmo(i,j,k,QRHO   ) = rrnewl
             qmo(i,j,k,QU     ) = runewl/rrnewl
             qmo(i,j,k,QV     ) = rvnewl/rrnewl
             qmo(i,j,k,QW     ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j,k,QREINT ) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo(i,j,k,QREINT ) = qm(i,j,k,QREINT) &
                        - cdtdy*(fyz(i-1,j+1,k,UEINT) - fyz(i-1,j,k,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i-1,j  ,k+1,UEINT) - fzy(i-1,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i,j,k,QPRES) - pynew - pznew
                   qmo(i,j,k,QPRES  ) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QREINT)*(qmo(i,j,k,QGAME)-ONE)
                end if
             else
                qmo(i,j,k,QPRES  ) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif

             qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo, qmo_lo, qmo_hi, i, j, k)

#ifdef RADIATION
             qmo(i,j,k,qrad:qrad-1+ngroups) = ernewl(:)
             qmo(i,j,k,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qrad-1+ngroups)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transyz
