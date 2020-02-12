
  subroutine transx_on_zstates(lo, hi, &
                               qzm, qzm_lo, qzm_hi, &
                               qzmo, qzmo_lo, qzmo_hi, &
                               qzp, qzp_lo, qzp_hi, &
                               qzpo, qzpo_lo, qzpo_hi, &
                               qaux, qa_lo, qa_hi, &
                               fx, fx_lo, fx_hi, &
#ifdef RADIATION
                               rfx, rfx_lo, rfx_hi, &
#endif
                               qx, qx_lo, qx_hi, &
                               hdt, cdtdx) bind(C, name="transx_on_zstates")

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
#if AMREX_SPACEDIM == 2
  use prob_params_module, only : mom_flux_has_p
#endif

    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
    integer, intent(in) :: qzmo_lo(3), qzmo_hi(3)
    integer, intent(in) :: qzpo_lo(3), qzpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rfx_lo(3), rfx_hi(3)
#endif
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: lo(3), hi(3)

#ifdef RADIATION
    real(rt) :: rfx(rfx_lo(1):rfx_hi(1),rfx_lo(2):rfx_hi(2),rfx_lo(3):rfx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in), value :: hdt, cdtdx

    real(rt), intent(in) :: qzm(qzm_lo(1):qzm_hi(1),qzm_lo(2):qzm_hi(2),qzm_lo(3):qzm_hi(3),NQ)
    real(rt), intent(in) :: qzp(qzp_lo(1):qzp_hi(1),qzp_lo(2):qzp_hi(2),qzp_lo(3):qzp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

    real(rt), intent(out) :: qzmo(qzmo_lo(1):qzmo_hi(1),qzmo_lo(2):qzmo_hi(2),qzmo_lo(3):qzmo_hi(3),NQ)
    real(rt), intent(out) :: qzpo(qzpo_lo(1):qzpo_hi(1),qzpo_lo(2):qzpo_hi(2),qzpo_lo(3):qzpo_hi(3),NQ)

    integer i, j, k, n, nqp, ipassive

    real(rt)         rhoinv
    real(rt)         rrnew
    real(rt)         rrrz, rrlz
    real(rt)         rurz, rulz
    real(rt)         rvrz, rvlz
    real(rt)         rwrz, rwlz
    real(rt)         ekenrz, ekenlz
    real(rt)         rerz, relz
    real(rt)         rrnewrz, rrnewlz
    real(rt)         runewrz, runewlz
    real(rt)         rvnewrz, rvnewlz
    real(rt)         rwnewrz, rwnewlz
    real(rt)         renewrz, renewlz
    real(rt)         pnewrz, pnewlz
    real(rt)         rhoekenrz, rhoekenlz
    real(rt)         pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         compu
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom, divu
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

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qzp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qzp(i,j,k,QRHO)*qzp(i,j,k,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qzpo(i,j,k,nqp) = compu/rrnew

                rrnew = qzm(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k-1,URHO) - fx(i,j,k-1,URHO))
                compu = qzm(i,j,k,QRHO)*qzm(i,j,k,nqp) - cdtdx*(fx(i+1,j,k-1,n) - fx(i,j,k-1,n))
                qzmo(i,j,k,nqp) = compu/rrnew
             end do
          end do
       end do

    end do


    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to z-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qzpo state
             !-------------------------------------------------------------------

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

             ! Convert to conservation form
             rrrz = qzp(i,j,k,QRHO)
             rurz = rrrz*qzp(i,j,k,QU)
             rvrz = rrrz*qzp(i,j,k,QV)
             rwrz = rrrz*qzp(i,j,k,QW)
             ekenrz = HALF*rrrz*sum(qzp(i,j,k,QU:QW)**2)
             rerz = qzp(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
             err  = qzp(i,j,k,qrad:qrad-1+ngroups)
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

             call reset_edge_state_thermo(qzpo, qzpo_lo, qzpo_hi, i, j, k)

#ifdef RADIATION
             qzpo(i,j,k,qrad:qrad-1+ngroups) = ernewr(:)
             qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
             qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qrad-1+ngroups)) + qzpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             pgp  = qx(i+1,j,k-1,GDPRES)
             pgm  = qx(i  ,j,k-1,GDPRES)
             ugp  = qx(i+1,j,k-1,GDU   )
             ugm  = qx(i  ,j,k-1,GDU   )
             gegp = qx(i+1,j,k-1,GDGAME)
             gegm = qx(i  ,j,k-1,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j,k-1,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j,k-1,GDERADS:GDERADS-1+ngroups)
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
             gamc = qaux(i,j,k-1,QGAMCG)
#else
             gamc = qaux(i,j,k-1,QGAMC)
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

             ! Convert to conservation form
             rrlz = qzm(i,j,k,QRHO)
             rulz = rrlz*qzm(i,j,k,QU)
             rvlz = rrlz*qzm(i,j,k,QV)
             rwlz = rrlz*qzm(i,j,k,QW)
             ekenlz = HALF*rrlz*sum(qzm(i,j,k,QU:QW)**2)
             relz = qzm(i,j,k,QREINT) + ekenlz
#ifdef RADIATION
             erl  = qzm(i,j,k,qrad:qrad-1+ngroups)
#endif

             ! Add transverse predictor
             rrnewlz = rrlz - cdtdx*(fx(i+1,j,k-1,URHO) - fx(i,j,k-1,URHO))
             runewlz = rulz - cdtdx*(fx(i+1,j,k-1,UMX) - fx(i,j,k-1,UMX))
             rvnewlz = rvlz - cdtdx*(fx(i+1,j,k-1,UMY) - fx(i,j,k-1,UMY))
             rwnewlz = rwlz - cdtdx*(fx(i+1,j,k-1,UMZ) - fx(i,j,k-1,UMZ))
             renewlz = relz - cdtdx*(fx(i+1,j,k-1,UEDEN) - fx(i,j,k-1,UEDEN))
#ifdef RADIATION
             runewlz = runewlz + dmom
             renewlz = renewlz + dre
             ernewl  = erl(:) - cdtdx*(rfx(i+1,j,k-1,:) - rfx(i,j,k-1,:)) + der(:)
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
             qzmo(i,j,k,QRHO) = rrnewlz
             rhoinv = ONE/rrnewlz
             qzmo(i,j,k,QU) = runewlz*rhoinv
             qzmo(i,j,k,QV) = rvnewlz*rhoinv
             qzmo(i,j,k,QW) = rwnewlz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
             qzmo(i,j,k,QREINT) = renewlz - rhoekenlz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzmo(i,j,k,QREINT) = qzm(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k-1,UEINT) - fx(i,j,k-1,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qzmo(i,j,k,QPRES) = max(pnewlz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzmo(i,j,k,QPRES) = qzmo(i,j,k,QREINT)*(qzmo(i,j,k,QGAME)-ONE)
                   qzmo(i,j,k,QPRES) = max(qzmo(i,j,k,QPRES), small_pres)
                endif
             else
                qzmo(i,j,k,QPRES) = qzm(i,j,k,QPRES)
                qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzmo, qzmo_lo, qzmo_hi, i, j, k)

#ifdef RADIATION
             qzmo(i,j,k,qrad:qrad-1+ngroups) = ernewl(:)
             qzmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k,QPRES)
             qzmo(i,j,k,qreitot) = sum(qzmo(i,j,k,qrad:qrad-1+ngroups)) + qzmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transx_on_zstates

