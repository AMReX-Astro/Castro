module transverse_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

contains


  ! add the transverse flux difference in direction idir_t to the
  ! interface states in direction idir_n

  subroutine trans_single(lo, hi, &
                          idir_t, idir_n,
                          qm, qm_lo, qm_hi, &
                          qmo, qmo_lo, qmo_hi, &
                          qp, qp_lo, qp_hi, &
                          qpo, qpo_lo, qpo_hi, &
                          qaux, qa_lo, qa_hi, &
                          flux_t, f_lo, f_hi, &
#ifdef RADIATION
                          rflux_t, rf_lo, rf_hi, &
#endif
                          q_t, q_lo, q_hi, &
#if AMREX_SPACEDIM == 2
                          area_t, a_lo, a_hi, &
                          vol, vol_lo, vol_hi, &
#endif
                          hdt, cdtdx) bind(C, name="trans_single")

    ! here, lo and hi are the bounds of the interfaces we are looping over

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

    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qmo_lo(3), qmo_hi(3)
    integer, intent(in) :: qpo_lo(3), qpo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rf_lo(3), rf_hi(3)
#endif
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: lo(3), hi(3)
#if AMREX_SPACEDIM == 2
    integer, intent(in) :: a_lo(3), a_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#endif

#ifdef RADIATION
    real(rt) :: rfx(rfx_lo(1):rfx_hi(1),rfx_lo(2):rfx_hi(2),rfx_lo(3):rfx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in), value :: hdt, cdtdx

    real(rt), intent(in) :: qm(qym_lo(1):qym_hi(1),qym_lo(2):qym_hi(2),qym_lo(3):qym_hi(3),NQ)
    real(rt), intent(in) :: qp(qyp_lo(1):qyp_hi(1),qyp_lo(2):qyp_hi(2),qyp_lo(3):qyp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)
    real(rt), intent(in) :: q_t(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: a_t(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
#endif

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
    real(rt)         rhoekenry, rhoekenly, rhoekenrz, rhoekenlz
    real(rt)         pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt)         compu
    real(rt)         :: gamc

#ifdef RADIATION
    real(rt)         :: dre, dmom, divu
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt)         eddf, f1
    integer :: g
#endif

    logical :: reset_state

    real(rt) :: volinv

    !$gpu

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transverse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------


    !       qm|qp
    !         |
    ! --------+--------
    !   i-1       i
    !        i-1/2
    !
    ! the qm state will see the transverse flux in zone i-1

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             do ipassive = 1, npassive
                n  = upass_map(ipassive)
                nqp = qpass_map(ipassive)

#if AMREX_SPACEDIM == 2
                volinv = ONE/vol(i,j,k)

                rrnew = qp(i,j,k,QRHO) - hdt*(area_t(i+1,j,k)*fx(i+1,j,k,URHO) - &
                                               area_t(i,j,k)*fx(i,j,k,URHO)) * volinv
                compu = qp(i,j,k,QRHO)*qp(i,j,k,nqp) - &
                     hdt*(area_t(i+1,j,k)*fx(i+1,j,k,n) - &
                          area_t(i,j,k)*fx(i,j,k,n)) * volinv
                qpo(i,j,k,nqp) = compu/rrnew
#else
                rrnew = qp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qp(i,j,k,QRHO)*qp(i,j,k,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qpo(i,j,k,nqp) = compu/rrnew
#endif

#if AMREX_SPACEDIM == 2
                volinv = ONE/vol(i,j-1,k)

                rrnew = qm(i,j,k,QRHO) - hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,URHO) - &
                                               area_t(i,j-1,k)*fx(i,j-1,k,URHO)) * volinv
                compu = qm(i,j,k,QRHO)*qm(i,j,k,nqp) - &
                     hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,n) - &
                          area_t(i,j-1,k)*fx(i,j-1,k,n)) * volinv
                qmo(i,j,k,nqp) = compu/rrnew
#else
                rrnew = qm(i,j,k,QRHO) - cdtdx*(fx(i+1,j-1,k,URHO) - fx(i,j-1,k,URHO))
                compu = qm(i,j,k,QRHO)*qm(i,j,k,nqp) - cdtdx*(fx(i+1,j-1,k,n) - fx(i,j-1,k,n))
                qmo(i,j,k,nqp) = compu/rrnew
#endif
             end do


             ! add the transverse flux difference in the x-direction to y-states
             ! for the fluid variables


             !----------------------------------------------------------------
             ! qpo state
             !----------------------------------------------------------------

#if AMREX_SPACEDIM == 2
             volinv = ONE/vol(i,j,k)
#endif

             pgp  = q_t(i+1,j,k,GDPRES)
             pgm  = q_t(i  ,j,k,GDPRES)
             ugp  = q_t(i+1,j,k,GDU   )
             ugm  = q_t(i  ,j,k,GDU   )
             gegp = q_t(i+1,j,k,GDGAME)
             gegm = q_t(i  ,j,k,GDGAME)

#ifdef RADIATION
             lambda = qaux(i,j,k,QLAMS:QLAMS+ngroups-1)
             ergp = q_t(i+1,j,k,GDERADS:GDERADS-1+ngroups)
             ergm = q_t(i  ,j,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
             dup = area_t(i+1,j,k)*pgp*ugp - area_t(i,j,k)*pgm*ugm
             du = area_t(i+1,j,k)*ugp-area_t(i,j,k)*ugm
#else
             dup = pgp*ugp - pgm*ugm
             du = ugp-ugm
#endif
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
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
             luge = uav * lamge(:)
             dre = -cdtdx*sum(luge)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdx * uav * f1 * (ergp(g) - ergm(g))
                end do
             else if (fspace_type .eq. 2) then
#if AMREX_SPACEDIM == 2
                divu = (area_t(i+1,j,k)*ugp-area_t(i,j,k)*ugm) * volinv
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = 0.5e0_rt*(1.e0_rt-eddf)
                   der(g) = -hdt * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * divu
                end do
#else
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
                end do
#endif
             else ! mixed frame
                der(:) = cdtdx * luge
             end if
#endif

             ! Convert to conservation form
             rrry = qp(i,j,k,QRHO)
             rury = rrry*qp(i,j,k,QU)
             rvry = rrry*qp(i,j,k,QV)
             rwry = rrry*qp(i,j,k,QW)
             ekenry = HALF*rrry*sum(qp(i,j,k,QU:QW)**2)
             rery = qp(i,j,k,QREINT) + ekenry
#ifdef RADIATION
             err  = qp(i,j,k,qrad:qrad-1+ngroups)
#endif

#if AMREX_SPACEDIM == 2
             ! Add transverse predictor
             rrnewry = rrry - hdt*(area_t(i+1,j,k)*fx(i+1,j,k,URHO) -  &
                                   area_t(i,j,k)*fx(i,j,k,URHO)) * volinv

             ! Note that pressure may be treated specially here, depending on
             ! the geometry.  Our y-interface equation for (rho u) is:
             !
             !  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
             !
             ! in cylindrical coords -- note that the p term is not in
             ! a divergence, so there are no area factors.  For this
             ! geometry, we do not include p in our definition of the
             ! flux in the x-direction, for we need to fix this now.
             runewry = rury - hdt*(area_t(i+1,j,k)*fx(i+1,j,k,UMX)  -  &
                                   area_t(i,j,k)*fx(i,j,k,UMX)) * volinv
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewry = runewry - cdtdx *(pgp-pgm)
             endif
             rvnewry = rvry - hdt*(area_t(i+1,j,k)*fx(i+1,j,k,UMY)  -  &
                                   area_t(i,j,k)*fx(i,j,k,UMY)) * volinv
             rwnewry = rwry - hdt*(area_t(i+1,j,k)*fx(i+1,j,k,UMZ)  -  &
                                   area_t(i,j,k)*fx(i,j,k,UMZ)) * volinv
             renewry = rery - hdt*(area_t(i+1,j,k)*fx(i+1,j,k,UEDEN)-  &
                                  area_t(i,j,k)*fx(i,j,k,UEDEN)) * volinv

#ifdef RADIATION
             runewry = runewry - HALF*hdt*(area_t(i+1,j,k)+area_t(i,j,k))*sum(lamge) * volinv
             renewry = renewry + dre
             ernewr(:) = err(:) - hdt*(area_t(i+1,j,k)*rfx(i+1,j,k,:)-  &
                                       area_t(i,j,k)*rfx(i,j,k,:)) * volinv + der(:)
#endif

#else
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
             qpo(i,j,k,QRHO) = rrnewry
             rhoinv = ONE/rrnewry
             qpo(i,j,k,QU) = runewry*rhoinv
             qpo(i,j,k,QV) = rvnewry*rhoinv
             qpo(i,j,k,QW) = rwnewry*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
             qpo(i,j,k,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) - &
                        hdt*(area_t(i+1,j,k)*fx(i+1,j,k,UEINT)-  &
                             area_t(i,j,k)*fx(i,j,k,UEINT) + pav*du) * volinv
#else
                   qpo(i,j,k,QREINT) = qp(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   ! the divergences here, dup and du, already have area factors
                   pnewry = qp(i,j,k,QPRES) - hdt*(dup + pav*du*(gamc - ONE)) * volinv
#else
                   pnewry = qp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   qpo(i,j,k,QPRES) = max(pnewry, small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du) * volinv - cdtdx*uav*dge
#else
                   qpo(i,j,k,QGAME) = qp(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif
                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,k,QPRES) = qpo(i,j,k,QREINT)*(qpo(i,j,k,QGAME)-ONE)
                   qpo(i,j,k,QPRES) = max(qpo(i,j,k,QPRES),small_pres)
                end if
             else
                qpo(i,j,k,QPRES) = qp(i,j,k,QPRES)
                qpo(i,j,k,QGAME) = qp(i,j,k,QGAME)
             endif


#ifdef RADIATION
             qpo(i,j,k,qrad:qrad-1+ngroups) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qrad-1+ngroups)) + qpo(i,j,k,QREINT)
#endif


             !-------------------------------------------------------------------
             ! qmo state
             !-------------------------------------------------------------------

#if AMREX_SPACEDIM == 2
             volinv = ONE/vol(i,j-1,k)
#endif

             pgp  = q_t(i+1,j-1,k,GDPRES)
             pgm  = q_t(i  ,j-1,k,GDPRES)
             ugp  = q_t(i+1,j-1,k,GDU   )
             ugm  = q_t(i  ,j-1,k,GDU   )
             gegp = q_t(i+1,j-1,k,GDGAME)
             gegm = q_t(i  ,j-1,k,GDGAME)

#ifdef RADIATION
             lambda = qaux(i,j-1,k,QLAMS:QLAMS+ngroups-1)
             ergp = q_t(i+1,j-1,k,GDERADS:GDERADS-1+ngroups)
             ergm = q_t(i  ,j-1,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
             dup = area_t(i+1,j-1,k)*pgp*ugp - area_t(i,j-1,k)*pgm*ugm
             du = area_t(i+1,j-1,k)*ugp-area_t(i,j-1,k)*ugm
#else
             dup = pgp*ugp - pgm*ugm
             du = ugp-ugm
#endif
             pav = HALF*(pgp+pgm)
             uav = HALF*(ugp+ugm)
             geav = HALF*(gegp+gegm)
             dge = gegp-gegm

             ! this is the gas gamma_1
#ifdef RADIATION
             gamc = qaux(i,j-1,k,QGAMCG)
#else
             gamc = qaux(i,j-1,k,QGAMC)
#endif

#ifdef RADIATION
             lamge = lambda(:) * (ergp(:)-ergm(:))
             dmom = - cdtdx*sum(lamge(:))
             luge = uav * lamge(:)
             dre = -cdtdx*sum(luge)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdx * uav * f1 * (ergp(g) - ergm(g))
                end do
             else if (fspace_type .eq. 2) then
#if AMREX_SPACEDIM == 2
                divu = (area_t(i+1,j-1,k)*ugp-area_t(i,j-1,k)*ugm) * volinv
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = 0.5e0_rt*(1.e0_rt-eddf)
                   der(g) = -hdt * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * divu
                end do
#else
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdx * f1 * HALF*(ergp(g)+ergm(g)) * (ugm-ugp)
                end do
#endif
             else ! mixed frame
                der(:) = cdtdx * luge
             end if
#endif

             ! Convert to conservation form
             rrly = qm(i,j,k,QRHO)
             ruly = rrly*qm(i,j,k,QU)
             rvly = rrly*qm(i,j,k,QV)
             rwly = rrly*qm(i,j,k,QW)
             ekenly = HALF*rrly*sum(qm(i,j,k,QU:QW)**2)
             rely = qm(i,j,k,QREINT) + ekenly
#ifdef RADIATION
             erl  = qm(i,j,k,qrad:qrad-1+ngroups)
#endif

#if AMREX_SPACEDIM == 2
             rrnewly = rrly - hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,URHO) -  &
                                   area_t(i,j-1,k)*fx(i,j-1,k,URHO)) * volinv
             runewly = ruly - hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,UMX)  -  &
                                   area_t(i,j-1,k)*fx(i,j-1,k,UMX)) * volinv
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewly = runewly - cdtdx *(pgp-pgm)
             endif
             rvnewly = rvly - hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,UMY)  -  &
                                   area_t(i,j-1,k)*fx(i,j-1,k,UMY)) * volinv
             rwnewly = rwly - hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,UMZ)  -  &
                                   area_t(i,j-1,k)*fx(i,j-1,k,UMZ)) * volinv
             renewly = rely - hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,UEDEN)-  &
                                   area_t(i,j-1,k)*fx(i,j-1,k,UEDEN)) * volinv

#ifdef RADIATION
             runewly = runewly - HALF*hdt*(area_t(i+1,j-1,k)+area_t(i,j-1,k))*sum(lamge) * volinv
             renewly = renewly + dre
             ernewl(:) = erl(:) - hdt*(area_t(i+1,j-1,k)*rfx(i+1,j-1,k,:)-  &
                                       area_t(i,j-1,k)*rfx(i,j-1,k,:)) * volinv + der(:)
#endif

#else
             ! Add transverse predictor
             rrnewly = rrly - cdtdx*(fx(i+1,j-1,k,URHO) - fx(i,j-1,k,URHO))
             runewly = ruly - cdtdx*(fx(i+1,j-1,k,UMX) - fx(i,j-1,k,UMX))
             rvnewly = rvly - cdtdx*(fx(i+1,j-1,k,UMY) - fx(i,j-1,k,UMY))
             rwnewly = rwly - cdtdx*(fx(i+1,j-1,k,UMZ) - fx(i,j-1,k,UMZ))
             renewly = rely - cdtdx*(fx(i+1,j-1,k,UEDEN) - fx(i,j-1,k,UEDEN))
#ifdef RADIATION
             runewly = runewly + dmom
             renewly = renewly + dre
             ernewl  = erl(:) - cdtdx*(rfx(i+1,j-1,k,:) - rfx(i,j-1,k,:)) + der(:)
#endif
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
             qmo(i,j,k,QRHO) = rrnewly
             rhoinv = ONE/rrnewly
             qmo(i,j,k,QU) = runewly*rhoinv
             qmo(i,j,k,QV) = rvnewly*rhoinv
             qmo(i,j,k,QW) = rwnewly*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
             qmo(i,j,k,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qmo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) - &
                        hdt*(area_t(i+1,j-1,k)*fx(i+1,j-1,k,UEINT)-  &
                             area_t(i,j-1,k)*fx(i,j-1,k,UEINT) + pav*du) * volinv
#else
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j-1,k,UEINT) - fx(i,j-1,k,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   pnewly = qm(i,j,k,QPRES) - hdt*(dup + pav*du*(gamc - ONE)) * volinv
#else
                   pnewly = qm(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   qmo(i,j,k,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du) * volinv - cdtdx*uav*dge
#else
                   qmo(i,j,k,QGAME) = qm(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j,k,QPRES) = qmo(i,j,k,QREINT)*(qmo(i,j,k,QGAME)-ONE)
                   qmo(i,j,k,QPRES) = max(qmo(i,j,k,QPRES), small_pres)
                end if
             else
                qmo(i,j,k,QPRES) = qm(i,j,k,QPRES)
                qmo(i,j,k,QGAME) = qm(i,j,k,QGAME)
             endif


#ifdef RADIATION
             qmo(i,j,k,qrad:qrad-1+ngroups) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qrad-1+ngroups)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transx_on_ystates

end module transverse_module
