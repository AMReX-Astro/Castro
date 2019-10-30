module transverse_module

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

contains

  !===========================================================================
  ! transx routines
  !===========================================================================
  subroutine transx_on_ystates(lo, hi, &
                               qym, qym_lo, qym_hi, &
                               qymo, qymo_lo, qymo_hi, &
                               qyp, qyp_lo, qyp_hi, &
                               qypo, qypo_lo, qypo_hi, &
                               qaux, qa_lo, qa_hi, &
                               fx, fx_lo, fx_hi, &
#ifdef RADIATION
                               rfx, rfx_lo, rfx_hi, &
#endif
                               qx, qx_lo, qx_hi, &
#if AMREX_SPACEDIM == 2
                               area1, area1_lo, area1_hi, &
                               vol, vol_lo, vol_hi, &
#endif
                               hdt, cdtdx) bind(C, name="transx_on_ystates")

    ! here, lo and hi are the bounds of the interfaces we are looping over

    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t
#if AMREX_SPACEDIM == 2
    use prob_params_module, only : mom_flux_has_p
#endif

    integer, intent(in) :: qym_lo(3), qym_hi(3)
    integer, intent(in) :: qyp_lo(3), qyp_hi(3)
    integer, intent(in) :: qymo_lo(3), qymo_hi(3)
    integer, intent(in) :: qypo_lo(3), qypo_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
#ifdef RADIATION
    integer, intent(in) :: rfx_lo(3), rfx_hi(3)
#endif
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: lo(3), hi(3)
#if AMREX_SPACEDIM == 2
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
#endif

#ifdef RADIATION
    real(rt) :: rfx(rfx_lo(1):rfx_hi(1),rfx_lo(2):rfx_hi(2),rfx_lo(3):rfx_hi(3),0:ngroups-1)
#endif

    real(rt), intent(in), value :: hdt, cdtdx

    real(rt), intent(in) :: qym(qym_lo(1):qym_hi(1),qym_lo(2):qym_hi(2),qym_lo(3):qym_hi(3),NQ)
    real(rt), intent(in) :: qyp(qyp_lo(1):qyp_hi(1),qyp_lo(2):qyp_hi(2),qyp_lo(3):qyp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

    real(rt), intent(out) :: qymo(qymo_lo(1):qymo_hi(1),qymo_lo(2):qymo_hi(2),qymo_lo(3):qymo_hi(3),NQ)
    real(rt), intent(out) :: qypo(qypo_lo(1):qypo_hi(1),qypo_lo(2):qypo_hi(2),qypo_lo(3):qypo_hi(3),NQ)
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
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
    real(rt)         eddf, f1, ugc
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

    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
#if AMREX_SPACEDIM == 2
                volinv = ONE/vol(i,j,k)

                rrnew = qyp(i,j,k,QRHO) - hdt*(area1(i+1,j,k)*fx(i+1,j,k,URHO) - &
                                               area1(i,j,k)*fx(i,j,k,URHO)) * volinv
                compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nqp) - &
                     hdt*(area1(i+1,j,k)*fx(i+1,j,k,n) - &
                             area1(i,j,k)*fx(i,j,k,n)) * volinv
                qypo(i,j,k,nqp) = compu/rrnew
#else
                rrnew = qyp(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qyp(i,j,k,QRHO)*qyp(i,j,k,nqp) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qypo(i,j,k,nqp) = compu/rrnew
#endif

#if AMREX_SPACEDIM == 2
                volinv = ONE/vol(i,j-1,k)

                rrnew = qym(i,j,k,QRHO) - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,URHO) - &
                                               area1(i,j-1,k)*fx(i,j-1,k,URHO)) * volinv
                compu = qym(i,j,k,QRHO)*qym(i,j,k,nqp) - &
                     hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,n) - &
                          area1(i,j-1,k)*fx(i,j-1,k,n)) * volinv
                qymo(i,j,k,nqp) = compu/rrnew
#else
                rrnew = qym(i,j,k,QRHO) - cdtdx*(fx(i+1,j-1,k,URHO) - fx(i,j-1,k,URHO))
                compu = qym(i,j,k,QRHO)*qym(i,j,k,nqp) - cdtdx*(fx(i+1,j-1,k,n) - fx(i,j-1,k,n))
                qymo(i,j,k,nqp) = compu/rrnew
#endif
             end do
          end do
       end do
    end do

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !----------------------------------------------------------------
             ! qypo state
             !----------------------------------------------------------------

             volinv = ONE/vol(i,j,k)

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

#if AMREX_SPACEDIM == 2
             dup = area1(i+1,j,k)*pgp*ugp - area1(i,j,k)*pgm*ugm
             du = area1(i+1,j,k)*ugp-area1(i,j,k)*ugm
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
             luge = ugc * lamge(:)
             dre = -cdtdx*sum(luge)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
                end do
             else if (fspace_type .eq. 2) then
#if AMREX_SPACEDIM == 2
                divu = (area1(i+1,j,k)*ugp-area1(i,j,k)*ugm) * volinv
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
             rrry = qyp(i,j,k,QRHO)
             rury = rrry*qyp(i,j,k,QU)
             rvry = rrry*qyp(i,j,k,QV)
             rwry = rrry*qyp(i,j,k,QW)
             ekenry = HALF*rrry*sum(qyp(i,j,k,QU:QW)**2)
             rery = qyp(i,j,k,QREINT) + ekenry
#ifdef RADIATION
             err  = qyp(i,j,k,qrad:qradhi)
#endif

#if AMREX_SPACEDIM == 2
             ! Add transverse predictor
             rrnewry = rrry - hdt*(area1(i+1,j,k)*fx(i+1,j,k,URHO) -  &
                                   area1(i,j,k)*fx(i,j,k,URHO)) * volinv

             ! Note that pressure may be treated specially here, depending on
             ! the geometry.  Our y-interface equation for (rho u) is:
             !
             !  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
             !
             ! in cylindrical coords -- note that the p term is not in
             ! a divergence, so there are no area factors.  For this
             ! geometry, we do not include p in our definition of the
             ! flux in the x-direction, for we need to fix this now.
             runewry = rury - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UMX)  -  &
                                   area1(i,j,k)*fx(i,j,k,UMX)) * volinv
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewry = runewry - cdtdx *(pgp-pgm)
             endif
             rvnewry = rvry - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UMY)  -  &
                                   area1(i,j,k)*fx(i,j,k,UMY)) * volinv
             rwnewry = rwry - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UMZ)  -  &
                                   area1(i,j,k)*fx(i,j,k,UMZ)) * volinv
             renewry = rery - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UEDEN)-  &
                                  area1(i,j,k)*fx(i,j,k,UEDEN)) * volinv

#ifdef RADIATION
             runewry = runewry - HALF*hdt*(area1(i+1,j,k)+area1(i,j,k))*sum(lamge) * volinv
             renewry = renewry + dre
             ernewr(:) = err(:) - hdt*(area1(i+1,j,k)*rfx(i+1,j,k,:)-  &
                                       area1(i,j,k)*rfx(i,j,k,:)) * volinv + der(:)
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
#if AMREX_SPACEDIM == 2
                   qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                        hdt*(area1(i+1,j,k)*fx(i+1,j,k,UEINT)-  &
                             area1(i,j,k)*fx(i,j,k,UEINT) + pav*du) * volinv
#else
                   qypo(i,j,k,QREINT) = qyp(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   ! the divergences here, dup and du, already have area factors
                   pnewry = qyp(i,j,k,QPRES) - hdt*(dup + pav*du*(gamc - ONE)) * volinv
#else
                   pnewry = qyp(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   qypo(i,j,k,QPRES) = max(pnewry, small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du) * volinv - cdtdx*uav*dge
#else
                   qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif
                   ! and compute the p edge state from this and (rho e)
                   qypo(i,j,k,QPRES) = qypo(i,j,k,QREINT)*(qypo(i,j,k,QGAME)-ONE)
                   qypo(i,j,k,QPRES) = max(qypo(i,j,k,QPRES),small_pres)
                end if
             else
                qypo(i,j,k,QPRES) = qyp(i,j,k,QPRES)
                qypo(i,j,k,QGAME) = qyp(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qypo, qypo_lo, qypo_hi, i, j, k)

#ifdef RADIATION
             qypo(i,j,k,qrad:qradhi) = ernewr(:)
             qypo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo(i,j,k,QPRES)
             qypo(i,j,k,qreitot) = sum(qypo(i,j,k,qrad:qradhi)) + qypo(i,j,k,QREINT)
#endif


             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------
             volinv = ONE/vol(i,j-1,k)

             pgp  = qx(i+1,j-1,k,GDPRES)
             pgm  = qx(i  ,j-1,k,GDPRES)
             ugp  = qx(i+1,j-1,k,GDU   )
             ugm  = qx(i  ,j-1,k,GDU   )
             gegp = qx(i+1,j-1,k,GDGAME)
             gegm = qx(i  ,j-1,k,GDGAME)

#ifdef RADIATION
             lambda = qaux(i,j-1,k,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qx(i+1,j-1,k,GDERADS:GDERADS-1+ngroups)
             ergm = qx(i  ,j-1,k,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
             dup = area1(i+1,j-1,k)*pgp*ugp - area1(i,j-1,k)*pgm*ugm
             du = area1(i+1,j-1,k)*ugp-area1(i,j-1,k)*ugm
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
             luge = ugc * lamge(:)
             dre = -cdtdx*sum(luge)

             if (fspace_type .eq. 1 .and. comoving) then
                do g=0, ngroups-1
                   eddf = Edd_factor(lambda(g))
                   f1 = HALF*(ONE-eddf)
                   der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
                end do
             else if (fspace_type .eq. 2) then
#if AMREX_SPACEDIM == 2
                divu = (area1(i+1,j-1,k)*ugp-area1(i,j-1,k)*ugm) * volinv
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
             rrly = qym(i,j,k,QRHO)
             ruly = rrly*qym(i,j,k,QU)
             rvly = rrly*qym(i,j,k,QV)
             rwly = rrly*qym(i,j,k,QW)
             ekenly = HALF*rrly*sum(qym(i,j,k,QU:QW)**2)
             rely = qym(i,j,k,QREINT) + ekenly
#ifdef RADIATION
             erl  = qym(i,j,k,qrad:qradhi)
#endif

#if AMREX_SPACEDIM == 2
             rrnewly = rrly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,URHO) -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,URHO)) * volinv
             runewly = ruly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UMX)  -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UMX)) * volinv
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewly = runewly - cdtdx *(pgp-pgm)
             endif
             rvnewly = rvly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UMY)  -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UMY)) * volinv
             rwnewly = rwly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UMZ)  -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UMZ)) * volinv
             renewly = rely - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UEDEN)-  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UEDEN)) * volinv

#ifdef RADIATION
             runewly = runewly - HALF*hdt*(area1(i+1,j-1,k)+area1(i,j-1,k))*sum(lamge) * volinv
             renewly = renewly + dre
             ernewl(:) = erl(:) - hdt*(area1(i+1,j-1,k)*rfx(i+1,j-1,k,:)-  &
                                       area1(i,j-1,k)*rfx(i,j-1,k,:)) * volinv + der(:)
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
#if AMREX_SPACEDIM == 2
                   qymo(i,j,k,QREINT) = qym(i,j,k,QREINT) - &
                        hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UEINT)-  &
                             area1(i,j-1,k)*fx(i,j-1,k,UEINT) + pav*du) * volinv
#else
                   qymo(i,j,k,QREINT) = qym(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j-1,k,UEINT) - fx(i,j-1,k,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   pnewly = qym(i,j,k,QPRES) - hdt*(dup + pav*du*(gamc - ONE)) * volinv
#else
                   pnewly = qym(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   qymo(i,j,k,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   qymo(i,j,k,QGAME) = qym(i,j,k,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du) * volinv - cdtdx*uav*dge
#else
                   qymo(i,j,k,QGAME) = qym(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif

                   ! and compute the p edge state from this and (rho e)
                   qymo(i,j,k,QPRES) = qymo(i,j,k,QREINT)*(qymo(i,j,k,QGAME)-ONE)
                   qymo(i,j,k,QPRES) = max(qymo(i,j,k,QPRES), small_pres)
                end if
             else
                qymo(i,j,k,QPRES) = qym(i,j,k,QPRES)
                qymo(i,j,k,QGAME) = qym(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qymo, qymo_lo, qymo_hi, i, j, k)

#ifdef RADIATION
             qymo(i,j,k,qrad:qradhi) = ernewl(:)
             qymo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo(i,j,k,QPRES)
             qymo(i,j,k,qreitot) = sum(qymo(i,j,k,qrad:qradhi)) + qymo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transx_on_ystates


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
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
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
  use fluxlimiter_module, only : Edd_factor
#endif
  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_t
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

             call reset_edge_state_thermo(qzpo, qzpo_lo, qzpo_hi, i, j, k)

#ifdef RADIATION
             qzpo(i,j,k,qrad:qradhi) = ernewr(:)
             qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
             qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qradhi)) + qzpo(i,j,k,QREINT)
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
             erl  = qzm(i,j,k,qrad:qradhi)
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
             qzmo(i,j,k,qrad:qradhi) = ernewl(:)
             qzmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k,QPRES)
             qzmo(i,j,k,qreitot) = sum(qzmo(i,j,k,qrad:qradhi)) + qzmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transx_on_zstates


  !===========================================================================
  ! transy routines
  !===========================================================================
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
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


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

             call reset_edge_state_thermo(qxpo, qxpo_lo, qxpo_hi, i, j, k)

#ifdef RADIATION
             qxpo(i,j,k,qrad:qradhi) = ernewr(:)
             qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo(i,j,k,QPRES)
             qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qradhi)) + qxpo(i,j,k,QREINT)
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
             erl  = qxm(i,j,k,qrad:qradhi)
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
             qxmo(i,j,k,qrad:qradhi) = ernewl(:)
             qxmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo(i,j,k,QPRES)
             qxmo(i,j,k,qreitot) = sum(qxmo(i,j,k,qrad:qradhi)) + qxmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transy_on_xstates

  subroutine transy_on_zstates(lo, hi, &
                               qzm, qzm_lo, qzm_hi, &
                               qzmo, qzmo_lo, qzmo_hi, &
                               qzp, qzp_lo, qzp_hi, &
                               qzpo, qzpo_lo, qzpo_hi, &
                               qaux, qa_lo, qa_hi, &
                               fy, fy_lo, fy_hi, &
#ifdef RADIATION
                               rfy, rfy_lo, rfy_hi, &
#endif
                               qy, qy_lo, qy_hi, &
                               cdtdy) bind(C, name="transy_on_zstates")

    ! here, lo and hi are the bounds of edges we are looping over


    use amrex_constants_module, only : ZERO, ONE, HALF

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t

    integer, intent(in) :: qzm_lo(3), qzm_hi(3)
    integer, intent(in) :: qzp_lo(3), qzp_hi(3)
    integer, intent(in) :: qzmo_lo(3), qzmo_hi(3)
    integer, intent(in) :: qzpo_lo(3), qzpo_hi(3)
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

    real(rt), intent(in) :: qzm(qzm_lo(1):qzm_hi(1),qzm_lo(2):qzm_hi(2),qzm_lo(3):qzm_hi(3),NQ)
    real(rt), intent(in) :: qzp(qzp_lo(1):qzp_hi(1),qzp_lo(2):qzp_hi(2),qzp_lo(3):qzp_hi(3),NQ)
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

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
                rrnew = qzp(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = qzp(i,j,k,QRHO)*qzp(i,j,k,nqp) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                qzpo(i,j,k,nqp) = compu/rrnew

                rrnew = qzm(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k-1,URHO) - fy(i,j,k-1,URHO))
                compu = qzm(i,j,k,QRHO)*qzm(i,j,k,nqp) - cdtdy*(fy(i,j+1,k-1,n) - fy(i,j,k-1,n))
                qzmo(i,j,k,nqp) = compu/rrnew
             end do
          end do
       end do
    end do


    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to z-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             !-------------------------------------------------------------------
             ! qzpo states
             !-------------------------------------------------------------------

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

             call reset_edge_state_thermo(qzpo, qzpo_lo, qzpo_hi, i, j, k)

#ifdef RADIATION
             qzpo(i,j,k,qrad:qradhi) = ernewr(:)
             qzpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo(i,j,k,QPRES)
             qzpo(i,j,k,qreitot) = sum(qzpo(i,j,k,qrad:qradhi)) + qzpo(i,j,k,QREINT)
#endif

             !-------------------------------------------------------------------
             ! qzmo state
             !-------------------------------------------------------------------

             pgp  = qy(i,j+1,k-1,GDPRES)
             pgm  = qy(i,j  ,k-1,GDPRES)
             ugp  = qy(i,j+1,k-1,GDV   )
             ugm  = qy(i,j  ,k-1,GDV   )
             gegp = qy(i,j+1,k-1,GDGAME)
             gegm = qy(i,j  ,k-1,GDGAME)
#ifdef RADIATION
             lambda(:) = qaux(i,j,k-1,QLAMS:QLAMS+ngroups-1)
             ugc = HALF*(ugp+ugm)
             ergp = qy(i,j+1,k-1,GDERADS:GDERADS-1+ngroups)
             ergm = qy(i,j  ,k-1,GDERADS:GDERADS-1+ngroups)
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
             rrlz = qzm(i,j,k,QRHO)
             rulz = rrlz*qzm(i,j,k,QU)
             rvlz = rrlz*qzm(i,j,k,QV)
             rwlz = rrlz*qzm(i,j,k,QW)
             ekenlz = HALF*rrlz*sum(qzm(i,j,k,QU:QW)**2)
             relz = qzm(i,j,k,QREINT) + ekenlz
#ifdef RADIATION
             erl  = qzm(i,j,k,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewlz = rrlz - cdtdy*(fy(i,j+1,k-1,URHO) - fy(i,j,k-1,URHO))
             runewlz = rulz - cdtdy*(fy(i,j+1,k-1,UMX) - fy(i,j,k-1,UMX))
             rvnewlz = rvlz - cdtdy*(fy(i,j+1,k-1,UMY) - fy(i,j,k-1,UMY))
             rwnewlz = rwlz - cdtdy*(fy(i,j+1,k-1,UMZ) - fy(i,j,k-1,UMZ))
             renewlz = relz - cdtdy*(fy(i,j+1,k-1,UEDEN) - fy(i,j,k-1,UEDEN))
#ifdef RADIATION
             rvnewlz = rvnewlz + dmom
             renewlz = renewlz + dre
             ernewl  = erl(:) - cdtdy*(rfy(i,j+1,k-1,:) - rfy(i,j,k-1,:)) + der
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
                        cdtdy*(fy(i,j+1,k-1,UEINT) - fy(i,j,k-1,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qzmo(i,j,k,QPRES) = max(pnewlz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzmo(i,j,k,QGAME) = qzm(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

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
             qzmo(i,j,k,qrad:qradhi) = ernewl(:)
             qzmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo(i,j,k,QPRES)
             qzmo(i,j,k,qreitot) = sum(qzmo(i,j,k,qrad:qradhi)) + qzmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transy_on_zstates


#if AMREX_SPACEDIM == 3
  !===========================================================================
  ! transz routines
  !===========================================================================

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
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


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

             call reset_edge_state_thermo(qxpo, qxpo_lo, qxpo_hi, i, j, k)

#ifdef RADIATION
             qxpo(i,j,k,qrad:qradhi) = ernewrx(:)
             qxpo(i,j,k,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo(i,j,k,QPRES)
             qxpo(i,j,k,qreitot) = sum(qxpo(i,j,k,qrad:qradhi)) + qxpo(i,j,k,QREINT)
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
             erlx = qxm(i,j,k,qrad:qradhi)
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
             qxmo(i,j,k,qrad:qradhi) = ernewlx(:)
             qxmo(i,j,k,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo(i,j,k,QPRES)
             qxmo(i,j,k,qreitot) = sum(qxmo(i,j,k,qrad:qradhi)) + qxmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transz_on_xstates


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
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


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

             call reset_edge_state_thermo(qypo, qypo_lo, qypo_hi, i, j, k)

#ifdef RADIATION
             qypo(i,j,k,qrad:qradhi) = ernewry(:)
             qypo(i,j,k,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo(i,j,k,QPRES)
             qypo(i,j,k,qreitot) = sum(qypo(i,j,k,qrad:qradhi)) + qypo(i,j,k,QREINT)
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
             erly = qym(i,j,k,qrad:qradhi)
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
             qymo(i,j,k,qrad:qradhi) = ernewly(:)
             qymo(i,j,k,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo(i,j,k,QPRES)
             qymo(i,j,k,qreitot) = sum(qymo(i,j,k,qrad:qradhi)) + qymo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transz_on_ystates

  !===========================================================================
  ! transyz
  !===========================================================================
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
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


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
             qpo(i,j,k,qrad:qradhi) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
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
             erl = qm(i,j,k,qrad:qradhi)
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
             qmo(i,j,k,qrad:qradhi) = ernewl(:)
             qmo(i,j,k,qptot) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qradhi)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transyz

  !===========================================================================
  ! transxz
  !===========================================================================
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
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


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
             qpo(i,j,k,qrad:qradhi) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
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
             erl = qm(i,j,k,qrad:qradhi)
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
             qmo(i,j,k,qrad:qradhi) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qradhi)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxz

  !===========================================================================
  ! transxy
  !===========================================================================
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
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   qrad, qradhi, qptot, qreitot, &
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
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_t


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
             qpo(i,j,k,qrad:qradhi) = ernewr(:)
             qpo(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo(i,j,k,QPRES)
             qpo(i,j,k,qreitot) = sum(qpo(i,j,k,qrad:qradhi)) + qpo(i,j,k,QREINT)
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
             erl = qm(i,j,k,qrad:qradhi)
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
             qmo(i,j,k,qrad:qradhi) = ernewl(:)
             qmo(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo(i,j,k,QPRES)
             qmo(i,j,k,qreitot) = sum(qmo(i,j,k,qrad:qradhi)) + qmo(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxy
#endif

  subroutine reset_edge_state_thermo(qedge, qd_lo, qd_hi, ii, jj, kk)

  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQ, NVAR, NQAUX, QRHO, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 small_pres, small_temp, &
                                 ppm_predict_gammae, &
                                 transverse_use_eos, transverse_reset_rhoe

  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer, intent(in) :: ii, jj, kk
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    real(rt)        , intent(inout) :: qedge(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    logical :: reset
    type (eos_t) :: eos_state

    !$gpu

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
