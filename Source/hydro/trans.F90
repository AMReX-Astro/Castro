module transverse_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg

  implicit none

contains

  !===========================================================================
  ! transx routines
  !===========================================================================
  subroutine transx_on_ystates(lo, hi, &
                               qym_core, qycm_lo, qycm_hi, &
                               qymo_core, qycmo_lo, qycmo_hi, &
                               qym_pass, qypm_lo, qypm_hi, &
                               qymo_pass, qypmo_lo, qypmo_hi, &
#ifdef RADIATION
                               qym_rad, qyrm_lo, qyrm_hi, &
                               qymo_rad, qyrmo_lo, qyrmo_hi, &
#endif
                               qyp_core, qycp_lo, qycp_hi, &
                               qypo_core, qycpo_lo, qycpo_hi, &
                               qyp_pass, qypp_lo, qypp_hi, &
                               qypo_pass, qyppo_lo, qyppo_hi, &
#ifdef RADIATION
                               qyp_rad, qyrp_lo, qyrp_hi, &
                               qypo_rad, qyrpo_lo, qyrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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

    integer, intent(in) :: qycm_lo(3), qycm_hi(3)
    integer, intent(in) :: qypm_lo(3), qypm_hi(3)
    integer, intent(in) :: qycmo_lo(3), qycmo_hi(3)
    integer, intent(in) :: qypmo_lo(3), qypmo_hi(3)
    integer, intent(in) :: qycp_lo(3), qycp_hi(3)
    integer, intent(in) :: qypp_lo(3), qypp_hi(3)
    integer, intent(in) :: qycpo_lo(3), qycpo_hi(3)
    integer, intent(in) :: qyppo_lo(3), qyppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qyrm_lo(3), qyrm_hi(3)
    integer, intent(in) :: qyrmo_lo(3), qyrmo_hi(3)
    integer, intent(in) :: qyrp_lo(3), qyrp_hi(3)
    integer, intent(in) :: qyrpo_lo(3), qyrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qym_core(qycm_lo(1):qycm_hi(1),qycm_lo(2):qycm_hi(2),qycm_lo(3):qycm_hi(3),NQC)
    real(rt), intent(in) :: qyp_core(qycp_lo(1):qycp_hi(1),qycp_lo(2):qycp_hi(2),qycp_lo(3):qycp_hi(3),NQC)
    real(rt), intent(in) :: qym_pass(qypm_lo(1):qypm_hi(1),qypm_lo(2):qypm_hi(2),qypm_lo(3):qypm_hi(3),NQP)
    real(rt), intent(in) :: qyp_pass(qypp_lo(1):qypp_hi(1),qypp_lo(2):qypp_hi(2),qypp_lo(3):qypp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qym_rad(qyrm_lo(1):qyrm_hi(1),qyrm_lo(2):qyrm_hi(2),qyrm_lo(3):qyrm_hi(3),NQR)
    real(rt), intent(in) ::  qyp_rad(qyrp_lo(1):qyrp_hi(1),qyrp_lo(2):qyrp_hi(2),qyrp_lo(3):qyrp_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

    real(rt), intent(out) :: qymo_core(qycmo_lo(1):qycmo_hi(1),qycmo_lo(2):qycmo_hi(2),qycmo_lo(3):qycmo_hi(3),NQC)
    real(rt), intent(out) :: qypo_core(qycpo_lo(1):qycpo_hi(1),qycpo_lo(2):qycpo_hi(2),qycpo_lo(3):qycpo_hi(3),NQC)
    real(rt), intent(out) :: qymo_pass(qypmo_lo(1):qypmo_hi(1),qypmo_lo(2):qypmo_hi(2),qypmo_lo(3):qypmo_hi(3),NQP)
    real(rt), intent(out) :: qypo_pass(qyppo_lo(1):qyppo_hi(1),qyppo_lo(2):qyppo_hi(2),qyppo_lo(3):qyppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qymo_rad(qyrmo_lo(1):qyrmo_hi(1),qyrmo_lo(2):qyrmo_hi(2),qyrmo_lo(3):qyrmo_hi(3),NQR)
    real(rt), intent(out) ::  qypo_rad(qyrpo_lo(1):qyrpo_hi(1),qyrpo_lo(2):qyrpo_hi(2),qyrpo_lo(3):qyrpo_hi(3),NQR)
#endif
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
#endif

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
#if AMREX_SPACEDIM == 2
                rrnew = qyp_core(i,j,k,QRHO) - hdt*(area1(i+1,j,k)*fx(i+1,j,k,URHO) - &
                                                    area1(i,j,k)*fx(i,j,k,URHO))/vol(i,j,k)
                compu = qyp_core(i,j,k,QRHO)*qyp_pass(i,j,k,nqs) - &
                     hdt*(area1(i+1,j,k)*fx(i+1,j,k,n) - &
                             area1(i,j,k)*fx(i,j,k,n))/vol(i,j,k)
                qypo_pass(i,j,k,nqs) = compu/rrnew
#else
                rrnew = qyp_core(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qyp_core(i,j,k,QRHO)*qyp_pass(i,j,k,nqs) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qypo_pass(i,j,k,nqs) = compu/rrnew
#endif

#if AMREX_SPACEDIM == 2
                rrnew = qym_core(i,j,k,QRHO) - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,URHO) - &
                                                    area1(i,j-1,k)*fx(i,j-1,k,URHO))/vol(i,j-1,k)
                compu = qym_core(i,j,k,QRHO)*qym_pass(i,j,k,nqs) - &
                     hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,n) - &
                          area1(i,j-1,k)*fx(i,j-1,k,n))/vol(i,j-1,k)
                qymo_pass(i,j,k,nqs) = compu/rrnew
#else
                rrnew = qym_core(i,j,k,QRHO) - cdtdx*(fx(i+1,j-1,k,URHO) - fx(i,j-1,k,URHO))
                compu = qym_core(i,j,k,QRHO)*qym_pass(i,j,k,nqs) - cdtdx*(fx(i+1,j-1,k,n) - fx(i,j-1,k,n))
                qymo_pass(i,j,k,nqs) = compu/rrnew
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
                divu = (area1(i+1,j,k)*ugp-area1(i,j,k)*ugm)/vol(i,j,k)
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
             rrry = qyp_core(i,j,k,QRHO)
             rury = rrry*qyp_core(i,j,k,QU)
             rvry = rrry*qyp_core(i,j,k,QV)
             rwry = rrry*qyp_core(i,j,k,QW)
             ekenry = HALF*rrry*sum(qyp_core(i,j,k,QU:QW)**2)
             rery = qyp_core(i,j,k,QREINT) + ekenry
#ifdef RADIATION
             err  = qyp_rad(i,j,k,qrad:qradhi)
#endif

#if AMREX_SPACEDIM == 2
             ! Add transverse predictor
             rrnewry = rrry - hdt*(area1(i+1,j,k)*fx(i+1,j,k,URHO) -  &
                                   area1(i,j,k)*fx(i,j,k,URHO))/vol(i,j,k)

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
                                   area1(i,j,k)*fx(i,j,k,UMX))/vol(i,j,k)
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewry = runewry - cdtdx *(pgp-pgm)
             endif
             rvnewry = rvry - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UMY)  -  &
                                   area1(i,j,k)*fx(i,j,k,UMY))/vol(i,j,k)
             rwnewry = rwry - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UMZ)  -  &
                                   area1(i,j,k)*fx(i,j,k,UMZ))/vol(i,j,k)
             renewry = rery - hdt*(area1(i+1,j,k)*fx(i+1,j,k,UEDEN)-  &
                                  area1(i,j,k)*fx(i,j,k,UEDEN))/vol(i,j,k)

#ifdef RADIATION
             runewry = runewry - HALF*hdt*(area1(i+1,j,k)+area1(i,j,k))*sum(lamge)/vol(i,j,k)
             renewry = renewry + dre
             ernewr(:) = err(:) - hdt*(area1(i+1,j,k)*rfx(i+1,j,k,:)-  &
                                       area1(i,j,k)*rfx(i,j,k,:))/vol(i,j,k) + der(:)
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
             qypo_core(i,j,k,QRHO) = rrnewry
             rhoinv = ONE/rrnewry
             qypo_core(i,j,k,QU) = runewry*rhoinv
             qypo_core(i,j,k,QV) = rvnewry*rhoinv
             qypo_core(i,j,k,QW) = rwnewry*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
             qypo_core(i,j,k,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qypo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                   qypo_core(i,j,k,QREINT) = qyp_core(i,j,k,QREINT) - &
                        hdt*(area1(i+1,j,k)*fx(i+1,j,k,UEINT)-  &
                             area1(i,j,k)*fx(i,j,k,UEINT) + pav*du)/vol(i,j,k)
#else
                   qypo_core(i,j,k,QREINT) = qyp_core(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   ! the divergences here, dup and du, already have area factors
                   pnewry = qyp_core(i,j,k,QPRES) - hdt*(dup + pav*du*(gamc - ONE))/vol(i,j,k)
#else
                   pnewry = qyp_core(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   qypo_core(i,j,k,QPRES) = max(pnewry, small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   qypo_core(i,j,k,QGAME) = qyp_core(i,j,k,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du)/vol(i,j,k) - cdtdx*uav*dge
#else
                   qypo_core(i,j,k,QGAME) = qyp_core(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif
                   ! and compute the p edge state from this and (rho e)
                   qypo_core(i,j,k,QPRES) = qypo_core(i,j,k,QREINT)*(qypo_core(i,j,k,QGAME)-ONE)
                   qypo_core(i,j,k,QPRES) = max(qypo_core(i,j,k,QPRES),small_pres)
                end if
             else
                qypo_core(i,j,k,QPRES) = qyp_core(i,j,k,QPRES)
                qypo_core(i,j,k,QGAME) = qyp_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qypo_core, qycpo_lo, qycpo_hi, &
                                          qypo_pass, qyppo_lo, qyppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qypo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qypo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qypo_core(i,j,k,QPRES)
             qypo_rad(i,j,k,qreitot) = sum(qypo_rad(i,j,k,qrad:qradhi)) + qypo_core(i,j,k,QREINT)
#endif


             !-------------------------------------------------------------------
             ! qymo state
             !-------------------------------------------------------------------
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
                divu = (area1(i+1,j-1,k)*ugp-area1(i,j-1,k)*ugm)/vol(i,j-1,k)
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
             rrly = qym_core(i,j,k,QRHO)
             ruly = rrly*qym_core(i,j,k,QU)
             rvly = rrly*qym_core(i,j,k,QV)
             rwly = rrly*qym_core(i,j,k,QW)
             ekenly = HALF*rrly*sum(qym_core(i,j,k,QU:QW)**2)
             rely = qym_core(i,j,k,QREINT) + ekenly
#ifdef RADIATION
             erl  = qym_rad(i,j,k,qrad:qradhi)
#endif

#if AMREX_SPACEDIM == 2
             rrnewly = rrly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,URHO) -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,URHO))/vol(i,j-1,k)
             runewly = ruly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UMX)  -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UMX))/vol(i,j-1,k)
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewly = runewly - cdtdx *(pgp-pgm)
             endif
             rvnewly = rvly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UMY)  -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UMY))/vol(i,j-1,k)
             rwnewly = rwly - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UMZ)  -  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UMZ))/vol(i,j-1,k)
             renewly = rely - hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UEDEN)-  &
                                   area1(i,j-1,k)*fx(i,j-1,k,UEDEN))/vol(i,j-1,k)

#ifdef RADIATION
             runewly = runewly - HALF*hdt*(area1(i+1,j-1,k)+area1(i,j-1,k))*sum(lamge)/vol(i,j-1,k)
             renewly = renewly + dre
             ernewl(:) = erl(:) - hdt*(area1(i+1,j-1,k)*rfx(i+1,j-1,k,:)-  &
                                       area1(i,j-1,k)*rfx(i,j-1,k,:))/vol(i,j-1,k) + der(:)
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
             qymo_core(i,j,k,QRHO) = rrnewly
             rhoinv = ONE/rrnewly
             qymo_core(i,j,k,QU) = runewly*rhoinv
             qymo_core(i,j,k,QV) = rvnewly*rhoinv
             qymo_core(i,j,k,QW) = rwnewly*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
             qymo_core(i,j,k,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qymo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                   qymo_core(i,j,k,QREINT) = qym_core(i,j,k,QREINT) - &
                        hdt*(area1(i+1,j-1,k)*fx(i+1,j-1,k,UEINT)-  &
                             area1(i,j-1,k)*fx(i,j-1,k,UEINT) + pav*du)/vol(i,j-1,k)
#else
                   qymo_core(i,j,k,QREINT) = qym_core(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j-1,k,UEINT) - fx(i,j-1,k,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   pnewly = qym_core(i,j,k,QPRES) - hdt*(dup + pav*du*(gamc - ONE))/vol(i,j-1,k)
#else
                   pnewly = qym_core(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   qymo_core(i,j,k,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   qymo_core(i,j,k,QGAME) = qym_core(i,j,k,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du)/vol(i,j-1,k) - cdtdx*uav*dge
#else
                   qymo_core(i,j,k,QGAME) = qym_core(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif

                   ! and compute the p edge state from this and (rho e)
                   qymo_core(i,j,k,QPRES) = qymo_core(i,j,k,QREINT)*(qymo_core(i,j,k,QGAME)-ONE)
                   qymo_core(i,j,k,QPRES) = max(qymo_core(i,j,k,QPRES), small_pres)
                end if
             else
                qymo_core(i,j,k,QPRES) = qym_core(i,j,k,QPRES)
                qymo_core(i,j,k,QGAME) = qym_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qymo_core, qycmo_lo, qycmo_hi, &
                                          qymo_pass, qypmo_lo, qypmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qymo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qymo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qymo_core(i,j,k,QPRES)
             qymo_rad(i,j,k,qreitot) = sum(qymo_rad(i,j,k,qrad:qradhi)) + qymo_pass(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transx_on_ystates


  subroutine transx_on_zstates(lo, hi, &
                               qzm_core, qzcm_lo, qzcm_hi, &
                               qzmo_core, qzcmo_lo, qzcmo_hi, &
                               qzm_pass, qzpm_lo, qzpm_hi, &
                               qzmo_pass, qzpmo_lo, qzpmo_hi, &
#ifdef RADIATION
                               qzm_rad, qzrm_lo, qzrm_hi, &
                               qzmo_rad, qzrmo_lo, qzrmo_hi, &
#endif
                               qzp_core, qzcp_lo, qzcp_hi, &
                               qzpo_core, qzcpo_lo, qzcpo_hi, &
                               qzp_pass, qzpp_lo, qzpp_hi, &
                               qzpo_pass, qzppo_lo, qzppo_hi, &
#ifdef RADIATION
                               qzp_rad, qzrp_lo, qzrp_hi, &
                               qzpo_rad, qzrpo_lo, qzrpo_hi, &
#endif
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
  use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 QC, QGAMC, &
#ifdef RADIATION
                                 NQR, qrad, qradhi, qptot, qreitot, &
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

    integer, intent(in) :: qzcm_lo(3), qzcm_hi(3)
    integer, intent(in) :: qzcp_lo(3), qzcp_hi(3)
    integer, intent(in) :: qzcmo_lo(3), qzcmo_hi(3)
    integer, intent(in) :: qzcpo_lo(3), qzcpo_hi(3)
    integer, intent(in) :: qzpm_lo(3), qzpm_hi(3)
    integer, intent(in) :: qzpp_lo(3), qzpp_hi(3)
    integer, intent(in) :: qzpmo_lo(3), qzpmo_hi(3)
    integer, intent(in) :: qzppo_lo(3), qzppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qzrm_lo(3), qzrm_hi(3)
    integer, intent(in) :: qzrmo_lo(3), qzrmo_hi(3)
    integer, intent(in) :: qzrp_lo(3), qzrp_hi(3)
    integer, intent(in) :: qzrpo_lo(3), qzrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qzm_core(qzcm_lo(1):qzcm_hi(1),qzcm_lo(2):qzcm_hi(2),qzcm_lo(3):qzcm_hi(3),NQC)
    real(rt), intent(in) :: qzp_core(qzcp_lo(1):qzcp_hi(1),qzcp_lo(2):qzcp_hi(2),qzcp_lo(3):qzcp_hi(3),NQC)
    real(rt), intent(in) :: qzm_pass(qzpm_lo(1):qzpm_hi(1),qzpm_lo(2):qzpm_hi(2),qzpm_lo(3):qzpm_hi(3),NQP)
    real(rt), intent(in) :: qzp_pass(qzpp_lo(1):qzpp_hi(1),qzpp_lo(2):qzpp_hi(2),qzpp_lo(3):qzpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qzm_rad(qzrm_lo(1):qzrm_hi(1),qzrm_lo(2):qzrm_hi(2),qzrm_lo(3):qzrm_hi(3),NQR)
    real(rt), intent(in) ::  qzp_rad(qzrp_lo(1):qzrp_hi(1),qzrp_lo(2):qzrp_hi(2),qzrp_lo(3):qzrp_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

    real(rt), intent(out) :: qzmo_core(qzcmo_lo(1):qzcmo_hi(1),qzcmo_lo(2):qzcmo_hi(2),qzcmo_lo(3):qzcmo_hi(3),NQC)
    real(rt), intent(out) :: qzpo_core(qzcpo_lo(1):qzcpo_hi(1),qzcpo_lo(2):qzcpo_hi(2),qzcpo_lo(3):qzcpo_hi(3),NQC)
    real(rt), intent(out) :: qzmo_pass(qzpmo_lo(1):qzpmo_hi(1),qzpmo_lo(2):qzpmo_hi(2),qzpmo_lo(3):qzpmo_hi(3),NQP)
    real(rt), intent(out) :: qzpo_pass(qzppo_lo(1):qzppo_hi(1),qzppo_lo(2):qzppo_hi(2),qzppo_lo(3):qzppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qzmo_rad(qzrmo_lo(1):qzrmo_hi(1),qzrmo_lo(2):qzrmo_hi(2),qzrmo_lo(3):qzrmo_hi(3),NQR)
    real(rt), intent(out) ::  qzpo_rad(qzrpo_lo(1):qzrpo_hi(1),qzrpo_lo(2):qzrpo_hi(2),qzrpo_lo(3):qzrpo_hi(3),NQR)
#endif

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qzp_core(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k,URHO) - fx(i,j,k,URHO))
                compu = qzp_core(i,j,k,QRHO)*qzp_pass(i,j,k,nqs) - cdtdx*(fx(i+1,j,k,n) - fx(i,j,k,n))
                qzpo_pass(i,j,k,nqs) = compu/rrnew

                rrnew = qzm_core(i,j,k,QRHO) - cdtdx*(fx(i+1,j,k-1,URHO) - fx(i,j,k-1,URHO))
                compu = qzm_core(i,j,k,QRHO)*qzm_pass(i,j,k,nqs) - cdtdx*(fx(i+1,j,k-1,n) - fx(i,j,k-1,n))
                qzmo_pass(i,j,k,nqs) = compu/rrnew
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
             rrrz = qzp_core(i,j,k,QRHO)
             rurz = rrrz*qzp_core(i,j,k,QU)
             rvrz = rrrz*qzp_core(i,j,k,QV)
             rwrz = rrrz*qzp_core(i,j,k,QW)
             ekenrz = HALF*rrrz*sum(qzp_core(i,j,k,QU:QW)**2)
             rerz = qzp_core(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
             err  = qzp_rad(i,j,k,qrad:qradhi)
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
             qzpo_core(i,j,k,QRHO) = rrnewrz
             rhoinv = ONE/rrnewrz
             qzpo_core(i,j,k,QU) = runewrz*rhoinv
             qzpo_core(i,j,k,QV) = rvnewrz*rhoinv
             qzpo_core(i,j,k,QW) = rwnewrz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
             qzpo_core(i,j,k,QREINT) = renewrz - rhoekenrz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzpo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qzpo_core(i,j,k,QREINT) = qzp_core(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k,UEINT) - fx(i,j,k,UEINT) + pav*du)
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrz = qzp_core(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qzpo_core(i,j,k,QPRES) = max(pnewrz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzpo_core(i,j,k,QGAME) = qzp_core(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzpo_core(i,j,k,QPRES) = qzpo_core(i,j,k,QREINT)*(qzpo_core(i,j,k,QGAME)-ONE)
                   qzpo_core(i,j,k,QPRES) = max(qzpo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qzpo_core(i,j,k,QPRES) = qzp_core(i,j,k,QPRES)
                qzpo_core(i,j,k,QGAME) = qzp_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzpo_core, qzcpo_lo, qzcpo_hi, &
                                          qzpo_pass, qzppo_lo, qzppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qzpo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qzpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo_core(i,j,k,QPRES)
             qzpo_rad(i,j,k,qreitot) = sum(qzpo_rad(i,j,k,qrad:qradhi)) + qzpo_core(i,j,k,QREINT)
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
             rrlz = qzm_core(i,j,k,QRHO)
             rulz = rrlz*qzm_core(i,j,k,QU)
             rvlz = rrlz*qzm_core(i,j,k,QV)
             rwlz = rrlz*qzm_core(i,j,k,QW)
             ekenlz = HALF*rrlz*sum(qzm_core(i,j,k,QU:QW)**2)
             relz = qzm_core(i,j,k,QREINT) + ekenlz
#ifdef RADIATION
             erl  = qzm_rad(i,j,k,qrad:qradhi)
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
             qzmo_core(i,j,k,QRHO) = rrnewlz
             rhoinv = ONE/rrnewlz
             qzmo_core(i,j,k,QU) = runewlz*rhoinv
             qzmo_core(i,j,k,QV) = rvnewlz*rhoinv
             qzmo_core(i,j,k,QW) = rwnewlz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
             qzmo_core(i,j,k,QREINT) = renewlz - rhoekenlz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzmo_core(i,j,k,QREINT) = qzm_core(i,j,k,QREINT) - &
                        cdtdx*(fx(i+1,j,k-1,UEINT) - fx(i,j,k-1,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm_core(i,j,k,QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
                   qzmo_core(i,j,k,QPRES) = max(pnewlz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzmo_core(i,j,k,QGAME) = qzm_core(i,j,k,QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzmo_core(i,j,k,QPRES) = qzmo_core(i,j,k,QREINT)*(qzmo_core(i,j,k,QGAME)-ONE)
                   qzmo_core(i,j,k,QPRES) = max(qzmo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qzmo_core(i,j,k,QPRES) = qzm_core(i,j,k,QPRES)
                qzmo_core(i,j,k,QGAME) = qzm_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzmo_core, qzcmo_lo, qzcmo_hi, &
                                          qzmo_pass, qzpmo_lo, qzpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qzmo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qzmo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo_core(i,j,k,QPRES)
             qzmo_rad(i,j,k,qreitot) = sum(qzmo_rad(i,j,k,qrad:qradhi)) + qzmo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transx_on_zstates


  !===========================================================================
  ! transy routines
  !===========================================================================
  subroutine transy_on_xstates(lo, hi, &
                               qxm_core, qxcm_lo, qxcm_hi, &
                               qxmo_core, qxcmo_lo, qxcmo_hi, &
                               qxm_pass, qxpm_lo, qxpm_hi, &
                               qxmo_pass, qxpmo_lo, qxpmo_hi, &
#ifdef RADIATION
                               qxm_rad, qxrm_lo, qxrm_hi, &
                               qxmo_rad, qxrmo_lo, qxrmo_hi, &
#endif
                               qxp_core, qxcp_lo, qxcp_hi, &
                               qxpo_core, qxcpo_lo, qxcpo_hi, &
                               qxp_pass, qxpp_lo, qxpp_hi, &
                               qxpo_pass, qxppo_lo, qxppo_hi, &
#ifdef RADIATION
                               qxp_rad, qxrp_lo, qxrp_hi, &
                               qxpo_rad, qxrpo_lo, qxrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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


    integer, intent(in) :: qxcm_lo(3), qxcm_hi(3)
    integer, intent(in) :: qxcp_lo(3), qxcp_hi(3)
    integer, intent(in) :: qxpm_lo(3), qxpm_hi(3)
    integer, intent(in) :: qxpp_lo(3), qxpp_hi(3)
    integer, intent(in) :: qxcmo_lo(3), qxcmo_hi(3)
    integer, intent(in) :: qxpmo_lo(3), qxpmo_hi(3)
    integer, intent(in) :: qxcpo_lo(3), qxcpo_hi(3)
    integer, intent(in) :: qxppo_lo(3), qxppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qxrm_lo(3), qxrm_hi(3)
    integer, intent(in) :: qxrp_lo(3), qxrp_hi(3)
    integer, intent(in) :: qxrmo_lo(3), qxrmo_hi(3)
    integer, intent(in) :: qxrpo_lo(3), qxrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qxm_core(qxcm_lo(1):qxcm_hi(1),qxcm_lo(2):qxcm_hi(2),qxcm_lo(3):qxcm_hi(3),NQC)
    real(rt), intent(in) :: qxp_core(qxcp_lo(1):qxcp_hi(1),qxcp_lo(2):qxcp_hi(2),qxcp_lo(3):qxcp_hi(3),NQC)
    real(rt), intent(in) :: qxm_pass(qxpm_lo(1):qxpm_hi(1),qxpm_lo(2):qxpm_hi(2),qxpm_lo(3):qxpm_hi(3),NQP)
    real(rt), intent(in) :: qxp_pass(qxpp_lo(1):qxpp_hi(1),qxpp_lo(2):qxpp_hi(2),qxpp_lo(3):qxpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qxm_rad(qxrm_lo(1):qxrm_hi(1),qxrm_lo(2):qxrm_hi(2),qxrm_lo(3):qxrm_hi(3),NQR)
    real(rt), intent(in) ::  qxp_rad(qxrp_lo(1):qxrp_hi(1),qxrp_lo(2):qxrp_hi(2),qxrp_lo(3):qxrp_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    real(rt), intent(out) :: qxmo_core(qxcmo_lo(1):qxcmo_hi(1),qxcmo_lo(2):qxcmo_hi(2),qxcmo_lo(3):qxcmo_hi(3),NQC)
    real(rt), intent(out) :: qxpo_core(qxcpo_lo(1):qxcpo_hi(1),qxcpo_lo(2):qxcpo_hi(2),qxcpo_lo(3):qxcpo_hi(3),NQC)
    real(rt), intent(out) :: qxmo_pass(qxpmo_lo(1):qxpmo_hi(1),qxpmo_lo(2):qxpmo_hi(2),qxpmo_lo(3):qxpmo_hi(3),NQP)
    real(rt), intent(out) :: qxpo_pass(qxppo_lo(1):qxppo_hi(1),qxppo_lo(2):qxppo_hi(2),qxppo_lo(3):qxppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qxmo_rad(qxrmo_lo(1):qxrmo_hi(1),qxrmo_lo(2):qxrmo_hi(2),qxrmo_lo(3):qxrmo_hi(3),NQR)
    real(rt), intent(out) ::  qxpo_rad(qxrpo_lo(1):qxrpo_hi(1),qxrpo_lo(2):qxrpo_hi(2),qxrpo_lo(3):qxrpo_hi(3),NQR)
#endif

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qxp_core(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = qxp_core(i,j,k,QRHO)*qxp_pass(i,j,k,nqs) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                qxpo_pass(i,j,k,nqs) = compu/rrnew

                rrnew = qxm_core(i,j,k,QRHO) - cdtdy*(fy(i-1,j+1,k,URHO) - fy(i-1,j,k,URHO))
                compu = qxm_core(i,j,k,QRHO)*qxm_pass(i,j,k,nqs) - cdtdy*(fy(i-1,j+1,k,n) - fy(i-1,j,k,n))
                qxmo_pass(i,j,k,nqs) = compu/rrnew
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
             rrrx = qxp_core(i,j,k,QRHO)
             rurx = rrrx*qxp_core(i,j,k,QU)
             rvrx = rrrx*qxp_core(i,j,k,QV)
             rwrx = rrrx*qxp_core(i,j,k,QW)
             ekenrx = HALF*rrrx*sum(qxp_core(i,j,k,QU:QW)**2)
             rerx = qxp_core(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
             err  = qxp_rad(i,j,k,qrad:qradhi)
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
             qxpo_core(i,j,k,QRHO) = rrnewrx
             rhoinv = ONE/rrnewrx
             qxpo_core(i,j,k,QU) = runewrx*rhoinv
             qxpo_core(i,j,k,QV) = rvnewrx*rhoinv
             qxpo_core(i,j,k,QW) = rwnewrx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
             qxpo_core(i,j,k,QREINT) = renewrx - rhoekenrx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxpo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qxpo_Core(i,j,k,QREINT) = qxp_core(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp_core(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qxpo_core(i,j,k,QPRES) = max(pnewrx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxpo_core(i,j,k,QGAME) = qxp_core(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxpo_core(i,j,k,QPRES) = qxpo_core(i,j,k,QREINT)*(qxpo_core(i,j,k,QGAME)-ONE)
                   qxpo_core(i,j,k,QPRES) = max(qxpo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qxpo_core(i,j,k,QPRES) = qxp_core(i,j,k,QPRES)
                qxpo_core(i,j,k,QGAME) = qxp_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qxpo_core, qxcpo_lo, qxcpo_hi, &
                                          qxpo_pass, qxppo_lo, qxppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qxpo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qxpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qxpo_core(i,j,k,QPRES)
             qxpo_rad(i,j,k,qreitot) = sum(qxpo_rad(i,j,k,qrad:qradhi)) + qxpo_core(i,j,k,QREINT)
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
             rrlx = qxm_core(i,j,k,QRHO)
             rulx = rrlx*qxm_core(i,j,k,QU)
             rvlx = rrlx*qxm_core(i,j,k,QV)
             rwlx = rrlx*qxm_core(i,j,k,QW)
             ekenlx = HALF*rrlx*sum(qxm_core(i,j,k,QU:QW)**2)
             relx = qxm_core(i,j,k,QREINT) + ekenlx
#ifdef RADIATION
             erl  = qxm_rad(i,j,k,qrad:qradhi)
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

             qxmo_core(i,j,k,QRHO) = rrnewlx
             rhoinv = ONE/rrnewlx
             qxmo_core(i,j,k,QU) = runewlx*rhoinv
             qxmo_core(i,j,k,QV) = rvnewlx*rhoinv
             qxmo_core(i,j,k,QW) = rwnewlx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
             qxmo_core(i,j,k,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo_core(i,j,k,QREINT) = qxm_core(i,j,k,QREINT) - &
                        cdtdy*(fy(i-1,j+1,k,UEINT) - fy(i-1,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm_core(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qxmo_core(i,j,k,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo_core(i,j,k,QGAME) = qxm_core(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo_core(i,j,k,QPRES) = qxmo_core(i,j,k,QREINT)*(qxmo_core(i,j,k,QGAME)-ONE)
                   qxmo_core(i,j,k,QPRES) = max(qxmo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qxmo_core(i,j,k,QPRES) = qxm_core(i,j,k,QPRES)
                qxmo_core(i,j,k,QGAME) = qxm_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qxmo_core, qxcmo_lo, qxcmo_hi, &
                                          qxmo_pass, qxpmo_lo, qxpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qxmo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qxmo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qxmo_core(i,j,k,QPRES)
             qxmo_rad(i,j,k,qreitot) = sum(qxmo_rad(i,j,k,qrad:qradhi)) + qxmo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transy_on_xstates

  subroutine transy_on_zstates(lo, hi, &
                               qzm_core, qzcm_lo, qzcm_hi, &
                               qzmo_core, qzcmo_lo, qzcmo_hi, &
                               qzm_pass, qzpm_lo, qzpm_hi, &
                               qzmo_pass, qzpmo_lo, qzpmo_hi, &
#ifdef RADIATION
                               qzm_rad, qzrm_lo, qzrm_hi, &
                               qzmo_rad, qzrmo_lo, qzrmo_hi, &
#endif
                               qzp_core, qzcp_lo, qzcp_hi, &
                               qzpo_core, qzcpo_lo, qzcpo_hi, &
                               qzp_pass, qzpp_lo, qzpp_hi, &
                               qzpo_pass, qzppo_lo, qzppo_hi, &
#ifdef RADIATION
                               qzp_rad, qzrp_lo, qzrp_hi, &
                               qzpo_rad, qzrpo_lo, qzrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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

    integer, intent(in) :: qzcm_lo(3), qzcm_hi(3)
    integer, intent(in) :: qzpm_lo(3), qzpm_hi(3)
    integer, intent(in) :: qzcp_lo(3), qzcp_hi(3)
    integer, intent(in) :: qzpp_lo(3), qzpp_hi(3)
    integer, intent(in) :: qzcmo_lo(3), qzcmo_hi(3)
    integer, intent(in) :: qzcpo_lo(3), qzcpo_hi(3)
    integer, intent(in) :: qzpmo_lo(3), qzpmo_hi(3)
    integer, intent(in) :: qzppo_lo(3), qzppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qzrm_lo(3), qzrm_hi(3)
    integer, intent(in) :: qzrp_lo(3), qzrp_hi(3)
    integer, intent(in) :: qzrmo_lo(3), qzrmo_hi(3)
    integer, intent(in) :: qzrpo_lo(3), qzrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qzm_core(qzcm_lo(1):qzcm_hi(1),qzcm_lo(2):qzcm_hi(2),qzcm_lo(3):qzcm_hi(3),NQC)
    real(rt), intent(in) :: qzp_core(qzcp_lo(1):qzcp_hi(1),qzcp_lo(2):qzcp_hi(2),qzcp_lo(3):qzcp_hi(3),NQC)
    real(rt), intent(in) :: qzm_pass(qzpm_lo(1):qzpm_hi(1),qzpm_lo(2):qzpm_hi(2),qzpm_lo(3):qzpm_hi(3),NQP)
    real(rt), intent(in) :: qzp_pass(qzpp_lo(1):qzpp_hi(1),qzpp_lo(2):qzpp_hi(2),qzpp_lo(3):qzpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qzm_rad(qzrm_lo(1):qzrm_hi(1),qzrm_lo(2):qzrm_hi(2),qzrm_lo(3):qzrm_hi(3),NQR)
    real(rt), intent(in) ::  qzp_rad(qzrp_lo(1):qzrp_hi(1),qzrp_lo(2):qzrp_hi(2),qzrp_lo(3):qzrp_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    real(rt), intent(out) :: qzmo_core(qzcmo_lo(1):qzcmo_hi(1),qzcmo_lo(2):qzcmo_hi(2),qzcmo_lo(3):qzcmo_hi(3),NQC)
    real(rt), intent(out) :: qzpo_core(qzcpo_lo(1):qzcpo_hi(1),qzcpo_lo(2):qzcpo_hi(2),qzcpo_lo(3):qzcpo_hi(3),NQC)
    real(rt), intent(out) :: qzmo_pass(qzpmo_lo(1):qzpmo_hi(1),qzpmo_lo(2):qzpmo_hi(2),qzpmo_lo(3):qzpmo_hi(3),NQP)
    real(rt), intent(out) :: qzpo_pass(qzppo_lo(1):qzppo_hi(1),qzppo_lo(2):qzppo_hi(2),qzppo_lo(3):qzppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qzmo_rad(qzrmo_lo(1):qzrmo_hi(1),qzrmo_lo(2):qzrmo_hi(2),qzrmo_lo(3):qzrmo_hi(3),NQR)
    real(rt), intent(out) ::  qzpo_rad(qzrpo_lo(1):qzrpo_hi(1),qzrpo_lo(2):qzrpo_hi(2),qzrpo_lo(3):qzrpo_hi(3),NQR)
#endif

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qzp_core(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k,URHO) - fy(i,j,k,URHO))
                compu = qzp_core(i,j,k,QRHO)*qzp_pass(i,j,k,nqs) - cdtdy*(fy(i,j+1,k,n) - fy(i,j,k,n))
                qzpo_pass(i,j,k,nqs) = compu/rrnew

                rrnew = qzm_core(i,j,k,QRHO) - cdtdy*(fy(i,j+1,k-1,URHO) - fy(i,j,k-1,URHO))
                compu = qzm_core(i,j,k,QRHO)*qzm_pass(i,j,k,nqs) - cdtdy*(fy(i,j+1,k-1,n) - fy(i,j,k-1,n))
                qzmo_pass(i,j,k,nqs) = compu/rrnew
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
             rrrz = qzp_core(i,j,k,QRHO)
             rurz = rrrz*qzp_core(i,j,k,QU)
             rvrz = rrrz*qzp_core(i,j,k,QV)
             rwrz = rrrz*qzp_core(i,j,k,QW)
             ekenrz = HALF*rrrz*sum(qzp_core(i,j,k,QU:QW)**2)
             rerz = qzp_core(i,j,k,QREINT) + ekenrz
#ifdef RADIATION
             err  = qzp_rad(i,j,k,qrad:qradhi)
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
             qzpo_core(i,j,k,QRHO) = rrnewrz
             rhoinv = ONE/rrnewrz
             qzpo_core(i,j,k,QU) = runewrz*rhoinv
             qzpo_core(i,j,k,QV) = rvnewrz*rhoinv
             qzpo_core(i,j,k,QW) = rwnewrz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)*rhoinv
             qzpo_core(i,j,k,QREINT) = renewrz - rhoekenrz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzpo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzpo_core(i,j,k,QREINT) = qzp_core(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k,UEINT) - fy(i,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrz = qzp_core(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qzpo_core(i,j,k,QPRES) = max(pnewrz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzpo_core(i,j,k,QGAME) = qzp_core(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzpo_core(i,j,k,QPRES) = qzpo_core(i,j,k,QREINT)*(qzpo_core(i,j,k,QGAME)-ONE)
                   qzpo_core(i,j,k,QPRES) = max(qzpo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qzpo_core(i,j,k,QPRES) = qzp_core(i,j,k,QPRES)
                qzpo_core(i,j,k,QGAME) = qzp_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzpo_core, qzcpo_lo, qzcpo_hi, &
                                          qzpo_pass, qzppo_lo, qzppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qzpo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qzpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qzpo_core(i,j,k,QPRES)
             qzpo_rad(i,j,k,qreitot) = sum(qzpo_rad(i,j,k,qrad:qradhi)) + qzpo_core(i,j,k,QREINT)
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
             rrlz = qzm_core(i,j,k,QRHO)
             rulz = rrlz*qzm_core(i,j,k,QU)
             rvlz = rrlz*qzm_core(i,j,k,QV)
             rwlz = rrlz*qzm_core(i,j,k,QW)
             ekenlz = HALF*rrlz*sum(qzm_core(i,j,k,QU:QW)**2)
             relz = qzm_core(i,j,k,QREINT) + ekenlz
#ifdef RADIATION
             erl  = qzm_rad(i,j,k,qrad:qradhi)
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
             qzmo_core(i,j,k,QRHO) = rrnewlz
             rhoinv = ONE/rrnewlz
             qzmo_core(i,j,k,QU) = runewlz*rhoinv
             qzmo_core(i,j,k,QV) = rvnewlz*rhoinv
             qzmo_core(i,j,k,QW) = rwnewlz*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)*rhoinv
             qzmo_core(i,j,k,QREINT) = renewlz - rhoekenlz

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qzmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qzmo_core(i,j,k,QREINT) = qzm_core(i,j,k,QREINT) - &
                        cdtdy*(fy(i,j+1,k-1,UEINT) - fy(i,j,k-1,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlz = qzm_core(i,j,k,QPRES) - cdtdy*(dup + pav*du*(gamc - ONE))
                   qzmo_core(i,j,k,QPRES) = max(pnewlz,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qzmo_core(i,j,k,QGAME) = qzm_core(i,j,k,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qzmo_core(i,j,k,QPRES) = qzmo_core(i,j,k,QREINT)*(qzmo_core(i,j,k,QGAME)-ONE)
                   qzmo_core(i,j,k,QPRES) = max(qzmo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qzmo_core(i,j,k,QPRES) = qzm_core(i,j,k,QPRES)
                qzmo_core(i,j,k,QGAME) = qzm_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qzmo_core, qzcmo_lo, qzcmo_hi, &
                                          qzmo_pass, qzpmo_lo, qzpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qzmo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qzmo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qzmo_core(i,j,k,QPRES)
             qzmo_rad(i,j,k,qreitot) = sum(qzmo_rad(i,j,k,qrad:qradhi)) + qzmo_core(i,j,k,QREINT)
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
                               qxm_core, qxcm_lo, qxcm_hi, &
                               qxmo_core, qxcmo_lo, qxcmo_hi, &
                               qxm_pass, qxpm_lo, qxpm_hi, &
                               qxmo_pass, qxpmo_lo, qxpmo_hi, &
#ifdef RADIATION
                               qxm_rad, qxrm_lo, qxrm_hi, &
                               qxmo_rad, qxrmo_lo, qxrmo_hi, &
#endif
                               qxp_core, qxcp_lo, qxcp_hi, &
                               qxpo_core, qxcpo_lo, qxcpo_hi, &
                               qxp_pass, qxpp_lo, qxpp_hi, &
                               qxpo_pass, qxppo_lo, qxppo_hi, &
#ifdef RADIATION
                               qxp_rad, qxrp_lo, qxrp_hi, &
                               qxpo_rad, qxrpo_lo, qxrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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


    integer, intent(in) :: qxcm_lo(3), qxcm_hi(3)
    integer, intent(in) :: qxcmo_lo(3), qxcmo_hi(3)
    integer, intent(in) :: qxpm_lo(3), qxpm_hi(3)
    integer, intent(in) :: qxpmo_lo(3), qxpmo_hi(3)
    integer, intent(in) :: qxcp_lo(3), qxcp_hi(3)
    integer, intent(in) :: qxcpo_lo(3), qxcpo_hi(3)
    integer, intent(in) :: qxpp_lo(3), qxpp_hi(3)
    integer, intent(in) :: qxppo_lo(3), qxppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qxrm_lo(3), qxrm_hi(3)
    integer, intent(in) :: qxrmo_lo(3), qxrmo_hi(3)
    integer, intent(in) :: qxrp_lo(3), qxrp_hi(3)
    integer, intent(in) :: qxrpo_lo(3), qxrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qxm_core(qxcm_lo(1):qxcm_hi(1),qxcm_lo(2):qxcm_hi(2),qxcm_lo(3):qxcm_hi(3),NQC)
    real(rt), intent(in) :: qxp_core(qxcp_lo(1):qxcp_hi(1),qxcp_lo(2):qxcp_hi(2),qxcp_lo(3):qxcp_hi(3),NQC)
    real(rt), intent(in) :: qxm_pass(qxpm_lo(1):qxpm_hi(1),qxpm_lo(2):qxpm_hi(2),qxpm_lo(3):qxpm_hi(3),NQP)
    real(rt), intent(in) :: qxp_pass(qxpp_lo(1):qxpp_hi(1),qxpp_lo(2):qxpp_hi(2),qxpp_lo(3):qxpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qxm_rad(qxrm_lo(1):qxrm_hi(1),qxrm_lo(2):qxrm_hi(2),qxrm_lo(3):qxrm_hi(3),NQR)
    real(rt), intent(in) ::  qxp_rad(qxrp_lo(1):qxrp_hi(1),qxrp_lo(2):qxrp_hi(2),qxrp_lo(3):qxrp_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    real(rt), intent(out) :: qxmo_core(qxcmo_lo(1):qxcmo_hi(1),qxcmo_lo(2):qxcmo_hi(2),qxcmo_lo(3):qxcmo_hi(3),NQC)
    real(rt), intent(out) :: qxpo_core(qxcpo_lo(1):qxcpo_hi(1),qxcpo_lo(2):qxcpo_hi(2),qxcpo_lo(3):qxcpo_hi(3),NQC)
    real(rt), intent(out) :: qxmo_pass(qxpmo_lo(1):qxpmo_hi(1),qxpmo_lo(2):qxpmo_hi(2),qxpmo_lo(3):qxpmo_hi(3),NQP)
    real(rt), intent(out) :: qxpo_pass(qxppo_lo(1):qxppo_hi(1),qxppo_lo(2):qxppo_hi(2),qxppo_lo(3):qxppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qxmo_rad(qxrmo_lo(1):qxrmo_hi(1),qxrmo_lo(2):qxrmo_hi(2),qxrmo_lo(3):qxrmo_hi(3),NQR)
    real(rt), intent(out) ::  qxpo_rad(qxrpo_lo(1):qxrpo_hi(1),qxrpo_lo(2):qxrpo_hi(2),qxrpo_lo(3):qxrpo_hi(3),NQR)
#endif

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

        do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qxp_core(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                compu = qxp_core(i,j,k,QRHO)*qxp_pass(i,j,k,nqs) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                qxpo_pass(i,j,k,nqs) = compu/rrnew

                rrnew = qxm_core(i,j,k,QRHO) - cdtdz*(fz(i-1,j,k+1,URHO) - fz(i-1,j,k,URHO))
                compu = qxm_core(i,j,k,QRHO)*qxm_pass(i,j,k,nqs) - cdtdz*(fz(i-1,j,k+1,n) - fz(i-1,j,k,n))
                qxmo_pass(i,j,k,nqs) = compu/rrnew
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
             rrrx = qxp_core(i,j,k,QRHO)
             rurx = rrrx*qxp_core(i,j,k,QU)
             rvrx = rrrx*qxp_core(i,j,k,QV)
             rwrx = rrrx*qxp_core(i,j,k,QW)
             ekenrx = HALF*rrrx*sum(qxp_core(i,j,k,QU:QW)**2)
             rerx = qxp_core(i,j,k,QREINT) + ekenrx
#ifdef RADIATION
             errx = qxp_rad(i,j,k,qrad:qradhi)
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
             qxpo_core(i,j,k,QRHO) = rrnewrx
             rhoinv = ONE/rrnewrx
             qxpo_core(i,j,k,QU) = runewrx*rhoinv
             qxpo_core(i,j,k,QV) = rvnewrx*rhoinv
             qxpo_core(i,j,k,QW) = rwnewrx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)*rhoinv
             qxpo_core(i,j,k,QREINT) = renewrx - rhoekenrx

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qxpo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qxpo_core(i,j,k,QREINT) = qxp_core(i,j,k,QREINT) - &
                        cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp_core(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qxpo_core(i,j,k,QPRES) = max(pnewrx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxpo_core(i,j,k,QGAME) = qxp_core(i,j,k,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxpo_core(i,j,k,QPRES) = qxpo_core(i,j,k,QREINT)*(qxpo_core(i,j,k,QGAME)-ONE)
                   qxpo_core(i,j,k,QPRES) = max(qxpo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qxpo_core(i,j,k,QPRES) = qxp_core(i,j,k,QPRES)
                qxpo_core(i,j,k,QGAME) = qxp_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qxpo_core, qxcpo_lo, qxcpo_hi, &
                                          qxpo_pass, qxppo_lo, qxppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qxpo_rad(i,j,k,qrad:qradhi) = ernewrx(:)
             qxpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewrx(:)) + qxpo_core(i,j,k,QPRES)
             qxpo_pass(i,j,k,qreitot) = sum(qxpo_pass(i,j,k,qrad:qradhi)) + qxpo_core(i,j,k,QREINT)
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
             rrlx = qxm_core(i,j,k,QRHO)
             rulx = rrlx*qxm_core(i,j,k,QU)
             rvlx = rrlx*qxm_core(i,j,k,QV)
             rwlx = rrlx*qxm_core(i,j,k,QW)
             ekenlx = HALF*rrlx*sum(qxm_core(i,j,k,QU:QW)**2)
             relx = qxm_core(i,j,k,QREINT) + ekenlx
#ifdef RADIATION
             erlx = qxm_rad(i,j,k,qrad:qradhi)
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
             qxmo_core(i,j,k,QRHO) = rrnewlx
             rhoinv = ONE/rrnewlx
             qxmo_core(i,j,k,QU) = runewlx*rhoinv
             qxmo_core(i,j,k,QV) = rvnewlx*rhoinv
             qxmo_core(i,j,k,QW) = rwnewlx*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)*rhoinv
             qxmo_core(i,j,k,QREINT) = renewlx - rhoekenlx

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qxmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qxmo_core(i,j,k,QREINT) = qxm_core(i,j,k,QREINT) - &
                        cdtdz*(fz(i-1,j,k+1,UEINT) - fz(i-1,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm_core(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qxmo_core(i,j,k,QPRES) = max(pnewlx,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qxmo_core(i,j,k,QGAME) = qxm_core(i,j,k,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qxmo_core(i,j,k,QPRES) = qxmo_core(i,j,k,QREINT)*(qxmo_core(i,j,k,QGAME)-ONE)
                   qxmo_core(i,j,k,QPRES) = max(qxmo_core(i,j,k,QPRES), small_pres)
                end if
             else
                qxmo_core(i,j,k,QPRES) = qxm_core(i,j,k,QPRES)
                qxmo_core(i,j,k,QGAME) = qxm_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qxmo_core, qxcmo_lo, qxcmo_hi, &
                                          qxmo_pass, qxpmo_lo, qxpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qxmo_rad(i,j,k,qrad:qradhi) = ernewlx(:)
             qxmo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewlx(:)) + qxmo_core(i,j,k,QPRES)
             qxmo_rad(i,j,k,qreitot) = sum(qxmo_rad(i,j,k,qrad:qradhi)) + qxmo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transz_on_xstates


  subroutine transz_on_ystates(lo, hi, &
                               qym_core, qycm_lo, qycm_hi, &
                               qymo_core, qycmo_lo, qycmo_hi, &
                               qym_pass, qypm_lo, qypm_hi, &
                               qymo_pass, qypmo_lo, qypmo_hi, &
#ifdef RADIATION
                               qym_rad, qyrm_lo, qyrm_hi, &
                               qymo_rad, qyrmo_lo, qyrmo_hi, &
#endif
                               qyp_core, qycp_lo, qycp_hi, &
                               qypo_core, qycpo_lo, qycpo_hi, &
                               qyp_pass, qypp_lo, qypp_hi, &
                               qypo_pass, qyppo_lo, qyppo_hi, &
#ifdef RADIATION
                               qyp_rad, qyrp_lo, qyrp_hi, &
                               qypo_rad, qyrpo_lo, qyrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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


    integer, intent(in) :: qycm_lo(3), qycm_hi(3)
    integer, intent(in) :: qycmo_lo(3), qycmo_hi(3)
    integer, intent(in) :: qycp_lo(3), qycp_hi(3)
    integer, intent(in) :: qycpo_lo(3), qycpo_hi(3)
    integer, intent(in) :: qypm_lo(3), qypm_hi(3)
    integer, intent(in) :: qypmo_lo(3), qypmo_hi(3)
    integer, intent(in) :: qypp_lo(3), qypp_hi(3)
    integer, intent(in) :: qyppo_lo(3), qyppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qyrm_lo(3), qyrm_hi(3)
    integer, intent(in) :: qyrmo_lo(3), qyrmo_hi(3)
    integer, intent(in) :: qyrp_lo(3), qyrp_hi(3)
    integer, intent(in) :: qyrpo_lo(3), qyrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qym_core(qycm_lo(1):qycm_hi(1),qycm_lo(2):qycm_hi(2),qycm_lo(3):qycm_hi(3),NQC)
    real(rt), intent(in) :: qyp_core(qycp_lo(1):qycp_hi(1),qycp_lo(2):qycp_hi(2),qycp_lo(3):qycp_hi(3),NQC)
    real(rt), intent(in) :: qym_pass(qypm_lo(1):qypm_hi(1),qypm_lo(2):qypm_hi(2),qypm_lo(3):qypm_hi(3),NQP)
    real(rt), intent(in) :: qyp_pass(qypp_lo(1):qypp_hi(1),qypp_lo(2):qypp_hi(2),qypp_lo(3):qypp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qym_rad(qyrm_lo(1):qyrm_hi(1),qyrm_lo(2):qyrm_hi(2),qyrm_lo(3):qyrm_hi(3),NQR)
    real(rt), intent(in) ::  qyp_rad(qyrp_lo(1):qyrp_hi(1),qyrp_lo(2):qyrp_hi(2),qyrp_lo(3):qyrp_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)
    real(rt), intent(in) :: fz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),NVAR)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    real(rt), intent(out) :: qymo_core(qycmo_lo(1):qycmo_hi(1),qycmo_lo(2):qycmo_hi(2),qycmo_lo(3):qycmo_hi(3),NQC)
    real(rt), intent(out) :: qypo_core(qycpo_lo(1):qycpo_hi(1),qycpo_lo(2):qycpo_hi(2),qycpo_lo(3):qycpo_hi(3),NQC)
    real(rt), intent(out) :: qymo_pass(qypmo_lo(1):qypmo_hi(1),qypmo_lo(2):qypmo_hi(2),qypmo_lo(3):qypmo_hi(3),NQP)
    real(rt), intent(out) :: qypo_pass(qyppo_lo(1):qyppo_hi(1),qyppo_lo(2):qyppo_hi(2),qyppo_lo(3):qyppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qymo_rad(qyrmo_lo(1):qyrmo_hi(1),qyrmo_lo(2):qyrmo_hi(2),qyrmo_lo(3):qyrmo_hi(3),NQR)
    real(rt), intent(out) ::  qypo_rad(qyrpo_lo(1):qyrpo_hi(1),qyrpo_lo(2):qyrpo_hi(2),qyrpo_lo(3):qyrpo_hi(3),NQR)
#endif

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrnew = qyp_core(i,j,k,QRHO) - cdtdz*(fz(i,j,k+1,URHO) - fz(i,j,k,URHO))
                compu = qyp_core(i,j,k,QRHO)*qyp_pass(i,j,k,nqs) - cdtdz*(fz(i,j,k+1,n) - fz(i,j,k,n))
                qypo_pass(i,j,k,nqs) = compu/rrnew

                rrnew = qym_core(i,j,k,QRHO) - cdtdz*(fz(i,j-1,k+1,URHO) - fz(i,j-1,k,URHO))
                compu = qym_core(i,j,k,QRHO)*qym_pass(i,j,k,nqs) - cdtdz*(fz(i,j-1,k+1,n) - fz(i,j-1,k,n))
                qymo_pass(i,j,k,nqs) = compu/rrnew
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
             rrry = qyp_core(i,j,k,QRHO)
             rury = rrry*qyp_core(i,j,k,QU)
             rvry = rrry*qyp_core(i,j,k,QV)
             rwry = rrry*qyp_core(i,j,k,QW)
             ekenry = HALF*rrry*sum(qyp_core(i,j,k,QU:QW)**2)
             rery = qyp_core(i,j,k,QREINT) + ekenry
#ifdef RADIATION
             erry = qyp_rad(i,j,k,qrad:qradhi)
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
             qypo_core(i,j,k,QRHO) = rrnewry
             rhoinv = ONE/rrnewry
             qypo_core(i,j,k,QU) = runewry*rhoinv
             qypo_core(i,j,k,QV) = rvnewry*rhoinv
             qypo_core(i,j,k,QW) = rwnewry*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)*rhoinv
             qypo_core(i,j,k,QREINT) = renewry - rhoekenry

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qypo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qypo_core(i,j,k,QREINT) = qyp_core(i,j,k,QREINT) - &
                        cdtdz*(fz(i,j,k+1,UEINT) - fz(i,j,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp_core(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qypo_core(i,j,k,QPRES) = max(pnewry,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qypo_core(i,j,k,QGAME) = qyp_core(i,j,k,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qypo_core(i,j,k,QPRES) = qypo_core(i,j,k,QREINT)*(qypo_core(i,j,k,QGAME)-ONE)
                   qypo_core(i,j,k,QPRES) = max(qypo_core(i,j,k,QPRES), small_pres)
                endif
             else
                qypo_core(i,j,k,QPRES) = qyp_core(i,j,k,QPRES)
                qypo_core(i,j,k,QGAME) = qyp_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qypo_core, qycpo_lo, qycpo_hi, &
                                          qypo_pass, qyppo_lo, qyppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qypo_rad(i,j,k,qrad:qradhi) = ernewry(:)
             qypo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewry(:)) + qypo_core(i,j,k,QPRES)
             qypo_rad(i,j,k,qreitot) = sum(qypo_rad(i,j,k,qrad:qradhi)) + qypo_core(i,j,k,QREINT)
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
             rrly = qym_core(i,j,k,QRHO)
             ruly = rrly*qym_core(i,j,k,QU)
             rvly = rrly*qym_core(i,j,k,QV)
             rwly = rrly*qym_core(i,j,k,QW)
             ekenly = HALF*rrly*sum(qym_core(i,j,k,QU:QW)**2)
             rely = qym_core(i,j,k,QREINT) + ekenly
#ifdef RADIATION
             erly = qym_rad(i,j,k,qrad:qradhi)
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
             qymo_core(i,j,k,QRHO) = rrnewly
             rhoinv = ONE/rrnewly
             qymo_core(i,j,k,QU) = runewly*rhoinv
             qymo_core(i,j,k,QV) = rvnewly*rhoinv
             qymo_core(i,j,k,QW) = rwnewly*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)*rhoinv
             qymo_core(i,j,k,QREINT) = renewly - rhoekenly

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qymo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qymo_core(i,j,k,QREINT) = qym_core(i,j,k,QREINT) - &
                        cdtdz*(fz(i,j-1,k+1,UEINT) - fz(i,j-1,k,UEINT) + pav*du)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym_core(i,j,k,QPRES) - cdtdz*(dup + pav*du*(gamc - ONE))
                   qymo_core(i,j,k,QPRES) = max(pnewly,small_pres)
                else
                   ! Update gammae with its transverse terms
                   qymo_core(i,j,k,QGAME) = qym_core(i,j,k,QGAME) + &
                        cdtdz*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qymo_core(i,j,k,QPRES) = qymo_core(i,j,k,QREINT)*(qymo_core(i,j,k,QGAME)-ONE)
                   qymo_core(i,j,k,QPRES) = max(qymo_Core(i,j,k,QPRES), small_pres)
                endif
             else
                qymo_core(i,j,k,QPRES) = qym_core(i,j,k,QPRES)
                qymo_core(i,j,k,QGAME) = qym_core(i,j,k,QGAME)
             endif

             call reset_edge_state_thermo(qymo_core, qycmo_lo, qycmo_hi, &
                                          qymo_pass, qypmo_lo, qypmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qymo_rad(i,j,k,qrad:qradhi) = ernewly(:)
             qymo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewly(:)) + qymo_core(i,j,k,QPRES)
             qymo_rad(i,j,k,qreitot) = sum(qymo_rad(i,j,k,qrad:qradhi)) + qymo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transz_on_ystates

  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz(lo, hi, &
                     qm_core, qcm_lo, qcm_hi, &
                     qmo_core, qcmo_lo, qcmo_hi, &
                     qm_pass, qpm_lo, qpm_hi, &
                     qmo_pass, qpmo_lo, qpmo_hi, &
#ifdef RADIATION
                     qm_rad, qrm_lo, qrm_hi, &
                     qmo_rad, qrmo_lo, qrmo_hi, &
#endif
                     qp_core, qcp_lo, qcp_hi, &
                     qpo_core, qcpo_lo, qcpo_hi, &
                     qp_pass, qpp_lo, qpp_hi, &
                     qpo_pass, qppo_lo, qppo_hi, &
#ifdef RADIATION
                     qp_rad, qrp_lo, qrp_hi, &
                     qpo_rad, qrpo_lo, qrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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


    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcmo_lo(3), qcmo_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)
    integer, intent(in) :: qcpo_lo(3), qcpo_hi(3)
    integer, intent(in) :: qpm_lo(3), qpm_hi(3)
    integer, intent(in) :: qpmo_lo(3), qpmo_hi(3)
    integer, intent(in) :: qpp_lo(3), qpp_hi(3)
    integer, intent(in) :: qppo_lo(3), qppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qrm_lo(3), qrm_hi(3)
    integer, intent(in) :: qrmo_lo(3), qrmo_hi(3)
    integer, intent(in) :: qrp_lo(3), qrp_hi(3)
    integer, intent(in) :: qrpo_lo(3), qrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(in) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)
    real(rt), intent(in) :: qm_pass(qpm_lo(1):qpm_hi(1),qpm_lo(2):qpm_hi(2),qpm_lo(3):qpm_hi(3),NQP)
    real(rt), intent(in) :: qp_pass(qpp_lo(1):qpp_hi(1),qpp_lo(2):qpp_hi(2),qpp_lo(3):qpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qm_rad(qrm_lo(1):qrm_hi(1),qrm_lo(2):qrm_hi(2),qrm_lo(3):qrm_hi(3),NQR)
    real(rt), intent(in) ::  qp_rad(qrp_lo(1):qrp_hi(1),qrp_lo(2):qrp_hi(2),qrp_lo(3):qrp_hi(3),NQR)
#endif
    real(rt), intent(out) :: qmo_core(qcmo_lo(1):qcmo_hi(1),qcmo_lo(2):qcmo_hi(2),qcmo_lo(3):qcmo_hi(3),NQC)
    real(rt), intent(out) :: qpo_core(qcpo_lo(1):qcpo_hi(1),qcpo_lo(2):qcpo_hi(2),qcpo_lo(3):qcpo_hi(3),NQC)
    real(rt), intent(out) :: qmo_pass(qpmo_lo(1):qpmo_hi(1),qpmo_lo(2):qpmo_hi(2),qpmo_lo(3):qpmo_hi(3),NQP)
    real(rt), intent(out) :: qpo_pass(qppo_lo(1):qppo_hi(1),qppo_lo(2):qppo_hi(2),qppo_lo(3):qppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qmo_rad(qrmo_lo(1):qrmo_hi(1),qrmo_lo(2):qrmo_hi(2),qrmo_lo(3):qrmo_hi(3),NQR)
    real(rt), intent(out) ::  qpo_rad(qrpo_lo(1):qrpo_hi(1),qrpo_lo(2):qrpo_hi(2),qrpo_lo(3):qrpo_hi(3),NQR)
#endif

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fyz(fyz_lo(1):fyz_hi(1),fyz_lo(2):fyz_hi(2),fyz_lo(3):fyz_hi(3),NVAR)
    real(rt), intent(in) :: fzy(fzy_lo(1):fzy_hi(1),fzy_lo(2):fzy_hi(2),fzy_lo(3):fzy_hi(3),NVAR)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                rrr = qp_core(i,j,k,QRHO)
                compr = rrr*qp_pass(i,j,k,nqs)
                rrnewr = rrr - cdtdy*(fyz(i,j+1,k,URHO) - fyz(i,j,k,URHO)) &
                             - cdtdz*(fzy(i,j  ,k+1,URHO) - fzy(i,j,k,URHO))
                compnr = compr - cdtdy*(fyz(i,j+1,k,n   ) - fyz(i,j,k,n)) &
                               - cdtdz*(fzy(i,j  ,k+1,n   ) - fzy(i,j,k,n))

                qpo_pass(i  ,j,k,nqs) = compnr/rrnewr

                rrl = qm_core(i,j,k,QRHO)
                compl = rrl*qm_pass(i,j,k,nqs)
                rrnewl = rrl - cdtdy*(fyz(i-1,j+1,k,URHO) - fyz(i-1,j,k,URHO)) &
                             - cdtdz*(fzy(i-1,j  ,k+1,URHO) - fzy(i-1,j,k,URHO))
                compnl = compl - cdtdy*(fyz(i-1,j+1,k,n   ) - fyz(i-1,j,k,n)) &
                               - cdtdz*(fzy(i-1,j  ,k+1,n   ) - fzy(i-1,j,k,n))

                qmo_pass(i,j,k,nqs) = compnl/rrnewl
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
             rrr = qp_core(i,j,k,QRHO)
             rur = rrr*qp_core(i,j,k,QU)
             rvr = rrr*qp_core(i,j,k,QV)
             rwr = rrr*qp_core(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp_core(i,j,k,QU:QW)**2)
             rer = qp_core(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp_rad(i,j,k,qrad:qradhi)
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

             qpo_core(i,j,k,QRHO  ) = rrnewr
             qpo_core(i,j,k,QU    ) = runewr/rrnewr
             qpo_core(i,j,k,QV    ) = rvnewr/rrnewr
             qpo_core(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo_core(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo_Core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) &
                        - cdtdy*(fyz(i,j+1,k,UEINT) - fyz(i,j,k,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,k+1,UEINT) - fzy(i,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp_core(i,j,k,QPRES) - pynew - pznew
                   qpo_core(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo_core(i,j,k,QGAME) = qp_core(i,j,k,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo_core(i,j,k,QPRES) = qpo_core(i,j,k,QREINT)*(qpo_core(i,j,k,QGAME)-ONE)
                end if
             else
                qpo_core(i,j,k,QPRES) = qp_core(i,j,k,QPRES)
                qpo_core(i,j,k,QGAME) = qp_core(i,j,k,QGAME)
             endif

             qpo_Core(i,j,k,QPRES) = max(qpo_core(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo_core, qcpo_lo, qcpo_hi, &
                                          qpo_pass, qppo_lo, qppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qpo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo_core(i,j,k,QPRES)
             qpo_rad(i,j,k,qreitot) = sum(qpo_rad(i,j,k,qrad:qradhi)) + qpo_core(i,j,k,QREINT)
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
             rrl = qm_core(i,j,k,QRHO)
             rul = rrl*qm_core(i,j,k,QU)
             rvl = rrl*qm_core(i,j,k,QV)
             rwl = rrl*qm_core(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm_core(i,j,k,QU:QW)**2)
             rel = qm_core(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm_rad(i,j,k,qrad:qradhi)
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

             qmo_core(i,j,k,QRHO   ) = rrnewl
             qmo_core(i,j,k,QU     ) = runewl/rrnewl
             qmo_core(i,j,k,QV     ) = rvnewl/rrnewl
             qmo_core(i,j,k,QW     ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo_core(i,j,k,QREINT ) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo_core(i,j,k,QREINT ) = qm_core(i,j,k,QREINT) &
                        - cdtdy*(fyz(i-1,j+1,k,UEINT) - fyz(i-1,j,k,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i-1,j  ,k+1,UEINT) - fzy(i-1,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm_core(i,j,k,QPRES) - pynew - pznew
                   qmo_core(i,j,k,QPRES  ) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo_core(i,j,k,QGAME) = qm_core(i,j,k,QGAME) + geynew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo_core(i,j,k,QPRES) = qmo_core(i,j,k,QREINT)*(qmo_core(i,j,k,QGAME)-ONE)
                end if
             else
                qmo_core(i,j,k,QPRES) = qm_core(i,j,k,QPRES)
                qmo_core(i,j,k,QGAME) = qm_core(i,j,k,QGAME)
             endif

             qmo_core(i,j,k,QPRES) = max(qmo_core(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo_core, qcmo_lo, qcmo_hi, &
                                          qmo_pass, qpmo_lo, qpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qmo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qmo_rad(i,j,k,qptot) = sum(lambda(:)*ernewl(:)) + qmo_core(i,j,k,QPRES)
             qmo_rad(i,j,k,qreitot) = sum(qmo_rad(i,j,k,qrad:qradhi)) + qmo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transyz

  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz(lo, hi, &
                     qm_core, qcm_lo, qcm_hi, &
                     qmo_core, qcmo_lo, qcmo_hi, &
                     qm_pass, qpm_lo, qpm_hi, &
                     qmo_pass, qpmo_lo, qpmo_hi, &
#ifdef RADIATION
                     qm_rad, qrm_lo, qrm_hi, &
                     qmo_rad, qrmo_lo, qrmo_hi, &
#endif
                     qp_core, qcp_lo, qcp_hi, &
                     qpo_core, qcpo_lo, qcpo_hi, &
                     qp_pass, qpp_lo, qpp_hi, &
                     qpo_pass, qppo_lo, qppo_hi, &
#ifdef RADIATION
                     qp_rad, qrp_lo, qrp_hi, &
                     qpo_rad, qrpo_lo, qrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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


    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcmo_lo(3), qcmo_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)
    integer, intent(in) :: qcpo_lo(3), qcpo_hi(3)
    integer, intent(in) :: qpm_lo(3), qpm_hi(3)
    integer, intent(in) :: qpmo_lo(3), qpmo_hi(3)
    integer, intent(in) :: qpp_lo(3), qpp_hi(3)
    integer, intent(in) :: qppo_lo(3), qppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qrm_lo(3), qrm_hi(3)
    integer, intent(in) :: qrmo_lo(3), qrmo_hi(3)
    integer, intent(in) :: qrp_lo(3), qrp_hi(3)
    integer, intent(in) :: qrpo_lo(3), qrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(in) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)
    real(rt), intent(in) :: qm_pass(qpm_lo(1):qpm_hi(1),qpm_lo(2):qpm_hi(2),qpm_lo(3):qpm_hi(3),NQP)
    real(rt), intent(in) :: qp_pass(qpp_lo(1):qpp_hi(1),qpp_lo(2):qpp_hi(2),qpp_lo(3):qpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qm_rad(qrm_lo(1):qrm_hi(1),qrm_lo(2):qrm_hi(2),qrm_lo(3):qrm_hi(3),NQR)
    real(rt), intent(in) ::  qp_rad(qrp_lo(1):qrp_hi(1),qrp_lo(2):qrp_hi(2),qrp_lo(3):qrp_hi(3),NQR)
#endif
    real(rt), intent(out) :: qmo_core(qcmo_lo(1):qcmo_hi(1),qcmo_lo(2):qcmo_hi(2),qcmo_lo(3):qcmo_hi(3),NQC)
    real(rt), intent(out) :: qpo_core(qcpo_lo(1):qcpo_hi(1),qcpo_lo(2):qcpo_hi(2),qcpo_lo(3):qcpo_hi(3),NQC)
    real(rt), intent(out) :: qmo_pass(qpmo_lo(1):qpmo_hi(1),qpmo_lo(2):qpmo_hi(2),qpmo_lo(3):qpmo_hi(3),NQP)
    real(rt), intent(out) :: qpo_pass(qppo_lo(1):qppo_hi(1),qppo_lo(2):qppo_hi(2),qppo_lo(3):qppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qmo_rad(qrmo_lo(1):qrmo_hi(1),qrmo_lo(2):qrmo_hi(2),qrmo_lo(3):qrmo_hi(3),NQR)
    real(rt), intent(out) ::  qpo_rad(qrpo_lo(1):qrpo_hi(1),qrpo_lo(2):qrpo_hi(2),qrpo_lo(3):qrpo_hi(3),NQR)
#endif

    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxz(fxz_lo(1):fxz_hi(1),fxz_lo(2):fxz_hi(2),fxz_lo(3):fxz_hi(3),NVAR)
    real(rt), intent(in) :: fzx(fzx_lo(1):fzx_hi(1),fzx_lo(2):fzx_hi(2),fzx_lo(3):fzx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrr = qp_core(i,j,k,QRHO)
                compr = rrr*qp_pass(i,j,k,nqs)
                rrnewr = rrr - cdtdx*(fxz(i+1,j,k,URHO) - fxz(i,j,k,URHO)) &
                             - cdtdz*(fzx(i  ,j,k+1,URHO) - fzx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxz(i+1,j,k,n) - fxz(i,j,k,n)) &
                               - cdtdz*(fzx(i  ,j,k+1,n) - fzx(i,j,k,n))

                qpo_pass(i,j  ,k,nqs) = compnr/rrnewr

                rrl = qm_core(i,j,k,QRHO)
                compl = rrl*qm_pass(i,j,k,nqs)
                rrnewl = rrl - cdtdx*(fxz(i+1,j-1,k,URHO) - fxz(i,j-1,k,URHO)) &
                             - cdtdz*(fzx(i  ,j-1,k+1,URHO) - fzx(i,j-1,k,URHO))
                compnl = compl - cdtdx*(fxz(i+1,j-1,k,n) - fxz(i,j-1,k,n)) &
                               - cdtdz*(fzx(i  ,j-1,k+1,n) - fzx(i,j-1,k,n))

                qmo_pass(i,j,k,nqs) = compnl/rrnewl
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
             rrr = qp_core(i,j,k,QRHO)
             rur = rrr*qp_core(i,j,k,QU)
             rvr = rrr*qp_core(i,j,k,QV)
             rwr = rrr*qp_core(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp_core(i,j,k,QU:QW)**2)
             rer = qp_core(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp_rad(i,j,k,qrad:qradhi)
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
             qpo_core(i,j,k,QRHO  ) = rrnewr
             qpo_core(i,j,k,QU    ) = runewr/rrnewr
             qpo_core(i,j,k,QV    ) = rvnewr/rrnewr
             qpo_core(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo_core(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) &
                        - cdtdx*(fxz(i+1,j,k,UEINT) - fxz(i,j,k,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i  ,j,k+1,UEINT) - fzx(i,j,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp_core(i,j,k,QPRES) - pxnew - pznew
                   qpo_core(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo_core(i,j,k,QGAME) = qp_core(i,j,k,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qpo_core(i,j,k,QPRES) = qpo_core(i,j,k,QREINT)*(qpo_core(i,j,k,QGAME)-ONE)
                endif
             else
                qpo_core(i,j,k,QPRES) = qp_core(i,j,k,QPRES)
                qpo_core(i,j,k,QGAME) = qp_core(i,j,k,QGAME)
             endif

             qpo_core(i,j,k,QPRES) = max(qpo_core(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo_core, qcpo_lo, qcpo_hi, &
                                          qpo_pass, qppo_lo, qppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qpo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo_core(i,j,k,QPRES)
             qpo_rad(i,j,k,qreitot) = sum(qpo_rad(i,j,k,qrad:qradhi)) + qpo_core(i,j,k,QREINT)
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
             rrl = qm_core(i,j,k,QRHO)
             rul = rrl*qm_core(i,j,k,QU)
             rvl = rrl*qm_core(i,j,k,QV)
             rwl = rrl*qm_core(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm_core(i,j,k,QU:QW)**2)
             rel = qm_core(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm_rad(i,j,k,qrad:qradhi)
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
             qmo_core(i,j,k,QRHO  ) = rrnewl
             qmo_core(i,j,k,QU    ) = runewl/rrnewl
             qmo_core(i,j,k,QV    ) = rvnewl/rrnewl
             qmo_core(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo_core(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo_core(i,j,k,QREINT) = qm_core(i,j,k,QREINT) &
                        - cdtdx*(fxz(i+1,j-1,k,UEINT) - fxz(i,j-1,k,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i,j-1,k+1,UEINT) - fzx(i,j-1,k,UEINT) + pzav*duz)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm_core(i,j,k,QPRES) - pxnew - pznew
                   qmo_core(i,j,k,QPRES) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo_core(i,j,k,QGAME) = qm_core(i,j,k,QGAME) + gexnew + geznew

                   ! and compute the p edge state from this and (rho e)
                   qmo_core(i,j,k,QPRES) = qmo_core(i,j,k,QREINT)*(qmo_core(i,j,k,QGAME)-ONE)
                endif
             else
                qmo_core(i,j,k,QPRES) = qm_core(i,j,k,QPRES)
                qmo_core(i,j,k,QGAME) = qm_core(i,j,k,QGAME)
             endif

             qmo_core(i,j,k,QPRES) = max(qmo_core(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo_core, qcmo_lo, qcmo_hi, &
                                          qmo_pass, qpmo_lo, qpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qmo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qmo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo_core(i,j,k,QPRES)
             qmo_rad(i,j,k,qreitot) = sum(qmo_rad(i,j,k,qrad:qradhi)) + qmo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxz

  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy(lo, hi, &
                     qm_core, qcm_lo, qcm_hi, &
                     qmo_core, qcmo_lo, qcmo_hi, &
                     qm_pass, qpm_lo, qpm_hi, &
                     qmo_pass, qpmo_lo, qpmo_hi, &
#ifdef RADIATION
                     qm_rad, qrm_lo, qrm_hi, &
                     qmo_rad, qrmo_lo, qrmo_hi, &
#endif
                     qp_core, qcp_lo, qcp_hi, &
                     qpo_core, qcpo_lo, qcpo_hi, &
                     qp_pass, qpp_lo, qpp_hi, &
                     qpo_pass, qppo_lo, qppo_hi, &
#ifdef RADIATION
                     qp_rad, qrp_lo, qrp_hi, &
                     qpo_rad, qrpo_lo, qrpo_hi, &
#endif
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
    use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QFS, QFX, &
                                   QC, QGAMC, &
#ifdef RADIATION
                                   NQR, qrad, qradhi, qptot, qreitot, &
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


    integer, intent(in) :: qcm_lo(3), qcm_hi(3)
    integer, intent(in) :: qcmo_lo(3), qcmo_hi(3)
    integer, intent(in) :: qcp_lo(3), qcp_hi(3)
    integer, intent(in) :: qcpo_lo(3), qcpo_hi(3)
    integer, intent(in) :: qpm_lo(3), qpm_hi(3)
    integer, intent(in) :: qpmo_lo(3), qpmo_hi(3)
    integer, intent(in) :: qpp_lo(3), qpp_hi(3)
    integer, intent(in) :: qppo_lo(3), qppo_hi(3)
#ifdef RADIATION
    integer, intent(in) :: qrm_lo(3), qrm_hi(3)
    integer, intent(in) :: qrmo_lo(3), qrmo_hi(3)
    integer, intent(in) :: qrp_lo(3), qrp_hi(3)
    integer, intent(in) :: qrpo_lo(3), qrpo_hi(3)
#endif
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

    real(rt), intent(in) :: qm_core(qcm_lo(1):qcm_hi(1),qcm_lo(2):qcm_hi(2),qcm_lo(3):qcm_hi(3),NQC)
    real(rt), intent(in) :: qp_core(qcp_lo(1):qcp_hi(1),qcp_lo(2):qcp_hi(2),qcp_lo(3):qcp_hi(3),NQC)
    real(rt), intent(in) :: qm_pass(qpm_lo(1):qpm_hi(1),qpm_lo(2):qpm_hi(2),qpm_lo(3):qpm_hi(3),NQP)
    real(rt), intent(in) :: qp_pass(qpp_lo(1):qpp_hi(1),qpp_lo(2):qpp_hi(2),qpp_lo(3):qpp_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(in) ::  qm_rad(qrm_lo(1):qrm_hi(1),qrm_lo(2):qrm_hi(2),qrm_lo(3):qrm_hi(3),NQR)
    real(rt), intent(in) ::  qp_rad(qrp_lo(1):qrp_hi(1),qrp_lo(2):qrp_hi(2),qrp_lo(3):qrp_hi(3),NQR)
#endif
    real(rt), intent(out) :: qmo_core(qcmo_lo(1):qcmo_hi(1),qcmo_lo(2):qcmo_hi(2),qcmo_lo(3):qcmo_hi(3),NQC)
    real(rt), intent(out) :: qpo_core(qcpo_lo(1):qcpo_hi(1),qcpo_lo(2):qcpo_hi(2),qcpo_lo(3):qcpo_hi(3),NQC)
    real(rt), intent(out) :: qmo_pass(qpmo_lo(1):qpmo_hi(1),qpmo_lo(2):qpmo_hi(2),qpmo_lo(3):qpmo_hi(3),NQP)
    real(rt), intent(out) :: qpo_pass(qppo_lo(1):qppo_hi(1),qppo_lo(2):qppo_hi(2),qppo_lo(3):qppo_hi(3),NQP)
#ifdef RADIATION
    real(rt), intent(out) ::  qmo_rad(qrmo_lo(1):qrmo_hi(1),qrmo_lo(2):qrmo_hi(2),qrmo_lo(3):qrmo_hi(3),NQR)
    real(rt), intent(out) ::  qpo_rad(qrpo_lo(1):qrpo_hi(1),qrpo_lo(2):qrpo_hi(2),qrpo_lo(3):qrpo_hi(3),NQR)
#endif
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(in) :: fxy(fxy_lo(1):fxy_hi(1),fxy_lo(2):fxy_hi(2),fxy_lo(3):fxy_hi(3),NVAR)
    real(rt), intent(in) :: fyx(fyx_lo(1):fyx_hi(1),fyx_lo(2):fyx_hi(2),fyx_lo(3):fyx_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)

    integer i, j, k, n, nqs, ipassive

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
       nqs = qpass_map(ipassive)

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                rrr = qp_core(i,j,k,QRHO)
                compr = rrr*qp_pass(i,j,k,nqs)
                rrnewr = rrr - cdtdx*(fxy(i+1,j,k,URHO) - fxy(i,j,k,URHO)) &
                             - cdtdy*(fyx(i,j+1,k,URHO) - fyx(i,j,k,URHO))
                compnr = compr - cdtdx*(fxy(i+1,j,k,n) - fxy(i,j,k,n)) &
                               - cdtdy*(fyx(i,j+1,k,n) - fyx(i,j,k,n))

                qpo_pass(i,j,k,nqs) = compnr/rrnewr

                rrl = qm_core(i,j,k,QRHO)
                compl = rrl*qm_pass(i,j,k,nqs)
                rrnewl = rrl - cdtdx*(fxy(i+1,j,k-1,URHO) - fxy(i,j,k-1,URHO)) &
                             - cdtdy*(fyx(i,j+1,k-1,URHO) - fyx(i,j,k-1,URHO))
                compnl = compl - cdtdx*(fxy(i+1,j,k-1,n) - fxy(i,j,k-1,n)) &
                               - cdtdy*(fyx(i,j+1,k-1,n) - fyx(i,j,k-1,n))

                qmo_pass(i,j,k,nqs) = compnl/rrnewl
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
             rrr = qp_core(i,j,k,QRHO)
             rur = rrr*qp_core(i,j,k,QU)
             rvr = rrr*qp_core(i,j,k,QV)
             rwr = rrr*qp_core(i,j,k,QW)
             ekenr = HALF*rrr*sum(qp_core(i,j,k,QU:QW)**2)
             rer = qp_core(i,j,k,QREINT) + ekenr
#ifdef RADIATION
             err = qp_rad(i,j,k,qrad:qradhi)
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
             qpo_core(i,j,k,QRHO  ) = rrnewr
             qpo_core(i,j,k,QU    ) = runewr/rrnewr
             qpo_core(i,j,k,QV    ) = rvnewr/rrnewr
             qpo_core(i,j,k,QW    ) = rwnewr/rrnewr

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo_core(i,j,k,QREINT) = renewr - rhoekenr

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating
                   ! (rho e).
                   qpo_core(i,j,k,QREINT) = qp_core(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k,UEINT) - fxy(i,j,k,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k,UEINT) - fyx(i,j,k,UEINT) + pyav*duy)
                endif


                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp_core(i,j,k,QPRES) - pxnew - pynew
                   qpo_core(i,j,k,QPRES) = pnewr
                else
                   ! Update gammae with its transverse terms
                   qpo_core(i,j,k,QGAME) = qp_core(i,j,k,QGAME) + gexnew + geynew

                   ! and compute the p edge state from this and (rho e)
                   qpo_core(i,j,k,QPRES) = qpo_core(i,j,k,QREINT)*(qpo_core(i,j,k,QGAME)-ONE)
                endif
             else
                qpo_core(i,j,k,QPRES) = qp_core(i,j,k,QPRES)
                qpo_core(i,j,k,QGAME) = qp_core(i,j,k,QGAME)
             endif

             qpo_core(i,j,k,QPRES) = max(qpo_core(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qpo_core, qcpo_lo, qcpo_hi, &
                                          qpo_pass, qppo_lo, qppo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qpo_rad(i,j,k,qrad:qradhi) = ernewr(:)
             qpo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewr(:)) + qpo_core(i,j,k,QPRES)
             qpo_rad(i,j,k,qreitot) = sum(qpo_rad(i,j,k,qrad:qradhi)) + qpo_core(i,j,k,QREINT)
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
             rrl = qm_core(i,j,k,QRHO)
             rul = rrl*qm_core(i,j,k,QU)
             rvl = rrl*qm_core(i,j,k,QV)
             rwl = rrl*qm_core(i,j,k,QW)
             ekenl = HALF*rrl*sum(qm_core(i,j,k,QU:QW)**2)
             rel = qm_core(i,j,k,QREINT) + ekenl
#ifdef RADIATION
             erl = qm_rad(i,j,k,qrad:qradhi)
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
             qmo_core(i,j,k,QRHO  ) = rrnewl
             qmo_core(i,j,k,QU    ) = runewl/rrnewl
             qmo_core(i,j,k,QV    ) = rvnewl/rrnewl
             qmo_core(i,j,k,QW    ) = rwnewl/rrnewl

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo_core(i,j,k,QREINT) = renewl - rhoekenl

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo_core(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).
                   qmo_core(i,j,k,QREINT) = qm_core(i,j,k,QREINT) &
                        - cdtdx*(fxy(i+1,j,k-1,UEINT) - fxy(i,j,k-1,UEINT) + pxav*dux) &
                        - cdtdy*(fyx(i,j+1,k-1,UEINT) - fyx(i,j,k-1,UEINT) + pyav*duy)
                endif

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm_core(i,j,k,QPRES) - pxnew - pynew
                   qmo_core(i,j,k,QPRES) = pnewl
                else
                   ! Update gammae with its transverse terms
                   qmo_core(i,j,k,QGAME) = qm_core(i,j,k,QGAME) + gexnew + geynew

                   ! and compute the p edge state from this and (rho e)
                   qmo_core(i,j,k,QPRES) = qmo_core(i,j,k,QREINT)*(qmo_core(i,j,k,QGAME)-ONE)
                endif
             else
                qmo_core(i,j,k,QPRES) = qm_core(i,j,k,QPRES)
                qmo_core(i,j,k,QGAME) = qm_core(i,j,k,QGAME)
             endif

             qmo_core(i,j,k,QPRES) = max(qmo_core(i,j,k,QPRES), small_pres)

             call reset_edge_state_thermo(qmo_core, qcmo_lo, qcmo_hi, &
                                          qmo_pass, qpmo_lo, qpmo_hi, &
                                          i, j, k)

#ifdef RADIATION
             qmo_rad(i,j,k,qrad:qradhi) = ernewl(:)
             qmo_rad(i,j,k,qptot  ) = sum(lambda(:)*ernewl(:)) + qmo_core(i,j,k,QPRES)
             qmo_rad(i,j,k,qreitot) = sum(qmo_rad(i,j,k,qrad:qradhi)) + qmo_core(i,j,k,QREINT)
#endif

          end do
       end do
    end do

  end subroutine transxy
#endif

  subroutine reset_edge_state_thermo(qedge_core, qc_lo, qc_hi, &
                                     qedge_pass, qp_lo, qp_hi, &
                                     ii, jj, kk)

  use amrex_constants_module, only : ZERO, ONE, HALF

  use network, only : nspec, naux
  use meth_params_module, only : NQC, NQP, NVAR, NQAUX, QRHO, &
                                 QPRES, QREINT, QGAME, QFS, QFX, &
                                 small_pres, small_temp, &
                                 ppm_predict_gammae, &
                                 transverse_use_eos, transverse_reset_rhoe

  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer, intent(in) :: ii, jj, kk
    integer, intent(in) :: qc_lo(3), qc_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    real(rt), intent(inout) :: qedge_core(qc_lo(1):qc_hi(1),qc_lo(2):qc_hi(2),qc_lo(3):qc_hi(3),NQC)
    real(rt), intent(inout) :: qedge_pass(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQP)

    logical :: reset
    type (eos_t) :: eos_state

    !$gpu

    reset = .false.

    if (transverse_reset_rhoe == 1) then
       ! if we are still negative, then we need to reset
       if (qedge_core(ii,jj,kk,QREINT) < ZERO) then
          reset = .true.

          eos_state % rho = qedge_core(ii,jj,kk,QRHO)
          eos_state % T = small_temp
          eos_state % xn(:) = qedge_pass(ii,jj,kk,QFS:QFS-1+nspec)
          eos_state % aux(:) = qedge_pass(ii,jj,kk,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          qedge_core(ii,jj,kk,QREINT) = qedge_core(ii,jj,kk,QRHO)*eos_state % e
          qedge_core(ii,jj,kk,QPRES) = eos_state % p
       endif

    end if

    if (ppm_predict_gammae == 0 ) then

       if (transverse_use_eos == 1) then
          eos_state % rho = qedge_core(ii,jj,kk,QRHO)
          eos_state % e   = qedge_core(ii,jj,kk,QREINT) / qedge_core(ii,jj,kk,QRHO)
          eos_state % T   = small_temp
          eos_state % xn  = qedge_pass(ii,jj,kk,QFS:QFS+nspec-1)
          eos_state % aux = qedge_pass(ii,jj,kk,QFX:QFX+naux-1)

          call eos(eos_input_re, eos_state)

          qedge_core(ii,jj,kk,QREINT) = eos_state % e * eos_state % rho
          qedge_core(ii,jj,kk,QPRES) = max(eos_state % p, small_pres)
       end if

    else
       if (reset) then
          ! recompute the p edge state from this and (rho e), since we reset
          ! qreint  (actually, is this code even necessary?)
          qedge_core(ii,jj,kk,QPRES) = qedge_core(ii,jj,kk,QREINT)*(qedge_core(ii,jj,kk,QGAME)-ONE)
          qedge_core(ii,jj,kk,QPRES) = max(qedge_core(ii,jj,kk,QPRES), small_pres)
       end if
    end if

  end subroutine reset_edge_state_thermo

end module transverse_module
