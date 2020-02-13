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
    real(rt), intent(in) :: flux_t(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)
    real(rt), intent(in) :: q_t(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    real(rt), intent(out) :: qmo(qmo_lo(1):qmo_hi(1),qmo_lo(2):qmo_hi(2),qmo_lo(3):qmo_hi(3),NQ)
    real(rt), intent(out) :: qpo(qpo_lo(1):qpo_hi(1),qpo_lo(2):qpo_hi(2),qpo_lo(3):qpo_hi(3),NQ)
#if AMREX_SPACEDIM == 2
    real(rt), intent(in) :: a_t(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
#endif

    integer i, j, k, n, nqp, ipassive

    real(rt) :: rhoinv
    real(rt) :: rrnew
    real(rt) :: rrry, rrly
    real(rt) :: rury, ruly
    real(rt) :: rvry, rvly
    real(rt) :: rwry, rwly
    real(rt) :: ekenry, ekenly
    real(rt) :: rery, rely
    real(rt) :: rrnewry, rrnewly
    real(rt) :: runewry, runewly
    real(rt) :: rvnewry, rvnewly
    real(rt) :: rwnewry, rwnewly
    real(rt) :: renewry, renewly
    real(rt) :: pnewry, pnewly
    real(rt) :: rhoekenry, rhoekenly, rhoekenrz, rhoekenlz
    real(rt) :: pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    real(rt) :: compu
    real(rt) :: :: gamc

#ifdef RADIATION
    real(rt) :: :: dre, dmom, divu
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, ernewr, ernewl, &
         lamge, luge, der
    real(rt) :: eddf, f1
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

             !----------------------------------------------------------------
             ! qpo state
             !----------------------------------------------------------------

             ! We are handling the states at the interface of
             ! (i, i+1) in the x-direction, and similarly for
             ! the y- and z- directions.

             il = i
             jl = j
             kl = k

             if (idir_t == 1) then
                ir = i+1
                jr = j
                kr = k
             else if (idir_t == 2) then
                ir = i
                jr = j+1
                kr = k
             else
                ir = i
                jr = j
                kr = k+1
             end if

             ! the "p" interface state in a local array
             lqn(:) = qp(i,j,k,:)

             do ipassive = 1, npassive
                n  = upass_map(ipassive)
                nqp = qpass_map(ipassive)

#if AMREX_SPACEDIM == 2
                volinv = ONE/vol(i,j,k)

                rrnew = lqn(QRHO) - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,URHO) - &
                                         area_t(il,jl,kl)*flux_t(il,jl,kl,URHO)) * volinv
                compu = lqn(QRHO)*lqn(nqp) - &
                     hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,n) - &
                          area_t(il,jl,kl)*flux_t(il,jl,kl,n)) * volinv
                lqno(nqp) = compu/rrnew
#else
                rrnew = lqn(QRHO) - cdtdx*(flux_t(ir,jr,kr,URHO) - flux_t(il,jl,kl,URHO))
                compu = lqn(QRHO)*lqn(nqp) - cdtdx*(flux_t(ir,jr,kr,n) - flux_t(il,jl,kl,n))
                lqno(nqp) = compu/rrnew
#endif
             end do

#if AMREX_SPACEDIM == 2
             volinv = ONE/vol(i,j,k)
#endif

             pgp  = q_t(ir,jr,kr,GDPRES)
             pgm  = q_t(il,jl,kl,GDPRES)
             ugp  = q_t(ir,jr,kr,GDU-1+idir_t)
             ugm  = q_t(il,jl,kl,GDU-1+dir_t)
             gegp = q_t(ir,jr,kr,GDGAME)
             gegm = q_t(il,jl,kl,GDGAME)

#ifdef RADIATION
             lambda(:) = qaux(il,jl,kl,QLAMS:QLAMS+ngroups-1)
             ergp(:) = q_t(ir,jr,kr,GDERADS:GDERADS-1+ngroups)
             ergm(:) = q_t(il,jl,kl,GDERADS:GDERADS-1+ngroups)
#endif

             ! we need to augment our conserved system with either a p
             ! equation or gammae (if we have ppm_predict_gammae = 1) to
             ! be able to deal with the general EOS

#if AMREX_SPACEDIM == 2
             dup = area_t(ir,jr,kr)*pgp*ugp - area_t(il,jl,kl)*pgm*ugm
             du = area_t(ir,jr,kr)*ugp-area_t(il,jl,kl)*ugm
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
             gamc = qaux(il,jl,kl,QGAMCG)
#else
             gamc = qaux(il,jl,kl,QGAMC)
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
                divu = (area_t(ir,jr,kr)*ugp - area_t(il,jl,kl)*ugm) * volinv
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
             rrn = lqn(QRHO)
             run = rrrn*lqn(QU)
             rvn = rrrn*lqn(QV)
             rwn = rrrn*lqn(QW)
             ekenn = HALF*rrrn*sum(lqn(QU:QW)**2)
             ren = lqn(QREINT) + ekenrn
#ifdef RADIATION
             ern = lqn(qrad:qrad-1+ngroups)
#endif

#if AMREX_SPACEDIM == 2
             ! Add transverse predictor
             rrnewn = rrn - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,URHO) -  &
                                 area_t(i,j,k)*flux_t(i,j,k,URHO)) * volinv

             ! Note that pressure may be treated specially here, depending on
             ! the geometry.  Our y-interface equation for (rho u) is:
             !
             !  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
             !
             ! in cylindrical coords -- note that the p term is not in
             ! a divergence, so there are no area factors.  For this
             ! geometry, we do not include p in our definition of the
             ! flux in the x-direction, for we need to fix this now.
             runewn = run - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UMX) -  &
                                 area_t(il,jl,kl)*flux_t(il,jl,kl,UMX)) * volinv
             if (.not. mom_flux_has_p(idir_t)%comp(UMX)) then
                runewn = runewn - cdtdx * (pgp-pgm)
             endif
             rvnewn = rvn - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UMY) -  &
                                 area_t(il,jl,kl)*flux_t(il,jl,kl,UMY)) * volinv
             rwnewn = rwn - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UMZ) -  &
                                 area_t(il,jl,kl)*flux_t(il,jl,kl,UMZ)) * volinv
             renewn = ren - hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UEDEN) -  &
                                 area_t(il,jl,kl)*flux_t(il,jl,kl,UEDEN)) * volinv

#ifdef RADIATION
             runewn = runewn - HALF*hdt*(area_t(ir,jr,kr)+area_t(il,jl,kl))*sum(lamge) * volinv
             renewn = renewn + dre
             ernewn(:) = err(:) - hdt*(area_t(ir,jr,kr)*rflux_t(ir,jr,kr,:) -  &
                                       area_t(il,jl,kl)*rflux_t(il,jl,kl,:)) * volinv + der(:)
#endif

#else
             ! Add transverse predictor
             rrnewn = rrn - cdtdx*(flux_t(ir,jr,kr,URHO) - flux_t(il,jl,kl,URHO))
             runewn = run - cdtdx*(flux_t(ir,jr,kr,UMX) - flux_t(il,jl,kl,UMX))
             rvnewn = rvn - cdtdx*(flux_t(ir,jr,kr,UMY) - flux_t(il,jl,kl,UMY))
             rwnewn = rwn - cdtdx*(flux_t(ir,jr,kr,UMZ) - flux_t(il,jl,kl,UMZ))
             renewn = ren - cdtdx*(flux_t(ir,jr,kr,UEDEN) - flux_t(il,jl,kl,UEDEN))
#ifdef RADIATION
             runewn = runewn + dmom
             renewn = renewn + dre
             ernewn  = err(:) - cdtdx*(rflux_t(ir,jr,kr,:) - rflux_t(il,jl,kl,:)) + der(:)
#endif
#endif

             ! Reset to original value if adding transverse terms made density negative
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewry < ZERO) then
                rrnewn = rrn
                runewn = run
                rvnewn = rvn
                rwnewn = rwn
                renewn = ren
#ifdef RADIATION
                ernewn = ern(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             lqno(QRHO) = rrnewn
             rhoinv = ONE/rrnewn
             lqno(QU) = runewn*rhoinv
             lqno(QV) = rvnewn*rhoinv
             lqno(QW) = rwnewn*rhoinv

             ! note: we run the risk of (rho e) being negative here
             rhoekenn = HALF*(runewn**2 + rvnewn**2 + rwnewn**2)*rhoinv
             lqno(QREINT) = renewn - rhoekenn

             if (.not. reset_state) then
                ! do the transverse terms for p, gamma, and rhoe, as necessary

                if (transverse_reset_rhoe == 1 .and. qpo(i,j,k,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
#if AMREX_SPACEDIM == 2
                   lqno(QREINT) = lqn(QREINT) - &
                        hdt*(area_t(ir,jr,kr)*flux_t(ir,jr,kr,UEINT) - &
                             area_t(il,jl,kl)*flux_t(il,jl,kl,UEINT) + pav*du) * volinv
#else
                   lqno(QREINT) = lqn(QREINT) - &
                        cdtdx*(flux_t(ir,jr,kr,UEINT) - flux_t(il,jl,kl,UEINT) + pav*du)
#endif
                end if

                ! Pretend QREINT has been fixed and transverse_use_eos .ne. 1.
                ! If we are wrong, we will fix it later

                if (ppm_predict_gammae == 0) then
                   ! add the transverse term to the p evolution eq here
#if AMREX_SPACEDIM == 2
                   ! the divergences here, dup and du, already have area factors
                   pnewn = lqn(QPRES) - hdt*(dup + pav*du*(gamc - ONE)) * volinv
#else
                   pnewn = lqn(QPRES) - cdtdx*(dup + pav*du*(gamc - ONE))
#endif
                   lqno(QPRES) = max(pnewn, small_pres)
                else
                   ! Update gammae with its transverse terms
#if AMREX_SPACEDIM == 2
                   lqno(QGAME) = lqn(QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*du) * volinv - cdtdx*uav*dge
#else
                   lqno(QGAME) = lqn(QGAME) + &
                        cdtdx*( (geav-ONE)*(geav - gamc)*du - uav*dge )
#endif
                   ! and compute the p edge state from this and (rho e)
                   lqno(QPRES) = lqno(QREINT)*(lqno(QGAME) - ONE)
                   lqno(QPRES) = max(lqno(QPRES), small_pres)
                end if
             else
                lqno(QPRES) = lqn(QPRES)
                lqno(QGAME) = lqn(QGAME)
             endif

#ifdef RADIATION
             lqno(qrad:qrad-1+ngroups) = ernewn(:)
             lqno(qptot  ) = sum(lambda(:)*ernewn(:)) + lqno(QPRES)
             lqno(qreitot) = sum(lqno(qrad:qrad-1+ngroups)) + lqno(QREINT)
#endif

             ! store the state
             qpo(i,j,k,:) = lqno(:)


             !-------------------------------------------------------------------
             ! qmo state
             !-------------------------------------------------------------------

             do ipassive = 1, npassive
                n  = upass_map(ipassive)
                nqp = qpass_map(ipassive)

#if AMREX_SPACEDIM == 2
                volinv = ONE/vol(i,j-1,k)

                rrnew = qm(i,j,k,QRHO) - hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,URHO) - &
                                               area_t(i,j-1,k)*flux_t(i,j-1,k,URHO)) * volinv
                compu = qm(i,j,k,QRHO)*qm(i,j,k,nqp) - &
                     hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,n) - &
                          area_t(i,j-1,k)*flux_t(i,j-1,k,n)) * volinv
                qmo(i,j,k,nqp) = compu/rrnew
#else
                rrnew = qm(i,j,k,QRHO) - cdtdx*(flux_t(i+1,j-1,k,URHO) - flux_t(i,j-1,k,URHO))
                compu = qm(i,j,k,QRHO)*qm(i,j,k,nqp) - cdtdx*(flux_t(i+1,j-1,k,n) - flux_t(i,j-1,k,n))
                qmo(i,j,k,nqp) = compu/rrnew
#endif
             end do


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
             rrnewly = rrly - hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,URHO) -  &
                                   area_t(i,j-1,k)*flux_t(i,j-1,k,URHO)) * volinv
             runewly = ruly - hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,UMX)  -  &
                                   area_t(i,j-1,k)*flux_t(i,j-1,k,UMX)) * volinv
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewly = runewly - cdtdx *(pgp-pgm)
             endif
             rvnewly = rvly - hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,UMY)  -  &
                                   area_t(i,j-1,k)*flux_t(i,j-1,k,UMY)) * volinv
             rwnewly = rwly - hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,UMZ)  -  &
                                   area_t(i,j-1,k)*flux_t(i,j-1,k,UMZ)) * volinv
             renewly = rely - hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,UEDEN)-  &
                                   area_t(i,j-1,k)*flux_t(i,j-1,k,UEDEN)) * volinv

#ifdef RADIATION
             runewly = runewly - HALF*hdt*(area_t(i+1,j-1,k)+area_t(i,j-1,k))*sum(lamge) * volinv
             renewly = renewly + dre
             ernewl(:) = erl(:) - hdt*(area_t(i+1,j-1,k)*rflux_t(i+1,j-1,k,:)-  &
                                       area_t(i,j-1,k)*rflux_t(i,j-1,k,:)) * volinv + der(:)
#endif

#else
             ! Add transverse predictor
             rrnewly = rrly - cdtdx*(flux_t(i+1,j-1,k,URHO) - flux_t(i,j-1,k,URHO))
             runewly = ruly - cdtdx*(flux_t(i+1,j-1,k,UMX) - flux_t(i,j-1,k,UMX))
             rvnewly = rvly - cdtdx*(flux_t(i+1,j-1,k,UMY) - flux_t(i,j-1,k,UMY))
             rwnewly = rwly - cdtdx*(flux_t(i+1,j-1,k,UMZ) - flux_t(i,j-1,k,UMZ))
             renewly = rely - cdtdx*(flux_t(i+1,j-1,k,UEDEN) - flux_t(i,j-1,k,UEDEN))
#ifdef RADIATION
             runewly = runewly + dmom
             renewly = renewly + dre
             ernewl  = erl(:) - cdtdx*(rflux_t(i+1,j-1,k,:) - rflux_t(i,j-1,k,:)) + der(:)
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
                        hdt*(area_t(i+1,j-1,k)*flux_t(i+1,j-1,k,UEINT)-  &
                             area_t(i,j-1,k)*flux_t(i,j-1,k,UEINT) + pav*du) * volinv
#else
                   qmo(i,j,k,QREINT) = qm(i,j,k,QREINT) - &
                        cdtdx*(flux_t(i+1,j-1,k,UEINT) - flux_t(i,j-1,k,UEINT) + pav*du)
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
