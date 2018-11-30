module transverse_module

  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH, TWO
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public transx, transy

contains

  subroutine transx(qm, qmo, qp, qpo, qd_lo, qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    fx, fx_lo, fx_hi, &
#ifdef RADIATION
                    rfx, rfx_lo, rfx_hi, &
#endif
                    qgdx, qgdx_lo, qgdx_hi, &
                    srcQ, src_lo, src_hi, &
                    hdt, cdtdx,  &
                    area1, area1_lo, area1_hi, &
                    vol, vol_lo, vol_hi, &
                    ilo, ihi, jlo, jhi)

    use amrex_constants_module
    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NQAUX, &
                                 NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QGAME, &
                                 URHO, UMX, UMY, UEDEN, UEINT, QFS, QFX, &
                                 GDU, GDV, GDPRES, GDGAME, &
                                 NGDNV, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 GDERADS, QGAMCG, QLAMS, &
                                 fspace_type, comoving, &
#endif
                                 small_pres, small_temp, &
                                 npassive, qpass_map, upass_map, &
                                 transverse_use_eos, ppm_type, &
                                 transverse_reset_density, transverse_reset_rhoe, &
                                 ppm_predict_gammae
    use prob_params_module, only : mom_flux_has_p
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer qd_lo(3), qd_hi(3)
    integer qa_lo(3), qa_hi(3)
    integer fx_lo(3), fx_hi(3)
    integer qgdx_lo(3), qgdx_hi(3)
    integer src_lo(3), src_hi(3)
    integer area1_lo(3), area1_hi(3)
    integer vol_lo(3), vol_hi(3)
    integer ilo, ihi, jlo, jhi

#ifdef RADIATION
    integer rfx_lo(3), rfx_hi(3)
    real(rt)         rfx(rfx_lo(1):rfx_hi(1),rfx_lo(2):rfx_hi(2),0:ngroups-1)
#endif

    real(rt)         qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)

    real(rt)         fx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),NVAR)
    real(rt)         qgdx(qgdx_lo(1):qgdx_hi(1),qgdx_lo(2):qgdx_hi(2),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),QVAR)
    real(rt)         area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
    real(rt)         vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))

    real(rt)         hdt, cdtdx

    integer          :: i, j, g
    integer          :: n, nqp, ipassive

    real(rt)         :: rr, rrnew, compo, compn
    real(rt)         :: rrr, rur, rvr, rer, ekinr, rhoekinr
    real(rt)         :: rrnewr, runewr, rvnewr, renewr
    real(rt)         :: rrl, rul, rvl, rel, ekinl, rhoekinl
    real(rt)         :: rrnewl, runewl, rvnewl, renewl

    ! here, pggp/pggm is the Godunov gas pressure (not radiation contribution)
    real(rt)         :: pggp, pggm, ugp, ugm, dAup, pav, uav, dAu, pnewl,pnewr
    real(rt)         :: geav, dge, gegp, gegm, gamc
    real(rt)         :: rhotmp

#ifdef RADIATION
    real(rt)         dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    real(rt)         eddf, f1, ugc, divu
#endif

    type (eos_t) :: eos_state

    logical :: reset_state

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the x-direction to y-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    ! here cdtdx = 0.5 dt / dx

    ! update all of the passively-advected quantities in a single loop
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do j = jlo, jhi
          do i = ilo, ihi

             compn = hdt*(area1(i+1,j)*fx(i+1,j,n) - &
                          area1(i  ,j)*fx(i  ,j,n))/vol(i,j)

             if (j >= jlo+1) then
                rr = qp(i,j,  QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                  area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)
                compo = rr*qp(i,j,nqp) - compn
                qpo(i,j,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if

             if (j <= jhi-1) then
                rr = qm(i,j+1,QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                  area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)
                compo = rr*qm(i,j+1,nqp) - compn
                qmo(i,j+1,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if
          enddo
       enddo
    enddo

    ! hydro variables
    do j = jlo, jhi
       do i = ilo, ihi

          pggp = qgdx(i+1,j,GDPRES)
          pggm = qgdx(i,  j,GDPRES)
          ugp = qgdx(i+1,j,GDU)
          ugm = qgdx(i,  j,GDU)
          gegp = qgdx(i+1,j,GDGAME)
          gegm = qgdx(i,  j,GDGAME)

#ifdef RADIATION
          lambda(:) = qaux(i,j,QLAMS:QLAMS+ngroups-1)
          ugc = 0.5e0_rt*(ugp+ugm)
          ergp(:) = qgdx(i+1,j,GDERADS:GDERADS-1+ngroups)
          ergm(:) = qgdx(i  ,j,GDERADS:GDERADS-1+ngroups)
#endif

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dAup = area1(i+1,j)*pggp*ugp - area1(i,j)*pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          dAu = area1(i+1,j)*ugp-area1(i,j)*ugm
          geav = HALF*(gegp+gegm)
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,QGAMCG)
#else
          gamc = qaux(i,j,QGAMC)
#endif

#ifdef RADIATION
          lamge(:) = lambda(:) * (ergp(:)-ergm(:))
          luge(:) = ugc * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             divu = (area1(i+1,j)*ugp-area1(i,j)*ugm)/vol(i,j)
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = -hdt * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * divu
             end do
          else ! mixed frame
             der(:) = cdtdx * luge(:)
          end if
#endif

          !-------------------------------------------------------------------
          ! qp state
          !-------------------------------------------------------------------

          ! "right" state on the j-1/2 interface
          if (j >= jlo+1) then

             ! Convert to conservation form
             rrr = qp(i,j,QRHO)
             rur = rrr*qp(i,j,QU)
             rvr = rrr*qp(i,j,QV)
             ekinr = HALF*rrr*sum(qp(i,j,QU:QW)**2)
             rer = qp(i,j,QREINT) + ekinr
#ifdef RADIATION
             err(:) = qp(i,j,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                 area1(i,j)*fx(i,j,URHO))/vol(i,j)

             ! Note that pressure may be treated specially here, depending on 
             ! the geometry.  Our y-interface equation for (rho u) is:
             !
             !  d(rho u)/dt + d(rho u v)/dy = - 1/r d(r rho u u)/dr - dp/dr
             !
             ! in cylindrical coords -- note that the p term is not in
             ! a divergence, so there are no area factors.  For this
             ! geometry, we do not include p in our definition of the
             ! flux in the x-direction, for we need to fix this now.
             runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                 area1(i,j)*fx(i,j,UMX))/vol(i,j)
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewr = runewr -cdtdx *(pggp-pggm)
             endif
             rvnewr = rvr - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                                 area1(i,j)*fx(i,j,UMY))/vol(i,j)
             renewr = rer - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                                 area1(i,j)*fx(i,j,UEDEN))/vol(i,j)

#ifdef RADIATION
             runewr = runewr - HALF*hdt*(area1(i+1,j)+area1(i,j))*sum(lamge)/vol(i,j)
             renewr = renewr + dre
             ernewr(:) = err(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
                  area1(i,j)*rfx(i,j,:))/vol(i,j) &
                  + der(:)
#endif

             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr < ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                renewr = rer
#ifdef RADIATION
                ernewr(:) = err(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             rhotmp = rrnewr
             qpo(i,j,QRHO) = rhotmp
             qpo(i,j,QU  ) = runewr/rhotmp
             qpo(i,j,QV  ) = rvnewr/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qpo(i,j,QRHO) = qpo(i,j,QRHO) + hdt*srcQ(i,j,QRHO)
                qpo(i,j,QU:QV) = qpo(i,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekinr = HALF*(runewr**2+rvnewr**2+(rhotmp*qpo(i,j,QW))**2)/rhotmp
             qpo(i,j,QREINT) = renewr - rhoekinr

             if (ppm_type == 0) then
                qpo(i,j,QREINT) = qpo(i,j,QREINT) + hdt*srcQ(i,j,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by
                   ! using the discretized expression for updating (rho e).
                   qpo(i,j,QREINT) = qp(i,j,QREINT) - &
                        hdt*(area1(i+1,j)*fx(i+1,j,UEINT)-  &
                             area1(i,j)*fx(i,j,UEINT) + pav*dAu)/vol(i,j)

                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,QREINT) < ZERO .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qpo(i,j,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qpo(i,j,QREINT) = qpo(i,j,QRHO)*eos_state % e
                      qpo(i,j,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % e   = qpo(i,j,QREINT) / qpo(i,j,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qpo(i,j,QFS:QFS+nspec-1)
                      eos_state % aux = qpo(i,j,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qpo(i,j,QPRES ) = pnewr
                      qpo(i,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      ! we are expressing the pressure evolution as:
                      !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
                      ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
                      ! but these are divergences, so we need area factors
                      pnewr = qp(i,j,QPRES) - hdt*(dAup + pav*dAu*(gamc - ONE))/vol(i,j)
                      qpo(i,j,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qpo(i,j,QGAME) = qp(i,j,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*dAu)/vol(i,j) - cdtdx*uav*dge

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,QPRES) = qpo(i,j,QREINT)*(qpo(i,j,QGAME)-ONE)

                endif
             else
                qpo(i,j,QPRES) = qp(i,j,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qpo(i,j,QGAME) = qp(i,j,QGAME)
             endif

#ifdef RADIATION
             qpo(i,j,qrad:qradhi) = ernewr(:)
             qpo(i,j,qptot)   = sum(lambda*ernewr) + qpo(i,j,QPRES)
             qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qm state
          !-------------------------------------------------------------------

          ! "left" state on the j+1/2 interface
          if (j <= jhi-1) then

             rrl = qm(i,j+1,QRHO)
             rul = rrl*qm(i,j+1,QU)
             rvl = rrl*qm(i,j+1,QV)
             ekinl = HALF*rrl*sum(qm(i,j+1,QU:QW)**2)
             rel = qm(i,j+1,QREINT) + ekinl
#ifdef RADIATION
             erl(:) = qm(i,j+1,qrad:qradhi)
#endif

             rrnewl = rrl - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                 area1(i,j)*fx(i,j,URHO))/vol(i,j)
             runewl = rul - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                 area1(i,j)*fx(i,j,UMX))/vol(i,j) 
             if (.not. mom_flux_has_p(1)%comp(UMX)) then
                runewl = runewl -cdtdx *(pggp-pggm)
             endif
             rvnewl = rvl - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                                 area1(i,j)*fx(i,j,UMY))/vol(i,j)
             renewl = rel - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                                 area1(i,j)*fx(i,j,UEDEN))/vol(i,j)


#ifdef RADIATION
             runewl = runewl - HALF*hdt*(area1(i+1,j)+area1(i,j))*sum(lamge)/vol(i,j)
             renewl = renewl + dre
             ernewl(:) = erl(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
                  area1(i,j)*rfx(i,j,:))/vol(i,j) &
                  + der(:)
#endif
             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl < ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                renewl = rel
#ifdef RADIATION
                ernewl(:) = erl(:)
#endif
                reset_state = .true.
             endif

             ! Convert back to primitive form
             rhotmp = rrnewl
             qmo(i,j+1,QRHO) = rhotmp
             qmo(i,j+1,QU  ) = runewl/rhotmp
             qmo(i,j+1,QV  ) = rvnewl/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qmo(i,j+1,QRHO)  = qmo(i,j+1,QRHO) + hdt*srcQ(i,j,QRHO)
                qmo(i,j+1,QU:QV) = qmo(i,j+1,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinl = HALF*(runewl**2+rvnewl**2+(rhotmp*qmo(i,j+1,QW))**2)/rhotmp
             qmo(i,j+1,QREINT) = renewl - rhoekinl
             if (ppm_type == 0) then
                qmo(i,j+1,QREINT) = qmo(i,j+1,QREINT) + hdt*srcQ(i,j,QREINT)
             endif


             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i,j+1,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using
                   ! the discretized expression for updating (rho e).

                   qmo(i,j+1,QREINT) = qm(i,j+1,QREINT) - &
                        hdt*(area1(i+1,j)*fx(i+1,j,UEINT)-  &
                             area1(i,j)*fx(i,j,UEINT) + pav*dAu)/vol(i,j)

                   ! if we are still negative, then we need to reset
                   if (qmo(i,j+1,QREINT) < ZERO .and. qmo(i,j+1,QRHO) > ZERO) then
                      eos_state % rho = qmo(i,j+1,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i,j+1,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qmo(i,j+1,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qmo(i,j+1,QREINT) = qmo(i,j+1,QRHO)*eos_state % e
                      qmo(i,j+1,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qmo(i,j+1,QRHO) > ZERO) then
                      eos_state % rho = qmo(i,j+1,QRHO)
                      eos_state % e   = qmo(i,j+1,QREINT) / qmo(i,j+1,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qmo(i,j+1,QFS:QFS+nspec-1)
                      eos_state % aux = qmo(i,j+1,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qmo(i,j+1,QPRES ) = pnewr
                      qmo(i,j+1,QREINT) = eos_state % e * eos_state % rho
                   else
                      ! we are expressing the pressure evolution as:
                      !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
                      ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
                      ! but these are divergences, so we need area factors
                      pnewl = qm(i,j+1,QPRES) - hdt*(dAup + pav*dAu*(gamc - ONE))/vol(i,j)
                      qmo(i,j+1,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i,j+1,QPRES) = qmo(i,j+1,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qmo(i,j+1,QPRES) = max(qmo(i,j+1,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qmo(i,j+1,QGAME) = qm(i,j+1,QGAME) + &
                        hdt*( (geav-ONE)*(geav - gamc)*dAu)/vol(i,j) - cdtdx*uav*dge

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j+1,QPRES) = qmo(i,j+1,QREINT)*(qmo(i,j+1,QGAME)-ONE)

                endif
             else
                qmo(i,j+1,QPRES) = qm(i,j+1,QPRES)
                if (ppm_type == 0) then
                   qmo(i,j+1,QPRES) = qmo(i,j+1,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qmo(i,j+1,QGAME) = qm(i,j+1,QGAME)
             endif

#ifdef RADIATION
             qmo(i,j+1,qrad:qradhi) = ernewl(:)
             qmo(i,j+1,qptot)   = sum(lambda*ernewl) + qmo(i,j+1,QPRES)
             qmo(i,j+1,qreitot) = sum(qmo(i,j+1,qrad:qradhi)) + qmo(i,j+1,QREINT)
#endif

          end if

       enddo
    enddo

  end subroutine transx


  subroutine transy(qm, qmo, qp, qpo, qd_lo, qd_hi, &
                    qaux, qa_lo, qa_hi, &
                    fy, fy_lo, fy_hi, &
#ifdef RADIATION
                    rfy, rfy_lo, rfy_hi, &
#endif
                    qgdy, qgdy_lo, qgdy_hi, &
                    srcQ, src_lo, src_hi, &
                    hdt, cdtdy, ilo, ihi, jlo, jhi)

    use amrex_constants_module
    use network, only : nspec, naux
    use meth_params_module, only : NQ, QVAR, NQAUX, &
                                 NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QGAME, &
                                 URHO, UMX, UMY, UEDEN, UEINT, QFS, QFX, &
                                 GDU, GDV, GDPRES, GDGAME, &
                                 NGDNV, QGAMC, &
#ifdef RADIATION
                                 qrad, qradhi, qptot, qreitot, &
                                 GDERADS, QGAMCG, QLAMS, &
                                 fspace_type, comoving, &
#endif
                                 small_pres, small_temp, &
                                 npassive, qpass_map, upass_map, &
                                 transverse_use_eos, ppm_type, &
                                 transverse_reset_density, transverse_reset_rhoe, &
                                 ppm_predict_gammae
    use prob_params_module, only : mom_flux_has_p
#ifdef RADIATION
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t

    integer qd_lo(3), qd_hi(3)
    integer qa_lo(3), qa_hi(3)
    integer fy_lo(3), fy_hi(3)
    integer qgdy_lo(3), qgdy_hi(3)
    integer src_lo(3), src_hi(3)
    integer ilo, ihi, jlo, jhi

#ifdef RADIATION
    integer rfy_lo(3), rfy_hi(3)
    real(rt)         rfy(rfy_lo(1):rfy_hi(1),rfy_lo(2):rfy_hi(2),0:ngroups-1)
#endif

    real(rt)         qm(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qmo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qp(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)
    real(rt)         qpo(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),NQ)

    real(rt)         qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)

    real(rt)         fy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),NVAR)
    real(rt)         qgdy(qgdy_lo(1):qgdy_hi(1),qgdy_lo(2):qgdy_hi(2),NGDNV)
    real(rt)         srcQ(src_lo(1):src_hi(1),src_lo(2):src_hi(2),QVAR)

    real(rt)         hdt, cdtdy

    integer          :: i, j, g
    integer          :: n, nqp, ipassive

    real(rt)         :: rr,rrnew
    real(rt)         :: pggp, pggm, ugp, ugm, dup, pav, uav, du, pnewr,pnewl
    real(rt)         :: gegp, gegm, geav, dge, gamc
    real(rt)         :: rrr, rur, rvr, rer, ekinr, rhoekinr
    real(rt)         :: rrnewr, runewr, rvnewr, renewr
    real(rt)         :: rrl, rul, rvl, rel, ekinl, rhoekinl
    real(rt)         :: rrnewl, runewl, rvnewl, renewl
    real(rt)         :: rhotmp
    real(rt)         :: compo, compn

#ifdef RADIATION
    real(rt)         :: dre, dmom
    real(rt)        , dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    real(rt)         :: eddf, f1, ugc
#endif

    type (eos_t) :: eos_state

    logical :: reset_state

    !-------------------------------------------------------------------
    ! add the transverse flux difference in the y-direction to x-states
    ! for the fluid variables
    !-------------------------------------------------------------------

    ! update all of the passively-advected quantities in a single loop
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nqp = qpass_map(ipassive)

       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdy*(fy(i,j+1,n)-fy(i,j,n))

             if (i >= ilo+1) then
                rr = qp(i,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
                compo = rr*qp(i,j,nqp) - compn
                qpo(i,j,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if

             if (i <= ihi-1) then
                rr = qm(i+1,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
                compo = rr*qm(i+1,j,nqp) - compn
                qmo(i+1,j,nqp) = compo/rrnew + hdt*srcQ(i,j,nqp)
             end if
          enddo
       enddo
    enddo

    ! hydro variables
    do j = jlo, jhi
       do i = ilo, ihi

          pggp = qgdy(i,j+1,GDPRES)
          pggm = qgdy(i,j  ,GDPRES)
          ugp = qgdy(i,j+1,GDV)
          ugm = qgdy(i,j  ,GDV)
          gegp = qgdy(i,j+1,GDGAME)
          gegm = qgdy(i,j  ,GDGAME)

#ifdef RADIATION
          lambda(:) = qaux(i,j,QLAMS:QLAMS+ngroups-1)
          ugc = 0.5e0_rt*(ugp+ugm)
          ergp(:) = qgdy(i,j+1,GDERADS:GDERADS-1+ngroups)
          ergm(:) = qgdy(i,j  ,GDERADS:GDERADS-1+ngroups)
#endif
          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pggp*ugp - pggm*ugm
          pav = HALF*(pggp+pggm)
          uav = HALF*(ugp+ugm)
          du = ugp-ugm
          geav = HALF*(gegp+gegm)
          dge = gegp-gegm

          ! this is the gas gamma_1
#ifdef RADIATION
          gamc = qaux(i,j,QGAMCG)
#else
          gamc = qaux(i,j,QGAMC)
#endif


#ifdef RADIATION
          lamge(:) = lambda(:) * (ergp(:)-ergm(:))
          luge(:) = ugc * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                der(g) = cdtdy * f1 * 0.5e0_rt*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! qp state
          !-------------------------------------------------------------------

          ! right state on the i-1/2 interface
          if (i >= ilo+1) then

             ! Convert to conservation form
             rrr = qp(i,j,QRHO)
             rur = rrr*qp(i,j,QU)
             rvr = rrr*qp(i,j,QV)
             ekinr = HALF*rrr*sum(qp(i,j,QU:QW)**2)
             rer = qp(i,j,QREINT) + ekinr
#ifdef RADIATION
             err(:) = qp(i,j,qrad:qradhi)
#endif

             ! Add transverse predictor
             rrnewr = rrr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))

             runewr = rur - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
             ! note: we are always Cartesian in the y-direction, so the
             ! pressure term is already accounted for in the flux
             rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) 
             renewr = rer - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN))

#ifdef RADIATION
             rvnewr = rvnewr - cdtdy*sum(lamge)
             renewr = renewr + dre
             ernewr(:) = err(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
                  + der(:)
#endif

             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewr <= ZERO) then
                rrnewr = rrr
                runewr = rur
                rvnewr = rvr
                renewr = rer
#ifdef RADIATION
                ernewr(:) = err(:)
#endif
                reset_state = .true.
             end if

             ! convert back to non-conservation form
             rhotmp =  rrnewr
             qpo(i,j,QRHO  ) = rhotmp
             qpo(i,j,QU    ) = runewr/rhotmp
             qpo(i,j,QV    ) = rvnewr/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qpo(i,j,QRHO  ) = qpo(i,j,QRHO  ) + hdt*srcQ(i,j,QRHO)
                qpo(i,j,QU:QV) = qpo(i,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinr = HALF*(runewr**2+rvnewr**2+(rhotmp*qpo(i,j,QW))**2)/rhotmp
             qpo(i,j,QREINT) = renewr - rhoekinr
             if (ppm_type == 0) then
                qpo(i,j,QREINT) = qpo(i,j,QREINT) + hdt*srcQ(i,j,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,QREINT) <= ZERO) then
                   qpo(i,j,QREINT) = qp(i,j,QREINT) - &
                        cdtdy*(fy(i,j+1,UEINT)- fy(i,j,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,QREINT) < ZERO .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qpo(i,j,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qpo(i,j,QREINT) = qpo(i,j,QRHO) * eos_state % e
                      qpo(i,j,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % e   = qpo(i,j,QREINT) / qpo(i,j,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qpo(i,j,QFS:QFS+nspec-1)
                      eos_state % aux = qpo(i,j,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qpo(i,j,QPRES ) = pnewr
                      qpo(i,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      pnewr = qp(i  ,j,QPRES)-cdtdy*(dup + pav*du*(gamc - ONE))
                      qpo(i,j,QPRES) = pnewr
                      if (ppm_type == 0) then
                         qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qpo(i,j,QGAME) = qp(i,j,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,QPRES) = qpo(i,j,QREINT)*(qpo(i,j,QGAME)-ONE)

                endif
             else
                qpo(i,j,QPRES) = qp(i,j,QPRES)
                if (ppm_type == 0) then
                   qpo(i,j,QPRES) = qpo(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qpo(i,j,QGAME) = qp(i,j,QGAME)
             endif

#ifdef RADIATION
             qpo(i,j,qrad:qradhi) = ernewr(:)
             qpo(i,j,qptot  ) = sum(lambda*ernewr) + qpo(i,j,QPRES)
             qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)
#endif

          end if

          !-------------------------------------------------------------------
          ! qm state
          !-------------------------------------------------------------------

          ! left state on the i+1/2 interface
          if (i <= ihi-1) then

             rrl = qm(i+1,j,QRHO)
             rul = rrl*qm(i+1,j,QU)
             rvl = rrl*qm(i+1,j,QV)
             ekinl = HALF*rrl*sum(qm(i+1,j,QU:QW)**2)
             rel = qm(i+1,j,QREINT) + ekinl
#ifdef RADIATION
             erl(:) = qm(i+1,j,qrad:qradhi)
#endif

             rrnewl = rrl - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))
             runewl = rul - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
             rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) 
             renewl = rel - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN))

#ifdef RADIATION
             rvnewl = rvnewl - cdtdy*sum(lamge)
             renewl = renewl + dre
             ernewl(:) = erl(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
                  + der(:)
#endif

             reset_state = .false.
             if (transverse_reset_density == 1 .and. rrnewl <= ZERO) then
                rrnewl = rrl
                runewl = rul
                rvnewl = rvl
                renewl = rel
#ifdef RADIATION
                ernewl(:) = erl(:)
#endif
                reset_state = .true.
             endif

             rhotmp =  rrnewl
             qmo(i+1,j,QRHO  ) = rhotmp            
             qmo(i+1,j,QU    ) = runewl/rhotmp
             qmo(i+1,j,QV    ) = rvnewl/rhotmp

             ! for ppm_type > 0 we already added the piecewise
             ! parabolic traced source terms to the normal edge states
             if (ppm_type == 0) then
                qmo(i+1,j,QRHO  ) = qmo(i+1,j,QRHO  ) + hdt*srcQ(i,j,QRHO)
                qmo(i+1,j,QU:QV) = qmo(i+1,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinl = HALF*(runewl**2+rvnewl**2+(rhotmp*qmo(i+1,j,QW))**2)/rhotmp
             qmo(i+1,j,QREINT) = renewl - rhoekinl
             if (ppm_type == 0) then
                qmo(i+1,j,QREINT) = qmo(i+1,j,QREINT) + hdt*srcQ(i,j,QREINT)
             endif

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qmo(i+1,j,QREINT) <= ZERO) then
                   ! If it is negative, reset the internal energy by using the discretized
                   ! expression for updating (rho e).

                   qmo(i+1,j,QREINT) = qm(i+1,j,QREINT) - &
                        cdtdy*(fy(i,j+1,UEINT) - fy(i,j,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qmo(i+1,j,QREINT) < ZERO .and. qmo(i+1,j,QRHO) > ZERO) then
                      eos_state % rho = qmo(i+1,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i+1,j,QFS:QFS-1+nspec)
                      eos_state % aux(:) = qmo(i+1,j,QFX:QFX-1+naux)

                      call eos(eos_input_rt, eos_state)

                      qmo(i+1,j,QREINT) = qmo(i+1,j,QRHO)*eos_state % e
                      qmo(i+1,j,QPRES) = eos_state % p
                   endif
                endif

                if (ppm_predict_gammae == 0) then

                   ! Optionally, use the EOS to calculate the pressure.

                   if (transverse_use_eos .eq. 1 .and. qmo(i+1,j,QRHO) > ZERO) then
                      eos_state % rho = qmo(i+1,j,QRHO)
                      eos_state % e   = qmo(i+1,j,QREINT) / qmo(i+1,j,QRHO)
                      eos_state % T   = small_temp
                      eos_state % xn  = qmo(i+1,j,QFS:QFS+nspec-1)
                      eos_state % aux = qmo(i+1,j,QFX:QFX+naux-1)

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qmo(i+1,j,QPRES ) = pnewr
                      qmo(i+1,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      pnewl = qm(i+1,j,QPRES)-cdtdy*(dup + pav*du*(gamc - ONE))
                      qmo(i+1,j,QPRES) = pnewl
                      if (ppm_type == 0) then
                         qmo(i+1,j,QPRES) = qmo(i+1,j,QPRES) + hdt*srcQ(i,j,QPRES)
                      endif
                   endif

                   qmo(i+1,j,QPRES) = max(qmo(i+1,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qmo(i+1,j,QGAME) = qm(i+1,j,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav - gamc)*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qmo(i+1,j,QPRES) = qmo(i+1,j,QREINT)*(qmo(i+1,j,QGAME)-ONE)

                endif
             else
                qmo(i+1,j,QPRES) = qm(i+1,j,QPRES)
                if (ppm_type == 0) then
                   qmo(i+1,j,QPRES) = qmo(i+1,j,QPRES) + hdt*srcQ(i,j,QPRES)
                endif
                qmo(i+1,j,QGAME) = qm(i+1,j,QGAME)
             endif

#ifdef RADIATION
             qmo(i+1,j,qrad:qradhi) = ernewl(:)
             qmo(i+1,j,qptot  ) = sum(lambda*ernewl) + qmo(i+1,j,QPRES)
             qmo(i+1,j,qreitot) = sum(qmo(i+1,j,qrad:qradhi)) + qmo(i+1,j,QREINT)
#endif

          end if

       enddo
    enddo

  end subroutine transy

end module transverse_module
