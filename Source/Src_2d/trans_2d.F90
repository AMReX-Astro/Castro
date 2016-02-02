module transverse_module

  use bl_constants_module
  use network, only : nspec
  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QGAME, &
                                 URHO, UMX, UMY, UEDEN, UEINT, QFS, &
                                 GDU, GDV, GDPRES, GDGAME, ngdnv, &
                                 small_pres, small_temp, &
                                 npassive, qpass_map, upass_map, &
                                 transverse_use_eos, ppm_type, ppm_trace_sources, &
                                 transverse_reset_density, transverse_reset_rhoe, &
                                 ppm_predict_gammae
#ifdef RADIATION
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
                                     fspace_type, comoving
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
#endif

  use eos_module

  implicit none

  private

  public :: transx, transy

contains

  subroutine transx( &
#ifdef RADIATION
                    lam, lam_l1, lam_l2, lam_h1, lam_h2, &
#endif
                    qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                    fx, fx_l1, fx_l2, fx_h1, fx_h2, &
#ifdef RADIATION
                    rfx, rfx_l1, rfx_l2, rfx_h1, rfx_h2, &
#endif
                    qgdx, qgdx_l1, qgdx_l2, qgdx_h1, qgdx_h2, &
                    gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                    srcQ, src_l1, src_l2, src_h1, src_h2, &
                    hdt, cdtdx,  &
                    area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                    vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                    ilo, ihi, jlo, jhi)

#ifdef RADIATION
    integer lam_l1,lam_l2,lam_h1,lam_h2
    integer rfx_l1, rfx_l2, rfx_h1, rfx_h2
#endif
    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer gc_l1, gc_l2, gc_h1, gc_h2
    integer fx_l1, fx_l2, fx_h1, fx_h2
    integer qgdx_l1, qgdx_l2, qgdx_h1, qgdx_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer area1_l1, area1_l2, area1_h1, area1_h2
    integer vol_l1, vol_l2, vol_h1, vol_h2
    integer ilo, ihi, jlo, jhi

#ifdef RADIATION
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision rfx(rfx_l1:rfx_h1,rfx_l2:rfx_h2,0:ngroups-1)
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
#else
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
#endif

    double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,NVAR)
    double precision qgdx(qgdx_l1:qgdx_h1,qgdx_l2:qgdx_h2,ngdnv)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
    double precision hdt, cdtdx

    integer          :: i, j, g
    integer          :: n, nq, ipassive

    double precision :: rr, rrnew, compo, compn
    double precision :: rrr, rur, rvr, rer, ekinr, rhoekinr
    double precision :: rrnewr, runewr, rvnewr, renewr
    double precision :: rrl, rul, rvl, rel, ekinl, rhoekinl
    double precision :: rrnewl, runewl, rvnewl, renewl

    ! here, pggp/pggm is the Godunov gas pressure (not radiation contribution)
    double precision :: pggp, pggm, ugp, ugm, dAup, pav, uav, dAu, pnewl,pnewr
    double precision :: geav, dge, gegp, gegm
    double precision :: rhotmp

#ifdef RADIATION
    double precision dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    double precision eddf, f1, ugc, divu
#endif

    type (eos_t) :: eos_state

    logical :: reset_state

    eos_state % check_small = .false.

    ! update all of the passively-advected quantities in a single loop
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do j = jlo, jhi
          do i = ilo, ihi

             compn = hdt*(area1(i+1,j)*fx(i+1,j,n) - &
                          area1(i  ,j)*fx(i  ,j,n))/vol(i,j)

             if (j >= jlo+1) then
                rr = qp(i,j,  QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                  area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)
                compo = rr*qp(i,j,nq) - compn
                qpo(i,j,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
             end if

             if (j <= jhi-1) then
                rr = qm(i,j+1,QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                  area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)
                compo = rr*qm(i,j+1,nq) - compn
                qmo(i,j+1,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
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
          lambda(:) = lam(i,j,:)
          ugc = 0.5d0*(ugp+ugm)
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

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to y-states
          ! for the fluid variables
          !-------------------------------------------------------------------

#ifdef RADIATION
          lamge(:) = lambda(:) * (ergp(:)-ergm(:))
          luge(:) = ugc * lamge(:)
          dre = -cdtdx*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             divu = (area1(i+1,j)*ugp-area1(i,j)*ugm)/vol(i,j)
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = -hdt * f1 * 0.5d0*(ergp(g)+ergm(g)) * divu
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
             runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                 area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                      -HALF*hdt*(area1(i+1,j)+area1(i,j))*(pggp-pggm)/vol(i,j)
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
             qpo(i,j,QRHO) = rhotmp        + hdt*srcQ(i,j,QRHO)
             qpo(i,j,QU  ) = runewr/rhotmp
             qpo(i,j,QV  ) = rvnewr/rhotmp

             ! if ppm_trace_sources == 1, then we already added the
             ! piecewise parabolic traced source terms to the normal edge
             ! states
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qpo(i,j,QU:QV) = qpo(i,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             ! note: we run the risk of (rho e) being negative here
             rhoekinr = HALF*(runewr**2+rvnewr**2+(rhotmp*qpo(i,j,QW))**2)/rhotmp
             qpo(i,j,QREINT) = renewr - rhoekinr + hdt*srcQ(i,j,QREINT)

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

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qpo(i,j,QPRES ) = pnewr
                      qpo(i,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      ! we are expressing the pressure evolution as:
                      !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
                      ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
                      ! but these are divergences, so we need area factors
                      pnewr = qp(i,j,QPRES) - hdt*(dAup + pav*dAu*(gamc(i,j)-ONE))/vol(i,j)
                      qpo(i,j,QPRES) = pnewr + hdt*srcQ(i,j,QPRES)
                   endif

                   qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qpo(i,j,QGAME) = qp(i,j,QGAME) + &
                        hdt*( (geav-ONE)*(geav-gamc(i,j))*dAu)/vol(i,j) - cdtdx*uav*dge

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,QPRES) = qpo(i,j,QREINT)*(qpo(i,j,QGAME)-ONE)

                endif
             else
                qpo(i,j,QPRES) = qp(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
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
                                 area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                      -HALF*hdt*(area1(i+1,j)+area1(i,j))*(pggp-pggm)/vol(i,j)
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
             qmo(i,j+1,QRHO) = rhotmp         + hdt*srcQ(i,j,QRHO)
             qmo(i,j+1,QU  ) = runewl/rhotmp
             qmo(i,j+1,QV  ) = rvnewl/rhotmp

             ! if ppm_trace_sources == 1, then we already added the
             ! piecewise parabolic traced source terms to the normal edge
             ! states
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qmo(i,j+1,QU:QV) = qmo(i,j+1,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinl = HALF*(runewl**2+rvnewl**2+(rhotmp*qmo(i,j+1,QW))**2)/rhotmp
             qmo(i,j+1,QREINT)= renewl - rhoekinl +hdt*srcQ(i,j,QREINT)

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

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qmo(i,j+1,QPRES ) = pnewr
                      qmo(i,j+1,QREINT) = eos_state % e * eos_state % rho
                   else
                      ! we are expressing the pressure evolution as:
                      !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
                      ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
                      ! but these are divergences, so we need area factors
                      pnewl = qm(i,j+1,QPRES) - hdt*(dAup + pav*dAu*(gamc(i,j)-ONE))/vol(i,j)
                      qmo(i,j+1,QPRES) = pnewl + hdt*srcQ(i,j,QPRES)
                   endif

                   qmo(i,j+1,QPRES) = max(qmo(i,j+1,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qmo(i,j+1,QGAME) = qm(i,j+1,QGAME) + &
                        hdt*( (geav-ONE)*(geav-gamc(i,j))*dAu)/vol(i,j) - cdtdx*uav*dge

                   ! and compute the p edge state from this and (rho e)
                   qmo(i,j+1,QPRES) = qmo(i,j+1,QREINT)*(qmo(i,j+1,QGAME)-ONE)

                endif
             else
                qmo(i,j+1,QPRES) = qm(i,j+1,QPRES) + hdt*srcQ(i,j,QPRES)
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


  subroutine transy( &
#ifdef RADIATION
                    lam, lam_l1, lam_l2, lam_h1, lam_h2, &
#endif
                    qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                    fy,fy_l1,fy_l2,fy_h1,fy_h2, &
#ifdef RADIATION
                    rfy,rfy_l1,rfy_l2,rfy_h1,rfy_h2, &
#endif
                    qgdy, qgdy_l1, qgdy_l2, qgdy_h1, qgdy_h2, &
                    gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                    srcQ, src_l1, src_l2, src_h1, src_h2, &
                    hdt, cdtdy, ilo, ihi, jlo, jhi)

#ifdef RADIATION
    integer lam_l1,lam_l2,lam_h1,lam_h2
    integer rfy_l1, rfy_l2, rfy_h1, rfy_h2
#endif

    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer gc_l1, gc_l2, gc_h1, gc_h2
    integer fy_l1, fy_l2, fy_h1, fy_h2
    integer qgdy_l1, qgdy_l2, qgdy_h1, qgdy_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer ilo, ihi, jlo, jhi

#ifdef RADIATION
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
#else
    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision rfy(rfy_l1:rfy_h1,rfy_l2:rfy_h2,0:ngroups-1)
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
#endif
    double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,NVAR)
    double precision qgdy(qgdy_l1:qgdy_h1,qgdy_l2:qgdy_h2,ngdnv)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision hdt, cdtdy

    integer          :: i, j, g
    integer          :: n, nq, ipassive

    double precision :: rr,rrnew
    double precision :: pggp, pggm, ugp, ugm, dup, pav, uav, du, pnewr,pnewl
    double precision :: gegp, gegm, geav, dge
    double precision :: rrr, rur, rvr, rer, ekinr, rhoekinr
    double precision :: rrnewr, runewr, rvnewr, renewr
    double precision :: rrl, rul, rvl, rel, ekinl, rhoekinl
    double precision :: rrnewl, runewl, rvnewl, renewl
    double precision :: rhotmp
    double precision :: compo, compn

#ifdef RADIATION
    double precision :: dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    double precision :: eddf, f1, ugc
#endif

    type (eos_t) :: eos_state

    logical :: reset_state

    eos_state % check_small = .false.


    ! update all of the passively-advected quantities in a single loop
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       do j = jlo, jhi
          do i = ilo, ihi

             compn = cdtdy*(fy(i,j+1,n)-fy(i,j,n))

             if (i >= ilo+1) then
                rr = qp(i,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
                compo = rr*qp(i,j,nq) - compn
                qpo(i,j,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
             end if

             if (i <= ihi-1) then
                rr = qm(i+1,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
                compo = rr*qm(i+1,j,nq) - compn
                qmo(i+1,j,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
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
          lambda(:) = lam(i,j,:)
          ugc = 0.5d0*(ugp+ugm)
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

#ifdef RADIATION
          lamge(:) = lambda(:) * (ergp(:)-ergm(:))
          luge(:) = ugc * lamge(:)
          dre = -cdtdy*sum(luge)

          if (fspace_type .eq. 1 .and. comoving) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
             end do
          else if (fspace_type .eq. 2) then
             do g=0, ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
             end do
          else ! mixed frame
             der(:) = cdtdy * luge
          end if
#endif

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to x-states
          ! for the fluid variables
          !-------------------------------------------------------------------

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
             rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                          -cdtdy*(pggp-pggm)
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
             qpo(i,j,QRHO  ) = rhotmp           + hdt*srcQ(i,j,QRHO)
             qpo(i,j,QU    ) = runewr/rhotmp
             qpo(i,j,QV    ) = rvnewr/rhotmp

             ! if ppm_trace_sources == 1, then we already added the
             ! piecewise parabolic traced source terms to the normal edge
             ! states
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qpo(i,j,QU:QV) = qpo(i,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinr = HALF*(runewr**2+rvnewr**2+(rhotmp*qpo(i,j,QW))**2)/rhotmp
             qpo(i,j,QREINT) = renewr - rhoekinr + hdt*srcQ(i,j,QREINT)

             if (.not. reset_state) then
                if (transverse_reset_rhoe == 1 .and. qpo(i,j,QREINT) <= ZERO) then
                   qpo(i,j,QREINT) = qp(i,j,QREINT) - &
                        cdtdy*(fy(i,j+1,UEINT)- fy(i,j,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,QREINT) < ZERO .and. qpo(i,j,QRHO) > ZERO) then
                      eos_state % rho = qpo(i,j,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,QFS:QFS-1+nspec)

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

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qpo(i,j,QPRES ) = pnewr
                      qpo(i,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      pnewr = qp(i  ,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-ONE))
                      qpo(i,j,QPRES) = pnewr + hdt*srcQ(i,j,QPRES)
                   endif

                   qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qpo(i,j,QGAME) = qp(i,j,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav-gamc(i,j))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qpo(i,j,QPRES) = qpo(i,j,QREINT)*(qpo(i,j,QGAME)-ONE)

                endif
             else
                qpo(i,j,QPRES) = qp(i,j,QPRES) + hdt*srcQ(i,j,QPRES)
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
             rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                           -cdtdy*(pggp-pggm)
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
             qmo(i+1,j,QRHO  ) = rhotmp            + hdt*srcQ(i,j,QRHO)
             qmo(i+1,j,QU    ) = runewl/rhotmp
             qmo(i+1,j,QV    ) = rvnewl/rhotmp

             ! if ppm_trace_sources == 1, then we already added the
             ! piecewise parabolic traced source terms to the normal edge
             ! states
             if (ppm_trace_sources == 0 .or. ppm_type == 0) then
                qmo(i+1,j,QU:QV) = qmo(i+1,j,QU:QV) + hdt*srcQ(i,j,QU:QV)
             endif

             rhoekinl = HALF*(runewl**2+rvnewl**2+(rhotmp*qmo(i+1,j,QW))**2)/rhotmp
             qmo(i+1,j,QREINT) = renewl - rhoekinl + hdt*srcQ(i,j,QREINT)

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

                      call eos(eos_input_re, eos_state)

                      pnewr = eos_state % p
                      qmo(i+1,j,QPRES ) = pnewr
                      qmo(i+1,j,QREINT) = eos_state % e * eos_state % rho
                   else
                      pnewl = qm(i+1,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-ONE))
                      qmo(i+1,j,QPRES) = pnewl + hdt*srcQ(i,j,QPRES)
                   endif

                   qmo(i+1,j,QPRES) = max(qmo(i+1,j,QPRES),small_pres)

                else

                   ! Update gammae with its transverse terms
                   qmo(i+1,j,QGAME) = qm(i+1,j,QGAME) + &
                        cdtdy*( (geav-ONE)*(geav-gamc(i,j))*du - uav*dge )

                   ! and compute the p edge state from this and (rho e)
                   qmo(i+1,j,QPRES) = qmo(i+1,j,QREINT)*(qmo(i+1,j,QGAME)-ONE)

                endif
             else
                qmo(i+1,j,QPRES) = qm(i+1,j,QPRES) + hdt*srcQ(i,j,QPRES)
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
