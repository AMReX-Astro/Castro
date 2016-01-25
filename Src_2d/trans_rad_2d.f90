module transverse_rad_module

  use bl_constants_module
  use network, only : nspec
  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                 URHO, UMX, UMY, UEDEN, &
                                 GDPRES, GDU, GDV, GDERADS, ngdnv, &
                                 npassive, upass_map, qpass_map, ppm_trace_sources, &
                                 ppm_type, small_pres, transverse_reset_density
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
       fspace_type, comoving
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor

  implicit none

  private

  public :: transx_rad, transy_rad

contains
  
  subroutine transx_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &
                        qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                        fx,  fx_l1,  fx_l2,  fx_h1,  fx_h2, &
                        rfx, rfx_l1, rfx_l2, rfx_h1, rfx_h2, &
                        qgdx, qgdx_l1, qgdx_l2, qgdx_h1, qgdx_h2, &
                        gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                        srcQ, src_l1, src_l2, src_h1, src_h2, &
                        hdt, cdtdx,  &
                        area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                        vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                        ilo, ihi, jlo, jhi)

    integer lam_l1,lam_l2,lam_h1,lam_h2
    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer gc_l1, gc_l2, gc_h1, gc_h2
    integer fx_l1, fx_l2, fx_h1, fx_h2
    integer rfx_l1, rfx_l2, rfx_h1, rfx_h2
    integer qgdx_l1, qgdx_l2, qgdx_h1, qgdx_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer area1_l1, area1_l2, area1_h1, area1_h2
    integer vol_l1, vol_l2, vol_h1, vol_h2
    integer ilo, ihi, jlo, jhi

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision qm (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qp (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision  fx( fx_l1: fx_h1, fx_l2: fx_h2,NVAR)
    double precision rfx(rfx_l1:rfx_h1,rfx_l2:rfx_h2,0:ngroups-1)
    double precision  qgdx(qgdx_l1:qgdx_h1, qgdx_l2:qgdx_h2,ngdnv)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
    double precision hdt, cdtdx

    integer :: i, j, g
    integer :: n, nq, ipassive

    double precision :: rr, rrnew, compo, compn
    double precision :: rrr, rur, rvr, rer, ekinr, rhoekinr
    double precision :: rrnewr, runewr, rvnewr, renewr
    double precision :: rrl, rul, rvl, rel, ekinl, rhoekinl
    double precision :: rrnewl, runewl, rvnewl, renewl
    double precision :: ugp, ugm, dAup, pav, dAu, pnewl,pnewr
    double precision :: rhotmp

    double precision pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    double precision eddf, f1, ugc, divu


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

                compo = rr*qp(i,j  ,nq) - compn
                qpo(i,j  ,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
             endif

             if (j <= jhi-1) then
                rr = qm(i,j+1,QRHO)
                rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                     area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j)

                compo = rr*qm(i,j+1,nq) - compn
                qmo(i,j+1,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
             endif

          enddo
       enddo
    enddo

    ! hydro variables
    do j = jlo, jhi
       do i = ilo, ihi

          lambda(:) = lam(i,j,:)
          pggp = qgdx(i+1,j,GDPRES)
          pggm = qgdx(i  ,j,GDPRES)
          ugp = qgdx(i+1,j,GDU)
          ugm = qgdx(i  ,j,GDU)
          ugc = 0.5d0*(ugp+ugm)
          ergp(:) = qgdx(i+1,j,GDERADS:GDERADS-1+ngroups)
          ergm(:) = qgdx(i  ,j,GDERADS:GDERADS-1+ngroups)

          ! we need to augment our conserved system with either a p
          ! equation to be able to deal with the general EOS
          dAup = area1(i+1,j)*pggp*ugp - area1(i,j)*pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          dAu = area1(i+1,j)*ugp-area1(i,j)*ugm

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to y-states
          ! for the fluid variables
          !-------------------------------------------------------------------

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
             err(:) = qp(i,j,qrad:qradhi)


             ! Add transverse predictor
             rrnewr = rrr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                  area1(i,j)*fx(i,j,URHO))/vol(i,j)
             runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                  area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                  -HALF*hdt*(area1(i+1,j)+area1(i,j))*((pggp-pggm)+sum(lamge))/vol(i,j)
             rvnewr = rvr - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                  area1(i,j)*fx(i,j,UMY))/vol(i,j)
             renewr = rer - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                  area1(i,j)*fx(i,j,UEDEN))/vol(i,j) &
                  + dre

             ernewr(:) = err(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
                  area1(i,j)*rfx(i,j,:))/vol(i,j) &
                  + der(:)

             if (transverse_reset_density == 1) then
                if (rrnewr < ZERO) then
                   rrnewr = rrr
                   runewr = rur
                   rvnewr = rvr
                   renewr = rer
                   ernewr(:) = err(:)
                endif
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
             qpo(i,j,QREINT)= renewr - rhoekinr + hdt*srcQ(i,j,QREINT)

             ! we are expressing the pressure evolution as:
             !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
             ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
             ! but these are divergences, so we need area factors
             pnewr = qp(i,j  ,QPRES) - hdt*(dAup + pav*dAu*(gamc(i,j)-ONE))/vol(i,j)
             qpo(i,j,QPRES) =  pnewr            + hdt*srcQ(i,j,QPRES)
             qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

             qpo(i,j,qrad:qradhi) = ernewr(:)
             qpo(i,j,qptot)   = sum(lambda*ernewr) + qpo(i,j,QPRES)
             qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)

          endif


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
             erl(:) = qm(i,j+1,qrad:qradhi)

             rrnewl = rrl - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                  area1(i,j)*fx(i,j,URHO))/vol(i,j)
             runewl = rul - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                  area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                  -HALF*hdt*(area1(i+1,j)+area1(i,j))*((pggp-pggm)+sum(lamge))/vol(i,j)
             rvnewl = rvl - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                  area1(i,j)*fx(i,j,UMY))/vol(i,j)
             renewl = rel - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                  area1(i,j)*fx(i,j,UEDEN))/vol(i,j) &
                  + dre
             ernewl(:) = erl(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
                  area1(i,j)*rfx(i,j,:))/vol(i,j) &
                  + der(:)

             if (transverse_reset_density == 1) then
                if (rrnewl < ZERO) then
                   rrnewl = rrl
                   runewl = rul
                   rvnewl = rvl
                   renewl = rel
                   ernewl(:) = erl(:)
                endif
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

             ! we are expressing the pressure evolution as:
             !   p_t + div{Up} + (gamma_1 - 1)p div{U} = 0
             ! The transverse term is d(up)/dx + (gamma_1 - 1)p du/dx,
             ! but these are divergences, so we need area factors
             pnewl = qm(i,j+1,QPRES) - hdt*(dAup + pav*dAu*(gamc(i,j)-ONE))/vol(i,j)
             qmo(i,j+1,QPRES) = pnewl +hdt*srcQ(i,j,QPRES)
             qmo(i,j+1,QPRES) = max(qmo(i,j+1,QPRES),small_pres)

             qmo(i,j+1,qrad:qradhi) = ernewl(:)
             qmo(i,j+1,qptot)   = sum(lambda*ernewl) + qmo(i,j+1,QPRES)
             qmo(i,j+1,qreitot) = sum(qmo(i,j+1,qrad:qradhi)) + qmo(i,j+1,QREINT)
             
          endif

       enddo
    enddo

  end subroutine transx_rad


  subroutine transy_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &
                        qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                        fy, fy_l1, fy_l2, fy_h1, fy_h2, &
                        rfy,rfy_l1,rfy_l2,rfy_h1,rfy_h2, &
                        qgdy, qgdy_l1, qgdy_l2, qgdy_h1, qgdy_h2, &
                        gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                        srcQ, src_l1, src_l2, src_h1, src_h2, &
                        hdt, cdtdy, ilo, ihi, jlo, jhi)

    integer lam_l1,lam_l2,lam_h1,lam_h2
    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer gc_l1, gc_l2, gc_h1, gc_h2
    integer fy_l1, fy_l2, fy_h1, fy_h2
    integer rfy_l1, rfy_l2, rfy_h1, rfy_h2
    integer qgdy_l1, qgdy_l2, qgdy_h1, qgdy_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer ilo, ihi, jlo, jhi

    double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
    double precision qm (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qp (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
    double precision  fy( fy_l1: fy_h1, fy_l2: fy_h2,NVAR)
    double precision rfy(rfy_l1:rfy_h1,rfy_l2:rfy_h2,0:ngroups-1)
    double precision  qgdy(qgdy_l1:qgdy_h1,qgdy_l2:qgdy_h2,ngdnv)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision hdt, cdtdy

    integer :: i, j, g
    integer :: n, nq, ipassive

    double precision :: rr,rrnew
    double precision :: ugp, ugm, dup, pav, du, pnewr,pnewl
    double precision :: rrr, rur, rvr, rer, ekinr, rhoekinr
    double precision :: rrnewr, runewr, rvnewr, renewr
    double precision :: rrl, rul, rvl, rel, ekinl, rhoekinl
    double precision :: rrnewl, runewl, rvnewl, renewl
    double precision :: rhotmp
    double precision :: compo, compn

    double precision :: pggp, pggm, dre, dmom
    double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
         der, ernewr, ernewl
    double precision :: eddf, f1, ugc


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
             endif

             if (i <= ihi-1) then
                rr = qm(i+1,j,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))

                compo = rr*qm(i+1,j,nq) - compn
                qmo(i+1,j,nq) = compo/rrnew + hdt*srcQ(i,j,nq)
             endif
          enddo
       enddo
    enddo

    ! hydro variables
    do j = jlo, jhi
       do i = ilo, ihi

          lambda(:) = lam(i,j,:)
          pggp = qgdy(i,j+1,GDPRES)
          pggm = qgdy(i,j  ,GDPRES)
          ugp = qgdy(i,j+1,GDV)
          ugm = qgdy(i,j  ,GDV)
          ugc = 0.5d0*(ugp+ugm)
          ergp(:) = qgdy(i,j+1,GDERADS:GDERADS-1+ngroups)
          ergm(:) = qgdy(i,j  ,GDERADS:GDERADS-1+ngroups)

          ! we need to augment our conserved system with a p
          ! equation to be able to deal with the general EOS
          dup = pggp*ugp - pggm*ugm
          pav = 0.5d0*(pggp+pggm)
          du = ugp-ugm

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
             err(:) = qp(i,j,qrad:qradhi)

             ! Add transverse predictor
             rrnewr = rrr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))
             
             runewr = rur - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
             rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                  -cdtdy*((pggp-pggm) + sum(lamge))
             renewr = rer - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) &
                  + dre

             ernewr(:) = err(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
                  + der(:)

             if (transverse_reset_density == 1) then
                if (rrnewr <= ZERO) then
                   rrnewr = rrr
                   runewr = rur
                   rvnewr = rvr
                   renewr = rer
                   ernewr(:) = err(:)
                endif
             endif

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

             pnewr = qp(i  ,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-ONE))
             qpo(i,j,QPRES) =  pnewr            + hdt*srcQ(i,j,QPRES)
             qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

             qpo(i,j,qrad:qradhi) = ernewr(:)
             qpo(i,j,qptot  ) = sum(lambda*ernewr) + qpo(i,j,QPRES)
             qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)

          endif


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
             erl(:) = qm(i+1,j,qrad:qradhi)

             rrnewl = rrl - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))
             runewl = rul - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
             rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                  -cdtdy*((pggp-pggm) + sum(lamge))
             renewl = rel - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) &
                  + dre

             ernewl(:) = erl(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
                  + der(:)

             if (transverse_reset_density == 1) then
                if (rrnewl <= ZERO) then
                   rrnewl = rrl
                   runewl = rul
                   rvnewl = rvl
                   renewl = rel
                   ernewl(:) = erl(:)
                endif
             endif

             ! convert back to non-conservation form
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

             pnewl = qm(i+1,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-ONE))
             qmo(i+1,j,QPRES) = pnewl             + hdt*srcQ(i,j,QPRES)
             qmo(i+1,j,QPRES) = max(qmo(i+1,j,QPRES),small_pres)
             
             qmo(i+1,j,qrad:qradhi) = ernewl(:)
             qmo(i+1,j,qptot  ) = sum(lambda*ernewl) + qmo(i+1,j,QPRES)
             qmo(i+1,j,qreitot) = sum(qmo(i+1,j,qrad:qradhi)) + qmo(i+1,j,QREINT)
             
          endif

       enddo
    enddo

  end subroutine transy_rad

end module transverse_rad_module


