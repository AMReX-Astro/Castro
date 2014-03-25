module transverse_module

  implicit none

contains

  subroutine transx(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                    fx, fx_l1, fx_l2, fx_h1, fx_h2, &
                    pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                    ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                    gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                    srcQ, src_l1, src_l2, src_h1, src_h2, &
                    grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                    hdt, cdtdx,  &
                    area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                    vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                    ilo, ihi, jlo, jhi)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, &
                                   QFA, QFS, QFX, &
                                   URHO, UMX, UMY, UEDEN, UEINT, &
                                   UFA, UFS, UFX, &
                                   nadv, small_pres, small_temp, &
                                   transverse_use_eos, ppm_type, ppm_trace_grav
    use eos_module

    implicit none

    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer gc_l1, gc_l2, gc_h1, gc_h2
    integer fx_l1, fx_l2, fx_h1, fx_h2
    integer pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2
    integer ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer gv_l1, gv_l2, gv_h1, gv_h2
    integer area1_l1, area1_l2, area1_h1, area1_h2
    integer vol_l1, vol_l2, vol_h1, vol_h2
    integer ilo, ihi, jlo, jhi
    
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,NVAR)
    double precision pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
    double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
    double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
    double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
    double precision hdt, cdtdx
    
    integer i, j
    integer n, nq, iadv
    integer ispec, iaux
    
    double precision rr, rrnew, compo, compn
    double precision rrr, rur, rvr, rer, ekinr, rhoekinr
    double precision rrnewr, runewr, rvnewr, renewr
    double precision rrl, rul, rvl, rel, ekinl, rhoekinl
    double precision rrnewl, runewl, rvnewl, renewl
    double precision pgp, pgm, ugp, ugm, dup, pav, du, pnewl,pnewr
    double precision rhotmp

    type (eos_t) :: eos_state
    
    ! NOTE: it is better *not* to protect against small density in this routine

    do iadv = 1, nadv
       n  = UFA + iadv - 1
       nq = QFA + iadv - 1
       
       do j = jlo, jhi 
          do i = ilo, ihi 

             rr = qp(i,j,  QRHO)
             rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                               area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

             compo = rr*qp(i,j  ,nq)
             compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                  area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

             qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

             rr = qm(i,j+1,QRHO)
             rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                               area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

             compo = rr*qm(i,j+1,nq)
             compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                  area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 
           
             qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)
             
          enddo
       enddo
    enddo

    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       
       do j = jlo, jhi 
          do i = ilo, ihi 

             rr = qp(i  ,j,QRHO)
             rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                               area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

             compo = rr*qp(i,j  ,nq)
             compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                  area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

             qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

             rr = qm(i,j+1,QRHO)
             rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                               area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

             compo  = rr*qm(i,j+1,nq)
             compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                  area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

             qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

          enddo
       enddo
    enddo

    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rr = qp(i  ,j,QRHO)
             rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                               area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

             compo = rr*qp(i,j  ,nq)
             compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                  area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

             qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

             rr = qm(i,j+1,QRHO)
             rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                               area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

             compo = rr*qm(i,j+1,nq)
             compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                  area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

             qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          pgp = pgdx(i+1,j)
          pgm = pgdx(i,j)
          ugp = ugdx(i+1,j)
          ugm = ugdx(i,j)

          ! Convert to conservation form
          rrr = qp(i,j,QRHO)
          rur = rrr*qp(i,j,QU)
          rvr = rrr*qp(i,j,QV)
          ekinr = HALF*rrr*(qp(i,j,QU)**2 + qp(i,j,QV)**2)
          rer = qp(i,j,QREINT) + ekinr
          
          rrl = qm(i,j+1,QRHO)
          rul = rrl*qm(i,j+1,QU)
          rvl = rrl*qm(i,j+1,QV)
          ekinl = HALF*rrl*(qm(i,j+1,QU)**2 + qm(i,j+1,QV)**2)
          rel = qm(i,j+1,QREINT) + ekinl
          
          ! Add transverse predictor
          rrnewr = rrr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                              area1(i,j)*fx(i,j,URHO))/vol(i,j) 
          runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                              area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                              -0.5*hdt*(area1(i+1,j)+area1(i,j))*(pgp-pgm)/vol(i,j) 
          rvnewr = rvr - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                              area1(i,j)*fx(i,j,UMY))/vol(i,j) 
          renewr = rer - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                              area1(i,j)*fx(i,j,UEDEN))/vol(i,j) 

          rrnewl = rrl - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                              area1(i,j)*fx(i,j,URHO))/vol(i,j) 
          runewl = rul - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                              area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                              -0.5*hdt*(area1(i+1,j)+area1(i,j))*(pgp-pgm)/vol(i,j) 
          rvnewl = rvl - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                              area1(i,j)*fx(i,j,UMY))/vol(i,j) 
          renewl = rel - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                              area1(i,j)*fx(i,j,UEDEN))/vol(i,j) 
          
          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          du = ugp-ugm
        
          ! Convert back to non-conservation form
          rhotmp = rrnewr
          qpo(i,j,QRHO) = rhotmp        + hdt*srcQ(i,j,QRHO)
          qpo(i,j,QU  ) = runewr/rhotmp + hdt*srcQ(i,j,QU)  
          qpo(i,j,QV  ) = rvnewr/rhotmp + hdt*srcQ(i,j,QV)  
          rhoekinr = HALF*(runewr**2+rvnewr**2)/rhotmp
          qpo(i,j,QREINT)= renewr - rhoekinr + hdt*srcQ(i,j,QREINT)

          ! Optionally, use the EOS to calculate the pressure.

          if (transverse_use_eos .eq. 1) then
             eos_state % rho = qpo(i,j,QRHO)
             eos_state % e   = qpo(i,j,QREINT) / qpo(i,j,QRHO)
             eos_state % T   = small_temp
             eos_state % xn  = qpo(i,j,QFS:QFS+nspec-1)

             call eos(eos_input_re, eos_state)

             pnewr = eos_state % p
             qpo(i,j,QPRES ) = pnewr
             qpo(i,j,QREINT) = eos_state % e * eos_state % rho
          else           
             pnewr = qp(i,j,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j)-ONE))
             qpo(i,j,QPRES) = pnewr + hdt*srcQ(i,j,QPRES)
          endif
          
          qpo(i,j,QPRES) = max(qpo(i,j,QPRES),small_pres)

          ! Convert back to non-conservation form
          rhotmp = rrnewl
          qmo(i,j+1,QRHO) = rhotmp         + hdt*srcQ(i,j,QRHO)
          qmo(i,j+1,QU  ) = runewl/rhotmp  + hdt*srcQ(i,j,QU) 
          qmo(i,j+1,QV  ) = rvnewl/rhotmp  + hdt*srcQ(i,j,QV) 
          rhoekinl = HALF*(runewl**2+rvnewl**2)/rhotmp
          qmo(i,j+1,QREINT)= renewl - rhoekinl +hdt*srcQ(i,j,QREINT)

          ! Optionally, use the EOS to calculate the pressure.

          if (transverse_use_eos .eq. 1) then
             eos_state % rho = qmo(i,j+1,QRHO)
             eos_state % e   = qmo(i,j+1,QREINT) / qmo(i,j+1,QRHO)
             eos_state % T   = small_temp
             eos_state % xn  = qmo(i,j+1,QFS:QFS+nspec-1)

             call eos(eos_input_re, eos_state)

             pnewr = eos_state % p
             qmo(i,j+1,QPRES ) = pnewr
             qmo(i,j+1,QREINT) = eos_state % e * eos_state % rho
          else           
             pnewl = qm(i,j+1,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j)-ONE))
             qmo(i,j+1,QPRES) = pnewl + hdt*srcQ(i,j,QPRES)
          endif
          
          qmo(i,j+1,QPRES) = max(qmo(i,j+1,QPRES),small_pres)
        
          ! if ppm_trace_grav == 1, then we already added the
          ! piecewise parabolic traced gravity to the normal edge
          ! states
          if (ppm_trace_grav == 0 .or. ppm_type == 0) then
             qpo(i,j,QU  ) = qpo(i,j,QU  ) + hdt*grav(i,j,1)
             qpo(i,j,QV  ) = qpo(i,j,QV  ) + hdt*grav(i,j,2)
             
             qmo(i,j+1,QU  ) = qmo(i,j+1,QU  ) + hdt*grav(i,j,1)
             qmo(i,j+1,QV  ) = qmo(i,j+1,QV  ) + hdt*grav(i,j,2)
          endif
          
       enddo
    enddo
    
  end subroutine transx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine transy(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                    fy,fy_l1,fy_l2,fy_h1,fy_h2, &
                    pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                    ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                    gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                    srcQ, src_l1, src_l2, src_h1, src_h2, &
                    grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                    hdt, cdtdy, ilo, ihi, jlo, jhi)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, &
                                   QFA, QFS, QFX, &
                                   URHO, UMX, UMY, UEDEN, UEINT, &
                                   UFA, UFS, UFX, &
                                   nadv, small_pres, small_temp, &
                                   transverse_use_eos, ppm_type, ppm_trace_grav
    use eos_module

    implicit none

    integer qd_l1, qd_l2, qd_h1, qd_h2
    integer gc_l1, gc_l2, gc_h1, gc_h2
    integer fy_l1, fy_l2, fy_h1, fy_h2
    integer pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2
    integer ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2
    integer src_l1, src_l2, src_h1, src_h2
    integer gv_l1, gv_l2, gv_h1, gv_h2
    integer ilo, ihi, jlo, jhi

    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,NVAR)
    double precision pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
    double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
    double precision hdt, cdtdy
    
    integer i, j
    integer n, nq, iadv, ispec, iaux
  
    double precision rr,rrnew
    double precision pgp, pgm, ugp, ugm, dup, pav, du, pnewr,pnewl
    double precision rrr, rur, rvr, rer, ekinr, rhoekinr
    double precision rrnewr, runewr, rvnewr, renewr
    double precision rrl, rul, rvl, rel, ekinl, rhoekinl
    double precision rrnewl, runewl, rvnewl, renewl
    double precision rhotmp
    double precision compo, compn

    type (eos_t) :: eos_state
    
    ! NOTE: it is better *not* to protect against small density in this routine

    do iadv = 1, nadv
       n  = UFA + iadv - 1
       nq = QFA + iadv - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rr = qp(i,j,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

             compo = rr*qp(i,j,nq)
             compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

             qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)
           
             rr = qm(i+1,j,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

             compo = rr*qm(i+1,j,nq)
             compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

             qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

          enddo
       enddo
    enddo

    do ispec = 1, nspec 
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1

       do j = jlo, jhi 
          do i = ilo, ihi 

             rr = qp(i,j,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

             compo = rr*qp(i,j,nq)
             compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

             qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)
             
             rr = qm(i+1,j,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 
           
             compo = rr*qm(i+1,j,nq)
             compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 
           
             qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)
           
          enddo
       enddo
    enddo

    do iaux = 1, naux 
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       
       do j = jlo, jhi 
          do i = ilo, ihi 
             
             rr = qp(i,j,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 
             
             compo = rr*qp(i,j,nq)
             compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 
             
             qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)
             
             rr = qm(i+1,j,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 
             
             compo = rr*qm(i+1,j,nq)
             compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 
             
             qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)
             
          enddo
       enddo
    enddo
    
    do j = jlo, jhi 
       do i = ilo, ihi 

          pgp = pgdy(i,j+1)
          pgm = pgdy(i,j)
          ugp = ugdy(i,j+1)
          ugm = ugdy(i,j)
          
          ! Convert to conservation form
          rrr = qp(i,j,QRHO)
          rur = rrr*qp(i,j,QU)
          rvr = rrr*qp(i,j,QV)
          ekinr = HALF*rrr*(qp(i,j,QU)**2 + qp(i,j,QV)**2)
          rer = qp(i,j,QREINT) + ekinr
          
          rrl = qm(i+1,j,QRHO)
          rul = rrl*qm(i+1,j,QU)
          rvl = rrl*qm(i+1,j,QV)
          ekinl = HALF*rrl*(qm(i+1,j,QU)**2 + qm(i+1,j,QV)**2)
          rel = qm(i+1,j,QREINT) + ekinl
          
          ! Add transverse predictor
          rrnewr = rrr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO)) 
          
          runewr = rur - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX)) 
          rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                        -cdtdy*(pgp-pgm) 
          renewr = rer - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) 

          rrnewl = rrl - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO)) 
          runewl = rul - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX)) 
          rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                        -cdtdy*(pgp-pgm) 
          renewl = rel - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) 

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          du = ugp-ugm

          ! convert back to non-conservation form
          rhotmp =  rrnewr
          qpo(i,j,QRHO  ) = rhotmp           + hdt*srcQ(i,j,QRHO)
          qpo(i,j,QU    ) = runewr/rhotmp    + hdt*srcQ(i,j,QU) 
          qpo(i,j,QV    ) = rvnewr/rhotmp    + hdt*srcQ(i,j,QV) 
          rhoekinr = HALF*(runewr**2+rvnewr**2)/rhotmp
          qpo(i,j,QREINT) = renewr - rhoekinr + hdt*srcQ(i,j,QREINT)

          ! Optionally, use the EOS to calculate the pressure.

          if (transverse_use_eos .eq. 1) then
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

          
          rhotmp =  rrnewl
          qmo(i+1,j,QRHO  ) = rhotmp            + hdt*srcQ(i,j,QRHO)
          qmo(i+1,j,QU    ) = runewl/rhotmp     + hdt*srcQ(i,j,QU) 
          qmo(i+1,j,QV    ) = rvnewl/rhotmp     + hdt*srcQ(i,j,QV) 
          rhoekinl = HALF*(runewl**2+rvnewl**2)/rhotmp
          qmo(i+1,j,QREINT) = renewl - rhoekinl + hdt*srcQ(i,j,QREINT)

          ! Optionally, use the EOS to calculate the pressure.

          if (transverse_use_eos .eq. 1) then
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
          
          ! if ppm_trace_grav == 1, then we already added the
          ! piecewise parabolic traced gravity to the normal edge
          ! states
          if (ppm_trace_grav == 0 .or. ppm_type == 0) then
             qpo(i,j,QU    ) = qpo(i,j,QU    ) + hdt*grav(i,j,1)
             qpo(i,j,QV    ) = qpo(i,j,QV    ) + hdt*grav(i,j,2)
             
             qmo(i+1,j,QU    ) = qmo(i+1,j,QU    ) + hdt*grav(i,j,1)
             qmo(i+1,j,QV    ) = qmo(i+1,j,QV    ) + hdt*grav(i,j,2)
          endif
          
       enddo
    enddo
    
  end subroutine transy

end module transverse_module

