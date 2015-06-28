module transverse_module

  use bl_constants_module

  implicit none

contains

  !===========================================================================
  ! transx1
  !===========================================================================
  subroutine transx1(qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                     ugdnvx,pgdnvx,gegdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                     gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                     cdtdx,ilo,ihi,jlo,jhi,kc,k3d)
    
    ! Note that what we call ilo here is ilo = lo(1)
    ! Note that what we call ihi here is ihi = hi(1)
    ! Note that what we call jlo here is jlo = lo(2) - 1
    ! Note that what we call jhi here is jhi = hi(2) + 1
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module
    
    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,k3d
    
    double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision gegdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdx
    
    integer i, j
    integer n, nq

    double precision rrnew, rr
    double precision rrry, rrly
    double precision rury, ruly
    double precision rvry, rvly
    double precision rwry, rwly
    double precision ekenry, ekenly
    double precision rery, rely
    double precision rrnewry, rrnewly
    double precision runewry, runewly
    double precision rvnewry, rvnewly
    double precision rwnewry, rwnewly
    double precision renewry, renewly
    double precision pnewry, pnewly
    double precision rhoekenry, rhoekenly
    double precision compn, compu
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))
             
             if (j.ge.jlo+1) then
                rr = qyp(i,j,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qyp(i,j,kc,nq) - compn
                qypo(i,j,kc,nq) = compu/rrnew
             end if
             
             if (j.le.jhi-1) then
                rr = qym(i,j+1,kc,QRHO)
                rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
                compu = rr*qym(i,j+1,kc,nq) - compn
                qymo(i,j+1,kc,nq) = compu/rrnew
             end if
             
          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to y-states
          ! for the fluid variables
          !-------------------------------------------------------------------
          
          pgp = pgdnvx(i+1,j,kc)
          pgm = pgdnvx(i,j,kc)
          ugp = ugdnvx(i+1,j,kc)
          ugm = ugdnvx(i,j,kc)
          gegp = gegdnvx(i+1,j,kc)
          gegm = gegdnvx(i,j,kc)

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm
                    
          !----------------------------------------------------------------
          ! qypo state
          !----------------------------------------------------------------
                    
          ! Convert back to primitive form
          if (j.ge.jlo+1) then

             ! Convert to conservation form
             rrry = qyp(i,j,kc,QRHO)
             rury = rrry*qyp(i,j,kc,QU)
             rvry = rrry*qyp(i,j,kc,QV)
             rwry = rrry*qyp(i,j,kc,QW)
             ekenry = HALF*rrry* &
                  (qyp(i,j,kc,QU)**2 + qyp(i,j,kc,QV)**2 + qyp(i,j,kc,QW)**2)
             rery = qyp(i,j,kc,QREINT) + ekenry
             
             ! Add transverse predictor
             rrnewry = rrry - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewry = rury - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))          
             rvnewry = rvry - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewry = rwry - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewry = rery - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewry .lt. ZERO) then
                   rrnewry = rrry
                   runewry = rury
                   rvnewry = rvry
                   rwnewry = rwry
                   renewry = rery
                endif
             endif

             qypo(i,j,kc,QRHO) = rrnewry
             qypo(i,j,kc,QU) = runewry/qypo(i,j,kc,QRHO)
             qypo(i,j,kc,QV) = rvnewry/qypo(i,j,kc,QRHO)
             qypo(i,j,kc,QW) = rwnewry/qypo(i,j,kc,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,kc,QRHO)
             qypo(i,j,kc,QREINT) = renewry - rhoekenry

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by
                ! using the discretized expression for updating (rho e).

                if (qypo(i,j,kc,QREINT) .le. ZERO) then
                   qypo(i,j,kc,QREINT) = qyp(i,j,kc,QREINT) - &
                        cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qypo(i,j,kc,QREINT) < ZERO) then
                      eos_state % rho = qypo(i,j,kc,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qypo(i,j,kc,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qypo(i,j,kc,QREINT) = qypo(i,j,kc,QRHO)*eos_state % e
                      qypo(i,j,kc,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qypo(i,j,kc,QRHO)
                   eos_state % e   = qypo(i,j,kc,QREINT) / qypo(i,j,kc,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qypo(i,j,kc,QFS:QFS+nspec-1)
                
                   call eos(eos_input_re, eos_state)
                
                   pnewry = eos_state % p
                   qypo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-ONE))
                endif

                qypo(i,j,kc,QPRES) = max(pnewry,small_pres)

             else

                ! Update gammae with its transverse terms
                qypo(i,j,kc,QGAME) = qyp(i,j,kc,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qypo(i,j,kc,QPRES) = qypo(i,j,kc,QREINT)*(qypo(i,j,kc,QGAME)-ONE)
                qypo(i,j,kc,QPRES) = max(qypo(i,j,kc,QPRES),small_pres)
             endif

          endif

          !-------------------------------------------------------------------   
          ! qymo state
          !-------------------------------------------------------------------

          ! Convert back to primitive form
          if (j.le.jhi-1) then

             ! Convert to conservation form
             rrly = qym(i,j+1,kc,QRHO)
             ruly = rrly*qym(i,j+1,kc,QU)
             rvly = rrly*qym(i,j+1,kc,QV)
             rwly = rrly*qym(i,j+1,kc,QW)
             ekenly = HALF*rrly* &
                  (qym(i,j+1,kc,QU)**2 + qym(i,j+1,kc,QV)**2 + qym(i,j+1,kc,QW)**2)
             rely = qym(i,j+1,kc,QREINT) + ekenly
             
             ! Add transverse predictor
             rrnewly = rrly - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             runewly = ruly - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
             rvnewly = rvly - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
             rwnewly = rwly - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
             renewly = rely - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewly .lt. ZERO) then
                   rrnewly = rrly 
                   runewly = ruly 
                   rvnewly = rvly 
                   rwnewly = rwly 
                   renewly = rely 
                endif
             endif

             qymo(i,j+1,kc,QRHO) = rrnewly
             qymo(i,j+1,kc,QU) = runewly/qymo(i,j+1,kc,QRHO)
             qymo(i,j+1,kc,QV) = rvnewly/qymo(i,j+1,kc,QRHO)
             qymo(i,j+1,kc,QW) = rwnewly/qymo(i,j+1,kc,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,kc,QRHO)
             qymo(i,j+1,kc,QREINT) = renewly - rhoekenly
             
             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qymo(i,j+1,kc,QREINT) .le. ZERO) then
                   qymo(i,j+1,kc,QREINT) = qym(i,j+1,kc,QREINT) - &
                        cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qymo(i,j+1,kc,QREINT) < ZERO) then
                      eos_state % rho = qymo(i,j+1,kc,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qymo(i,j+1,kc,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qymo(i,j+1,kc,QREINT) = qymo(i,j+1,kc,QRHO)*eos_state % e
                      qymo(i,j+1,kc,QPRES) = eos_state % p
                   endif
                endif
             endif


             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qymo(i,j+1,kc,QRHO)
                   eos_state % e   = qymo(i,j+1,kc,QREINT) / qymo(i,j+1,kc,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qymo(i,j+1,kc,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)
  
                   pnewly = eos_state % p
                   qymo(i,j+1,kc,QREINT) = eos_state % e * eos_state % rho
                else
                   pnewly = qym(i,j+1,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-ONE))
                endif
                
                qymo(i,j+1,kc,QPRES) = max(pnewly,small_pres)

             else
                
                ! Update gammae with its transverse terms
                qymo(i,j+1,kc,QGAME) = qym(i,j+1,kc,QGAME) + &
                     cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qymo(i,j+1,kc,QPRES) = qymo(i,j+1,kc,QREINT)*(qymo(i,j+1,kc,QGAME)-ONE)
                qymo(i,j+1,kc,QPRES) = max(qymo(i,j+1,kc,QPRES), small_pres)

             end if
          
          endif
       enddo
    enddo
    
  end subroutine transx1


  !===========================================================================
  ! transx2
  !===========================================================================
  subroutine transx2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fx,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                     ugdnvx,pgdnvx,gegdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                     gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                     cdtdx,ilo,ihi,jlo,jhi,kc,km,k3d)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,km,k3d
    
    double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision gegdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdx
    
    integer i, j
    integer n, nq
    
    double precision rrnew, rr
    double precision rrrz, rrlz
    double precision rurz, rulz
    double precision rvrz, rvlz
    double precision rwrz, rwlz
    double precision ekenrz, ekenlz
    double precision rerz, relz
    double precision rrnewrz, rrnewlz
    double precision runewrz, runewlz
    double precision rvnewrz, rvnewlz
    double precision rwnewrz, rwnewlz
    double precision renewrz, renewlz
    double precision pnewrz, pnewlz
    double precision rhoekenrz, rhoekenlz
    double precision compn, compu
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             compn = cdtdx*(fx(i+1,j,kc,n) - fx(i,j,kc,n))
             
             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nq) - compn
             qzpo(i,j,kc,nq) = compu/rrnew
             
             compn = cdtdx*(fx(i+1,j,km,n) - fx(i,j,km,n))
             
             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nq) - compn
             qzmo(i,j,kc,nq) = compu/rrnew
             
          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 
          
          !-------------------------------------------------------------------
          ! add the transverse flux difference in the x-direction to z-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------

          pgp = pgdnvx(i+1,j,kc)
          pgm = pgdnvx(i,j,kc)
          ugp = ugdnvx(i+1,j,kc)
          ugm = ugdnvx(i,j,kc)
          gegp = gegdnvx(i+1,j,kc)
          gegm = gegdnvx(i,j,kc)

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm
          
          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
          
          ! Add transverse predictor
          rrnewrz = rrrz - cdtdx*(fx(i+1,j,kc,URHO) - fx(i,j,kc,URHO))
          runewrz = rurz - cdtdx*(fx(i+1,j,kc,UMX) - fx(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdx*(fx(i+1,j,kc,UMY) - fx(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdx*(fx(i+1,j,kc,UMZ) - fx(i,j,kc,UMZ))
          renewrz = rerz - cdtdx*(fx(i+1,j,kc,UEDEN) - fx(i,j,kc,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1 .and. rrnewrz .lt. ZERO) then
             rrnewrz = rrrz 
             runewrz = rurz 
             rvnewrz = rvrz 
             rwnewrz = rwrz 
             renewrz = rerz 
          endif          
                 
          ! Convert back to primitive form
          qzpo(i,j,kc,QRHO) = rrnewrz
          qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)

          ! note: we run the risk of (rho e) being negative here
          rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).

             if (qzpo(i,j,kc,QREINT) .le. ZERO) then
                qzpo(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) - &
                     cdtdx*(fx(i+1,j,kc,UEINT) - fx(i,j,kc,UEINT) + pav*du)
                
                ! if we are still negative, then we need to reset
                if (qzpo(i,j,kc,QREINT) < ZERO) then
                   eos_state % rho = qzpo(i,j,kc,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qzpo(i,j,kc,QFS:QFS-1+nspec)
                   
                   call eos(eos_input_rt, eos_state)
                
                   qzpo(i,j,kc,QREINT) = qzpo(i,j,kc,QRHO)*eos_state % e
                   qzpo(i,j,kc,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then

             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qzpo(i,j,kc,QRHO)
                eos_state % e   = qzpo(i,j,kc,QREINT) / qzpo(i,j,kc,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qzpo(i,j,kc,QFS:QFS+nspec-1)

                call eos(eos_input_re, eos_state)

                pnewrz = eos_state % p
                qzpo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             else
                ! add the transverse term to the p evolution eq here
                pnewrz = qzp(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d)-ONE))
             endif

             qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)
             
          else

             ! Update gammae with its transverse terms
             qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME) + &
                  cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

             ! and compute the p edge state from this and (rho e)
             qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT)*(qzpo(i,j,kc,QGAME)-ONE)
             qzpo(i,j,kc,QPRES) = max(qzpo(i,j,kc,QPRES), small_pres)
             
          endif


          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------
          
          pgp = pgdnvx(i+1,j,km)
          pgm = pgdnvx(i,j,km)
          ugp = ugdnvx(i+1,j,km)
          ugm = ugdnvx(i,j,km)
          gegp = gegdnvx(i+1,j,km)
          gegm = gegdnvx(i,j,km)
          
          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm         
          dge = gegp-gegm

          ! Convert to conservation form
          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
          
          ! Add transverse predictor
          rrnewlz = rrlz - cdtdx*(fx(i+1,j,km,URHO) - fx(i,j,km,URHO))
          runewlz = rulz - cdtdx*(fx(i+1,j,km,UMX) - fx(i,j,km,UMX))
          rvnewlz = rvlz - cdtdx*(fx(i+1,j,km,UMY) - fx(i,j,km,UMY))
          rwnewlz = rwlz - cdtdx*(fx(i+1,j,km,UMZ) - fx(i,j,km,UMZ))
          renewlz = relz - cdtdx*(fx(i+1,j,km,UEDEN) - fx(i,j,km,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1 .and. rrnewlz .lt. ZERO) then
             rrnewlz = rrlz
             runewlz = rulz
             rvnewlz = rvlz
             rwnewlz = rwlz
             renewlz = relz
          endif

          ! Convert back to primitive form
          qzmo(i,j,kc,QRHO) = rrnewlz
          qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)

          ! note: we run the risk of (rho e) being negative here
          rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).

             if (qzmo(i,j,kc,QREINT) .le. ZERO) then
                qzmo(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) - &
                     cdtdx*(fx(i+1,j,km,UEINT) - fx(i,j,km,UEINT) + pav*du)
                
                ! if we are still negative, then we need to reset
                if (qzmo(i,j,kc,QREINT) < ZERO) then
                   eos_state % rho = qzmo(i,j,kc,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qzmo(i,j,kc,QFS:QFS-1+nspec)
                
                   call eos(eos_input_rt, eos_state)
                   
                   qzmo(i,j,kc,QREINT) = qzmo(i,j,kc,QRHO)*eos_state % e
                   qzmo(i,j,kc,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then

             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qzmo(i,j,kc,QRHO)
                eos_state % e   = qzmo(i,j,kc,QREINT) / qzmo(i,j,kc,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qzmo(i,j,kc,QFS:QFS+nspec-1)
                
                call eos(eos_input_re, eos_state)

                pnewlz = eos_state % p
                qzmo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             else
                pnewlz = qzm(i,j,kc,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j,k3d-1)-ONE))
             endif
             
             qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)
          
          else

             ! Update gammae with its transverse terms
             qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME) + &
                  cdtdx*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

             ! and compute the p edge state from this and (rho e)
             qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT)*(qzmo(i,j,kc,QGAME)-ONE)
             qzmo(i,j,kc,QPRES) = max(qzmo(i,j,kc,QPRES), small_pres)

          endif

       enddo
    enddo
    
  end subroutine transx2


  !===========================================================================
  ! transy1
  !===========================================================================
  subroutine transy1(qxm,qxmo,qxp,qxpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                     ugdnvy,pgdnvy,gegdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                     gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                     cdtdy,ilo,ihi,jlo,jhi,kc,k3d)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module

    implicit none
      
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,k3d
    
    double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gegdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdy
    
    integer i, j
    integer n, nq
    
    double precision rrnew, rr
    double precision compn, compu
    double precision rrrx, rrlx
    double precision rurx, rulx
    double precision rvrx, rvlx
    double precision rwrx, rwlx
    double precision ekenrx, ekenlx
    double precision rerx, relx
    double precision rrnewrx, rrnewlx
    double precision runewrx, runewlx
    double precision rvnewrx, rvnewlx
    double precision rwnewrx, rwnewlx
    double precision renewrx, renewlx
    double precision pnewrx, pnewlx
    double precision rhoekenrx, rhoekenlx
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))             

             if (i.ge.ilo+1) then
                rr = qxp(i,j,kc,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                compu = rr*qxp(i,j,kc,nq) - compn
                qxpo(i,j,kc,nq) = compu/rrnew
             end if

             if (i.le.ihi-1) then
                rr = qxm(i+1,j,kc,QRHO)
                rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
                compu = rr*qxm(i+1,j,kc,nq) - compn
                qxmo(i+1,j,kc,nq) = compu/rrnew
             end if
             
          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to x-states
          ! for the fluid variables
          !-------------------------------------------------------------------
          
          pgp = pgdnvy(i,j+1,kc)
          pgm = pgdnvy(i,j,kc)
          ugp = ugdnvy(i,j+1,kc)
          ugm = ugdnvy(i,j,kc)
          gegp = gegdnvy(i,j+1,kc)
          gegm = gegdnvy(i,j,kc)

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm
          
          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------    

          ! Convert back to primitive form
          if (i.ge.ilo+1) then
             ! Convert to conservation form
             rrrx = qxp(i,j,kc,QRHO)
             rurx = rrrx*qxp(i,j,kc,QU)
             rvrx = rrrx*qxp(i,j,kc,QV)
             rwrx = rrrx*qxp(i,j,kc,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,kc,QU)**2 + qxp(i,j,kc,QV)**2 &
                  + qxp(i,j,kc,QW)**2)
             rerx = qxp(i,j,kc,QREINT) + ekenrx
             
             ! Add transverse predictor
             rrnewrx = rrrx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             runewrx = rurx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
             rvnewrx = rvrx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
             rwnewrx = rwrx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
             renewrx = rerx - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewrx .lt. ZERO) then
                   rrnewrx = rrrx 
                   runewrx = rurx 
                   rvnewrx = rvrx 
                   rwnewrx = rwrx 
                   renewrx = rerx 
                endif
             endif
             
             qxpo(i,j,kc,QRHO) = rrnewrx
             qxpo(i,j,kc,QU) = runewrx/qxpo(i,j,kc,QRHO)
             qxpo(i,j,kc,QV) = rvnewrx/qxpo(i,j,kc,QRHO)
             qxpo(i,j,kc,QW) = rwnewrx/qxpo(i,j,kc,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,kc,QRHO)
             qxpo(i,j,kc,QREINT) = renewrx - rhoekenrx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by
                ! using the discretized expression for updating (rho e).
             
                if (qxpo(i,j,kc,QREINT) .le. ZERO) then
                   qxpo(i,j,kc,QREINT) = qxp(i,j,kc,QREINT) - &
                        cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qxpo(i,j,kc,QREINT) < ZERO) then
                      eos_state % rho = qxpo(i,j,kc,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxpo(i,j,kc,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qxpo(i,j,kc,QREINT) = qxpo(i,j,kc,QRHO) * eos_state % e
                      qxpo(i,j,kc,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.             

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxpo(i,j,kc,QRHO)
                   eos_state % e   = qxpo(i,j,kc,QREINT) / qxpo(i,j,kc,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxpo(i,j,kc,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewrx = eos_state % p
                   qxpo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - ONE))
                endif

                qxpo(i,j,kc,QPRES) = max(pnewrx,small_pres)

             else

                ! Update gammae with its transverse terms
                qxpo(i,j,kc,QGAME) = qxp(i,j,kc,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qxpo(i,j,kc,QPRES) = qxpo(i,j,kc,QREINT)*(qxpo(i,j,kc,QGAME)-ONE)
                qxpo(i,j,kc,QPRES) = max(qxpo(i,j,kc,QPRES), small_pres)
                
             endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------
          
          if (i.le.ihi-1) then
             ! Convert to conservation form
             rrlx = qxm(i+1,j,kc,QRHO)
             rulx = rrlx*qxm(i+1,j,kc,QU)
             rvlx = rrlx*qxm(i+1,j,kc,QV)
             rwlx = rrlx*qxm(i+1,j,kc,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,kc,QU)**2 + qxm(i+1,j,kc,QV)**2 &
                  + qxm(i+1,j,kc,QW)**2)
             relx = qxm(i+1,j,kc,QREINT) + ekenlx
             
             ! Add transverse predictor
             rrnewlx = rrlx - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             runewlx = rulx - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
             rvnewlx = rvlx - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
             rwnewlx = rwlx - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
             renewlx = relx - cdtdy*(fy(i,j+1,kc,UEDEN)- fy(i,j,kc,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewlx .lt. ZERO) then
                   rrnewlx = rrlx 
                   runewlx = rulx 
                   rvnewlx = rvlx 
                   rwnewlx = rwlx 
                   renewlx = relx 
                endif
             endif

             qxmo(i+1,j,kc,QRHO) = rrnewlx
             qxmo(i+1,j,kc,QU) = runewlx/qxmo(i+1,j,kc,QRHO)
             qxmo(i+1,j,kc,QV) = rvnewlx/qxmo(i+1,j,kc,QRHO)
             qxmo(i+1,j,kc,QW) = rwnewlx/qxmo(i+1,j,kc,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,kc,QRHO)
             qxmo(i+1,j,kc,QREINT) = renewlx - rhoekenlx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qxmo(i+1,j,kc,QREINT) .le. ZERO) then
                   qxmo(i+1,j,kc,QREINT) = qxm(i+1,j,kc,QREINT) - &
                        cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qxmo(i+1,j,kc,QREINT) < ZERO) then
                      eos_state % rho = qxmo(i+1,j,kc,QRHO) 
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxmo(i+1,j,kc,QFS:QFS-1+nspec) 
                   
                      call eos(eos_input_rt, eos_state)
                      
                      qxmo(i+1,j,kc,QREINT) = qxmo(i+1,j,kc,QRHO)*eos_state % e
                      qxmo(i+1,j,kc,QPRES) = eos_state % p
                   endif
                endif
             endif


             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.             

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxmo(i+1,j,kc,QRHO)
                   eos_state % e   = qxmo(i+1,j,kc,QREINT) / qxmo(i+1,j,kc,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxmo(i+1,j,kc,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)

                   pnewlx = eos_state % p
                   qxmo(i+1,j,kc,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - ONE))
                endif

                qxmo(i+1,j,kc,QPRES) = max(pnewlx,small_pres)

             else

                ! Update gammae with its transverse terms
                qxmo(i+1,j,kc,QGAME) = qxm(i+1,j,kc,QGAME) + &
                     cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qxmo(i+1,j,kc,QPRES) = qxmo(i+1,j,kc,QREINT)*(qxmo(i+1,j,kc,QGAME)-ONE)
                qxmo(i+1,j,kc,QPRES) = max(qxmo(i+1,j,kc,QPRES), small_pres)

             endif

          endif
          
       enddo
    enddo

  end subroutine transy1


  !===========================================================================
  ! transy2
  !===========================================================================
  subroutine transy2(qzm,qzmo,qzp,qzpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fy,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                     ugdnvy,pgdnvy,gegdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                     gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                     cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module

    implicit none
    
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,kc,km,k3d
    
    double precision  qzm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qzp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qzmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qzpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gegdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdy
    
    integer i, j
    integer n, nq
    
    double precision rrnew, rr
    double precision compn, compu
    double precision rrrz, rrlz
    double precision rurz, rulz
    double precision rvrz, rvlz
    double precision rwrz, rwlz
    double precision ekenrz, ekenlz
    double precision rerz, relz
    double precision rrnewrz, rrnewlz
    double precision runewrz, runewlz
    double precision rvnewrz, rvnewlz
    double precision rwnewrz, rwnewlz
    double precision renewrz, renewlz
    double precision pnewrz, pnewlz
    double precision rhoekenrz, rhoekenlz
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             compn = cdtdy*(fy(i,j+1,kc,n) - fy(i,j,kc,n))
             
             rr = qzp(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
             compu = rr*qzp(i,j,kc,nq) - compn
             qzpo(i,j,kc,nq) = compu/rrnew
             
             compn = cdtdy*(fy(i,j+1,km,n) - fy(i,j,km,n))
             
             rr = qzm(i,j,kc,QRHO)
             rrnew = rr - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
             compu = rr*qzm(i,j,kc,nq) - compn
             qzmo(i,j,kc,nq) = compu/rrnew
             
          enddo
       enddo
    enddo

    do j = jlo, jhi
       do i = ilo, ihi

          !-------------------------------------------------------------------
          ! add the transverse flux difference in the y-direction to z-states
          ! for the fluid variables
          !-------------------------------------------------------------------

          !-------------------------------------------------------------------
          ! qzpo states
          !-------------------------------------------------------------------          
          
          pgp = pgdnvy(i,j+1,kc)
          pgm = pgdnvy(i,j,kc)
          ugp = ugdnvy(i,j+1,kc)
          ugm = ugdnvy(i,j,kc)
          gegp = gegdnvy(i,j+1,kc)
          gegm = gegdnvy(i,j,kc)

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm
                    
          ! Convert to conservation form
          rrrz = qzp(i,j,kc,QRHO)
          rurz = rrrz*qzp(i,j,kc,QU)
          rvrz = rrrz*qzp(i,j,kc,QV)
          rwrz = rrrz*qzp(i,j,kc,QW)
          ekenrz = HALF*rrrz*(qzp(i,j,kc,QU)**2 + qzp(i,j,kc,QV)**2 &
               + qzp(i,j,kc,QW)**2)
          rerz = qzp(i,j,kc,QREINT) + ekenrz
          
          ! Add transverse predictor
          rrnewrz = rrrz - cdtdy*(fy(i,j+1,kc,URHO) - fy(i,j,kc,URHO))
          runewrz = rurz - cdtdy*(fy(i,j+1,kc,UMX) - fy(i,j,kc,UMX))
          rvnewrz = rvrz - cdtdy*(fy(i,j+1,kc,UMY) - fy(i,j,kc,UMY))
          rwnewrz = rwrz - cdtdy*(fy(i,j+1,kc,UMZ) - fy(i,j,kc,UMZ))
          renewrz = rerz - cdtdy*(fy(i,j+1,kc,UEDEN) - fy(i,j,kc,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1 .and. rrnewrz .lt. ZERO) then
             rrnewrz = rrrz 
             runewrz = rurz 
             rvnewrz = rvrz 
             rwnewrz = rwrz 
             renewrz = rerz 
          endif

          ! Convert back to primitive form
          qzpo(i,j,kc,QRHO) = rrnewrz
          qzpo(i,j,kc,QU) = runewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QV) = rvnewrz/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QW) = rwnewrz/qzpo(i,j,kc,QRHO)

          ! note: we run the risk of (rho e) being negative here
          rhoekenrz = HALF*(runewrz**2 + rvnewrz**2 + rwnewrz**2)/qzpo(i,j,kc,QRHO)
          qzpo(i,j,kc,QREINT) = renewrz - rhoekenrz

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).
          
             if (qzpo(i,j,kc,QREINT) .le. ZERO) then
                qzpo(i,j,kc,QREINT) = qzp(i,j,kc,QREINT) - &
                     cdtdy*(fy(i,j+1,kc,UEINT) - fy(i,j,kc,UEINT) + pav*du)
                
                ! if we are still negative, then we need to reset
                if (qzpo(i,j,kc,QREINT) < ZERO) then
                   eos_state % rho = qzpo(i,j,kc,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qzpo(i,j,kc,QFS:QFS-1+nspec)
                   
                   call eos(eos_input_rt, eos_state)
                   
                   qzpo(i,j,kc,QREINT) = qzpo(i,j,kc,QRHO)*eos_state % e
                   qzpo(i,j,kc,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then

             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qzpo(i,j,kc,QRHO)
                eos_state % e   = qzpo(i,j,kc,QREINT) / qzpo(i,j,kc,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qzpo(i,j,kc,QFS:QFS+nspec-1)

                call eos(eos_input_re, eos_state)
                
                pnewrz = eos_state % p
                qzpo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             else
                ! add the transverse term to the p evolution eq here
                pnewrz = qzp(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d) - ONE))
             endif

             qzpo(i,j,kc,QPRES) = max(pnewrz,small_pres)
          
          else

             ! Update gammae with its transverse terms
             qzpo(i,j,kc,QGAME) = qzp(i,j,kc,QGAME) + &
                  cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d))*du - uav*dge )

             ! and compute the p edge state from this and (rho e)
             qzpo(i,j,kc,QPRES) = qzpo(i,j,kc,QREINT)*(qzpo(i,j,kc,QGAME)-ONE)
             qzpo(i,j,kc,QPRES) = max(qzpo(i,j,kc,QPRES), small_pres)

          endif


          !-------------------------------------------------------------------
          ! qzmo states
          !-------------------------------------------------------------------

          pgp  =  pgdnvy(i,j+1,km)
          pgm  =  pgdnvy(i,j,km)
          ugp  =  ugdnvy(i,j+1,km)
          ugm  =  ugdnvy(i,j,km)
          gegp = gegdnvy(i,j+1,km)
          gegm = gegdnvy(i,j,km)

          ! we need to augment our conserved system with either a p
          ! equation or gammae (if we have ppm_predict_gammae = 1) to
          ! be able to deal with the general EOS
          
          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm

          ! Convert to conservation form          
          rrlz = qzm(i,j,kc,QRHO)
          rulz = rrlz*qzm(i,j,kc,QU)
          rvlz = rrlz*qzm(i,j,kc,QV)
          rwlz = rrlz*qzm(i,j,kc,QW)
          ekenlz = HALF*rrlz*(qzm(i,j,kc,QU)**2 + qzm(i,j,kc,QV)**2 &
               + qzm(i,j,kc,QW)**2)
          relz = qzm(i,j,kc,QREINT) + ekenlz
          
          ! Add transverse predictor
          rrnewlz = rrlz - cdtdy*(fy(i,j+1,km,URHO) - fy(i,j,km,URHO))
          runewlz = rulz - cdtdy*(fy(i,j+1,km,UMX) - fy(i,j,km,UMX))
          rvnewlz = rvlz - cdtdy*(fy(i,j+1,km,UMY) - fy(i,j,km,UMY))
          rwnewlz = rwlz - cdtdy*(fy(i,j+1,km,UMZ) - fy(i,j,km,UMZ))
          renewlz = relz - cdtdy*(fy(i,j+1,km,UEDEN)- fy(i,j,km,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1 .and. rrnewlz .lt. ZERO) then
             rrnewlz = rrlz
             runewlz = rulz
             rvnewlz = rvlz
             rwnewlz = rwlz
             renewlz = relz
          endif

          ! Convert back to primitive form
          qzmo(i,j,kc,QRHO) = rrnewlz
          qzmo(i,j,kc,QU) = runewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QV) = rvnewlz/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QW) = rwnewlz/qzmo(i,j,kc,QRHO)

          ! note: we run the risk of (rho e) being negative here
          rhoekenlz = HALF*(runewlz**2 + rvnewlz**2 + rwnewlz**2)/qzmo(i,j,kc,QRHO)
          qzmo(i,j,kc,QREINT) = renewlz - rhoekenlz

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).
             
             if (qzmo(i,j,kc,QREINT) .le. ZERO) then
                qzmo(i,j,kc,QREINT) = qzm(i,j,kc,QREINT) - &
                     cdtdy*(fy(i,j+1,km,UEINT) - fy(i,j,km,UEINT) + pav*du)
                
                ! if we are still negative, then we need to reset
                if (qzmo(i,j,kc,QREINT) < ZERO) then
                   eos_state % rho = qzmo(i,j,kc,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qzmo(i,j,kc,QFS:QFS-1+nspec)

                   call eos(eos_input_rt, eos_state)
                   
                   qzmo(i,j,kc,QREINT) = qzmo(i,j,kc,QRHO)*eos_state % e
                   qzmo(i,j,kc,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then

             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qzmo(i,j,kc,QRHO)
                eos_state % e   = qzmo(i,j,kc,QREINT) / qzmo(i,j,kc,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qzmo(i,j,kc,QFS:QFS+nspec-1)
                
                call eos(eos_input_re, eos_state)
                
                pnewlz = eos_state % p
                qzmo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             else
                ! add the transverse term to the p evolution eq here
                pnewlz = qzm(i,j,kc,QPRES) - cdtdy*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
             endif
             
             qzmo(i,j,kc,QPRES) = max(pnewlz,small_pres)
          
          else

             ! Update gammae with its transverse terms
             qzmo(i,j,kc,QGAME) = qzm(i,j,kc,QGAME) + &
                  cdtdy*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

             ! and compute the p edge state from this and (rho e)
             qzmo(i,j,kc,QPRES) = qzmo(i,j,kc,QREINT)*(qzmo(i,j,kc,QGAME)-ONE)
             qzmo(i,j,kc,QPRES) = max(qzmo(i,j,kc,QPRES), small_pres)

          endif

       enddo
    enddo
    
  end subroutine transy2


  !===========================================================================
  ! transz
  !===========================================================================
  subroutine transz(qxm,qxmo,qxp,qxpo, &
                    qym,qymo,qyp,qypo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                    fz,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                    ugdnvz,pgdnvz,gegdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                    gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                    cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)

    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe
    use eos_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d
    
    double precision  qxm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qxp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qym(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qyp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qxmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qxpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qymo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qypo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
    double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision gegdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision cdtdz
    
    integer n, nq
    integer i, j
    
    double precision rrnew, rr
    double precision compn, compu
    double precision rrrx, rrry, rrlx, rrly
    double precision rurx, rury, rulx, ruly
    double precision rvrx, rvry, rvlx, rvly
    double precision rwrx, rwry, rwlx, rwly
    double precision ekenrx, ekenry, ekenlx, ekenly
    double precision rerx, rery, relx, rely
    double precision rrnewrx, rrnewry, rrnewlx, rrnewly
    double precision runewrx, runewry, runewlx, runewly
    double precision rvnewrx, rvnewry, rvnewlx, rvnewly
    double precision rwnewrx, rwnewry, rwnewlx, rwnewly
    double precision renewrx, renewry, renewlx, renewly
    double precision pnewrx, pnewry, pnewlx, pnewly
    double precision rhoekenrx, rhoekenry, rhoekenlx, rhoekenly
    double precision pgp, pgm, ugp, ugm, gegp, gegm, dup, pav, du, dge, uav, geav

    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             compn = cdtdz*(fz(i,j,kc,n) - fz(i,j,km,n))
             
             if (i.ge.ilo+1) then
                rr = qxp(i,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qxp(i,j,km,nq) - compn
                qxpo(i,j,km,nq) = compu/rrnew
             end if
             
             if (j.ge.jlo+1) then
                rr = qyp(i,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qyp(i,j,km,nq) - compn
                qypo(i,j,km,nq) = compu/rrnew
             end if
             
             if (i.le.ihi-1) then
                rr = qxm(i+1,j,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qxm(i+1,j,km,nq) - compn
                qxmo(i+1,j,km,nq) = compu/rrnew
             end if
             
             if (j.le.jhi-1) then
                rr = qym(i,j+1,km,QRHO)
                rrnew = rr - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
                compu = rr*qym(i,j+1,km,nq) - compn
                qymo(i,j+1,km,nq) = compu/rrnew
             end if
             
          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          !-------------------------------------------------------------------
          ! add transverse flux difference in the z-direction to the x- and
          ! y-states for the fluid variables
          !-------------------------------------------------------------------    

          pgp  =  pgdnvz(i,j,kc)
          pgm  =  pgdnvz(i,j,km)
          ugp  =  ugdnvz(i,j,kc)
          ugm  =  ugdnvz(i,j,km)
          gegp = gegdnvz(i,j,kc)
          gegm = gegdnvz(i,j,km)

          dup = pgp*ugp - pgm*ugm
          pav = HALF*(pgp+pgm)
          uav = HALF*(ugp+ugm)
          geav = HALF*(gegp+gegm)
          du = ugp-ugm
          dge = gegp-gegm
          
          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------
          
          ! Convert back to primitive form
          if (i.ge.ilo+1) then
             ! Convert to conservation form
             rrrx = qxp(i,j,km,QRHO)
             rurx = rrrx*qxp(i,j,km,QU)
             rvrx = rrrx*qxp(i,j,km,QV)
             rwrx = rrrx*qxp(i,j,km,QW)
             ekenrx = HALF*rrrx*(qxp(i,j,km,QU)**2 + qxp(i,j,km,QV)**2 &
                  + qxp(i,j,km,QW)**2)
             rerx = qxp(i,j,km,QREINT) + ekenrx
             
             ! Add transverse predictor
             rrnewrx = rrrx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewrx = rurx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewrx = rvrx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewrx = rwrx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewrx = rerx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewrx .lt. ZERO) then
                   rrnewrx = rrrx 
                   runewrx = rurx 
                   rvnewrx = rvrx 
                   rwnewrx = rwrx 
                   renewrx = rerx 
                endif
             end if
             
             qxpo(i,j,km,QRHO) = rrnewrx
             qxpo(i,j,km,QU) = runewrx/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QV) = rvnewrx/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QW) = rwnewrx/qxpo(i,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenrx = HALF*(runewrx**2 + rvnewrx**2 + rwnewrx**2)/qxpo(i,j,km,QRHO)
             qxpo(i,j,km,QREINT) = renewrx - rhoekenrx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qxpo(i,j,km,QREINT) .le. ZERO) then
                   qxpo(i,j,km,QREINT) = qxp(i,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)

                   ! if we are still negative, then we need to reset
                   if (qxpo(i,j,km,QREINT) < ZERO) then
                      eos_state % rho = qxpo(i,j,km,QRHO) 
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxpo(i,j,km,QFS:QFS-1+nspec) 
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qxpo(i,j,km,QREINT) = qxpo(i,j,km,QRHO)*eos_state % e
                      qxpo(i,j,km,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
             
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxpo(i,j,km,QRHO)
                   eos_state % e   = qxpo(i,j,km,QREINT) / qxpo(i,j,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxpo(i,j,km,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewrx = eos_state % p
                   qxpo(i,j,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewrx = qxp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))

                endif

                qxpo(i,j,km,QPRES) = max(pnewrx,small_pres)

             else

                ! Update gammae with its transverse terms
                qxpo(i,j,km,QGAME) = qxp(i,j,km,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qxpo(i,j,km,QPRES) = qxpo(i,j,km,QREINT)*(qxpo(i,j,km,QGAME)-ONE)
                qxpo(i,j,km,QPRES) = max(qxpo(i,j,km,QPRES), small_pres)

             endif

          end if

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------
          
          if (j.ge.jlo+1) then
             ! Convert to conservation form
             rrry = qyp(i,j,km,QRHO)
             rury = rrry*qyp(i,j,km,QU)
             rvry = rrry*qyp(i,j,km,QV)
             rwry = rrry*qyp(i,j,km,QW)
             ekenry = HALF*rrry*(qyp(i,j,km,QU)**2 + qyp(i,j,km,QV)**2 &
                  + qyp(i,j,km,QW)**2)
             rery = qyp(i,j,km,QREINT) + ekenry
             
             ! Add transverse predictor
             rrnewry = rrry - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewry = rury - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewry = rvry - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewry = rwry - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewry = rery - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewry .lt. ZERO) then
                   rrnewry = rrry 
                   runewry = rury 
                   rvnewry = rvry 
                   rwnewry = rwry 
                   renewry = rery 
                endif
             end if
             
             qypo(i,j,km,QRHO) = rrnewry
             qypo(i,j,km,QU) = runewry/qypo(i,j,km,QRHO)
             qypo(i,j,km,QV) = rvnewry/qypo(i,j,km,QRHO)
             qypo(i,j,km,QW) = rwnewry/qypo(i,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenry = HALF*(runewry**2 + rvnewry**2 + rwnewry**2)/qypo(i,j,km,QRHO)
             qypo(i,j,km,QREINT) = renewry - rhoekenry

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qypo(i,j,km,QREINT) .le. ZERO) then
                   qypo(i,j,km,QREINT) = qyp(i,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qypo(i,j,km,QREINT) < ZERO) then
                      eos_state % rho = qypo(i,j,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qypo(i,j,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qypo(i,j,km,QREINT) = qypo(i,j,km,QRHO)*eos_state % e
                      qypo(i,j,km,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qypo(i,j,km,QRHO)
                   eos_state % e   = qypo(i,j,km,QREINT) / qypo(i,j,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qypo(i,j,km,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)
                   
                   pnewry = eos_state % p
                   qypo(i,j,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewry = qyp(i,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                endif

                qypo(i,j,km,QPRES) = max(pnewry,small_pres)

             else

                ! Update gammae with its transverse terms
                qypo(i,j,km,QGAME) = qyp(i,j,km,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qypo(i,j,km,QPRES) = qypo(i,j,km,QREINT)*(qypo(i,j,km,QGAME)-ONE)
                qypo(i,j,km,QPRES) = max(qypo(i,j,km,QPRES), small_pres)

             endif

          end if

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------
          
          if (i.le.ihi-1) then

             ! Convert to conservation form
             rrlx = qxm(i+1,j,km,QRHO)
             rulx = rrlx*qxm(i+1,j,km,QU)
             rvlx = rrlx*qxm(i+1,j,km,QV)
             rwlx = rrlx*qxm(i+1,j,km,QW)
             ekenlx = HALF*rrlx*(qxm(i+1,j,km,QU)**2 + qxm(i+1,j,km,QV)**2 &
                  + qxm(i+1,j,km,QW)**2)
             relx = qxm(i+1,j,km,QREINT) + ekenlx
             
             ! Add transverse predictor
             rrnewlx = rrlx - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewlx = rulx - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewlx = rvlx - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewlx = rwlx - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewlx = relx - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewlx .lt. ZERO) then
                   rrnewlx = rrlx 
                   runewlx = rulx 
                   rvnewlx = rvlx 
                   rwnewlx = rwlx 
                   renewlx = relx 
                endif
             end if
             
             qxmo(i+1,j,km,QRHO) = rrnewlx
             qxmo(i+1,j,km,QU) = runewlx/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QV) = rvnewlx/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QW) = rwnewlx/qxmo(i+1,j,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenlx = HALF*(runewlx**2 + rvnewlx**2 + rwnewlx**2)/qxmo(i+1,j,km,QRHO)
             qxmo(i+1,j,km,QREINT) = renewlx - rhoekenlx

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qxmo(i+1,j,km,QREINT) .le. ZERO) then
                   qxmo(i+1,j,km,QREINT) = qxm(i+1,j,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qxmo(i+1,j,km,QREINT) < ZERO) then
                      eos_state % rho = qxmo(i+1,j,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qxmo(i+1,j,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qxmo(i+1,j,km,QREINT) = qxmo(i+1,j,km,QRHO)*eos_state % e 
                      qxmo(i+1,j,km,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
             
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qxmo(i+1,j,km,QRHO)
                   eos_state % e   = qxmo(i+1,j,km,QREINT) / qxmo(i+1,j,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qxmo(i+1,j,km,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)

                   pnewlx = eos_state % p
                   qxmo(i+1,j,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewlx = qxm(i+1,j,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                endif
                
                qxmo(i+1,j,km,QPRES) = max(pnewlx,small_pres)

             else
                
                ! Update gammae with its transverse terms             
                qxmo(i+1,j,km,QGAME) = qxm(i+1,j,km,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )
                
                ! and compute the p edge state from this and (rho e)
                qxmo(i+1,j,km,QPRES) = qxmo(i+1,j,km,QREINT)*(qxmo(i+1,j,km,QGAME)-ONE)
                qxmo(i+1,j,km,QPRES) = max(qxmo(i+1,j,km,QPRES), small_pres)
                
             end if

          endif


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------
          
          if (j.le.jhi-1) then
             ! Convert to conservation form
             rrly = qym(i,j+1,km,QRHO)
             ruly = rrly*qym(i,j+1,km,QU)
             rvly = rrly*qym(i,j+1,km,QV)
             rwly = rrly*qym(i,j+1,km,QW)
             ekenly = HALF*rrly*(qym(i,j+1,km,QU)**2 + qym(i,j+1,km,QV)**2 &
                  + qym(i,j+1,km,QW)**2)
             rely = qym(i,j+1,km,QREINT) + ekenly
             
             ! Add transverse predictor
             rrnewly = rrly - cdtdz*(fz(i,j,kc,URHO) - fz(i,j,km,URHO))
             runewly = ruly - cdtdz*(fz(i,j,kc,UMX) - fz(i,j,km,UMX))
             rvnewly = rvly - cdtdz*(fz(i,j,kc,UMY) - fz(i,j,km,UMY))
             rwnewly = rwly - cdtdz*(fz(i,j,kc,UMZ) - fz(i,j,km,UMZ))
             renewly = rely - cdtdz*(fz(i,j,kc,UEDEN) - fz(i,j,km,UEDEN))
             
             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewly .lt. ZERO) then
                   rrnewly = rrly 
                   runewly = ruly 
                   rvnewly = rvly 
                   rwnewly = rwly 
                   renewly = rely 
                endif
             endif

             qymo(i,j+1,km,QRHO) = rrnewly
             qymo(i,j+1,km,QU) = runewly/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QV) = rvnewly/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QW) = rwnewly/qymo(i,j+1,km,QRHO)

             ! note: we run the risk of (rho e) being negative here
             rhoekenly = HALF*(runewly**2 + rvnewly**2 + rwnewly**2)/qymo(i,j+1,km,QRHO)
             qymo(i,j+1,km,QREINT) = renewly - rhoekenly

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qymo(i,j+1,km,QREINT) .le. ZERO) then
                   qymo(i,j+1,km,QREINT) = qym(i,j+1,km,QREINT) - &
                        cdtdz*(fz(i,j,kc,UEINT) - fz(i,j,km,UEINT) + pav*du)
                   
                   ! if we are still negative, then we need to reset
                   if (qymo(i,j+1,km,QREINT) < ZERO) then
                      eos_state % rho = qymo(i,j+1,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qymo(i,j+1,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qymo(i,j+1,km,QREINT) =  qymo(i,j+1,km,QRHO)*eos_state % e
                      qymo(i,j+1,km,QPRES) =  eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.             

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qymo(i,j+1,km,QRHO)
                   eos_state % e   = qymo(i,j+1,km,QREINT) / qymo(i,j+1,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qymo(i,j+1,km,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)
                
                   pnewly = eos_state % p
                   qymo(i,j+1,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewly = qym(i,j+1,km,QPRES) - cdtdz*(dup + pav*du*(gamc(i,j,k3d-1) - ONE))
                endif

                qymo(i,j+1,km,QPRES) = max(pnewly,small_pres)

             else

                ! Update gammae with its transverse terms
                qymo(i,j+1,km,QGAME) = qym(i,j+1,km,QGAME) + &
                     cdtdz*( (geav-ONE)*(geav-gamc(i,j,k3d-1))*du - uav*dge )

                ! and compute the p edge state from this and (rho e)
                qymo(i,j+1,km,QPRES) = qymo(i,j+1,km,QREINT)*(qymo(i,j+1,km,QGAME)-ONE)
                qymo(i,j+1,km,QPRES) = max(qymo(i,j+1,km,QPRES), small_pres)
                
             endif

          endif

       enddo
    enddo

  end subroutine transz


  !===========================================================================
  ! transxy
  !===========================================================================
  subroutine transxy(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fxy,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                     fyx,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                     ugdnvx,pgdnvx,gegdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                     ugdnvy,pgdnvy,gegdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                     gamc,gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3, &
                     srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                     grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                     rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3, &
                     hdt,cdtdx,cdtdy,ilo,ihi,jlo,jhi,kc,km,k3d)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &            
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_type, ppm_trace_grav, rot_period, ppm_trace_rot, &
                                   do_grav, do_rotation
    use eos_module

    implicit none
    
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer gd_l1,gd_l2,gd_l3,gd_h1,gd_h2,gd_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
    integer rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d
    
    double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fxy(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision fyx(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision gegdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gegdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gamc(gd_l1:gd_h1,gd_l2:gd_h2,gd_l3:gd_h3)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
    double precision rot(rt_l1:rt_h1,rt_l2:rt_h2,rt_l3:rt_h3,3)
    double precision hdt,cdtdx,cdtdy
    
    integer i, j
    integer n , nq
    
    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    double precision pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    double precision uxav, gexav, dgex, uyav, geyav, dgey
    double precision pgxpm, pgxmm, ugxpm, ugxmm, gegxpm, gegxmm, duxpm, pxavm, duxm, pxnewm, gexnewm
    double precision pgypm, pgymm, ugypm, ugymm, gegypm, gegymm, duypm, pyavm, duym, pynewm, geynewm
    double precision uxavm, gexavm, dgexm, uyavm, geyavm, dgeym
    double precision compr, compl, compnr, compnl
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             rrr = qp(i,j,kc,QRHO)
             rrl = qm(i,j,kc,QRHO)
             
             compr = rrr*qp(i,j,kc,nq)
             compl = rrl*qm(i,j,kc,nq)
             
             rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                          - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
             rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                          - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))

             compnr = compr - cdtdx*(fxy(i+1,j,kc,n) - fxy(i,j,kc,n)) &
                            - cdtdy*(fyx(i,j+1,kc,n) - fyx(i,j,kc,n))
             compnl = compl - cdtdx*(fxy(i+1,j,km,n) - fxy(i,j,km,n)) &
                            - cdtdy*(fyx(i,j+1,km,n) - fyx(i,j,km,n))
             
             qpo(i,j,kc,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d  ,nq)
             qmo(i,j,kc,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d-1,nq)
             
          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          !-------------------------------------------------------------------
          ! add the transverse xy and yx differences to the z-states for the
          ! fluid variables
          !-------------------------------------------------------------------          

          pgxp = pgdnvx(i+1,j,kc)
          pgxm = pgdnvx(i,j,kc)
          ugxp = ugdnvx(i+1,j,kc)          
          ugxm = ugdnvx(i,j,kc)
          gegxp = gegdnvx(i+1,j,kc)          
          gegxm = gegdnvx(i,j,kc)
          
          pgyp = pgdnvy(i,j+1,kc)
          pgym = pgdnvy(i,j,kc)
          ugyp = ugdnvy(i,j+1,kc)
          ugym = ugdnvy(i,j,kc)
          gegyp = gegdnvy(i,j+1,kc)
          gegym = gegdnvy(i,j,kc)
          
          pgxpm = pgdnvx(i+1,j,km)
          pgxmm = pgdnvx(i,j,km)
          ugxpm = ugdnvx(i+1,j,km)
          ugxmm = ugdnvx(i,j,km)
          gegxpm = gegdnvx(i+1,j,km)
          gegxmm = gegdnvx(i,j,km)
          
          pgypm = pgdnvy(i,j+1,km)
          pgymm = pgdnvy(i,j,km)
          ugypm = ugdnvy(i,j+1,km)
          ugymm = ugdnvy(i,j,km)
          gegypm = gegdnvy(i,j+1,km)
          gegymm = gegdnvy(i,j,km)

          duxp = pgxp*ugxp - pgxm*ugxm
          pxav = HALF*(pgxp+pgxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav-gamc(i,j,k3d))*dux - uxav*dgex )

          duxpm = pgxpm*ugxpm - pgxmm*ugxmm
          pxavm = HALF*(pgxpm+pgxmm)
          uxavm = HALF*(ugxpm+ugxmm)
          gexavm = HALF*(gegxpm+gegxmm)
          duxm = ugxpm-ugxmm
          dgexm = gegxpm-gegxmm
          pxnewm = cdtdx*(duxpm + pxavm*duxm*(gamc(i,j,k3d-1)-ONE))
          gexnewm = cdtdx*( (gexavm-ONE)*(gexavm-gamc(i,j,k3d-1))*duxm - uxavm*dgexm )
          
          duyp = pgyp*ugyp - pgym*ugym
          pyav = HALF*(pgyp+pgym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav-gamc(i,j,k3d))*duy - uyav*dgey )

          duypm = pgypm*ugypm - pgymm*ugymm
          pyavm = HALF*(pgypm+pgymm)
          uyavm = HALF*(ugypm+ugymm)
          geyavm = HALF*(gegypm+gegymm)
          duym = ugypm-ugymm
          dgeym = gegypm-gegymm
          pynewm = cdtdy*(duypm + pyavm*duym*(gamc(i,j,k3d-1)-ONE))
          geynewm = cdtdy*( (geyavm-ONE)*(geyavm-gamc(i,j,k3d-1))*duym - uyavm*dgeym )
          
          ! Convert to conservation form
          rrr = qp(i,j,kc,QRHO)
          rur = rrr*qp(i,j,kc,QU)
          rvr = rrr*qp(i,j,kc,QV)
          rwr = rrr*qp(i,j,kc,QW)
          ekenr = HALF*rrr*(qp(i,j,kc,QU)**2 + qp(i,j,kc,QV)**2 + &
               qp(i,j,kc,QW)**2)
          rer = qp(i,j,kc,QREINT) + ekenr
          
          rrl = qm(i,j,kc,QRHO)
          rul = rrl*qm(i,j,kc,QU)
          rvl = rrl*qm(i,j,kc,QV)
          rwl = rrl*qm(i,j,kc,QW)
          ekenl = HALF*rrl*(qm(i,j,kc,QU)**2 + qm(i,j,kc,QV)**2 + &
               qm(i,j,kc,QW)**2)
          rel = qm(i,j,kc,QREINT) + ekenl
          
          ! Add transverse predictor
          rrnewr = rrr - cdtdx*(fxy(i+1,j,kc,URHO) - fxy(i,j,kc,URHO)) &
                       - cdtdy*(fyx(i,j+1,kc,URHO) - fyx(i,j,kc,URHO))
          runewr = rur - cdtdx*(fxy(i+1,j,kc,UMX) - fxy(i,j,kc,UMX)) &
                       - cdtdy*(fyx(i,j+1,kc,UMX) - fyx(i,j,kc,UMX))
          rvnewr = rvr - cdtdx*(fxy(i+1,j,kc,UMY) - fxy(i,j,kc,UMY)) &
                       - cdtdy*(fyx(i,j+1,kc,UMY) - fyx(i,j,kc,UMY))
          rwnewr = rwr - cdtdx*(fxy(i+1,j,kc,UMZ) - fxy(i,j,kc,UMZ)) &
                       - cdtdy*(fyx(i,j+1,kc,UMZ) - fyx(i,j,kc,UMZ))
          renewr = rer - cdtdx*(fxy(i+1,j,kc,UEDEN) - fxy(i,j,kc,UEDEN)) &
                       - cdtdy*(fyx(i,j+1,kc,UEDEN) - fyx(i,j,kc,UEDEN))


          rrnewl = rrl - cdtdx*(fxy(i+1,j,km,URHO) - fxy(i,j,km,URHO)) &
                       - cdtdy*(fyx(i,j+1,km,URHO) - fyx(i,j,km,URHO))
          runewl = rul - cdtdx*(fxy(i+1,j,km,UMX) - fxy(i,j,km,UMX)) &
                       - cdtdy*(fyx(i,j+1,km,UMX) - fyx(i,j,km,UMX))
          rvnewl = rvl - cdtdx*(fxy(i+1,j,km,UMY) - fxy(i,j,km,UMY)) &
                       - cdtdy*(fyx(i,j+1,km,UMY) - fyx(i,j,km,UMY))
          rwnewl = rwl - cdtdx*(fxy(i+1,j,km,UMZ) - fxy(i,j,km,UMZ)) &
                       - cdtdy*(fyx(i,j+1,km,UMZ) - fyx(i,j,km,UMZ))
          renewl = rel - cdtdx*(fxy(i+1,j,km,UEDEN) - fxy(i,j,km,UEDEN)) &
                       - cdtdy*(fyx(i,j+1,km,UEDEN) - fyx(i,j,km,UEDEN))

          ! Reset to original value if adding transverse terms made density negative
          if (transverse_reset_density == 1) then
             if (rrnewr .lt. ZERO) then
                rrnewr = rrr 
                runewr = rur 
                rvnewr = rvr 
                rwnewr = rwr 
                renewr = rer 
             endif
             if (rrnewl .lt. ZERO) then
                rrnewl = rrl 
                runewl = rul 
                rvnewl = rvl 
                rwnewl = rwl 
                renewl = rel 
             endif
          endif

          rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
          rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl

          !-------------------------------------------------------------------
          ! qzpo state
          !-------------------------------------------------------------------          
          
          ! Convert back to primitive form
          qpo(i,j,kc,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
          qpo(i,j,kc,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k3d,QU) 
          qpo(i,j,kc,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k3d,QV) 
          qpo(i,j,kc,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k3d,QW) 

          ! note: we run the risk of (rho e) being negative here
          qpo(i,j,kc,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).
             
             if (qpo(i,j,kc,QREINT) .le. ZERO) then
                qpo(i,j,kc,QREINT) = qp(i,j,kc,QREINT) &
                     - cdtdx*(fxy(i+1,j,kc,UEINT) - fxy(i,j,kc,UEINT) + pxav*dux) &
                     - cdtdy*(fyx(i,j+1,kc,UEINT) - fyx(i,j,kc,UEINT) + pyav*duy) &
                     + hdt*srcQ(i,j,k3d,QREINT)
                
                ! if we are still negative, then we need to reset
                if (qpo(i,j,kc,QREINT) < ZERO) then
                   eos_state % rho = qpo(i,j,kc,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qpo(i,j,kc,QFS:QFS-1+nspec)
                   
                   call eos(eos_input_rt, eos_state)
                   
                   qpo(i,j,kc,QREINT) = qpo(i,j,kc,QRHO)*eos_state % e
                   qpo(i,j,kc,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then
             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qpo(i,j,kc,QRHO)
                eos_state % e   = qpo(i,j,kc,QREINT) / qpo(i,j,kc,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qpo(i,j,kc,QFS:QFS+nspec-1)
                
                call eos(eos_input_re, eos_state)

                pnewr = eos_state % p
                qpo(i,j,kc,QPRES ) = pnewr
                qpo(i,j,kc,QREINT) = eos_state % e * eos_state % rho    
             else
                ! add the transverse term to the p evolution eq here
                pnewr = qp(i,j,kc,QPRES) - pxnew - pynew
                qpo(i,j,kc,QPRES) = pnewr + hdt*srcQ(i,j,k3d,QPRES)
             endif

             qpo(i,j,kc,QPRES) = max(qpo(i,j,kc,QPRES),small_pres)

          else
             
             ! Update gammae with its transverse terms
             qpo(i,j,kc,QGAME) = qp(i,j,kc,QGAME) + gexnew + geynew

             ! and compute the p edge state from this and (rho e)
             qpo(i,j,kc,QPRES) = qpo(i,j,kc,QREINT)*(qpo(i,j,kc,QGAME)-ONE)
             qpo(i,j,kc,QPRES) = max(qpo(i,j,kc,QPRES), small_pres)

          endif


          !-------------------------------------------------------------------
          ! qzmo state
          !-------------------------------------------------------------------          

          qmo(i,j,kc,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k3d-1,QRHO)
          qmo(i,j,kc,QU    ) = runewl/rrnewl + hdt*srcQ(i,j,k3d-1,QU)
          qmo(i,j,kc,QV    ) = rvnewl/rrnewl + hdt*srcQ(i,j,k3d-1,QV)
          qmo(i,j,kc,QW    ) = rwnewl/rrnewl + hdt*srcQ(i,j,k3d-1,QW)

          ! note: we run the risk of (rho e) being negative here
          qmo(i,j,kc,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d-1,QREINT)

          if (transverse_reset_rhoe == 1) then
             ! If it is negative, reset the internal energy by using the discretized
             ! expression for updating (rho e).
             
             if (qmo(i,j,kc,QREINT) .le. ZERO) then
                qmo(i,j,kc,QREINT) = qm(i,j,kc,QREINT) &
                     - cdtdx*(fxy(i+1,j,km,UEINT) - fxy(i,j,km,UEINT) + pxavm*duxm) &
                     - cdtdy*(fyx(i,j+1,km,UEINT) - fyx(i,j,km,UEINT) + pyavm*duym) &
                     + hdt*srcQ(i,j,k3d-1,QREINT)
                
                ! if we are still negative, then we need to reset
                if (qmo(i,j,kc,QREINT) < ZERO) then
                   eos_state % rho = qmo(i,j,kc,QRHO)
                   eos_state % T = small_temp
                   eos_state % xn(:) = qmo(i,j,kc,QFS:QFS-1+nspec)
                   
                   call eos(eos_input_rt, eos_state)
                   
                   qmo(i,j,kc,QREINT) = qmo(i,j,kc,QRHO)*eos_state % e
                   qmo(i,j,kc,QPRES) = eos_state % p
                endif
             endif
          endif

          if (ppm_predict_gammae == 0) then

             ! Optionally, use the EOS to calculate the pressure.

             if (transverse_use_eos .eq. 1) then
                eos_state % rho = qmo(i,j,kc,QRHO)
                eos_state % e   = qmo(i,j,kc,QREINT) / qmo(i,j,kc,QRHO)
                eos_state % T   = small_temp
                eos_state % xn  = qmo(i,j,kc,QFS:QFS+nspec-1)
                
                call eos(eos_input_re, eos_state)

                pnewl = eos_state % p
                qmo(i,j,kc,QPRES ) = pnewl
                qmo(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             else
                ! add the transverse term to the p evolution eq here
                pnewl = qm(i,j,kc,QPRES) - pxnewm - pynewm
                qmo(i,j,kc,QPRES) = pnewl + hdt*srcQ(i,j,k3d-1,QPRES)
             endif

             qmo(i,j,kc,QPRES) = max(qmo(i,j,kc,QPRES),small_pres)
          
          else

             ! Update gammae with its transverse terms
             qmo(i,j,kc,QGAME) = qm(i,j,kc,QGAME) + gexnewm + geynewm
                 
             ! and compute the p edge state from this and (rho e)
             qmo(i,j,kc,QPRES) = qmo(i,j,kc,QREINT)*(qmo(i,j,kc,QGAME)-ONE)
             qmo(i,j,kc,QPRES) = max(qmo(i,j,kc,QPRES), small_pres)

          endif

       enddo
    enddo


    ! if ppm_trace_grav == 1, then we already added the piecewise parabolic traced
    ! gravity to the normal edge states
    if ((do_grav .eq. 1) .and. (ppm_trace_grav == 0 .or. ppm_type == 0)) then
       do j = jlo, jhi 
          do i = ilo, ihi          
             qpo(i,j,kc,QU    ) = qpo(i,j,kc,QU    ) + hdt*grav(i,j,k3d,1)
             qpo(i,j,kc,QV    ) = qpo(i,j,kc,QV    ) + hdt*grav(i,j,k3d,2)
             qpo(i,j,kc,QW    ) = qpo(i,j,kc,QW    ) + hdt*grav(i,j,k3d,3)
             
             qmo(i,j,kc,QU    ) = qmo(i,j,kc,QU    ) + hdt*grav(i,j,k3d-1,1)
             qmo(i,j,kc,QV    ) = qmo(i,j,kc,QV    ) + hdt*grav(i,j,k3d-1,2)
             qmo(i,j,kc,QW    ) = qmo(i,j,kc,QW    ) + hdt*grav(i,j,k3d-1,3)
          enddo
       enddo
    endif

    ! if ppm_trace_rot == 1, then we already added the piecewise parabolic traced
    ! rotation to the normal edge states
    if ((do_rotation .eq. 1) .and. (ppm_trace_rot == 0 .or. ppm_type == 0)) then
       do j = jlo, jhi
          do i = ilo, ihi
             qpo(i,j,kc,QU    ) = qpo(i,j,kc,QU    ) + hdt*rot(i,j,k3d,1)
             qpo(i,j,kc,QV    ) = qpo(i,j,kc,QV    ) + hdt*rot(i,j,k3d,2)
             qpo(i,j,kc,QW    ) = qpo(i,j,kc,QW    ) + hdt*rot(i,j,k3d,3)

             qmo(i,j,kc,QU    ) = qmo(i,j,kc,QU    ) + hdt*rot(i,j,k3d-1,1)
             qmo(i,j,kc,QV    ) = qmo(i,j,kc,QV    ) + hdt*rot(i,j,k3d-1,2)
             qmo(i,j,kc,QW    ) = qmo(i,j,kc,QW    ) + hdt*rot(i,j,k3d-1,3)
          enddo
       enddo
    endif
    
  end subroutine transxy

  
  !===========================================================================
  ! transxz
  !===========================================================================
  subroutine transxz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fxz,fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3, &
                     fzx,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                     ugdnvx,pgdnvx,gegdnvx,pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3, &
                     ugdnvz,pgdnvz,gegdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                     gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                     srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                     grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                     rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3, &
                     hdt,cdtdx,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_type, ppm_trace_grav, rot_period, ppm_trace_rot, &
                                   do_grav, do_rotation
    use eos_module

    implicit none      
    
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fx_l1,fx_l2,fx_l3,fx_h1,fx_h2,fx_h3
    integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
    integer pgdx_l1,pgdx_l2,pgdx_l3,pgdx_h1,pgdx_h2,pgdx_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
    integer rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d
    
    double precision  qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision  qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fxz(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
    double precision fzx(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
    double precision ugdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision pgdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision gegdnvx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2,pgdx_l3:pgdx_h3)
    double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision gegdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
    double precision rot(rt_l1:rt_h1,rt_l2:rt_h2,rt_l3:rt_h3,3)
    double precision hdt,cdtdx,cdtdz
    
    integer i, j
    integer n, nq
    
    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pgxp, pgxm, ugxp, ugxm, gegxp, gegxm, duxp, pxav, dux, pxnew, gexnew
    double precision pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    double precision uxav, gexav, dgex, uzav, gezav, dgez
    double precision compr, compl, compnr, compnl, drr, dcompn
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------    
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             drr    = - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                      - cdtdz*(fzx(i  ,j,kc,URHO) - fzx(i,j,km,URHO))
             dcompn = - cdtdx*(fxz(i+1,j,km,n   ) - fxz(i,j,km,n)) &
                      - cdtdz*(fzx(i  ,j,kc,n   ) - fzx(i,j,km,n))

             if (j.ge.jlo+1) then
                rrr = qp(i,j,km,QRHO)
                compr = rrr*qp(i,j,km,nq)

                rrnewr = rrr + drr
                compnr = compr + dcompn

                qpo(i,j  ,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             end if

             if (j.le.jhi-1) then
                rrl = qm(i,j+1,km,QRHO)
                compl = rrl*qm(i,j+1,km,nq)

                rrnewr = rrl + drr
                compnl = compl + dcompn

                qmo(i,j+1,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)
             end if

          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          !-------------------------------------------------------------------
          ! add the transverse xz and zx differences to the y-states for the
          ! fluid variables
          !-------------------------------------------------------------------            
          pgxp = pgdnvx(i+1,j,km)
          pgxm = pgdnvx(i,j,km)
          ugxp = ugdnvx(i+1,j,km)
          ugxm = ugdnvx(i,j,km)
          gegxp = gegdnvx(i+1,j,km)
          gegxm = gegdnvx(i,j,km)
          
          pgzp = pgdnvz(i,j,kc)
          pgzm = pgdnvz(i,j,km)
          ugzp = ugdnvz(i,j,kc)
          ugzm = ugdnvz(i,j,km)
          gegzp = gegdnvz(i,j,kc)
          gegzm = gegdnvz(i,j,km)

          duxp = pgxp*ugxp - pgxm*ugxm
          pxav = HALF*(pgxp+pgxm)
          uxav = HALF*(ugxp+ugxm)
          gexav = HALF*(gegxp+gegxm)
          dux = ugxp-ugxm
          dgex = gegxp-gegxm
          pxnew = cdtdx*(duxp + pxav*dux*(gamc(i,j,k3d)-ONE))
          gexnew = cdtdx*( (gexav-ONE)*(gexav-gamc(i,j,k3d))*dux - uxav*dgex )

          duzp = pgzp*ugzp - pgzm*ugzm
          pzav = HALF*(pgzp+pgzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav-gamc(i,j,k3d))*duz - uzav*dgez )

          !-------------------------------------------------------------------
          ! qypo state
          !-------------------------------------------------------------------

          ! Convert back to primitive form
          if (j.ge.jlo+1) then
             ! Convert to conservation form
             rrr = qp(i,j,km,QRHO)
             rur = rrr*qp(i,j,km,QU)
             rvr = rrr*qp(i,j,km,QV)
             rwr = rrr*qp(i,j,km,QW)
             ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + qp(i,j,km,QW)**2)
             rer = qp(i,j,km,QREINT) + ekenr
             
             ! Add transverse predictor
             rrnewr = rrr - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                          - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
             runewr = rur - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                          - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
             rvnewr = rvr - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                          - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
             rwnewr = rwr - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                          - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
             renewr = rer - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                          - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewr .lt. ZERO) then
                   rrnewr = rrr 
                   runewr = rur
                   rvnewr = rvr 
                   rwnewr = rwr 
                   renewr = rer 
                endif
             end if

             qpo(i,j,km,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
             qpo(i,j,km,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k3d,QU)
             qpo(i,j,km,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k3d,QV)
             qpo(i,j,km,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k3d,QW)

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,km,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).

                if (qpo(i,j,km,QREINT) .le. ZERO) then
                   qpo(i,j,km,QREINT) = qp(i,j,km,QREINT) &
                        - cdtdx*(fxz(i+1,j,km,UEINT) - fxz(i,j,km,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i  ,j,kc,UEINT) - fzx(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                   
                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,km,QREINT) < ZERO) then
                      eos_state % rho = qpo(i,j,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qpo(i,j,km,QREINT) = qpo(i,j,km,QRHO)*eos_state % e
                      qpo(i,j,km,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then
                
                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then

                   eos_state % rho = qpo(i,j,km,QRHO)
                   eos_state % e   = qpo(i,j,km,QREINT) / qpo(i,j,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qpo(i,j,km,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewr = eos_state % p
                   qpo(i,j,km,QPRES ) = pnewr
                   qpo(i,j,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,km,QPRES) - pxnew - pznew
                   qpo(i,j,km,QPRES) = pnewr + hdt*srcQ(i,j,k3d,QPRES)
                endif
                
                qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES),small_pres)

             else
                
                ! Update gammae with its transverse terms
                qpo(i,j,km,QGAME) = qp(i,j,km,QGAME) + gexnew + geznew

                ! and compute the p edge state from this and (rho e)
                qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT)*(qpo(i,j,km,QGAME)-ONE)
                qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES), small_pres)

             endif

          end if


          !-------------------------------------------------------------------
          ! qymo state
          !-------------------------------------------------------------------          

          if (j.le.jhi-1) then
             ! Convert to conservation form
             rrl = qm(i,j+1,km,QRHO)
             rul = rrl*qm(i,j+1,km,QU)
             rvl = rrl*qm(i,j+1,km,QV)
             rwl = rrl*qm(i,j+1,km,QW)
             ekenl = HALF*rrl*(qm(i,j+1,km,QU)**2 + qm(i,j+1,km,QV)**2 + qm(i,j+1,km,QW)**2)
             rel = qm(i,j+1,km,QREINT) + ekenl
             
             ! Add transverse predictor
             rrnewl = rrl - cdtdx*(fxz(i+1,j,km,URHO) - fxz(i,j,km,URHO)) &
                          - cdtdz*(fzx(i,j,kc,URHO) - fzx(i,j,km,URHO))
             runewl = rul - cdtdx*(fxz(i+1,j,km,UMX) - fxz(i,j,km,UMX)) &
                          - cdtdz*(fzx(i,j,kc,UMX) - fzx(i,j,km,UMX))
             rvnewl = rvl - cdtdx*(fxz(i+1,j,km,UMY) - fxz(i,j,km,UMY)) &
                          - cdtdz*(fzx(i,j,kc,UMY) - fzx(i,j,km,UMY))
             rwnewl = rwl - cdtdx*(fxz(i+1,j,km,UMZ) - fxz(i,j,km,UMZ)) &
                          - cdtdz*(fzx(i,j,kc,UMZ) - fzx(i,j,km,UMZ))
             renewl = rel - cdtdx*(fxz(i+1,j,km,UEDEN) - fxz(i,j,km,UEDEN)) &
                          - cdtdz*(fzx(i,j,kc,UEDEN) - fzx(i,j,km,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewl .lt. ZERO) then
                   rrnewl = rrl 
                   runewl = rul 
                   rvnewl = rvl 
                   rwnewl = rwl 
                   renewl = rel 
                endif
             endif
             
             qmo(i,j+1,km,QRHO  ) = rrnewl        + hdt*srcQ(i,j,k3d,QRHO)
             qmo(i,j+1,km,QU    ) = runewl/rrnewl + hdt*srcQ(i,j,k3d,QU)
             qmo(i,j+1,km,QV    ) = rvnewl/rrnewl + hdt*srcQ(i,j,k3d,QV)
             qmo(i,j+1,km,QW    ) = rwnewl/rrnewl + hdt*srcQ(i,j,k3d,QW)

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i,j+1,km,QREINT) = renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qmo(i,j+1,km,QREINT) .le. ZERO) then
                   qmo(i,j+1,km,QREINT) = qm(i,j+1,km,QREINT) &
                        - cdtdx*(fxz(i+1,j,km,UEINT) - fxz(i,j,km,UEINT) + pxav*dux) &
                        - cdtdz*(fzx(i,j,kc,UEINT) - fzx(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                   
                   ! if we are still negative, then we need to reset
                   if (qmo(i,j+1,km,QREINT) < ZERO) then
                      eos_state % rho = qmo(i,j+1,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i,j+1,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)

                      qmo(i,j+1,km,QREINT) = qmo(i,j+1,km,QRHO)*eos_state % e
                      qmo(i,j+1,km,QPRES) = eos_state % p
                   endif
                endif
             endif
             
             
             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qmo(i,j+1,km,QRHO)
                   eos_state % e   = qmo(i,j+1,km,QREINT) / qmo(i,j+1,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qmo(i,j+1,km,QFS:QFS+nspec-1)
                   
                   call eos(eos_input_re, eos_state)

                   pnewl = eos_state % p
                   qmo(i,j+1,km,QPRES ) = pnewl
                   qmo(i,j+1,km,QREINT) = eos_state % e * eos_state % rho
                else
                   pnewl = qm(i,j+1,km,QPRES) - pxnew - pznew
                   qmo(i,j+1,km,QPRES) = pnewl + hdt*srcQ(i,j,k3d,QPRES)
                endif

                qmo(i,j+1,km,QPRES) = max(qmo(i,j+1,km,QPRES),small_pres)
             
             else

                ! Update gammae with its transverse terms
                qmo(i,j+1,km,QGAME) = qm(i,j+1,km,QGAME) + gexnew + geznew

                ! and compute the p edge state from this and (rho e)
                qmo(i,j+1,km,QPRES) = qmo(i,j+1,km,QREINT)*(qmo(i,j+1,km,QGAME)-ONE)
                qmo(i,j+1,km,QPRES) = max(qmo(i,j+1,km,QPRES), small_pres)

             endif

          endif
          
       enddo
    enddo

    ! if ppm_trace_grav == 1, then we already added the piecewise parabolic traced
    ! gravity to the normal edge states
    if ((do_grav .eq. 1) .and. (ppm_trace_grav == 0 .or. ppm_type == 0)) then
       do j = jlo+1, jhi 
          do i = ilo, ihi 
             qpo(i,j,km,QU    ) = qpo(i,j,km,QU    ) + hdt*grav(i,j,k3d,1)
             qpo(i,j,km,QV    ) = qpo(i,j,km,QV    ) + hdt*grav(i,j,k3d,2)
             qpo(i,j,km,QW    ) = qpo(i,j,km,QW    ) + hdt*grav(i,j,k3d,3)
          end do
       end do
       do j = jlo, jhi-1 
          do i = ilo, ihi 
             qmo(i,j+1,km,QU    ) = qmo(i,j+1,km,QU    ) + hdt*grav(i,j,k3d,1)
             qmo(i,j+1,km,QV    ) = qmo(i,j+1,km,QV    ) + hdt*grav(i,j,k3d,2)
             qmo(i,j+1,km,QW    ) = qmo(i,j+1,km,QW    ) + hdt*grav(i,j,k3d,3)
          enddo
       enddo
    endif

    ! if ppm_trace_rot == 1, then we already added the piecewise parabolic traced
    ! rotation to the normal edge states
    if ((do_rotation .eq. 1) .and. (ppm_trace_rot == 0 .or. ppm_type == 0)) then
       do j = jlo+1, jhi 
          do i = ilo, ihi 
             qpo(i,j,km,QU    ) = qpo(i,j,km,QU    ) + hdt*rot(i,j,k3d,1)
             qpo(i,j,km,QV    ) = qpo(i,j,km,QV    ) + hdt*rot(i,j,k3d,2)
             qpo(i,j,km,QW    ) = qpo(i,j,km,QW    ) + hdt*rot(i,j,k3d,3)
          end do
       end do
       do j = jlo, jhi-1 
          do i = ilo, ihi              
             qmo(i,j+1,km,QU    ) = qmo(i,j+1,km,QU    ) + hdt*rot(i,j,k3d,1)
             qmo(i,j+1,km,QV    ) = qmo(i,j+1,km,QV    ) + hdt*rot(i,j,k3d,2)
             qmo(i,j+1,km,QW    ) = qmo(i,j+1,km,QW    ) + hdt*rot(i,j,k3d,3)
          enddo
       enddo
    endif
    
  end subroutine transxz


  !===========================================================================
  ! transyz
  !===========================================================================
  subroutine transyz(qm,qmo,qp,qpo,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     fyz,fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3, &
                     fzy,fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3, &
                     ugdnvy,pgdnvy,gegdnvy,pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3, &
                     ugdnvz,pgdnvz,gegdnvz,pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3, &
                     gamc,gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3, &
                     srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,&
                     grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                     rot,rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3, &
                     hdt,cdtdy,cdtdz,ilo,ihi,jlo,jhi,km,kc,k3d)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QGAME, QESGS, QFA, QFS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, &
                                   nadv, small_pres, small_temp, &
                                   npassive, upass_map, qpass_map, &
                                   ppm_predict_gammae, &
                                   transverse_use_eos, transverse_reset_density, transverse_reset_rhoe, &
                                   ppm_type, ppm_trace_grav, rot_period, ppm_trace_rot, &
                                   do_grav, do_rotation
    use eos_module

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer fy_l1,fy_l2,fy_l3,fy_h1,fy_h2,fy_h3
    integer fz_l1,fz_l2,fz_l3,fz_h1,fz_h2,fz_h3
    integer pgdy_l1,pgdy_l2,pgdy_l3,pgdy_h1,pgdy_h2,pgdy_h3
    integer pgdz_l1,pgdz_l2,pgdz_l3,pgdz_h1,pgdz_h2,pgdz_h3
    integer gc_l1,gc_l2,gc_l3,gc_h1,gc_h2,gc_h3
    integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
    integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
    integer rt_l1,rt_l2,rt_l3,rt_h1,rt_h2,rt_h3
    integer ilo,ihi,jlo,jhi,km,kc,k3d
    
    double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision fyz(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
    double precision fzy(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
    double precision ugdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision pgdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision gegdnvy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2,pgdy_l3:pgdy_h3)
    double precision ugdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision pgdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision gegdnvz(pgdz_l1:pgdz_h1,pgdz_l2:pgdz_h2,pgdz_l3:pgdz_h3)
    double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2,gc_l3:gc_h3)
    double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
    double precision rot(rt_l1:rt_h1,rt_l2:rt_h2,rt_l3:rt_h3,3)
    double precision hdt,cdtdy,cdtdz
    
    integer i, j
    integer n, nq
    
    double precision rrr, rur, rvr, rwr, rer, ekenr, rhoekenr
    double precision rrl, rul, rvl, rwl, rel, ekenl, rhoekenl
    double precision rrnewr, runewr, rvnewr, rwnewr, renewr
    double precision rrnewl, runewl, rvnewl, rwnewl, renewl
    double precision pnewr, pnewl
    double precision pgyp, pgym, ugyp, ugym, gegyp, gegym, duyp, pyav, duy, pynew, geynew
    double precision pgzp, pgzm, ugzp, ugzm, gegzp, gegzm, duzp, pzav, duz, pznew, geznew
    double precision uyav, geyav, dgey, uzav, gezav, dgez
    double precision compr, compl, compnr, compnl
    double precision drr, dcompn
    
    integer ipassive

    type (eos_t) :: eos_state

    !-------------------------------------------------------------------------
    ! update all of the passively-advected quantities with the
    ! transerse term and convert back to the primitive quantity
    !-------------------------------------------------------------------------

    do ipassive = 1,npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)
       do j = jlo, jhi
          do i = ilo, ihi
             
             drr    = - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                      - cdtdz*(fzy(i,j  ,kc,URHO) - fzy(i,j,km,URHO))
             dcompn = - cdtdy*(fyz(i,j+1,km,n   ) - fyz(i,j,km,n)) &
                      - cdtdz*(fzy(i,j  ,kc,n   ) - fzy(i,j,km,n))

             if (i.ge.ilo+1) then
                rrr = qp(i,j,km,QRHO)
                compr = rrr*qp(i,j,km,nq)

                rrnewr = rrr +drr
                compnr = compr +dcompn

                qpo(i  ,j,km,nq) = compnr/rrnewr + hdt*srcQ(i,j,k3d,nq)
             end if

             if (i.le.ihi-1) then
                rrl = qm(i+1,j,km,QRHO)
                compl = rrl*qm(i+1,j,km,nq)

                rrnewl = rrl + drr
                compnl = compl +dcompn

                qmo(i+1,j,km,nq) = compnl/rrnewl + hdt*srcQ(i,j,k3d,nq)
             end if
          enddo
       enddo
    enddo

    do j = jlo, jhi 
       do i = ilo, ihi 

          !-------------------------------------------------------------------
          ! add the transverse yz and zy differences to the x-states for the 
          ! fluid variables
          !-------------------------------------------------------------------          

          pgyp = pgdnvy(i,j+1,km)
          pgym = pgdnvy(i,j,km)
          ugyp = ugdnvy(i,j+1,km)
          ugym = ugdnvy(i,j,km)
          gegyp = gegdnvy(i,j+1,km)
          gegym = gegdnvy(i,j,km)
          
          pgzp = pgdnvz(i,j,kc)
          pgzm = pgdnvz(i,j,km)
          ugzp = ugdnvz(i,j,kc)
          ugzm = ugdnvz(i,j,km)
          gegzp = gegdnvz(i,j,kc)
          gegzm = gegdnvz(i,j,km)

          duyp = pgyp*ugyp - pgym*ugym
          pyav = HALF*(pgyp+pgym)
          uyav = HALF*(ugyp+ugym)
          geyav = HALF*(gegyp+gegym)
          duy = ugyp-ugym
          dgey = gegyp-gegym
          pynew = cdtdy*(duyp + pyav*duy*(gamc(i,j,k3d)-ONE))
          geynew = cdtdy*( (geyav-ONE)*(geyav-gamc(i,j,k3d))*duy - uyav*dgey )
          
          duzp = pgzp*ugzp - pgzm*ugzm
          pzav = HALF*(pgzp+pgzm)
          uzav = HALF*(ugzp+ugzm)
          gezav = HALF*(gegzp+gegzm)
          duz = ugzp-ugzm
          dgez = gegzp-gegzm
          pznew = cdtdz*(duzp + pzav*duz*(gamc(i,j,k3d)-ONE))
          geznew = cdtdz*( (gezav-ONE)*(gezav-gamc(i,j,k3d))*duz - uzav*dgez )
          

          !-------------------------------------------------------------------
          ! qxpo state
          !-------------------------------------------------------------------
          
          ! Convert back to primitive form
          if (i.ge.ilo+1) then
             ! Convert to conservation form
             rrr = qp(i,j,km,QRHO)
             rur = rrr*qp(i,j,km,QU)
             rvr = rrr*qp(i,j,km,QV)
             rwr = rrr*qp(i,j,km,QW)
             ekenr = HALF*rrr*(qp(i,j,km,QU)**2 + qp(i,j,km,QV)**2 + &
                  qp(i,j,km,QW)**2)
             rer = qp(i,j,km,QREINT) + ekenr
             
             ! Add transverse predictor
             rrnewr = rrr - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                          - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
             runewr = rur - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                          - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
             rvnewr = rvr - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                          - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
             rwnewr = rwr - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                          - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
             renewr = rer - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                          - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewr .lt. ZERO) then
                   rrnewr = rrr 
                   runewr = rur 
                   rvnewr = rvr 
                   rwnewr = rwr 
                   renewr = rer 
                endif
             end if
                
             qpo(i,j,km,QRHO  ) = rrnewr        + hdt*srcQ(i,j,k3d,QRHO)
             qpo(i,j,km,QU    ) = runewr/rrnewr + hdt*srcQ(i,j,k3d,QU)
             qpo(i,j,km,QV    ) = rvnewr/rrnewr + hdt*srcQ(i,j,k3d,QV)
             qpo(i,j,km,QW    ) = rwnewr/rrnewr + hdt*srcQ(i,j,k3d,QW)

             ! note: we run the risk of (rho e) being negative here
             rhoekenr = HALF*(runewr**2 + rvnewr**2 + rwnewr**2)/rrnewr
             qpo(i,j,km,QREINT) = renewr - rhoekenr + hdt*srcQ(i,j,k3d,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qpo(i,j,km,QREINT) .le. ZERO) then
                   qpo(i,j,km,QREINT) = qp(i,j,km,QREINT) &
                        - cdtdy*(fyz(i,j+1,km,UEINT) - fyz(i,j,km,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,kc,UEINT) - fzy(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)
                   
                   ! if we are still negative, then we need to reset
                   if (qpo(i,j,km,QREINT) .le. ZERO) then
                      eos_state % rho = qpo(i,j,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qpo(i,j,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)
                      
                      qpo(i,j,km,QREINT) = qpo(i,j,km,QRHO)*eos_state % e
                      qpo(i,j,km,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS to calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qpo(i,j,km,QRHO)
                   eos_state % e   = qpo(i,j,km,QREINT) / qpo(i,j,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qpo(i,j,km,QFS:QFS+nspec-1)

                   call eos(eos_input_re, eos_state)

                   pnewr = eos_state % p
                   qpo(i,j,km,QPRES ) = pnewr
                   qpo(i,j,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewr = qp(i,j,km,QPRES) - pynew - pznew
                   qpo(i,j,km,QPRES) = pnewr + hdt*srcQ(i,j,k3d,QPRES)
                endif

                qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES),small_pres)

             else

                ! Update gammae with its transverse terms
                qpo(i,j,km,QGAME) = qp(i,j,km,QGAME) + geynew + geznew

                ! and compute the p edge state from this and (rho e)
                qpo(i,j,km,QPRES) = qpo(i,j,km,QREINT)*(qpo(i,j,km,QGAME)-ONE)
                qpo(i,j,km,QPRES) = max(qpo(i,j,km,QPRES), small_pres)

             end if

          endif

          !-------------------------------------------------------------------
          ! qxmo state
          !-------------------------------------------------------------------

          if (i.le.ihi-1) then
             ! Convert to conservation form
             rrl = qm(i+1,j,km,QRHO)
             rul = rrl*qm(i+1,j,km,QU)
             rvl = rrl*qm(i+1,j,km,QV)
             rwl = rrl*qm(i+1,j,km,QW)
             ekenl = HALF*rrl*(qm(i+1,j,km,QU)**2 + qm(i+1,j,km,QV)**2 + &
                  qm(i+1,j,km,QW)**2)
             rel = qm(i+1,j,km,QREINT) + ekenl
             
             ! Add transverse predictor
             rrnewl = rrl - cdtdy*(fyz(i,j+1,km,URHO) - fyz(i,j,km,URHO)) &
                          - cdtdz*(fzy(i,j,kc,URHO) - fzy(i,j,km,URHO))
             runewl = rul - cdtdy*(fyz(i,j+1,km,UMX) - fyz(i,j,km,UMX)) &
                          - cdtdz*(fzy(i,j,kc,UMX) - fzy(i,j,km,UMX))
             rvnewl = rvl - cdtdy*(fyz(i,j+1,km,UMY) - fyz(i,j,km,UMY)) &
                          - cdtdz*(fzy(i,j,kc,UMY) - fzy(i,j,km,UMY))
             rwnewl = rwl - cdtdy*(fyz(i,j+1,km,UMZ) - fyz(i,j,km,UMZ)) &
                          - cdtdz*(fzy(i,j,kc,UMZ) - fzy(i,j,km,UMZ))
             renewl = rel - cdtdy*(fyz(i,j+1,km,UEDEN) - fyz(i,j,km,UEDEN)) &
                          - cdtdz*(fzy(i,j,kc,UEDEN) - fzy(i,j,km,UEDEN))

             ! Reset to original value if adding transverse terms made density negative
             if (transverse_reset_density == 1) then
                if (rrnewl .lt. ZERO) then
                   rrnewl = rrl 
                   runewl = rul 
                   rvnewl = rvl 
                   rwnewl = rwl 
                   renewl = rel 
                endif
             endif

             qmo(i+1,j,km,QRHO   ) = rrnewl        + hdt*srcQ(i,j,k3d,QRHO)
             qmo(i+1,j,km,QU     ) = runewl/rrnewl + hdt*srcQ(i,j,k3d,QU)
             qmo(i+1,j,km,QV     ) = rvnewl/rrnewl + hdt*srcQ(i,j,k3d,QV)
             qmo(i+1,j,km,QW     ) = rwnewl/rrnewl + hdt*srcQ(i,j,k3d,QW)

             ! note: we run the risk of (rho e) being negative here
             rhoekenl = HALF*(runewl**2 + rvnewl**2 + rwnewl**2)/rrnewl
             qmo(i+1,j,km,QREINT ) = renewl - rhoekenl + hdt*srcQ(i,j,k3d,QREINT)

             if (transverse_reset_rhoe == 1) then
                ! If it is negative, reset the internal energy by using the discretized
                ! expression for updating (rho e).
                
                if (qmo(i+1,j,km,QREINT) .le. ZERO) then
                   qmo(i+1,j,km,QREINT ) = qm(i+1,j,km,QREINT) &
                        - cdtdy*(fyz(i,j+1,km,UEINT) - fyz(i,j,km,UEINT) + pyav*duy) &
                        - cdtdz*(fzy(i,j  ,kc,UEINT) - fzy(i,j,km,UEINT) + pzav*duz) &
                        + hdt*srcQ(i,j,k3d,QREINT)

                   ! if we are still negative, then we need to reset
                   if (qmo(i+1,j,km,QREINT) < ZERO) then
                      eos_state % rho = qmo(i+1,j,km,QRHO)
                      eos_state % T = small_temp
                      eos_state % xn(:) = qmo(i+1,j,km,QFS:QFS-1+nspec)
                      
                      call eos(eos_input_rt, eos_state)

                      qmo(i+1,j,km,QREINT) = qmo(i+1,j,km,QRHO)*eos_state % e
                      qmo(i+1,j,km,QPRES) = eos_state % p
                   endif
                endif
             endif

             if (ppm_predict_gammae == 0) then

                ! Optionally, use the EOS To calculate the pressure.

                if (transverse_use_eos .eq. 1) then
                   eos_state % rho = qmo(i+1,j,km,QRHO)
                   eos_state % e   = qmo(i+1,j,km,QREINT) / qmo(i+1,j,km,QRHO)
                   eos_state % T   = small_temp
                   eos_state % xn  = qmo(i+1,j,km,QFS:QFS+nspec-1)
                
                   call eos(eos_input_re, eos_state)
                   
                   pnewl = eos_state % p
                   qmo(i+1,j,km,QPRES ) = pnewl
                   qmo(i+1,j,km,QREINT) = eos_state % e * eos_state % rho
                else
                   ! add the transverse term to the p evolution eq here
                   pnewl = qm(i+1,j,km,QPRES) - pynew - pznew
                   qmo(i+1,j,km,QPRES  ) = pnewl + hdt*srcQ(i,j,k3d,QPRES)
                endif

                qmo(i+1,j,km,QPRES  ) = max(qmo(i+1,j,km,QPRES),small_pres)

             else
                
                ! Update gammae with its transverse terms
                qmo(i+1,j,km,QGAME) = qm(i+1,j,km,QGAME) + geynew + geznew

                ! and compute the p edge state from this and (rho e)
                qmo(i+1,j,km,QPRES) = qmo(i+1,j,km,QREINT)*(qmo(i+1,j,km,QGAME)-ONE)
                qmo(i+1,j,km,QPRES) = max(qmo(i+1,j,km,QPRES), small_pres)

             end if

          endif

       enddo
    enddo
    
    ! if ppm_trace_grav == 1, then we already added the piecewise parabolic traced
    ! gravity to the normal edge states
    if ((do_grav .eq. 1) .and. (ppm_trace_grav == 0 .or. ppm_type == 0)) then
       do j = jlo, jhi 
          do i = ilo+1, ihi 
             qpo(i,j,km,QU    ) = qpo(i,j,km,QU    ) + hdt*grav(i,j,k3d,1)
             qpo(i,j,km,QV    ) = qpo(i,j,km,QV    ) + hdt*grav(i,j,k3d,2)
             qpo(i,j,km,QW    ) = qpo(i,j,km,QW    ) + hdt*grav(i,j,k3d,3)
          end do
          do i = ilo, ihi-1
             qmo(i+1,j,km,QU     ) = qmo(i+1,j,km,QU     ) + hdt*grav(i,j,k3d,1)
             qmo(i+1,j,km,QV     ) = qmo(i+1,j,km,QV     ) + hdt*grav(i,j,k3d,2)
             qmo(i+1,j,km,QW     ) = qmo(i+1,j,km,QW     ) + hdt*grav(i,j,k3d,3)
          enddo
       enddo
    endif

    ! if ppm_trace_rot == 1, then we already added the piecewise parabolic traced
    ! rotation to the normal edge states
    if ((do_rotation .eq. 1) .and. (ppm_trace_rot == 0 .or. ppm_type == 0)) then
       do j = jlo, jhi 
          do i = ilo+1, ihi 
             qpo(i,j,km,QU    ) = qpo(i,j,km,QU    ) + hdt*rot(i,j,k3d,1)
             qpo(i,j,km,QV    ) = qpo(i,j,km,QV    ) + hdt*rot(i,j,k3d,2)
             qpo(i,j,km,QW    ) = qpo(i,j,km,QW    ) + hdt*rot(i,j,k3d,3)
          end do
          do i = ilo, ihi-1             
             qmo(i+1,j,km,QU     ) = qmo(i+1,j,km,QU     ) + hdt*rot(i,j,k3d,1)
             qmo(i+1,j,km,QV     ) = qmo(i+1,j,km,QV     ) + hdt*rot(i,j,k3d,2)
             qmo(i+1,j,km,QW     ) = qmo(i+1,j,km,QW     ) + hdt*rot(i,j,k3d,3)
          enddo
       enddo
    endif
    
  end subroutine transyz

end module transverse_module
