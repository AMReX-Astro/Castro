module riemann_module

  implicit none

  private

  public riemannus, riemanncg

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_h1,uflx_h2, &
                       pgdnv,pg_l1,pg_l2,pg_h1,pg_h2, &
                       ugdnv,ug_l1,ug_l2,ug_h1,ug_h2, &
                       idir,ilo1,ihi1,ilo2,ihi2)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo,Symmetry
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, &
                                   QPRES, QREINT, QFA, QFS, &
                                   QFX, URHO, UMX, UMY, UEDEN, UEINT, &
                                   UFA, UFS, UFX, &
                                   nadv, small_dens, small_pres

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_h1,uflx_h2
    integer :: ug_l1,ug_l2,ug_h1,ug_h2
    integer :: pg_l1,pg_l2,pg_h1,pg_h2
    integer :: idir,ilo1,ihi1,ilo2,ihi2
    integer :: i,j,ilo,jlo,ihi,jhi

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision ::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: ugdnv(ug_l1:ug_h1,ug_l2:ug_h2)
    double precision :: pgdnv(pg_l1:pg_h1,pg_l2:pg_h2)

    integer :: n, nq
    integer :: iadv, ispec, iaux
    
    double precision :: rgdnv,vgdnv,ustar,gamgdnv
    double precision :: rl, ul, vl, pl, rel
    double precision :: rr, ur, vr, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, gamco, entho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

    double precision :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    double precision :: zm, zp
    double precision :: denom, dpditer, dpjmp
    double precision :: gamc_bar, game_bar
    double precision :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter
    integer, parameter :: iter_max= 8
    double precision, parameter :: tol = 1.d-5
    double precision :: err

    logical :: converged

    double precision :: pstnm1
    double precision :: taul, taur, tauo
    double precision :: ustarm, ustarp, ustnm1, ustnp1

    double precision, parameter :: weakwv = 1.d-3

    !************************************************************
    !  set min/max based on normal direction
    if(idir.eq.1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1
    endif


    do j = jlo, jhi
       do i = ilo, ihi

          ! left state
          rl = max(ql(i,j,QRHO),small_dens)
          
          ! pick left velocities based on direction
          if(idir.eq.1) then
             ul = ql(i,j,QU)
             vl = ql(i,j,QV)
          else
             ul = ql(i,j,QV)
             vl = ql(i,j,QU)
          endif
          
          pl  = max(ql(i,j,QPRES ),small_pres)
          rel =     ql(i,j,QREINT)


          ! right state
          rr = max(qr(i,j,QRHO),small_dens)
          
          ! pick right velocities based on direction
          if(idir.eq.1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
          endif
          
          pr  = max(qr(i,j,QPRES),small_pres)
          rer =     qr(i,j,QREINT)

            
          ! common quantities
          taul = 1.d0/rl
          taur = 1.d0/rr
          
          ! lagrangian sound speeds
          clsql = gamcl(i,j)*pl*rl
          clsqr = gamcr(i,j)*pr*rr
          

          ! Note: in the original Colella & Glaz paper, they predicted
          ! gamma_e to the interfaces using a special (non-hyperbolic)
          ! evolution equation.  In Castro, we instead bring (rho e)
          ! to the edges, so we construct the necessary gamma_e here from
          ! what we have on the interfaces.
          gamel = pl/rel + 1
          gamer = pr/rer + 1
          
          ! these should consider a wider average of the cell-centered
          ! gammas
          gmin = min(gamel, gamer, 4.d0/3.d0)
          gmax = max(gamel, gamer, 5.d0/3.d0)
          
          game_bar = 0.5d0*(gamel + gamer)
          gamc_bar = 0.5d0*(gamcl(i,j) + gamcr(i,j))
          
          gdot = 2.d0*(1.d0 - game_bar/gamc_bar)*(game_bar - 1.0)
          
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(clsql)))
          wr = max(wsmall,sqrt(abs(clsqr)))
          
          ! make an initial guess for pstar -- this is a two-shock 
          ! approximation
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = max(pstar,small_pres)

          ! get the shock speeds -- this computes W_s from CG Eq. 34
          call wsqge(pl,taul,gamel,gdot,  &
                     gamstar,pstar,wlsq,clsql,gmin,gmax)

          call wsqge(pr,taur,gamer,gdot,  &
                     gamstar,pstar,wrsq,clsqr,gmin,gmax)

          pstnm1 = pstar

          wl = sqrt(wlsq)
          wr = sqrt(wrsq)

          ! R-H jump conditions give ustar across each wave -- these should
          ! be equal when we are done iterating
          ustarp = ul - (pstar-pl)/wl
          ustarm = ur + (pstar-pr)/wr

          ! revise our pstar guess
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = max(pstar,small_pres)

          ! sectant iteration
          converged = .false.
          iter = 1
          do while (iter < iter_max .and. .not. converged)
               
             call wsqge(pl,taul,gamel,gdot,  &
                        gamstar,pstar,wlsq,clsql,gmin,gmax)

             call wsqge(pr,taur,gamer,gdot,  &
                        gamstar,pstar,wrsq,clsqr,gmin,gmax)

             wl = 1.d0 / sqrt(wlsq)
             wr = 1.d0 / sqrt(wrsq)
             
             ustnm1 = ustarm
             ustnp1 = ustarp
             
             ustarm = ur-(pr-pstar)*wr
             ustarp = ul+(pl-pstar)*wl
             
             dpditer=abs(pstnm1-pstar)
             
             zp=abs(ustarp-ustnp1)
             if(zp-weakwv*cav(i,j) <= 0.d0)then
                zp = dpditer*wl
             endif
             
             zm=abs(ustarm-ustnm1)
             if(zm-weakwv*cav(i,j) <= 0.d0)then
                zm = dpditer*wr
             endif
             
             ! the new pstar is found via CG Eq. 18
             denom=dpditer/max(zp+zm,small*(cav(i,j)))
             pstnm1 = pstar
             pstar=pstar-denom*(ustarm-ustarp)
             pstar=max(pstar,small_pres)

             err = abs(pstar - pstnm1)
             if (err < tol*pstar) converged = .true.

             iter = iter + 1
             
          enddo

          if (.not. converged) then
             print *, iter
             print *, pstar, pstnm1
             print *, ustarm, ustarp
             print *, err, tol*pstar
             call bl_error("ERROR: non-convergence in the Riemann solver")
          endif
          
          
          ! we converged!  construct the single ustar for the region
          ! between the left and right waves
          ustar = 0.5d0* ( ustarp + ustarm )

          
          ! sample the solution -- here we look first at the direction
          ! that the contact is moving.  This tells us if we need to
          ! worry about the L/L* states or the R*/R states.  
          if (ustar .gt. 0.d0) then
             ro = rl
             uo = ul
             po = pl
             tauo = taul
             reo = rel
             gamco = gamcl(i,j)
             gameo = gamel
             
          else if (ustar .lt. 0.d0) then
             ro = rr
             uo = ur
             po = pr
             tauo = taur
             reo = rer
             gamco = gamcr(i,j)
             gameo = gamer
          else
             ro = 0.5d0*(rl+rr)
             uo = 0.5d0*(ul+ur)
             po = 0.5d0*(pl+pr)
             tauo = 1.d0/ro
             reo = 0.5d0*(rel+rer)
             gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))
             gameo = 0.5d0*(gamel + gamer)
          endif

          ro = max(small_dens,ro)
         
          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          clsq = (co*ro)**2

          call wsqge(po,tauo,gameo,gdot,   &
                     gamstar,pstar,wosq,clsq,gmin,gmax)

          sgnm = sign(1.d0,ustar)
          
          wo = sqrt(wosq)
          dpjmp = pstar - po

          ! is this max really necessary?
          rstar=max(1.d0-ro*dpjmp/wosq, (gameo-1.d0)/(gameo+1.d0))
          rstar=ro/rstar
          rstar = max(small_dens,rstar)

          entho = (reo/ro + po/ro)/co**2
          estar = reo + (pstar - po)*entho
          
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)
          
          
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          
          ushock = 0.5d0*(spin + spout)
          
          if (pstar-po .ge. 0.d0) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. 0.d0) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (1.d0 + (spout + spin)/scr)*0.5d0
          frac = max(0.d0,min(1.d0,frac))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar .gt. 0.d0) then
             vgdnv = vl
          else if (ustar .lt. 0.d0) then
             vgdnv = vr
          else
             vgdnv = 0.5d0*(vl+vr)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          rgdnv = frac*rstar + (1.d0 - frac)*ro          
          ugdnv(i,j) = frac*ustar + (1.d0 - frac)*uo
          pgdnv(i,j) = frac*pstar + (1.d0 - frac)*po
          gamgdnv =  frac*gamstar + (1.d0-frac)*gameo          

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. 0.d0) then
             rgdnv = ro
             ugdnv(i,j) = uo
             pgdnv(i,j) = po
             gamgdnv = gameo
          endif
          if (spin .ge. 0.d0) then
             rgdnv = rstar
             ugdnv(i,j) = ustar
             pgdnv(i,j) = pstar
             gamgdnv = gamstar
          endif

          
          pgdnv(i,j) = max(pgdnv(i,j),small_pres)

          ! Enforce that fluxes through a symmetry plane are hard zero.
          if (i.eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) &
               ugdnv(i,j) = 0.d0
          if (j.eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) &
               ugdnv(i,j) = 0.d0
          
          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgdnv*ugdnv(i,j)
          
          if(idir.eq.1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j) + pgdnv(i,j)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgdnv
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgdnv
             uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j) + pgdnv(i,j)
          endif

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = pgdnv(i,j)/(gamgdnv - 1.0d0) + &
               0.5d0*rgdnv*(ugdnv(i,j)**2 + vgdnv**2)

          uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv(i,j))
          uflx(i,j,UEINT) = ugdnv(i,j)*pgdnv(i,j)/(gamgdnv - 1.d0)


          ! advected quantities -- only the contact matters
          do iadv = 1, nadv
             n  = UFA + iadv - 1
             nq = QFA + iadv - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else
                qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
          
          ! species -- only the contact matters
          do ispec = 1, nspec
             n  = UFS + ispec - 1
             nq = QFS + ispec - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else
                qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
          
          ! auxillary quantities -- only the contact matters
          do iaux = 1, naux
             n  = UFX + iaux - 1
             nq = QFX + iaux - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else
                qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
          
       enddo
    enddo

  end subroutine riemanncg

  subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

    double precision p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax
    double precision smlp1,small,divide,temp

    data smlp1,small/.001d0,1.d-07/

    ! First predict a value of game across the shock

    ! CG Eq. 31
    gstar=(pstar-p)*gdot/(pstar+p) + gam
    gstar=max(gmin,min(gmax,gstar))

    ! Now use that predicted value of game with the R-H jump conditions
    ! to compute the wave speed.

    ! CG Eq. 34
    wsq = (0.5d0*(gstar-1.0d0)*(pstar+p)+pstar)
    temp = ((gstar-gam)/(gam-1.0d0))

    if (pstar-p.eq.0.0d0) then
       divide=small
    else
       divide=pstar-p
    endif
    
    temp=temp/divide
    wsq = wsq/(v - temp*p*v)
    if (abs(pstar/p-1.d0)-smlp1 .lt. 0.0d0 ) then
       wsq = csq
    endif
    wsq=max(wsq,(.5d0*(gam-1.d0)/gam)*csq)
    
    return
  end subroutine wsqge

! ::: 
! ::: ------------------------------------------------------------------
! ::: 


  subroutine riemannus(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                       gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                       uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                       pgdnv, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                       ugdnv, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                       idir, ilo1, ihi1, ilo2, ihi2)

    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo,Symmetry
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
                                   URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv, &
                                   small_dens, small_pres

    implicit none

    double precision, parameter:: small = 1.d-8

    integer qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer gd_l1, gd_l2, gd_h1, gd_h2
    integer uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer pgd_l1, pgd_l2, pgd_h1, pgd_h2
    integer ugd_l1, ugd_l2, ugd_h1, ugd_h2
    integer idir, ilo1, ihi1, ilo2, ihi2
    integer ilo,ihi,jlo,jhi

    double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision pgdnv(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
    double precision ugdnv(ugd_l1:ugd_h1,ugd_l2:ugd_h2)
    
    integer iadv, ispec, iaux, n, nq
    integer i, j

    double precision rgd, vgd, regd, ustar
    double precision rl, ul, vl, pl, rel
    double precision rr, ur, vr, pr, rer
    double precision wl, wr, rhoetot, scr
    double precision rstar, cstar, estar, pstar
    double precision ro, uo, po, reo, co, gamco, entho
    double precision sgnm, spin, spout, ushock, frac
    double precision wsmall, csmall,qavg

    !************************************************************
    !  set min/max based on normal direction
    if(idir.eq.1) then
       ilo = ilo1
       ihi = ihi1 + 1
       jlo = ilo2
       jhi = ihi2
    else
       ilo = ilo1
       ihi = ihi1
       jlo = ilo2
       jhi = ihi2+1
    endif

    !     Solve Riemann Problem
    !     NOTE: The calling routine will order velocity unknowns so that
    !     for the purposes of this routine, the normal component is always
    !     loaded in the QU slot.
    do j = jlo, jhi
       do i = ilo, ihi

          rl = ql(i,j,QRHO)

          !  pick left velocities based on direction
          if(idir.eq.1) then
             ul = ql(i,j,QU)
             vl = ql(i,j,QV)
          else
             ul = ql(i,j,QV)
             vl = ql(i,j,QU)
          endif
          
          pl = ql(i,j,QPRES)
          rel = ql(i,j,QREINT)
          
          rr = qr(i,j,QRHO)

          !  pick right velocities based on direction
          if(idir.eq.1) then
             ur = qr(i,j,QU)
             vr = qr(i,j,QV)
          else
             ur = qr(i,j,QV)
             vr = qr(i,j,QU)
          endif
          
          pr = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))
          
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = max(pstar,small_pres)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
          
          if (ustar .gt. 0.d0) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl(i,j)
          else if (ustar .lt. 0.d0) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
          else
             ro = 0.5d0*(rl+rr)
             uo = 0.5d0*(ul+ur)
             po = 0.5d0*(pl+pr)
             reo = 0.5d0*(rel+rer)
             gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))               
          endif
          ro = max(small_dens,ro)
          
          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          entho = (reo/ro + po/ro)/co**2
          rstar = ro + (pstar - po)/co**2
          rstar = max(small_dens,rstar)
          estar = reo + (pstar - po)*entho
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)
          
          sgnm = sign(1.d0,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = 0.5d0*(spin + spout)
          if (pstar-po .ge. 0.d0) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. 0.d0) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (1.d0 + (spout + spin)/scr)*0.5d0
          frac = max(0.d0,min(1.d0,frac))
          
          if (ustar .gt. 0.d0) then
             vgd = vl
          else if (ustar .lt. 0.d0) then
             vgd = vr
          else
             vgd = 0.5d0*(vl+vr)
          endif
          rgd = frac*rstar + (1.d0 - frac)*ro
          
          ugdnv(i,j) = frac*ustar + (1.d0 - frac)*uo
          pgdnv(i,j) = frac*pstar + (1.d0 - frac)*po

          regd = frac*estar + (1.d0 - frac)*reo
          if (spout .lt. 0.d0) then
             rgd = ro
             ugdnv(i,j) = uo
             pgdnv(i,j) = po
             regd = reo
          endif
          if (spin .ge. 0.d0) then
             rgd = rstar
             ugdnv(i,j) = ustar
             pgdnv(i,j) = pstar
             regd = estar
          endif

          ! Enforce that fluxes through a symmetry plane are hard zero.
          if (i.eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) ugdnv(i,j) = 0.d0
          if (j.eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) ugdnv(i,j) = 0.d0
          
          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgd*ugdnv(i,j)
          if(idir.eq.1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j)
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgd
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgd
             uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j)
          endif
          
          rhoetot = regd + 0.5d0*rgd*(ugdnv(i,j)**2 + vgd**2)
          uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv(i,j))
          uflx(i,j,UEINT) = ugdnv(i,j)*regd
          
          do iadv = 1, nadv
             n = UFA + iadv - 1
             nq = QFA + iadv - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else 
                qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
          
          do ispec = 1, nspec
             n  = UFS + ispec - 1
             nq = QFS + ispec - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else 
                qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
          
          do iaux = 1, naux
             n  = UFX + iaux - 1
             nq = QFX + iaux - 1
             if (ustar .gt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. 0.d0) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else 
                qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
          
       enddo
    enddo
  end subroutine riemannus

end module riemann_module
