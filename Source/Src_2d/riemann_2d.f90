module riemann_module

  implicit none

  private

  public cmpflx

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                    flx,flx_l1,flx_l2,flx_h1,flx_h2, &
                    pgd,pgd_l1,pgd_l2,pgd_h1,pgd_h2, &
                    ugd,ugd_l1,ugd_l2,ugd_h1,ugd_h2, &
                    gamc,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
                    idir,ilo,ihi,jlo,jhi,domlo,domhi)

    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QFS, QPRES, QREINT, &
                                   use_colglaz, ppm_temp_fix

    implicit none

    integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer flx_l1,flx_l2,flx_h1,flx_h2
    integer pgd_l1,pgd_l2,pgd_h1,pgd_h2
    integer ugd_l1,ugd_l2,ugd_h1,ugd_h2
    integer qd_l1,qd_l2,qd_h1,qd_h2
    integer idir,ilo,ihi,jlo,jhi
    integer domlo(2),domhi(2)
    
    double precision    qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision    qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision   flx(flx_l1:flx_h1,flx_l2:flx_h2,NVAR)
    double precision pgd(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
    double precision ugd(ugd_l1:ugd_h1,ugd_l2:ugd_h2)
    double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2)
    
    ! Local variables
    integer i, j
    
    double precision, allocatable :: smallc(:,:), cavg(:,:)
    double precision, allocatable :: gamcm(:,:), gamcp(:,:)
    
    integer :: imin, imax, jmin, jmax

    type (eos_t) :: eos_state

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )

    if(idir.eq.1) then
       do j = jlo, jhi
          do i = ilo, ihi+1
             smallc(i,j) = max( csml(i,j), csml(i-1,j) )
             cavg(i,j) = 0.5d0*( c(i,j) + c(i-1,j) )
             gamcm(i,j) = gamc(i-1,j)
             gamcp(i,j) = gamc(i,j)
          enddo
       enddo
    else
       do j = jlo, jhi+1
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j), csml(i,j-1) )
             cavg(i,j) = 0.5d0*( c(i,j) + c(i,j-1) )
             gamcm(i,j) = gamc(i,j-1)
             gamcp(i,j) = gamc(i,j)
          enddo
       enddo
    endif

    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent
       
       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are 
       ! thermodynamically consistent.
       
       if (idir == 1) then
          imin = ilo
          imax = ihi+1
          jmin = jlo
          jmax = jhi
       else
          imin = ilo
          imax = ihi
          jmin = jlo
          jmax = jhi+1
       endif
       
       do j = jmin, jmax
          do i = imin, imax
             
             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 10000.0d0   
             
             ! minus state
             eos_state%rho = qm(i,j,QRHO)
             eos_state%p   = qm(i,j,QPRES)
             eos_state%e   = qm(i,j,QREINT)/qm(i,j,QRHO)
             eos_state%xn  = qm(i,j,QFS:QFS-1+nspec)

             call eos(eos_input_re, eos_state, .false.)

             qm(i,j,QREINT) = qm(i,j,QRHO)*eos_state%e
             qm(i,j,QPRES) = eos_state%p
             gamcm(i,j) = eos_state%gam1


             ! plus state
             eos_state%rho = qp(i,j,QRHO)
             eos_state%p   = qp(i,j,QPRES)
             eos_state%e   = qp(i,j,QREINT)/qp(i,j,QRHO)
             eos_state%xn  = qp(i,j,QFS:QFS-1+nspec)

             call eos(eos_input_re, eos_state, .false.)

             qp(i,j,QREINT) = qp(i,j,QRHO)*eos_state%e
             qp(i,j,QPRES) = eos_state%p  
             gamcp(i,j) = eos_state%gam1

          enddo
       enddo
         
    endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved
    if (use_colglaz == 0) then
       call riemannus(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                      ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)
    else
       call riemanncg(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                      ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)
    endif
    
    deallocate(smallc,cavg,gamcm,gamcp)
    
  end subroutine cmpflx


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_h1,uflx_h2, &
                       pgdnv,pg_l1,pg_l2,pg_h1,pg_h2, &
                       ugdnv,ug_l1,ug_l2,ug_h1,ug_h2, &
                       idir,ilo1,ihi1,ilo2,ihi2,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module 
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall 
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, &
                                   QPRES, QREINT, QFA, QFS, &
                                   QFX, URHO, UMX, UMY, UEDEN, UEINT, &
                                   UFA, UFS, UFX, &
                                   nadv, small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_h1,uflx_h2
    integer :: ug_l1,ug_l2,ug_h1,ug_h2
    integer :: pg_l1,pg_l2,pg_h1,pg_h2
    integer :: idir,ilo1,ihi1,ilo2,ihi2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision ::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: ugdnv(ug_l1:ug_h1,ug_l2:ug_h2)
    double precision :: pgdnv(pg_l1:pg_h1,pg_l2:pg_h2)

    integer :: i,j,ilo,jlo,ihi,jhi
    integer :: n, nq
    integer :: iadv, ispec, iaux
    
    double precision :: rgdnv,vgdnv,ustar,gamgdnv
    double precision :: rl, ul, vl, pl, rel
    double precision :: rr, ur, vr, pr, rer
    double precision :: wl, wr, rhoetot
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

    double precision :: gcl, gcr
    double precision :: clsq, clsql, clsqr, wlsq, wosq, wrsq, wo
    double precision :: zm, zp
    double precision :: denom, dpditer, dpjmp
    double precision :: gamc_bar, game_bar
    double precision :: gamel, gamer, gameo, gamstar, gmin, gmax, gdot

    integer :: iter, iter_max
    double precision :: tol
    double precision :: err

    logical :: converged

    double precision :: pstnm1
    double precision :: taul, taur, tauo
    double precision :: ustarm, ustarp, ustnm1, ustnp1

    double precision, parameter :: weakwv = 1.d-3

    double precision, allocatable :: pstar_hist(:)

    type (eos_t) :: eos_state

    tol = cg_tol
    iter_max = cg_maxiter


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

    allocate (pstar_hist(iter_max))

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
          
          pl  = ql(i,j,QPRES )
          rel = ql(i,j,QREINT)
          gcl = gamcl(i,j)

          ! sometimes we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= 0.0d0 .or. pl <= small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
             eos_state%T = small_temp
             eos_state%rho = rl
             eos_state%xn(:) = ql(i,j,QFS:QFS-1+nspec)

             call eos(eos_input_rt, eos_state, .false.)

             rel = rl*eos_state%e
             pl = eos_state%p
             gcl = eos_state%gam1
          endif

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
          
          pr  = qr(i,j,QPRES)
          rer = qr(i,j,QREINT)
          gcr = gamcr(i,j)

          if (rer <= 0.0d0 .or. pr <= small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
             eos_state%T = small_temp
             eos_state%rho = rr
             eos_state%xn(:) = qr(i,j,QFS:QFS-1+nspec)

             call eos(eos_input_rt, eos_state, .false.)

             rer = rr*eos_state%e
             pr = eos_state%p
             gcr = eos_state%gam1
          endif
            
          ! common quantities
          taul = 1.d0/rl
          taur = 1.d0/rr
          
          ! lagrangian sound speeds
          clsql = gcl*pl*rl
          clsqr = gcr*pr*rr
          

          ! Note: in the original Colella & Glaz paper, they predicted
          ! gamma_e to the interfaces using a special (non-hyperbolic)
          ! evolution equation.  In Castro, we instead bring (rho e)
          ! to the edges, so we construct the necessary gamma_e here from
          ! what we have on the interfaces.
          gamel = pl/rel + 1.d0
          gamer = pr/rer + 1.d0
          
          ! these should consider a wider average of the cell-centered
          ! gammas
          gmin = min(gamel, gamer, 1.0, 4.d0/3.d0)
          gmax = max(gamel, gamer, 2.0, 5.d0/3.d0)
          
          game_bar = 0.5d0*(gamel + gamer)
          gamc_bar = 0.5d0*(gcl + gcr)
          
          gdot = 2.d0*(1.d0 - game_bar/gamc_bar)*(game_bar - 1.0d0)
          
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(clsql)))
          wr = max(wsmall,sqrt(abs(clsqr)))
          
          ! make an initial guess for pstar -- this is a two-shock 
          ! approximation
          !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
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
          !pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr)
          pstar = max(pstar,small_pres)

          ! sectant iteration
          converged = .false.
          iter = 1
          do while ((iter <= iter_max .and. .not. converged) .or. iter <= 2)
               
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

             pstar_hist(iter) = pstar

             iter = iter + 1
             
          enddo

          if (.not. converged) then
             print *, 'pstar history: '
             do iter = 1, iter_max
                print *, iter, pstar_hist(iter)
             enddo

             print *, ' '
             print *, 'left state  (r,u,p,re,gc): ', rl, ul, pl, rel, gcl
             print *, 'right state (r,u,p,re,gc): ', rr, ur, pr, rer, gcr
             call bl_error("ERROR: non-convergence in the Riemann solver")
          endif
          
          
          ! we converged!  construct the single ustar for the region
          ! between the left and right waves, using the updated wave speeds
          ustarm = ur-(pr-pstar)*wr  ! careful -- here wl, wr are 1/W
          ustarp = ul+(pl-pstar)*wl

          ustar = 0.5d0* ( ustarp + ustarm )

          
          ! sample the solution -- here we look first at the direction
          ! that the contact is moving.  This tells us if we need to
          ! worry about the L/L* states or the R*/R states.  
          if (ustar .gt. 0.d0) then
             ro = rl
             uo = ul
             po = pl
             tauo = taul
             !reo = rel
             gamco = gcl
             gameo = gamel
             
          else if (ustar .lt. 0.d0) then
             ro = rr
             uo = ur
             po = pr
             tauo = taur
             !reo = rer
             gamco = gcr
             gameo = gamer
          else
             ro = 0.5d0*(rl+rr)
             uo = 0.5d0*(ul+ur)
             po = 0.5d0*(pl+pr)
             tauo = 0.5d0*(taul+taur)
             !reo = 0.5d0*(rel+rer)
             gamco = 0.5d0*(gcl+gcr)
             gameo = 0.5d0*(gamel + gamer)
          endif

          ! use tau = 1/rho as the independent variable here
          ro = max(small_dens,1.0d0/tauo)
          tauo = 1.0d0/ro

          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          clsq = (co*ro)**2

          ! now that we know which state (left or right) we need to worry
          ! about, get the value of gamstar and wosq across the wave we
          ! are dealing with.
          call wsqge(po,tauo,gameo,gdot,   &
                     gamstar,pstar,wosq,clsq,gmin,gmax)

          sgnm = sign(1.d0,ustar)
          
          wo = sqrt(wosq)
          dpjmp = pstar - po

          ! is this max really necessary?
          !rstar=max(1.d0-ro*dpjmp/wosq, (gameo-1.d0)/(gameo+1.d0))
          rstar=1.d0-ro*dpjmp/wosq
          rstar=ro/rstar
          rstar = max(small_dens,rstar)

          !entho = (reo/ro + po/ro)/co**2
          !estar = reo + (pstar - po)*entho
          
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)
          
          
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          
          !ushock = 0.5d0*(spin + spout)
          ushock = wo/ro - sgnm*uo

          if (pstar-po .ge. 0.d0) then
             spin = ushock
             spout = ushock
          endif
          ! if (spout-spin .eq. 0.d0) then
          !    scr = small*cav(i,j)
          ! else
          !    scr = spout-spin
          ! endif
          ! frac = (1.d0 + (spout + spin)/scr)*0.5d0
          ! frac = max(0.d0,min(1.d0,frac))

          frac = 0.5d0*(1.0d0 + (spin + spout)/max(spout-spin,spin+spout, small*cav(i,j)))

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

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                 (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
             if (i.eq.domhi(1)+1 .and. &
                 (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                 (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
             if (j.eq.domhi(2)+1 .and. &
                 (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
          end if

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgdnv*ugdnv(i,j)
          
          ! note: here we do not include the pressure, since in 2-d,
          ! for some geometries, div{F} + grad{p} cannot be written
          ! in a flux difference form
          if(idir.eq.1) then
             uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j) 
             uflx(i,j,UMY) = uflx(i,j,URHO)*vgdnv
          else
             uflx(i,j,UMX) = uflx(i,j,URHO)*vgdnv
             uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j) 
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

    double precision, parameter :: smlp1 = 1.d-10
    double precision, parameter :: small = 1.d-7

    double precision :: alpha, beta

    ! First predict a value of game across the shock

    ! CG Eq. 31
    gstar=(pstar-p)*gdot/(pstar+p) + gam
    gstar=max(gmin,min(gmax,gstar))

    ! Now use that predicted value of game with the R-H jump conditions
    ! to compute the wave speed.

    ! CG Eq. 34
    ! wsq = (0.5d0*(gstar-1.0d0)*(pstar+p)+pstar)
    ! temp = ((gstar-gam)/(gam-1.0d0))

    ! if (pstar-p == 0.0d0) then
    !    divide=small
    ! else
    !    divide=pstar-p
    ! endif
    
    ! temp=temp/divide
    ! wsq = wsq/(v - temp*p*v)

    alpha = pstar-(gstar-1.0d0)*p/(gam-1.0d0)
    if (alpha == 0.0d0) alpha = smlp1*(pstar + p)

    beta = pstar + 0.5d0*(gstar-1.0)*(pstar+p)

    wsq = (pstar-p)*beta/(v*alpha)

    if (abs(pstar - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq=max(wsq,(0.5d0*(gam-1.d0)/gam)*csq)
    
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
                       idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall 
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
                                   URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv, &
                                   small_dens, small_pres

    implicit none

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: pgd_l1, pgd_l2, pgd_h1, pgd_h2
    integer :: ugd_l1, ugd_l2, ugd_h1, ugd_h2
    integer :: idir, ilo1, ihi1, ilo2, ihi2
    integer :: domlo(2),domhi(2)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision :: gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
    double precision :: pgdnv(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
    double precision :: ugdnv(ugd_l1:ugd_h1,ugd_l2:ugd_h2)
    
    integer :: ilo,ihi,jlo,jhi
    integer :: iadv, ispec, iaux, n, nq
    integer :: i, j

    double precision :: rgd, vgd, regd, ustar
    double precision :: rl, ul, vl, pl, rel
    double precision :: rr, ur, vr, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, gamco, entho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg

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

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                 (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
             if (i.eq.domhi(1)+1 .and. &
                 (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                 (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
             if (j.eq.domhi(2)+1 .and. &
                 (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = 0.d0
          end if
          
          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,URHO) = rgd*ugdnv(i,j)

          ! note: here we do not include the pressure, since in 2-d,
          ! for some geometries, div{F} + grad{p} cannot be written
          ! in a flux difference form
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
