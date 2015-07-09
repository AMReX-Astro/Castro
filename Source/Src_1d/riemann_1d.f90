module riemann_module

  implicit none

  private

  public cmpflx

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine cmpflx(lo,hi,domlo,domhi, &
                    qm,qp,qpd_l1,qpd_h1, &
                    flx,flx_l1,flx_h1, &
                    pgdnv,pg_l1,pg_h1, &
                    ugdnv,ug_l1,ug_h1, &
                    gamc,csml,c,qd_l1,qd_h1,ilo,ihi)
    
    use meth_params_module, only : QVAR, NVAR, &
                                   use_colglaz
    
    implicit none
    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer ilo,ihi
    integer qpd_l1,qpd_h1
    integer flx_l1, flx_h1
    integer  pg_l1, pg_h1
    integer  ug_l1, ug_h1
    integer  qd_l1,  qd_h1
    double precision    qm(qpd_l1:qpd_h1, QVAR)
    double precision    qp(qpd_l1:qpd_h1, QVAR)
    double precision   flx(flx_l1:flx_h1, NVAR)
    double precision pgdnv( pg_l1: pg_h1)
    double precision ugdnv( ug_l1: ug_h1)
    double precision  gamc( qd_l1: qd_h1)
    double precision     c( qd_l1: qd_h1)
    double precision  csml( qd_l1: qd_h1)
    
    ! Local variables
    integer i
    double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:)
    
    allocate ( smallc(ilo:ihi+1) )
    allocate ( cavg(ilo:ihi+1) )
    allocate ( gamcp(ilo:ihi+1) )
    allocate ( gamcm(ilo:ihi+1) )
    
    do i = ilo, ihi+1 
       smallc(i) = max( csml(i), csml(i-1) )
       cavg(i) = 0.5d0*( c(i) + c(i-1) )
       gamcm(i) = gamc(i-1)
       gamcp(i) = gamc(i)
    enddo
    
    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (use_colglaz == 0) then
       call riemannus(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                      gamcm, gamcp, flx, flx_l1, flx_h1, &
                      pgdnv, pg_l1, pg_h1, &
                      ugdnv, ug_l1, ug_h1, ilo, ihi, domlo, domhi )
    else
       call riemanncg(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                      gamcm, gamcp, flx, flx_l1, flx_h1, &
                      pgdnv, pg_l1, pg_h1, &
                      ugdnv, ug_l1, ug_h1, ilo, ihi)
    endif

    deallocate (smallc,cavg,gamcm,gamcp)
    
  end subroutine cmpflx


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemanncg(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                       gamcl,gamcr,uflx,uflx_l1,uflx_h1, &
                       pgdnv, pg_l1, pg_h1, &
                       ugdnv, ug_l1, ug_h1,ilo,ihi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)
        
    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, &
                                   QPRES, QREINT, QFA, QFS, &
                                   QFX, URHO, UMX, UEDEN, UEINT, &
                                   UFA, UFS, UFX, &
                                   nadv, small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol


    double precision, parameter :: small = 1.d-8

    integer :: qpd_l1, qpd_h1
    integer :: uflx_l1, uflx_h1
    integer :: ug_l1, ug_h1
    integer :: pg_l1, pg_h1


    double precision :: ql(qpd_l1:qpd_h1, QVAR+3)
    double precision :: qr(qpd_l1:qpd_h1, QVAR+3)
    double precision ::  gamcl(ilo:ihi+1)
    double precision ::  gamcr(ilo:ihi+1)
    double precision ::    cav(ilo:ihi+1)
    double precision :: smallc(ilo:ihi+1)
    double precision :: uflx(uflx_l1:uflx_h1, NVAR)
    double precision :: ugdnv(ug_l1:ug_h1)
    double precision :: pgdnv(pg_l1:pg_h1)

    integer :: ilo,ihi
    integer :: n, nq
    integer :: iadv, ispec, iaux

    double precision :: rgdnv,ustar,gamgdnv
    double precision :: rl, ul, pl, rel
    double precision :: rr, ur, pr, rer
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

    integer :: k
    
    tol = cg_tol
    iter_max = cg_maxiter

    allocate (pstar_hist(iter_max))

    do k = ilo, ihi+1

       ! left state
       rl  = max( ql(k,QRHO), small_dens)

       ul  = ql(k,QU)

       pl  = ql(k,QPRES)
       rel = ql(k,QREINT)
       gcl = gamcl(k)

       ! sometimes we come in here with negative energy or pressure
       ! note: reset both in either case, to remain thermo 
       ! consistent
       if (rel <= 0.d0 .or. pl <= small_pres) then
          print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
          eos_state % T   = small_temp
          eos_state % rho = rl
          eos_state % xn  = ql(k,QFS:QFS-1+nspec)
          eos_state % aux = ql(k,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          rel = rl*eos_state%e
          pl  = eos_state%p
          gcl = eos_state%gam1
       endif

       ! right state
       rr  = max( qr(k,QRHO), small_dens)

       ur  = qr(k,QU)

       pr  = qr(k,QPRES)
       rer = qr(k,QREINT)
       gcr = gamcr(k)

       if (rer <= 0.d0 .or. pr <= small_pres) then
          print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
          eos_state % T   = small_temp
          eos_state % rho = rr
          eos_state % xn  = qr(k,QFS:QFS-1+nspec)
          eos_state % aux = qr(k,QFX:QFX-1+naux)

          call eos(eos_input_rt, eos_state)

          rer = rr*eos_state%e
          pr  = eos_state%p
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

       csmall = smallc(k)
       wsmall = small_dens*csmall
       wl = max(wsmall,sqrt(abs(clsql)))
       wr = max(wsmall,sqrt(abs(clsqr)))

       ! make an initial guess for pstar -- this is a two-shock
       ! approximation
       !  pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
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

       ! secant iteration
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
          if(zp-weakwv*cav(k) <= 0.d0)then
             zp = dpditer*wl
          endif
          
          zm=abs(ustarm-ustnm1)
          if(zm-weakwv*cav(k) <= 0.d0)then
             zm = dpditer*wr
          endif
          
          ! the new pstar is found via CG Eq. 18
          denom=dpditer/max(zp+zm,small*(cav(k)))
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
       

       ! we converged! construct the single ustar for the region
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
       !rstar=max(1.-ro*dpjmp/wosq,(gameo-1.)/(gameo+1.))
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
       !   scr = small*cav(k)
       !  else
       !     scr = spout-spin
       !  endif
       !  frac = (1.d0 + (spout + spin)/scr)*0.5d0
       !  frac = max(0.d0,min(1.d0,frac))

       frac = 0.5d0*(1.0d0 + (spin + spout)/max(spout-spin,spin+spout, small*cav(k)))

       ! linearly interpolate between the star and normal state -- this covers the
       ! case where we are inside the rarefaction fan.
       rgdnv = frac*rstar + (1.d0 - frac)*ro
       ugdnv(k) = frac*ustar + (1.d0 - frac)*uo
       pgdnv(k) = frac*pstar + (1.d0 - frac)*po
       gamgdnv =  frac*gamstar + (1.d0-frac)*gameo

       ! now handle the cases where instead we are fully in the
       ! star or fully in the original (l/r) state
       if (spout .lt. 0.d0) then
          rgdnv = ro
          ugdnv(k) = uo
          pgdnv(k) = po
          gamgdnv = gameo
       endif
       if (spin .ge. 0.d0) then
          rgdnv = rstar
          ugdnv(k) = ustar
          pgdnv(k) = pstar
          gamgdnv = gamstar
       endif

       pgdnv(k) = max(pgdnv(k),small_pres)

       ! Compute fluxes, order as conserved state (not q)
       uflx(k,URHO) = rgdnv*ugdnv(k)

       ! note: here we do not include the pressure, since in 1-d,
       ! for some geometries, div{F} + grad{p} cannot be written
       ! in a flux difference form
       uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 

       ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
       rhoetot = pgdnv(k)/(gamgdnv - 1.0d0) + &
            0.5d0*rgdnv*ugdnv(k)**2

       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv(k))
       uflx(k,UEINT) = ugdnv(k)*pgdnv(k)/(gamgdnv - 1.d0)

       ! advected quantities -- only the contact matters
       do iadv = 1, nadv
          n = UFA + iadv - 1
          nq = QFA + iadv - 1
          if (ustar .gt. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else if (ustar .lt. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          else
             qavg = 0.5d0 * (ql(k,nq) + qr(k,nq))
             uflx(k,n) = uflx(k,URHO)*qavg
          endif
       enddo

       ! species -- only the contact matters
       do ispec = 1, nspec
          n = UFS + ispec - 1
          nq = QFS + ispec - 1
          if (ustar .gt. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else if (ustar .lt. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          else
             qavg = 0.5d0 * (ql(k,nq) + qr(k,nq))
             uflx(k,n) = uflx(k,URHO)*qavg
          endif
       enddo

       ! auxillary quantities -- only the contact matters
       do iaux = 1, naux
          n = UFX + iaux - 1
          nq = QFX + iaux - 1
          if (ustar .gt. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else if (ustar .lt. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          else
             qavg = 0.5d0 * (ql(k,nq) + qr(k,nq))
             uflx(k,n) = uflx(k,URHO)*qavg
          endif
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
    ! temp = ((gstar-gam)/(gam-1.))
    
    !if (pstar-p.eq.0.0d0) then
    !   divide=small
    !else
    !   divide=pstar-p
    !endif
    
    !temp=temp/divide
    !wsq = wsq/(v - temp*p*v)
    
    alpha = pstar-(gstar-1.0d0)*p/(gam-1.0d0)
    if (alpha == 0.0d0) alpha = smlp1*(pstar + p)
    
    beta = pstar + 0.5d0*(gstar-1.0d0)*(pstar+p)
    
    wsq = (pstar-p)*beta/(v*alpha)
    
    if (abs(pstar  - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq=max(wsq,(0.5d0*(gam-1.d0)/gam)*csq)
    
    return
  end subroutine wsqge
      
! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemannus(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                       gamcl,gamcr,uflx,uflx_l1,uflx_h1,&
                       pgdnv,pg_l1,pg_h1,ugdnv,ug_l1,ug_h1, &
                       ilo,ihi,domlo,domhi)
    
    use network, only : nspec, naux
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QPRES, QREINT, QFA, QFS, QFX, &
                                   URHO, UMX, UEDEN, UEINT, &
                                   UFA, UFS, UFX, nadv, small_dens, small_pres, &
                                   fix_mass_flux
    use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry

    implicit none

    double precision, parameter:: small = 1.d-8

    integer ilo,ihi
    integer domlo(1),domhi(1)
    integer  qpd_l1,  qpd_h1
    integer   pg_l1,   pg_h1
    integer   ug_l1,   ug_h1
    integer uflx_l1, uflx_h1
    double precision ql(qpd_l1:qpd_h1, QVAR)
    double precision qr(qpd_l1:qpd_h1, QVAR)
    double precision   cav(ilo:ihi+1), smallc(ilo:ihi+1)
    double precision gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
    double precision  uflx(uflx_l1:uflx_h1, NVAR)
    double precision pgdnv( pg_l1: pg_h1)
    double precision ugdnv( ug_l1: ug_h1)
    
    double precision rgdnv, regdnv, ustar
    double precision rl, ul, pl, rel
    double precision rr, ur, pr, rer
    double precision wl, wr, rhoetot, scr
    double precision rstar, cstar, estar, pstar
    double precision ro, uo, po, reo, co, gamco, entho
    double precision sgnm, spin, spout, ushock, frac
    
    double precision wsmall, csmall
    integer iadv, n, nq
    integer k,ispec, iaux
    logical :: fix_mass_flux_lo, fix_mass_flux_hi
    
    ! Solve Riemann Problem

    fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
    fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

    do k = ilo, ihi+1
       rl  = ql(k,QRHO)
       ul  = ql(k,QU)
       pl  = ql(k,QPRES)
       rel = ql(k,QREINT)
       rr  = qr(k,QRHO)
       ur  = qr(k,QU)
       pr  = qr(k,QPRES)
       rer = qr(k,QREINT)
       
       csmall = smallc(k)
       wsmall = small_dens*csmall
       wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
       wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))
       pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
       pstar = max(pstar,small_pres)
       ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
       
       if (ustar .gt. 0.d0) then
          ro = rl
          uo = ul
          po = pl
          reo = rel
          gamco = gamcl(k)
       else if (ustar .lt. 0.d0) then
          ro = rr
          uo = ur
          po = pr
          reo = rer
          gamco = gamcr(k)
       else
          ro = 0.5d0*(rl+rr)
          uo = 0.5d0*(ul+ur)
          po = 0.5d0*(pl+pr)
          reo = 0.5d0*(rel+rer)
          gamco = 0.5d0*(gamcl(k)+gamcr(k))
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
          scr = small*cav(k)
       else
          scr = spout-spin
       endif
       frac = (1.d0 + (spout + spin)/scr)*0.5d0
       frac = max(0.d0,min(1.d0,frac))
       
       rgdnv = frac*rstar + (1.d0 - frac)*ro
       ugdnv(k) = frac*ustar + (1.d0 - frac)*uo
       pgdnv(k) = frac*pstar + (1.d0 - frac)*po
       regdnv = frac*estar + (1.d0 - frac)*reo
       
       if (spout .lt. 0.d0) then
          rgdnv = ro
          ugdnv(k) = uo
          pgdnv(k) = po
          regdnv = reo
       endif
       if (spin .ge. 0.d0) then
          rgdnv = rstar
          ugdnv(k) = ustar
          pgdnv(k) = pstar
          regdnv = estar
       endif
       
       if (k.eq.0 .and. physbc_lo(1) .eq. Symmetry) ugdnv(k) = 0.d0
       
       if (fix_mass_flux_lo .and. k.eq.domlo(1) .and. ugdnv(k) .ge. 0.d0) then
          rgdnv    = ql(k,QRHO)
          ugdnv(k) = ql(k,QU)
          regdnv   = ql(k,QREINT)
       end if
       if (fix_mass_flux_hi .and. k.eq.domhi(1)+1 .and. ugdnv(k) .le. 0.d0) then
          rgdnv    = qr(k,QRHO)
          ugdnv(k) = qr(k,QU)
          regdnv   = qr(k,QREINT)
       end if

       ! Compute fluxes, order as conserved state (not q)
       uflx(k,URHO) = rgdnv*ugdnv(k)
       uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 
       
       rhoetot = regdnv + 0.5d0*rgdnv*ugdnv(k)**2 
       uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv(k))
       uflx(k,UEINT) = ugdnv(k)*regdnv
       
       do iadv = 1, nadv
          n = UFA + iadv - 1
          nq = QFA + iadv - 1
          if (ustar .ge. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif
       enddo
       
       do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          if (ustar .ge. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif
       enddo
       
       do iaux = 1, naux
          n  = UFX + iaux - 1
          nq = QFX + iaux - 1
          if (ustar .ge. 0.d0) then
             uflx(k,n) = uflx(k,URHO)*ql(k,nq)
          else
             uflx(k,n) = uflx(k,URHO)*qr(k,nq)
          endif
       enddo
       
    enddo
  end subroutine riemannus

end module riemann_module
