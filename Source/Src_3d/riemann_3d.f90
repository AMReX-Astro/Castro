module riemann_module

  implicit none

  private

  public cmpflx

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                    flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                    ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                    gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                    idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QFS, QFX, QPRES, QREINT, &
                                   use_colglaz, ppm_temp_fix
    use bl_constants_module

    implicit none

    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
    integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer idir,ilo,ihi,jlo,jhi
    integer i,j,kc,kflux,k3d
    integer domlo(3),domhi(3)

    double precision qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)
    double precision ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision gegdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision gamc(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision    c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision, allocatable :: smallc(:,:),cavg(:,:)
    double precision, allocatable :: gamcm(:,:),gamcp(:,:)

    type (eos_t) :: eos_state

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
      
    if(idir.eq.1) then
       do j = jlo, jhi
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i-1,j,k3d) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i-1,j,k3d) )
             gamcm(i,j) = gamc(i-1,j,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
          enddo
       enddo
    elseif(idir.eq.2) then
       do j = jlo, jhi
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i,j-1,k3d) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i,j-1,k3d) )
             gamcm(i,j) = gamc(i,j-1,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
          enddo
       enddo
    else
       do j = jlo, jhi
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i,j,k3d-1) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i,j,k3d-1) )
             gamcm(i,j) = gamc(i,j,k3d-1)
             gamcp(i,j) = gamc(i,j,k3d)
          enddo
       enddo
    endif
    
    if (ppm_temp_fix == 2) then
       ! recompute the thermodynamics on the interface to make it
       ! all consistent

       ! we want to take the edge states of rho, p, and X, and get
       ! new values for gamc and (rho e) on the edges that are 
       ! thermodynamically consistent.
       do j = jlo, jhi
          do i = ilo, ihi

             ! this is an initial guess for iterations, since we
             ! can't be certain that temp is on interfaces
             eos_state%T = 10000.0d0   
                              
             ! minus state
             eos_state%rho = qm(i,j,kc,QRHO)
             eos_state%p   = qm(i,j,kc,QPRES)
             eos_state%e   = qm(i,j,kc,QREINT)/qm(i,j,kc,QRHO)
             eos_state%xn  = qm(i,j,kc,QFS:QFS-1+nspec)
             eos_state%aux = qm(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_re, eos_state, .false.)

             qm(i,j,kc,QREINT) = qm(i,j,kc,QRHO)*eos_state%e
             qm(i,j,kc,QPRES) = eos_state%p
             gamcm(i,j) = eos_state%gam1


             ! plus state
             eos_state%rho = qp(i,j,kc,QRHO)
             eos_state%p   = qp(i,j,kc,QPRES)
             eos_state%e   = qp(i,j,kc,QREINT)/qp(i,j,kc,QRHO)
             eos_state%xn  = qp(i,j,kc,QFS:QFS-1+nspec)
             eos_state%aux = qp(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_re, eos_state, .false.)

             qp(i,j,kc,QREINT) = qp(i,j,kc,QRHO)*eos_state%e
             qp(i,j,kc,QPRES) = eos_state%p
             gamcp(i,j) = eos_state%gam1

          enddo
       enddo
         
    endif

    ! Solve Riemann problem
    if (use_colglaz == 1) then
       call riemanncg(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                      gamcm,gamcp,cavg,smallc,ilo-1,jlo-1,ihi+1,jhi+1, &
                      flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                      ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                      idir,ilo,ihi,jlo,jhi,kc,kflux,domlo,domhi)

    else
       call riemannus(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                      gamcm,gamcp,cavg,smallc,ilo-1,jlo-1,ihi+1,jhi+1, &
                      flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                      ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                      idir,ilo,ihi,jlo,jhi,kc,kflux,domlo,domhi)
    endif

    deallocate(smallc,cavg,gamcm,gamcp)

  end subroutine cmpflx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                       ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                       idir,ilo,ihi,jlo,jhi,kc,kflux,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QESGS, QFA, QFS, &
                                   QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UESGS, UFA, UFS, UFX, &
                                   nadv, small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol
    use bl_constants_module

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3
    integer :: pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision ::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,uflx_l3:uflx_h3,NVAR)
    double precision :: ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision :: pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision :: gegdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)

    integer :: i,j,kc,kflux
    integer :: n, nq
    integer :: iadv, ispec, iaux
    
    double precision :: rgdnv,v1gdnv,v2gdnv,regdnv,ustar,gamgdnv
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, pstar
    double precision :: ro, uo, po, co, gamco
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg
    double precision :: rho_K_contrib

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

    allocate (pstar_hist(iter_max))

    !$OMP PARALLEL DO PRIVATE(i,j) &
    !$OMP PRIVATE(rl,ul,v1l,v2l,pl,rel,rr,ur,v1r,v2r,pr,rer,gcl,gcr) &
    !$OMP PRIVATE(taul,taur,clsql,clsqr,gamel,gamer,gmin,gmax) &
    !$OMP PRIVATE(game_bar,gamc_bar,gdot,csmall,wsmall,wl,wr) &
    !$OMP PRIVATE(pstar,gamstar,wlsq,wrsq,pstnm1) &
    !$OMP PRIVATE(ustarp,ustarm,converged,iter,ustnm1,ustnp1) &
    !$OMP PRIVATE(dpditer,zp,zm,denom,err,ustar) &
    !$OMP PRIVATE(ro,uo,po,tauo,gamco,gameo,co,clsq,wosq,sgnm,wo,dpjmp) &
    !$OMP PRIVATE(rstar,cstar,spout,spin,ushock,scr,frac) &
    !$OMP PRIVATE(v1gdnv,v2gdnv,rgdnv,gamgdnv) &
    !$OMP PRIVATE(rhoetot,n,nq,qavg,rho_K_contrib,iadv,ispec,iaux) &
    !$OMP PRIVATE(pstar_hist) &
    !$OMP PRIVATE(eos_state)
    do j = jlo, jhi
       do i = ilo, ihi

          ! left state
          rl = max(ql(i,j,kc,QRHO),small_dens)
          
          ! pick left velocities based on direction
          if(idir.eq.1) then
             ul  = ql(i,j,kc,QU)
             v1l = ql(i,j,kc,QV)
             v2l = ql(i,j,kc,QW)
          elseif(idir.eq.2) then
             ul  = ql(i,j,kc,QV)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QW)
          else
             ul  = ql(i,j,kc,QW)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QV)
          endif
          
          pl  = ql(i,j,kc,QPRES)
          rel = ql(i,j,kc,QREINT)
          gcl = gamcl(i,j)

          ! sometime we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= ZERO .or. pl < small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres

             eos_state%T   = small_temp
             eos_state%rho = rl
             eos_state%xn  = ql(i,j,kc,QFS:QFS-1+nspec)
             eos_state%aux = ql(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state, .false.)

             rel = rl*eos_state%e
             pl  = eos_state%p
             gcl = eos_state%gam1
          endif

          ! right state
          rr = max(qr(i,j,kc,QRHO),small_dens)
          
          ! pick right velocities based on direction
          if(idir.eq.1) then
             ur  = qr(i,j,kc,QU)
             v1r = qr(i,j,kc,QV)
             v2r = qr(i,j,kc,QW)
          elseif(idir.eq.2) then
             ur  = qr(i,j,kc,QV)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QW)
          else
             ur  = qr(i,j,kc,QW)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QV)
          endif
          
          pr  = qr(i,j,kc,QPRES)
          rer = qr(i,j,kc,QREINT)
          gcr = gamcr(i,j)

          if (rer <= ZERO .or. pr < small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres

             eos_state % T   = small_temp
             eos_state % rho = rr
             eos_state % xn  = qr(i,j,kc,QFS:QFS-1+nspec)
             eos_state % aux = qr(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state, .false.)

             rer = rr*eos_state%e
             pr  = eos_state%p
             gcr = eos_state%gam1
          endif

            
          ! common quantities
          taul = ONE/rl
          taur = ONE/rr
          
          ! lagrangian sound speeds
          clsql = gcl*pl*rl
          clsqr = gcr*pr*rr
          

          ! Note: in the original Colella & Glaz paper, they predicted
          ! gamma_e to the interfaces using a special (non-hyperbolic)
          ! evolution equation.  In Castro, we instead bring (rho e)
          ! to the edges, so we construct the necessary gamma_e here from
          ! what we have on the interfaces.
          gamel = pl/rel + ONE
          gamer = pr/rer + ONE
          
          ! these should consider a wider average of the cell-centered
          ! gammas
          gmin = min(gamel, gamer, ONE, FOUR3RD)
          gmax = max(gamel, gamer, TWO, FIVE3RD)
          
          game_bar = HALF*(gamel + gamer)
          gamc_bar = HALF*(gcl + gcr)
          
          gdot = TWO*(ONE - game_bar/gamc_bar)*(game_bar - ONE)
          
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

          ! R-H jump conditions give ustar across each wave -- these
          ! should be equal when we are done iterating.  Our notation
          ! here is a little funny, comparing to CG, ustarp = u*_L and
          ! ustarm = u*_R.
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

             ! NOTE: these are really the inverses of the wave speeds!
             wl = ONE / sqrt(wlsq)
             wr = ONE / sqrt(wrsq)
             
             ustnm1 = ustarm
             ustnp1 = ustarp
             
             ustarm = ur-(pr-pstar)*wr
             ustarp = ul+(pl-pstar)*wl
             
             dpditer=abs(pstnm1-pstar)
             
             ! Here we are going to do the Secant iteration version in
             ! CG.  Note that what we call zp and zm here are not
             ! actually the Z_p = |dp*/du*_p| defined in CG, by rather
             ! simply |du*_p| (or something that looks like dp/Z!).
             zp=abs(ustarp-ustnp1)
             if(zp-weakwv*cav(i,j) <= ZERO)then
                zp = dpditer*wl
             endif
             
             zm=abs(ustarm-ustnm1)
             if(zm-weakwv*cav(i,j) <= ZERO)then
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

          ustar = HALF* ( ustarp + ustarm )

          
          ! sample the solution -- here we look first at the direction
          ! that the contact is moving.  This tells us if we need to
          ! worry about the L/L* states or the R*/R states.  
          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             tauo = taul
             gamco = gcl
             gameo = gamel
             
          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             tauo = taur
             gamco = gcr
             gameo = gamer
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             tauo = HALF*(taul+taur)
             gamco = HALF*(gcl+gcr)
             gameo = HALF*(gamel + gamer)
          endif

          ! use tau = 1/rho as the independent variable here
          ro = max(small_dens,ONE/tauo)
          tauo = ONE/ro
         
          co = sqrt(abs(gamco*po/ro))
          co = max(csmall,co)
          clsq = (co*ro)**2

          ! now that we know which state (left or right) we need to worry 
          ! about, get the value of gamstar and wosq across the wave we
          ! are dealing with.
          call wsqge(po,tauo,gameo,gdot,   &
                     gamstar,pstar,wosq,clsq,gmin,gmax)

          sgnm = sign(ONE,ustar)
          
          wo = sqrt(wosq)
          dpjmp = pstar - po

          ! is this max really necessary?
          !rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
          rstar=ONE-ro*dpjmp/wosq
          rstar=ro/rstar
          rstar = max(small_dens,rstar)

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)
          
          
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          
          !ushock = HALF*(spin + spout)
          ushock = wo/ro - sgnm*uo
          
          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          !if (spout-spin .eq. ZERO) then
          !   scr = small*cav(i,j)
          !else
          !   scr = spout-spin
          !endif
          !frac = (ONE + (spout + spin)/scr)*HALF
          !frac = max(ZERO,min(ONE,frac))

          frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav(i,j)))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar .gt. ZERO) then
             v1gdnv = v1l
             v2gdnv = v2l
          else if (ustar .lt. ZERO) then
             v1gdnv = v1r
             v2gdnv = v2r
          else
             v1gdnv = HALF*(v1l+v1r)
             v2gdnv = HALF*(v2l+v2r)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          rgdnv = frac*rstar + (ONE - frac)*ro          
          ugdnv(i,j,kc) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j,kc) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo          

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j,kc) = uo
             pgdnv(i,j,kc) = po
             gamgdnv = gameo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j,kc) = ustar
             pgdnv(i,j,kc) = pstar
             gamgdnv = gamstar
          endif
          
          pgdnv(i,j,kc) = max(pgdnv(i,j,kc),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                 (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                 (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                 (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                 (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          if (idir .eq. 3) then
             if (kflux.eq.domlo(3) .and. &
                 (physbc_lo(3) .eq. Symmetry .or.  physbc_lo(3) .eq. SlipWall .or. &
                  physbc_lo(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (kflux.eq.domhi(3)+1 .and. &
                 (physbc_hi(3) .eq. Symmetry .or.  physbc_hi(3) .eq. SlipWall .or. &
                  physbc_hi(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)
          
          if(idir.eq.1) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
          elseif(idir.eq.2) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
          else
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v2gdnv
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
          endif

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = pgdnv(i,j,kc)/(gamgdnv - ONE) + &
               HALF*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

          uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv(i,j,kc))
          uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*pgdnv(i,j,kc)/(gamgdnv - ONE)


          ! Treat K as a passively advected quantity but allow it to
          ! affect fluxes of (rho E) and momenta.
          if (UESGS .gt. -1) then
             n  = UESGS
             nq = QESGS
             if (ustar .gt. ZERO) then
                qavg = ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                qavg = qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
             endif
             
             uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             
             rho_K_contrib =  TWO3RD * rgdnv * qavg
             
             if(idir.eq.1) then
                uflx(i,j,kflux,UMX) = uflx(i,j,kflux,UMX) + rho_K_contrib
             elseif(idir.eq.2) then
                uflx(i,j,kflux,UMY) = uflx(i,j,kflux,UMY) + rho_K_contrib
             elseif(idir.eq.3) then
                uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,UMZ) + rho_K_contrib
             endif
             
             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + ugdnv(i,j,kc) * rho_K_contrib
          end if

          ! advected quantities -- only the contact matters
          do iadv = 1, nadv
             n  = UFA + iadv - 1
             nq = QFA + iadv - 1
             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
          ! species -- only the contact matters
          do ispec = 1, nspec
             n  = UFS + ispec - 1
             nq = QFS + ispec - 1
             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
          ! auxillary quantities -- only the contact matters
          do iaux = 1, naux
             n  = UFX + iaux - 1
             nq = QFX + iaux - 1
             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
       enddo
    enddo

  end subroutine riemanncg

  subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

    use bl_constants_module
    
    implicit none

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
    ! wsq = (HALF*(gstar-ONE)*(pstar+p)+pstar)
    ! temp = ((gstar-gam)/(gam-ONE))

    ! if (pstar-p == ZERO) then
    !    divide=small
    ! else
    !    divide=pstar-p
    ! endif
    
    ! temp=temp/divide
    ! wsq = wsq/(v - temp*p*v)

    alpha = pstar-(gstar-ONE)*p/(gam-ONE)
    if (alpha == ZERO) alpha = smlp1*(pstar + p)

    beta = pstar + HALF*(gstar-ONE)*(pstar+p)

    wsq = (pstar-p)*beta/(v*alpha)

    if (abs(pstar - p) < smlp1*(pstar + p)) then
       wsq = csq
    endif
    wsq=max(wsq,(HALF*(gam-ONE)/gam)*csq)
    
    return
  end subroutine wsqge

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemannus(ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                       ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                       idir,ilo,ihi,jlo,jhi,kc,kflux,domlo,domhi)

    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QESGS, QFA, QFS, &
                                     QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, UFA, UFS, UFX, &
                                     nadv, small_dens, small_pres
    use bl_constants_module

    implicit none
    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3
    integer :: pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)
    
    double precision :: ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision :: qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision ::  gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::  gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision ::    cav(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: smallc(gd_l1:gd_h1,gd_l2:gd_h2)
    double precision :: uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,uflx_l3:uflx_h3,NVAR)
    double precision :: ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision :: pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision :: gegdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    
    integer :: i,j,kc,kflux
    integer :: n, nq
    integer :: iadv, ispec, iaux
    
    double precision :: rgdnv,v1gdnv,v2gdnv,regdnv,ustar
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar
    double precision :: ro, uo, po, reo, co, gamco, entho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg
    double precision :: rho_K_contrib
    
    !$OMP PARALLEL DO PRIVATE(i,j,rl,ul,v1l,v2l,pl,rel,rr,ur,v1r,v2r,pr,rer,csmall,wsmall,wl,wr,pstar,ustar,ro,uo) &
    !$OMP PRIVATE(po,reo,gamco,co,entho,rstar,estar,cstar,sgnm,spout,spin,ushock,scr,frac,v1gdnv,v2gdnv,rgdnv,regdnv) &
    !$OMP PRIVATE(rhoetot,iadv,n,nq,qavg,ispec,iaux,rho_K_contrib)
    do j = jlo, jhi
       do i = ilo, ihi

          rl = max(ql(i,j,kc,QRHO),small_dens)

          ! pick left velocities based on direction
          if(idir.eq.1) then
             ul  = ql(i,j,kc,QU)
             v1l = ql(i,j,kc,QV)
             v2l = ql(i,j,kc,QW)
          elseif(idir.eq.2) then
             ul  = ql(i,j,kc,QV)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QW)
          else
             ul  = ql(i,j,kc,QW)
             v1l = ql(i,j,kc,QU)
             v2l = ql(i,j,kc,QV)
          endif
          
          pl  = max(ql(i,j,kc,QPRES ),small_pres)
          rel =     ql(i,j,kc,QREINT)
          
          rr = max(qr(i,j,kc,QRHO),small_dens)
          
          ! pick right velocities based on direction
          if(idir.eq.1) then
             ur  = qr(i,j,kc,QU)
             v1r = qr(i,j,kc,QV)
             v2r = qr(i,j,kc,QW)
          elseif(idir.eq.2) then
             ur  = qr(i,j,kc,QV)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QW)
          else
             ur  = qr(i,j,kc,QW)
             v1r = qr(i,j,kc,QU)
             v2r = qr(i,j,kc,QV)
          endif
          
          pr  = max(qr(i,j,kc,QPRES),small_pres)
          rer =     qr(i,j,kc,QREINT)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))
          
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)
          pstar = max(pstar,small_pres)
          
          if (ustar .gt. ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl(i,j)
          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
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
          
          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)
          if (pstar-po .ge. ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin .eq. ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))
          
          if (ustar .gt. ZERO) then
             v1gdnv = v1l
             v2gdnv = v2l
          else if (ustar .lt. ZERO) then
             v1gdnv = v1r
             v2gdnv = v2r
          else
             v1gdnv = HALF*(v1l+v1r)
             v2gdnv = HALF*(v2l+v2r)
          endif
          rgdnv = frac*rstar + (ONE - frac)*ro
          
          ugdnv(i,j,kc) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j,kc) = frac*pstar + (ONE - frac)*po
          
          regdnv = frac*estar + (ONE - frac)*reo
          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j,kc) = uo
             pgdnv(i,j,kc) = po
             regdnv = reo
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j,kc) = ustar
             pgdnv(i,j,kc) = pstar
             regdnv = estar
          endif
          
          pgdnv(i,j,kc) = max(pgdnv(i,j,kc),small_pres)
            
          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                  (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                  (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                  (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                  (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          if (idir .eq. 3) then
             if (kflux.eq.domlo(3) .and. &
                  (physbc_lo(3) .eq. Symmetry .or.  physbc_lo(3) .eq. SlipWall .or. &
                  physbc_lo(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
             if (kflux.eq.domhi(3)+1 .and. &
                  (physbc_hi(3) .eq. Symmetry .or.  physbc_hi(3) .eq. SlipWall .or. &
                  physbc_hi(3) .eq. NoSlipWall) ) &
                  ugdnv(i,j,kc) = ZERO
          end if
          
          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)
          
          if(idir.eq.1) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
          elseif(idir.eq.2) then
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*v2gdnv
          else
             uflx(i,j,kflux,UMX) = uflx(i,j,kflux,URHO)*v1gdnv
             uflx(i,j,kflux,UMY) = uflx(i,j,kflux,URHO)*v2gdnv
             uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
          endif
          
          rhoetot = regdnv + HALF*rgdnv*(ugdnv(i,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

          uflx(i,j,kflux,UEDEN) = ugdnv(i,j,kc)*(rhoetot + pgdnv(i,j,kc))
          uflx(i,j,kflux,UEINT) = ugdnv(i,j,kc)*regdnv
          
          ! Treat K as a passively advected quantity but allow it to affect fluxes of (rho E) and momenta.
          if (UESGS .gt. -1) then
             n  = UESGS
             nq = QESGS
             if (ustar .gt. ZERO) then
                qavg = ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                qavg = qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
             endif
             
             uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             
             rho_K_contrib =  TWO3RD * rgdnv * qavg
             
             if(idir.eq.1) then
                uflx(i,j,kflux,UMX) = uflx(i,j,kflux,UMX) + rho_K_contrib
             elseif(idir.eq.2) then
                uflx(i,j,kflux,UMY) = uflx(i,j,kflux,UMY) + rho_K_contrib
             elseif(idir.eq.3) then
                uflx(i,j,kflux,UMZ) = uflx(i,j,kflux,UMZ) + rho_K_contrib
             endif
             
             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + ugdnv(i,j,kc) * rho_K_contrib
          end if
          
          do iadv = 1, nadv
             n  = UFA + iadv - 1
             nq = QFA + iadv - 1
             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
          do ispec = 1, nspec
             n  = UFS + ispec - 1
             nq = QFS + ispec - 1
             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
          do iaux = 1, naux
             n  = UFX + iaux - 1
             nq = QFX + iaux - 1
             if (ustar .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo
          
       enddo
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine riemannus

end module riemann_module
