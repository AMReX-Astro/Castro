module riemann_module

  use bl_constants_module

  implicit none

  private

  public cmpflx, shock

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                    flx,flx_l1,flx_l2,flx_h1,flx_h2, &
                    pgd,pgd_l1,pgd_l2,pgd_h1,pgd_h2, &
                    ugd,ugd_l1,ugd_l2,ugd_h1,ugd_h2, &
                    gegd,ggd_l1,ggd_l2,ggd_h1,ggd_h2, &
                    gamc,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
                    shk,s_l1,s_l2,s_h1,s_h2, &
                    idir,ilo,ihi,jlo,jhi,domlo,domhi)

    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QFS, QFX, QPRES, QREINT, &
                                   use_colglaz, ppm_temp_fix, hybrid_riemann, &
                                   small_temp, allow_negative_energy


    implicit none

    integer, intent(in) :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer, intent(in) :: flx_l1,flx_l2,flx_h1,flx_h2
    integer, intent(in) :: pgd_l1,pgd_l2,pgd_h1,pgd_h2
    integer, intent(in) :: ugd_l1,ugd_l2,ugd_h1,ugd_h2
    integer, intent(in) :: ggd_l1,ggd_l2,ggd_h1,ggd_h2
    integer, intent(in) :: qd_l1,qd_l2,qd_h1,qd_h2
    integer, intent(in) :: s_l1,s_l2,s_h1,s_h2
    integer, intent(in) :: idir,ilo,ihi,jlo,jhi
    integer, intent(in) :: domlo(2),domhi(2)
    
    double precision, intent(inout) ::  qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision, intent(inout) ::  qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
    double precision, intent(inout) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,NVAR)
    double precision, intent(inout) :: pgd(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
    double precision, intent(inout) :: ugd(ugd_l1:ugd_h1,ugd_l2:ugd_h2)
    double precision, intent(inout) ::gegd(ggd_l1:ggd_h1,ggd_l2:ggd_h2)

    double precision, intent(in) :: gamc(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) ::    c(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) :: csml(qd_l1:qd_h1,qd_l2:qd_h2)
    double precision, intent(in) ::  shk( s_l1: s_h1, s_l2: s_h2)
    
    ! Local variables
    integer i, j
    
    double precision, allocatable :: smallc(:,:), cavg(:,:)
    double precision, allocatable :: gamcm(:,:), gamcp(:,:)
    
    integer :: imin, imax, jmin, jmax
    integer :: is_shock
    double precision :: cl, cr
    type (eos_t) :: eos_state

    allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
    allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )

    if (idir == 1) then
       do j = jlo, jhi
          do i = ilo, ihi+1
             smallc(i,j) = max( csml(i,j), csml(i-1,j) )
             cavg(i,j) = HALF*( c(i,j) + c(i-1,j) )
             gamcm(i,j) = gamc(i-1,j)
             gamcp(i,j) = gamc(i,j)
          enddo
       enddo

    else
       do j = jlo, jhi+1
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j), csml(i,j-1) )
             cavg(i,j) = HALF*( c(i,j) + c(i,j-1) )
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
             eos_state % rho = qm(i,j,QRHO)
             eos_state % p   = qm(i,j,QPRES)
             eos_state % e   = qm(i,j,QREINT)/qm(i,j,QRHO)
             eos_state % xn  = qm(i,j,QFS:QFS-1+nspec)
             eos_state % aux = qm(i,j,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state)                
             else
                call eos(eos_input_re, eos_state)
             endif

             qm(i,j,QREINT) = qm(i,j,QRHO)*eos_state%e
             qm(i,j,QPRES) = eos_state%p
             gamcm(i,j) = eos_state%gam1


             ! plus state
             eos_state % rho = qp(i,j,QRHO)
             eos_state % p   = qp(i,j,QPRES)
             eos_state % e   = qp(i,j,QREINT)/qp(i,j,QRHO)
             eos_state % xn  = qp(i,j,QFS:QFS-1+nspec)
             eos_state % aux = qp(i,j,QFX:QFX-1+naux)

             ! Protect against negative energies

             if (allow_negative_energy .eq. 0 .and. eos_state % e < ZERO) then
                eos_state % T = small_temp
                call eos(eos_input_rt, eos_state)                
             else
                call eos(eos_input_re, eos_state)
             endif

             qp(i,j,QREINT) = qp(i,j,QRHO)*eos_state%e
             qp(i,j,QPRES) = eos_state%p  
             gamcp(i,j) = eos_state%gam1

          enddo
       enddo
         
    endif

    ! Solve Riemann problem (godunov state passed back, but only (u,p) saved)
    if (use_colglaz == 0) then
       call riemannus(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                      ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                      gegd, ggd_l1, ggd_l2, ggd_h1, ggd_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)
    else
       call riemanncg(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                      gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                      pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                      ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                      gegd, ggd_l1, ggd_l2, ggd_h1, ggd_h2, &
                      idir, ilo, ihi, jlo, jhi, domlo, domhi)
    endif

    if (hybrid_riemann == 1) then
       ! correct the fluxes using an HLL scheme if we are in a shock
       ! and doing the hybrid approach
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
         
             if (idir == 1) then
                is_shock = shk(i-1,j) + shk(i,j)
             else
                is_shock = shk(i,j-1) + shk(i,j)
             endif

             if (is_shock >= 1) then

                if (idir == 1) then
                   cl = c(i-1,j)
                   cr = c(i,j)
                else
                   cl = c(i,j-1)
                   cr = c(i,j)
                endif

                call HLL(qm(i,j,:), qp(i,j,:), cl, cr, &
                         idir, flx(i,j,:))
                
             endif
    
          enddo
       enddo
    
    endif

    deallocate(smallc,cavg,gamcm,gamcp)
    
  end subroutine cmpflx


  subroutine shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                   shk,s_l1,s_l2,s_h1,s_h2, &
                   ilo1,ilo2,ihi1,ihi2,dx,dy)  

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV, QPRES, QVAR

    integer, intent(in) :: qd_l1, qd_l2, qd_h1, qd_h2
    integer, intent(in) :: s_l1, s_l2, s_h1, s_h2
    integer, intent(in) :: ilo1, ilo2, ihi1, ihi2
    double precision, intent(in) :: dx, dy
    double precision, intent(in) :: q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
    double precision, intent(inout) :: shk(s_l1:s_h1,s_l2:s_h2)

    integer :: i, j

    double precision :: divU
    double precision :: px_pre, px_post, py_pre, py_post
    double precision :: e_x, e_y, d
    double precision :: p_pre, p_post, pjump

    double precision, parameter :: small = 1.d-10
    double precision, parameter :: eps = 0.33d0

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    if (coord_type /= 0) then
       call bl_error("ERROR: invalid geometry in shock()")
    endif

    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          ! construct div{U}
          divU = HALF*(q(i+1,j,QU) - q(i-1,j,QU))/dx + &
                 HALF*(q(i,j+1,QV) - q(i,j-1,QV))/dy
          
          ! find the pre- and post-shock pressures in each direction
          if (q(i+1,j,QPRES) - q(i-1,j,QPRES) < ZERO) then
             px_pre  = q(i+1,j,QPRES)
             px_post = q(i-1,j,QPRES)
          else
             px_pre  = q(i-1,j,QPRES)
             px_post = q(i+1,j,QPRES)
          endif

          if (q(i,j+1,QPRES) - q(i,j-1,QPRES) < ZERO) then
             py_pre  = q(i,j+1,QPRES)
             py_post = q(i,j-1,QPRES)
          else
             py_pre  = q(i,j-1,QPRES)
             py_post = q(i,j+1,QPRES)
          endif

          ! use compression to create unit vectors for the shock direction
          e_x = (q(i+1,j,QU) - q(i-1,j,QU))**2
          e_y = (q(i,j+1,QV) - q(i,j-1,QV))**2
          d = ONE/(e_x + e_y + small)

          e_x = e_x*d
          e_y = e_y*d 

          ! project the pressures onto the shock direction 
          p_pre  = e_x*px_pre + e_y*py_pre
          p_post = e_x*px_post + e_y*py_post

          ! test for compression + pressure jump to flag a shock
          pjump = eps - (p_post - p_pre)/p_pre

          if (pjump < ZERO .and. divU < ZERO) then
             shk(i,j) = ONE
          else
             shk(i,j) = ZERO
          endif

       enddo
    enddo

  end subroutine shock


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_h1,uflx_h2, &
                       pgdnv,pg_l1,pg_l2,pg_h1,pg_h2, &
                       ugdnv,ug_l1,ug_l2,ug_h1,ug_h2, &
                       gegdnv,gg_l1,gg_l2,gg_h1,gg_h2, &
                       idir,ilo1,ihi1,ilo2,ihi2,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module 
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall 
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, &
                                   QPRES, QREINT, QFS, &
                                   QFX, URHO, UMX, UMY, UEDEN, UEINT, &
                                   small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol, &
                                   npassive, upass_map, qpass_map
    

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1,qpd_l2,qpd_h1,qpd_h2
    integer :: gd_l1,gd_l2,gd_h1,gd_h2
    integer :: uflx_l1,uflx_l2,uflx_h1,uflx_h2
    integer :: ug_l1,ug_l2,ug_h1,ug_h2
    integer :: pg_l1,pg_l2,pg_h1,pg_h2
    integer :: gg_l1,gg_l2,gg_h1,gg_h2
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
    double precision :: gegdnv(gg_l1:gg_h1,gg_l2:gg_h2)

    integer :: i,j,ilo,jlo,ihi,jhi, ipassive
    integer :: n, nq
    
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
          if (rel <= ZERO .or. pl <= small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres
             eos_state % T   = small_temp
             eos_state % rho = rl
             eos_state % xn  = ql(i,j,QFS:QFS-1+nspec)
             eos_state % aux = ql(i,j,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             rel = rl*eos_state%e
             pl  = eos_state%p
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

          if (rer <= ZERO .or. pr <= small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres
             eos_state % T   = small_temp
             eos_state % rho = rr
             eos_state % xn  = qr(i,j,QFS:QFS-1+nspec)
             eos_state % aux = qr(i,j,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

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

             wl = ONE / sqrt(wlsq)
             wr = ONE / sqrt(wrsq)
             
             ustnm1 = ustarm
             ustnp1 = ustarp
             
             ustarm = ur-(pr-pstar)*wr
             ustarp = ul+(pl-pstar)*wl
             
             dpditer=abs(pstnm1-pstar)
             
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
             !reo = rel
             gamco = gcl
             gameo = gamel
             
          else if (ustar .lt. ZERO) then
             ro = rr
             uo = ur
             po = pr
             tauo = taur
             !reo = rer
             gamco = gcr
             gameo = gamer
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             tauo = HALF*(taul+taur)
             !reo = HALF*(rel+rer)
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

          !entho = (reo/ro + po/ro)/co**2
          !estar = reo + (pstar - po)*entho
          
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
          ! if (spout-spin .eq. ZERO) then
          !    scr = small*cav(i,j)
          ! else
          !    scr = spout-spin
          ! endif
          ! frac = (ONE + (spout + spin)/scr)*HALF
          ! frac = max(ZERO,min(ONE,frac))

          frac = HALF*(ONE + (spin + spout)/max(spout-spin,spin+spout, small*cav(i,j)))

          ! the transverse velocity states only depend on the
          ! direction that the contact moves
          if (ustar .gt. ZERO) then
             vgdnv = vl
          else if (ustar .lt. ZERO) then
             vgdnv = vr
          else
             vgdnv = HALF*(vl+vr)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          rgdnv = frac*rstar + (ONE - frac)*ro          
          ugdnv(i,j) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo          

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. ZERO) then
             rgdnv = ro
             ugdnv(i,j) = uo
             pgdnv(i,j) = po
             gamgdnv = gameo 
          endif
          if (spin .ge. ZERO) then
             rgdnv = rstar
             ugdnv(i,j) = ustar
             pgdnv(i,j) = pstar
             gamgdnv = gamstar
          endif

          gegdnv(i,j) = gamgdnv

          pgdnv(i,j) = max(pgdnv(i,j),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                 (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                 (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                 (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                 (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
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
          rhoetot = pgdnv(i,j)/(gamgdnv - ONE) + &
               HALF*rgdnv*(ugdnv(i,j)**2 + vgdnv**2)

          uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv(i,j))
          uflx(i,j,UEINT) = ugdnv(i,j)*pgdnv(i,j)/(gamgdnv - ONE)

          ! advected quantities -- only the contact matters
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nq = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else
                qavg = HALF * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
                   
       enddo
    enddo

  end subroutine riemanncg

  subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)

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


  subroutine riemannus(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                       gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                       uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                       pgdnv, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                       ugdnv, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                       gegdnv, ggd_l1, ggd_l2, ggd_h1, ggd_h2, &
                       idir, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall 
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, &
                                   URHO, UMX, UMY, UEDEN, UEINT, &
                                   small_dens, small_pres, &
                                   npassive, upass_map, qpass_map    

    implicit none

    double precision, parameter:: small = 1.d-8

    integer :: qpd_l1, qpd_l2, qpd_h1, qpd_h2
    integer :: gd_l1, gd_l2, gd_h1, gd_h2
    integer :: uflx_l1, uflx_l2, uflx_h1, uflx_h2
    integer :: pgd_l1, pgd_l2, pgd_h1, pgd_h2
    integer :: ugd_l1, ugd_l2, ugd_h1, ugd_h2
    integer :: ggd_l1, ggd_l2, ggd_h1, ggd_h2
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
    double precision :: gegdnv(ggd_l1:ggd_h1,ggd_l2:ggd_h2)
    
    integer :: ilo,ihi,jlo,jhi
    integer :: n, nq
    integer :: i, j, ipassive

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
             vgd = vl
          else if (ustar .lt. ZERO) then
             vgd = vr
          else
             vgd = HALF*(vl+vr)
          endif
          rgd = frac*rstar + (ONE - frac)*ro
          
          ugdnv(i,j) = frac*ustar + (ONE - frac)*uo
          pgdnv(i,j) = frac*pstar + (ONE - frac)*po

          regd = frac*estar + (ONE - frac)*reo
          if (spout .lt. ZERO) then
             rgd = ro
             ugdnv(i,j) = uo
             pgdnv(i,j) = po
             regd = reo
          endif
          if (spin .ge. ZERO) then
             rgd = rstar
             ugdnv(i,j) = ustar
             pgdnv(i,j) = pstar
             regd = estar
          endif

          gegdnv(i,j) = pgdnv(i,j)/regd + 1.0d0


          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if (idir .eq. 1) then
             if (i.eq.domlo(1) .and. &
                 (physbc_lo(1) .eq. Symmetry .or.  physbc_lo(1) .eq. SlipWall .or. &
                  physbc_lo(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
             if (i.eq.domhi(1)+1 .and. &
                 (physbc_hi(1) .eq. Symmetry .or.  physbc_hi(1) .eq. SlipWall .or. &
                  physbc_hi(1) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
          end if
          if (idir .eq. 2) then
             if (j.eq.domlo(2) .and. &
                 (physbc_lo(2) .eq. Symmetry .or.  physbc_lo(2) .eq. SlipWall .or. &
                  physbc_lo(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
             if (j.eq.domhi(2)+1 .and. &
                 (physbc_hi(2) .eq. Symmetry .or.  physbc_hi(2) .eq. SlipWall .or. &
                  physbc_hi(2) .eq. NoSlipWall) ) &
                  ugdnv(i,j) = ZERO
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
          
          rhoetot = regd + HALF*rgd*(ugdnv(i,j)**2 + vgd**2)
          uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv(i,j))
          uflx(i,j,UEINT) = ugdnv(i,j)*regd
          
          do ipassive = 1, npassive
             n  = upass_map(ipassive)
             nq = qpass_map(ipassive)

             if (ustar .gt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
             else if (ustar .lt. ZERO) then
                uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
             else 
                qavg = HALF * (ql(i,j,nq) + qr(i,j,nq))
                uflx(i,j,n) = uflx(i,j,URHO)*qavg
             endif
          enddo
                    
       enddo
    enddo
  end subroutine riemannus

  subroutine HLL(ql, qr, cl, cr, idir, f)
 
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, &
                                   URHO, UMX, UMY, UEDEN, UEINT, &
                                   npassive, upass_map, qpass_map

    use network, only : nspec, naux

    double precision, intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    double precision, intent(inout) :: f(NVAR)
    integer, intent(in) :: idir

    integer :: ivel, ivelt, imom, imomt
    double precision :: a1, a4, bd, bl, bm, bp, br
    double precision :: cavg, uavg
    double precision :: fl_tmp, fr_tmp
    double precision :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    double precision, parameter :: small = 1.d-10

    if (idir == 1) then
       ivel = QU
       ivelt = QV

       imom = UMX
       imomt = UMY
    else
       ivel = QV
       ivelt = QU

       imom = UMY
       imomt = UMX
    endif

    rhol_sqrt = sqrt(ql(QRHO))
    rhor_sqrt = sqrt(qr(QRHO))

    rhod = ONE/(rhol_sqrt + rhor_sqrt)

    

    ! compute the average sound speed. This uses an approximation from
    ! E88, eq. 5.6, 5.7 that assumes gamma falls between 1
    ! and 5/3
    cavg = sqrt( (rhol_sqrt*cl**2 + rhor_sqrt*cr**2)*rhod + &
         HALF*rhol_sqrt*rhor_sqrt*rhod**2*(qr(ivel) - ql(ivel))**2 )


    ! Roe eigenvalues (E91, eq. 5.3b)
    uavg = (rhol_sqrt*ql(ivel) + rhor_sqrt*qr(ivel))*rhod

    a1 = uavg - cavg
    a4 = uavg + cavg
    

    ! signal speeds (E91, eq. 4.5)
    bl = min(a1, ql(ivel) - cl)
    br = max(a4, qr(ivel) + cr)

    bm = min(ZERO, bl)
    bp = max(ZERO, br)

    bd = bp - bm

    if (abs(bd) < small*max(abs(bm),abs(bp))) return

    bd = ONE/bd
    

    ! compute the fluxes according to E91, eq. 4.4b -- note that the
    ! min/max above picks the correct flux if we are not in the star
    ! region

    ! density flux
    fl_tmp = ql(QRHO)*ql(ivel)
    fr_tmp = qr(QRHO)*qr(ivel)

    f(URHO) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO) - ql(QRHO))


    ! normal momentum flux -- leave out the pressure term -- we handle
    ! that separately
    fl_tmp = ql(QRHO)*ql(ivel)**2
    fr_tmp = qr(QRHO)*qr(ivel)**2

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2)
    fr_tmp = qr(ivel)*(rhoEr + qr(QPRES))

    f(UEDEN) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl)


    ! eint flux
    fl_tmp = ql(QREINT)*ql(ivel)
    fr_tmp = qr(QREINT)*qr(ivel)

    f(UEINT) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QREINT) - ql(QREINT))


    ! passively-advected scalar fluxes
    do ipassive = 1, npassive
       n  = upass_map(ipassive)
       nq = qpass_map(ipassive)

       fl_tmp = ql(QRHO)*ql(nq)*ql(ivel)
       fr_tmp = qr(QRHO)*qr(nq)*qr(ivel)

       f(n) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(nq) - ql(QRHO)*ql(nq))
    enddo

  end subroutine HLL

end module riemann_module
