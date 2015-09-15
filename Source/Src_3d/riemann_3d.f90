module riemann_module

  use bl_constants_module

  implicit none

  private

  public cmpflx, shock

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                    flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                    ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                    gamc,csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                    shk,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                    idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QFS, QFX, QPRES, QREINT, &
                                   use_colglaz, ppm_temp_fix, hybrid_riemann, &
                                   small_temp, allow_negative_energy
    use bl_constants_module

    implicit none

    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
    integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
    integer idir,ilo,ihi,jlo,jhi
    integer i,j,kc,kflux,k3d
    integer domlo(3),domhi(3)

    ! note: qm, qp, ugdnv, pgdnv, gegdnv come in as planes (all
    ! of x,y zones but only 2 elements in the z dir) instead of being
    ! dimensioned the same as the full box.  We index these with kc
    double precision qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
    double precision gegdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)

    ! flux either comes in as planes (like qm, qp, ... above), or
    ! comes in dimensioned as the full box.  We index the flux
    ! with kflux -- this will be set correctly for the different
    ! cases.
    double precision flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)
    
    ! gamc, csml, c, shk come in dimensioned as the full box, so we
    ! use k3d here to index it in z
    double precision gamc(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision    c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  shk(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3)
    
    double precision, pointer :: smallc(:,:),cavg(:,:)
    double precision, pointer :: gamcm(:,:),gamcp(:,:)

    type (eos_t_3D) :: eos_state

    integer :: n
    integer :: is_shock
    double precision :: cl, cr

    double precision :: rhoInv
    
    call eos_allocate(eos_state, (/ ilo, jlo, 1 /), (/ ihi, jhi, 1 /) )
    
    call bl_allocate ( smallc, ilo,ihi,jlo,jhi)
    call bl_allocate (   cavg, ilo,ihi,jlo,jhi)
    call bl_allocate (  gamcm, ilo,ihi,jlo,jhi)
    call bl_allocate (  gamcp, ilo,ihi,jlo,jhi)

    if(idir.eq.1) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i-1,j,k3d) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i-1,j,k3d) )
             gamcm(i,j) = gamc(i-1,j,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
          enddo
       enddo
    elseif(idir.eq.2) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i,j-1,k3d) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i,j-1,k3d) )
             gamcm(i,j) = gamc(i,j-1,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
          enddo
       enddo
    else
       do j = jlo, jhi
          !dir$ ivdep
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
             eos_state % T(i,j,1) = 10000.0d0

             ! minus state
             eos_state % rho(i,j,1)   = qm(i,j,kc,QRHO)
             eos_state % p(i,j,1)     = qm(i,j,kc,QPRES)
             eos_state % xn(i,j,1,:)  = qm(i,j,kc,QFS:QFS+nspec-1)
             eos_state % aux(i,j,1,:) = qm(i,j,kc,QFX:QFX+naux-1)
             
          enddo
       enddo

       call eos(eos_input_rp, eos_state)

       do j = jlo, jhi
          do i = ilo, ihi

             qm(i,j,kc,QREINT) = eos_state % e(i,j,1) * eos_state % rho(i,j,1)
             qm(i,j,kc,QPRES)  = eos_state % p(i,j,1)
             gamcm(i,j)        = eos_state % gam1(i,j,1)

          enddo
       enddo

       ! plus state
       do j = jlo, jhi
          do i = ilo, ihi
             rhoInv = ONE / qp(i,j,kc,QRHO)
             
             eos_state % rho(i,j,1)   = qp(i,j,kc,QRHO)
             eos_state % p(i,j,1)     = qp(i,j,kc,QPRES)
             eos_state % xn(i,j,1,:)  = qp(i,j,kc,QFS:QFS+nspec-1) * rhoInv
             eos_state % aux(i,j,1,:) = qp(i,j,kc,QFX:QFX+naux-1) * rhoInv

          enddo
       enddo

       call eos(eos_input_rp, eos_state)

       do j = jlo, jhi
          do i = ilo, ihi

             qp(i,j,kc,QREINT) = eos_state % e(i,j,1) * eos_state % rho(i,j,1)
             qp(i,j,kc,QPRES)  = eos_state % p(i,j,1)
             gamcp(i,j)        = eos_state % gam1(i,j,1)

          enddo
       enddo

    endif

    ! Solve Riemann problem
    if (use_colglaz == 1) then
       call riemanncg(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                      gamcm,gamcp,cavg,smallc,ilo,jlo,ihi,jhi, &
                      flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                      ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                      idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    else
       call riemannus(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                      gamcm,gamcp,cavg,smallc,ilo,jlo,ihi,jhi, &
                      flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                      ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                      idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)
    endif


    if (hybrid_riemann == 1) then
       ! correct the fluxes using an HLL scheme if we are in a shock
       ! and doing the hybrid approach
       do j = jlo, jhi
          do i = ilo, ihi

             select case (idir)
             case (1)
                is_shock = shk(i-1,j,k3d) + shk(i,j,k3d)
             case (2)
                is_shock = shk(i,j-1,k3d) + shk(i,j,k3d)
             case (3)
                is_shock = shk(i,j,k3d-1) + shk(i,j,k3d)
             end select

             if (is_shock >= 1) then

                select case (idir)
                case (1)
                   cl = c(i-1,j,k3d)
                   cr = c(i,j,k3d)
                case (2)
                   cl = c(i,j-1,k3d)
                   cr = c(i,j,k3d)
                case (3)
                   cl = c(i,j,k3d-1)
                   cr = c(i,j,k3d)
                end select

                call HLL(qm(i,j,kc,:), qp(i,j,kc,:), cl, cr, &
                         idir, flx(i,j,kflux,:))

             endif

          enddo
       enddo

    endif

    call bl_deallocate(smallc)
    call bl_deallocate(  cavg)
    call bl_deallocate( gamcm)
    call bl_deallocate( gamcp)

    call eos_deallocate(eos_state)

  end subroutine cmpflx


  subroutine shock(q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   shk,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                   ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV, QW, QPRES, QVAR
    use bl_constants_module

    integer, intent(in) :: qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
    integer, intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer, intent(in) :: ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
    double precision, intent(in) :: dx, dy, dz
    double precision, intent(in) :: q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision, intent(inout) :: shk(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3)

    integer :: i, j, k

    double precision :: dxinv, dyinv, dzinv
    double precision :: divU
    double precision :: px_pre, px_post, py_pre, py_post, pz_pre, pz_post
    double precision :: e_x, e_y, e_z, d
    double precision :: p_pre, p_post, pjump

    double precision, parameter :: small = 1.d-10
    double precision, parameter :: eps = 0.33d0

    ! This is a basic multi-dimensional shock detection algorithm.
    ! This implementation follows Flash, which in turn follows
    ! AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
    !
    ! The spirit of this follows the shock detection in Colella &
    ! Woodward (1984)

    dxinv = ONE/dx
    dyinv = ONE/dy
    dzinv = ONE/dz

    if (coord_type /= 0) then
       call bl_error("ERROR: invalid geometry in shock()")
    endif

    do k = ilo3-1, ihi3+1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! construct div{U}
             divU = HALF*(q(i+1,j,k,QU) - q(i-1,j,k,QU))*dxinv + &
                    HALF*(q(i,j+1,k,QV) - q(i,j-1,k,QV))*dyinv + &
                    HALF*(q(i,j,k+1,QW) - q(i,j,k-1,QW))*dzinv

             ! find the pre- and post-shock pressures in each direction
             if (q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES) < ZERO) then
                px_pre  = q(i+1,j,k,QPRES)
                px_post = q(i-1,j,k,QPRES)
             else
                px_pre  = q(i-1,j,k,QPRES)
                px_post = q(i+1,j,k,QPRES)
             endif

             if (q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES) < ZERO) then
                py_pre  = q(i,j+1,k,QPRES)
                py_post = q(i,j-1,k,QPRES)
             else
                py_pre  = q(i,j-1,k,QPRES)
                py_post = q(i,j+1,k,QPRES)
             endif

             if (q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES) < ZERO) then
                pz_pre  = q(i,j,k+1,QPRES)
                pz_post = q(i,j,k-1,QPRES)
             else
                pz_pre  = q(i,j,k-1,QPRES)
                pz_post = q(i,j,k+1,QPRES)
             endif

             ! use compression to create unit vectors for the shock direction
             e_x = (q(i+1,j,k,QU) - q(i-1,j,k,QU))**2
             e_y = (q(i,j+1,k,QV) - q(i,j-1,k,QV))**2
             e_z = (q(i,j,k+1,QW) - q(i,j,k-1,QW))**2
             d = ONE/(e_x + e_y + e_z + small)

             e_x = e_x*d
             e_y = e_y*d
             e_z = e_z*d

             ! project the pressures onto the shock direction
             p_pre  = e_x*px_pre + e_y*py_pre + e_z*pz_pre
             p_post = e_x*px_post + e_y*py_post + e_z*pz_post

             ! test for compression + pressure jump to flag a shock
             pjump = eps - (p_post - p_pre)/p_pre

             if (pjump < ZERO .and. divU < ZERO) then
                shk(i,j,k) = ONE
             else
                shk(i,j,k) = ZERO
             endif

          enddo
       enddo
    enddo

  end subroutine shock



! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemanncg(ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                       gamcl,gamcr,cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                       uflx,uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                       ugdnv,pgdnv,gegdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use mempool_module, only : bl_allocate, bl_deallocate
    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QREINT, QESGS, QFS, &
                                   QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UESGS, UFS, UFX, &
                                   small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol, &
                                   npassive, upass_map, qpass_map
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

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so in this routine
    !             it takes values only of  1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nq, ipassive

    double precision :: rgdnv,v1gdnv,v2gdnv,ustar,gamgdnv
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

    double precision, pointer :: pstar_hist(:)

    type (eos_t) :: eos_state

    double precision, pointer :: us1d(:)

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    double precision :: bnd_fac_x, bnd_fac_y, bnd_fac_z

    if (idir .eq. 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ
    else if (idir .eq. 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       im1 = UMY
       im2 = UMX
       im3 = UMZ
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       im1 = UMZ
       im2 = UMX
       im3 = UMY
    end if

    special_bnd_lo = (physbc_lo(idir) .eq. Symmetry &
         .or.         physbc_lo(idir) .eq. SlipWall &
         .or.         physbc_lo(idir) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) .eq. Symmetry &
         .or.         physbc_hi(idir) .eq. SlipWall &
         .or.         physbc_hi(idir) .eq. NoSlipWall)

    if (idir .eq. 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    bnd_fac_z = ONE
    if (idir.eq.3) then
       if ( k3d .eq. domlo(3)   .and. special_bnd_lo .or. &
            k3d .eq. domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = ZERO
       end if
    end if

    tol = cg_tol
    iter_max = cg_maxiter

    call bl_allocate(pstar_hist, 1,iter_max)
    call bl_allocate(us1d, ilo,ihi)

    do j = jlo, jhi

       bnd_fac_y = ONE
       if (idir .eq. 2) then
          if ( j .eq. domlo(2)   .and. special_bnd_lo .or. &
               j .eq. domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = ZERO
          end if
       end if

       do i = ilo, ihi

          ! left state
          rl = max(ql(i,j,kc,QRHO),small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = ql(i,j,kc,QPRES)
          rel = ql(i,j,kc,QREINT)
          gcl = gamcl(i,j)

          ! sometime we come in here with negative energy or pressure
          ! note: reset both in either case, to remain thermo
          ! consistent
          if (rel <= ZERO .or. pl < small_pres) then
             print *, "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: ", rel, pl, small_pres

             eos_state % T   = small_temp
             eos_state % rho = rl
             eos_state % xn  = ql(i,j,kc,QFS:QFS-1+nspec)
             eos_state % aux = ql(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             rel = rl*eos_state % e
             pl  = eos_state % p
             gcl = eos_state % gam1
          endif

          ! right state
          rr = max(qr(i,j,kc,QRHO),small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = qr(i,j,kc,QPRES)
          rer = qr(i,j,kc,QREINT)
          gcr = gamcr(i,j)

          if (rer <= ZERO .or. pr < small_pres) then
             print *, "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: ", rer, pr, small_pres

             eos_state % T   = small_temp
             eos_state % rho = rr
             eos_state % xn  = qr(i,j,kc,QFS:QFS-1+nspec)
             eos_state % aux = qr(i,j,kc,QFX:QFX-1+naux)

             call eos(eos_input_rt, eos_state)

             rer = rr*eos_state % e
             pr  = eos_state % p
             gcr = eos_state % gam1
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

          gegdnv(i,j,kc) = gamgdnv

          pgdnv(i,j,kc) = max(pgdnv(i,j,kc),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          ugdnv(i,j,kc) = ugdnv(i,j,kc) * bnd_fac_x*bnd_fac_y*bnd_fac_z

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*v1gdnv
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*v2gdnv

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

             uflx(i,j,kflux,im1) = uflx(i,j,kflux,im1) + rho_K_contrib

             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + ugdnv(i,j,kc) * rho_K_contrib
          end if

          us1d(i) = ustar
       end do

       ! advected quantities -- only the contact matters
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)

          do i = ilo, ihi
             if (us1d(i) .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (us1d(i) .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    call bl_deallocate(pstar_hist)
    call bl_deallocate(us1d)

    call eos_deallocate(eos_state)

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
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use network, only : nspec, naux
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QESGS, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, &
                                     small_dens, small_pres, npassive, upass_map, qpass_map
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

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nq, ipassive

    double precision :: rgdnv,v1gdnv,v2gdnv,regdnv
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar, ustar
    double precision :: ro, uo, po, reo, co, gamco, entho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg
    double precision :: rho_K_contrib

    double precision, pointer :: us1d(:)

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    double precision :: bnd_fac_x, bnd_fac_y, bnd_fac_z
    double precision :: wwinv, roinv, co2inv

    call bl_allocate(us1d,ilo,ihi)

    if (idir .eq. 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       im1 = UMX
       im2 = UMY
       im3 = UMZ
    else if (idir .eq. 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       im1 = UMY
       im2 = UMX
       im3 = UMZ
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       im1 = UMZ
       im2 = UMX
       im3 = UMY
    end if

    special_bnd_lo = (physbc_lo(idir) .eq. Symmetry &
         .or.         physbc_lo(idir) .eq. SlipWall &
         .or.         physbc_lo(idir) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) .eq. Symmetry &
         .or.         physbc_hi(idir) .eq. SlipWall &
         .or.         physbc_hi(idir) .eq. NoSlipWall)

    if (idir .eq. 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    bnd_fac_z = ONE
    if (idir.eq.3) then
       if ( k3d .eq. domlo(3)   .and. special_bnd_lo .or. &
            k3d .eq. domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = ZERO
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = ONE
       if (idir .eq. 2) then
          if ( j .eq. domlo(2)   .and. special_bnd_lo .or. &
               j .eq. domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = ZERO
          end if
       end if

       !dir$ ivdep
       do i = ilo, ihi

          rl = max(ql(i,j,kc,QRHO),small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

          pl  = max(ql(i,j,kc,QPRES ),small_pres)
          rel =     ql(i,j,kc,QREINT)

          rr = max(qr(i,j,kc,QRHO),small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

          pr  = max(qr(i,j,kc,QPRES),small_pres)
          rer =     qr(i,j,kc,QREINT)

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv
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

          roinv = ONE/ro
          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall,co)
          co2inv = ONE/(co*co)
          entho = (reo + po)*roinv*co2inv
          rstar = ro + (pstar - po)*co2inv
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

          gegdnv(i,j,kc) = pgdnv(i,j,kc)/regdnv + 1.0d0

          pgdnv(i,j,kc) = max(pgdnv(i,j,kc),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          ugdnv(i,j,kc) = ugdnv(i,j,kc) * bnd_fac_x*bnd_fac_y*bnd_fac_z

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = rgdnv*ugdnv(i,j,kc)

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*ugdnv(i,j,kc) + pgdnv(i,j,kc)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*v1gdnv
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*v2gdnv

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

             uflx(i,j,kflux,im1) = uflx(i,j,kflux,im1) + rho_K_contrib

             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + ugdnv(i,j,kc) * rho_K_contrib
          end if

          us1d(i) = ustar
       end do

       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)

          !dir$ ivdep
          do i = ilo, ihi
             if (us1d(i) .gt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
             else if (us1d(i) .lt. ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg
             endif
          enddo

       enddo
    enddo

    call bl_deallocate(us1d)

  end subroutine riemannus


  subroutine HLL(ql, qr, cl, cr, idir, f)

    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   npassive, upass_map, qpass_map

    use network, only : nspec, naux

    double precision, intent(in) :: ql(QVAR), qr(QVAR), cl, cr
    double precision, intent(inout) :: f(NVAR)
    integer, intent(in) :: idir

    integer :: ivel, ivelt, iveltt, imom, imomt, imomtt
    double precision :: a1, a4, bd, bl, bm, bp, br
    double precision :: cavg, uavg
    double precision :: fl_tmp, fr_tmp
    double precision :: rhod, rhoEl, rhoEr, rhol_sqrt, rhor_sqrt
    integer :: n, nq

    integer :: ipassive

    double precision, parameter :: small = 1.d-10

    select case (idir)
    case (1)
       ivel = QU
       ivelt = QV
       iveltt = QW

       imom = UMX
       imomt = UMY
       imomtt = UMZ

    case (2)
       ivel = QV
       ivelt = QU
       iveltt = QW

       imom = UMY
       imomt = UMX
       imomtt = UMZ

    case (3)
       ivel = QW
       ivelt = QU
       iveltt = QV

       imom = UMZ
       imomt = UMX
       imomtt = UMY

    end select

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


    ! normal momentum flux -- we handle that separately
    fl_tmp = ql(QRHO)*ql(ivel)**2 + ql(QPRES)
    fr_tmp = qr(QRHO)*qr(ivel)**2 + qr(QPRES)

    f(imom) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivel) - ql(QRHO)*ql(ivel))


    ! transverse momentum flux
    fl_tmp = ql(QRHO)*ql(ivel)*ql(ivelt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(ivelt)

    f(imomt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(ivelt) - ql(QRHO)*ql(ivelt))


    fl_tmp = ql(QRHO)*ql(ivel)*ql(iveltt)
    fr_tmp = qr(QRHO)*qr(ivel)*qr(iveltt)

    f(imomtt) = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr(QRHO)*qr(iveltt) - ql(QRHO)*ql(iveltt))


    ! total energy flux
    rhoEl = ql(QREINT) + HALF*ql(QRHO)*(ql(ivel)**2 + ql(ivelt)**2 + ql(iveltt)**2)
    fl_tmp = ql(ivel)*(rhoEl + ql(QPRES))

    rhoEr = qr(QREINT) + HALF*qr(QRHO)*(qr(ivel)**2 + qr(ivelt)**2 + qr(iveltt)**2)
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
