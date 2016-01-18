module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module

  implicit none

  private

  public cmpflx, shock

  real (kind=dp_t), parameter :: smallu = 1.e-12_dp_t

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, q_lo, q_hi, &
#ifdef RADIATION
                    lam, lam_lo, lam_hi, &
                    rflx, rflx_lo, rflx_hi, &
                    gamcg, &
#endif
                    gamc,csml,c,qd_lo,qd_hi, &
                    shk,s_lo,s_hi, &
                    idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QFS, QFX, QPRES, QREINT, NGDNV, &
                                   riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                   small_temp, allow_negative_energy
#ifdef RADIATION
    use radhydro_params_module, only : QRADVAR
    use rad_params_module, only : ngroups
#endif

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)

#ifdef RADIATION
    integer lam_lo(3), lam_hi(3)
    integer rflx_lo(3), rflx_hi(3)
#endif

    integer, intent(in) :: idir,ilo,ihi,jlo,jhi,kc,kflux,k3d
    integer, intent(in) :: domlo(3),domhi(3)

    ! note: qm, qp, q come in as planes (all of x,y
    ! zones but only 2 elements in the z dir) instead of being
    ! dimensioned the same as the full box.  We index these with kc
    ! flux either comes in as planes (like qm, qp, ... above), or
    ! comes in dimensioned as the full box.  We index the flux with
    ! kflux -- this will be set correctly for the different cases.

#ifdef RADIATION
    double precision, intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QRADVAR)
    double precision, intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QRADVAR)
#else
    double precision, intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
    double precision, intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
#endif

    double precision, intent(inout) ::    flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    double precision, intent(inout) ::   qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision rflx(rflx_lo(1):rflx_hi(1), rflx_lo(2):rflx_hi(2), rflx_lo(3):rflx_hi(3),0:ngroups-1)
    double precision gamcg(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
#endif

    ! gamc, csml, c, shk come in dimensioned as the full box, so we
    ! use k3d here to index it in z
    double precision, intent(inout) :: gamc(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(inout) ::    c(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(inout) :: csml(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
    double precision, intent(inout) ::  shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    ! local variables
    integer i, j
    integer :: gd_lo(2), gd_hi(2)
    double precision, pointer :: smallc(:,:), cavg(:,:)
    double precision, pointer :: gamcm(:,:), gamcp(:,:)
#ifdef RADIATION
    double precision, pointer :: gamcgm(:,:), gamcgp(:,:)
#endif

    integer :: is_shock
    double precision :: cl, cr
    type (eos_t) :: eos_state

    double precision :: rhoInv

    gd_lo = (/ ilo, jlo /)
    gd_hi = (/ ihi, jhi /)
    
    call bl_allocate ( smallc, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (   cavg, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (  gamcm, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (  gamcp, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
#ifdef RADIATION
    call bl_allocate ( gamcgm, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate ( gamcgp, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
#endif

    if (idir == 1) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i-1,j,k3d) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i-1,j,k3d) )
             gamcm(i,j) = gamc(i-1,j,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
#ifdef RADIATION
             gamcgm(i,j) = gamcg(i-1,j,k3d)
             gamcgp(i,j) = gamcg(i,j,k3d)
#endif
          enddo
       enddo
    elseif (idir == 2) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( csml(i,j,k3d), csml(i,j-1,k3d) )
             cavg(i,j) = HALF*( c(i,j,k3d) + c(i,j-1,k3d) )
             gamcm(i,j) = gamc(i,j-1,k3d)
             gamcp(i,j) = gamc(i,j,k3d)
#ifdef RADIATION
             gamcgm(i,j) = gamcg(i,j-1,k3d)
             gamcgp(i,j) = gamcg(i,j,k3d)
#endif
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
#ifdef RADIATION
             gamcgm(i,j) = gamcg(i,j,k3d-1)
             gamcgp(i,j) = gamcg(i,j,k3d)
#endif
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
             eos_state % T = 10000.0d0

             ! minus state
             eos_state % rho = qm(i,j,kc,QRHO)
             eos_state % p   = qm(i,j,kc,QPRES)
             eos_state % e   = qm(i,j,kc,QREINT)/qm(i,j,kc,QRHO) 
             eos_state % xn  = qm(i,j,kc,QFS:QFS+nspec-1)
             eos_state % aux = qm(i,j,kc,QFX:QFX+naux-1)

             call eos(eos_input_re, eos_state)

             qm(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             qm(i,j,kc,QPRES)  = eos_state % p
             gamcm(i,j)        = eos_state % gam1

          enddo
       enddo

       ! plus state
       do j = jlo, jhi
          do i = ilo, ihi
             rhoInv = ONE / qp(i,j,kc,QRHO)

             eos_state % rho = qp(i,j,kc,QRHO)
             eos_state % p   = qp(i,j,kc,QPRES)
             eos_state % e   = qp(i,j,kc,QREINT)/qp(i,j,kc,QRHO)
             eos_state % xn  = qp(i,j,kc,QFS:QFS+nspec-1) * rhoInv
             eos_state % aux = qp(i,j,kc,QFX:QFX+naux-1) * rhoInv

             call eos(eos_input_re, eos_state)

             qp(i,j,kc,QREINT) = eos_state % e * eos_state % rho
             qp(i,j,kc,QPRES)  = eos_state % p
             gamcp(i,j)        = eos_state % gam1

          enddo
       enddo

    endif

    ! Solve Riemann problem
    if (riemann_solver == 0) then
       ! Colella, Glaz, & Ferguson solver
       call riemannus(qm, qp, qpd_lo, qpd_hi, &
                      gamcm, gamcp, cavg, smallc, gd_lo, gd_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, q_lo, q_hi, &
#ifdef RADIATION
                      lam, lam_lo, lam_hi, &
                      gamcgm, gamcgp, &
                      rflx, rflx_lo, rflx_hi, &
#endif
                      idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    elseif (riemann_solver == 1) then
       ! Colella & Glaz solver
       call riemanncg(qm, qp, qpd_lo, qpd_hi, &
                      gamcm, gamcp, cavg, smallc, gd_lo, gd_hi, &
                      flx, flx_lo, flx_hi, &
                      qint, q_lo, q_hi, &
                      idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    elseif (riemann_solver == 2) then
       ! HLLC
       call HLLC(qm, qp, qpd_lo, qpd_hi, &
                 gamcm, gamcp, cavg, smallc, gd_lo, gd_hi, &
                 flx, flx_lo, flx_hi, &
                 qint, q_lo, q_hi, &
                 idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)
    else
       call bl_error("ERROR: invalid value of riemann_solver")
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
                         idir, 3, flx(i,j,kflux,:))

             endif

          enddo
       enddo

    endif

    call bl_deallocate(smallc)
    call bl_deallocate(  cavg)
    call bl_deallocate( gamcm)
    call bl_deallocate( gamcp)
#ifdef RADIATION
    call bl_deallocate(gamcgm)
    call bl_deallocate(gamcgp)
#endif

  end subroutine cmpflx


  subroutine shock(q,qd_lo,qd_hi,shk,s_lo,s_hi,lo,hi,dx)

    use prob_params_module, only : coord_type
    use meth_params_module, only : QU, QV, QW, QPRES, QVAR
    use bl_constants_module

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    double precision, intent(in) :: dx(3)
    double precision, intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)
    double precision, intent(inout) :: shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

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

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    if (coord_type /= 0) then
       call bl_error("ERROR: invalid geometry in shock()")
    endif

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

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

  subroutine riemanncg(ql,qr,qpd_lo,qpd_hi, &
                       gamcl,gamcr,cav,smallc,gd_lo,gd_hi, &
                       uflx,uflx_lo,uflx_hi, &
                       qint,q_lo,q_hi, &
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    ! this implements the approximate Riemann solver of Colella & Glaz (1985)

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use bl_error_module
    use network, only : nspec, naux
    use eos_type_module
    use eos_module
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, &
                                   QPRES, QGAME, QREINT, QESGS, QFS, &
                                   QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   UESGS, UFS, UFX, &
                                   NGDNV, GDRHO, GDPRES, GDGAME, &
                                   small_dens, small_pres, small_temp, &
                                   cg_maxiter, cg_tol, &
                                   npassive, upass_map, qpass_map

    double precision, parameter:: small = 1.d-8

    integer :: qpd_lo(3),qpd_hi(3)
    integer :: gd_lo(2),gd_hi(2)
    integer :: uflx_lo(3),uflx_hi(3)
    integer :: q_lo(3),q_hi(3)
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

    double precision :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
    double precision :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
    double precision ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    double precision :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so in this routine
    !             it takes values only of  1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nq, ipassive

    double precision :: ustar,gamgdnv
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

          ! secant iteration
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

          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

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
             qint(i,j,kc,iv1) = v1l
             qint(i,j,kc,iv2) = v2l
          else if (ustar .lt. ZERO) then
             qint(i,j,kc,iv1) = v1r
             qint(i,j,kc,iv2) = v2r
          else
             qint(i,j,kc,iv1) = HALF*(v1l+v1r)
             qint(i,j,kc,iv2) = HALF*(v2l+v2r)
          endif

          ! linearly interpolate between the star and normal state -- this covers the
          ! case where we are inside the rarefaction fan.
          qint(i,j,kc,GDRHO ) = frac*rstar + (ONE - frac)*ro
          qint(i,j,kc,iu   ) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          gamgdnv =  frac*gamstar + (ONE-frac)*gameo

          ! now handle the cases where instead we are fully in the
          ! star or fully in the original (l/r) state
          if (spout .lt. ZERO) then
             qint(i,j,kc,GDRHO ) = ro
             qint(i,j,kc,iu   ) = uo
             qint(i,j,kc,GDPRES) = po
             gamgdnv = gameo
          endif
          if (spin .ge. ZERO) then
             qint(i,j,kc,GDRHO ) = rstar
             qint(i,j,kc,iu   ) = ustar
             qint(i,j,kc,GDPRES) = pstar
             gamgdnv = gamstar
          endif

          qint(i,j,kc,GDGAME) = gamgdnv

          qint(i,j,kc,GDPRES) = max(qint(i,j,kc,GDPRES),small_pres)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          qint(i,j,kc,iu) = qint(i,j,kc,iu) * bnd_fac_x*bnd_fac_y*bnd_fac_z

          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = qint(i,j,kc,GDRHO)*qint(i,j,kc,iu)

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iu) + qint(i,j,kc,GDPRES)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv2)

          ! compute the total energy from the internal, p/(gamma - 1), and the kinetic
          rhoetot = qint(i,j,kc,GDPRES)/(gamgdnv - ONE) + &
               HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = qint(i,j,kc,iu)*(rhoetot + qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = qint(i,j,kc,iu)*qint(i,j,kc,GDPRES)/(gamgdnv - ONE)


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

             rho_K_contrib =  TWO3RD * qint(i,j,kc,GDRHO) * qavg

             uflx(i,j,kflux,im1) = uflx(i,j,kflux,im1) + rho_K_contrib

             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + qint(i,j,kc,iu) * rho_K_contrib
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

  end subroutine riemanncg

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemannus(ql,qr,qpd_lo,qpd_hi, &
                       gamcl,gamcr,cav,smallc,gd_lo,gd_hi, &
                       uflx,uflx_lo,uflx_hi, &
                       qint,q_lo,q_hi, &
#ifdef RADIATION
                       lam,lam_lo,lam_hi, &
                       gamcgl,gamcgr, &
                       rflx, rflx_lo,rflx_hi, &
#endif
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QGAME, QREINT, QESGS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, &
                                   NGDNV, GDRHO, GDPRES, GDGAME, &
                                   small_dens, small_pres, npassive, upass_map, qpass_map
#ifdef RADIATION
    use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
    use rad_params_module, only : ngroups
    use fluxlimiter_module, only : Edd_factor
#endif

    double precision, parameter:: small = 1.d-8

    integer :: qpd_lo(3),qpd_hi(3)
    integer :: gd_lo(2),gd_hi(2)
    integer :: uflx_lo(3),uflx_hi(3)
    integer :: q_lo(3),q_hi(3)
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

#ifdef RADIATION
    integer lam_lo(3),lam_hi(3)
    integer rflx_lo(3),rflx_hi(3)
#endif

#ifdef RADIATION
    double precision :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QRADVAR)
    double precision :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QRADVAR)
#else
    double precision :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
    double precision :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
#endif
    double precision ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    double precision ::    qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

#ifdef RADIATION
    double precision lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    double precision gamcgl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision gamcgr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision rflx(rflx_lo(1):rflx_hi(1),rflx_lo(2):rflx_hi(2),rflx_lo(3):rflx_hi(3),0:ngroups-1)
#endif

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nq, ipassive

    double precision :: regdnv
    double precision :: rl, ul, v1l, v2l, pl, rel
    double precision :: rr, ur, v1r, v2r, pr, rer
    double precision :: wl, wr, rhoetot, scr
    double precision :: rstar, cstar, estar, pstar, ustar, v1g, v2g
    double precision :: ro, uo, po, reo, co, gamco, entho, drho
    double precision :: sgnm, spin, spout, ushock, frac
    double precision :: wsmall, csmall,qavg
    double precision :: rho_K_contrib

#ifdef RADIATION
    double precision, dimension(0:ngroups-1) :: erl, err
    double precision :: reo_g, po_g, co_g, gamco_g
    double precision :: pl_g, rel_g, pr_g, rer_g
    double precision :: regdnv_g, pgdnv_g, pgdnv_t
    double precision :: estar_g, pstar_g
    double precision, dimension(0:ngroups-1) :: lambda, laml, lamr, reo_r, po_r, estar_r, regdnv_r
    double precision :: eddf, f1
    integer :: g
#endif

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

#ifdef RADIATION
          if (idir == 1) then
             laml = lam(i-1,j,k3d,:)
          elseif (idir == 2) then
             laml = lam(i,j-1,k3d,:)
          else
             laml = lam(i,j,k3d-1,:)
          end if
          lamr = lam(i,j,k3d,:)
#endif

          rl = max(ql(i,j,kc,QRHO), small_dens)

          ! pick left velocities based on direction
          ul  = ql(i,j,kc,iu)
          v1l = ql(i,j,kc,iv1)
          v2l = ql(i,j,kc,iv2)

#ifdef RADIATION
          pl = ql(i,j,kc,qptot)
          rel = ql(i,j,kc,qreitot)
          erl(:) = ql(i,j,kc,qrad:qradhi)
          pl_g = ql(i,j,kc,QPRES)
          rel_g = ql(i,j,kc,QREINT)
#else
          pl  = max(ql(i,j,kc,QPRES ), small_pres)
          rel =     ql(i,j,kc,QREINT)
#endif

          rr = max(qr(i,j,kc,QRHO), small_dens)

          ! pick right velocities based on direction
          ur  = qr(i,j,kc,iu)
          v1r = qr(i,j,kc,iv1)
          v2r = qr(i,j,kc,iv2)

#ifdef RADIATION
          pr = qr(i,j,kc,qptot)
          rer = qr(i,j,kc,qreitot)
          err(:) = qr(i,j,kc,qrad:qradhi)
          pr_g = qr(i,j,kc,QPRES)
          rer_g = qr(i,j,kc,QREINT)
#else
          pr  = max(qr(i,j,kc,QPRES), small_pres)
          rer =     qr(i,j,kc,QREINT)
#endif

          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

          pstar = max(pstar,small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             reo = rel
             gamco = gamcl(i,j)
#ifdef RADIATION
             lambda = laml
             po_g = pl_g
             po_r(:) = erl(:) * lambda(:)
             reo_r(:) = erl(:)
             reo_g = rel_g
             gamco_g = gamcgl(i,j)
#endif

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             reo = rer
             gamco = gamcr(i,j)
#ifdef RADIATION
             lambda = lamr
             po_g = pr_g
             po_r(:) = err(:) * lambda(:)
             reo_r(:) = err(:)
             reo_g = rer_g
             gamco_g = gamcgr(i,j)
#endif

          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             reo = HALF*(rel+rer)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
#ifdef RADIATION
             do g=0, ngroups-1
                lambda(g) = 2.0d0*(laml(g)*lamr(g))/(laml(g)+lamr(g)+1.d-50)
             end do
             po_r(:) = lambda(:) * reo_r(:)
             po_g = 0.5*(pr_g+pl_g)
             reo_r(:) = 0.5d0*(erl(:)+err(:))
             reo_g = 0.5d0*(rel_g+rer_g)
             gamco_g = 0.5d0*(gamcgl(i,j)+gamcgr(i,j))
#endif

          endif

          ro = max(small_dens,ro)

          roinv = ONE/ro

          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall,co)
          co2inv = ONE/(co*co)

          drho = (pstar - po)*co2inv
          rstar = ro + drho
          rstar = max(small_dens,rstar)

#ifdef RADIATION
          estar_g = reo_g + drho*(reo_g + po_g)/ro
          co_g = sqrt(abs(gamco_g*po_g/ro))
          co_g = max(csmall,co_g)
          pstar_g = po_g + drho*co_g**2
          pstar_g = max(pstar_g,small_pres)
          estar_r = reo_r(:) + drho*(reo_r(:) + po_r(:))/ro
#else
          entho = (reo + po)*roinv*co2inv
          estar = reo + (pstar - po)*entho
#endif
          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif

          if (spout-spin == ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif

          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          if (ustar > ZERO) then
             qint(i,j,kc,iv1) = v1l
             qint(i,j,kc,iv2) = v2l
          else if (ustar < ZERO) then
             qint(i,j,kc,iv1) = v1r
             qint(i,j,kc,iv2) = v2r
          else
             qint(i,j,kc,iv1) = HALF*(v1l+v1r)
             qint(i,j,kc,iv2) = HALF*(v2l+v2r)
          endif
          qint(i,j,kc,GDRHO) = frac*rstar + (ONE - frac)*ro
          qint(i,j,kc,iu  ) = frac*ustar + (ONE - frac)*uo

#ifdef RADIATION
          pgdnv_t       = frac*pstar + (1.d0 - frac)*po
          pgdnv_g       = frac*pstar_g + (1.d0 - frac)*po_g
          regdnv_g = frac*estar_g + (1.d0 - frac)*reo_g
          regdnv_r(:) = frac*estar_r(:) + (1.d0 - frac)*reo_r(:)
#else
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          regdnv = frac*estar + (ONE - frac)*reo
#endif
          
          if (spout < ZERO) then
             qint(i,j,kc,GDRHO) = ro
             qint(i,j,kc,iu  ) = uo
#ifdef RADIATION
             pgdnv_t = po
             pgdnv_g = po_g
             regdnv_g = reo_g
             regdnv_r = reo_r(:)
#else
             qint(i,j,kc,GDPRES) = po
             regdnv = reo
#endif
          endif

          if (spin >= ZERO) then
             qint(i,j,kc,GDRHO) = rstar
             qint(i,j,kc,iu  ) = ustar
#ifdef RADIATION
             pgdnv_t = pstar
             pgdnv_g = pstar_g
             regdnv_g = estar_g
             regdnv_r = estar_r(:)
#else
             qint(i,j,kc,GDPRES) = pstar
             regdnv = estar
#endif
          endif


#ifdef RADIATION
          do g=0, ngroups-1
             qint(i,j,kc,GDERADS+g) = max(regdnv_r(g), 0.d0)
          end do
          
          qint(i,j,kc,GDPRES) = pgdnv_g          
          qint(i,j,kc,GDLAMS:GDLAMS-1+ngroups) = lambda(:)

          qint(i,j,kc,GDGAME) = pgdnv_g/regdnv_g + ONE
#else
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE
          qint(i,j,kc,GDPRES) = max(qint(i,j,kc,GDPRES),small_pres)
#endif

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          qint(i,j,kc,iu) = qint(i,j,kc,iu) * bnd_fac_x*bnd_fac_y*bnd_fac_z


          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = qint(i,j,kc,GDRHO)*qint(i,j,kc,iu)

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iu ) + qint(i,j,kc,GDPRES)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv2)

#ifdef RADIATION
          rhoetot = regdnv_g + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)
          
          uflx(i,j,kflux,UEDEN) = qint(i,j,kc,iu)*(rhoetot + pgdnv_g)
          
          uflx(i,j,kflux,UEINT) = qint(i,j,kc,iu)*regdnv_g
#else
          rhoetot = regdnv + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = qint(i,j,kc,iu)*(rhoetot + qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = qint(i,j,kc,iu)*regdnv
#endif

          ! Treat K as a passively advected quantity but allow it to
          ! affect fluxes of (rho E) and momenta.
          if (UESGS > -1) then
             n  = UESGS
             nq = QESGS

             if (ustar > ZERO) then
                qavg = ql(i,j,kc,nq)

             else if (ustar < ZERO) then
                qavg = qr(i,j,kc,nq)

             else
                qavg = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
             endif

             uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg

             rho_K_contrib =  TWO3RD * qint(i,j,kc,GDRHO) * qavg

             uflx(i,j,kflux,im1) = uflx(i,j,kflux,im1) + rho_K_contrib

             uflx(i,j,kflux,UEDEN) = uflx(i,j,kflux,UEDEN) + qint(i,j,kc,iu) * rho_K_contrib
          end if

          ! store this for vectorization
          us1d(i) = ustar

#ifdef RADIATION
          if (fspace_type.eq.1) then
             do g=0,ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5d0*(1.d0-eddf)
                rflx(i,j,kflux,g) = (1.d0+f1) * qint(i,j,kc,GDERADS+g) * qint(i,j,kc,iu)
             end do
          else ! type 2
             do g=0,ngroups-1
                rflx(i,j,kflux,g) = qint(i,j,kc,GDERADS+g) * qint(i,j,kc,iu)
             end do
          end if
#endif
       end do

       ! passively advected quantities
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)

          !dir$ ivdep
          do i = ilo, ihi
             if (us1d(i) > ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)

             else if (us1d(i) < ZERO) then
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


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine HLLC(ql,qr,qpd_lo,qpd_hi, &
                  gamcl,gamcr,cav,smallc,gd_lo,gd_hi, &
                  uflx,uflx_lo,uflx_hi, &
                  qint,q_lo,q_hi, &
                  idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)


    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QESGS, &
                                   URHO, UMX, UMY, UMZ, UEDEN, UEINT, UESGS, &
                                   NGDNV, GDRHO, GDPRES, GDGAME, &
                                   small_dens, small_pres, npassive, upass_map, qpass_map

    double precision, parameter:: small = 1.d-8

    integer :: qpd_lo(3),qpd_hi(3)
    integer :: gd_lo(2),gd_hi(2)
    integer :: uflx_lo(3),uflx_hi(3)
    integer :: q_lo(3),q_hi(3)
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

    double precision :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
    double precision :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),QVAR)
    double precision ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    double precision :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    double precision :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nq, ipassive

    double precision :: rgdnv,regdnv
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
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    double precision :: wwinv, roinv, co2inv

    double precision :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    double precision :: S_l, S_r, S_c

    call bl_allocate(us1d,ilo,ihi)

    if (UESGS > 0) then
       call bl_error("ERROR: HLLC doesn't support SGS")
    endif

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

    bnd_fac_z = 1
    if (idir.eq.3) then
       if ( k3d .eq. domlo(3)   .and. special_bnd_lo .or. &
            k3d .eq. domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac_z = 0
       end if
    end if

    do j = jlo, jhi

       bnd_fac_y = 1
       if (idir .eq. 2) then
          if ( j .eq. domlo(2)   .and. special_bnd_lo .or. &
               j .eq. domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac_y = 0
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

          ! now we essentially do the CGF solver to get p and u on the
          ! interface, but we won't use these in any flux construction.
          csmall = smallc(i,j)
          wsmall = small_dens*csmall
          wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
          wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

          wwinv = ONE/(wl + wr)
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
          ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

          pstar = max(pstar,small_pres)
          ! for symmetry preservation, if ustar is really small, then we
          ! set it to zero
          if (abs(ustar) < smallu*HALF*(abs(ul) + abs(ur))) then
             ustar = ZERO
          endif

          if (ustar > ZERO) then
             ro = rl
             uo = ul
             po = pl
             gamco = gamcl(i,j)

          else if (ustar < ZERO) then
             ro = rr
             uo = ur
             po = pr
             gamco = gamcr(i,j)
          else
             ro = HALF*(rl+rr)
             uo = HALF*(ul+ur)
             po = HALF*(pl+pr)
             gamco = HALF*(gamcl(i,j)+gamcr(i,j))
          endif
          ro = max(small_dens,ro)

          roinv = ONE/ro
          co = sqrt(abs(gamco*po*roinv))
          co = max(csmall,co)
          co2inv = ONE/(co*co)

          rstar = ro + (pstar - po)*co2inv
          rstar = max(small_dens,rstar)

          cstar = sqrt(abs(gamco*pstar/rstar))
          cstar = max(cstar,csmall)

          sgnm = sign(ONE,ustar)
          spout = co - sgnm*uo
          spin = cstar - sgnm*ustar
          ushock = HALF*(spin + spout)

          if (pstar-po > ZERO) then
             spin = ushock
             spout = ushock
          endif
          if (spout-spin == ZERO) then
             scr = small*cav(i,j)
          else
             scr = spout-spin
          endif
          frac = (ONE + (spout + spin)/scr)*HALF
          frac = max(ZERO,min(ONE,frac))

          rgdnv = frac*rstar + (ONE - frac)*ro

          qint(i,j,kc,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po

          ! TODO
          !qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE


          ! now we do the HLLC construction


          ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = 0
          else
             bnd_fac_x = 1
          end if

          bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z
          
          ! use the simplest estimates of the wave speeds
          S_l = min(ul - sqrt(gamcl(i,j)*pl/rl), ur - sqrt(gamcr(i,j)*pr/rr))
          S_r = max(ul + sqrt(gamcl(i,j)*pl/rl), ur + sqrt(gamcr(i,j)*pr/rr))

          ! estimate of the contact speed -- this is Toro Eq. 10.8
          S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
               (rl*(S_l - ul) - rr*(S_r - ur))

          if (S_r <= ZERO) then
             ! R region
             call cons_state(qr(i,j,kc,:), U_state)
             call compute_flux(idir, 3, bnd_fac, U_state, pr, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call cons_state(qr(i,j,kc,:), U_state)
             call compute_flux(idir, 3, bnd_fac, U_state, pr, F_state)
             
             call HLLC_state(idir, S_r, S_c, qr(i,j,kc,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call cons_state(ql(i,j,kc,:), U_state)
             call compute_flux(idir, 3, bnd_fac, U_state, pl, F_state)

             call HLLC_state(idir, S_l, S_c, ql(i,j,kc,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call cons_state(ql(i,j,kc,:), U_state)
             call compute_flux(idir, 3, bnd_fac, U_state, pl, F_state)

          endif

          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

    call bl_deallocate(us1d)

  end subroutine HLLC

end module riemann_module
