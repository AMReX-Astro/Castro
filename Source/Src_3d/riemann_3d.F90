module riemann_module

  use bl_types
  use bl_constants_module
  use riemann_util_module
  use meth_params_module, only : NQ, NQAUX, NVAR, QRHO, QU, QV, QW, &
                                 QPRES, QGAME, QREINT, QFS, &
                                 QFX, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                 UFS, UFX, &
                                 NGDNV, GDRHO, GDPRES, GDGAME, &
                                 QC, QCSML, QGAMC, &
#ifdef RADIATION
                                 GDERADS, GDLAMS, QGAMCG, QLAMS, &
                                 qrad, qradhi, qptot, qreitot, fspace_type, &
#endif
                                 small_dens, small_pres, small_temp, &
                                 cg_maxiter, cg_tol, cg_blend, &
                                 npassive, upass_map, qpass_map, &
                                 riemann_solver, ppm_temp_fix, hybrid_riemann, &
                                 allow_negative_energy
#ifdef RADIATION
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor
  use rad_params_module, only : ngroups
#endif

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  private

  public cmpflx

  real(rt), parameter :: smallu = 1.e-12_rt

contains

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine cmpflx(qm, qp, qpd_lo, qpd_hi, &
                    flx, flx_lo, flx_hi, &
                    qint, q_lo, q_hi, &
#ifdef RADIATION
                    rflx, rflx_lo, rflx_hi, &
#endif
                    qaux, qa_lo, qa_hi, &
                    shk, s_lo, s_hi, &
                    idir, ilo, ihi, jlo, jhi, kc, kflux, k3d, domlo, domhi)

    use actual_riemann_module
    use mempool_module, only : bl_allocate, bl_deallocate
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use network, only: nspec, naux
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: qpd_lo(3), qpd_hi(3)
    integer, intent(in) :: flx_lo(3), flx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)

#ifdef RADIATION
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

    real(rt)        , intent(inout) :: qm(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt)        , intent(inout) :: qp(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    real(rt)        , intent(inout) ::    flx(flx_lo(1):flx_hi(1),flx_lo(2):flx_hi(2),flx_lo(3):flx_hi(3),NVAR)
    real(rt)        , intent(inout) ::   qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

#ifdef RADIATION
    real(rt)        , intent(inout) :: rflx(rflx_lo(1):rflx_hi(1), rflx_lo(2):rflx_hi(2), rflx_lo(3):rflx_hi(3),0:ngroups-1)
#endif

    ! qaux come in dimensioned as the full box, so we use k3d here to
    ! index it in z

    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt)        , intent(in) ::  shk(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    ! local variables

    integer i, j
    integer :: gd_lo(2), gd_hi(2)
    real(rt)        , pointer :: smallc(:,:), cavg(:,:)
    real(rt)        , pointer :: gamcm(:,:), gamcp(:,:)
#ifdef RADIATION
    real(rt)        , pointer :: gamcgm(:,:), gamcgp(:,:)
    real(rt)        , pointer :: lam(:,:,:,:)
#endif

    integer :: is_shock
    real(rt)         :: cl, cr
    type (eos_t) :: eos_state

    gd_lo = [ ilo, jlo ]
    gd_hi = [ ihi, jhi ]

    call bl_allocate ( smallc, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (   cavg, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (  gamcm, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate (  gamcp, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
#ifdef RADIATION
    call bl_allocate ( gamcgm, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate ( gamcgp, gd_lo(1),gd_hi(1),gd_lo(2),gd_hi(2))
    call bl_allocate ( lam, qa_lo(1), qa_hi(1), qa_lo(2), qa_hi(2), qa_lo(3), qa_hi(3), 0, ngroups-1)
#endif

    if (idir == 1) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i-1,j,k3d,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i-1,j,k3d,QC) )
             gamcm(i,j) = qaux(i-1,j,k3d,QGAMC)
             gamcp(i,j) = qaux(i,j,k3d,QGAMC)
#ifdef RADIATION
             gamcgm(i,j) = qaux(i-1,j,k3d,QGAMCG)
             gamcgp(i,j) = qaux(i,j,k3d,QGAMCG)
#endif
          enddo
       enddo
    elseif (idir == 2) then
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i,j-1,k3d,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i,j-1,k3d,QC) )
             gamcm(i,j) = qaux(i,j-1,k3d,QGAMC)
             gamcp(i,j) = qaux(i,j,k3d,QGAMC)
#ifdef RADIATION
             gamcgm(i,j) = qaux(i,j-1,k3d,QGAMCG)
             gamcgp(i,j) = qaux(i,j,k3d,QGAMCG)
#endif
          enddo
       enddo
    else
       do j = jlo, jhi
          !dir$ ivdep
          do i = ilo, ihi
             smallc(i,j) = max( qaux(i,j,k3d,QCSML), qaux(i,j,k3d-1,QCSML) )
             cavg(i,j) = HALF*( qaux(i,j,k3d,QC) + qaux(i,j,k3d-1,QC) )
             gamcm(i,j) = qaux(i,j,k3d-1,QGAMC)
             gamcp(i,j) = qaux(i,j,k3d,QGAMC)
#ifdef RADIATION
             gamcgm(i,j) = qaux(i,j,k3d-1,QGAMCG)
             gamcgp(i,j) = qaux(i,j,k3d,QGAMCG)
#endif
          enddo
       enddo
    endif

#ifdef RADIATION
    lam(:,:,:,0:ngroups-1) = qaux(:,:,:,QLAMS:QLAMS+ngroups-1)
#endif

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
             eos_state % T = 10000.0e0_rt

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

             eos_state % rho = qp(i,j,kc,QRHO)
             eos_state % p   = qp(i,j,kc,QPRES)
             eos_state % e   = qp(i,j,kc,QREINT)/qp(i,j,kc,QRHO)
             eos_state % xn  = qp(i,j,kc,QFS:QFS+nspec-1)
             eos_state % aux = qp(i,j,kc,QFX:QFX+naux-1)

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
                      lam, qa_lo, qa_hi, &
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
                   cl = qaux(i-1,j,k3d,QC)
                   cr = qaux(i,j,k3d,QC)
                case (2)
                   cl = qaux(i,j-1,k3d,QC)
                   cr = qaux(i,j,k3d,QC)
                case (3)
                   cl = qaux(i,j,k3d-1,QC)
                   cr = qaux(i,j,k3d,QC)
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
#ifdef RADIATION
    call bl_deallocate(gamcgm)
    call bl_deallocate(gamcgp)
    call bl_deallocate(lam)
#endif

  end subroutine cmpflx


! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine riemannus(ql,qr,qpd_lo,qpd_hi, &
                       gamcl,gamcr,cav,smallc,gd_lo,gd_hi, &
                       uflx,uflx_lo,uflx_hi, &
                       qint,q_lo,q_hi, &
#ifdef RADIATION
                       lam, lam_lo, lam_hi, &
                       gamcgl, gamcgr, &
                       rflx, rflx_lo,rflx_hi, &
#endif
                       idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,domlo,domhi)

    use mempool_module, only : bl_allocate, bl_deallocate
    use prob_params_module, only : physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : compute_hybrid_flux
#endif

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , parameter:: small = 1.e-8_rt

    integer :: qpd_lo(3),qpd_hi(3)
    integer :: gd_lo(2),gd_hi(2)
    integer :: uflx_lo(3),uflx_hi(3)
    integer :: q_lo(3),q_hi(3)
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

#ifdef RADIATION
    integer lam_lo(3), lam_hi(3)
    integer rflx_lo(3),rflx_hi(3)
#endif

    real(rt)         :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt)         :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)

    real(rt)         ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt)         ::    qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

#ifdef RADIATION
    real(rt)         lam(lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2),lam_lo(3):lam_hi(3),0:ngroups-1)
    real(rt)         gamcgl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         gamcgr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         rflx(rflx_lo(1):rflx_hi(1),rflx_lo(2):rflx_hi(2),rflx_lo(3):rflx_hi(3),0:ngroups-1)
#endif

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d
    integer :: n, nqp, ipassive

    real(rt)         :: regdnv
    real(rt)         :: rl, ul, v1l, v2l, pl, rel
    real(rt)         :: rr, ur, v1r, v2r, pr, rer
    real(rt)         :: wl, wr, rhoetot, scr
    real(rt)         :: rstar, cstar, estar, pstar, ustar
    real(rt)         :: ro, uo, po, reo, co, gamco, entho, drho
    real(rt)         :: sgnm, spin, spout, ushock, frac
    real(rt)         :: wsmall, csmall,qavg

#ifdef RADIATION
    real(rt)        , dimension(0:ngroups-1) :: erl, err
    real(rt)         :: reo_g, po_g, co_g, gamco_g
    real(rt)         :: pl_g, rel_g, pr_g, rer_g
    real(rt)         :: regdnv_g, pgdnv_g, pgdnv_t
    real(rt)         :: estar_g, pstar_g
    real(rt)        , dimension(0:ngroups-1) :: lambda, laml, lamr, reo_r, po_r, estar_r, regdnv_r
    real(rt)         :: eddf, f1
    integer :: g
#endif

    real(rt)        , pointer :: us1d(:)

#ifdef ROTATION
    real(rt)         :: vel(3)
#endif

    real(rt)         :: u_adv

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    real(rt)         :: bnd_fac_x, bnd_fac_y, bnd_fac_z
    real(rt)         :: wwinv, roinv, co2inv

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
                lambda(g) = 2.0e0_rt*(laml(g)*lamr(g))/(laml(g)+lamr(g)+1.e-50_rt)
             end do

             reo_r(:) = 0.5e0_rt*(erl(:)+err(:))
             reo_g = 0.5e0_rt*(rel_g+rer_g)
             po_r(:) = lambda(:) * reo_r(:)
             gamco_g = 0.5e0_rt*(gamcgl(i,j)+gamcgr(i,j))
             po_g = 0.5*(pr_g+pl_g)
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
          pgdnv_t       = frac*pstar + (1.e0_rt - frac)*po
          pgdnv_g       = frac*pstar_g + (1.e0_rt - frac)*po_g
          regdnv_g = frac*estar_g + (1.e0_rt - frac)*reo_g
          regdnv_r(:) = frac*estar_r(:) + (1.e0_rt - frac)*reo_r(:)
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
             qint(i,j,kc,GDERADS+g) = max(regdnv_r(g), 0.e0_rt)
          end do

          qint(i,j,kc,GDPRES) = pgdnv_g
          qint(i,j,kc,GDLAMS:GDLAMS-1+ngroups) = lambda(:)

          qint(i,j,kc,GDGAME) = pgdnv_g/regdnv_g + ONE
#else
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE
          qint(i,j,kc,GDPRES) = max(qint(i,j,kc,GDPRES),small_pres)
#endif

          u_adv = qint(i,j,kc,iu)

          ! Enforce that fluxes through a symmetry plane or wall are hard zero.
          if ( special_bnd_lo_x .and. i.eq.domlo(1) .or. &
               special_bnd_hi_x .and. i.eq.domhi(1)+1 ) then
             bnd_fac_x = ZERO
          else
             bnd_fac_x = ONE
          end if
          u_adv = u_adv * bnd_fac_x*bnd_fac_y*bnd_fac_z


          ! Compute fluxes, order as conserved state (not q)
          uflx(i,j,kflux,URHO) = qint(i,j,kc,GDRHO)*u_adv

          uflx(i,j,kflux,im1) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iu ) + qint(i,j,kc,GDPRES)
          uflx(i,j,kflux,im2) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv1)
          uflx(i,j,kflux,im3) = uflx(i,j,kflux,URHO)*qint(i,j,kc,iv2)

#ifdef HYBRID_MOMENTUM
          call compute_hybrid_flux(qint(i,j,kc,:), uflx(i,j,kflux,:), idir, [i, j, k3d])
#endif

#ifdef RADIATION
          rhoetot = regdnv_g + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + pgdnv_g)

          uflx(i,j,kflux,UEINT) = u_adv*regdnv_g
#else
          rhoetot = regdnv + HALF*qint(i,j,kc,GDRHO)*(qint(i,j,kc,iu)**2 + qint(i,j,kc,iv1)**2 + qint(i,j,kc,iv2)**2)

          uflx(i,j,kflux,UEDEN) = u_adv*(rhoetot + qint(i,j,kc,GDPRES))
          uflx(i,j,kflux,UEINT) = u_adv*regdnv
#endif

          ! store this for vectorization
          us1d(i) = ustar

#ifdef RADIATION
          if (fspace_type.eq.1) then
             do g=0,ngroups-1
                eddf = Edd_factor(lambda(g))
                f1 = 0.5e0_rt*(1.e0_rt-eddf)
                rflx(i,j,kflux,g) = (1.e0_rt+f1) * qint(i,j,kc,GDERADS+g) * u_adv
             end do
          else ! type 2
             do g=0,ngroups-1
                rflx(i,j,kflux,g) = qint(i,j,kc,GDERADS+g) * u_adv
             end do
          end if
#endif
       end do

       ! passively advected quantities
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nqp = qpass_map(ipassive)

          !dir$ ivdep
          do i = ilo, ihi
             if (us1d(i) > ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nqp)

             else if (us1d(i) < ZERO) then
                uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nqp)

             else
                qavg = HALF * (ql(i,j,kc,nqp) + qr(i,j,kc,nqp))
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

    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , parameter:: small = 1.e-8_rt

    integer :: qpd_lo(3),qpd_hi(3)
    integer :: gd_lo(2),gd_hi(2)
    integer :: uflx_lo(3),uflx_hi(3)
    integer :: q_lo(3),q_hi(3)
    integer :: idir,ilo,ihi,jlo,jhi
    integer :: domlo(3),domhi(3)

    real(rt)         :: ql(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt)         :: qr(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),NQ)
    real(rt)         ::  gamcl(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         ::  gamcr(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         ::    cav(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         :: smallc(gd_lo(1):gd_hi(1),gd_lo(2):gd_hi(2))
    real(rt)         :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt)         :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NGDNV)

    ! Note:  Here k3d is the k corresponding to the full 3d array --
    !         it should be used for print statements or tests against domlo, domhi, etc
    !         kc  is the k corresponding to the 2-wide slab of k-planes, so takes values
    !             only of 1 or 2
    !         kflux is used for indexing into the uflx array -- in the initial calls to
    !             cmpflx when uflx = {fx,fy,fxy,fyx,fz,fxz,fzx,fyz,fzy}, kflux = kc,
    !             but in later calls, when uflx = {flux1,flux2,flux3}  , kflux = k3d
    integer :: i,j,kc,kflux,k3d

    real(rt)         :: rgdnv,regdnv
    real(rt)         :: rl, ul, v1l, v2l, pl, rel
    real(rt)         :: rr, ur, v1r, v2r, pr, rer
    real(rt)         :: wl, wr, scr
    real(rt)         :: rstar, cstar, estar, pstar, ustar
    real(rt)         :: ro, uo, po, reo, co, gamco, entho
    real(rt)         :: sgnm, spin, spout, ushock, frac
    real(rt)         :: wsmall, csmall

    integer :: iu, iv1, iv2, im1, im2, im3
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt)         :: wwinv, roinv, co2inv

    real(rt)         :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt)         :: S_l, S_r, S_c

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
             reo = rel
             gamco = gamcl(i,j)

          else if (ustar < ZERO) then
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

          rstar = ro + (pstar - po)*co2inv
          rstar = max(small_dens,rstar)

          entho = (reo + po)*co2inv/ro
          estar = reo + (pstar - po)*entho

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
          regdnv = frac*estar + (ONE - frac)*reo

          qint(i,j,kc,iu) = frac*ustar + (ONE - frac)*uo
          qint(i,j,kc,GDPRES) = frac*pstar + (ONE - frac)*po
          qint(i,j,kc,GDGAME) = qint(i,j,kc,GDPRES)/regdnv + ONE


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
             call compute_flux(idir, bnd_fac, U_state, pr, F_state)

          else if (S_r > ZERO .and. S_c <= ZERO) then
             ! R* region
             call cons_state(qr(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pr, F_state)

             call HLLC_state(idir, S_r, S_c, qr(i,j,kc,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

          else if (S_c > ZERO .and. S_l < ZERO) then
             ! L* region
             call cons_state(ql(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pl, F_state)

             call HLLC_state(idir, S_l, S_c, ql(i,j,kc,:), U_hllc_state)

             ! correct the flux
             F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

          else
             ! L region
             call cons_state(ql(i,j,kc,:), U_state)
             call compute_flux(idir, bnd_fac, U_state, pl, F_state)

          endif

          uflx(i,j,kflux,:) = F_state(:)
       enddo
    enddo

  end subroutine HLLC

end module riemann_module
