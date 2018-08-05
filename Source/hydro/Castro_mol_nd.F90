! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(lo, hi, time, &
                               domlo, domhi, &
                               stage_weight, &
                               uin, uin_lo, uin_hi, &
                               uout, uout_lo, uout_hi, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               srcU, srU_lo, srU_hi, &
                               update, updt_lo, updt_hi, &
                               update_flux, uf_lo, uf_hi, &
                               dx, dt, &
                               flux1, flux1_lo, flux1_hi, &
#if BL_SPACEDIM >= 2
                               flux2, flux2_lo, flux2_hi, &
#endif
#if BL_SPACEDIM == 3
                               flux3, flux3_lo, flux3_hi, &
#endif
                               area1, area1_lo, area1_hi, &
#if BL_SPACEDIM >= 2
                               area2, area2_lo, area2_hi, &
#endif
#if BL_SPACEDIM == 3
                               area3, area3_lo, area3_hi, &
#endif
#if BL_SPACEDIM < 3
                               pradial, p_lo, p_hi, &
                               dloga, dloga_lo, dloga_hi, &
#endif
                               vol, vol_lo, vol_hi, &
                               verbose) bind(C, name="ca_mol_single_stage")

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                 UTEMP, UEINT, USHK, GDU, GDV, GDW, UMX, &
                                 use_flattening, QPRES, NQAUX, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens, ppm_type, ppm_temp_fix
  use advection_util_module, only : limit_hydro_fluxes_on_small_dens, shock, &
                                    divu, normalize_species_fluxes, calc_pdivu
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: uflatten
  use riemann_module, only: cmpflx
  use ppm_module, only : ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
#ifdef RADIATION  
  use rad_params_module, only : ngroups
#endif
#ifdef HYBRID_MOMENTUM
  use hybrid_advection_module, only : add_hybrid_advection_source
#endif
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec, naux
  use prob_params_module, only : dg, coord_type

    
  implicit none

  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) ::  domlo(3), domhi(3)
  real(rt), intent(in) :: stage_weight
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
#if BL_SPACEDIM >= 2
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if BL_SPACEDIM == 3
  integer, intent(in) :: flux3_lo(3), flux3_hi(3)
  integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
#if BL_SPACEDIM <= 2
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
#endif
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
  real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
  real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
  real(rt), intent(inout) :: update_flux(uf_lo(1):uf_hi(1), uf_lo(2):uf_hi(2), uf_lo(3):uf_hi(3), NVAR)
  real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
  real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if BL_SPACEDIM >= 2
  real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if BL_SPACEDIM == 3
  real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
  real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
#if BL_SPACEDIM <= 2
  real(rt), intent(inout) :: pradial(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
  real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1), dloga_lo(2):dloga_hi(2), dloga_lo(3):dloga_hi(3))
#endif
  real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt), intent(in) :: dx(3), dt, time

  ! Automatic arrays for workspace
  real(rt)        , pointer:: flatn(:,:,:)
  real(rt)        , pointer:: div(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , pointer:: q1(:,:,:,:)
  real(rt)        , pointer:: q2(:,:,:,:)
  real(rt)        , pointer:: q3(:,:,:,:)
  real(rt)        , pointer:: qint(:,:,:,:)

#ifdef RADIATION
  ! radiation fluxes (need these to get things to compile)
  real(rt)        , pointer:: rflx(:,:,:,:)
  real(rt)        , pointer:: rfly(:,:,:,:)
  real(rt)        , pointer:: rflz(:,:,:,:)
#endif

  real(rt)        , pointer:: shk(:,:,:)

  ! temporary interface values of the parabola
  real(rt)        , pointer :: sxm(:,:,:), sym(:,:,:), szm(:,:,:)
  real(rt)        , pointer :: sxp(:,:,:), syp(:,:,:), szp(:,:,:)

  real(rt)        , pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
  real(rt)        , pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

  integer :: ngf
  integer :: It_lo(3), It_hi(3)
  integer :: st_lo(3), st_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: div1
  integer :: i, j, k, n
  integer :: kc, km, kt, k3d

  type (eos_t) :: eos_state
  
  ngf = 1

  It_lo = [lo(1) - 1, lo(2) - dg(2), dg(3)]
  It_hi = [hi(1) + 1, hi(2) + dg(2), 2*dg(3)]

  st_lo = [lo(1) - 2, lo(2) - 2*dg(2), dg(3)]
  st_hi = [hi(1) + 2, hi(2) + 2*dg(2), 2*dg(3)]

  shk_lo(:) = lo(:) - dg(:)
  shk_hi(:) = hi(:) + dg(:)

  call bl_allocate(   div, lo(1), hi(1)+1, lo(2), hi(2)+dg(2), lo(3), hi(3)+dg(3))

  call bl_allocate(q1, flux1_lo, flux1_hi, NGDNV)
#if BL_SPACEDIM >= 2
  call bl_allocate(q2, flux2_lo, flux2_hi, NGDNV)
#endif
#if BL_SPACEDIM == 3
  call bl_allocate(q3, flux3_lo, flux3_hi, NGDNV)
#endif

#ifdef RADIATION
  ! when we do radiation, these would be passed out
  call bl_allocate(rflx, flux1_lo, flux1_hi, ngroups)
#if BL_SPACEDIM >= 2
  call bl_allocate(rfly, flux2_lo, flux2_hi, ngroups)
#endif
#if BL_SPACEDIM == 3
  call bl_allocate(rflz, flux3_lo, flux3_hi, ngroups)
#endif
#endif

  call bl_allocate(sxm, st_lo, st_hi)
  call bl_allocate(sxp, st_lo, st_hi)
  call bl_allocate(qxm, It_lo, It_hi, NQ)
  call bl_allocate(qxp, It_lo, It_hi, NQ)

#if BL_SPACEDIM >= 2
  call bl_allocate(sym, st_lo, st_hi)
  call bl_allocate(syp, st_lo, st_hi)
  call bl_allocate(qym, It_lo, It_hi, NQ)
  call bl_allocate(qyp, It_lo, It_hi, NQ)
#endif
#if BL_SPACEDIM == 3
  call bl_allocate(szm, st_lo, st_hi)
  call bl_allocate(szp, st_lo, st_hi)
  call bl_allocate(qzm, It_lo, It_hi, NQ)
  call bl_allocate(qzp, It_lo, It_hi, NQ)
#endif

  call bl_allocate(qint, It_lo, It_hi, NGDNV)

  call bl_allocate(shk, shk_lo, shk_hi)

#ifndef AMREX_USE_CUDA  
  if (ppm_type == 0) then
     call amrex_error("ERROR: method of lines integration does not support ppm_type = 0")
  endif
#endif

#ifdef SHOCK_VAR
    uout(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), USHK) = ZERO

    call shock(q, q_lo, q_hi, shk, shk_lo, shk_hi, lo, hi, dx)

    ! Store the shock data for future use in the burning step.

    do k3d = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uout(i,j,k3d,USHK) = shk(i,j,k3d)
          enddo
       enddo
    enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:,:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, q_lo, q_hi, shk, shk_lo, shk_hi, lo, hi, dx)
    else
       shk(:,:,:) = ZERO
    endif
#endif

  ! Compute flattening coefficient for slope calculations.
  call bl_allocate(flatn, q_lo, q_hi)

  if (first_order_hydro == 1) then
     flatn = ZERO
  elseif (use_flattening == 1) then
     call uflatten(lo - ngf*dg, hi + ngf*dg, &
                   q, flatn, q_lo, q_hi, QPRES)
  else
     flatn = ONE
  endif

  ! We come into this routine with a 3-d box of data, but we operate
  ! on it locally by considering 2 planes that encompass all of the
  ! x, y indices of the original box, but each plane corresponds to
  ! a single z index.
  !
  ! In the notation below, k3d will always been the index into the
  ! original 3-d box.  kc will be the z-index in the local "planar"
  ! data and km will be the previously used index in the local
  ! planar data.
  !
  ! With each loop in the k direction, we will overwrite the old
  ! data in the planar arrays.

  ! Initialize kc (current k-level) and km (previous k-level)
#if BL_SPACEDIM == 3
  kc = 1
  km = 2
#else
  kc = 0
  km = 0
#endif

  do k3d = lo(3)-dg(3), hi(3)+dg(3)

#if BL_SPACEDIM == 3
     ! Swap pointers to levels
     kt = km
     km = kc
     kc = kt
#endif

     do n = 1, NQ
        call ppm_reconstruct(q, q_lo, q_hi, NQ, n, &
                             flatn, q_lo, q_hi, &
                             sxm, sxp, &
#if BL_SPACEDIM >= 2
                             sym, syp, &
#endif
#if BL_SPACEDIM == 3
                             szm, szp, &
#endif
                             st_lo, st_hi, &
                             lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

        ! Construct the interface states -- this is essentially just a
        ! reshuffling of interface states from zone-center indexing to
        ! edge-centered indexing
        do j = lo(2)-dg(2), hi(2)+dg(2)
           do i = lo(1)-1, hi(1)+1

              ! x-edges

              ! left state at i-1/2 interface
              qxm(i,j,kc,n) = sxp(i-1,j,kc)

              ! right state at i-1/2 interface
              qxp(i,j,kc,n) = sxm(i,j,kc)

#if BL_SPACEDIM >= 2
              ! y-edges

              ! left state at j-1/2 interface
              qym(i,j,kc,n) = syp(i,j-1,kc)

              ! right state at j-1/2 interface
              qyp(i,j,kc,n) = sym(i,j,kc)
#endif

#if BL_SPACEDIM == 3
              ! z-edges

              ! left state at k3d-1/2 interface
              qzm(i,j,km,n) = szp(i,j,kc)

              ! right state at k3d-1/2 interface
              qzp(i,j,kc,n) = szm(i,j,kc)
#endif

           enddo
        enddo

     enddo

     ! use T to define p
     if (ppm_temp_fix == 1) then
        do j = lo(2)-dg(2), hi(2)+dg(2)
           do i = lo(1)-1, hi(1)+1

              eos_state%rho    = qxp(i,j,kc,QRHO)
              eos_state%T      = qxp(i,j,kc,QTEMP)
              eos_state%xn(:)  = qxp(i,j,kc,QFS:QFS-1+nspec)
              eos_state%aux(:) = qxp(i,j,kc,QFX:QFX-1+naux)

              call eos(eos_input_rt, eos_state)

              qxp(i,j,kc,QPRES) = eos_state%p
              qxp(i,j,kc,QREINT) = qxp(i,j,kc,QRHO)*eos_state%e
              ! should we try to do something about Gamma_! on interface?

              eos_state%rho    = qxm(i,j,kc,QRHO)
              eos_state%T      = qxm(i,j,kc,QTEMP)
              eos_state%xn(:)  = qxm(i,j,kc,QFS:QFS-1+nspec)
              eos_state%aux(:) = qxm(i,j,kc,QFX:QFX-1+naux)

              call eos(eos_input_rt, eos_state)

              qxm(i,j,kc,QPRES) = eos_state%p
              qxm(i,j,kc,QREINT) = qxm(i,j,kc,QRHO)*eos_state%e
              ! should we try to do something about Gamma_! on interface?

#if BL_SPACEDIM >= 2
              eos_state%rho    = qyp(i,j,kc,QRHO)
              eos_state%T      = qyp(i,j,kc,QTEMP)
              eos_state%xn(:)  = qyp(i,j,kc,QFS:QFS-1+nspec)
              eos_state%aux(:) = qyp(i,j,kc,QFX:QFX-1+naux)

              call eos(eos_input_rt, eos_state)

              qyp(i,j,kc,QPRES) = eos_state%p
              qyp(i,j,kc,QREINT) = qyp(i,j,kc,QRHO)*eos_state%e
              ! should we try to do something about Gamma_! on interface?

              eos_state%rho    = qym(i,j,kc,QRHO)
              eos_state%T      = qym(i,j,kc,QTEMP)
              eos_state%xn(:)  = qym(i,j,kc,QFS:QFS-1+nspec)
              eos_state%aux(:) = qym(i,j,kc,QFX:QFX-1+naux)

              call eos(eos_input_rt, eos_state)

              qym(i,j,kc,QPRES) = eos_state%p
              qym(i,j,kc,QREINT) = qym(i,j,kc,QRHO)*eos_state%e
              ! should we try to do something about Gamma_! on interface?
#endif

#if BL_SPACEDIM == 3
              eos_state%rho    = qzp(i,j,kc,QRHO)
              eos_state%T      = qzp(i,j,kc,QTEMP)
              eos_state%xn(:)  = qzp(i,j,kc,QFS:QFS-1+nspec)
              eos_state%aux(:) = qzp(i,j,kc,QFX:QFX-1+naux)

              call eos(eos_input_rt, eos_state)

              qzp(i,j,kc,QPRES) = eos_state%p
              qzp(i,j,kc,QREINT) = qzp(i,j,kc,QRHO)*eos_state%e
              ! should we try to do something about Gamma_! on interface?

              eos_state%rho    = qzm(i,j,kc,QRHO)
              eos_state%T      = qzm(i,j,kc,QTEMP)
              eos_state%xn(:)  = qzm(i,j,kc,QFS:QFS-1+nspec)
              eos_state%aux(:) = qzm(i,j,kc,QFX:QFX-1+naux)

              call eos(eos_input_rt, eos_state)

              qzm(i,j,kc,QPRES) = eos_state%p
              qzm(i,j,kc,QREINT) = qzm(i,j,kc,QRHO)*eos_state%e
              ! should we try to do something about Gamma_! on interface?
#endif

           enddo
        enddo
     endif


     if (k3d >= lo(3)) then

        ! Compute F^x at kc (k3d)
        if (k3d <= hi(3)) then
           call cmpflx(qxm, qxp, It_lo, It_hi, &
                       flux1, flux1_lo, flux1_hi, &
                       qint, It_lo, It_hi, &  ! temporary
#ifdef RADIATION
                       rflx, flux1_lo, flux1_hi, &
#endif
                       qaux, qa_lo, qa_hi, &
                       shk, shk_lo, shk_hi, &
                       1, lo(1), hi(1)+1, lo(2), hi(2), kc, k3d, k3d, domlo, domhi)

           do j = lo(2), hi(2)
              do i = lo(1), hi(1)+1
                 q1(i,j,k3d,:) = qint(i,j,kc,:)
              enddo
           enddo

#if BL_SPACEDIM >= 2
           ! Compute F^y at kc (k3d)
           call cmpflx(qym, qyp, It_lo, It_hi, &
                       flux2, flux2_lo, flux2_hi, &
                       qint, It_lo, It_hi, &  ! temporary
#ifdef RADIATION
                       rfly, flux2_lo, flux2_hi, &
#endif
                       qaux, qa_lo, qa_hi, &
                       shk, shk_lo, shk_hi, &
                       2, lo(1), hi(1), lo(2), hi(2)+1, kc, k3d, k3d, domlo, domhi)

           do j = lo(2), hi(2)+1
              do i = lo(1), hi(1)
                 q2(i,j,k3d,:) = qint(i,j,kc,:)
              enddo
           enddo
#endif
        endif  ! hi(3) check

#if BL_SPACEDIM == 3
        ! Compute F^z at kc (k3d)

        call cmpflx(qzm, qzp, It_lo, It_hi, &
                    flux3, flux3_lo, flux3_hi, &
                    qint, It_lo, It_hi, &
#ifdef RADIATION
                    rflz, flux3_lo, flux3_hi, &
#endif
                    qaux, qa_lo, qa_hi, &
                    shk, shk_lo, shk_hi, &
                    3, lo(1), hi(1), lo(2), hi(2), kc, k3d, k3d, domlo, domhi)

        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              q3(i,j,k3d,:) = qint(i,j,kc,:)
           enddo
        enddo
#endif

     endif

  enddo

  call bl_deallocate(flatn)

  call bl_deallocate(sxm)
  call bl_deallocate(sxp)
  call bl_deallocate(qxm)
  call bl_deallocate(qxp)

#if BL_SPACEDIM >= 2
  call bl_deallocate(sym)
  call bl_deallocate(syp)
  call bl_deallocate(qym)
  call bl_deallocate(qyp)
#endif

#if BL_SPACEDIM == 3
  call bl_deallocate(szm)
  call bl_deallocate(szp)
  call bl_deallocate(qzm)
  call bl_deallocate(qzp)
#endif

  call bl_deallocate(qint)
  call bl_deallocate(shk)


  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo, hi, q, q_lo, q_hi, &
            dx, div, lo, hi+dg)

  do n = 1, NVAR

     if ( n == UTEMP ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if BL_SPACEDIM >= 2
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if BL_SPACEDIM == 3
        flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

#ifdef SHOCK_VAR
     else if ( n == USHK ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if BL_SPACEDIM >= 2
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if BL_SPACEDIM == 3
        flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif
#endif

     else
        ! do the artificial viscosity
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)+1

                 div1 = FOURTH*(div(i,j,k) + div(i,j+dg(2),k) + &
                                div(i,j,k+dg(3)) + div(i,j+dg(2),k+dg(3)))
                 div1 = difmag*min(ZERO, div1)

                 flux1(i,j,k,n) = flux1(i,j,k,n) + &
                      dx(1) * div1 * (uin(i,j,k,n) - uin(i-1,j,k,n))
              enddo
           enddo
        enddo
#if BL_SPACEDIM >= 2
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)+1
              do i = lo(1), hi(1)
                 div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                div(i,j,k+dg(3)) + div(i+1,j,k+dg(3)))
                 div1 = difmag*min(ZERO, div1)

                 flux2(i,j,k,n) = flux2(i,j,k,n) + &
                      dx(2) * div1 * (uin(i,j,k,n) - uin(i,j-1,k,n))
              enddo
           enddo
        enddo
#endif
#if BL_SPACEDIM == 3
        do k = lo(3), hi(3)+1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                div(i,j+1,k) + div(i+1,j+1,k))
                 div1 = difmag*min(ZERO, div1)

                 flux3(i,j,k,n) = flux3(i,j,k,n) + &
                      dx(3) * div1 * (uin(i,j,k,n) - uin(i,j,k-1,n))
              enddo
           enddo
        enddo
#endif
     endif

  enddo

  if (limit_fluxes_on_small_dens == 1) then
     call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
                                           q,q_lo,q_hi, &
                                           vol,vol_lo,vol_hi, &
                                           flux1,flux1_lo,flux1_hi, &
                                           area1,area1_lo,area1_hi, &
#if BL_SPACEDIM >= 2
                                           flux2,flux2_lo,flux2_hi, &
                                           area2,area2_lo,area2_hi, &
#endif
#if BL_SPACEDIM == 3
                                           flux3,flux3_lo,flux3_hi, &
                                           area3,area3_lo,area3_hi, &
#endif
                                           lo,hi,dt,dx)

  endif

  call normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
#if BL_SPACEDIM >= 2
                                flux2,flux2_lo,flux2_hi, &
#endif
#if BL_SPACEDIM == 3
                                flux3,flux3_lo,flux3_hi, &
#endif
                                lo,hi)

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

#if BL_SPACEDIM == 1
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) ) / vol(i,j,k)

#elif BL_SPACEDIM == 2
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                    flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)

#else
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                    flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                    flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)
#endif

#if BL_SPACEDIM == 1
              if (n == UMX) then
                 update(i,j,k,UMX) = update(i,j,k,UMX) - ( q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES) ) / dx(1)
              endif
#endif

#if BL_SPACEDIM == 2
              if (n == UMX) then
                 ! add the pressure source term for axisymmetry
                 if (coord_type > 0) then
                    update(i,j,k,n) = update(i,j,k,n) - (q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES))/ dx(1)
                 endif
              endif
#endif

              ! for storage
              update_flux(i,j,k,n) = update_flux(i,j,k,n) + &
                   stage_weight * update(i,j,k,n)

              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

           enddo
        enddo
     enddo
  enddo

#if BL_SPACEDIM == 3
#ifdef HYBRID_MOMENTUM
  call add_hybrid_advection_source(lo, hi, dt, &
                                   update, uout_lo, uout_hi, &
                                   q1, flux1_lo, flux1_hi, &
                                   q2, flux2_lo, flux2_hi, &
                                   q3, flux3_lo, flux3_hi)
#endif
#endif



  ! Scale the fluxes for the form we expect later in refluxing.

  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1
              flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)

#if BL_SPACEDIM == 1
              if (coord_type .eq. 0 .and. n == UMX) then
                 flux1(i,j,k,n) = flux1(i,j,k,n) + dt * area1(i,j,k) * q1(i,j,k,GDPRES)
              endif
#endif

           enddo
        enddo
     enddo
  enddo

#if BL_SPACEDIM >= 2
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2) + 1
           do i = lo(1), hi(1)
              flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
           enddo
        enddo
     enddo
  enddo
#endif

#if BL_SPACEDIM == 3
  do n = 1, NVAR
     do k = lo(3), hi(3) + 1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
           enddo
        enddo
     enddo
  enddo
#endif

#if BL_SPACEDIM < 3
  if (coord_type > 0) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
  end if
#endif

  call bl_deallocate(   div)

  call bl_deallocate(q1)
#if BL_SPACEDIM >= 2
  call bl_deallocate(q2)
#endif
#if BL_SPACEDIM == 3
  call bl_deallocate(q3)
#endif

#ifdef RADIATION
  call bl_deallocate(rflx)
#if BL_SPACEDIM >= 2
  call bl_deallocate(rfly)
#endif
#if BL_SPACEDIM == 3
  call bl_deallocate(rflz)
#endif
#endif

end subroutine ca_mol_single_stage
