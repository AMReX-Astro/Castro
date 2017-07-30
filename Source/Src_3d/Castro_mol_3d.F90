! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(time, &
                               lo, hi, domlo, domhi, &
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
                               flux2, flux2_lo, flux2_hi, &
                               flux3, flux3_lo, flux3_hi, &
                               area1, area1_lo, area1_hi, &
                               area2, area2_lo, area2_hi, &
                               area3, area3_lo, area3_hi, &
                               vol, vol_lo, vol_hi, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                 UTEMP, UEINT, USHK, GDU, GDV, GDW, &
                                 use_flattening, QPRES, NQAUX, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens, ppm_type
  use advection_util_3d_module, only : divu, normalize_species_fluxes
  use advection_util_module, only : compute_cfl, limit_hydro_fluxes_on_small_dens
  use bl_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: uflatten
  use riemann_module, only: cmpflx, shock
  use ppm_module, only : ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
#ifdef RADIATION  
  use rad_params_module, only : ngroups
#endif
#ifdef HYBRID_MOMENTUM
    use hybrid_advection_module, only : add_hybrid_advection_source
#endif

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
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: flux3_lo(3), flux3_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
  integer, intent(in) :: area3_lo(3), area3_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt)        , intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
  real(rt)        , intent(inout) :: update_flux(uf_lo(1):uf_hi(1), uf_lo(2):uf_hi(2), uf_lo(3):uf_hi(3), NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt)        , intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
  real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
  real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
  real(rt)        , intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
  real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt)        , intent(in) :: dx(3), dt, time
  real(rt)        , intent(inout) :: courno

  ! Automatic arrays for workspace
  real(rt)        , pointer:: flatn(:,:,:)
  real(rt)        , pointer:: div(:,:,:)
  real(rt)        , pointer:: pdivu(:,:,:)

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
  real(rt)        , pointer :: sxm(:,:,:,:), sym(:,:,:,:), szm(:,:,:,:)
  real(rt)        , pointer :: sxp(:,:,:,:), syp(:,:,:,:), szp(:,:,:,:)

  real(rt)        , pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
  real(rt)        , pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

  integer :: ngf
  integer :: It_lo(3), It_hi(3)
  integer :: st_lo(3), st_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: div1
  integer :: i, j, k, n
  integer :: kc, km, kt, k3d

  ngf = 1

  It_lo = [lo(1) - 1, lo(2) - 1, 1]
  It_hi = [hi(1) + 1, hi(2) + 1, 2]

  st_lo = [lo(1) - 2, lo(2) - 2, 1]
  st_hi = [hi(1) + 2, hi(2) + 2, 2]

  shk_lo(:) = lo(:) - 1
  shk_hi(:) = hi(:) + 1

  call bl_allocate(   div, lo(1), hi(1)+1, lo(2), hi(2)+1, lo(3), hi(3)+1)
  call bl_allocate( pdivu, lo(1), hi(1)  , lo(2), hi(2)  , lo(3), hi(3)  )

  call bl_allocate(q1, flux1_lo, flux1_hi, NGDNV)
  call bl_allocate(q2, flux2_lo, flux2_hi, NGDNV)
  call bl_allocate(q3, flux3_lo, flux3_hi, NGDNV)

#ifdef RADIATION
  ! when we do radiation, these would be passed out
  call bl_allocate(rflx, flux1_lo, flux1_hi, ngroups)
  call bl_allocate(rfly, flux2_lo, flux2_hi, ngroups)
  call bl_allocate(rflz, flux3_lo, flux3_hi, ngroups)
#endif

  call bl_allocate(sxm, st_lo, st_hi, NQ)
  call bl_allocate(sxp, st_lo, st_hi, NQ)
  call bl_allocate(sym, st_lo, st_hi, NQ)
  call bl_allocate(syp, st_lo, st_hi, NQ)
  call bl_allocate(szm, st_lo, st_hi, NQ)
  call bl_allocate(szp, st_lo, st_hi, NQ)

  call bl_allocate ( qxm, It_lo, It_hi, NQ)
  call bl_allocate ( qxp, It_lo, It_hi, NQ)

  call bl_allocate ( qym, It_lo, It_hi, NQ)
  call bl_allocate ( qyp, It_lo, It_hi, NQ)

  call bl_allocate ( qzm, It_lo, It_hi, NQ)
  call bl_allocate ( qzp, It_lo, It_hi, NQ)

  call bl_allocate(qint, It_lo, It_hi, NGDNV)

  call bl_allocate(shk, shk_lo, shk_hi)

  if (ppm_type == 0) then
     call bl_error("ERROR: method of lines integration does not support ppm_type = 0")
  endif
  
#ifdef SHOCK_VAR
    uout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),USHK) = ZERO

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

  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo, hi, dt, dx, courno)

  ! Compute flattening coefficient for slope calculations.
  call bl_allocate( flatn, q_lo, q_hi)

  if (first_order_hydro == 1) then
     flatn = ZERO
  elseif (use_flattening == 1) then
     call uflatten(lo - ngf, hi + ngf, &
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
  kc = 1
  km = 2

  do k3d = lo(3)-1, hi(3)+1

     ! Swap pointers to levels
     kt = km
     km = kc
     kc = kt

     do n = 1, NQ
        call ppm_reconstruct(q(:,:,:,n  ), q_lo, q_hi, &
                             flatn, q_lo, q_hi, &
                             sxm(:,:,:,n), sxp(:,:,:,n), sym(:,:,:,n), &
                             syp(:,:,:,n), szm(:,:,:,n), szp(:,:,:,n), st_lo, st_hi, &
                             lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

        ! Construct the interface states -- this is essentially just a
        ! reshuffling of interface states from zone-center indexing to
        ! edge-centered indexing
        do j = lo(2)-1, hi(2)+1
           do i = lo(1)-1, hi(1)+1

              ! x-edges

              ! left state at i-1/2 interface
              qxm(i,j,kc,n) = sxp(i-1,j,kc,n)

              ! right state at i-1/2 interface
              qxp(i,j,kc,n) = sxm(i,j,kc,n)

              ! y-edges

              ! left state at j-1/2 interface
              qym(i,j,kc,n) = syp(i,j-1,kc,n)

              ! right state at j-1/2 interface
              qyp(i,j,kc,n) = sym(i,j,kc,n)

              ! z-edges

              ! left state at k3d-1/2 interface
              qzm(i,j,kc,n) = szp(i,j,km,n)

              ! right state at k3d-1/2 interface
              qzp(i,j,kc,n) = szm(i,j,kc,n)

           enddo
        enddo

     enddo

     ! use T to define p
     if (ppm_temp_fix == 1) then
        do j = lo(2)-1, hi(2)+1
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
        endif  ! hi(3) check

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

     endif
     

  enddo

  call bl_deallocate(flatn)

  call bl_deallocate(sxm)
  call bl_deallocate(sxp)
  call bl_deallocate(sym)
  call bl_deallocate(syp)
  call bl_deallocate(szm)
  call bl_deallocate(szp)

  call bl_deallocate(qxm)
  call bl_deallocate(qxp)

  call bl_deallocate(qym)
  call bl_deallocate(qyp)

  call bl_deallocate(qzm)
  call bl_deallocate(qzp)

  call bl_deallocate(qint)
  call bl_deallocate(shk)


  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo,hi,q,q_lo,q_hi,dx,div,lo,hi+1)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           pdivu(i,j,k) = &
                HALF*(q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                     (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1) + &
                HALF*(q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                     (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2) + &
                HALF*(q3(i,j,k+1,GDPRES) + q3(i,j,k,GDPRES)) * &
                     (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)
        enddo
     enddo
  enddo

  do n = 1, NVAR

     if ( n == UTEMP ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
        flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO

#ifdef SHOCK_VAR
     else if ( n == USHK ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
        flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

     else
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)+1
                 div1 = FOURTH*(div(i,j,k) + div(i,j+1,k) + &
                                div(i,j,k+1) + div(i,j+1,k+1))
                 div1 = difmag*min(ZERO,div1)

                 flux1(i,j,k,n) = flux1(i,j,k,n) + &
                      dx(1) * div1 * (uin(i,j,k,n)-uin(i-1,j,k,n))
              enddo
           enddo
        enddo

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)+1
              do i = lo(1), hi(1)
                 div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                div(i,j,k+1) + div(i+1,j,k+1))
                 div1 = difmag*min(ZERO,div1)

                 flux2(i,j,k,n) = flux2(i,j,k,n) + &
                      dx(2) * div1 * (uin(i,j,k,n)-uin(i,j-1,k,n))
              enddo
           enddo
        enddo

        do k = lo(3), hi(3)+1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                div(i,j+1,k) + div(i+1,j+1,k))
                 div1 = difmag*min(ZERO,div1)

                 flux3(i,j,k,n) = flux3(i,j,k,n) + &
                      dx(3) * div1 * (uin(i,j,k,n)-uin(i,j,k-1,n))
              enddo
           enddo
        enddo

     endif

  enddo

  if (limit_fluxes_on_small_dens == 1) then
     call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
                                           q,q_lo,q_hi, &
                                           vol,vol_lo,vol_hi, &
                                           flux1,flux1_lo,flux1_hi, &
                                           area1,area1_lo,area1_hi, &
                                           flux2,flux2_lo,flux2_hi, &
                                           area2,area2_lo,area2_hi, &
                                           flux3,flux3_lo,flux3_hi, &
                                           area3,area3_lo,area3_hi, &
                                           lo,hi,dt,dx)

  endif

  call normalize_species_fluxes(flux1,flux1_lo,flux1_hi, &
                                flux2,flux2_lo,flux2_hi, &
                                flux3,flux3_lo,flux3_hi, &
                                lo,hi)

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                    flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                    flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)

              ! Add the p div(u) source term to (rho e).
              if (n .eq. UEINT) then
                 update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
              endif

              ! for storage
              update_flux(i,j,k,n) = update_flux(i,j,k,n) + &
                   stage_weight * update(i,j,k,n)

              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

           enddo
        enddo
     enddo
  enddo

#ifdef HYBRID_MOMENTUM
  call add_hybrid_advection_source(lo, hi, dt, &
                                   update, uout_lo, uout_hi, &
                                   q1, flux1_lo, flux1_hi, &
                                   q2, flux2_lo, flux2_hi, &
                                   q3, flux3_lo, flux3_hi)
#endif

  ! Scale the fluxes for the form we expect later in refluxing.

  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1
              flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
           enddo
        enddo
     enddo
  enddo

  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2) + 1
           do i = lo(1), hi(1)
              flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
           enddo
        enddo
     enddo
  enddo

  do n = 1, NVAR
     do k = lo(3), hi(3) + 1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              flux3(i,j,k,n) = dt * flux3(i,j,k,n) * area3(i,j,k)
           enddo
        enddo
     enddo
  enddo


  call bl_deallocate(   div)
  call bl_deallocate( pdivu)

  call bl_deallocate(    q1)
  call bl_deallocate(    q2)
  call bl_deallocate(    q3)

#ifdef RADIATION
  call bl_deallocate(rflx)
  call bl_deallocate(rfly)
  call bl_deallocate(rflz)
#endif

end subroutine ca_mol_single_stage
