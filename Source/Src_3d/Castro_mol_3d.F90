! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(time, &
                               lo, hi, domlo, domhi, &
                               uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                               uout, uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
                               q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                               qaux, qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3, &
                               srcU, srU_l1, srU_l2, srU_l3, srU_h1, srU_h2, srU_h3, &
                               update, updt_l1, updt_l2, updt_l3, updt_h1, updt_h2, updt_h3, &
                               dx, dt, &
                               flux1, flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3, &
                               flux2, flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3, &
                               flux3, flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3, &
                               area1, area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3, &
                               area2, area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3, &
                               area3, area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3, &
                               vol, vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                 UTEMP, UEINT, USHK, UMX, GDU, GDV, GDW, &
                                 use_flattening, QU, QV, QW, QPRES, NQAUX, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens
  use advection_util_3d_module, only : divu, normalize_species_fluxes
  use advection_util_module, only : compute_cfl, limit_hydro_fluxes_on_small_dens
  use bl_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: uflaten
  use riemann_module, only: cmpflx, shock
  use ppm_module, only : ppm_reconstruct
  use bl_fort_module, only : rt => c_real

  implicit none

  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) ::  domlo(3), domhi(3)
  integer, intent(in) :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
  integer, intent(in) :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
  integer, intent(in) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer, intent(in) :: qa_l1, qa_l2, qa_l3, qa_h1, qa_h2, qa_h3
  integer, intent(in) :: srU_l1, srU_l2, srU_l3, srU_h1, srU_h2, srU_h3
  integer, intent(in) :: updt_l1, updt_l2, updt_l3, updt_h1, updt_h2, updt_h3
  integer, intent(in) :: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
  integer, intent(in) :: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
  integer, intent(in) :: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
  integer, intent(in) :: area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3
  integer, intent(in) :: area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3
  integer, intent(in) :: area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3
  integer, intent(in) :: vol_l1, vol_l2, vol_l3, vol_h1, vol_h2, vol_h3

  real(rt)        , intent(in) :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
  real(rt)        , intent(inout) :: uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NVAR)
  real(rt)        , intent(inout) :: q(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, NQ)
  real(rt)        , intent(inout) :: qaux(qa_l1:qa_h1, qa_l2:qa_h2, qa_l3:qa_h3, NQAUX)
  real(rt)        , intent(in) :: srcU(srU_l1:srU_h1, srU_l2:srU_h2, srU_l3:srU_h3, QVAR)
  real(rt)        , intent(inout) :: update(updt_l1:updt_h1, updt_l2:updt_h2, updt_l3:updt_h3, NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, flux1_l3:flux1_h3, NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, flux2_l3:flux2_h3, NVAR)
  real(rt)        , intent(inout) :: flux3(flux3_l1:flux3_h1, flux3_l2:flux3_h2, flux3_l3:flux3_h3, NVAR)
  real(rt)        , intent(in) :: area1(area1_l1:area1_h1, area1_l2:area1_h2, area1_l3:area1_h3)
  real(rt)        , intent(in) :: area2(area2_l1:area2_h1, area2_l2:area2_h2, area2_l3:area2_h3)
  real(rt)        , intent(in) :: area3(area3_l1:area3_h1, area3_l2:area3_h2, area3_l3:area3_h3)
  real(rt)        , intent(in) :: vol(vol_l1:vol_h1, vol_l2:vol_h2, vol_l3:vol_h3)
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

  ! radiation fluxes (need these to get things to compile)
  real(rt)        , pointer:: rflx(:,:,:,:)
  real(rt)        , pointer:: rfly(:,:,:,:)
  real(rt)        , pointer:: rflz(:,:,:,:)

  real(rt)        , pointer:: shk(:,:,:)

  ! temporary interface values of the parabola
  real(rt)        , pointer :: sxm(:,:,:), sym(:,:,:), szm(:,:,:)
  real(rt)        , pointer :: sxp(:,:,:), syp(:,:,:), szp(:,:,:)

  real(rt)        , pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
  real(rt)        , pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

  integer :: ngq, ngf
  integer :: uin_lo(3), uin_hi(3)
  integer :: uout_lo(3), uout_hi(3)
  integer :: updt_lo(3), updt_hi(3)
  integer :: flux1_lo(3), flux1_hi(3)
  integer :: flux2_lo(3), flux2_hi(3)
  integer :: flux3_lo(3), flux3_hi(3)
  integer :: area1_lo(3), area1_hi(3)
  integer :: area2_lo(3), area2_hi(3)
  integer :: area3_lo(3), area3_hi(3)
  integer :: vol_lo(3), vol_hi(3)
  integer :: q_lo(3), q_hi(3)
  integer :: qa_lo(3), qa_hi(3)
  integer :: srU_lo(3), srU_hi(3)
  integer :: It_lo(3), It_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: div1
  integer :: i, j, k, n
  integer :: kc, km, kt, k3d

  ngf = 1

  q_lo = [ q_l1, q_l2, q_l3 ]
  q_hi = [ q_h1, q_h2, q_h3 ]

  qa_lo = [ qa_l1, qa_l2, qa_l3 ]
  qa_hi = [ qa_h1, qa_h2, qa_h3 ]

  srU_lo = [ srU_l1, srU_l2, srU_l3 ]
  srU_hi = [ srU_h1, srU_h2, srU_h3 ]

  uin_lo = [ uin_l1, uin_l2, uin_l3 ]
  uin_hi = [ uin_h1, uin_h2, uin_h3 ]

  uout_lo = [ uout_l1, uout_l2, uout_l3 ]
  uout_hi = [ uout_h1, uout_h2, uout_h3 ]

  updt_lo = [ updt_l1, updt_l2, updt_l3 ]
  updt_hi = [ updt_h1, updt_h2, updt_h3 ]

  flux1_lo = [ flux1_l1, flux1_l2, flux1_l3 ]
  flux1_hi = [ flux1_h1, flux1_h2, flux1_h3 ]

  flux2_lo = [ flux2_l1, flux2_l2, flux2_l3 ]
  flux2_hi = [ flux2_h1, flux2_h2, flux2_h3 ]

  flux3_lo = [ flux3_l1, flux3_l2, flux3_l3 ]
  flux3_hi = [ flux3_h1, flux3_h2, flux3_h3 ]

  area1_lo = [ area1_l1, area1_l2, area1_l3 ]
  area1_hi = [ area1_h1, area1_h2, area1_h3 ]

  area2_lo = [ area2_l1, area2_l2, area2_l3 ]
  area2_hi = [ area2_h1, area2_h2, area2_h3 ]

  area3_lo = [ area3_l1, area3_l2, area3_l3 ]
  area3_hi = [ area3_h1, area3_h2, area3_h3 ]

  vol_lo = [ vol_l1, vol_l2, vol_l3 ]
  vol_hi = [ vol_h1, vol_h2, vol_h3 ]

  It_lo = [lo(1) - 1, lo(2) - 1, 1]
  It_hi = [hi(1) + 1, hi(2) + 1, 2]

  shk_lo(:) = lo(:) - 1
  shk_hi(:) = hi(:) + 1

  call bl_allocate(   div, lo(1), hi(1)+1, lo(2), hi(2)+1, lo(3), hi(3)+1)
  call bl_allocate( pdivu, lo(1), hi(1)  , lo(2), hi(2)  , lo(3), hi(3)  )

  call bl_allocate(q1, flux1_lo, flux1_hi, NGDNV)
  call bl_allocate(q2, flux2_lo, flux2_hi, NGDNV)
  call bl_allocate(q3, flux3_lo, flux3_hi, NGDNV)

  call bl_allocate(sxm, It_lo, It_hi)
  call bl_allocate(sxp, It_lo, It_hi)
  call bl_allocate(sym, It_lo, It_hi)
  call bl_allocate(syp, It_lo, It_hi)
  call bl_allocate(szm, It_lo, It_hi)
  call bl_allocate(szp, It_lo, It_hi)

  call bl_allocate ( qxm, It_lo, It_hi, NQ)
  call bl_allocate ( qxp, It_lo, It_hi, NQ)

  call bl_allocate ( qym, It_lo, It_hi, NQ)
  call bl_allocate ( qyp, It_lo, It_hi, NQ)

  call bl_allocate ( qzm, It_lo, It_hi, NQ)
  call bl_allocate ( qzp, It_lo, It_hi, NQ)

  call bl_allocate(qint, It_lo, It_hi, NQ)

  call bl_allocate(shk, shk_lo, shk_hi)
  
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
     call uflaten(lo - ngf, hi + ngf, &
                  q(:,:,:,QPRES), q(:,:,:,QU), q(:,:,:,QV), q(:,:,:,QW), &
                  flatn, q_lo, q_hi)
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
                             sxm, sxp, sym, syp, szm, szp, It_lo, It_hi, &
                             lo(1), lo(2), hi(1), hi(2), dx, k3d, kc)

        ! Construct the interface states -- this is essentially just a
        ! reshuffling of interface states from zone-center indexing to
        ! edge-centered indexing
        if (k3d >= lo(3)) then
           do j = lo(2)-1, hi(2)+1
              do i = lo(1)-1, hi(1)+1

                 ! x-edges

                 ! left state at i-1/2 interface
                 qxm(i,j,kc,n) = sxp(i-1,j,kc)

                 ! right state at i-1/2 interface
                 qxp(i,j,kc,n) = sxm(i,j,kc)

                 ! y-edges

                 ! left state at j-1/2 interface
                 qym(i,j,kc,n) = syp(i,j-1,kc)

                 ! right state at j-1/2 interface
                 qyp(i,j,kc,n) = sym(i,j,kc)

                 ! z-edges

                 ! left state at k3d-1/2 interface
                 qzm(i,j,kc,n) = szp(i,j,km)

                 ! right state at k3d-1/2 interface
                 qzp(i,j,kc,n) = szm(i,j,kc)

              enddo
           enddo

           ! Compute F^x at kc (k3d)
           call cmpflx(qxm, qxp, It_lo, It_hi, &
                       flux1, flux1_lo, flux1_hi, &
                       qint, It_lo, It_hi, &  ! temporary
#ifdef RADIATION
                       rflux1, rfd1_lo, rfd1_hi, &
#endif
                       qaux, qa_lo, qa_hi, &
                       shk, shk_lo, shk_hi, &
                       1, lo(1), hi(1)+1, lo(2), hi(2), kc, k3d, k3d, domlo, domhi)

           do j = lo(2)-1, hi(2)+1
              do i = lo(1)-1, hi(1)+2
                 q1(i,j,k3d,:) = qint(i,j,kc,:)
              enddo
           enddo


           ! Compute F^y at kc (k3d)
           call cmpflx(qym, qyp, It_lo, It_hi, &
                       flux2, flux2_lo, flux2_hi, &
                       qint, It_lo, It_hi, &  ! temporary
#ifdef RADIATION
                       rflux2, rfd2_lo, rfd2_hi, &
#endif
                       qaux, qa_lo, qa_hi, &
                       shk, shk_lo, shk_hi, &
                       2, lo(1), hi(1), lo(2), hi(2)+1, kc, k3d, k3d, domlo, domhi)

           do j = lo(2)-1, hi(2)+2
              do i = lo(1)-1, hi(1)+1
                 q2(i,j,k3d,:) = qint(i,j,kc,:)
              enddo
           enddo
           
           ! Compute F^z at kc (k3d)
           call cmpflx(qzm, qzp, It_lo, It_hi, &
                       flux3, flux3_lo, flux3_hi, &
                       qint, It_lo, It_hi, &
#ifdef RADIATION
                       rflux3, rfd3_lo, rfd3_hi, &
#endif
                       qaux, qa_lo, qa_hi, &
                       shk, shk_lo, shk_hi, &
                       3, lo(1), hi(1), lo(2), hi(2), kc, k3d, k3d, domlo, domhi)

           do j=lo(2)-1,hi(2)+1
              do i=lo(1)-1,hi(1)+1
                 q3(i,j,k3d,:) = qint(i,j,kc,:)
              enddo
           enddo

        endif

     enddo
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

              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

           enddo
        enddo
     enddo
  enddo

#ifdef HYBRID_MOMENTUM
  call add_hybrid_advection_source(lo, hi, dt, &
                                   update, uout_lo, uout_hi, &
                                   qx, qx_lo, qx_hi, &
                                   qy, qy_lo, qy_hi, &
                                   qz, qz_lo, qz_hi)
#endif

  call bl_deallocate(   div)
  call bl_deallocate( pdivu)

  call bl_deallocate(    q1)
  call bl_deallocate(    q2)
  call bl_deallocate(    q3)

end subroutine ca_mol_single_stage
