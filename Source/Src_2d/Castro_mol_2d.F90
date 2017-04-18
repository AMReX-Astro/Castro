! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(time, &
                               lo, hi, domlo, domhi, &
                               uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                               uout, uout_l1, uout_l2, uout_h1, uout_h2, &
                               q, qd_l1, qd_l2, qd_h1, qd_h2, &
                               qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
                               srcU, srU_l1, srU_l2, srU_h1, srU_h2, &
                               update, updt_l1, updt_l2, updt_h1, updt_h2, &
                               delta, dt, &
                               flux1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
                               flux2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
                               pradial, p_l1, p_l2, p_h1, p_h2, &
                               area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                               area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                               dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                               vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                 UTEMP, UEINT, USHK, UMX, GDU, GDV, &
                                 use_flattening, QU, QV, QW, QPRES, NQAUX, &
                                 first_order_hydro, difmag, hybrid_riemann
  use advection_util_2d_module, only : divu, normalize_species_fluxes
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, HALF, ONE
  use flatten_module, only : uflaten
  use prob_params_module, only : coord_type
  use riemann_module, only: cmpflx, shock
  use ppm_module, only : ppm_reconstruct
  use bl_fort_module, only : rt => c_real
#ifdef RADIATION  
  use rad_params_module, only : ngroups
#endif
  implicit none

  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
  integer, intent(in) :: uin_l1, uin_l2, uin_h1, uin_h2
  integer, intent(in) :: uout_l1,uout_l2,uout_h1,uout_h2
  integer, intent(in) :: qd_l1, qd_l2, qd_h1, qd_h2
  integer, intent(in) :: qa_l1, qa_l2, qa_h1, qa_h2
  integer, intent(in) :: srU_l1, srU_l2, srU_h1, srU_h2
  integer, intent(in) :: updt_l1,updt_l2,updt_h1,updt_h2
  integer, intent(in) :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
  integer, intent(in) :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
  integer, intent(in) :: p_l1,p_l2,p_h1,p_h2
  integer, intent(in) :: area1_l1,area1_l2,area1_h1,area1_h2
  integer, intent(in) :: area2_l1,area2_l2,area2_h1,area2_h2
  integer, intent(in) :: dloga_l1,dloga_l2,dloga_h1,dloga_h2
  integer, intent(in) :: vol_l1,vol_l2,vol_h1,vol_h2

  real(rt)        , intent(in) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  real(rt)        , intent(inout) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
  real(rt)        , intent(inout) :: q(qd_l1:qd_h1,qd_l2:qd_h2,NQ)
  real(rt)        , intent(inout) :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  real(rt)        , intent(in) :: srcU(srU_l1:srU_h1,srU_l2:srU_h2,NVAR)
  real(rt)        , intent(inout) :: update(updt_l1:updt_h1,updt_l2:updt_h2,NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
  real(rt)        , intent(inout) :: pradial(p_l1:p_h1,p_l2:p_h2)
  real(rt)        , intent(in) :: area1(area1_l1:area1_h1,area1_l2:area1_h2)
  real(rt)        , intent(in) :: area2(area2_l1:area2_h1,area2_l2:area2_h2)
  real(rt)        , intent(in) :: dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  real(rt)        , intent(in) :: vol(vol_l1:vol_h1,vol_l2:vol_h2)
  real(rt)        , intent(in) :: delta(2), dt, time
  real(rt)        , intent(inout) :: courno

  ! Automatic arrays for workspace
  real(rt)        , allocatable :: flatn(:,:)
  real(rt)        , allocatable :: div(:,:)
  real(rt)        , allocatable :: pdivu(:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , allocatable :: q1(:,:,:)
  real(rt)        , allocatable :: q2(:,:,:)

  ! radiation fluxes (need these to get things to compile)
  real(rt)        , allocatable :: rflx(:,:,:)
  real(rt)        , allocatable :: rfly(:,:,:)

  real(rt)        , allocatable :: shk(:,:)

  ! temporary interface values of the parabola
  real(rt), allocatable :: sxm(:,:), sxp(:,:), sym(:,:), syp(:,:)

  real(rt)        , allocatable:: qxm(:,:,:), qym(:,:,:)
  real(rt)        , allocatable:: qxp(:,:,:), qyp(:,:,:)

  integer ngf
  real(rt) :: dx, dy

  integer ::  uin_lo(2),  uin_hi(2)
  integer :: uout_lo(2), uout_hi(2)
  integer :: q_lo(2), q_hi(2)

  integer :: lo_3D(3), hi_3D(3)
  integer :: q_lo_3D(3), q_hi_3D(3)
  integer :: uin_lo_3D(3), uin_hi_3D(3)
  real(rt) :: dx_3D(3)

  real(rt) :: div1

  integer :: i, j, n

  uin_lo  = [uin_l1, uin_l2]
  uin_hi  = [uin_h1, uin_h2]

  uout_lo = [uout_l1, uout_l2]
  uout_hi = [uout_h1, uout_h2]
  
  ngf = 1

  q_lo = [qd_l1, qd_l2]
  q_hi = [qd_h1, qd_h2]

  lo_3D   = [lo(1), lo(2), 0]
  hi_3D   = [hi(1), hi(2), 0]

  q_lo_3D = [qd_l1, qd_l2, 0]
  q_hi_3D = [qd_h1, qd_h2, 0]

  uin_lo_3D = [uin_l1, uin_l2, 0]
  uin_hi_3D = [uin_h1, uin_h2, 0]

  dx_3D   = [delta(1), delta(2), ZERO]

  allocate(flatn(qd_l1:qd_h1, qd_l2:qd_h2))

  allocate(q1(flux1_l1:flux1_h1, flux1_l2:flux1_h2, NGDNV))
  allocate(q2(flux2_l1:flux2_h1, flux2_l2:flux2_h2, NGDNV))

#ifdef RADIATION
  ! when we do radiation, these would be passed out
  allocate(rflx(flux1_l1:flux1_h1, flux1_l2:flux1_h2, ngroups))
  allocate(rfly(flux2_l1:flux2_h1, flux2_l2:flux2_h2, ngroups))
#endif

  allocate(shk(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1))
  allocate(pdivu(lo(1):hi(1), lo(2):hi(2)))

  dx = delta(1)
  dy = delta(2)


#ifdef SHOCK_VAR
    uout(lo(1):hi(1),lo(2):hi(2),USHK) = ZERO

    call shock(q, qd_l1, qd_l2, qd_h1, qd_h2, &
               shk, lo(1)-1, lo(2)-1, hi(1)+1, hi(2)+1, &
               lo(1), lo(2), hi(1), hi(2), dx, dy)

    ! Store the shock data for future use in the burning step.
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          uout(i,j,USHK) = shk(i,j)
       enddo
    enddo

    ! Discard it locally if we don't need it in the hydro update.

    if (hybrid_riemann /= 1) then
       shk(:,:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, qd_l1, qd_l2, qd_h1, qd_h2, &
                  shk, lo(1)-1, lo(2)-1, hi(1)+1, hi(2)+1, &
                  lo(1), lo(2), hi(1), hi(2), dx, dy)
    else
       shk(:,:) = ZERO
    endif
#endif


  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo_3D, q_hi_3D, &
                   qaux, [qa_l1, qa_l2, 0], [qa_h1, qa_h2, 0], &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  elseif (use_flattening == 1) then
     call uflaten([lo(1) - ngf, lo(2) - ngf, 0], [hi(1) + ngf, hi(2) + ngf, 0], &
                  q(:,:,QPRES), q(:,:,QU), q(:,:,QV), q(:,:,QW), &
                  flatn, [qd_l1, qd_l2, 0], [qd_h1, qd_h2, 0])
  else
     flatn = ONE
  endif


  ! sm and sp are the minus and plus parts of the parabola -- they are
  ! defined for a single zone, so for zone i, sm is the left value of
  ! the parabola and sp is the right value of the parabola
  allocate(sxm(qd_l1:qd_h1, qd_l2:qd_h2))
  allocate(sxp(qd_l1:qd_h1, qd_l2:qd_h2))
  allocate(sym(qd_l1:qd_h1, qd_l2:qd_h2))
  allocate(syp(qd_l1:qd_h1, qd_l2:qd_h2))

  ! qm and qp are the left and right states for an interface -- they
  ! are defined for a particular interface, with the convention that
  ! qm(i) and qp(i) correspond to the i-1/2 interface
  allocate ( qxm(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )
  allocate ( qxp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )
  allocate ( qym(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )
  allocate ( qyp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )

  ! Do PPM reconstruction
  do n = 1, QVAR
     call ppm_reconstruct(q(:,:,n), qd_l1, qd_l2, qd_h1, qd_h2, &
                          flatn, qd_l1, qd_l2, qd_h1, qd_h2, &
                          sxm, sxp, sym, syp, &
                          lo(1), lo(2), hi(1), hi(2), dx, dy)

     ! Construct the interface states -- this is essentially just a
     ! reshuffling of interface states from zone-center indexing to
     ! edge-centered indexing
     do j = lo(2)-1, hi(2)+1
        do i = lo(1)-1, hi(1)+1

           ! x-edges

           ! left state at i-1/2 interface
           qxm(i,j,n) = sxp(i-1,j)

           ! right state at i-1/2 interface
           qxp(i,j,n) = sxm(i,j)

           ! y-edges

           ! left state at j-1/2 interface
           qym(i,j,n) = syp(i,j-1)

           ! right state at j-1/2 interface
           qyp(i,j,n) = sym(i,j)

        enddo
     enddo
  enddo

  deallocate(sxm, sxp, sym, syp)


  ! Get the fluxes from the Riemann solver
  call cmpflx(qxm, qxp, lo(1)-1, lo(2)-1, hi(1)+2, hi(2)+2, &
              flux1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
              q1, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
#ifdef RADIATION
              rflx, flux1_l1, flux1_l2, flux1_h1, flux1_h2, &
#endif
              qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
              shk, lo(1)-1, lo(2)-1, hi(1)+1, hi(2)+1, &
              1, lo(1), hi(1), lo(2), hi(2), domlo, domhi)


  call cmpflx(qym, qyp, lo(1)-1, lo(2)-1, hi(1)+2, hi(2)+2, &
              flux2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
              q2, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
#ifdef RADIATION
              rfly, flux2_l1, flux2_l2, flux2_h1, flux2_h2, &
#endif
              qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
              shk, lo(1)-1, lo(2)-1, hi(1)+1, hi(2)+1, &
              2, lo(1), hi(1), lo(2), hi(2), domlo, domhi)

  deallocate(qxm, qxp, qym, qyp)


  ! construct p div{U}
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        pdivu(i,j) = HALF*( &
             (q1(i+1,j,GDPRES) + q1(i,j,GDPRES)) * &
             (q1(i+1,j,GDU)*area1(i+1,j) - q1(i,j,GDU)*area1(i,j)) + &
             (q2(i,j+1,GDPRES) + q2(i,j,GDPRES)) * &
             (q2(i,j+1,GDV)*area2(i,j+1) - q2(i,j,GDV)*area2(i,j)) ) / vol(i,j)
     enddo
  enddo

  ! Compute the artifical viscosity

  ! Compute divergence of velocity field (on surrounding nodes(lo,hi))
  ! this is used for the artifical viscosity
  allocate(div(lo(1):hi(1)+1 ,lo(2):hi(2)+1))

  call divu(lo, hi, q, qd_l1, qd_l2, qd_h1, qd_h2, &
            delta, div, lo(1), lo(2), hi(1)+1, hi(2)+1)

  do n = 1, NVAR
     if (n == UTEMP) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),n) = ZERO
        flux2(lo(1):hi(1),lo(2):hi(2)+1,n) = ZERO
#ifdef SHOCK_VAR
     else if (n == USHK) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),n) = ZERO
        flux2(lo(1):hi(1),lo(2):hi(2)+1,n) = ZERO
#endif
     else
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              div1 = HALF*(div(i,j) + div(i,j+1))
              div1 = difmag*min(ZERO, div1)

              flux1(i,j,n) = flux1(i,j,n) + &
                   dx*div1*(uin(i,j,n) - uin(i-1,j,n))
           enddo
        enddo

        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              div1 = HALF*(div(i,j) + div(i+1,j))
              div1 = difmag*min(ZERO,div1)

              flux2(i,j,n) = flux2(i,j,n) + &
                   dy*div1*(uin(i,j,n) - uin(i,j-1,n))
           enddo
        enddo

     endif
  enddo


  ! Normalize the species fluxes
  call normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                                flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                                lo,hi)


  ! Make the update for this state

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update
  
  do n = 1, NVAR
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           update(i,j,n) = update(i,j,n) + &
                ( flux1(i,j,n) * area1(i,j) - flux1(i+1,j,n) * area1(i+1,j) + &
                  flux2(i,j,n) * area2(i,j) - flux2(i,j+1,n) * area2(i,j+1) ) / vol(i,j)

           if (n == UEINT) then
              ! Add p div(u) source term to (rho e)
              update(i,j,n) = update(i,j,n) - pdivu(i,j)

           else if (n == UMX) then
              ! add the pressure source term for axisummetry 
              if (coord_type == 1) then
                 update(i,j,n) = update(i,j,n) - (q1(i+1,j,GDPRES) - q1(i,j,GDPRES))/ dx
              endif
           endif

           ! include source terms
           update(i,j,n) = update(i,j,n) + srcU(i,j,n)

        enddo
     enddo
  enddo

  ! Scale the fluxes for the form we expect later in refluxing.
  do n = 1, NVAR
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
           flux1(i,j,n) = dt * flux1(i,j,n) * area1(i,j)
        enddo
     enddo
  enddo
  
  do n = 1, NVAR
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           flux2(i,j,n) = dt * flux2(i,j,n) * area2(i,j)
        enddo
     enddo
  enddo

  ! Store fluxes for flux correction
  if (coord_type .eq. 1) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2)) = q1(lo(1):hi(1)+1,lo(2):hi(2),GDPRES) * dt
  end if

  deallocate(flatn,div,q1,q2,pdivu)

#ifdef RADIATION
  deallocate(rflx, rfly)
#endif

end subroutine ca_mol_single_stage
