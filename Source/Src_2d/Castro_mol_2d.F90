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
                               delta, dt, &
                               flux1, flux1_lo, flux1_hi, &
                               flux2, flux2_lo, flux2_hi, &
                               pradial, p_lo, p_hi, &
                               area1, area1_lo, area1_hi, &
                               area2, area2_lo, area2_hi, &
                               dloga, dloga_lo, dloga_hi, &
                               vol, vol_lo, vol_hi, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                 UTEMP, UEINT, USHK, UMX, GDU, GDV, &
                                 use_flattening, QU, QV, QW, QPRES, NQAUX, &
                                 first_order_hydro, difmag, hybrid_riemann
  use advection_util_2d_module, only : divu, normalize_species_fluxes
  use advection_util_module, only : compute_cfl
  use bl_constants_module, only : ZERO, HALF, ONE
  use flatten_module, only : uflatten
  use prob_params_module, only : coord_type
  use riemann_module, only: cmpflx, shock
  use ppm_module, only : ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
#ifdef RADIATION  
  use rad_params_module, only : ngroups
#endif
  implicit none

  integer, intent(in) :: lo(2), hi(2), verbose
  integer, intent(in) :: domlo(2), domhi(2)
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
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),NVAR)
  real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),NVAR)
  real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),NQ)
  real(rt)        , intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),NQAUX)
  real(rt)        , intent(in) :: srcU(srU_lo(1):srU_hi(1),srU_lo(2):srU_hi(2),NVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),NVAR)
  real(rt)        , intent(inout) :: update_flux(uf_lo(1):uf_hi(1),uf_lo(2):uf_hi(2),NVAR)
  real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),NVAR)
  real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),NVAR)
  real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2))
  real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2))
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2))
  real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2))
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

  integer :: lo_3D(3), hi_3D(3)
  real(rt) :: dx_3D(3)

  real(rt) :: div1

  integer :: i, j, n

  ngf = 1

  lo_3D   = [lo(1), lo(2), 0]
  hi_3D   = [hi(1), hi(2), 0]

  dx_3D   = [delta(1), delta(2), ZERO]

  allocate(flatn(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))

  allocate(q1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), NGDNV))
  allocate(q2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), NGDNV))

#ifdef RADIATION
  ! when we do radiation, these would be passed out
  allocate(rflx(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), ngroups))
  allocate(rfly(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), ngroups))
#endif

  allocate(shk(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1))
  allocate(pdivu(lo(1):hi(1), lo(2):hi(2)))

  dx = delta(1)
  dy = delta(2)


#ifdef SHOCK_VAR
    uout(lo(1):hi(1),lo(2):hi(2),USHK) = ZERO

    call shock(q, q_lo(1), q_lo(2), q_hi(1), q_hi(2), &
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
       call shock(q, q_lo(1), q_lo(2), q_hi(1), q_hi(2), &
                  shk, lo(1)-1, lo(2)-1, hi(1)+1, hi(2)+1, &
                  lo(1), lo(2), hi(1), hi(2), dx, dy)
    else
       shk(:,:) = ZERO
    endif
#endif


  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  elseif (use_flattening == 1) then
     call uflatten([lo(1) - ngf, lo(2) - ngf, 0], [hi(1) + ngf, hi(2) + ngf, 0], &
                   q, flatn, q_lo, q_hi)
  else
     flatn = ONE
  endif


  ! sm and sp are the minus and plus parts of the parabola -- they are
  ! defined for a single zone, so for zone i, sm is the left value of
  ! the parabola and sp is the right value of the parabola
  allocate(sxm(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
  allocate(sxp(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
  allocate(sym(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))
  allocate(syp(q_lo(1):q_hi(1), q_lo(2):q_hi(2)))

  ! qm and qp are the left and right states for an interface -- they
  ! are defined for a particular interface, with the convention that
  ! qm(i) and qp(i) correspond to the i-1/2 interface
  allocate ( qxm(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )
  allocate ( qxp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )
  allocate ( qym(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )
  allocate ( qyp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NQ) )

  ! Do PPM reconstruction
  do n = 1, QVAR
     call ppm_reconstruct(q(:,:,n), q_lo(1), q_lo(2), q_hi(1), q_hi(2), &
                          flatn, q_lo(1), q_lo(2), q_hi(1), q_hi(2), &
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
              flux1, flux1_lo(1), flux1_lo(2), flux1_hi(1), flux1_hi(2), &
              q1, flux1_lo(1), flux1_lo(2), flux1_hi(1), flux1_hi(2), &
#ifdef RADIATION
              rflx, flux1_lo(1), flux1_lo(2), flux1_hi(1), flux1_hi(2), &
#endif
              qaux, qa_lo(1), qa_lo(2), qa_hi(1), qa_hi(2), &
              shk, lo(1)-1, lo(2)-1, hi(1)+1, hi(2)+1, &
              1, lo(1), hi(1), lo(2), hi(2), domlo, domhi)


  call cmpflx(qym, qyp, lo(1)-1, lo(2)-1, hi(1)+2, hi(2)+2, &
              flux2, flux2_lo(1), flux2_lo(2), flux2_hi(1), flux2_hi(2), &
              q2, flux2_lo(1), flux2_lo(2), flux2_hi(1), flux2_hi(2), &
#ifdef RADIATION
              rfly, flux2_lo(1), flux2_lo(2), flux2_hi(1), flux2_hi(2), &
#endif
              qaux, qa_lo(1), qa_lo(2), qa_hi(1), qa_hi(2), &
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

  call divu(lo, hi, q, q_lo(1), q_lo(2), q_hi(1), q_hi(2), &
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
  call normalize_species_fluxes(flux1,flux1_lo(1),flux1_lo(2),flux1_hi(1),flux1_hi(2), &
                                flux2,flux2_lo(1),flux2_lo(2),flux2_hi(1),flux2_hi(2), &
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

           ! for storage
           update_flux(i,j,n) = update_flux(i,j,n) + stage_weight * update(i,j,n)

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
