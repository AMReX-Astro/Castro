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
                               flux, flux_lo, flux_hi, &
                               pradial, p_lo, p_hi, &
                               area, area_lo, area_hi, &
                               dloga, dloga_lo, dloga_hi, &
                               vol, vol_lo, vol_hi, &
                               courno, verbose) bind(C, name="ca_mol_single_stage")

  use meth_params_module, only : NQ, QVAR, QU, QPRES, &
                                 UMX, UMY, UMZ, UTEMP, USHK, UEINT, &
                                 NQAUX, NVAR, NHYP, use_flattening, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 NGDNV, GDU, GDPRES, first_order_hydro, difmag, &
                                 hybrid_riemann, ppm_temp_fix
  use advection_util_module, only : compute_cfl, shock, normalize_species_fluxes
  use bl_constants_module, only : ZERO, HALF, ONE
  use flatten_module, only : uflatten
  use prob_params_module, only : coord_type
  use riemann_module, only: cmpflx
  use ppm_module, only : ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec, naux

  implicit none

  integer, intent(in) :: lo(1), hi(1), verbose
  integer, intent(in) :: domlo(1), domhi(1)
  real(rt), intent(in) :: stage_weight
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flux_lo(3), flux_hi(3)
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: area_lo(3), area_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt)        , intent(in) ::      uin(  uin_lo(1):  uin_hi(1),NVAR)
  real(rt)        , intent(inout) ::  uout( uout_lo(1): uout_hi(1),NVAR)
  real(rt)        , intent(inout) ::     q(    q_lo(1):    q_hi(1),NQ)
  real(rt)        , intent(in) ::     qaux(   qa_lo(1):   qa_hi(1),NQAUX)
  real(rt)        , intent(in) ::     srcU(  srU_lo(1):  srU_hi(1),NVAR)
  real(rt)        , intent(inout) :: update(updt_lo(1): updt_hi(1),NVAR)
  real(rt)        , intent(inout) :: update_flux(uf_lo(1): uf_hi(1),NVAR)
  real(rt)        , intent(inout) ::  flux( flux_lo(1): flux_hi(1),NVAR)
  real(rt)        , intent(inout) :: pradial(  p_lo(1):   p_hi(1))
  real(rt)        , intent(in) :: area( area_lo(1): area_hi(1)     )
  real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1)     )
  real(rt)        , intent(in) ::   vol(  vol_lo(1): vol_hi(1)      )
  real(rt)        , intent(in) :: delta(1), dt, time
  real(rt)        , intent(inout) :: courno

  ! Automatic arrays for workspace
  real(rt)        , allocatable :: flatn(:)
  real(rt)        , allocatable :: div(:)
  real(rt)        , allocatable :: pdivu(:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , allocatable :: q1(:,:)

  ! radiation fluxes (needthese to get things to compile)
  real(rt)        , allocatable :: rflx(:,:)

  real(rt)        , allocatable :: shk(:)

  ! temporary interface values of the parabola
  real(rt), allocatable :: sxm(:), sxp(:)

  real(rt), allocatable :: qxm(:,:), qxp(:,:)

  real(rt)         :: dx

  integer i, n, ngf

  integer :: lo_3D(3), hi_3D(3)
  real(rt)         :: dx_3D(3)
  integer :: qp_lo(3), qp_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: div1

  type(eos_t) :: eos_state
  
  ngf = 1

  lo_3D   = [lo(1), 0, 0]
  hi_3D   = [hi(1), 0, 0]

  dx_3D   = [delta(1), ZERO, ZERO]

  shk_lo = [lo(1)-1, 0, 0]
  shk_hi = [hi(1)+1, 0, 0]

  qp_lo = [lo(1)-1, 0, 0]
  qp_hi = [hi(1)+2, 0, 0]


  allocate( flatn(q_lo(1):q_hi(1)))

  allocate(    q1(flux_lo(1):flux_hi(1), NGDNV))

  ! when we do radiation, these would be passed out
  allocate(rflx(flux_lo(1):flux_hi(1), NGDNV))

  allocate( pdivu(lo(1):hi(1)  ))
  allocate( shk(shk_lo(1):shk_hi(1)))

  dx = delta(1)


#ifdef SHOCK_VAR
    uout(lo(1):hi(1),USHK) = ZERO

    call shock(q, q_lo, q_hi, shk, shk_lo, shk_hi, lo_3D, hi_3D, dx_3D)

    ! Store the shock data for future use in the burning step.
    do i = lo(1), hi(1)
       uout(i,USHK) = shk(i)
    enddo

    ! Discard it locally if we don't need it in the hydro update.
    if (hybrid_riemann /= 1) then
       shk(:) = ZERO
    endif
#else
    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, q_lo, q_hi, shk, shk_lo, shk_hi, lo_3D, hi_3D, dx_3D)
    else
       shk(:) = ZERO
    endif
#endif


  ! Check if we have violated the CFL criterion.
  call compute_cfl(q, q_lo, q_hi, &
                   qaux, qa_lo, qa_hi, &
                   lo_3D, hi_3D, dt, dx_3D, courno)

  ! Compute flattening coefficient for slope calculations.
  if (first_order_hydro == 1) then
     flatn = ZERO

  else if (use_flattening == 1) then
     call uflatten([lo(1) - ngf, 0, 0], [hi(1) + ngf, 0, 0], &
                   q, flatn, q_lo, q_hi, QPRES)
  else
     flatn = ONE
  endif


  ! sm and sp are the minus and plus parts of the parabola -- they are
  ! defined for a single zone, so for zone i, sm is the left value of
  ! the parabola and sp is the right value of the parabola
  allocate(sxm(q_lo(1):q_hi(1)))
  allocate(sxp(q_lo(1):q_hi(1)))

  ! qm and qp are the left and right states for an interface -- they
  ! are defined for a particular interface, with the convention that
  ! qm(i) and qp(i) correspond to the i-1/2 interface
  allocate ( qxm(qp_lo(1):qp_hi(1),NQ) )
  allocate ( qxp(qp_lo(1):qp_hi(1),NQ) )

  ! Do PPM reconstruction
  do n = 1, QVAR
     call ppm_reconstruct(q(:,n), q_lo, q_hi, &
                          flatn, q_lo, q_hi, &
                          sxm, sxp, sxm, sxp, sxm, sxp, q_lo, q_hi, &  ! extras are dummy
                          lo(1), 0, hi(1), 0, [dx, ZERO, ZERO], 0, 0)

     ! Construct the interface states -- this is essentially just a
     ! reshuffling of interface states from zone-center indexing to
     ! edge-centered indexing
     do i = lo(1)-1, hi(1)+1
        qxm(i,n) = sxp(i-1)
        qxp(i,n) = sxm(i)
     enddo
  enddo

  deallocate(sxm, sxp)

  ! use T to define p
  if (ppm_temp_fix == 1) then
     do i = lo(1)-1, hi(1)+1

        eos_state%rho    = qxp(i,QRHO)
        eos_state%T      = qxp(i,QTEMP)
        eos_state%xn(:)  = qxp(i,QFS:QFS-1+nspec)
        eos_state%aux(:) = qxp(i,QFX:QFX-1+naux)

        call eos(eos_input_rt, eos_state)

        qxp(i,QPRES) = eos_state%p
        qxp(i,QREINT) = qxp(i,QRHO)*eos_state%e
        ! should we try to do something about Gamma_! on interface?

        eos_state%rho    = qxm(i,QRHO)
        eos_state%T      = qxm(i,QTEMP)
        eos_state%xn(:)  = qxm(i,QFS:QFS-1+nspec)
        eos_state%aux(:) = qxm(i,QFX:QFX-1+naux)
        
        call eos(eos_input_rt, eos_state)

        qxm(i,QPRES) = eos_state%p
        qxm(i,QREINT) = qxm(i,QRHO)*eos_state%e
        ! should we try to do something about Gamma_! on interface?

     enddo
  endif

  ! Get the fluxes from the Riemann solver
  call cmpflx(lo, hi, domlo, domhi, &
              qxm, qxp, qp_lo, qp_hi, &
              flux, flux_lo, flux_hi, &
              q1, flux_lo, flux_hi, &
#ifdef RADIATION
              rflx, flux_lo, flux_hi, &
#endif
              qaux, qa_lo, qa_hi, lo(1), hi(1))

  deallocate(qxm, qxp)

  ! construct p div(U)
  do i = lo(1), hi(1)
     pdivu(i) = HALF * &
          (q1(i+1,GDPRES) + q1(i,GDPRES))* &
          (q1(i+1,GDU)*area(i+1) - q1(i,GDU)*area(i)) / vol(i)
  end do

  ! Compute the artifical viscosity
  allocate(div(lo(1):hi(1)+1))

  do i = lo(1), hi(1)+1
     div(i) = (q(i,QU) - q(i-1,QU)) / dx
  enddo

  do n = 1, NVAR
     if ( n == UTEMP ) then
        flux(lo(1):hi(1)+1,n) = ZERO
#ifdef SHOCK_VAR
     else if ( n == USHK) then
        flux(lo(1):hi(1)+1,n) = ZERO
#endif
     else if ( n == UMY ) then
        flux(lo(1):hi(1)+1,n) = ZERO
     else if ( n == UMZ ) then
        flux(lo(1):hi(1)+1,n) = ZERO
     else
        ! add the artifical viscosity
        do i = lo(1),hi(1)+1
           div1 = difmag*min(ZERO,div(i))
           flux(i,n) = flux(i,n) + dx*div1*(uin(i,n) - uin(i-1,n))
        enddo
     endif
  enddo

  ! Normalize the species fluxes
  call normalize_species_fluxes(flux, flux_lo(1), flux_hi(1), lo_3D, hi_3D)


  ! Make the update for this state

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update

  do n = 1, NVAR
     do i = lo(1), hi(1)
        update(i,n) = update(i,n) + ( flux(i,n) * area(i) - &
                                      flux(i+1,n) * area(i+1) ) / vol(i)

        ! Add p div(u) source term to (rho e)
        if (n == UEINT) then
           update(i,n) = update(i,n) - pdivu(i)
        else if (n == UMX) then
           update(i,UMX) = update(i,UMX) - ( q1(i+1,GDPRES) - q1(i,GDPRES) ) / dx
        endif

        ! for storage
        update_flux(i,n) = update_flux(i,n) + stage_weight * update(i,n)

        ! include source terms
        update(i,n) = update(i,n) + srcU(i,n)

     enddo
  enddo

  ! Scale the fluxes for the form we expect later in refluxing.

  do n = 1, NVAR
     do i = lo(1), hi(1)+1

        flux(i,n) = dt * area(i) * flux(i,n)

        ! Correct the momentum flux with the grad p part.
        if (coord_type .eq. 0 .and. n == UMX) then
           flux(i,n) = flux(i,n) + dt * area(i) * q1(i,GDPRES)
        endif

     enddo
  enddo


  if (coord_type .gt. 0) then
     pradial(lo(1):hi(1)+1) = q1(lo(1):hi(1)+1,GDPRES) * dt
  end if

  deallocate(flatn,div,pdivu,q1)

end subroutine ca_mol_single_stage
