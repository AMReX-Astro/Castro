!> @brief advection routines in support of method of lines integration
!!
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
#if AMREX_SPACEDIM >= 2
     flux2, flux2_lo, flux2_hi, &
#endif
#if AMREX_SPACEDIM == 3
     flux3, flux3_lo, flux3_hi, &
#endif
     area1, area1_lo, area1_hi, &
#if AMREX_SPACEDIM >= 2
     area2, area2_lo, area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
     area3, area3_lo, area3_hi, &
#endif
#if AMREX_SPACEDIM < 3
     pradial, p_lo, p_hi, &
     dloga, dloga_lo, dloga_hi, &
#endif
     vol, vol_lo, vol_hi, &
     verbose) bind(C, name="ca_mol_single_stage")

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
       UTEMP, USHK, UMX, &
       use_flattening, QPRES, NQAUX, &
       QTEMP, QFS, QFX, QREINT, QRHO, &
       first_order_hydro, difmag, hybrid_riemann, &
       limit_fluxes_on_small_dens, ppm_type, ppm_temp_fix
  use advection_util_module, only : limit_hydro_fluxes_on_small_dens, ca_shock, &
       divu, normalize_species_fluxes, calc_pdivu, &
       scale_flux, apply_av
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: ca_uflatten
  use riemann_module, only: cmpflx

  use riemann_util_module, only : ca_store_godunov_state
  use ppm_module, only : ca_ppm_reconstruct
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
#if AMREX_SPACEDIM >= 2
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
  integer, intent(in) :: flux3_lo(3), flux3_hi(3)
  integer, intent(in) :: area3_lo(3), area3_hi(3)
#endif
#if AMREX_SPACEDIM <= 2
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
#if AMREX_SPACEDIM >= 2
  real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if AMREX_SPACEDIM == 3
  real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
  real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
#if AMREX_SPACEDIM <= 2
  real(rt), intent(inout) :: pradial(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
  real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1), dloga_lo(2):dloga_hi(2), dloga_lo(3):dloga_hi(3))
#endif
  real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt), intent(in) :: dx(3), dt, time

  ! Automatic arrays for workspace
  real(rt)        , pointer:: flatn(:,:,:)
  real(rt)        , pointer:: div(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , pointer:: q_int(:,:,:,:)
#ifdef RADIATION
  real(rt)        , pointer:: lambda_int(:,:,:,:)
#endif
  real(rt)        , pointer:: q1(:,:,:,:)
  real(rt)        , pointer:: q2(:,:,:,:)
  real(rt)        , pointer:: q3(:,:,:,:)

#ifdef RADIATION
  ! radiation fluxes (need these to get things to compile)
  real(rt)        , pointer:: rflx(:,:,:,:)
  real(rt)        , pointer:: rfly(:,:,:,:)
  real(rt)        , pointer:: rflz(:,:,:,:)
#endif

  real(rt)        , pointer:: shk(:,:,:)

  ! temporary interface values of the parabola
  real(rt)        , pointer :: qm(:,:,:,:,:), qp(:,:,:,:,:)

  integer :: ngf
  integer :: It_lo(3), It_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  integer :: idir, i, j, k, n

  type (eos_t) :: eos_state

  ngf = 1

  It_lo = lo(:) - dg(:)
  It_hi = hi(:) + 2*dg(:)

  shk_lo(:) = lo(:) - dg(:)
  shk_hi(:) = hi(:) + dg(:)

  call bl_allocate(   div, lo, hi+dg)

  call bl_allocate(q_int, It_lo, It_hi, NQ)
#ifdef RADIATION
  call bl_allocate(lambda_int, It_lo(1), It_hi(1), It_lo(2), It_hi(2), It_lo(3), It_hi(3), 0, ngroups-1)
#endif

  call bl_allocate(q1, flux1_lo, flux1_hi, NGDNV)
#if AMREX_SPACEDIM >= 2
  call bl_allocate(q2, flux2_lo, flux2_hi, NGDNV)
#endif
#if AMREX_SPACEDIM == 3
  call bl_allocate(q3, flux3_lo, flux3_hi, NGDNV)
#endif

#ifdef RADIATION
  ! when we do radiation, these would be passed out
  call bl_allocate(rflx, flux1_lo, flux1_hi, ngroups)
#if AMREX_SPACEDIM >= 2
  call bl_allocate(rfly, flux2_lo, flux2_hi, ngroups)
#endif
#if AMREX_SPACEDIM == 3
  call bl_allocate(rflz, flux3_lo, flux3_hi, ngroups)
#endif
#endif

  call bl_allocate(qm, It_lo(1),It_hi(1), It_lo(2),It_hi(2), It_lo(3),It_hi(3), 1,NQ, 1,AMREX_SPACEDIM)
  call bl_allocate(qp, It_lo(1),It_hi(1), It_lo(2),It_hi(2), It_lo(3),It_hi(3), 1,NQ, 1,AMREX_SPACEDIM)


  call bl_allocate(shk, shk_lo, shk_hi)

#ifndef AMREX_USE_CUDA
  if (ppm_type == 0) then
     call amrex_error("ERROR: method of lines integration does not support ppm_type = 0")
  endif
#endif

#ifdef SHOCK_VAR
  uout(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), USHK) = ZERO

  call ca_shock(lo-dg, hi+dg, &
                q, q_lo, q_hi, &
                shk, shk_lo, shk_hi, &
                dx)

  ! Store the shock data for future use in the burning step.

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uout(i,j,k,USHK) = shk(i,j,k)
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
     call ca_shock(lo-dg, hi+dg, &
                   q, q_lo, q_hi, &
                   shk, shk_lo, shk_hi, &
                   dx)
  else
     shk(:,:,:) = ZERO
  endif
#endif

  ! Compute flattening coefficient for slope calculations.
  call bl_allocate(flatn, q_lo, q_hi)

  if (first_order_hydro == 1) then
     flatn = ZERO
  elseif (use_flattening == 1) then
     call ca_uflatten(lo - ngf*dg, hi + ngf*dg, &
          q, q_lo, q_hi, &
          flatn, q_lo, q_hi, QPRES)
  else
     flatn = ONE
  endif


  call ca_ppm_reconstruct(lo-dg, hi+dg, 1, &
       q, q_lo, q_hi, NQ, 1, NQ, &
       flatn, q_lo, q_hi, &
       qm, It_lo, It_hi, &
       qp, It_lo, It_hi, NQ, 1, NQ)

  ! use T to define p
  if (ppm_temp_fix == 1) then
     do idir = 1, AMREX_SPACEDIM
        do k = lo(3)-dg(3), hi(3)+dg(3)
           do j = lo(2)-dg(2), hi(2)+dg(2)
              do i = lo(1)-1, hi(1)+1

                 eos_state%rho    = qp(i,j,k,QRHO,idir)
                 eos_state%T      = qp(i,j,k,QTEMP,idir)
                 eos_state%xn(:)  = qp(i,j,k,QFS:QFS-1+nspec,idir)
                 eos_state%aux(:) = qp(i,j,k,QFX:QFX-1+naux,idir)

                 call eos(eos_input_rt, eos_state)

                 qp(i,j,k,QPRES,idir) = eos_state%p
                 qp(i,j,k,QREINT,idir) = qp(i,j,k,QRHO,idir)*eos_state%e
                 ! should we try to do something about Gamma_! on interface?

                 eos_state%rho    = qm(i,j,k,QRHO,idir)
                 eos_state%T      = qm(i,j,k,QTEMP,idir)
                 eos_state%xn(:)  = qm(i,j,k,QFS:QFS-1+nspec,idir)
                 eos_state%aux(:) = qm(i,j,k,QFX:QFX-1+naux,idir)

                 call eos(eos_input_rt, eos_state)

                 qm(i,j,k,QPRES,idir) = eos_state%p
                 qm(i,j,k,QREINT,idir) = qm(i,j,k,QRHO,idir)*eos_state%e
                 ! should we try to do something about Gamma_! on interface?

              end do
           end do
        end do
     end do
  end if

  ! Compute F^x at kc (k3d)
  call cmpflx([lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], &
              qm, qp, It_lo, It_hi, AMREX_SPACEDIM, 1, &
              flux1, flux1_lo, flux1_hi, &
              q_int, It_lo, It_hi, &
#ifdef RADIATION
              rflx, flux1_lo, flux1_hi, &
              lambda_int, It_lo, It_hi, &
#endif
              qaux, qa_lo, qa_hi, &
              shk, shk_lo, shk_hi, &
              1, domlo, domhi)

  call ca_store_godunov_state(lo, [hi(1)+1, hi(2), hi(3)], &
                              q_int, It_lo, It_hi, &
#ifdef RADIATION
                              lambda_int, It_lo, It_hi, &
#endif
                              q1, flux1_lo, flux1_hi)

#if AMREX_SPACEDIM >= 2
  ! Compute F^y at kc (k3d)
  call cmpflx([lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], &
              qm, qp, It_lo, It_hi, AMREX_SPACEDIM, 2, &
              flux2, flux2_lo, flux2_hi, &
              q_int, It_lo, It_hi, &  ! temporary
#ifdef RADIATION
              rfly, flux2_lo, flux2_hi, &
              lambda_int, It_lo, It_hi, &
#endif
              qaux, qa_lo, qa_hi, &
              shk, shk_lo, shk_hi, &
              2, domlo, domhi)

  call ca_store_godunov_state(lo, [hi(1), hi(2)+1, hi(3)], &
                              q_int, It_lo, It_hi, &
#ifdef RADIATION
                              lambda_int, It_lo, It_hi, &
#endif
                              q2, flux2_lo, flux2_hi)
#endif


#if AMREX_SPACEDIM == 3
  ! Compute F^z at kc (k3d)

  call cmpflx([lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1], &
              qm, qp, It_lo, It_hi, AMREX_SPACEDIM, 3, &
              flux3, flux3_lo, flux3_hi, &
              q_int, It_lo, It_hi, &
#ifdef RADIATION
              rflz, flux3_lo, flux3_hi, &
              lambda_int, It_lo, It_hi, &
#endif
              qaux, qa_lo, qa_hi, &
              shk, shk_lo, shk_hi, &
              3, domlo, domhi)

  call ca_store_godunov_state(lo, [hi(1), hi(2)+1, hi(3)], &
                              q_int, It_lo, It_hi, &
#ifdef RADIATION
                              lambda_int, It_lo, It_hi, &
#endif
                              q3, flux3_lo, flux3_hi)

#endif


  call bl_deallocate(flatn)

  call bl_deallocate(qm)
  call bl_deallocate(qp)


  call bl_deallocate(q_int)
#ifdef RADIATION
  call bl_deallocate(lambda_int)
#endif
  call bl_deallocate(shk)


  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo, hi+dg, q, q_lo, q_hi, &
       dx, div, lo, hi+dg)

  do n = 1, NVAR

     if ( n == UTEMP ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
        flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

#ifdef SHOCK_VAR
     else if ( n == USHK ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
        flux3(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif
#endif

     end if

  end do

  call apply_av(flux1_lo, flux1_hi, 1, dx, &
       div, lo, hi+dg, &
       uin, uin_lo, uin_hi, &
       flux1, flux1_lo, flux1_hi)

#if AMREX_SPACEDIM >= 2
  call apply_av(flux2_lo, flux2_hi, 2, dx, &
       div, lo, hi+dg, &
       uin, uin_lo, uin_hi, &
       flux2, flux2_lo, flux2_hi)
#endif

#if AMREX_SPACEDIM == 3
  call apply_av(flux3_lo, flux3_hi, 3, dx, &
       div, lo, hi+dg, &
       uin, uin_lo, uin_hi, &
       flux3, flux3_lo, flux3_hi)
#endif

  if (limit_fluxes_on_small_dens == 1) then
     call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
          q,q_lo,q_hi, &
          vol,vol_lo,vol_hi, &
          flux1,flux1_lo,flux1_hi, &
          area1,area1_lo,area1_hi, &
#if AMREX_SPACEDIM >= 2
          flux2,flux2_lo,flux2_hi, &
          area2,area2_lo,area2_hi, &
#endif
#if AMREX_SPACEDIM == 3
          flux3,flux3_lo,flux3_hi, &
          area3,area3_lo,area3_hi, &
#endif
          lo,hi,dt,dx)

  endif

  call normalize_species_fluxes(flux1_lo, flux1_hi, flux1, flux1_lo, flux1_hi)
#if AMREX_SPACEDIM >= 2
  call normalize_species_fluxes(flux2_lo, flux2_hi, flux2, flux2_lo, flux2_hi)
#endif

#if AMREX_SPACEDIM == 3
  call normalize_species_fluxes(flux3_lo, flux3_hi, flux3, flux3_lo, flux3_hi)
#endif

  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 1
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) ) / vol(i,j,k)

#elif AMREX_SPACEDIM == 2
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                   flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)

#else
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                   flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) + &
                   flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)
#endif

#if AMREX_SPACEDIM == 1
              if (n == UMX) then
                 update(i,j,k,UMX) = update(i,j,k,UMX) - ( q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES) ) / dx(1)
              endif
#endif

#if AMREX_SPACEDIM == 2
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

#if AMREX_SPACEDIM == 3
#ifdef HYBRID_MOMENTUM
  call add_hybrid_advection_source(lo, hi, dt, &
       update, uout_lo, uout_hi, &
       q1, flux1_lo, flux1_hi, &
       q2, flux2_lo, flux2_hi, &
       q3, flux3_lo, flux3_hi)
#endif
#endif



  ! Scale the fluxes for the form we expect later in refluxing.

  call scale_flux(flux1_lo, flux1_hi, flux1, flux1_lo, flux1_hi, area1, area1_lo, area1_hi, dt)
#if AMREX_SPACEDIM >= 2
  call scale_flux(flux2_lo, flux2_hi, flux2, flux2_lo, flux2_hi, area2, area2_lo, area2_hi, dt)
#endif
#if AMREX_SPACEDIM == 3
  call scale_flux(flux3_lo, flux3_hi, flux3, flux3_lo, flux3_hi, area3, area3_lo, area3_hi, dt)
#endif

#if AMREX_SPACEDIM == 1
  if (coord_type .eq. 0) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1
              flux1(i,j,k,UMX) = flux1(i,j,k,UMX) + dt * area1(i,j,k) * q1(i,j,k,GDPRES)
           enddo
        enddo
     enddo
  endif
#endif

#if AMREX_SPACEDIM < 3
  if (coord_type > 0) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
  end if
#endif

  call bl_deallocate(   div)

  call bl_deallocate(q1)
#if AMREX_SPACEDIM >= 2
  call bl_deallocate(q2)
#endif
#if AMREX_SPACEDIM == 3
  call bl_deallocate(q3)
#endif

#ifdef RADIATION
  call bl_deallocate(rflx)
#if AMREX_SPACEDIM >= 2
  call bl_deallocate(rfly)
#endif
#if AMREX_SPACEDIM == 3
  call bl_deallocate(rflz)
#endif
#endif

end subroutine ca_mol_single_stage


#ifndef RADIATION
module mol_module_cuda

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_construct_flux_cuda(lo, hi, domlo, domhi, dx, dt, idir, &
                                    uin, uin_lo, uin_hi, &
                                    div, div_lo, div_hi, &
                                    qaux, qa_lo, qa_hi, &
                                    shk, sk_lo, sk_hi, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    qint, qe_lo, qe_hi, &
                                    flux, f_lo, f_hi, &
                                    area, a_lo, a_hi) &
                                    bind(c,name='ca_construct_flux_cuda')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX, NQ
    use advection_util_module, only: apply_av, normalize_species_fluxes, scale_flux
    use riemann_module, only: cmpflx

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: sk_lo(3), sk_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: shk(sk_lo(1):sk_hi(1), sk_lo(2):sk_hi(2), sk_lo(3):sk_hi(3))
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,AMREX_SPACEDIM)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,AMREX_SPACEDIM)
    real(rt), intent(inout) :: qint(qe_lo(1):qe_hi(1), qe_lo(2):qe_hi(2), qe_lo(3):qe_hi(3), NQ)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    !$gpu

    call cmpflx(lo, hi, &
                qm, qp, qm_lo, qm_hi, AMREX_SPACEDIM, idir, &
                flux, f_lo, f_hi, &
                qint, qe_lo, qe_hi, &
                qaux, qa_lo, qa_hi, &
                shk, sk_lo, sk_hi, &
                idir, domlo, domhi)

    call apply_av(lo, hi, idir, dx, div, div_lo, div_hi, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
    call normalize_species_fluxes(lo, hi, flux, f_lo, f_hi)
    call scale_flux(lo, hi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)

  end subroutine ca_construct_flux_cuda

end module mol_module_cuda
#endif
