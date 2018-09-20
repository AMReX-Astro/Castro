! advection routines in support of method of lines integration

subroutine ca_fourth_single_stage(lo, hi, time, domlo, domhi, &
                                  stage_weight, &
                                  uin, uin_lo, uin_hi, &
                                  uout, uout_lo, uout_hi, &
                                  q, q_lo, q_hi, &
                                  q_bar, q_bar_lo, q_bar_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  srcU, srU_lo, srU_hi, &
                                  update, updt_lo, updt_hi, &
                                  update_flux, uf_lo, uf_hi, &
                                  dx, dt, &
                                  flx, flx_lo, flx_hi, &
#if AMREX_SPACEDIM >= 2
                                  fly, fly_lo, fly_hi, &
#endif
#if AMREX_SPACEDIM == 3
                                  flz, flz_lo, flz_hi, &
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
                                  verbose) bind(C, name="ca_fourth_single_stage")

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, NQAUX, GDPRES, &
                                 UTEMP, UEINT, USHK, GDU, GDV, GDW, UMX, &
                                 use_flattening, QPRES, NQAUX, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens, ppm_temp_fix
  use advection_util_module, only : limit_hydro_fluxes_on_small_dens, shock, &
                                    divu, normalize_species_fluxes, calc_pdivu
  use amrex_error_module
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: uflatten
  use riemann_module, only: riemann_state
  use riemann_util_module, only: compute_flux_q
  use fourth_order
  use amrex_fort_module, only : rt => amrex_real
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
  integer, intent(in) :: q_bar_lo(3), q_bar_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flx_lo(3), flx_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
#if AMREX_SPACEDIM >= 2
  integer, intent(in) :: fly_lo(3), fly_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
  integer, intent(in) :: flz_lo(3), flz_hi(3)
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
  real(rt), intent(inout) :: q_bar(q_bar_lo(1):q_bar_hi(1), q_bar_lo(2):q_bar_hi(2), q_bar_lo(3):q_bar_hi(3), NQ)
  real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
  real(rt), intent(inout) :: update_flux(uf_lo(1):uf_hi(1), uf_lo(2):uf_hi(2), uf_lo(3):uf_hi(3), NVAR)
  real(rt), intent(inout) :: flx(flx_lo(1):flx_hi(1), flx_lo(2):flx_hi(2), flx_lo(3):flx_hi(3), NVAR)
  real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
#if AMREX_SPACEDIM >= 2
  real(rt), intent(inout) :: fly(fly_lo(1):fly_hi(1), fly_lo(2):fly_hi(2), fly_lo(3):fly_hi(3), NVAR)
  real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
#endif
#if AMREX_SPACEDIM == 3
  real(rt), intent(inout) :: flz(flz_lo(1):flz_hi(1), flz_lo(2):flz_hi(2), flz_lo(3):flz_hi(3), NVAR)
  real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
#endif
#if AMREX_SPACEDIM <= 2
  real(rt), intent(inout) :: pradial(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
  real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1), dloga_lo(2):dloga_hi(2), dloga_lo(3):dloga_hi(3))
#endif
  real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt), intent(in) :: dx(3), dt, time

#ifndef RADIATION
  ! Automatic arrays for workspace
  real(rt), pointer :: flatn(:,:,:)
  real(rt), pointer :: div(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt), pointer :: qx_avg(:,:,:,:)
  real(rt), pointer :: qy_avg(:,:,:,:)
  real(rt), pointer :: qz_avg(:,:,:,:)

  real(rt), pointer :: qx_fc(:,:,:,:)
  real(rt), pointer :: qy_fc(:,:,:,:)
  real(rt), pointer :: qz_fc(:,:,:,:)

  ! Temporaries (for now)
  real(rt), pointer :: qgdnvx(:,:,:,:)
  real(rt), pointer :: qgdnvy(:,:,:,:)
  real(rt), pointer :: qgdnvz(:,:,:,:)

  real(rt), pointer :: qgdnvx_avg(:,:,:,:)
  real(rt), pointer :: qgdnvy_avg(:,:,:,:)
  real(rt), pointer :: qgdnvz_avg(:,:,:,:)

  real(rt), pointer :: flx_avg(:,:,:,:)
  real(rt), pointer :: fly_avg(:,:,:,:)
  real(rt), pointer :: flz_avg(:,:,:,:)

  real(rt), pointer :: shk(:,:,:)

  real(rt), pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
  real(rt), pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

  integer :: ngf
  integer :: It_lo(3), It_hi(3)
  integer :: st_lo(3), st_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: div1, lap
  integer :: i, j, k, n, m

  type (eos_t) :: eos_state


  ! to do 4th order for axisymmetry, we need to derive the transformations between
  ! averages and cell-centers with the correct volume terms in the integral.
#ifndef AMREX_USE_CUDA
  if (coord_type > 0) then
     call amrex_error("Error: fourth order not implemented for axisymmetric")
  endif
#endif

  ngf = 1

  It_lo = lo(:) - dg(:)
  It_hi = hi(:) + dg(:)

  shk_lo(:) = lo(:) - dg(:)
  shk_hi(:) = hi(:) + dg(:)

  call bl_allocate(   div, lo(1), hi(1)+1, lo(2), hi(2)+dg(2), lo(3), hi(3)+dg(3))

  call bl_allocate(qx_avg, q_lo, q_hi, NQ)
  call bl_allocate(qx_fc, q_lo, q_hi, NQ)
  call bl_allocate(qgdnvx, flx_lo, flx_hi, NGDNV)
  call bl_allocate(qgdnvx_avg, q_lo, q_hi, NGDNV)
  call bl_allocate(flx_avg, q_lo, q_hi, NVAR)
#if AMREX_SPACEDIM >= 2
  call bl_allocate(qy_avg, q_lo, q_hi, NQ)
  call bl_allocate(qy_fc, q_lo, q_hi, NQ)
  call bl_allocate(qgdnvy, fly_lo, fly_hi, NGDNV)
  call bl_allocate(qgdnvy_avg, q_lo, q_hi, NGDNV)
  call bl_allocate(fly_avg, q_lo, q_hi, NVAR)
#endif
#if AMREX_SPACEDIM == 3
  call bl_allocate(qz_avg, q_lo, q_hi, NQ)
  call bl_allocate(qz_fc, q_lo, q_hi, NQ)
  call bl_allocate(qgdnvz, flz_lo, flz_hi, NGDNV)
  call bl_allocate(qgdnvz_avg, q_lo, q_hi, NGDNV)
  call bl_allocate(flz_avg, q_lo, q_hi, NVAR)
#endif

  call bl_allocate(qxm, q_lo, q_hi, NQ)
  call bl_allocate(qxp, q_lo, q_hi, NQ)

#if AMREX_SPACEDIM >= 2
  call bl_allocate(qym, q_lo, q_hi, NQ)
  call bl_allocate(qyp, q_lo, q_hi, NQ)
#endif
#if AMREX_SPACEDIM == 3
  call bl_allocate(qzm, q_lo, q_hi, NQ)
  call bl_allocate(qzp, q_lo, q_hi, NQ)
#endif

  call bl_allocate(shk, shk_lo, shk_hi)

#ifdef SHOCK_VAR
  uout(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), USHK) = ZERO

  call shock(q_bar, q_bar_lo, q_bar_hi, shk, shk_lo, shk_hi, lo, hi, dx)

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
     call shock(q_bar, q_bar_lo, q_bar_hi, shk, shk_lo, shk_hi, lo, hi, dx)
  else
     shk(:,:,:) = ZERO
  endif
#endif

  ! Compute flattening coefficient for slope calculations -- we do
  ! this with q_bar, since we need all of the ghost cells
  call bl_allocate(flatn, q_bar_lo, q_bar_hi)

  if (use_flattening == 1) then
     call uflatten(lo - ngf*dg, hi + ngf*dg, &
                   q_bar, flatn, q_bar_lo, q_bar_hi, QPRES)
  else
     flatn = ONE
  endif

  ! in contrast to the other solvers, we do not use 2-d slabs for 3-d,
  ! but we consider the full 3-d box at once.


  ! do the reconstruction here -- get the interface states

  do n = 1, NQ

     ! x-interfaces
     call states(1, &
                 q, q_lo, q_hi, NQ, n, &
                 flatn, q_bar_lo, q_bar_hi, &
                 qxm, qxp, q_lo, q_hi, &
                 lo, hi)

#if AMREX_SPACEDIM >= 2
     ! y-interfaces
     call states(2, &
                 q, q_lo, q_hi, NQ, n, &
                 flatn, q_bar_lo, q_bar_hi, &
                 qym, qyp, q_lo, q_hi, &
                 lo, hi)
#endif

#if AMREX_SPACEDIM == 3
     ! z-interfaces
     call states(3, &
                 q, q_lo, q_hi, NQ, n, &
                 flatn, q_bar_lo, q_bar_hi, &
                 qzm, qzp, q_lo, q_hi, &
                 lo, hi)
#endif

  enddo

  ! this is where we would implement ppm_temp_fix


  ! solve the Riemann problems -- we just require the interface state
  ! at this point

  ! note that the Riemann solver is written to work in slabs, so we
  ! need to pass the k index for both the state and flux separately.

  ! TODO: we should explicitly compute a gamma with this state, since
  ! we cannot get away with the first-order construction that we pull
  ! from qaux in the Riemann solver

  do k = lo(3)-dg(3), hi(3)+dg(3)

     call riemann_state(qxm, qxp, q_lo, q_hi, &
                        qx_avg, q_lo, q_hi, &
                        qaux, qa_lo, qa_hi, &
                        1, [lo(1), lo(2)-dg(2), k], [hi(1)+1, hi(2)+dg(2), k], domlo, domhi)

     call compute_flux_q(1, qx_avg, q_lo, q_hi, &
                         flx_avg, q_lo, q_hi, &
                         qgdnvx_avg, q_lo, q_hi, &
                         [lo(1), lo(2)-dg(2), k], [hi(1)+1, hi(2)+dg(2), k])
  enddo

#if AMREX_SPACEDIM >= 2
  do k = lo(3)-dg(3), hi(3)+dg(3)

     call riemann_state(qym, qyp, q_lo, q_hi, &
                        qy_avg, q_lo, q_hi, &
                        qaux, qa_lo, qa_hi, &
                        2, [lo(1)-1, lo(2), k], [hi(1)+1, hi(2)+1, k], domlo, domhi)

     call compute_flux_q(2, qy_avg, q_lo, q_hi, &
                         fly_avg, q_lo, q_hi, &
                         qgdnvy_avg, q_lo, q_hi, &
                         [lo(1)-1, lo(2), k], [hi(1)+1, hi(2)+1, k])
  enddo
#endif

#if AMREX_SPACEDIM == 3
  do k = lo(3), hi(3)+dg(3)

     call riemann_state(qzm, qzp, q_lo, q_hi, &
                        qz_avg, q_lo, q_hi, &
                        qaux, qa_lo, qa_hi, &
                        3, [lo(1)-1, lo(2)-1, k], [hi(1)+1, hi(2)+1, k], domlo, domhi)

     call compute_flux_q(3, qz_avg, q_lo, q_hi, &
                         flz_avg, q_lo, q_hi, &
                         qgdnvz_avg, q_lo, q_hi, &
                         [lo(1)-1, lo(2)-1, k], [hi(1)+1, hi(2)+1, k])
  enddo
#endif


  call bl_deallocate(flatn)

  call bl_deallocate(qxm)
  call bl_deallocate(qxp)

#if AMREX_SPACEDIM >= 2
  call bl_deallocate(qym)
  call bl_deallocate(qyp)
#endif

#if AMREX_SPACEDIM == 3
  call bl_deallocate(qzm)
  call bl_deallocate(qzp)
#endif

  call bl_deallocate(shk)

  ! we now have the face-average interface states and fluxes evaluated with these
  ! for 1-d, we are done


  ! construct the face-center interface states

#if AMREX_SPACEDIM >= 2
  ! x-interfaces
  do n = 1, NQ
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1

              ! note: need to consider axisymmetry in the future
              lap = qx_avg(i,j+1,k,n) - TWO*qx_avg(i,j,k,n) + qx_avg(i,j-1,k,n)
#if AMREX_SPACEDIM == 3
              lap = lap + qx_avg(i,j,k+1,n) - TWO*qx_avg(i,j,k,n) + qx_avg(i,j,k-1,n)
#endif
              qx_fc(i,j,k,n) = qx_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
           enddo
        enddo
     enddo
  enddo

  ! y-interfaces
  do n = 1, NQ
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)

              ! note: need to consider axisymmetry in the future
              lap = qy_avg(i+1,j,k,n) - TWO*qy_avg(i,j,k,n) + qy_avg(i-1,j,k,n)
#if AMREX_SPACEDIM == 3
              lap = lap + qy_avg(i,j,k+1,n) - TWO*qy_avg(i,j,k,n) + qy_avg(i,j,k-1,n)
#endif
              qy_fc(i,j,k,n) = qy_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
           enddo
        enddo
     enddo
  enddo

#if AMREX_SPACEDIM == 3
  ! z-interfaces
  do n = 1, NQ
     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              ! note: need to consider axisymmetry in the future
              lap = qz_avg(i+1,j,k,n) - TWO*qz_avg(i,j,k,n) + qz_avg(i-1,j,k,n)
              lap = lap + qz_avg(i,j+1,k,n) - TWO*qz_avg(i,j,k,n) + qz_avg(i,j-1,k,n)

              qz_fc(i,j,k,n) = qz_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
           enddo
        enddo
     enddo
  enddo

#endif


  ! compute face-centered fluxes
  ! these will be stored in flx, fly, flz
  call compute_flux_q(1, qx_fc, q_lo, q_hi, &
                      flx, flx_lo, flx_hi, &
                      qgdnvx, flx_lo, flx_hi, &
                      [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)])

#if AMREX_SPACEDIM >= 2
  call compute_flux_q(2, qy_fc, q_lo, q_hi, &
                      fly, fly_lo, fly_hi, &
                      qgdnvy, fly_lo, fly_hi, &
                      [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)])
#endif

#if AMREX_SPACEDIM == 3
  call compute_flux_q(3, qz_fc, q_lo, q_hi, &
                      flz, flz_lo, flz_hi, &
                      qgdnvz, flz_lo, flz_hi, &
                      [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1])
#endif


  ! compute the final fluxes include the transverse correction
  ! x-interfaces
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1

              lap = flx_avg(i,j+1,k,n) - TWO*flx_avg(i,j,k,n) + flx_avg(i,j-1,k,n)
#if AMREX_SPACEDIM == 3
              lap = lap + flx_avg(i,j,k+1,n) - TWO*flx_avg(i,j,k,n) + flx_avg(i,j,k-1,n)
#endif
              flx(i,j,k,n) = flx(i,j,k,n) + 1.0_rt/24.0_rt * lap

           enddo
        enddo
     enddo
  enddo

  do n = 1, NGDNV
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1

              lap = qgdnvx_avg(i,j+1,k,n) - TWO*qgdnvx_avg(i,j,k,n) + qgdnvx_avg(i,j-1,k,n)
#if AMREX_SPACEDIM == 3
              lap = lap + qgdnvx_avg(i,j,k+1,n) - TWO*qgdnvx_avg(i,j,k,n) + qgdnvx_avg(i,j,k-1,n)
#endif
              qgdnvx(i,j,k,n) = qgdnvx(i,j,k,n) + 1.0_rt/24.0_rt * lap

           enddo
        enddo
     enddo
  enddo


  ! y-interfaces
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)

              lap = fly_avg(i+1,j,k,n) - TWO*fly_avg(i,j,k,n) + fly_avg(i-1,j,k,n)
#if AMREX_SPACEDIM == 3
              lap = lap + fly_avg(i,j,k+1,n) - TWO*fly_avg(i,j,k,n) + fly_avg(i,j,k-1,n)
#endif
              fly(i,j,k,n) = fly(i,j,k,n) + 1.0_rt/24.0_rt * lap
           enddo
        enddo
     enddo
  enddo

  do n = 1, NGDNV
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)

              lap = qgdnvy_avg(i+1,j,k,n) - TWO*qgdnvy_avg(i,j,k,n) + qgdnvy_avg(i-1,j,k,n)
#if AMREX_SPACEDIM == 3
              lap = lap + qgdnvy_avg(i,j,k+1,n) - TWO*qgdnvy_avg(i,j,k,n) + qgdnvy_avg(i,j,k-1,n)
#endif
              qgdnvy(i,j,k,n) = qgdnvy(i,j,k,n) + 1.0_rt/24.0_rt * lap
           enddo
        enddo
     enddo
  enddo

#if AMREX_SPACEDIM == 3
  ! z-interfaces
  do n = 1, NVAR
     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              lap = flz_avg(i+1,j,k,n) - TWO*flz_avg(i,j,k,n) + flz_avg(i-1,j,k,n)
              lap = lap + flz_avg(i,j+1,k,n) - TWO*flz_avg(i,j,k,n) + flz_avg(i,j-1,k,n)
              flz(i,j,k,n) = flz(i,j,k,n) + 1.0_rt/24.0_rt * lap

           enddo
        enddo
     enddo
  enddo

  do n = 1, NGDNV
     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              lap = qgdnvz_avg(i+1,j,k,n) - TWO*qgdnvz_avg(i,j,k,n) + qgdnvz_avg(i-1,j,k,n)
              lap = lap + qgdnvz_avg(i,j+1,k,n) - TWO*qgdnvz_avg(i,j,k,n) + qgdnvz_avg(i,j-1,k,n)
              qgdnvz(i,j,k,n) = qgdnvz(i,j,k,n) + 1.0_rt/24.0_rt * lap

           enddo
        enddo
     enddo
  enddo
#endif

#else
  ! for 1-d, we just copy flx_avg -> flx, since there is no face averaging
  flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = flx_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
#endif

  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo, hi, q, q_lo, q_hi, &
            dx, div, lo, hi+dg)

  do n = 1, NVAR

     if ( n == UTEMP ) then
        flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
        fly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
        flz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif

#ifdef SHOCK_VAR
     else if ( n == USHK ) then
        flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
#if AMREX_SPACEDIM >= 2
        fly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO
#endif
#if AMREX_SPACEDIM == 3
        flz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,n) = ZERO
#endif
#endif

     else
        ! do the artificial viscosity
        continue
#ifdef THIS_IS_NOT_FOURTH_ORDER_ACCURATE
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)+1

                 div1 = FOURTH*(div(i,j,k) + div(i,j+dg(2),k) + &
                                div(i,j,k+dg(3)) + div(i,j+dg(2),k+dg(3)))
                 div1 = difmag*min(ZERO, div1)

                 flx(i,j,k,n) = flx(i,j,k,n) + &
                      dx(1) * div1 * (uin(i,j,k,n) - uin(i-1,j,k,n))
              enddo
           enddo
        enddo
#if AMREX_SPACEDIM >= 2
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)+1
              do i = lo(1), hi(1)
                 div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                div(i,j,k+dg(3)) + div(i+1,j,k+dg(3)))
                 div1 = difmag*min(ZERO, div1)

                 fly(i,j,k,n) = fly(i,j,k,n) + &
                      dx(2) * div1 * (uin(i,j,k,n) - uin(i,j-1,k,n))
              enddo
           enddo
        enddo
#endif
#if AMREX_SPACEDIM == 3
        do k = lo(3), hi(3)+1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                div(i,j+1,k) + div(i+1,j+1,k))
                 div1 = difmag*min(ZERO, div1)

                 flz(i,j,k,n) = flz(i,j,k,n) + &
                      dx(3) * div1 * (uin(i,j,k,n) - uin(i,j,k-1,n))
              enddo
           enddo
        enddo
#endif
#endif  
     endif

  enddo

  call normalize_species_fluxes(flx_lo, flx_hi, flx, flx_lo, flx_hi)
#if AMREX_SPACEDIM >= 2
  call normalize_species_fluxes(fly_lo, fly_hi, fly, fly_lo, fly_hi)
#endif
#if AMREX_SPACEDIM == 3
  call normalize_species_fluxes(flz_lo, flz_hi, flz, flz_lo, flz_hi)
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
                   (flx(i,j,k,n) * area1(i,j,k) - flx(i+1,j,k,n) * area1(i+1,j,k) ) / vol(i,j,k)

#elif AMREX_SPACEDIM == 2
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flx(i,j,k,n) * area1(i,j,k) - flx(i+1,j,k,n) * area1(i+1,j,k) + &
                    fly(i,j,k,n) * area2(i,j,k) - fly(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)

#else
              update(i,j,k,n) = update(i,j,k,n) + &
                   (flx(i,j,k,n) * area1(i,j,k) - flx(i+1,j,k,n) * area1(i+1,j,k) + &
                    fly(i,j,k,n) * area2(i,j,k) - fly(i,j+1,k,n) * area2(i,j+1,k) + &
                    flz(i,j,k,n) * area3(i,j,k) - flz(i,j,k+1,n) * area3(i,j,k+1) ) / vol(i,j,k)
#endif

#if AMREX_SPACEDIM == 1
              if (n == UMX) then
                 update(i,j,k,UMX) = update(i,j,k,UMX) - ( qgdnvx(i+1,j,k,GDPRES) - qgdnvx(i,j,k,GDPRES) ) / dx(1)
              endif
#endif

#if AMREX_SPACEDIM == 2
              if (n == UMX) then
                 ! add the pressure source term for axisymmetry
                 if (coord_type > 0) then
                    update(i,j,k,n) = update(i,j,k,n) - (qgdnvx(i+1,j,k,GDPRES) - qgdnvx(i,j,k,GDPRES))/ dx(1)
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
                                   qgdnvx, flx_lo, flx_hi, &
                                   qgdnvy, fly_lo, fly_hi, &
                                   qgdnvz, flz_lo, flz_hi)
#endif
#endif



  ! Scale the fluxes for the form we expect later in refluxing.

  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) + 1
              flx(i,j,k,n) = dt * flx(i,j,k,n) * area1(i,j,k)

#if AMREX_SPACEDIM == 1
              if (coord_type .eq. 0 .and. n == UMX) then
                 flx(i,j,k,n) = flx(i,j,k,n) + dt * area1(i,j,k) * qgdnvx(i,j,k,GDPRES)
              endif
#endif

           enddo
        enddo
     enddo
  enddo

#if AMREX_SPACEDIM >= 2
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2) + 1
           do i = lo(1), hi(1)
              fly(i,j,k,n) = dt * fly(i,j,k,n) * area2(i,j,k)
           enddo
        enddo
     enddo
  enddo
#endif

#if AMREX_SPACEDIM == 3
  do n = 1, NVAR
     do k = lo(3), hi(3) + 1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              flz(i,j,k,n) = dt * flz(i,j,k,n) * area3(i,j,k)
           enddo
        enddo
     enddo
  enddo
#endif

#if AMREX_SPACEDIM < 3
  if (coord_type > 0) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = qgdnvx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
  end if
#endif

  call bl_deallocate(   div)

  call bl_deallocate(qx_avg)
  call bl_deallocate(qx_fc)
  call bl_deallocate(qgdnvx)
  call bl_deallocate(qgdnvx_avg)
  call bl_deallocate(flx_avg)
#if AMREX_SPACEDIM >= 2
  call bl_deallocate(qy_avg)
  call bl_deallocate(qy_fc)
  call bl_deallocate(qgdnvy)
  call bl_deallocate(qgdnvy_avg)
  call bl_deallocate(fly_avg)
#endif
#if AMREX_SPACEDIM == 3
  call bl_deallocate(qz_avg)
  call bl_deallocate(qz_fc)
  call bl_deallocate(qgdnvz)
  call bl_deallocate(qgdnvz_avg)
  call bl_deallocate(flz_avg)
#endif
#else
#ifndef AMREX_USE_CUDA
  ! RADIATION check
  call amrex_error("ERROR: ca_fourth_single_stage does not support radiation")
#endif
#endif
end subroutine ca_fourth_single_stage
