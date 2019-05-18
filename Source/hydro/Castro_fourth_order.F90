! advection routines in support of method of lines integration
!
subroutine ca_fourth_single_stage(lo, hi, time, domlo, domhi, &
                                  uin, uin_lo, uin_hi, &
                                  uout, uout_lo, uout_hi, &
                                  q_core, qc_lo, qc_hi, &
                                  q_pass, qp_lo, qp_hi, &
                                  q_core_bar, qc_bar_lo, qc_bar_hi, &
                                  q_pass_bar, qp_bar_lo, qp_bar_hi, &
                                  qaux, qa_lo, qa_hi, &
                                  qaux_bar, qa_bar_lo, qa_bar_hi, &
                                  srcU, srU_lo, srU_hi, &
                                  update, updt_lo, updt_hi, &
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
  use meth_params_module, only : NQC, NQP, NVAR, NGDNV, NQAUX, GDPRES, &
                                 UTEMP, UEINT, USHK, GDU, GDV, GDW, UMX, &
                                 use_flattening, QPRES, NQAUX, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, QGAME, QGC, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens, ppm_temp_fix, do_hydro
  use advection_util_module, only : limit_hydro_fluxes_on_small_dens, ca_shock, &
                                    normalize_species_fluxes, avisc

  use amrex_error_module
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: ca_uflatten
  use riemann_module, only: riemann_state
  use riemann_util_module, only: compute_flux_q
  use fourth_order
  use amrex_fort_module, only : rt => amrex_real
#ifdef HYBRID_MOMENTUM
  use hybrid_advection_module, only : add_hybrid_advection_source
  use riemann_util_module, only : ca_store_godunov_state
#endif
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec, naux
  use prob_params_module, only : dg, coord_type

  implicit none

  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) ::  domlo(3), domhi(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: qc_lo(3), qc_hi(3)
  integer, intent(in) :: qp_lo(3), qp_hi(3)
  integer, intent(in) :: qc_bar_lo(3), qc_bar_hi(3)
  integer, intent(in) :: qp_bar_lo(3), qp_bar_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: qa_bar_lo(3), qa_bar_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
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
  real(rt), intent(inout) :: q_core(qc_lo(1):qc_hi(1), qc_lo(2):qc_hi(2), qc_lo(3):qc_hi(3), NQC)
  real(rt), intent(inout) :: q_pass(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQP)
  real(rt), intent(inout) :: q_core_bar(qc_bar_lo(1):qc_bar_hi(1), qc_bar_lo(2):qc_bar_hi(2), qc_bar_lo(3):qc_bar_hi(3), NQC)
  real(rt), intent(inout) :: q_pass_bar(qp_bar_lo(1):qp_bar_hi(1), qp_bar_lo(2):qp_bar_hi(2), qp_bar_lo(3):qp_bar_hi(3), NQP)
  real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt), intent(inout) :: qaux_bar(qa_bar_lo(1):qa_bar_hi(1), qa_bar_lo(2):qa_bar_hi(2), qa_bar_lo(3):qa_bar_hi(3), NQAUX)
  real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
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
  real(rt), pointer :: avisx(:,:,:), avisy(:,:,:), avisz(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt), pointer :: qx_core_avg(:,:,:,:)
  real(rt), pointer :: qy_core_avg(:,:,:,:)
  real(rt), pointer :: qz_core_avg(:,:,:,:)

  real(rt), pointer :: qx_pass_avg(:,:,:,:)
  real(rt), pointer :: qy_pass_avg(:,:,:,:)
  real(rt), pointer :: qz_pass_avg(:,:,:,:)

  real(rt), pointer :: qx_core(:,:,:,:)
  real(rt), pointer :: qy_core(:,:,:,:)
  real(rt), pointer :: qz_core(:,:,:,:)

  real(rt), pointer :: qx_pass(:,:,:,:)
  real(rt), pointer :: qy_pass(:,:,:,:)
  real(rt), pointer :: qz_pass(:,:,:,:)

#ifdef HYBRID_MOMENTUM
  real(rt), pointer :: qgdnvx(:,:,:,:)
  real(rt), pointer :: qgdnvy(:,:,:,:)
  real(rt), pointer :: qgdnvz(:,:,:,:)
#endif

  ! Temporaries (for now)
  real(rt), pointer :: flx_avg(:,:,:,:)
  real(rt), pointer :: fly_avg(:,:,:,:)
  real(rt), pointer :: flz_avg(:,:,:,:)

  real(rt), pointer :: shk(:,:,:)

  real(rt), pointer :: qxm_core(:,:,:,:), qym_core(:,:,:,:), qzm_core(:,:,:,:)
  real(rt), pointer :: qxp_core(:,:,:,:), qyp_core(:,:,:,:), qzp_core(:,:,:,:)

  real(rt), pointer :: qxm_pass(:,:,:,:), qym_pass(:,:,:,:), qzm_pass(:,:,:,:)
  real(rt), pointer :: qxp_pass(:,:,:,:), qyp_pass(:,:,:,:), qzp_pass(:,:,:,:)

  integer :: ngf
  integer :: It_lo(3), It_hi(3)
  integer :: st_lo(3), st_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: lap
  integer :: i, j, k, n, m

  type (eos_t) :: eos_state

  ! artifical viscosity strength
  real(rt), parameter :: alpha = 0.3_rt
  real(rt) :: avisc_coeff

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

  call bl_allocate(avisx, lo, hi+dg)
#if BL_SPACEDIM >= 2
  call bl_allocate(avisy, lo, hi+dg)
#endif
#if BL_SPACEDIM == 3
  call bl_allocate(avisz, lo, hi+dg)
#endif

  call bl_allocate(qx_core_avg, qc_lo, qc_hi, NQC)
  call bl_allocate(qx_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qx_pass_avg, qc_lo, qc_hi, NQP)
  call bl_allocate(qx_pass, qc_lo, qc_hi, NQP)

  call bl_allocate(flx_avg, qc_lo, qc_hi, NVAR)

#if AMREX_SPACEDIM >= 2
  call bl_allocate(qy_core_avg, qc_lo, qc_hi, NQC)
  call bl_allocate(qy_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qy_pass_avg, qc_lo, qc_hi, NQP)
  call bl_allocate(qy_pass, qc_lo, qc_hi, NQP)

  call bl_allocate(fly_avg, qc_lo, qc_hi, NVAR)
#endif

#if AMREX_SPACEDIM == 3
  call bl_allocate(qz_core_avg, qc_lo, qc_hi, NQC)
  call bl_allocate(qz_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qz_pass_avg, qc_lo, qc_hi, NQP)
  call bl_allocate(qz_pass, qc_lo, qc_hi, NQP)

  call bl_allocate(flz_avg, qc_lo, qc_hi, NVAR)
#endif

  call bl_allocate(qxm_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qxp_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qxm_pass, qc_lo, qc_hi, NQP)
  call bl_allocate(qxp_pass, qc_lo, qc_hi, NQP)

#if AMREX_SPACEDIM >= 2
  call bl_allocate(qym_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qyp_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qym_pass, qc_lo, qc_hi, NQP)
  call bl_allocate(qyp_pass, qc_lo, qc_hi, NQP)
#endif
#if AMREX_SPACEDIM == 3
  call bl_allocate(qzm_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qzp_core, qc_lo, qc_hi, NQC)
  call bl_allocate(qzm_pass, qc_lo, qc_hi, NQP)
  call bl_allocate(qzp_pass, qc_lo, qc_hi, NQP)
#endif

  call bl_allocate(shk, shk_lo, shk_hi)

  if (do_hydro == 1) then

#ifdef SHOCK_VAR
     uout(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), USHK) = ZERO

     call ca_shock(lo-dg, hi+dg, &
                   q_core_bar, qc_bar_lo, qc_bar_hi, &
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
                      q_core_bar, qc_bar_lo, qc_bar_hi, &
                      shk, shk_lo, shk_hi, &
                      dx)
     else
        shk(:,:,:) = ZERO
     endif
#endif

     ! Compute flattening coefficient for slope calculations -- we do
     ! this with q_bar, since we need all of the ghost cells
     call bl_allocate(flatn, qc_bar_lo, qc_bar_hi)

     if (use_flattening == 1) then
        call ca_uflatten(lo - ngf*dg, hi + ngf*dg, &
                         q_core_bar, qc_bar_lo, qc_bar_hi, &
                         q_core_bar, qc_bar_lo, qc_bar_hi, NQC, QPRES, &
                         flatn, qc_bar_lo, qc_bar_hi)
     else
        flatn = ONE
     endif

     ! do the reconstruction here -- get the interface states

     do n = 1, NQC

        ! x-interfaces
        call states(1, &
                    q_core, qc_lo, qc_hi, NQC, n, &
                    flatn, qc_bar_lo, qc_bar_hi, &
                    qxm_core, qxp_core, qc_lo, qc_hi, &
                    lo, hi)

#if AMREX_SPACEDIM >= 2
        ! y-interfaces
        call states(2, &
                    q_core, qc_lo, qc_hi, NQC, n, &
                    flatn, qc_bar_lo, qc_bar_hi, &
                    qym_core, qyp_core, qc_lo, qc_hi, &
                    lo, hi)
#endif

#if AMREX_SPACEDIM == 3
        ! z-interfaces
        call states(3, &
                    q_core, qc_lo, qc_hi, NQC, n, &
                    flatn, qc_bar_lo, qc_bar_hi, &
                    qzm_core, qzp_core, qc_lo, qc_hi, &
                    lo, hi)
#endif

     enddo

     do n = 1, NQP

        ! x-interfaces
        call states(1, &
                    q_pass, qp_lo, qp_hi, NQP, n, &
                    flatn, qc_bar_lo, qc_bar_hi, &
                    qxm_pass, qxp_pass, qc_lo, qc_hi, &
                    lo, hi)

#if AMREX_SPACEDIM >= 2
        ! y-interfaces
        call states(2, &
                    q_pass, qp_lo, qp_hi, NQP, n, &
                    flatn, qc_bar_lo, qc_bar_hi, &
                    qym_pass, qyp_pass, qc_lo, qc_hi, &
                    lo, hi)
#endif

#if AMREX_SPACEDIM == 3
        ! z-interfaces
        call states(3, &
                    q_pass, qp_lo, qp_hi, NQP, n, &
                    flatn, qc_bar_lo, qc_bar_hi, &
                    qzm_pass, qzp_pass, qc_lo, qc_hi, &
                    lo, hi)
#endif

     enddo

     ! this is where we would implement ppm_temp_fix


     ! solve the Riemann problems -- we just require the interface state
     ! at this point

     call riemann_state(qxm_core, qc_lo, qc_hi, &
                        qxp_core, qc_lo, qc_hi, 1, 1, &
                        qxm_pass, qc_lo, qc_hi, &
                        qxp_pass, qc_lo, qc_hi, &
                        qx_core_avg, qc_lo, qc_hi, &
                        qx_pass_avg, qc_lo, qc_hi, &
                        qaux, qa_lo, qa_hi, &
                        1, &
                        [lo(1), lo(2)-dg(2), lo(3)-dg(3)], &
                        [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                        domlo, domhi)

     call compute_flux_q([lo(1), lo(2)-dg(2), lo(3)-dg(3)], &
                         [hi(1)+1, hi(2)+dg(2), hi(3)+dg(3)], &
                         qx_core_avg, qc_lo, qc_hi, &
                         qx_pass_avg, qc_lo, qc_hi, &
                         flx_avg, qc_lo, qc_hi, &
                         1)


#if AMREX_SPACEDIM >= 2
     call riemann_state(qym_core, qc_lo, qc_hi, &
                        qyp_core, qc_lo, qc_hi, 1, 1, &
                        qym_pass, qc_lo, qc_hi, &
                        qyp_pass, qc_lo, qc_hi, &
                        qy_core_avg, qc_lo, qc_hi, &
                        qy_pass_avg, qc_lo, qc_hi, &
                        qaux, qa_lo, qa_hi, &
                        2, &
                        [lo(1)-1, lo(2), lo(3)-dg(3)], &
                        [hi(1)+1, hi(2)+1, hi(3)+dg(3)], &
                        domlo, domhi)

     call compute_flux_q([lo(1)-1, lo(2), lo(3)-dg(3)], &
                         [hi(1)+1, hi(2)+1, hi(3)+dg(3)], &
                         qy_core_avg, qc_lo, qc_hi, &
                         qy_pass_avg, qc_lo, qc_hi, &
                         fly_avg, qc_lo, qc_hi, &
                         2)
#endif

#if AMREX_SPACEDIM == 3
     call riemann_state(qzm_core, qc_lo, qc_hi, &
                        qzp_core, qc_lo, qc_hi, 1, 1, &
                        qzm_pass, qc_lo, qc_hi, &
                        qzp_pass, qc_lo, qc_hi, &
                        qz_core_avg, qc_lo, qc_hi, &
                        qz_pass_avg, qc_lo, qc_hi, &
                        qaux, qa_lo, qa_hi, &
                        3, &
                        [lo(1)-1, lo(2)-1, lo(3)], &
                        [hi(1)+1, hi(2)+1, hi(3)+1], &
                        domlo, domhi)

     call compute_flux_q([lo(1)-1, lo(2)-1, lo(3)], &
                         [hi(1)+1, hi(2)+1, hi(3)+1], &
                         qz_core_avg, qc_lo, qc_hi, &
                         qz_pass_avg, qc_lo, qc_hi, &
                         flz_avg, qc_lo, qc_hi, &
                         3)
#endif


     call bl_deallocate(flatn)

     call bl_deallocate(qxm_core)
     call bl_deallocate(qxp_core)
     call bl_deallocate(qxm_pass)
     call bl_deallocate(qxp_pass)

#if AMREX_SPACEDIM >= 2
     call bl_deallocate(qym_core)
     call bl_deallocate(qyp_core)
     call bl_deallocate(qym_pass)
     call bl_deallocate(qyp_pass)
#endif

#if AMREX_SPACEDIM == 3
     call bl_deallocate(qzm_core)
     call bl_deallocate(qzp_core)
     call bl_deallocate(qzm_pass)
     call bl_deallocate(qzp_pass)
#endif

     call bl_deallocate(shk)

     ! we now have the face-average interface states and fluxes evaluated with these
     ! for 1-d, we are done


     ! construct the face-center interface states

#if AMREX_SPACEDIM >= 2

     do n = 1, NQC
        if (n == QGAME .or. n == QGC .or. n == QTEMP) cycle

        ! x-interfaces
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)+1

                 ! note: need to consider axisymmetry in the future
                 lap = qx_core_avg(i,j+1,k,n) - TWO*qx_core_avg(i,j,k,n) + qx_core_avg(i,j-1,k,n)
#if AMREX_SPACEDIM == 3
                 lap = lap + qx_core_avg(i,j,k+1,n) - TWO*qx_core_avg(i,j,k,n) + qx_core_avg(i,j,k-1,n)
#endif
                 qx_core(i,j,k,n) = qx_core_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
              end do
           end do
        end do

        ! y-interfaces
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)+1
              do i = lo(1), hi(1)

                 ! note: need to consider axisymmetry in the future
                 lap = qy_core_avg(i+1,j,k,n) - TWO*qy_core_avg(i,j,k,n) + qy_core_avg(i-1,j,k,n)
#if AMREX_SPACEDIM == 3
                 lap = lap + qy_core_avg(i,j,k+1,n) - TWO*qy_core_avg(i,j,k,n) + qy_core_avg(i,j,k-1,n)
#endif
                 qy_core(i,j,k,n) = qy_core_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
              end do
           end do
        end do

#if AMREX_SPACEDIM == 3
        ! z-interfaces
        do k = lo(3), hi(3)+1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 ! note: need to consider axisymmetry in the future
                 lap = qz_core_avg(i+1,j,k,n) - TWO*qz_core_avg(i,j,k,n) + qz_core_avg(i-1,j,k,n)
                 lap = lap + qz_core_avg(i,j+1,k,n) - TWO*qz_core_avg(i,j,k,n) + qz_core_avg(i,j-1,k,n)

                 qz_core(i,j,k,n) = qz_core_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
              end do
           end do
        end do
#endif
     end do


     do n = 1, NQP

        ! x-interfaces
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)+1

                 ! note: need to consider axisymmetry in the future
                 lap = qx_pass_avg(i,j+1,k,n) - TWO*qx_pass_avg(i,j,k,n) + qx_pass_avg(i,j-1,k,n)
#if AMREX_SPACEDIM == 3
                 lap = lap + qx_pass_avg(i,j,k+1,n) - TWO*qx_pass_avg(i,j,k,n) + qx_pass_avg(i,j,k-1,n)
#endif
                 qx_pass(i,j,k,n) = qx_pass_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
              end do
           end do
        end do

        ! y-interfaces
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)+1
              do i = lo(1), hi(1)

                 ! note: need to consider axisymmetry in the future
                 lap = qy_pass_avg(i+1,j,k,n) - TWO*qy_pass_avg(i,j,k,n) + qy_pass_avg(i-1,j,k,n)
#if AMREX_SPACEDIM == 3
                 lap = lap + qy_pass_avg(i,j,k+1,n) - TWO*qy_pass_avg(i,j,k,n) + qy_pass_avg(i,j,k-1,n)
#endif
                 qy_pass(i,j,k,n) = qy_pass_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
              end do
           end do
        end do

#if AMREX_SPACEDIM == 3
        ! z-interfaces
        do k = lo(3), hi(3)+1
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)

                 ! note: need to consider axisymmetry in the future
                 lap = qz_pass_avg(i+1,j,k,n) - TWO*qz_pass_avg(i,j,k,n) + qz_pass_avg(i-1,j,k,n)
                 lap = lap + qz_pass_avg(i,j+1,k,n) - TWO*qz_pass_avg(i,j,k,n) + qz_pass_avg(i,j-1,k,n)

                 qz_pass(i,j,k,n) = qz_pass_avg(i,j,k,n) - 1.0_rt/24.0_rt * lap
              end do
           end do
        end do
#endif
     end do

     ! compute face-centered fluxes
     ! these will be stored in flx, fly, flz
     call compute_flux_q([lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], &
                         qx_core, qc_lo, qc_hi, &
                         qx_pass, qc_lo, qc_hi, &
                         flx, flx_lo, flx_hi, &
                         1)

     call compute_flux_q([lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], &
                         qy_core, qc_lo, qc_hi, &
                         qy_pass, qc_lo, qc_hi, &
                         fly, fly_lo, fly_hi, &
                         2)

#if AMREX_SPACEDIM == 3
     call compute_flux_q([lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1], &
                         qz_core, qc_lo, qc_hi, &
                         qz_pass, qc_lo, qc_hi, &
                         flz, flz_lo, flz_hi, &
                         3)
#endif

     call bl_deallocate(qx_core)
     call bl_deallocate(qx_pass)
     call bl_deallocate(qy_core)
     call bl_deallocate(qy_pass)
#if AMREX_SPACEDIM == 3
     call bl_deallocate(qz_core)
     call bl_deallocate(qz_pass)
#endif

     ! compute the final fluxes (as an average over the interface), this
     ! requires a transverse correction.  Note, we don't need to do anything
     ! to get the average of the Godunov states over the interface--this is
     ! essentially what qx_avg already is

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
#endif

#else
     ! for 1-d, we just copy flx_avg -> flx, since there is no face averaging
     flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = flx_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
#endif


     ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
     call avisc(lo, hi, &
                q_core_bar, qc_bar_lo, qc_bar_hi, &
                qaux_bar, qa_bar_lo, qa_bar_hi, &
                dx, avisx, lo, hi+dg, 1)

#if BL_SPACEDIM >= 2
     call avisc(lo, hi, &
                q_core_bar, qc_bar_lo, qc_bar_hi, &
                qaux_bar, qa_bar_lo, qa_bar_hi, &
                dx, avisy, lo, hi+dg, 2)
#endif

#if BL_SPACEDIM == 3
     call avisc(lo, hi, &
                q_core_bar, qc_bar_lo, qc_bar_hi, &
                qaux_bar, qa_bar_lo, qa_bar_hi, &
                dx, avisz, lo, hi+dg, 3)
#endif

     ! avisc_coefficient is the coefficent we use.  The McCorquodale &
     ! Colella paper suggest alpha = 0.3, but our other hydro solvers use
     ! a coefficient on the divergence that defaults to 0.1, so we
     ! normalize to that value, to allow for adjustments
     avisc_coeff = alpha * (difmag / 0.1_rt)

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

           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)+1

                    flx(i,j,k,n) = flx(i,j,k,n) + &
                         avisc_coeff * avisx(i,j,k) * (uin(i,j,k,n) - uin(i-1,j,k,n))
                 enddo
              enddo
           enddo
#if AMREX_SPACEDIM >= 2
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)+1
                 do i = lo(1), hi(1)

                    fly(i,j,k,n) = fly(i,j,k,n) + &
                         avisc_coeff * avisy(i,j,k) * (uin(i,j,k,n) - uin(i,j-1,k,n))
                 enddo
              enddo
           enddo
#endif
#if AMREX_SPACEDIM == 3
           do k = lo(3), hi(3)+1
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    flz(i,j,k,n) = flz(i,j,k,n) + &
                         avisc_coeff * avisz(i,j,k) * (uin(i,j,k,n) - uin(i,j,k-1,n))
                 enddo
              enddo
           enddo
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

  else
     ! do_hydro = 0
     flx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = ZERO
     qx_core_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = ZERO
     qx_pass_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = ZERO
#if AMREX_SPACEDIM >= 2
     fly(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),:) = ZERO
     qy_core_avg(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),:) = ZERO
     qy_pass_avg(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),:) = ZERO
#endif
#if AMREX_SPACEDIM == 3
     flz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,:) = ZERO
     qz_core_avg(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,:) = ZERO
     qz_pass_avg(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,:) = ZERO
#endif

  end if

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
                 update(i,j,k,UMX) = update(i,j,k,UMX) - &
                      ( qx_core_avg(i+1,j,k,QPRES) - qx_core_avg(i,j,k,QPRES) ) / dx(1)
              end if
#endif

#if AMREX_SPACEDIM == 2
              if (n == UMX) then
                 ! add the pressure source term for axisymmetry
                 if (coord_type > 0) then
                    update(i,j,k,n) = update(i,j,k,n) - (qx_core_avg(i+1,j,k,QPRES) - qx_core_avg(i,j,k,QPRES))/ dx(1)
                 end if
              end if
#endif

              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

           end do
        end do
     end do
  end do

#if AMREX_SPACEDIM == 3
#ifdef HYBRID_MOMENTUM
  call bl_allocate(qgdnvx, qc_lo, qc_hi, NGDNV)
  call bl_allocate(qgdnvy, qc_lo, qc_hi, NGDNV)
  call bl_allocate(qgdnvz, qc_lo, qc_hi, NGDNV)

  call ca_store_godunov_state(lo, hi+dg, &
                              qx_core_avg, qc_lo, qc_hi, &
                              qx_pass_avg, qc_lo, qc_hi, &
                              qgdnvx, qc_lo, qc_hi)

  call ca_store_godunov_state(lo, hi+dg, &
                              qy_core_avg, qc_lo, qc_hi, &
                              qy_pass_avg, qc_lo, qc_hi, &
                              qgdnvy, qc_lo, qc_hi)

  call ca_store_godunov_state(lo, hi+dg, &
                              qz_core_avg, qc_lo, qc_hi, &
                              qz_pass_avg, qc_lo, qc_hi, &
                              qgdnvz, qc_lo, qc_hi)

  call add_hybrid_advection_source(lo, hi, dt, &
                                   update, uout_lo, uout_hi, &
                                   qgdnvx, flx_lo, flx_hi, &
                                   qgdnvy, fly_lo, fly_hi, &
                                   qgdnvz, flz_lo, flz_hi)
  call bl_deallocate(qgdnvx)
  call bl_deallocate(qgdnvy)
  call bl_deallocate(qgdnvz)
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
                 flx(i,j,k,n) = flx(i,j,k,n) + &
                      dt * area1(i,j,k) * qx_core_avg(i,j,k,QPRES)
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
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = &
          qx_core_avg(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),QPRES) * dt
  end if
#endif
  call bl_deallocate(avisx)
#if BL_SPACEDIM >= 2
  call bl_deallocate(avisy)
#endif
#if BL_SPACEDIM == 3
  call bl_deallocate(avisz)
#endif

  call bl_deallocate(qx_core_avg)
  call bl_deallocate(qx_pass_avg)
  call bl_deallocate(flx_avg)
#if AMREX_SPACEDIM >= 2
  call bl_deallocate(qy_core_avg)
  call bl_deallocate(qy_pass_avg)
  call bl_deallocate(fly_avg)
#endif
#if AMREX_SPACEDIM == 3
  call bl_deallocate(qz_core_avg)
  call bl_deallocate(qz_pass_avg)
  call bl_deallocate(flz_avg)
#endif
#else
#ifndef AMREX_USE_CUDA
  ! RADIATION check
  call amrex_error("ERROR: ca_fourth_single_stage does not support radiation")
#endif
#endif
end subroutine ca_fourth_single_stage
