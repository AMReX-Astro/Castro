! advection routines in support of method of lines integration
!

subroutine ca_mol_plm_reconstruct(lo, hi, &
                                  q_core, qc_lo, qc_hi, &
                                  q_pass, qp_lo, qp_hi, &
                                  flatn, fl_lo, fl_hi, &
                                  shk, shk_lo, shk_hi, &
                                  dq_core, dqc_lo, dqc_hi, &
                                  dq_pass, dqp_lo, dqp_hi, &
                                  qm_core, qmc_lo, qmc_hi, &
                                  qm_pass, qmp_lo, qmp_hi, &
                                  qp_core, qpc_lo, qpc_hi, &
                                  qp_pass, qpp_lo, qpp_hi, &
                                  dx) bind(C, name="ca_mol_plm_reconstruct")

  use amrex_error_module
  use meth_params_module, only : NQC, NQP, NVAR, NGDNV, GDPRES, &
                                 UTEMP, USHK, UMX, &
                                 use_flattening, QPRES, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, hybrid_riemann, &
                                 ppm_temp_fix
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use slope_module, only : uslope
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg, coord_type
  use advection_util_module, only : ca_shock

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: qc_lo(3), qc_hi(3)
  integer, intent(in) :: qp_lo(3), qp_hi(3)
  integer, intent(in) :: fl_lo(3), fl_hi(3)
  integer, intent(in) :: shk_lo(3), shk_hi(3)
  integer, intent(in) :: dqc_lo(3), dqc_hi(3)
  integer, intent(in) :: dqp_lo(3), dqp_hi(3)
  integer, intent(in) :: qm_core_lo(3), qm_core_hi(3)
  integer, intent(in) :: qm_pass_lo(3), qm_pass_hi(3)
  integer, intent(in) :: qp_core_lo(3), qp_core_hi(3)
  integer, intent(in) :: qp_pass_lo(3), qp_pass_hi(3)

  real(rt), intent(inout) :: q_core(qc_lo(1):qc_hi(1), qc_lo(2):qc_hi(2), qc_lo(3):qc_hi(3), NQC)
  real(rt), intent(inout) :: q_pass(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQP)
  real(rt), intent(in) :: flatn(fl_lo(1):fl_hi(1), fl_lo(2):fl_hi(2), fl_lo(3):fl_hi(3))
  real(rt), intent(inout) :: shk(shk_lo(1):shk_hi(1), shk_lo(2):shk_hi(2), shk_lo(3):shk_hi(3))
  real(rt), intent(inout) :: dq_core(dqc_lo(1):dqc_hi(1), dqc_lo(2):dqc_hi(2), dqc_lo(3):dqc_hi(3), NQC)
  real(rt), intent(inout) :: dq_pass(dqp_lo(1):dqp_hi(1), dqp_lo(2):dqp_hi(2), dqp_lo(3):dqp_hi(3), NQP)
  real(rt), intent(inout) :: qm_core(qmc_lo(1):qmc_hi(1), qmc_lo(2):qmc_hi(2), qmc_lo(3):qmc_hi(3), NQC, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qm_pass(qmp_lo(1):qmp_hi(1), qmp_lo(2):qmp_hi(2), qmp_lo(3):qmp_hi(3), NQP, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qp_core(qpc_lo(1):qpc_hi(1), qpc_lo(2):qpc_hi(2), qpc_lo(3):qpc_hi(3), NQC, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qp_pass(qpp_lo(1):qpp_hi(1), qpp_lo(2):qpp_hi(2), qpp_lo(3):qpp_hi(3), NQP, AMREX_SPACEDIM)
  real(rt), intent(in) :: dx(3)


  integer :: idir, i, j, k, n
  logical :: compute_shock

  !$gpu

#ifdef SHOCK_VAR
  compute_shock = .true.
#else
  compute_shock = .false.
#endif

  if (hybrid_riemann == 1 .or. compute_shock) then
     call ca_shock(lo, hi, &
                   q_core, qc_lo, qc_hi, &
                   shk, shk_lo, shk_hi, &
                   dx)
  else
     shk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
  endif

  do idir = 1, AMREX_SPACEDIM


     ! core variables

     do n = 1, NQC
        ! piecewise linear slopes
        call uslope(lo, hi, idir, &
                    q_core, qc_lo, qc_hi, n, &
                    flatn, fl_lo, fl_hi, &
                    dq_core, dqc_lo, dqc_hi)

     end do

     do n = 1, NQC

        ! for each slope, fill the two adjacent edge states
        if (idir == 1) then
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    ! left state at i+1/2 interface
                    qm_core(i+1,j,k,n,1) = q_core(i,j,k,n) + HALF*dq_core(i,j,k,n)

                    ! right state at i-1/2 interface
                    qp_core(i,j,k,n,1) = q_core(i,j,k,n) - HALF*dq_core(i,j,k,n)

                 end do
              end do
           end do

#if BL_SPACEDIM >= 2
        else if (idir == 2) then
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    ! left state at j+1/2 interface
                    qm_core(i,j+1,k,n,2) = q_core(i,j,k,n) + HALF*dq_core(i,j,k,n)

                    ! right state at j-1/2 interface
                    qp_core(i,j,k,n,2) = q_core(i,j,k,n) - HALF*dq_core(i,j,k,n)

                 end do
              end do
           end do
#endif

#if BL_SPACEDIM == 3
        else

           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    ! left state at k+1/2 interface
                    qm_core(i,j,k+1,n,3) = q_core(i,j,k,n) + HALF*dq_core(i,j,k,n)

                    ! right state at k-1/2 interface
                    qp_core(i,j,k,n,3) = q_core(i,j,k,n) - HALF*dq_core(i,j,k,n)

                 end do
              end do
           end do
#endif

        end if
     end do

     ! passive variables

     do n = 1, NQP
        ! piecewise linear slopes
        call uslope(lo, hi, idir, &
                    q_pass, qp_lo, qp_hi, n, &
                    flatn, fl_lo, fl_hi, &
                    dq_pass, dqp_lo, dqp_hi)

     end do

     do n = 1, NQP

        ! for each slope, fill the two adjacent edge states
        if (idir == 1) then
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    ! left state at i+1/2 interface
                    qm_pass(i+1,j,k,n,1) = q_pass(i,j,k,n) + HALF*dq_pass(i,j,k,n)

                    ! right state at i-1/2 interface
                    qp_pass(i,j,k,n,1) = q_pass(i,j,k,n) - HALF*dq_pass(i,j,k,n)

                 end do
              end do
           end do

#if BL_SPACEDIM >= 2
        else if (idir == 2) then
           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    ! left state at j+1/2 interface
                    qm_pass(i,j+1,k,n,2) = q_pass(i,j,k,n) + HALF*dq_pass(i,j,k,n)

                    ! right state at j-1/2 interface
                    qp_pass(i,j,k,n,2) = q_pass(i,j,k,n) - HALF*dq_pass(i,j,k,n)

                 end do
              end do
           end do
#endif

#if BL_SPACEDIM == 3
        else

           do k = lo(3), hi(3)
              do j = lo(2), hi(2)
                 do i = lo(1), hi(1)

                    ! left state at k+1/2 interface
                    qm_pass(i,j,k+1,n,3) = q_pass(i,j,k,n) + HALF*dq_pass(i,j,k,n)

                    ! right state at k-1/2 interface
                    qp_pass(i,j,k,n,3) = q_pass(i,j,k,n) - HALF*dq_pass(i,j,k,n)

                 end do
              end do
           end do
#endif

        end if

     end do
  end do

end subroutine ca_mol_plm_reconstruct


subroutine ca_mol_thermo_states(lo, hi, &
                                qm_core, qmc_lo, qmc_hi, &
                                qm_pass, qmp_lo, qmp_hi, &
                                qp_core, qpc_lo, qpc_hi, &
                                qp_pass, qpp_lo, qpp_hi) &
                                bind(C, name="ca_mol_thermo_states")

  use amrex_error_module
  use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                 UTEMP, USHK, UMX, &
                                 use_flattening, QPRES, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, hybrid_riemann, &
                                 ppm_temp_fix
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec, naux
  use prob_params_module, only : dg, coord_type

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: qmc_lo(3), qmc_hi(3)
  integer, intent(in) :: qmp_lo(3), qmp_hi(3)
  integer, intent(in) :: qpc_lo(3), qpc_hi(3)
  integer, intent(in) :: qpp_lo(3), qpp_hi(3)

  real(rt), intent(inout) :: qm_core(qmc_lo(1):qmc_hi(1), qmc_lo(2):qmc_hi(2), qmc_lo(3):qmc_hi(3), NQC, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qm_pass(qmp_lo(1):qmp_hi(1), qmp_lo(2):qmp_hi(2), qmp_lo(3):qmp_hi(3), NQP, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qp_core(qpc_lo(1):qpc_hi(1), qpc_lo(2):qpc_hi(2), qpc_lo(3):qpc_hi(3), NQC, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qp_pass(qpp_lo(1):qpp_hi(1), qpp_lo(2):qpp_hi(2), qpp_lo(3):qpp_hi(3), NQP, AMREX_SPACEDIM)


  integer :: idir, i, j, k, n
   type (eos_t) :: eos_state

  !$gpu

  do idir = 1, AMREX_SPACEDIM
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

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

end subroutine ca_mol_thermo_states

subroutine ca_mol_ppm_reconstruct(lo, hi, &
                                  q_core, qc_lo, qc_hi, &
                                  q_pass, qp_lo, qp_hi, &
                                  flatn, fl_lo, fl_hi, &
                                  shk, shk_lo, shk_hi, &
                                  qm_core, qmc_lo, qmc_hi, &
                                  qm_pass, qmp_lo, qmp_hi, &
                                  qp_core, qpc_lo, qpc_hi, &
                                  qp_pass, qpp_lo, qpp_hi, &
                                  dx) bind(C, name="ca_mol_ppm_reconstruct")

  use amrex_error_module
  use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                 UTEMP, USHK, UMX, &
                                 use_flattening, QPRES, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, hybrid_riemann
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use ppm_module, only : ca_ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : dg, coord_type
  use advection_util_module, only : ca_shock

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: qc_lo(3), qc_hi(3)
  integer, intent(in) :: qp_lo(3), qp_hi(3)
  integer, intent(in) :: fl_lo(3), fl_hi(3)
  integer, intent(in) :: shk_lo(3), shk_hi(3)
  integer, intent(in) :: qmc_lo(3), qmc_hi(3)
  integer, intent(in) :: qmp_lo(3), qmp_hi(3)
  integer, intent(in) :: qpc_lo(3), qpc_hi(3)
  integer, intent(in) :: qpp_lo(3), qpp_hi(3)

  real(rt), intent(inout) :: q_core(qc_lo(1):qc_hi(1), qc_lo(2):qc_hi(2), qc_lo(3):qc_hi(3), NQC)
  real(rt), intent(inout) :: q_pass(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQP)
  real(rt), intent(in) :: flatn(fl_lo(1):fl_hi(1), fl_lo(2):fl_hi(2), fl_lo(3):fl_hi(3))
  real(rt), intent(inout) :: shk(shk_lo(1):shk_hi(1), shk_lo(2):shk_hi(2), shk_lo(3):shk_hi(3))
  real(rt), intent(inout) :: qm_core(qmc_lo(1):qmc_hi(1), qmc_lo(2):qmc_hi(2), qmc_lo(3):qmc_hi(3), NQC, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qm_pass(qmp_lo(1):qmp_hi(1), qmp_lo(2):qmp_hi(2), qmp_lo(3):qmp_hi(3), NQP, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qp_core(qpc_lo(1):qpc_hi(1), qpc_lo(2):qpc_hi(2), qpc_lo(3):qpc_hi(3), NQC, AMREX_SPACEDIM)
  real(rt), intent(inout) :: qp_pass(qpp_lo(1):qpp_hi(1), qpp_lo(2):qpp_hi(2), qpp_lo(3):qpp_hi(3), NQP, AMREX_SPACEDIM)
  real(rt), intent(in) :: dx(3)

  integer :: idir, i, j, k, n
  logical :: compute_shock

  !$gpu

#ifdef SHOCK_VAR
  compute_shock = .true.
#else
  compute_shock = .false.
#endif

  if (hybrid_riemann == 1 .or. compute_shock) then
     call ca_shock(lo, hi, &
                   q_core, qc_lo, qc_hi, &
                   shk, shk_lo, shk_hi, &
                   dx)
  else
     shk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ZERO
  endif


  do idir = 1, AMREX_SPACEDIM
     call ca_ppm_reconstruct(lo, hi, 1, idir, &
                             q_core, qc_lo, qc_hi, NQC, 1, NQC, &
                             flatn, fl_lo, fl_hi, &
                             qm_core, qmc_lo, qmc_hi, &
                             qp_core, qpc_lo, qpc_hi, NQC, 1, NQC)

     call ca_ppm_reconstruct(lo, hi, 1, idir, &
                             q_pass, qp_lo, qp_hi, NQP, 1, NQP, &
                             flatn, fl_lo, fl_hi, &
                             qm_pass, qmp_lo, qmp_hi, &
                             qp_pass, qpp_lo, qpp_hi, NQP, 1, NQP)
  end do

end subroutine ca_mol_ppm_reconstruct

subroutine ca_mol_consup(lo, hi, &
                         uin, uin_lo, uin_hi, &
                         uout, uout_lo, uout_hi, &
                         srcU, srU_lo, srU_hi, &
                         update, updt_lo, updt_hi, &
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
                         q1, q1_lo, q1_hi, &
#if AMREX_SPACEDIM >= 2
                         q2, q2_lo, q2_hi, &
#endif
#if AMREX_SPACEDIM == 3
                         q3, q3_lo, q3_hi, &
#endif
                         vol, vol_lo, vol_hi) bind(C, name="ca_mol_consup")

  use amrex_error_module
  use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                 UTEMP, USHK, UMX, &
                                 use_flattening, QPRES, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens, ppm_type, ppm_temp_fix
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use slope_module, only : uslope
  use amrex_fort_module, only : rt => amrex_real
#ifdef HYBRID_MOMENTUM
  use hybrid_advection_module, only : add_hybrid_advection_source
#endif
  use network, only : nspec, naux
  use prob_params_module, only : dg, coord_type


  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
  integer, intent(in) :: q1_lo(3), q1_hi(3)
#if AMREX_SPACEDIM >= 2
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
  integer, intent(in) :: q2_lo(3), q2_hi(3)
#endif
#if AMREX_SPACEDIM == 3
  integer, intent(in) :: flux3_lo(3), flux3_hi(3)
  integer, intent(in) :: area3_lo(3), area3_hi(3)
  integer, intent(in) :: q3_lo(3), q3_hi(3)
#endif

  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
  real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)

  real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
  real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
  real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1), q1_lo(2):q1_hi(2), q1_lo(3):q1_hi(3), NGDNV)

#if AMREX_SPACEDIM >= 2
  real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
  real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1), q2_lo(2):q2_hi(2), q2_lo(3):q2_hi(3), NGDNV)
#endif

#if AMREX_SPACEDIM == 3
  real(rt), intent(inout) :: flux3(flux3_lo(1):flux3_hi(1), flux3_lo(2):flux3_hi(2), flux3_lo(3):flux3_hi(3), NVAR)
  real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1), area3_lo(2):area3_hi(2), area3_lo(3):area3_hi(3))
  real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1), q3_lo(2):q3_hi(2), q3_lo(3):q3_hi(3), NGDNV)
#endif

  real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt), intent(in) :: dx(3)
  real(rt), intent(in), value :: dt

  integer :: i, j, k, n

  !$gpu

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
              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

           enddo
        enddo
     enddo
  enddo

#if AMREX_SPACEDIM == 3
#ifdef HYBRID_MOMENTUM
  call add_hybrid_advection_source(lo, hi, dt, &
                                   update, uout_lo, uout_hi, &
                                   q1, q1_lo, q1_hi, &
                                   q2, q2_lo, q2_hi, &
                                   q3, q3_lo, q3_hi)
#endif
#endif

end subroutine ca_mol_consup
