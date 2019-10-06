! advection routines in support of method of lines integration
!

module castro_mol_module

  implicit none

contains

  subroutine ca_mol_plm_reconstruct(lo, hi, &
                                    q, q_lo, q_hi, &
                                    flatn, fl_lo, fl_hi, &
                                    dq, dq_lo, dq_hi, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    dx, domlo, domhi) bind(C, name="ca_mol_plm_reconstruct")

    use castro_error_module
    use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                   UTEMP, UMX, &
                                   plm_well_balanced, QPRES, &
                                   QTEMP, QFS, QFX, QREINT, QRHO, QU, QV, QW, &
                                   first_order_hydro, hybrid_riemann, &
                                   ppm_temp_fix, const_grav
    use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
    use slope_module, only : uslope
    use amrex_fort_module, only : rt => amrex_real
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux
    use prob_params_module, only : dg, coord_type, Symmetry, physbc_lo, physbc_hi

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: fl_lo(3), fl_hi(3)
    integer, intent(in) :: dq_lo(3), dq_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(in) :: flatn(fl_lo(1):fl_hi(1), fl_lo(2):fl_hi(2), fl_lo(3):fl_hi(3))
    real(rt), intent(inout) :: dq(dq_lo(1):dq_hi(1), dq_lo(2):dq_hi(2), dq_lo(3):dq_hi(3), NQ)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ, AMREX_SPACEDIM)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ, AMREX_SPACEDIM)
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(in) :: dx(3)

    integer :: idir, i, j, k, n
    type (eos_t) :: eos_state

    !$gpu


    do idir = 1, AMREX_SPACEDIM

       do n = 1, NQ
          ! piecewise linear slopes
          call uslope(lo, hi, idir, &
                      q, q_lo, q_hi, n, &
                      flatn, fl_lo, fl_hi, &
                      dq, dq_lo, dq_hi, &
                      dx, domlo, domhi)

       end do

       do n = 1, NQ

          ! for each slope, fill the two adjacent edge states
          if (idir == 1) then

             if (plm_well_balanced == 1 .and. n == QPRES .and. idir == AMREX_SPACEDIM) then

                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      do i = lo(1), hi(1)

                         ! left state at i+1/2 interface
                         qm(i+1,j,k,n,1) = q(i,j,k,n) + HALF*dq(i,j,k,n) + &
                              HALF*dx(1)*q(i,j,k,QRHO)*const_grav

                         ! right state at i-1/2 interface
                         qp(i,j,k,n,1) = q(i,j,k,n) - HALF*dq(i,j,k,n) - &
                              HALF*dx(1)*q(i,j,k,QRHO)*const_grav
                      end do
                   end do
                end do

             else
                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      do i = lo(1), hi(1)

                         ! left state at i+1/2 interface
                         qm(i+1,j,k,n,1) = q(i,j,k,n) + HALF*dq(i,j,k,n)

                         ! right state at i-1/2 interface
                         qp(i,j,k,n,1) = q(i,j,k,n) - HALF*dq(i,j,k,n)

                      end do
                   end do
                end do
             end if

#if BL_SPACEDIM >= 2
          else if (idir == 2) then

             if (plm_well_balanced == 1 .and. n == QPRES .and. idir == AMREX_SPACEDIM) then

                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      do i = lo(1), hi(1)

                         ! left state at i+1/2 interface
                         qm(i,j+1,k,n,2) = q(i,j,k,n) + HALF*dq(i,j,k,n) + &
                              HALF*dx(2)*q(i,j,k,QRHO)*const_grav

                         ! right state at i-1/2 interface
                         qp(i,j,k,n,2) = q(i,j,k,n) - HALF*dq(i,j,k,n) - &
                              HALF*dx(2)*q(i,j,k,QRHO)*const_grav
                      end do
                   end do
                end do

             else

                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      do i = lo(1), hi(1)

                         ! left state at j+1/2 interface
                         qm(i,j+1,k,n,2) = q(i,j,k,n) + HALF*dq(i,j,k,n)

                         ! right state at j-1/2 interface
                         qp(i,j,k,n,2) = q(i,j,k,n) - HALF*dq(i,j,k,n)

                      end do
                   end do
                end do
             end if
#endif

#if BL_SPACEDIM == 3
          else

             if (plm_well_balanced == 1 .and. n == QPRES .and. idir == AMREX_SPACEDIM) then

                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      do i = lo(1), hi(1)

                         ! left state at i+1/2 interface
                         qm(i,j,k+1,n,3) = q(i,j,k,n) + HALF*dq(i,j,k,n) + &
                              HALF*dx(3)*q(i,j,k,QRHO)*const_grav

                         ! right state at i-1/2 interface
                         qp(i,j,k,n,3) = q(i,j,k,n) - HALF*dq(i,j,k,n) - &
                              HALF*dx(3)*q(i,j,k,QRHO)*const_grav
                      end do
                   end do
                end do

             else

                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      do i = lo(1), hi(1)

                         ! left state at k+1/2 interface
                         qm(i,j,k+1,n,3) = q(i,j,k,n) + HALF*dq(i,j,k,n)

                         ! right state at k-1/2 interface
                         qp(i,j,k,n,3) = q(i,j,k,n) - HALF*dq(i,j,k,n)

                      end do
                   end do
                end do

             end if
#endif

          end if

       end do ! component loop

       ! special care for reflecting BCs
       if (idir == 1) then
          if (lo(1) == domlo(1) .and. physbc_lo(1) == Symmetry) then
             ! reset the left state at domlo(1) if needed -- it is outside the domain

             do k = lo(3)-dg(3), hi(3)+dg(3)
                do j = lo(2)-dg(2), hi(2)+dg(2)

                   qm(domlo(1),j,k,:,1) = qp(domlo(1),j,k,:,1)
                   qm(domlo(1),j,k,QU,1) = -qp(domlo(1),j,k,QU,1)

                end do
             end do

          end if

          if (hi(1)+1 == domhi(1)+1 .and. physbc_hi(1) == Symmetry) then
             ! reset the right state at domhi(1)+1 if needed -- it is outside the domain

             do k = lo(3)-dg(3), hi(3)+dg(3)
                do j = lo(2)-dg(2), hi(2)+dg(2)

                   qp(domhi(1)+1,j,k,:,1) = qm(domhi(1)+1,j,k,:,1)
                   qp(domhi(1)+1,j,k,QU,1) = -qm(domhi(1)+1,j,k,QU,1)

                end do
             end do

          end if

#if AMREX_SPACEDIM >= 2
       else if (idir == 2) then

          if (lo(2) == domlo(2) .and. physbc_lo(2) == Symmetry) then
             ! reset the left state at domlo(2) if needed -- it is outside the domain

             do k = lo(3)-dg(3), hi(3)+dg(3)
                do i = lo(1)-1, hi(1)+1

                   qm(i,domlo(2),k,:,2) = qp(i,domlo(2),k,:,2)
                   qm(i,domlo(2),k,QV,2) = -qp(i,domlo(2),k,QV,2)

                end do
             end do

          end if

          if (hi(2)+1 == domhi(2)+1 .and. physbc_hi(2) == Symmetry) then
             ! reset the right state at domhi(2)+1 if needed -- it is outside the domain

             do k = lo(3)-dg(3), hi(3)+dg(3)
                do i = lo(1)-1, hi(1)+1

                   qp(i,domhi(2)+1,k,:,2) = qm(i,domhi(2)+1,k,:,2)
                   qp(i,domhi(2)+1,k,QV,2) = -qm(i,domhi(2)+1,k,QV,2)

                end do
             end do

          end if
#endif
#if AMREX_SPACEDIM == 3
       else

          if (lo(3) == domlo(3) .and. physbc_lo(3) == Symmetry) then
             ! reset the left state at domlo(3) if needed -- it is outside the domain

             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   qm(i,j,domlo(3),:,3) = qp(i,j,domlo(3),:,3)
                   qm(i,j,domlo(3),QW,3) = -qp(i,j,domlo(3),QW,3)

                end do
             end do

          end if

          if (hi(3)+1 == domhi(3)+1 .and. physbc_hi(3) == Symmetry) then
             ! reset the right state at domhi(3)+1 if needed -- it is outside the domain

             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1

                   qp(i,j,domhi(3)+1,:,3) = qm(i,j,domhi(3)+1,:,3)
                   qp(i,j,domhi(3)+1,QW,3) = -qm(i,j,domhi(3)+1,QW,3)

                end do
             end do

          end if
#endif

       end if

    end do ! dimension loop

    ! use T to define p
    if (ppm_temp_fix == 1) then
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
    end if


  end subroutine ca_mol_plm_reconstruct

  subroutine ca_mol_ppm_reconstruct(lo, hi, &
                                    q, q_lo, q_hi, &
                                    flatn, fl_lo, fl_hi, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    dx) bind(C, name="ca_mol_ppm_reconstruct")

    use castro_error_module
    use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                   UTEMP, UMX, &
                                   QPRES, &
                                   QTEMP, QFS, QFX, QREINT, QRHO, &
                                   first_order_hydro, hybrid_riemann, &
                                   ppm_temp_fix
    use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
    use ppm_module, only: ppm_reconstruct
    use amrex_fort_module, only : rt => amrex_real
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux
    use prob_params_module, only : dg, coord_type

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: fl_lo(3), fl_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(in) :: flatn(fl_lo(1):fl_hi(1), fl_lo(2):fl_hi(2), fl_lo(3):fl_hi(3))
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ, AMREX_SPACEDIM)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ, AMREX_SPACEDIM)
    real(rt), intent(in) :: dx(3)

    integer :: idir, i, j, k, n
    type (eos_t) :: eos_state

    real(rt) :: s(-2:2)
    real(rt) :: sm, sp

    !$gpu

    do idir = 1, AMREX_SPACEDIM

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                do n = 1, NQ

                   if (idir == 1) then
                      s(:) = q(i-2:i+2,j,k,n)
                   else if (idir == 2) then
                      s(:) = q(i,j-2:j+2,k,n)
                   else
                      s(:) = q(i,j,k-2:k+2,n)
                   end if

                   call ppm_reconstruct(s, flatn(i,j,k), sm, sp)

                   if (idir == 1) then
                      ! right state at i-1/2
                      qp(i,j,k,n,1) = sm

                      ! left state at i+1/2
                      qm(i+1,j,k,n,1) = sp

                   else if (idir == 2) then
                      ! right state at j-1/2
                      qp(i,j,k,n,2) = sm

                      ! left state at j+1/2
                      qm(i,j+1,k,n,2) = sp

                   else
                      ! right state at k-1/2
                      qp(i,j,k,n,3) = sm

                      ! left state at k+1/2
                      qm(i,j,k+1,n,3) = sp

                   end if

                end do

             end do
          end do
       end do

    end do

    ! use T to define p
    if (ppm_temp_fix == 1) then

       do idir = 1, AMREX_SPACEDIM

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   ! we just got the extremes corresponding to a particular cell-center, but now
                   ! we need to assign them to interfaces

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

    use castro_error_module
    use meth_params_module, only : NQ, NVAR, NGDNV, GDPRES, &
                                   UTEMP, UMX, &
                                   QPRES, &
                                   QTEMP, QFS, QFX, QREINT, QRHO, &
                                   first_order_hydro, difmag, hybrid_riemann, &
                                   limit_fluxes_on_small_dens, ppm_type, ppm_temp_fix
    use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
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

  subroutine ca_mol_diffusive_flux(lo, hi, &
                                   idir, &
                                   U, U_lo, U_hi, &
                                   cond, c_lo, c_hi, &
                                   flux, f_lo, f_hi, &
                                   dx) bind(C, name="ca_mol_diffusive_flux")

    use castro_error_module
    use meth_params_module, only : NVAR, UTEMP, UEINT, UEDEN
    use amrex_constants_module, only : HALF
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: U_lo(3), U_hi(3)
    integer, intent(in) :: c_lo(3), c_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in), value :: idir

    real(rt), intent(in) :: U(U_lo(1):U_hi(1), U_lo(2):U_hi(2), U_lo(3):U_hi(3), NVAR)
    real(rt), intent(in) :: cond(c_lo(1):c_hi(1), c_lo(2):c_hi(2), c_lo(3):c_hi(3))
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in) :: dx(3)

    integer :: i, j, k

    real(rt) :: diff_term, cond_int

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (idir == 1) then
                cond_int = HALF * (cond(i,j,k) + cond(i-1,j,k))
                diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i-1,j,k,UTEMP))/dx(1)

             else if (idir == 2) then
                cond_int = HALF * (cond(i,j,k) + cond(i,j-1,k))
                diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i,j-1,k,UTEMP))/dx(2)

             else
                cond_int = HALF * (cond(i,j,k) + cond(i,j,k-1))
                diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i,j,k-1,UTEMP))/dx(3)
             end if

             flux(i,j,k,UEINT) = flux(i,j,k,UEINT) + diff_term
             flux(i,j,k,UEDEN) = flux(i,j,k,UEDEN) + diff_term

          end do
       end do
    end do

  end subroutine ca_mol_diffusive_flux

end module castro_mol_module
