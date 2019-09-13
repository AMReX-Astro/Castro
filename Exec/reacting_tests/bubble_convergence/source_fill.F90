module source_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

contains

  subroutine source_single_fill(lo, hi, state, s_lo, s_hi, domlo, domhi, delta, xlo, bc) bind(C, name="source_single_fill")
    ! Used for a source fill of any StateData.

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here -- user's can override as needed
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. s_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. s_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, state, s_lo, s_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine source_single_fill



  subroutine source_multi_fill(lo, hi, state, s_lo, s_hi, domlo, domhi, delta, xlo, bc) bind(C, name="source_multi_fill")

    use meth_params_module, only : UMY, const_grav, sdc_order
    use amrex_filcc_module, only: amrex_filccn
    use model_parser_module
    use amrex_constants_module, only : HALF, ONE, TWO
    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: d, n
    integer :: bc_temp(AMREX_SPACEDIM,2,NVAR)
    real(rt) :: y, yp, ym, rhoc, rhop, rhom, rho_avg
    integer :: i, j, k
    integer :: jmin, jmax

    !$gpu


    ! handle an external BC via extrapolation here -- user's can override as needed
    bc_temp(:,:,:) = bc(:,:,:)

    do n = 1, NVAR
       do d = 1, AMREX_SPACEDIM
          if (bc(d,1,n) == EXT_DIR .and. s_lo(d) < domlo(d)) then
             bc_temp(d,1,n) = FOEXTRAP
          end if

          if (bc(d,2,n) == EXT_DIR .and. s_hi(d) > domhi(d)) then
             bc_temp(d,2,n) = FOEXTRAP
          end if
       end do
    end do

    call amrex_filccn(lo, hi, state, s_lo, s_hi, NVAR, domlo, domhi, delta, xlo, bc_temp)

    ! if we are UMY, we want to fill the gravity source term
    ! explicitly so we maintain HSE.
    ! YLO

    if (lo(2) < domlo(2) .and. bc(2,1,1) == EXT_DIR) then
       jmin = s_lo(2)
       jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
       ! For CUDA, this should only be one thread doing the work:
       ! we'll arbitrary choose the zone with index domlo(2) - 1.
       if (hi(2) /= jmax) then
          jmax = jmin - 1
       end if
#endif
       do j = jmin, jmax
          yp = problo(2) + delta(2)*(dble(j+1) + HALF)
          y = problo(2) + delta(2)*(dble(j) + HALF)
          ym = problo(2) + delta(2)*(dble(j-1) + HALF)

          call interpolate_sub(rhop, yp, idens_model)
          call interpolate_sub(rhoc, y, idens_model)
          call interpolate_sub(rhom, ym, idens_model)

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                if (sdc_order == 4) then
                   rho_avg = rhoc + (ONE/24.0_rt)*(rhop - TWO*rhoc + rhom)
                   state(i,j,k,UMY) = rho_avg * const_grav
                else
                   state(i,j,k,UMY) = rhoc * const_grav
                end if

             end do
          end do
       end do
    end if

    ! YHI
    if (hi(2) > domhi(2) .and. bc(2,2,1) == EXT_DIR) then
       jmin = domhi(2)+1
       jmax = s_hi(2)
#ifdef AMREX_USE_CUDA
       if (lo(2) /= jmin) then
          jmin = jmax + 1
       end if
#endif
       do j = jmin, jmax
          yp = problo(2) + delta(2)*(dble(j+1) + HALF)
          y = problo(2) + delta(2)*(dble(j) + HALF)
          ym = problo(2) + delta(2)*(dble(j-1) + HALF)

          call interpolate_sub(rhop, yp, idens_model)
          call interpolate_sub(rhoc, y, idens_model)
          call interpolate_sub(rhom, ym, idens_model)

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                if (sdc_order == 4) then
                   rho_avg = rhoc + (ONE/24.0_rt)*(rhop - TWO*rhoc + rhom)
                   state(i,j,k,UMY) = rho_avg * const_grav
                else
                   state(i,j,k,UMY) = rhoc * const_grav
                end if

             end do
          end do
       end do
    end if

  end subroutine source_multi_fill

end module source_fill_module
