module bc_ext_fill_module
    ! this module contains different routines for filling the
    ! hydrodynamics boundary conditions

    ! .. note::
    !    the hydrostatic boundary conditions here rely on
    !    constant gravity

  use amrex_constants_module, only: ZERO, HALF
#ifndef AMREX_USE_CUDA
  use castro_error_module, only: castro_error
#endif
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UFX, UTEMP
  use prob_params_module, only: dim

  implicit none

  include 'AMReX_bc_types.fi'

contains

#ifndef CXX_MODEL_PARSER
  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only: problo
    use model_parser_module, only: idens_model, interpolate_sub
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real(rt) :: x, y, z

    !$gpu

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    ! XLO
    if ( bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       imin = adv_lo(1)
       imax = domlo(1)-1
#ifdef AMREX_USE_CUDA
       ! For CUDA, this should only be one thread doing the work:
       ! we'll arbitrary choose the zone with index domlo(2) - 1.
       if (hi(1) /= imax) then
          imax = imin - 1
       end if
#endif
       do i = imin, imax
          x = problo(1) + delta(1)*(dble(i) + HALF)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                call interpolate_sub(adv(i,j,k), x, idens_model)
             end do
          end do
       end do
    end if

    ! XHI
    if ( bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       imin = domhi(1)+1
       imax = adv_hi(1)
#ifdef AMREX_USE_CUDA
       if (lo(1) /= imin) then
          imin = imax + 1
       end if
#endif
       do i = imin, imax
          x = problo(1) + delta(1)*(dble(i)+ HALF)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                call interpolate_sub(adv(i,j,k), x, idens_model)
             end do
          end do
       end do
    endif

#if AMREX_SPACEDIM >= 2
    ! YLO
    if ( bc(2,1) == EXT_DIR .and. lo(2) < domlo(2)) then
       jmin = adv_lo(2)
       jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
       ! For CUDA, this should only be one thread doing the work:
       ! we'll arbitrary choose the zone with index domlo(2) - 1.
       if (hi(2) /= jmax) then
          jmax = jmin - 1
       end if
#endif
       do j = jmin, jmax
          y = problo(2) + delta(2)*(dble(j) + HALF)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), y, idens_model)
             end do
          end do
       end do
    end if

    ! YHI
    if ( bc(2,2) == EXT_DIR .and. hi(2) > domhi(2)) then
       jmin = domhi(2)+1
       jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
       if (lo(2) /= jmin) then
          jmin = jmax + 1
       end if
#endif
       do j = jmin, jmax
          y = problo(2) + delta(2)*(dble(j)+ HALF)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), y, idens_model)
             end do
          end do
       end do
    end if
#endif

#if AMREX_SPACEDIM == 3
    ! ZLO
    if ( bc(3,1) == EXT_DIR .and. lo(3) < domlo(3)) then
       kmin = adv_lo(3)
       kmax = domlo(3)-1
#ifdef AMREX_USE_CUDA
       if (hi(3) /= kmax) then
          kmax = kmin - 1
       end if
#endif
       do k = kmin, kmax
          z = problo(3) + delta(3)*(dble(k) + HALF)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), z, idens_model)
             end do
          end do
       end do
    end if

    ! ZHI
    if ( bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then
       kmin = domhi(3)+1
       kmax = adv_hi(3)
#ifdef AMREX_USE_CUDA
       if (lo(3) /= kmin) then
          kmin = kmax + 1
       end if
#endif
       do k = kmin, kmax
          z = problo(3) + delta(3)*(dble(k)+ HALF)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                call interpolate_sub(adv(i,j,k), z, idens_model)
             end do
          end do
       end do
    end if
#endif

  end subroutine ext_denfill
#endif

end module bc_ext_fill_module
