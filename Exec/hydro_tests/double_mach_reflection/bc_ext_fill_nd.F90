module bc_ext_fill_module
    ! this module contains different routines for filling the
    ! hydrodynamics boundary conditions

    ! .. note::
    !    the hydrostatic boundary conditions here rely on
    !    constant gravity

  use amrex_constants_module, only: ZERO, HALF, THREE, M_PI, FOURTH, SIXTH
#ifndef AMREX_USE_CUDA
  use castro_error_module, only: castro_error
#endif
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, &
                                 hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                 xl_ext, xr_ext, yl_ext, yr_ext, zl_ext, zr_ext
  use prob_params_module, only: dim

  implicit none

  include 'AMReX_bc_types.fi'

contains


  subroutine ext_fill(lo, hi, adv, adv_lo, adv_hi, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo, dim
    use probdata_module
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec
    use model_parser_module, only: model_r, model_state, npts_model, idens_model, itemp_model, ispec_model
    use amrex_filcc_module, only: amrex_filccn

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k, ii, jj, n
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real(rt) :: x, y, z, xcen, ycen, shockfront

    real(rt), parameter :: pi_over_3 = M_PI / THREE
    real(rt), parameter :: ff = FOURTH

    !$gpu

    !-------------------------------------------------------------------------
    ! x boundaries
    !-------------------------------------------------------------------------

    if (bc(1,1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             ! integrate downward
             imin = adv_lo(1)
             imax = domlo(1)-1
#ifdef AMREX_USE_CUDA
             ! For CUDA, this should only be one thread doing the work:
             ! we'll arbitrary choose the zone with index domlo(2) - 1.
             if (hi(1) /= imax) then
                imax = imin - 1
             end if
#endif
             do i = imax, imin, -1
                adv(i,j,k,URHO) = rho_l
                adv(i,j,k,UMX) = rho_l*u_l
                adv(i,j,k,UMY) = rho_l*v_l
                adv(i,j,k,UMZ) = 0.0e0_rt
                adv(i,j,k,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*(u_l*u_l + v_l*v_l)
                adv(i,j,k,UEINT) = rhoe_l
                adv(i,j,k,UTEMP) = T_l
                adv(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
                adv(i,j,k,UFS) = adv(i,j,k,URHO)
             end do

          end do
       end do
    endif

#ifndef AMREX_USE_CUDA
    if (bc(1,2,1) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: special boundary not defined at +x")
    end if
#endif

#if AMREX_SPACEDIM >= 2
    !-------------------------------------------------------------------------
    ! y boundaries
    !-------------------------------------------------------------------------

    ! YLO
    if (bc(2,1,1) == EXT_DIR .and. lo(2) < domlo(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             x = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

             ! integrate downward
             jmin = adv_lo(2)
             jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
             ! For CUDA, this should only be one thread doing the work:
             ! we'll arbitrary choose the zone with index domlo(2) - 1.
             if (hi(2) /= jmax) then
                jmax = jmin - 1
             end if
#endif
             do j = jmax, jmin, -1

                if (x <  1.e0_rt/6.e0_rt) then
                   ! ICs
                   adv(i,j,k,URHO) = rho_l
                   adv(i,j,k,UMX) = rho_l*u_l
                   adv(i,j,k,UMY) = rho_l*v_l
                   adv(i,j,k,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*(u_l*u_l + v_l*v_l)
                   adv(i,j,k,UEINT) = rhoe_l
                   adv(i,j,k,UTEMP) = T_l
                else
                   ! reflect
                   adv(i,j,k,URHO) =  adv(i,domlo(2),k,URHO)
                   adv(i,j,k,UMX)  =  adv(i,domlo(2),k,UMX)
                   adv(i,j,k,UMY)  =  -1.e0_rt*adv(i,domlo(2),k,UMY)
                   adv(i,j,k,UEDEN) = adv(i,domlo(2),k,UEDEN)
                   adv(i,j,k,UEINT) = adv(i,domlo(2),k,UEINT)
                   adv(i,j,k,UTEMP) = adv(i,domlo(2),k,UTEMP)
                endif
                adv(i,j,k,UMZ) = 0.0e0_rt
                adv(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
                adv(i,j,k,UFS) = adv(i,j,k,URHO)

             end do
          end do
       end do
    endif

    ! YHI
    if (bc(2,2,1) == EXT_DIR .and. hi(2) > domhi(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             xcen = problo(1) + delta(1) * (dble(i) + HALF)

             jmin = domhi(2)+1
             jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
             if (lo(2) /= jmin) then
                jmin = jmax + 1
             end if
#endif
             do j = jmin, jmax
                ycen = problo(2) + delta(2) * (dble(j) + HALF)

                adv(i,j,k,URHO ) = 0.e0_rt
                adv(i,j,k,UMX  ) = 0.e0_rt
                adv(i,j,k,UMY  ) = 0.e0_rt
                adv(i,j,k,UEDEN) = 0.e0_rt
                adv(i,j,k,UEINT) = 0.e0_rt
                adv(i,j,k,UTEMP) = 0.e0_rt

                do jj = -1, 1
                   if (jj == 0) cycle
                   y = ycen + HALF * delta(2) * (jj / sqrt(THREE))

                   shockfront = sixth + y/tan(pi_over_3) + (10.e0_rt/sin(pi_over_3))*time

                   do ii = -1, 1
                      if (ii == 0) cycle
                      x = xcen + HALF * delta(1) * (ii / sqrt(THREE))

                      if (x < shockfront) then
                         ! Post shock ICs
                         adv(i,j,k,URHO ) = adv(i,j,k,URHO ) + ff*rho_l
                         adv(i,j,k,UMX  ) = adv(i,j,k,UMX  ) + ff*rho_l*u_l
                         adv(i,j,k,UMY  ) = adv(i,j,k,UMY  ) + ff*rho_l*v_l
                         adv(i,j,k,UEDEN) = adv(i,j,k,UEDEN) + ff*(rhoe_l + 0.5e0_rt*rho_l*(u_l*u_l + v_l*v_l))
                         adv(i,j,k,UEINT) = adv(i,j,k,UEINT) + ff*rhoe_l
                         adv(i,j,k,UTEMP) = adv(i,j,k,UTEMP) + ff*T_l
                      else
                         ! Pre Shock ICs
                         adv(i,j,k,URHO ) = adv(i,j,k,URHO ) + ff*rho_r
                         adv(i,j,k,UMX  ) = adv(i,j,k,UMX  ) + ff*rho_r*u_r
                         adv(i,j,k,UMY  ) = adv(i,j,k,UMY  ) + ff*rho_r*v_r
                         adv(i,j,k,UEDEN) = adv(i,j,k,UEDEN) + ff*(rhoe_r + 0.5e0_rt*rho_r*(u_r*u_r + v_r*v_r))
                         adv(i,j,k,UEINT) = adv(i,j,k,UEINT) + ff*rhoe_r
                         adv(i,j,k,UTEMP) = adv(i,j,k,UTEMP) + ff*T_r
                      end if

                   end do
                end do

                adv(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
                adv(i,j,k,UFS) = adv(i,j,k,URHO)
                adv(i,j,k,UMZ  ) = 0.0e0_rt

             end do
          end do
       end do
    end if

#endif

#ifndef AMREX_USE_CUDA
#if AMREX_SPACEDIM == 3

    !-------------------------------------------------------------------------
    ! z boundaries
    !-------------------------------------------------------------------------

    ! ZLO
    if (bc(3,1,1) == EXT_DIR .and. lo(3) < domlo(3)) then
       call castro_error("ERROR: -z special BCs not implemented")
    endif

    ! ZHI
    if (bc(3,2,1) == EXT_DIR .and. hi(3) > domhi(3)) then
       call castro_error("ERROR: +z special BCs not implemented")
    end if
#endif
#endif

  end subroutine ext_fill


  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) &
                         bind(C, name="ext_denfill")

    use probdata_module
    use prob_params_module, only: problo
    use model_parser_module, only: npts_model, model_r, model_state, idens_model
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use amrex_filcc_module, only: amrex_filccn

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k, ii, jj
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real(rt) :: x, y, z, xcen, ycen, shockfront

    real(rt), parameter :: pi_over_3 = M_PI / THREE
    real(rt), parameter :: ff = FOURTH

    !$gpu

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    !-------------------------------------------------------------------------
    ! x boundaries
    !-------------------------------------------------------------------------

    if (bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)

             ! integrate downward
             imin = adv_lo(1)
             imax = domlo(1)-1
#ifdef AMREX_USE_CUDA
             ! For CUDA, this should only be one thread doing the work:
             ! we'll arbitrary choose the zone with index domlo(2) - 1.
             if (hi(1) /= imax) then
                imax = imin - 1
             end if
#endif
             do i = imax, imin, -1
                adv(i,j,k) = rho_l
             end do

          end do
       end do
    endif

#ifndef AMREX_USE_CUDA
    if (bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: special boundary not defined at +x")
    end if
#endif

#if AMREX_SPACEDIM >= 2
    !-------------------------------------------------------------------------
    ! y boundaries
    !-------------------------------------------------------------------------

    ! YLO
    if (bc(2,1) == EXT_DIR .and. lo(2) < domlo(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             x = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

             ! integrate downward
             jmin = adv_lo(2)
             jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
             ! For CUDA, this should only be one thread doing the work:
             ! we'll arbitrary choose the zone with index domlo(2) - 1.
             if (hi(2) /= jmax) then
                jmax = jmin - 1
             end if
#endif
             do j = jmax, jmin, -1

                if (x <  1.e0_rt/6.e0_rt) then
                   ! ICs
                   adv(i,j,k) = rho_l
                else
                   ! reflect
                   adv(i,j,k) =  adv(i,domlo(2),k)
                endif
             end do
          end do
       end do
    endif

    ! YHI
    if (bc(2,2) == EXT_DIR .and. hi(2) > domhi(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             xcen = problo(1) + delta(1) * (dble(i) + HALF)

             jmin = domhi(2)+1
             jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
             if (lo(2) /= jmin) then
                jmin = jmax + 1
             end if
#endif
             do j = jmin, jmax
                ycen = problo(2) + delta(2) * (dble(j) + HALF)

                adv(i,j,k) = 0.e0_rt

                do jj = -1, 1
                   if (jj == 0) cycle
                   y = ycen + HALF * delta(2) * (jj / sqrt(THREE))

                   shockfront = sixth + y/tan(pi_over_3) + (10.e0_rt/sin(pi_over_3))*time

                   do ii = -1, 1
                      if (ii == 0) cycle
                      x = xcen + HALF * delta(1) * (ii / sqrt(THREE))

                      if (x < shockfront) then
                         ! Post shock ICs
                         adv(i,j,k) = adv(i,j,k) + ff*rho_l
                      else
                         ! Pre Shock ICs
                         adv(i,j,k) = adv(i,j,k) + ff*rho_r
                      end if

                   end do
                end do

             end do
          end do
       end do
    end if

#endif

#ifndef AMREX_USE_CUDA
#if AMREX_SPACEDIM == 3

    !-------------------------------------------------------------------------
    ! z boundaries
    !-------------------------------------------------------------------------

    ! ZLO
    if (bc(3,1) == EXT_DIR .and. lo(3) < domlo(3)) then
       call castro_error("ERROR: -z special BCs not implemented")
    endif

    ! ZHI
    if (bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then
       call castro_error("ERROR: +z special BCs not implemented")
    end if
#endif
#endif

  end subroutine ext_denfill

end module bc_ext_fill_module
