module bc_ext_fill_module
    ! this module contains different routines for filling the
    ! hydrodynamics boundary conditions

    ! .. note::
    !    the hydrostatic boundary conditions here rely on
    !    constant gravity

  use amrex_constants_module, only: ZERO, HALF, ONE, TWO
#ifndef AMREX_USE_CUDA
  use castro_error_module, only: castro_error
#endif
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, &
                                 hse_zero_vels, hse_interp_temp, hse_reflect_vels, &
                                 xl_ext, xr_ext, yl_ext, yr_ext, zl_ext,zr_ext, EXT_HSE, EXT_INTERP
  use prob_params_module, only: dim

  implicit none

  include 'AMReX_bc_types.fi'

contains


  subroutine ext_fill(lo, hi, adv, adv_lo, adv_hi, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo, dim
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec
    use model_parser_module, only: model_r, model_state, npts_model, idens_model, itemp_model, ispec_model, interpolate_sub
    use amrex_filcc_module, only: amrex_filccn
    use meth_params_module, only : sdc_order

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k, q, iter, m, d
    integer :: ioff, joff, koff
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real(rt) :: x, y, yp, ym, z

    real(rt) :: rhoc, rhop, rhom, rho_avg
    real(rt) :: Tc, Tp, Tm, T_avg
    real(rt) :: rhoec, rhoep, rhoem, rhoe_avg

    real(rt) :: X_zone(nspec)

    type (eos_t) :: eos_state

    !$gpu


    !-------------------------------------------------------------------------
    ! x boundaries
    !-------------------------------------------------------------------------

    ! XLO
    if (bc(1,1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       call castro_error("BCs not implemented for -X")
    endif


    ! XHI
    if (bc(1,2,1) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("BCs not implemented for +X")
    endif


#if AMREX_SPACEDIM >= 2
    !-------------------------------------------------------------------------
    ! y boundaries
    !-------------------------------------------------------------------------

    ! YLO
    if (bc(2,1,1) == EXT_DIR .and. lo(2) < domlo(2)) then

       jmin = adv_lo(2)
       jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
       if (hi(2) /= jmax) then
          jmax = jmin - 1
       end if
#endif
       do j = jmax, jmin, -1
          yp = problo(2) + delta(2)*(dble(j+1)+HALF)
          y = problo(2) + delta(2)*(dble(j)+HALF)
          ym = problo(2) + delta(2)*(dble(j-1)+HALF)

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                call interpolate_sub(rhop, yp, idens_model)
                call interpolate_sub(rhoc, y, idens_model)
                call interpolate_sub(rhom, ym, idens_model)

                call interpolate_sub(Tp, yp, itemp_model)
                call interpolate_sub(Tc, y, itemp_model)
                call interpolate_sub(Tm, ym, itemp_model)

                ! the composition in our initial model is uniform, so we just need
                ! it at one point to get the average
                do q = 1, nspec
                   call interpolate_sub(X_zone(q), y, ispec_model-1+q)
                enddo

                ! extrap normal momentum
                adv(i,j,k,UMY) = min(ZERO, adv(i,domlo(2),k,UMY))

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMZ) = ZERO

                ! get the thermodynamics for the ym, y, yp
                eos_state % rho = rhop
                eos_state % T = Tp
                eos_state % xn(:) = X_zone
                call eos(eos_input_rt, eos_state)
                rhoep = rhop * eos_state % e

                eos_state % rho = rhoc
                eos_state % T = Tc
                eos_state % xn(:) = X_zone
                call eos(eos_input_rt, eos_state)
                rhoec = rhoc * eos_state % e

                eos_state % rho = rhom
                eos_state % T = Tm
                eos_state % xn(:) = X_zone
                call eos(eos_input_rt, eos_state)
                rhoem = rhom * eos_state % e

                ! now make the averages
                if (sdc_order == 4) then
                   rho_avg = rhoc + (ONE/24.0_rt) * (rhop - TWO*rhoc + rhom)
                   T_avg = Tc + (ONE/24.0_rt) * (Tp - TWO*Tc + Tm)
                   rhoe_avg = rhoec + (ONE/24.0_rt) * (rhoep - TWO*rhoec + rhoem)
                else
                   rho_avg = rhoc
                   T_avg = Tc
                   rhoe_avg = rhoec
                end if

                adv(i,j,k,URHO) = rho_avg
                adv(i,j,k,UEINT) = rhoe_avg
                adv(i,j,k,UTEMP) = T_avg
                adv(i,j,k,UFS:UFS-1+nspec) = rho_avg*X_zone(:)

                ! this is an approximation -- the KE term here is only
                ! second order accurate.  We are relying on the fact
                ! that for HSE, it should be much smaller than the
                ! internal energy
                adv(i,j,k,UEDEN) = rhoe_avg + &
                     HALF*sum(adv(i,j,k,UMX:UMZ)**2)/rhoc

             end do
          end do
       end do

    endif


    ! YHI
    if (bc(2,2,1) == EXT_DIR .and. hi(2) > domhi(2)) then

       ! interpolate thermodynamics from initial model

       jmin = domhi(2)+1
       jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
       if (lo(2) /= jmin) then
          jmin = jmax + 1
       end if
#endif
       do j = jmin, jmax
          yp = problo(2) + delta(2)*(dble(j+1)+HALF)
          y = problo(2) + delta(2)*(dble(j)+HALF)
          ym = problo(2) + delta(2)*(dble(j-1)+HALF)

          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                call interpolate_sub(rhop, yp, idens_model)
                call interpolate_sub(rhoc, y, idens_model)
                call interpolate_sub(rhom, ym, idens_model)

                call interpolate_sub(Tp, yp, itemp_model)
                call interpolate_sub(Tc, y, itemp_model)
                call interpolate_sub(Tm, ym, itemp_model)

                ! the composition in our initial model is uniform, so we just need
                ! it at one point to get the average
                do q = 1, nspec
                   call interpolate_sub(X_zone(q), y, ispec_model-1+q)
                enddo

                ! extrap normal momentum
                adv(i,j,k,UMY) = max(ZERO, adv(i,domhi(2),k,UMY))

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMZ) = ZERO

                ! get the thermodynamics for the ym, y, yp
                eos_state % rho = rhop
                eos_state % T = Tp
                eos_state % xn(:) = X_zone
                call eos(eos_input_rt, eos_state)
                rhoep = rhop * eos_state % e

                eos_state % rho = rhoc
                eos_state % T = Tc
                eos_state % xn(:) = X_zone
                call eos(eos_input_rt, eos_state)
                rhoec = rhoc * eos_state % e

                eos_state % rho = rhom
                eos_state % T = Tm
                eos_state % xn(:) = X_zone
                call eos(eos_input_rt, eos_state)
                rhoem = rhom * eos_state % e

                ! now make the averages
                if (sdc_order == 4) then
                   rho_avg = rhoc + (ONE/24.0_rt) * (rhop - TWO*rhoc + rhom)
                   T_avg = Tc + (ONE/24.0_rt) * (Tp - TWO*Tc + Tm)
                   rhoe_avg = rhoec + (ONE/24.0_rt) * (rhoep - TWO*rhoec + rhoem)
                else
                   rho_avg = rhoc
                   T_avg = Tc
                   rhoe_avg = rhoec
                end if

                adv(i,j,k,URHO) = rho_avg
                adv(i,j,k,UEINT) = rhoe_avg
                adv(i,j,k,UTEMP) = T_avg
                adv(i,j,k,UFS:UFS-1+nspec) = rho_avg*X_zone(:)

                ! this is an approximation -- the KE term here is only
                ! second order accurate.  We are relying on the fact
                ! that for HSE, it should be much smaller than the
                ! internal energy
                adv(i,j,k,UEDEN) = rhoe_avg + &
                     HALF*sum(adv(i,j,k,UMX:UMZ)**2)/rhoc

             end do
          end do
       end do

    endif
#endif

#if AMREX_SPACEDIM == 3

    !-------------------------------------------------------------------------
    ! z boundaries
    !-------------------------------------------------------------------------

    ! ZLO
    if (bc(3,1,1) == EXT_DIR .and. lo(3) < domlo(3)) then
       call castro_error("BCs not implemented for -Z")
    endif

    ! ZHI
    if (bc(3,2,1) == EXT_DIR .and. hi(3) > domhi(3)) then
       call castro_error("BCs not implemented for +Z")
    end if
#endif

  end subroutine ext_fill


  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only: problo
    use model_parser_module, only: npts_model, model_r, model_state, idens_model, interpolate_sub
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

    integer :: i, j, k
    integer :: jmin, jmax, kmin, kmax
    real(rt) :: y, z

    !$gpu

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

#ifndef AMREX_USE_CUDA
    ! XLO
    if ( bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       call castro_error("We should not be here (xlo denfill)")
    end if

    ! XHI
    if ( bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("We should not be here (xhi denfill)")
    endif
#endif

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

    end if

    ! ZHI
    if ( bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then

    end if
#endif

  end subroutine ext_denfill

#ifdef GRAVITY
  subroutine ext_gravxfill(lo, hi, grav, grav_lo, grav_hi, &
                           domlo, domhi, delta, xlo, time, bc) &
                           bind(C, name="ext_gravxfill")

    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time


    !$gpu

    ! this is currently a stub

  end subroutine ext_gravxfill


  subroutine ext_gravyfill(lo, hi, grav, grav_lo, grav_hi, &
                           domlo, domhi, delta, xlo, time, bc) &
                           bind(C, name="ext_gravyfill")

    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time


    !$gpu

    ! this is currently a stub

  end subroutine ext_gravyfill


  subroutine ext_gravzfill(lo, hi, grav, grav_lo, grav_hi, &
                           domlo, domhi, delta, xlo, time, bc) &
                           bind(C, name="ext_gravzfill")

    use prob_params_module, only: problo

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: grav_lo(3), grav_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3))
    real(rt), intent(in   ), value :: time


    !$gpu

    ! this is currently a stub

  end subroutine ext_gravzfill
#endif

end module bc_ext_fill_module
