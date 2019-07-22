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
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, T_guess
  use prob_params_module, only: dim

  implicit none

  include 'AMReX_bc_types.fi'

contains


  subroutine ext_fill(lo, hi, adv, adv_lo, adv_hi, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo, dim, center
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use network, only: nspec
    use probdata_module
    use model_module

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k, n
    integer :: jmin, jmax, kmin, kmax

    integer :: lo_model, hi_model

    real(rt), allocatable :: r_model(:), rho_model(:), T_model(:), &
                             e_model(:), p_model(:)

    ! compute background state.  Here we assume that the background is
    ! in the vertical direction.  These BCs only make sense for the
    ! boundary along the background direction (dimension
    ! AMREX_SPACEDIM).

    ! we'll generate the initial model at the needed resolution
    call get_model_size(rmin, rmax, delta(AMREX_SPACEDIM), lo_model, hi_model)

    allocate(  r_model(lo_model:hi_model))
    allocate(rho_model(lo_model:hi_model))
    allocate(  T_model(lo_model:hi_model))
    allocate(  e_model(lo_model:hi_model))
    allocate(  p_model(lo_model:hi_model))

    call get_model(rmin, rmax, delta(AMREX_SPACEDIM), &
                   pres_base, dens_base, do_isentropic, &
                   xn_model, &
                   r_model, rho_model, T_model, e_model, p_model, &
                   lo_model, hi_model)

    ! we'll just look at the density when determining how to fill the BCs for all vars
    n = URHO

#ifndef AMREX_USE_CUDA
    if (bc(1,1,n) == EXT_DIR .and. lo(1) < domlo(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for -x BC")
    endif

    if (bc(1,2,n) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for +x BC, d")
    end if
#endif

#if AMREX_SPACEDIM >= 2
    !-------------------------------------------------------------------------
    ! y boundaries
    !-------------------------------------------------------------------------

    ! YLO
    if (bc(2,1,n) == EXT_DIR .and. lo(2) < domlo(2)) then

       ! we will fill all the variables when we consider URHO
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

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

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMZ) = ZERO

                if (boundary_type .eq. 1) then
                   ! extrapolate normal momentum
                   ! enforces pi=0 at boundary
                   adv(i,j,k,UMY) = adv(i,domlo(2),k,UMY)
                else
                   ! zero normal momentum
                   ! permits pi to pass through boundary
                   adv(i,j,k,UMY) = ZERO
                end if

                adv(i,j,k,URHO) = rho_model(j)
                adv(i,j,k,UFS:UFS-1+nspec) = adv(i,j,k,URHO)*xn_model(:)
                adv(i,j,k,UEINT) = e_model(j)*adv(i,j,k,URHO)
                adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                           + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                adv(i,j,k,UTEMP) = T_model(j)

             end do
          end do
       end do

    endif


    ! YHI
    if (bc(2,2,n) == EXT_DIR .and. hi(2) > domhi(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             jmin = domhi(2)+1
             jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
             if (lo(2) /= jmin) then
                jmin = jmax + 1
             end if
#endif
             do j = jmin, jmax

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMZ) = ZERO

                if (boundary_type .eq. 1) then
                   ! extrapolate normal momentum
                   ! enforces pi=0 at boundary
                   adv(i,j,k,UMY) = adv(i,domhi(2),k,UMY)
                else
                         ! zero normal momentum
                   ! permits pi to pass through boundary
                   adv(i,j,k,UMY) = ZERO
                end if

                adv(i,j,k,URHO) = rho_model(j)
                adv(i,j,k,UFS:UFS-1+nspec) = adv(i,j,k,URHO)*xn_model(:)
                adv(i,j,k,UEINT) = e_model(j)*adv(i,j,k,URHO)
                adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                     + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                adv(i,j,k,UTEMP) = T_model(j)

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
    if (bc(3,1,n) == EXT_DIR .and. lo(3) < domlo(3)) then

       do j= lo(2), hi(2)
          do i = lo(1), hi(1)

             ! integrate downward
             kmin = adv_lo(3)
             kmax = domlo(3) - 1
#ifdef AMREX_USE_CUDA
             if (hi(3) /= kmax) then
                kmax = kmin - 1
             end if
#endif
             do k = kmax, kmin, -1

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMY) = ZERO

                if (boundary_type .eq. 1) then
                   ! extrapolate normal momentum
                   ! enforces pi=0 at boundary
                   adv(i,j,k,UMZ) = adv(i,j,domlo(3),UMZ)
                else
                   ! zero normal momentum
                   ! permits pi to pass through boundary
                   adv(i,j,k,UMZ) = ZERO
                end if

                adv(i,j,k,URHO) = density(k)
                adv(i,j,k,UFS:UFS-1+nspec) = adv(i,j,k,URHO)*xn_model(:)
                adv(i,j,k,UEINT) = e_model(k)*adv(i,j,k,URHO)
                adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                     + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                adv(i,j,k,UTEMP) = T_model(k)

             end do
          end do
       end do

    endif

    ! ZHI
    if (bc(3,2,n) == EXT_DIR .and. hi(3) > domhi(3)) then

       do j= lo(2), hi(2)
          do i = lo(1), hi(1)

             kmin = domhi(3) + 1
             kmax = adv_hi(3)
#ifdef AMREX_USE_CUDA
             if (lo(3) /= kmin) then
                kmin = kmax + 1
             end if
#endif
             do k = kmin, kmax

                ! zero transverse momentum
                adv(i,j,k,UMX) = ZERO
                adv(i,j,k,UMY) = ZERO

                if (boundary_type .eq. 1) then
                   ! extrapolate normal momentum
                   ! enforces pi=0 at boundary
                   adv(i,j,k,UMZ) = adv(i,j,domhi(3),UMZ)
                else
                   ! zero normal momentum
                   ! permits pi to pass through boundary
                   adv(i,j,k,UMZ) = ZERO
                end if

                adv(i,j,k,URHO) = rho_model(k)
                adv(i,j,k,UFS:UFS-1+nspec) = adv(i,j,k,URHO)*xn_model(:)
                adv(i,j,k,UEINT) = e_model(k)*adv(i,j,k,URHO)
                adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                     + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                adv(i,j,k,UTEMP) = T_model(k)

             end do
          end do
       end do

    end if
#endif

  end subroutine ext_fill


  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only: problo, center
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use probdata_module
    use model_module

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

    real(rt), allocatable :: r_model(:), rho_model(:), T_model(:), &
                             e_model(:), p_model(:)

    integer :: lo_model, hi_model

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    ! compute background state.  Here we assume that the background is
    ! in the vertical direction.  These BCs only make sense for the
    ! boundary along the background direction (dimension
    ! AMREX_SPACEDIM).

    ! first make a 1D initial model for the entire domain + 4 ghost cells on either end
    call get_model_size(rmin, rmax, delta(AMREX_SPACEDIM), lo_model, hi_model)

    allocate(  r_model(lo_model:hi_model))
    allocate(rho_model(lo_model:hi_model))
    allocate(  T_model(lo_model:hi_model))
    allocate(  e_model(lo_model:hi_model))
    allocate(  p_model(lo_model:hi_model))

    call get_model(rmin, rmax, delta(AMREX_SPACEDIM), &
                   pres_base, dens_base, do_isentropic, &
                   xn_model, &
                   r_model, rho_model, T_model, e_model, p_model, &
                   lo_model, hi_model)

#ifndef AMREX_USE_CUDA
    if (bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for -x BC")
    endif

    if (bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for +x BC, d")
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
                adv(i,j,k) = rho_model(j)
             end do
          end do
       end do

    endif


    ! YHI
    if (bc(2,2) == EXT_DIR .and. hi(2) > domhi(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             jmin = domhi(2)+1
             jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
             if (lo(2) /= jmin) then
                jmin = jmax + 1
             end if
#endif
             do j = jmin, jmax
                adv(i,j,k) = rho_model(j)
             end do

          end do
       end do
    end if
#endif

#if AMREX_SPACEDIM == 3

    !-------------------------------------------------------------------------
    ! z boundaries
    !-------------------------------------------------------------------------

    ! ZLO
    if (bc(3,1) == EXT_DIR .and. lo(3) < domlo(3)) then

       do j= lo(2), hi(2)
          do i = lo(1), hi(1)

             ! integrate downward
             kmin = adv_lo(3)
             kmax = domlo(3) - 1
#ifdef AMREX_USE_CUDA
             if (hi(3) /= kmax) then
                kmax = kmin - 1
             end if
#endif
             do k = kmax, kmin, -1
                adv(i,j,k) = rho_model(k)
             end do
          end do
       end do

    endif

    ! ZHI
    if (bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then

       do j= lo(2), hi(2)
          do i = lo(1), hi(1)

             kmin = domhi(3) + 1
             kmax = adv_hi(3)
#ifdef AMREX_USE_CUDA
             if (lo(3) /= kmin) then
                kmin = kmax + 1
             end if
#endif
             do k = kmin, kmax
                adv(i,j,k) = rho_model(k)
             end do

          end do
       end do
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
