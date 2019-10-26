! Process a sedov problem to produce rho, u, p and e as a
! function of r, for comparison to the analytic solution.

subroutine microphysics_initialize(name, namlen) bind(C, name='microphysics_initialize')
    use eos_module, only: eos_init
    use network, only: network_init
    use actual_rhs_module, only: actual_rhs_init
    use amrex_fort_module, only : rt => amrex_real

    integer, intent(in) :: namlen
    integer, intent(in) :: name(namlen)

    real(rt) :: small_dens = 1.e-20_rt, small_temp = 1.e-20_rt

    call runtime_init(name, namlen)

    call eos_init(small_dens=small_dens, small_temp=small_temp)

    call network_init()
    call actual_rhs_init()

end subroutine microphysics_initialize

subroutine microphysics_finalize() bind(C, name='microphysics_finalize')

    use eos_module, only: eos_finalize
    use network, only: network_finalize

    call eos_finalize()
    call network_finalize()

end subroutine microphysics_finalize


subroutine find_timestep_limiter(lo, hi, state, slo, shi, nc_s, &
     dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp, spec_comp, time_integration_method, dx, dt, dt_loc) &
     bind(C, name='find_timestep_limiter')

  use amrex_fort_module, only : rt => amrex_real, amrex_min
  use amrex_constants_module
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_re
  use network, only: nspec

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: slo(3), shi(3), nc_s
  real(rt), intent(in) :: state(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),0:nc_s-1)
  integer, intent(in), value :: dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp, spec_comp, time_integration_method
  real(rt), intent(in) :: dx(3)
  real(rt), intent(inout) :: dt
  real(rt), intent(inout) :: dt_loc(3)

  integer :: i, j, k
  real(rt) :: x, y, z, cs, dt1, dt2, dt3, dt_tmp, ux, uy, uz

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = (dble(k) + HALF)*dx(3)
     do j = lo(2), hi(2)
        y = (dble(j) + HALF)*dx(2)

        do i = lo(1), hi(1)
            x = (dble(i) + HALF)*dx(1)

            eos_state % rho = state(i,j,k,dens_comp)
            eos_state % e = state(i,j,k,rhoe_comp) / eos_state % rho
            eos_state % xn = state(i,j,k,spec_comp:spec_comp+nspec-1)

            call eos(eos_input_re, eos_state)

            cs = eos_state % cs

            ! calculate CFL timestep using the velocity 
            ux = state(i,j,k,xmom_comp) / state(i,j,k,dens_comp)

            dt1 = dx(1) / (cs + abs(ux))
            dt2 = dt1
            dt3 = dt1
#if AMREX_SPACEDIM >= 2
            uy = state(i,j,k,ymom_comp) / state(i,j,k,dens_comp)
            dt2 = dx(2) / (cs + abs(uy))
#endif
#if AMREX_SPACEDIM == 3
            uz = state(i,j,k,zmom_comp) / state(i,j,k,dens_comp)
            dt3 = dx(3) / (cs + abs(uz))
#endif

            if (time_integration_method == 0) then 
                if (min(dt1, dt2, dt3) < dt) then 
                    dt = min(dt1, dt2, dt3)
                    dt_loc(1:3) = (/ x, y, z /)
                endif
            else 
                dt_tmp = ONE / dt1
#if AMREX_SPACEDIM >= 2
                dt_tmp = dt_tmp + ONE / dt2
#endif
#if AMREX_SPACEDIM == 3
                dt_tmp = dt_tmp + ONE / dt3
#endif
                if (ONE/dt_tmp < dt) then 
                    dt = ONE/dt_tmp
                    dt_loc(1:3) = (/ x, y, z /)
                endif

            endif
        enddo
     enddo
  enddo

end subroutine find_timestep_limiter


subroutine find_timestep_limiter_burning(lo, hi, state, slo, shi, nc_s, &
    dens_comp, temp_comp, rhoe_comp, spec_comp, dx, dt, dt_loc) &
    bind(C, name='find_timestep_limiter_burning')

    use amrex_fort_module, only : rt => amrex_real, amrex_min
    use amrex_constants_module
    use actual_rhs_module, only: actual_rhs
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use burner_module, only: ok_to_burn
    use burn_type_module, only : burn_t, net_ienuc, burn_to_eos, eos_to_burn
    use temperature_integration_module, only: self_heat
    use network, only: nspec, aion
    use meth_params_module, only : dtnuc_e, dtnuc_X, dtnuc_X_threshold
    use extern_probin_module, only: small_x

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: slo(3), shi(3), nc_s
    real(rt), intent(in) :: state(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),0:nc_s-1)
    integer, intent(in), value :: dens_comp, temp_comp, rhoe_comp, spec_comp
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt
    real(rt), intent(inout) :: dt_loc(3)

    integer :: i, j, k, n
    real(rt) :: x, y, z, e, Xn(nspec), dedt, dXdt(nspec)

    type (eos_t) :: eos_state

    type (burn_t) :: state_old, state_new

    ! Set a floor on the minimum size of a derivative. This floor
    ! is small enough such that it will result in no timestep limiting.

    real(rt), parameter :: derivative_floor = 1.e-50_rt

    do k = lo(3), hi(3)
        z = (dble(k) + HALF)*dx(3)
        do j = lo(2), hi(2)
            y = (dble(j) + HALF)*dx(2)
            do i = lo(1), hi(1)
                x = (dble(i) + HALF)*dx(1)

                state_old % rho = state(i,j,k,dens_comp)
                state_old % T = state(i,j,k,temp_comp) 
                state_old % e = state(i,j,k,rhoe_comp) / state(i,j,k,dens_comp)
                state_old % xn = state(i,j,k,spec_comp:spec_comp+nspec-1)

                state_new % rho = state(i,j,k,dens_comp)
                state_new % T = state(i,j,k,temp_comp) 
                state_new % e = state(i,j,k,rhoe_comp) / state(i,j,k,dens_comp)
                state_new % xn = state(i,j,k,spec_comp:spec_comp+nspec-1)

                ! if (.not. ok_to_burn(state_new)) cycle

                e    = state_new % e
                Xn    = max(state_new % xn, small_x)

                call burn_to_eos(state_new, eos_state)
                call eos(eos_input_rt, eos_state)
                call eos_to_burn(eos_state, state_new)

                state_new % dx = minval(dx(1:AMREX_SPACEDIM))

                state_new % self_heat = .false.

                call actual_rhs(state_new)

                dedt = state_new % ydot(net_ienuc)
                dXdt = state_new % ydot(1:nspec) * aion

                dedt = max(abs(dedt), derivative_floor)

                do n = 1, nspec
                    if (Xn(n) .ge. dtnuc_X_threshold) then 
                        dXdt(n) = max(abs(dXdt(n)), derivative_floor)
                    else 
                        dXdt(n) = derivative_floor
                    endif
                enddo

                if (dtnuc_e * e / dedt < dt) then 
                    dt = dtnuc_e * e / dedt
                    dt_loc = (/ x, y, z /)
                else if (dtnuc_X * minval(X / dXdt) < dt) then 
                    dt = dtnuc_X * minval(X / dXdt)
                    dt_loc = (/ x, y, z /)
                endif
                
            enddo
        enddo
    enddo

end subroutine find_timestep_limiter_burning

#ifdef DIFFUSION
subroutine find_timestep_limiter_diffusion(lo, hi, state, slo, shi, nc_s, &
    dens_comp, temp_comp, rhoe_comp, spec_comp, dx, dt, dt_loc) &
    bind(C, name='find_timestep_limiter_diffusion')

    use amrex_fort_module, only : rt => amrex_real, amrex_min
    use amrex_constants_module
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use conductivity_module, only: conductivity
    use meth_params_module, only : diffuse_cutoff_density

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: slo(3), shi(3), nc_s
    real(rt), intent(in) :: state(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),0:nc_s-1)
    integer, intent(in), value :: dens_comp, temp_comp, rhoe_comp, spec_comp
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt
    real(rt), intent(inout) :: dt_loc(3)

    integer :: i, j, k
    real(rt) :: x, y, z, dt1, dt2, dt3, D

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
        z = (dble(k) + HALF)*dx(3)
        do j = lo(2), hi(2)
            y = (dble(j) + HALF)*dx(2)
            do i = lo(1), hi(1)
                x = (dble(i) + HALF)*dx(1)

                if (state(i,j,k,dens_comp) > diffuse_cutoff_density) then

                    eos_state % rho = state(i,j,k,dens_comp)
                    eos_state % T = state(i,j,k,temp_comp) 
                    eos_state % e = state(i,j,k,rhoe_comp) / state(i,j,k,dens_comp)
                    state_old % xn = state(i,j,k,spec_comp:spec_comp+nspec-1)

                    call eos(eos_input_re, eos_state)

                    call conductivity(eos_state)

                    D = eos_state % conductivity / (eos_state % rho * eos_state % cv)

                    dt1 = HALF * dx(1)**2 / D 
                    dt2 = dt1 
                    dt3 = dt1

#if AMREX_SPACEDIM >= 2
                    dt2 = HALF * dx(2)**2 / D
#endif 
#if AMREX_SPACEDIM == 3
                    dt3 = HALF * dx(3)**2 / D
#endif
                    if (min(dt1, dt2, dt3) < dt) then 
                        dt = min(dt1, dt2, dt3)
                        dt_loc = (/ x, y, z /)
                    endif
                endif
            enddo
        enddo
    enddo

end subroutine find_timestep_limiter_diffusion
#endif