module model_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

  subroutine get_model_size(ymin, ymax, dy, lo, hi)
    
    use amrex_fort_module, only : rt => amrex_real
    real(rt)        , intent(in) :: ymin, ymax, dy
    integer, intent(out) :: lo, hi
    
    integer :: npts

    ! number of points in the domain
    npts = (ymax - ymin)/dy + 1
    
    ! we'll do some ghost cells, for the boundary conditions
    ! by design, the base of the model will be at zone 0
    lo = -4
    hi = npts + 4

  end subroutine get_model_size

  subroutine get_model(ymin, ymax, dy, &
                       pres_base, dens_base, do_isentropic, &
                       xn_model, &
                       r_model, rho_model, T_model, e_model, p_model, &
                       lo, hi)

    use eos_module
    use eos_type_module
    use meth_params_module, only: const_grav

    use amrex_fort_module, only : rt => amrex_real
    integer, intent(in) :: lo, hi
    real(rt)        , intent(in) :: ymin, ymax, dy
    real(rt)        , intent(in) :: pres_base, dens_base
    logical,          intent(in) :: do_isentropic
    real(rt)        , intent(in) :: xn_model(nspec)

    real(rt)        , intent(out) ::   r_model(lo:hi)
    real(rt)        , intent(out) :: rho_model(lo:hi)
    real(rt)        , intent(out) ::   T_model(lo:hi)
    real(rt)        , intent(out) ::   e_model(lo:hi)
    real(rt)        , intent(out) ::   p_model(lo:hi)

    real(rt)         :: H, gamma_const
    
    integer :: j

    type (eos_t) :: eos_state

    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(const_grav)

    ! create the constant if we are isentropic
    eos_state % rho = dens_base
    eos_state % p = pres_base
    eos_state % xn(:) = xn_model(:)

    ! initial guess
    eos_state % T = 1000.0e0_rt

    call eos(eos_input_rp, eos_state)

    gamma_const = pres_base/(dens_base * eos_state % e) + 1.0e0_rt


    rho_model(0) = dens_base
    p_model(0) = pres_base

    r_model(0) = ymin + 0.5e0_rt*dy

    ! integrate up from the base
    do j = 1, hi

       r_model(j) = ymin + (dble(j)+0.5e0_rt)*dy

       if (do_isentropic) then
          rho_model(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)* &
               (r_model(j)-r_model(0))/ &
               (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
       else
          rho_model(j) = dens_base * exp(-(r_model(j)-r_model(0))/H)
       endif

       p_model(j) = p_model(j-1) - &
             dy * 0.5e0_rt * (rho_model(j)+rho_model(j-1)) * abs(const_grav)
       
    enddo

    ! integrate down from the base
    do j = -1, lo, -1
       
       r_model(j) = ymin + (dble(j)+0.5e0_rt)*dy

       if (do_isentropic) then
          rho_model(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)* &
               (r_model(j)-r_model(0))/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
       else
          rho_model(j) = dens_base * exp(-(r_model(j)-r_model(0))/H)
       endif

       p_model(j) = p_model(j+1) + &
             dy * 0.5e0_rt * (rho_model(j)+rho_model(j+1)) * abs(const_grav)
       
    enddo
     
    ! thermodynamics
    do j = lo, hi
       eos_state % rho = rho_model(j)
       eos_state % p = p_model(j)
       eos_state % xn(:) = xn_model(:)

       ! initial guess
       eos_state % T = 1000.0e0_rt

       call eos(eos_input_rp, eos_state)

       e_model(j) = eos_state % e
       T_model(j) = eos_state % T
    end do

  end subroutine get_model

end module model_module
