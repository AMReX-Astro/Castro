!======= Sound speed calc for ideal MHD ================
  subroutine eos_soundspeed_mhd(c, rho, e, temp, bx, by, bz, bd, compo)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module
     use eos_module, only : eos
     use eos_type_module, only : eos_t, eos_input_re
     use network, only: nspec

     implicit none
     ! In/out variables
     real(rt), intent(in   ) :: rho, e, temp, bx , by, bz, bd !density, internal energy, magnetic fields, directional mag field
     real(rt), intent(in   ) :: compo(nspec)
     real(rt), intent(  out) :: c

     !Gas Pressure
     real(rt) :: P
     !Sound Speed, Alfven Speed
     real(rt) :: as, ca

     type(eos_t) :: eos_state

     eos_state % rho = rho
     eos_state % e   = e
     eos_state % T   = temp
     eos_state % xn  = compo
     
     call eos(eos_input_re, eos_state)
  
     P = eos_state % p
     as = eos_state % cs
     ca = (bx**2 + by**2 + bz**2)/rho
     !Fast Magneto-Sonic Wave
     ! TODO: double check this expression
     c = eos_state % gam1 * P + (bx**2 + by**2 + bx**2) + &
         sqrt((eos_state % gam1 * P + (bx**2+by**2+bz**2))**2 - 4.d0 * eos_state % gam1 *P*(bd))
     c = 0.5d0*c/rho
     c = sqrt(c)

  end subroutine eos_soundspeed_mhd

