!======= Sound speed calc for ideal MHD ================
  subroutine eos_soundspeed_mhd(c, rho, P, gam1, bx, by, bz, bd)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module
   

     implicit none
     ! In/out variables
     real(rt), intent(in   ) :: rho, P, gam1, bx , by, bz, bd !gas P, gam1, magnetic fields, directional mag field
    
     real(rt), intent(  out) :: c


     !Sound Speed, Alfven Speed
     real(rt) :: as, ca, cad

     as = gam1 * P/rho
     ca = (bx**2 + by**2 + bz**2)/rho
     cad = bd**2/rho
     !Fast Magneto-Sonic Wave
     c = 0.5d0*((as + ca) + sqrt((as + ca)**2 -4.0d0*as*cad))
     c = sqrt(c)

  end subroutine eos_soundspeed_mhd

