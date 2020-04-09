!======= Sound speed calc for ideal MHD ================
  subroutine eos_soundspeed_mhd(c, as, ca, bd)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module
   

     implicit none
     ! In/out variables
     real(rt), intent(in   ) :: as, ca, bd !P_g* gam1/rho, B^2/rho, B_direction^2/rho
    
     real(rt), intent(  out) :: c

     

     !Fast Magneto-Sonic Wave
     c = 0.5d0*((as + ca) + sqrt((as + ca)**2 -4.0d0*as*bd))
     c = sqrt(c)

  end subroutine eos_soundspeed_mhd

