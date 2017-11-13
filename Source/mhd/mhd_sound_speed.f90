!======= Sound speed calc for ideal MHD ================
  subroutine eos_soundspeed_mhd(c, R, e, bx, by, bz, bd)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only: gamma_const, gamma_minus_1

     ! In/out variables
     real(rt), intent(in   ) :: R, e, bx , by, bz, bd !density, internal energy, magnetic fields, directional mag field
     real(rt), intent(  out) :: c

     !Gas Pressure
     real(rt) :: P
     !Sound Speed, Alfven Speed
     real(rt) :: as, ca

     P = R * e * gamma_minus_1
     as = gamma_const * P/R
     ca = (bx**2 + by**2 + bz**2)/R
     !Fast Magneto-Sonic Wave
     c = gamma_const*P + (bx**2 + by**2 + bx**2) + sqrt((gamma_const*P + (bx**2+by**2+bz**2))**2 - 4.d0*gamma_const*P*(bd))
     c = 0.5d0*c/R
     c = sqrt(c)

  end subroutine eos_soundspeed_mhd

