module burner_module

  use bl_types
  use bl_constants_module
  use network

contains

  subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

    ! outputs:
    !   Xout are the mass fractions after burning through timestep dt
    !   rho_omegadot = rho dX/dt
    !   rho_Hnuc = energy generation in erg/s/cm^3

    implicit none
    
    real(kind=dp_t), intent(in) :: dens, temp, Xin(nspec), dt
    real(kind=dp_t), intent(out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
        
    Xout(:) = Xin(:)
    
    rho_omegadot(:) = ZERO
    rho_Hnuc = ZERO
  
  end subroutine burner

end module burner_module
