module castro_burner_module

  use bl_types
  use bl_constants_module
  use network

contains

  subroutine burner(dens, temp, Xin, ein, dt, time, Xout, eout)

    ! outputs:
    !   Xout are the mass fractions after burning through timestep dt
    !   eout is the updated specific internal energy

    implicit none
    
    real(kind=dp_t), intent(in) :: dens, temp, Xin(nspec), ein, dt, time
    real(kind=dp_t), intent(out) :: Xout(nspec), eout
        
    Xout(:) = Xin(:)
    eout = ein
  
  end subroutine burner

end module castro_burner_module
