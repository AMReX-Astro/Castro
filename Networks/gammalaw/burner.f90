module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

contains

  subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

    ! outputs:
    !   Xout are the mass fractions after burning through timestep dt
    !   rho_omegadot = rho dX/dt
    !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

    implicit none
    
    real(kind=dp_t), intent(in) :: dens, temp, Xin(nspec), dt
    real(kind=dp_t), intent(out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
    
    integer :: n
    real(kind=dp_t) :: enuc, dX
    
    Xout(:) = Xin(:)
    
    ! compute the energy release.  Our convention is that the binding
    ! energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    enuc = 0.0_dp_t
    do n = 1, nspec
       dX = Xout(n)-Xin(n) 
       
       enuc = enuc - ebin(n) * dX
       
       rho_omegadot(n) = dens * dX / dt
    enddo
  
    rho_Hnuc = dens*enuc/dt
  
  end subroutine burner

end module burner_module
