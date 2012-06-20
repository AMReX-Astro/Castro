module castro_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

contains

  subroutine burner(dens, temp, Xin, ein, dt, time, Xout, eout)

    implicit none
    
    real(kind=dp_t), intent(in) :: dens, temp, Xin(nspec), ein, dt, time
    real(kind=dp_t), intent(out) :: Xout(nspec), eout
    
    integer :: n
    real(kind=dp_t) :: enuc, dX
    
    Xout(:) = Xin(:)
    
    enuc = 0.0_dp_t
    do n = 1, nspec
       dX = Xout(n)-Xin(n) 
       enuc = enuc - ebin(n) * dX
    enddo
  
    eout = ein + enuc
  
  end subroutine burner

end module castro_burner_module
