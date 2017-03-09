module probdata_module

 ! Problem setup data
  use bl_fort_module, only : rt => c_real
  real(rt)         :: rho1, rho2, pressure

  ! Problem number
  integer :: problem

  ! Uniform flow speed
  real(rt)         :: bulk_velocity

end module probdata_module
