module eos_data_module

  use bl_types

  implicit none

  integer, parameter :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter :: eos_input_ps = 6  ! p, s are inputs
  integer, parameter :: eos_input_ph = 7  ! p, h are inputs
  integer, parameter :: eos_input_th = 8  ! T, h are inputs

  integer, parameter :: itemp = 1
  integer, parameter :: idens = 2
  integer, parameter :: iener = 3
  integer, parameter :: ienth = 4
  integer, parameter :: ientr = 5
  integer, parameter :: ipres = 6

  integer, parameter :: ierr_general         = 1
  integer, parameter :: ierr_input           = 2
  integer, parameter :: ierr_iter_conv       = 3
  integer, parameter :: ierr_neg_e           = 4
  integer, parameter :: ierr_neg_p           = 5
  integer, parameter :: ierr_neg_h           = 6
  integer, parameter :: ierr_neg_s           = 7
  integer, parameter :: ierr_iter_var        = 8
  integer, parameter :: ierr_init            = 9
  integer, parameter :: ierr_init_xn         = 10
  integer, parameter :: ierr_out_of_bounds   = 11
  integer, parameter :: ierr_not_implemented = 12

  real(kind=dp_t), save :: smallt
  real(kind=dp_t), save :: smalld

  logical, save :: initialized = .false.

  public eos_get_small_temp, eos_get_small_dens

contains

  subroutine eos_get_small_temp(small_temp_out)
 
    real(kind=dp_t), intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
 
    real(kind=dp_t), intent(out) :: small_dens_out
 
    small_dens_out = smalld
 
  end subroutine eos_get_small_dens

end module eos_data_module
