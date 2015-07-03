module eos_data_module

  implicit none

  integer :: eos_input_rt = 1  ! rho, T are inputs
  integer :: eos_input_rh = 2  ! rho, h are inputs
  integer :: eos_input_tp = 3  ! T, p are inputs
  integer :: eos_input_rp = 4  ! rho, p are inputs
  integer :: eos_input_re = 5  ! rho, e are inputs
  integer :: eos_input_ps = 6  ! p, s are inputs
  integer :: eos_input_ph = 7  ! p, h are inputs
  integer :: eos_input_th = 8  ! T, h are inputs

  integer :: itemp = 1
  integer :: idens = 2
  integer :: iener = 3
  integer :: ienth = 4
  integer :: ientr = 5
  integer :: ipres = 6

  integer :: ierr_general         = 1
  integer :: ierr_input           = 2
  integer :: ierr_iter_conv       = 3
  integer :: ierr_neg_e           = 4
  integer :: ierr_neg_p           = 5
  integer :: ierr_neg_h           = 6
  integer :: ierr_neg_s           = 7
  integer :: ierr_iter_var        = 8
  integer :: ierr_init            = 9
  integer :: ierr_init_xn         = 10
  integer :: ierr_out_of_bounds   = 11
  integer :: ierr_not_implemented = 12

  double precision, save :: smallt
  double precision, save :: smalld

  logical, save :: initialized = .false.

  public eos_get_small_temp, eos_get_small_dens

contains

  subroutine eos_get_small_temp(small_temp_out)
 
    double precision,  intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
 
    double precision, intent(out) :: small_dens_out
 
    small_dens_out = smalld
 
  end subroutine eos_get_small_dens

end module eos_data_module
