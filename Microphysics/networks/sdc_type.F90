module sdc_type_module

  use actual_network, only: nspec
  use bl_types, only: dp_t

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn
  ! in the SDC formalism.

  integer, parameter :: SRHO  = 1
  integer, parameter :: SMX   = 2
  integer, parameter :: SMY   = 3
  integer, parameter :: SMZ   = 4
  integer, parameter :: SEDEN = 5
  integer, parameter :: SEINT = 6
  integer, parameter :: SFS   = 7

  integer, parameter :: SVAR  = SFS + nspec - 1

  type :: sdc_t

     real(dp_t) :: y(SVAR)
     real(dp_t) :: ydot_a(SVAR)

     logical :: T_from_eden

     integer :: i
     integer :: j
     integer :: k

  end type sdc_t

end module sdc_type_module
