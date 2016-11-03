module sdc_type_module

  use actual_network, only: nspec
  use bl_types, only: dp_t

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn
  ! in the SDC formalism.

  integer, parameter :: SEDEN = 1
  integer, parameter :: SEINT = 2
  integer, parameter :: SFS   = 3
  integer, parameter :: SRHO  = SFS + nspec
  integer, parameter :: SMX   = SRHO + 1
  integer, parameter :: SMY   = SRHO + 2
  integer, parameter :: SMZ   = SRHO + 3

  integer, parameter :: SVAR  = SMZ
  integer, parameter :: SVAR_EVOLVE = SRHO - 1

  type :: sdc_t

     real(dp_t) :: y(SVAR)
     real(dp_t) :: ydot_a(SVAR)

     logical :: T_from_eden

     integer :: i
     integer :: j
     integer :: k

     integer :: n_rhs
     integer :: n_jac

     ! this is not actually needed for SDC, but included for compatibility
     ! with the non-SDC
     logical :: self_heat
  end type sdc_t

end module sdc_type_module
