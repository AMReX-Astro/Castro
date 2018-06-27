module sdc_type_module

  use amrex_fort_module, only : rt => amrex_real
  use actual_network, only: nspec

  implicit none

  ! A generic structure holding data necessary to do a nuclear burn
  ! in the SDC formalism.

  ! these indicies represent the order that the conserved state comes
  ! into the ODE integration from the hydro code.
  !
  ! they also represent the order of the advective sources
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

     real(rt) :: y(SVAR)
     real(rt) :: ydot_a(SVAR)

     logical :: T_from_eden

     integer :: i
     integer :: j
     integer :: k

     integer :: n_rhs
     integer :: n_jac

     integer :: sdc_iter
  end type sdc_t

end module sdc_type_module
