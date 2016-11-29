! the network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
!
! aion_inv   -- 1/aion
!
! spec_names -- the name of the isotope
! short_spec_names -- the abbreviated name of the isotope
!
! This module contains the following routines:
!
!  network_init          -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name

module network

  use bl_types
  use bl_constants_module, only : ONE
  use actual_network
  use bl_error_module

  implicit none

  logical :: network_initialized = .false.

  ! this will be computed here, not in the actual network
  real(kind=dp_t) :: aion_inv(nspec)

  !$acc declare create(aion_inv)

contains

  subroutine network_init

    implicit none

    ! First, we call the specific network initialization.
    ! This should set the number of species and number of
    ! aux variables, and the components of the species.
    ! Note that the network MUST define nspec and naux
    ! as parameters, or else the compiler will throw an error.

    call actual_network_init()

    ! Check to make sure, and if not, throw an error.

    if ( nspec .le. 0 ) then
       call bl_error("Network cannot have a negative number of species.")
    endif

    if ( naux .lt. 0 ) then
       call bl_error("Network cannot have a negative number of auxiliary variables.")
    endif

    aion_inv(:) = ONE/aion(:)

    !$acc update device(aion_inv)

    network_initialized = .true.

  end subroutine network_init


  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          r = n
          return
       endif
    enddo

  end function network_species_index

end module network
