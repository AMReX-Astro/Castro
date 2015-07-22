! the network module provides the information about the species we are
! advecting: 
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
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
  use specific_network
  
  implicit none

  double precision :: smallx = 1.d-199
  
  logical :: network_initialized = .false.

contains
  
  subroutine network_init

    implicit none
    
    ! First, we call the specific network initialization.
    ! This should set the number of species and number of
    ! aux variables, and the components of the species.
    ! Note that the network MUST define nspec and naux
    ! as parameters, or else the compiler will throw an error.
    
    call specific_network_init

    ! Check to make sure, and if not, throw an error.

    if ( nspec .le. 0 ) then
       call bl_error("Network cannot have a negative number of species.")
    endif

    if ( naux .lt. 0 ) then
       call bl_error("Network cannot have a negative number of auxiliary variables.")
    endif

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
