! the network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
! eion       -- nuclear binding energy (in erg/g)
!
! spec_names -- the name of the isotope
!
!
! This module contains two routines:
!
!  network_init()        -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!

module network

  use bl_types

  implicit none

  integer        , parameter :: nspec  = 1
  integer        , parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)

  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    spec_names(1) = "X"
    short_spec_names(1) = "X"

    aion(1) = 1.0_dp_t
    zion(1) = 1.0_dp_t
    ebin(1) = 0.0_dp_t

    network_initialized = .true.

  end subroutine network_init

  function network_species_index(name)

    character(len=*) :: name
    integer :: network_species_index, n

    network_species_index = -1

    do n = 1, nspec
       if (trim(name) == trim(spec_names(n))) then
          network_species_index = n
          exit
       endif
    enddo
    
    return
  end function network_species_index

end module network
