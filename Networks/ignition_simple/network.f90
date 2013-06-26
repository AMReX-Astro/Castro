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
! short_spec_names -- the abbreviated name of the isotope
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

  ! nspec = number of species this network carries
  ! nspec_advance = the number of species that are explicitly integrated
  !                 in the ODE solve (the others are solved for 
  !                 algebraically).
  integer, parameter :: nspec = 3
  integer, parameter :: nspec_advance = 1
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)

  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    integer :: ic12, io16, img24

    ! integer keys -- for convinence.  In all other places, we will find
    ! these by querying based on species name using network_species_index
    ic12  = 1
    io16  = 2
    img24 = 3

    spec_names(ic12)  = "carbon-12"
    spec_names(io16)  = "oxygen-16"
    spec_names(img24) = "magnesium-24"

    short_spec_names(ic12)  = "C12"
    short_spec_names(io16)  = "O16"
    short_spec_names(img24) = "Mg24"

    aion(ic12)  = 12.0_dp_t
    aion(io16)  = 16.0_dp_t
    aion(img24) = 24.0_dp_t
    
    zion(ic12)  = 6.0_dp_t
    zion(io16)  = 8.0_dp_t
    zion(img24) = 12.0_dp_t

    ! our convention is that the binding energies are negative.  We convert
    ! from the MeV values that are traditionally written in astrophysics 
    ! papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    ! MeV values are per nucleus, so we divide by aion to make it per
    ! nucleon and we multiple by Avogardo's # (6.0221415e23) to get the 
    ! value in erg/g
    ebin(ic12)  = -7.4103097e18_dp_t     !  92.16294 MeV
    ebin(io16)  = -7.6959672e18_dp_t     ! 127.62093 MeV
    ebin(img24) = -7.9704080e18_dp_t     ! 198.2579  MeV

    network_initialized = .true.

  end subroutine network_init

  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          r = n
          exit
       endif
    enddo
    return
  end function network_species_index

end module network
