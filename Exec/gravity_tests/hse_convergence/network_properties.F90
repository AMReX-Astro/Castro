! An automatically generated file of network properties.  This provides the properties
! of a set of non-reacting species.
!
! nspec            -- the number of species
! naux             -- the number of auxiliary variables
!
! aion             -- atomic number
! zion             -- proton number
! aion_inv         -- 1/atomic number
!
! spec_names       -- the name of the isotope
! short_spec_names -- an abbreviated form of the species name
!
! aux_names        -- the name of the auxiliary variable
! short_aux_names  -- an abbreviated form of the auxiliary variable


module network_properties

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: nspec = 4
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), aion_inv(:), zion(:), nion(:)

  !$acc declare create(aion, aion_inv, zion, nion)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, aion_inv, zion, nion
#endif

contains

  subroutine network_properties_init

    spec_names(1) = "helium-4"
    spec_names(2) = "carbon-12"
    spec_names(3) = "oxygen-16"
    spec_names(4) = "iron-56"

    short_spec_names(1) = "He4"
    short_spec_names(2) = "C12"
    short_spec_names(3) = "O16"
    short_spec_names(4) = "Fe56"

    allocate(aion(nspec))
    allocate(aion_inv(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))

    aion(1) = 4.0_rt
    aion(2) = 12.0_rt
    aion(3) = 16.0_rt
    aion(4) = 56.0_rt

    zion(1) = 2.0_rt
    zion(2) = 6.0_rt
    zion(3) = 8.0_rt
    zion(4) = 26.0_rt

    aion_inv(1) = 1.0_rt/4.0_rt
    aion_inv(2) = 1.0_rt/12.0_rt
    aion_inv(3) = 1.0_rt/16.0_rt
    aion_inv(4) = 1.0_rt/56.0_rt

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)




    !$acc update device(aion, zion)

  end subroutine network_properties_init



  subroutine network_properties_finalize

    implicit none

    if (allocated(aion)) then
       deallocate(aion)
    end if

    if (allocated(aion_inv)) then
       deallocate(aion_inv)
    end if

    if (allocated(zion)) then
       deallocate(zion)
    end if

    if (allocated(nion)) then
       deallocate(nion)
    end if

  end subroutine network_properties_finalize

end module network_properties
