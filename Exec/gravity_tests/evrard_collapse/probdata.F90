module probdata_module

  ! Probin file
  use amrex_fort_module, only : rt => amrex_real
  character (len=:), allocatable :: probin

  ! Determine if we are the I/O processor
  integer :: ioproc

  ! Mass and radius of the collapsing sphere
  real(rt) :: sphere_mass, sphere_radius

  ! Smallest allowed mass fraction
  real(rt) :: smallx

  ! Smallest allowed velocity on the grid
  real(rt) :: smallu

  ! Density of ambient gas around the star
  real(rt) :: ambient_density

contains

  ! This routine calls all of the other subroutines at the beginning
  ! of a simulation to fill in the basic problem data.

  subroutine initialize(name, namlen)

    use amrex_constants_module, only: ZERO
    use castro_error_module, only: castro_error

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: namlen, i
    integer :: name(namlen)

    ! Build "probin" filename -- the name of the file containing the fortin namelist.
    allocate(character(len=namlen) :: probin)
    do i = 1, namlen
       probin(i:i) = char(name(i))
    enddo

    ! Determine if we are the I/O processor, and save it to the ioproc variable

    call get_ioproc

    ! Read in the namelist to set problem parameters

    call read_namelist

    ! Finalize ambient state, and get small_pres

    call set_small

  end subroutine



  ! This routine reads in the namelist

  subroutine read_namelist

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer :: untin

    namelist /fortin/ sphere_mass, sphere_radius, ambient_density, smallx, smallu

    ! Set namelist defaults

    sphere_mass   = 1.0e0_rt
    sphere_radius = 9.0e8_rt

    ambient_density = 1.0e0_rt

    smallx = 1.e-10_rt
    smallu = 1.e-12_rt

    ! Read namelist to override the defaults

    open(newunit=untin, file=probin, form='formatted', status='old')
    read(untin,fortin)
    close(unit=untin)

  end subroutine read_namelist



  ! Determine if we are the I/O processor, and save it to the ioproc variable

  subroutine get_ioproc

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    ! For outputting -- determine if we are the IO processor
    call bl_pd_is_ioproc(ioproc)

  end subroutine get_ioproc



  ! Define the ambient state, and calculate small_pres

  subroutine set_small

    use network, only: nspec
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos
    use meth_params_module, only: small_temp, small_pres, small_dens, small_ener
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    type (eos_t) :: eos_state

    ! Given the inputs of small_dens and small_temp, figure out small_pres.

    if (small_dens > 0.0e0_rt .and. small_temp > 0.0e0_rt) then
       eos_state % rho = small_dens
       eos_state % T   = small_temp
       eos_state % xn  = 1.0e0_rt / nspec

       call eos(eos_input_rt, eos_state)

       small_pres = eos_state % p
       small_ener = eos_state % e
    endif

  end subroutine set_small

end module probdata_module
