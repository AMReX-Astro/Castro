subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C)

  use probdata_module
  use model_parser_module
  use castro_error_module
  use prob_params_module, only : center
  use amrex_constants_module, only : ZERO, HALF
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer untin, i

  namelist /fortin/ &
       model_name

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer, parameter :: maxlen = 127
  character probin*(maxlen)
  character model*(maxlen)

  if (namlen > maxlen) call castro_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults if applicable here

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  ! read initial model
  call read_model_file(model_name)

#if AMREX_SPACEDIM == 1
  ! assume spherical
  center(1) = ZERO
  center(2) = ZERO
  center(3) = ZERO
#elif AMREX_SPACEDIM == 2
  ! assume axisymmetric
  center(1) = ZERO
  center(2) = HALF*(problo(2) + probhi(2))
  center(3) = ZERO
#else
  center(1) = HALF*(problo(1) + probhi(1))
  center(2) = HALF*(problo(2) + probhi(2))
  center(3) = HALF*(problo(3) + probhi(3))
#endif

end subroutine amrex_probinit

! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.
! :::
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: level     => amr level of grid
! ::: time      => time at which to init data
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, state_lo, state_hi, &
                       delta, xlo, xhi)

  use probdata_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP,&
                                 UEDEN, UEINT, UFS
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center, problo, probhi
  use eos_type_module
  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: x, y, z, pres, r1, t0, zc, velr, rad_cyl, rad_sph
  integer :: i, j, k, n

  ! direction angles for radial -> cartesian projection
  real(rt) :: sin_theta, cos_theta, sin_phi, cos_phi

  type(eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF) - center(3)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF) - center(2)

        do i = lo(1), hi(1)
           x = problo(1) + delta(1)*(dble(i) + HALF) - center(1)

           rad_sph = sqrt(x**2 + y**2 + z**2)
           rad_cyl = sqrt(x**2 + y**2)

           call interpolate_sub(state(i,j,k,URHO), rad_sph, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), rad_sph, itemp_model)

           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), rad_sph, ispec_model-1+n)
           end do

           ! get interpolated radial velocity
           call interpolate_sub(velr, rad_sph, ivelr_model)

           ! get angle sin, cos for spherical coordinates
           sin_phi = rad_cyl/rad_sph
           cos_phi = z/rad_sph
           sin_theta = y/rad_cyl
           cos_theta = x/rad_cyl

           ! project radial velocity onto cartesian velocity components
           state(i,j,k,UMX) = state(i,j,k,URHO) * velr * sin_phi * cos_theta
           state(i,j,k,UMY) = state(i,j,k,URHO) * velr * sin_phi * sin_theta
           state(i,j,k,UMZ) = state(i,j,k,URHO) * velr * cos_phi
        end do
     end do
  end do

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do
        end do
     end do
  end do

end subroutine ca_initdata
