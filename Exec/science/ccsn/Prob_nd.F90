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
        model_name, min_density, min_temperature, fluff_ye, &
        tag_max_density_fraction

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

  ! Initial model mapping parameter defaults

  !! default the minimum model density to just above the Weaklib EOS table lower limit in density
  min_density = 2.0e3_rt

  !! default the fluff electron fraction to 0.5
  fluff_ye = 0.5e0_rt

  !! default the minimum model temperature to just above the Weaklib EOS table lower limit in temperature
  min_temperature = 2.0e9_rt

  ! Tagging parameter defaults

  !! default the tag_max_density_fraction to 0.1
  tag_max_density_fraction = 0.1_rt

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
                                 UEDEN, UEINT, UFS, UFX, small_dens
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center, problo, probhi
  use eos_type_module
  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use amrex_fort_module, only : rt => amrex_real
  use UnitsModule, only: Gram, AtomicMassUnit

! Geometry Notes:
! 1D is not yet supported
! For 2D, we use Cylindrical r-z coordinates
! For 3D, we use Cartesian coordinates

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: x, y, z, pres, r1, t0, zc, velr, rad_sph
  integer :: i, j, k, n

  ! direction angles for radial -> cartesian projection and cylindrical radius
  real(rt) :: sin_theta, cos_theta, sin_phi, cos_phi, rad_cyl

  type(eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF) - center(3)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF) - center(2)

        do i = lo(1), hi(1)
           x = problo(1) + delta(1)*(dble(i) + HALF) - center(1)

           rad_sph = sqrt(x**2 + y**2 + z**2)

           call interpolate_sub(state(i,j,k,UTEMP), rad_sph, itemp_model)

           ! enforce the minimum temperature in the problem parameters
           state(i,j,k,UTEMP) = max(state(i,j,k,UTEMP), min_temperature)

           call interpolate_sub(state(i,j,k,URHO), rad_sph, idens_model)

           ! enforce the minimum density in the problem parameters
           ! and set Ye below this density to fluff_ye in the problem parameters
           if (state(i,j,k,URHO) .le. min_density) then
               state(i,j,k,URHO) = min_density
               do n = 1, nspec
                  state(i,j,k,UFS-1+n) = fluff_ye
               end do
           else
               do n = 1, nspec
                  call interpolate_sub(state(i,j,k,UFS-1+n), rad_sph, ispec_model-1+n)
               end do
           end if

           ! get interpolated radial velocity
           call interpolate_sub(velr, rad_sph, ivelr_model)

#if (AMREX_SPACEDIM==2)
           ! project spherical radial velocity onto cylindrical r-z velocity components
           rad_cyl = x
           sin_theta = rad_cyl/rad_sph
           cos_theta = y/rad_sph

           ! Cylindrical radial velocity
           state(i,j,k,UMX) = state(i,j,k,URHO) * velr * sin_theta

           ! Cylindrical axial velocity
           state(i,j,k,UMY) = state(i,j,k,URHO) * velr * cos_theta

           ! Out-of-plane velocity
           state(i,j,k,UMZ) = ZERO
#elif (AMREX_SPACEDIM==3)
           ! get angle sin, cos for spherical coordinates
           rad_cyl = sqrt(x**2 + y**2)
           sin_theta = rad_cyl/rad_sph
           cos_theta = z/rad_sph
           sin_phi = y/rad_cyl
           cos_phi = x/rad_cyl

           ! project radial velocity onto cartesian velocity components
           state(i,j,k,UMX) = state(i,j,k,URHO) * velr * sin_theta * cos_phi
           state(i,j,k,UMY) = state(i,j,k,URHO) * velr * sin_theta * sin_phi
           state(i,j,k,UMZ) = state(i,j,k,URHO) * velr * cos_theta
#endif
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
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                HALF * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2)/state(i,j,k,URHO)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

           ! Store Ne as an auxiliary quantity
           state(i,j,k,UFX) = state(i,j,k,UFS) * Gram / AtomicMassUnit
        end do
     end do
  end do

end subroutine ca_initdata

! hardwired assuming 4 moments
subroutine get_rad_ncomp(rad_ncomp) bind(C,name="ca_get_rad_ncomp")

  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule, only : nE, nNodesE

  integer :: rad_ncomp
  integer :: n_moments = 4

  rad_ncomp =  nSpecies * n_moments * nE * nNodesE

end subroutine get_rad_ncomp

! hardwired assuming 4 moments
subroutine get_thornado_node_averages(n_node_avgs) bind(C,name="ca_get_thornado_node_averages")

  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule  , only : nE

  integer :: n_node_avgs
  integer :: n_moments = 4

  n_node_avgs =  nSpecies * n_moments * nE

end subroutine get_thornado_node_averages

! hardwired assuming 4 moments
! streaming sine wave, J = H_x = 1 + sin(2*pi*x)
subroutine ca_init_thornado_data(level, time, lo, hi, nrad_comp, &
                                 rad_state, rad_state_lo, rad_state_hi, &
                                 state, state_lo, state_hi, &
                                 delta, xlo, xhi) bind(C,name="ca_init_thornado_data")

  use probdata_module
  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule, only : nE, nNodesE
  use amrex_fort_module, only : rt => amrex_real
  use amrex_error_module
  use amrex_constants_module, only : M_PI, zero, one
  use meth_params_module, only: NVAR, URHO, UEINT, UTEMP, UFX
  use MeshModule, only: MeshE, NodeCoordinate
  use UnitsModule
  use wlEOSInversionModule, only: DescribeEOSInversionError
  use NeutrinoOpacitiesComputationModule, only: FermiDirac
  use EquationOfStateModule_TABLE, only: ComputeTemperatureFromSpecificInternalEnergy_TABLE, &
      ComputeElectronChemicalPotential_TABLE, &
      ComputeProtonChemicalPotential_TABLE, &
      ComputeNeutronChemicalPotential_TABLE

  implicit none

  integer,  intent(in   ) :: level, nrad_comp
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: rad_state_lo(3), rad_state_hi(3)
  real(rt), intent(inout) :: rad_state(rad_state_lo(1):rad_state_hi(1),&
                                       rad_state_lo(2):rad_state_hi(2),&
                                       rad_state_lo(3):rad_state_hi(3),&
                                       0:nrad_comp-1)
  integer,  intent(in   ) :: state_lo(3), state_hi(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),&
                                   state_lo(2):state_hi(2),&
                                   state_lo(3):state_hi(3),NVAR)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)

  ! Local parameter
  integer, parameter :: n_moments = 4
  real(rt), parameter :: Log1d100 = LOG( 1.0d100 )

  ! local variables
  integer :: i,j,k,ienode,Error
  integer :: ii,ii_0,is,im,ie
  real(rt) :: rho_in, T_in, Ye_in, Evol, Ne_loc, Em_in, M_e, M_p, M_n, M_nu, E

  ! zero it out, just in case
  rad_state = zero

  ! print *,'nrad_comp ',nrad_comp
  ! print *,'nSpecies  ',nSpecies
  ! print *,'n_moments ',n_moments
  ! print *,'nE        ',nE
  ! print *,'nNodesE   ',nNodesE
  ! print *,'MULT ', nSpecies * n_moments * nE * nNodesE

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ! Get Castro fluid variables unit convert to thornado units
           rho_in = state(i,j,k,URHO ) * Gram / Centimeter**3
           T_in   = state(i,j,k,UTEMP) * Kelvin
           Evol   = state(i,j,k,UEINT) * (Erg/Centimeter**3)
           Ne_loc = state(i,j,k,UFX  ) / Centimeter**3

           ! Calculate chemical potentials via thornado subroutines
           Em_in  = Evol / rho_in
           Ye_in  = Ne_loc / rho_in * AtomicMassUnit
           call ComputeTemperatureFromSpecificInternalEnergy_TABLE(rho_in,Em_in,Ye_in,T_in,Error_Option=Error)
           if ( Error > 0 ) then
             call DescribeEOSInversionError( Error )
             stop
           end if
           call ComputeElectronChemicalPotential_TABLE(rho_in,T_in,Ye_in,M_e)
           call ComputeProtonChemicalPotential_TABLE(rho_in,T_in,Ye_in,M_p)
           call ComputeNeutronChemicalPotential_TABLE(rho_in,T_in,Ye_in,M_n)

           M_nu = M_e + M_p - M_n

           do is = 1, nSpecies
           do im = 1, n_moments
           do ie = 1, nE

              ii_0 = (is-1)*(n_moments*nE*nNodesE) + (im-1)*(nE*nNodesE) + (ie-1)*nNodesE

              do ienode = 1, nNodesE

                 ! Calculate the indices
                 ii = ii_0 + (ienode-1)

                 ! Get energy at given node coordinate via thornado subroutine
                 E = NodeCoordinate( MeshE, ie, ienode)

                 ! J moment, im = 1, is = 1
                 if (im .eq. 1 .and. is .eq. 1) then
                    rad_state(i,j,k,ii) = FermiDirac( E, +M_nu, BoltzmannConstant * T_in )
                 end if

                 ! J moment, im = 1, is = 2
                 if (im .eq. 1 .and. is .eq. 2) then
                    rad_state(i,j,k,ii) = FermiDirac( E, -M_nu, BoltzmannConstant * T_in )
                 end if

                 ! H_x moment, im = 2
                 if (im .eq. 2) rad_state(i,j,k,ii) = zero

                 ! H_y moment, im = 3
                 if (im .eq. 3) rad_state(i,j,k,ii) = zero

                 ! H_z moment, im = 4
                 if (im .eq. 4) rad_state(i,j,k,ii) = zero

              end do

           end do
           end do
           end do

        enddo
     enddo
  enddo

end subroutine ca_init_thornado_data
