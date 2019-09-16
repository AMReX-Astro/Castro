subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use model_parser_module
  use amrex_error_module
  use prob_params_module, only : center
  use probdata_module, only : heating_factor, g0, rho0, p0, gamma1
  use model_util_module, only : integrate_model
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  namelist /fortin/ &
       heating_factor, g0, rho0, p0, gamma1, do_pert, ny

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 127
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! allocate probdata variables
  allocate(heating_factor, g0, rho0, p0, gamma1)
  allocate(do_pert)
  allocate(ny)

  ! set namelist defaults

  heating_factor = 1.e3_rt
  g0 = -9.021899571e8_rt
  rho0 = 1.82094e6_rt
  p0 = 2.7647358e23_rt
  gamma1 = 1.4e0_rt
  do_pert = .true.
  ny = 256

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

#if AMREX_SPACEDIM == 1
  center(1) = ZERO

#elif AMREX_SPACEDIM == 2
  ! assume axisymmetric
  center(1) = ZERO
  center(2) = HALF*(problo(2)+probhi(2))

#else
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))
  center(3) = HALF*(problo(3)+probhi(3))
#endif

  call integrate_model(ny, problo(2), probhi(2), rho0, p0)

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
subroutine ca_initdata(lo, hi, &
                       state, s_lo, s_hi, &
                       dx, problo) bind(C, name='ca_initdata')

  use amrex_constants_module
  use probdata_module
  use interpolate_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP,&
                                 UEDEN, UEINT, UFS, T_guess
  use network, only : nspec
  use model_parser_module
  use eos_type_module
  use eos_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  real(rt), intent(in   ) :: dx(3), problo(3)

  real(rt) :: x, y, z, fheat, rhopert
  integer :: i, j, k, n

  type(eos_t) :: eos_state

  real(rt) :: pres

  !$gpu

  do k = lo(3), hi(3)
    z = problo(3) + dx(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + HALF)

        if (y < 1.125e0_rt * 4.e8_rt) then
            fheat = sin(8.e0_rt * M_PI * (y/ 4.e8_rt - ONE))
        else
            fheat = ZERO
        endif

        do i = lo(1), hi(1)
           x = problo(1) + dx(1)*(dble(i) + HALF)

           rhopert = ZERO

           if (do_pert) then
              rhopert = 5.e-5_rt * rho0 * fheat * (sin(3.e0_rt * M_PI * x / 4.e8_rt) + &
                                                   cos(M_PI * x / 4.e8_rt)) * &
                                                   (sin(3 * M_PI * z/4.e8_rt) - cos(M_PI * z/4.e8_rt))
           end if

           call interpolate_sub(state(i,j,k,URHO), y, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), y, itemp_model)

           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), y, ispec_model-1+n)
           end do

           ! get temporary pressure
           call interpolate_sub(pres, y, ipres_model)

           eos_state % rho = state(i,j,k,URHO) + rhopert
           eos_state % p = pres
           eos_state % xn(:) = state(i,j,k,UFS:UFS-1+nspec)
           eos_state % T = T_guess

           call eos(eos_input_rp, eos_state)

           state(i,j,k,URHO) = eos_state % rho
           state(i,j,k,UTEMP) = eos_state % T

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state%e
           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state%e

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

           ! Initial velocities = 0
           state(i,j,k,UMX:UMZ) = ZERO

        enddo
     enddo
  enddo


end subroutine ca_initdata
