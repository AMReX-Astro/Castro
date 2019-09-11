subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use model_parser_module
  use amrex_error_module
  use prob_params_module, only : center
  use probdata_module, only : heating_factor, g0, rho0, p0, gamma1

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  namelist /fortin/ &
       heating_factor, g0, rho0, p0, gamma1, do_pert

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 127
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  heating_factor = 1.e3_rt
  g0 = -9.021899571e8_rt
  rho0 = 1.82094e6_rt
  p0 = 2.7647358e23_rt
  gamma1 = 1.4e0_rt
  do_pert = .true.

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

  use amrex_constants_module
  use probdata_module
  use interpolate_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP,&
                                 UEDEN, UEINT, UFS, sdc_order
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center, problo, probhi
  use eos_type_module
  use eos_module
  use prescribe_grav_module, only : grav_zone
  use amrex_fort_module, only : rt => amrex_real
  use model_util_module, only : set_species, fv, dUdy

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), NVAR)

  real(rt) :: x, y, z, fheat, rhopert, U_old(2), U_new(2), h
  real(rt) :: ystart, k1(2), k2(2), k3(2), k4(2)
  real(rt), allocatable :: pres(:), dens(:)
  real(rt) :: xn(nspec)
  integer :: i, j, k, n, iter, n_dy
  real(rt), parameter :: TOL = 1.e-10_rt
  integer, parameter :: MAX_ITER = 200
  logical :: converged_hse

  type(eos_t) :: eos_state

  allocate(pres(0:hi(2)))
  allocate(dens(0:hi(2)))

  U_old(1) = log(rho0)
  U_old(2) = log(p0)

  if (sdc_order /= 4) then

     ! do HSE using RK2
     do j = 0, hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        ! our integration starts at y - h
        if (j .eq. 0) then
           h = delta(2) * HALF
        else
           h = delta(2)
        endif

        k1(:) = dUdy(y - h, U_old)
        U_new(:) = U_old(:) + h * dUdy(y - HALF*h, U_old + HALF*h * k1)

        dens(j) = exp(U_new(1))
        pres(j) = exp(U_new(2))

        U_old(:) = U_new(:)

     end do

  else

     ! do HSE using RK4
     do j = 0, hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        ! our integration starts at y - h
        if (j .eq. 0) then
           h = delta(2) * HALF
        else
           h = delta(2)
        endif

        ystart = y - h

        k1(:) = dUdy(ystart, U_old)
        U_new(:) = U_old(:) + HALF*h * k1(:)

        k2(:) = dUdy(ystart + HALF*h, U_new)
        U_new(:) = U_old(:) + HALF*h * k2(:)

        k3(:) = dUdy(ystart + HALF*h, U_new)
        U_new(:) = U_old(:) + h * k3(:)

        k4(:) = dUdy(ystart + h, U_new)

        U_new = U_old(:) + (1.0_rt/6.0_rt) * h * (k1(:) + TWO*k2(:) + TWO*k3(:) + k4(:))

        dens(j) = exp(U_new(1))
        pres(j) = exp(U_new(2))

        U_old(:) = U_new(:)

     end do

  end if

  do k = lo(3), hi(3)
    z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF) - center(3)

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

        if (y < 1.125e0_rt * 4.e8_rt) then
            fheat = sin(8.e0_rt * M_PI * (y/ 4.e8_rt - ONE))
        else
            fheat = ZERO
        endif

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF) - center(1)

           rhopert = ZERO

           if (do_pert) then
              rhopert = 5.e-5_rt * rho0 * fheat * (sin(3.e0_rt * M_PI * x / 4.e8_rt) + &
                                                   cos(M_PI * x / 4.e8_rt)) * &
                                                   (sin(3 * M_PI * z/4.e8_rt) - cos(M_PI * z/4.e8_rt))
           end if

           ! do species
           state(i,j,k,UFS:UFS-1+nspec) = set_species(y)

           if (j < 0) then
                eos_state % rho = rho0 + rhopert
                eos_state % p = p0
           else
                eos_state%rho = dens(j) + rhopert
                eos_state%p = pres(j)
           endif

           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rp, eos_state)

           state(i,j,k,URHO) = eos_state % rho

           state(i,j,k,UTEMP) = eos_state%T

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state%e
           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state%e

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        enddo
     enddo
  enddo

  deallocate(pres)
  deallocate(dens)

  ! Initial velocities = 0
  state(:,:,:,UMX:UMZ) = ZERO

end subroutine ca_initdata
