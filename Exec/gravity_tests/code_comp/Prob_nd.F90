subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use model_parser_module
  use amrex_error_module
  use prob_params_module, only : center
  use probdata_module, only : heating_factor, g0, rho0, p0

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  namelist /fortin/ &
       heating_factor, g0, rho0, p0

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
  use eos_module
  use actual_eos_module, only : gamma_const
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP,&
                                 UEDEN, UEINT, UFS
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center, problo, probhi
  use eos_type_module
  use eos_module
  use prescribe_grav_module, only : grav_zone
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), NVAR)

  real(rt) :: x, y, z, p_want, dens_zone, sum_X, fheat, rhopert, temp_zone
  real(rt) :: pres_zone, dpd, drho, g_zone, fv, fg, gp, gm
  real(rt), allocatable :: pres(:), dens(:)
  real(rt) :: xn(nspec)
  integer :: i, j, k, n, iter, n_dy
  real(rt), parameter :: TOL = 1.e-10_rt
  integer, parameter :: MAX_ITER = 200
  logical :: converged_hse

  type(eos_t) :: eos_state

  allocate(pres(0:hi(2)))
  allocate(dens(0:hi(2)))

  ! do HSE
  do j = 0, hi(2)
    y = problo(2) + delta(2)*(dble(j) + HALF) 

    if (j .eq. 0) then 
        dens_zone = rho0
        pres_zone = p0
    else
        dens_zone = dens(j-1)
        pres_zone = pres(j-1)
    endif

    eos_state%rho = dens_zone
    eos_state%p = pres_zone

    ! do species
    xn(:) = ZERO
    if (y < 1.9375e0_rt * 4.e8_rt) then 
         xn(1) = ONE
         xn(2) = ZERO
    else if (y > 2.0625e0_rt * 4.e8_rt) then
         xn(1) = ZERO 
         xn(2) = ONE 
    else
         fv = HALF * (ONE + sin(8.e0_rt * M_PI * (y/4.e8_rt - 2.e0_rt)))
         xn(1) = ONE - fv
         xn(2) = fv
    endif
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rp, eos_state)

    temp_zone = eos_state % T

    if (j .eq. 0) then
       ! compute the gravitational acceleration halfway between lower boundary 
       ! and cell center
       g_zone = grav_zone(y-delta(2)*0.25e0_rt)

       ! compute the gravitational acceleration on the lower boundary 
       gm = grav_zone(y-delta(2)*HALF)

    else
       ! compute the gravitational acceleration on the interface between zones
       ! i and i+1
       g_zone = grav_zone(y-delta(2)*HALF)

       ! compute the gravitational acceleration in the cell below
       gm = grav_zone(y-delta(2))
    endif

    ! compute the gravitational acceleration at cell center
    gp = grav_zone(y)

    converged_hse = .FALSE.

    do iter = 1, MAX_ITER

        if (j .eq. 0) then 
            p_want = p0 + HALF * delta(2) * HALF * (dens_zone*gp + rho0*gm) 
        else
            p_want = pres(j-1) + delta(2) * HALF * (dens_zone*gp + dens(j-1)*gm)
        endif

        ! (t, rho) -> (p, s)
        eos_state%T     = temp_zone
        eos_state%rho   = dens_zone
        eos_state%xn(:) = xn(:)

        ! do species
        eos_state%xn(:) = ZERO
        if (y < 1.9375e0_rt * 4.e8_rt) then 
             eos_state%xn(1) = ONE
             eos_state%xn(2) = ZERO
        else if (y > 2.0625e0_rt * 4.e8_rt) then
             eos_state%xn(1) = ZERO 
             eos_state%xn(2) = ONE 
        else
             fv = HALF * (ONE + sin(8.e0_rt * M_PI * (y/4.e8_rt - 2.e0_rt)))
             eos_state%xn(1) = ONE - fv
             eos_state%xn(2) = fv
        endif

        call eos(eos_input_rt, eos_state)

        pres_zone = eos_state%p
        temp_zone = eos_state%T

        dpd = eos_state%dpdr
        drho = (p_want - pres_zone) / (dpd - HALF*delta(2)*g_zone)

        dens_zone = max(0.9e0_rt*dens_zone, &
                min(dens_zone + drho, 1.1e0_rt*dens_zone))

        if (abs(drho) < TOL*dens_zone) then
            converged_hse = .TRUE.
            exit
        endif

    enddo

    if (.NOT. converged_hse) then

        print *, 'Error zone', j, ' did not converge in init_1d', y/4.e8_rt
        print *, 'integrate up'
        print *, 'dens_zone, temp_zone = ', dens_zone, temp_zone
        print *, "p_want = ", p_want
        print *, "drho = ", drho
        call bl_error('Error: HSE non-convergence')
    endif

    ! call the EOS one more time for this zone and then go on to the next
        ! (t, rho) -> (p, s)

    eos_state%T     = temp_zone
    eos_state%rho   = dens_zone
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rt, eos_state)

    pres(j) = eos_state%p
    dens(j) = dens_zone

  enddo

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

           rhopert = 5.e-5_rt * rho0 * fheat * (sin(3.e0_rt * M_PI * x / 4.e8_rt) + &
                                                cos(M_PI * x / 4.e8_rt)) * &
                     (sin(3 * M_PI * z/4.e8_rt) - cos(M_PI * z/4.e8_rt))

           ! do species
           state(i,j,k,UFS:UFS-1+nspec) = ZERO
           if (y < 1.9375e0_rt * 4.e8_rt) then 
                state(i,j,k,UFS) = ONE
                state(i,j,k,UFS+1) = ZERO
           else if (y > 2.0625e0_rt * 4.e8_rt) then
                state(i,j,k,UFS) = ZERO 
                state(i,j,k,UFS+1) = ONE 
           else
                fv = HALF * (ONE + sin(8.e0_rt * M_PI * (y/4.e8_rt - 2.e0_rt)))
                state(i,j,k,UFS) = ONE - fv
                state(i,j,k,UFS+1) = fv
           endif

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

