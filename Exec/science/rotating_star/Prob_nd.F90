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
                                 UEDEN, UEINT, UFS, UFX
  use network, only : nspec, naux
  use model_parser_module
  use prob_params_module, only : center, problo, probhi
  use eos_type_module
  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use amrex_fort_module, only : rt => amrex_real
  use burn_type_module
#ifdef NSE_THERMO
  use nse_check_module
  use nse_module, only: nse_interp
#endif

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: x, y, z, dist, pres, r1, t0, zc
  integer :: i, j, k, n

  type(eos_t) :: eos_state
  real(rt) :: sumX
  real(rt) :: abar, dq, dyedt
  real(rt) :: xn(nspec)
  type(burn_t) :: burn_state
  integer :: nse_check

  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF) - center(3)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF) - center(2)

        do i = lo(1), hi(1)
           x = problo(1) + delta(1)*(dble(i) + HALF) - center(1)

           dist = sqrt(x**2 + y**2 + z**2)

           call interpolate_sub(state(i,j,k,URHO), dist, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), dist, itemp_model)

           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), dist, ispec_model-1+n)
           end do

           sumX = 0.0_rt
           do n = 1, nspec
              sumX = sumX + state(i,j,k,UFS-1+n)
           end do
           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec) / sumX

#ifdef NSE_THERMO
           ! set the aux quantities -- we need to do this if we are using the NSE network
           state(i,j,k,UFX-1+iye) = sum(state(i,j,k,UFS:UFS-1+nspec) * zion(:) * aion_inv(:))
           state(i,j,k,UFX-1+iabar) = 1.0_rt / sum(state(i,j,k,UFS:UFS-1+nspec) * aion_inv(:))
           state(i,j,k,UFX-1+ibea) = sum(state(i,j,k,UFS:UFS-1+nspec) * bion(:) * aion_inv(:))

           burn_state % rho = state(i,j,k,URHO)
           burn_state % T = state(i,j,k,UTEMP)
           burn_state % xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call in_nse(burn_state, nse_check)
           if (nse_check == 1) then
              call nse_interp(state(i,j,k,UTEMP), state(i,j,k,URHO), &
                              state(i,j,k,UFS-1+iye), abar, dq, dyedt, xn)

              state(i,j,k,UFX-1+iabar) = abar
              state(i,j,k,UFX-1+ibea) = dq
              state(i,j,k,UFS:UFS-1+nspec) = xn(:)

              sumX = 0.0_rt
              do n = 1, nspec
                 sumX = sumX + state(i,j,k,UFS-1+n)
              end do
              state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec) / sumX

           end if
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
           eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state % e

           do n = 1, nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

           do n = 1, naux
              state(i,j,k,UFX+n-1) = state(i,j,k,URHO) * state(i,j,k,UFX+n-1)
           end do

           ! initial velocities = 0
           state(i,j,k,UMX:UMZ) = ZERO

        end do
     end do
  end do


end subroutine ca_initdata
