subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use prob_params_module, only : center
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos_on_host
  use extern_probin_module, only : small_x
  use network, only : nspec

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type(eos_t) :: eos_state

  call probdata_init(name, namlen)

  ! set explosion center
  center(:) = ZERO
  center(1) = HALF*(problo(1) + probhi(1))
#if BL_SPACEDIM >= 2
  center(2) = HALF*(problo(2) + probhi(2))
#endif
#if BL_SPACEDIM == 3
  center(3) = HALF*(problo(3) + probhi(3))
#endif

  xn_zone(:) = small_x
  xn_zone(1) = ONE - (nspec-1)*small_x

  eos_state % rho = rho0
  eos_state % T = T0
  eos_state % xn(:) = xn_zone(:)

  call eos_on_host(eos_input_rt, eos_state)

  p0 = eos_state % p
  s0 = eos_state % s

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
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_lo,state_hi, &
                       delta,xlo,xhi)

  use probdata_module
  use amrex_constants_module, only: M_PI, FOUR3RD, ZERO, HALF, ONE
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only : center, coord_type, problo, probhi
  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec
  use eos_type_module, only : eos_t, eos_input_ps
  use eos_module, only : eos

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: xx, yy, zz
  real(rt) :: dist, p, eint

  integer :: i,j, k

  type(eos_t) :: eos_state

  do k = lo(3), hi(3)
     zz = problo(3) + delta(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        yy = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)
           xx = problo(1) + delta(1)*(dble(i) + HALF)

           dist = sqrt((center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2)

           if (dist <= center(1)) then
              p = p0*(ONE + dp_fact*exp(-(dist/L_pert)**2) * cos(M_PI*(dist/(probhi(1)-problo(1))))**6)
           else
              p = p0
           endif

           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! we are isentropic, so find rho
           eos_state % p =  p
           eos_state % T = 1.e4_rt ! initial guess
           eos_state % rho = rho0  ! initial guess
           eos_state % s = s0
           eos_state % xn(:) = xn_zone(:)

           call eos(eos_input_ps, eos_state)

           state(i,j,k,URHO) = eos_state % rho

           state(i,j,k,UEDEN) = eos_state % rho * eos_state % e
           state(i,j,k,UEINT) = eos_state % rho * eos_state % e

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO)*xn_zone(:)

           state(i,j,k,UTEMP) = eos_state % T
        enddo
     enddo
  enddo

end subroutine ca_initdata
