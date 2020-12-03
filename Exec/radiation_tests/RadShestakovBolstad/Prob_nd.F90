subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use probdata_module
  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

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

  use fundamental_constants_module, only : hplanck, k_B, ev2erg, c_light
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use prob_params_module, only : problo

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), nscal)

  real(rt) :: xcen
  integer :: i, j, k
  real(rt) :: Tcgs, B0, nu0, l0, x0, u0
  real(rt) :: pi, rhoe_0

  Tcgs = T_0 * 1.e3_rt * ev2erg / k_B

  pi = 4.0e0_rt*atan(1.0e0_rt)
  B0 = 8.0*pi*hplanck/c_light**3
  nu0 = k_B*Tcgs/hplanck
  l0 = nu0**3/kappa_0
  x0 = l0/sqrt(3.e0_rt)
  u0 = B0*nu0**3

  ! rhoe_0 = R * k_B * Tcgs * u0 / hplanck
  rhoe_0 = 99968636.6828e0_rt * Tcgs * rho_0
  !          cv            * T    * rho

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           xcen = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

           state(i,j,k,URHO) = rho_0
           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! set the composition to be all in the first species
           state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
           state(i,j,k,UFS) = state(i,j,k,URHO)
           if (naux > 0) then
              state(i,j,k,UFX) = state(i,j,k,URHO)
           end if

           if (abs(xcen)/x0 .lt. x_jump) then
              state(i,j,k,UEDEN) = rhoe_0
              state(i,j,k,UEINT) = rhoe_0
              state(i,j,k,UTEMP) = Tcgs
           else
              state(i,j,k,UEDEN) = rhoe_0 * 1.d-12
              state(i,j,k,UEINT) = rhoe_0 * 1.d-12
              state(i,j,k,UTEMP) = Tcgs * 1.d-12
           end if

        end do
     end do
  end do

end subroutine ca_initdata

! :::
! ::: -----------------------------------------------------------
! :::

subroutine ca_initrad(level, time, lo, hi, nrad, &
                      rad_state, rad_state_lo, rad_state_hi, &
                      delta, xlo, xhi)

  use probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nrad
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: rad_state_lo(3), rad_state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: rad_state(rad_state_lo(1):rad_state_hi(1), &
                                       rad_state_lo(2):rad_state_hi(2), &
                                       rad_state_lo(3):rad_state_hi(3), 0:nrad-1)

  ! local variables
  integer :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rad_state(i,j,k,:) = 0.e0_rt
        end do
     end do
  end do

end subroutine ca_initrad
