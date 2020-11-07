subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use probdata_module
  use network, only : network_init
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_re
  use network, only : nspec
  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module, only : castro_error
  use prob_params_module, only : center

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  real(rt) :: X_in(nspec), e_0
  type(eos_t) :: eos_state

  ! get T_0 corresponding to rhoe_0 and rho_0 through the EOS
  e_0 = rhoe_0/rho_0
  X_in(:) = 0.0
  X_in(1) = 1.e0_rt

  eos_state % rho = rho_0
  eos_state % e = e_0
  eos_state % xn = X_in

  call eos(eos_input_re, eos_state)

  T_0 = eos_state % T

  ! the 'center' variables are the location of the middle of the
  ! domain -- this is where we put the interface
  center(1) = 0.5e0_rt*(problo(1)+probhi(1))


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
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only : nspec
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : problo

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

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           xcen = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

           state(i,j,k,URHO) = rho_0
           state(i,j,k,UMX:UMZ) = 0.e0_rt
           state(i,j,k,UEDEN) = rhoe_0
           state(i,j,k,UEINT) = rhoe_0

           ! set the composition to be all in the first species
           state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
           state(i,j,k,UFS) = state(i,j,k,URHO)

           state(i,j,k,UTEMP) = T_0

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

  use probdata_module,   only: E_rad
  use rad_params_module, only: nugroup, dnugroup
  use rad_params_module, only: pi, clight, hplanck, kboltz, stefbol
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
  real(rt) :: trad, mgfac, nu, radtmp
  integer :: i, j, k, n

  if (nrad == 1) then

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rad_state(i,j,k,0) = E_rad
           end do
        end do
     end do

  else

     trad = (E_rad * clight * 0.25e0_rt / stefbol) ** (0.25e0_rt)

     mgfac = 8.e0_rt * pi * hplanck / clight**3

     do n = 0, nrad-1
        nu     = nugroup(n)
        radtmp = exp(hplanck * nu / (kboltz * trad)) - 1.e0_rt
        radtmp = (mgfac * nu**3 / radtmp) * dnugroup(n)

        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 rad_state(i,j,k,n) = radtmp
              end do
           end do
        end do

     end do

  end if

end subroutine ca_initrad
