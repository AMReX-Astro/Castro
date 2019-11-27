subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(C, name="amrex_probinit")

  use probdata_module
  use amrex_fort_module, only : rt => amrex_real
  use castro_error_module, only: castro_error
  use eos_module
  use eos_type_module

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type(eos_t) :: eos_state

  call probdata_init(name, namlen)

  ! set some defaults we will need later in the BCs
  eos_state % rho = rho0
  eos_state % T   =   T0
  eos_state % xn  = 1.e0_rt

  call eos_on_host(eos_input_rt, eos_state)

  eint0 = rho0 * eos_state % e
  etot0 = eint0 + 0.5*rho0*v0**2

  eos_state % rho = rho1
  eos_state % T   =   T1

  call eos_on_host(eos_input_rt, eos_state)

  eint1 = rho1 * eos_state % e
  etot1 = eint1 + 0.5*rho1*v1**2

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
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use castro_error_module

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

  integer :: i, j, k
  real(rt) :: length_cell, rhoInv
  type(eos_t) :: eos_state

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           !Provides the simulation to be run in the x,y,or z direction
           !where length direction is the length side in a square prism
           if (idir == 1) then
              length_cell = problo(1) + delta(1) * (dble(i) + 0.5e0_rt)
           else if (idir == 2) then
              length_cell = problo(2) + delta(2) * (dble(j) + 0.5e0_rt)
           else if (idir == 3) then
              length_cell = problo(3) + delta(3) * (dble(k) + 0.5e0_rt)
           else
              call castro_error("Invalid direction please input idir = [1,3]")
           endif



           if (length_cell < 0.e0_rt) then
              state(i,j,k,URHO) = rho0

              ! set the composition to be all in the first species
              state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
              state(i,j,k,UFS  ) = state(i,j,k,URHO)

              state(i,j,k,UTEMP) = T0
              if (idir == 1) then
                 state(i,j,k,UMX) = rho0*v0
                 state(i,j,k,UMY) = 0.0e0_rt
                 state(i,j,k,UMZ) = 0.0e0_rt
              else if (idir == 2) then
                 state(i,j,k,UMX) = 0.0e0_rt
                 state(i,j,k,UMY) = rho0*v0
                 state(i,j,k,UMZ) = 0.0e0_rt
              else if (idir == 3) then
                 state(i,j,k,UMX) = 0.0e0_rt
                 state(i,j,k,UMY) = 0.0e0_rt
                 state(i,j,k,UMZ) = rho0*v0
              end if
           else
              state(i,j,k,URHO) = rho1

              ! set the composition to be all in the first species
              state(i,j,k,UFS:UFS-1+nspec) = 0.e0_rt
              state(i,j,k,UFS  ) = state(i,j,k,URHO)
              state(i,j,k,UTEMP) = T1

              if (idir == 1) then
                 state(i,j,k,UMX) = rho1*v1
                 state(i,j,k,UMY) = 0.0e0_rt
                 state(i,j,k,UMZ) = 0.0e0_rt
              else if (idir == 2) then
                 state(i,j,k,UMX) = 0.0e0_rt
                 state(i,j,k,UMY) = rho1*v1
                 state(i,j,k,UMZ) = 0.0e0_rt
              else if (idir == 3) then
                 state(i,j,k,UMX) = 0.0e0_rt
                 state(i,j,k,UMY) = 0.0e0_rt
                 state(i,j,k,UMZ) = rho1*v1
              end if
           end if


           if (naux > 0) then
              state(i,j,k,UFX) = state(i,j,k,URHO)
           end if

           eos_state % rho = state(i,j,k,URHO)
           eos_state % T   = state(i,j,k,UTEMP)

           rhoInv = 1.e0_rt / state(i,j,k,URHO)
           eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
           eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                0.5e0_rt*(state(i,j,k,UMX)**2 + &
                          state(i,j,k,UMY)**2 + &
                          state(i,j,k,UMZ)**2)/state(i,j,k,URHO)
        enddo
     enddo
  enddo

end subroutine ca_initdata


! :::
! ::: -----------------------------------------------------------
! :::

subroutine ca_initrad(level, time, lo, hi, nrad, &
                      rad_state, rad_state_lo, rad_state_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use fundamental_constants_module, only: a_rad
  use rad_params_module, only : xnu
  use blackbody_module, only : BGroup
  use castro_error_module
  use prob_params_module, only : problo
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
  integer :: i, j, k, igroup
  real(rt)         length_cell, t

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           !Provides the simulation to be run in the x,y,or z direction
           !where length direction is the length side in a square prism

           if (idir == 1) then
              length_cell = problo(1) + delta(1) * (dble(i) + 0.5e0_rt)
           else if (idir == 2) then
              length_cell = problo(2) + delta(2) * (dble(j) + 0.5e0_rt)
           else if (idir == 3) then
              length_cell = problo(3) + delta(3) * (dble(k) + 0.5e0_rt)
           else
              call castro_error("Invalid direction please input idir = [1,3]")
           endif

           if (length_cell < 0.e0_rt) then
              T = T0
           else
              T = T1
           end if

           if (nrad == 1) then
              rad_state(i,j,k,:) = a_rad*T**4
           else
              do igroup=0,nrad-1
                 rad_state(i,j,k,igroup) = BGroup(T, xnu(igroup), xnu(igroup+1))
              end do
           end if

        enddo
     enddo
  enddo

end subroutine ca_initrad
