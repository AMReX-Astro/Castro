subroutine amrex_probinit(init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use amrex_constants_module
  use castro_error_module
  use fundamental_constants_module
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin
  integer :: i

  namelist /fortin/ density, diameter, ambient_dens, problem

  integer, parameter :: maxlen=127
  character :: probin*(maxlen)
  character :: model*(maxlen)

  ! Build "probin" filename -- the name of file containing fortin namelist.
  if (namlen > maxlen) then
     call castro_error("ERROR: probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  density  = 1.0e0_rt
  diameter = 1.0e0_rt

  ambient_dens = 1.0e-8_rt

  problem = 1

  ! Read namelists -- override the defaults

  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

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
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
                                 UEDEN, UEINT, UFS, UFA
  use network, only : nspec
  use prob_params_module, only: problo, probhi, center

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

  real(rt) :: xx, yy, zz

  integer :: i, j, k

  !$OMP PARALLEL DO PRIVATE(i, j, k, xx, yy, zz)
  do k = lo(3), hi(3)
     zz = xlo(3) + delta(3) * (dble(k-lo(3))+HALF) - center(3)

     do j = lo(2), hi(2)
        yy = xlo(2) + delta(2) * (dble(j-lo(2))+HALF) - center(2)

        do i = lo(1), hi(1)
           xx = xlo(1) + delta(1) * (dble(i-lo(1))+HALF) - center(1)

           ! Establish the cube or sphere

           if (problem .eq. 1 .or. problem .eq. 2) then

              if ((xx**2 + yy**2 + zz**2)**0.5 < diameter / 2) then
                 state(i,j,k,URHO) = density
              else
                 state(i,j,k,URHO) = ambient_dens
              endif

           else if (problem .eq. 3) then

              if (abs(xx) < diameter/2 .and. abs(yy) < diameter/2 .and. abs(zz) < diameter/2) then
                 state(i,j,k,URHO) = density
              else
                 state(i,j,k,URHO) = ambient_dens
              endif

           else

              call castro_error("Problem not defined.")

           endif

           ! Establish the thermodynamic quantities. They don't have to be
           ! valid because this test will never do a hydro step.

           state(i,j,k,UTEMP) = 1.0e0_rt
           state(i,j,k,UEINT) = 1.0e0_rt
           state(i,j,k,UEDEN) = 1.0e0_rt

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO) / nspec

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
