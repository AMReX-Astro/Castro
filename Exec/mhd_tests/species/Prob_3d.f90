subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use prob_params_module, only : center
  use amrex_error_module
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_rp
  use eos_module, only: eos
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  type(eos_t) :: eos_state

  namelist /fortin/ p, rho, u_x, u_y, u_z, frac, &
       B_x, B_y, B_z, idir

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  p    = 1.0             ! pressure (erg/cc)
  u_x  = 1.0             ! velocity (cm/s)
  u_y  = 1.0
  u_z  = 1.0
  rho  = 1.0             ! density (g/cc)
  B_x  = 0.
  B_y  = 0.
  B_z  = 0.


  idir = 1                ! direction across which to jump
  frac = 0.5              ! fraction of the domain for the interface

  ! set local variable defaults
  center(1) = HALF*(problo(1) + probhi(1))
  center(2) = HALF*(problo(2) + probhi(2))
  center(3) = HALF*(problo(3) + probhi(3))

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

  xn_zone(:) = ZERO
  xn_zone(1) = ONE

  ! compute the internal energy (erg/cc) 
  eos_state%rho = rho
  eos_state%p = p
  eos_state%T = 100000.e0_rt  ! initial guess
  eos_state%xn(:) = xn_zone(:)

  call eos(eos_input_rp, eos_state)

  rhoe = rho*eos_state%e
  T  = eos_state % T

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
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)

  use probdata_module
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec
  use amrex_constants_module, only : ZERO
  use prob_params_module, only : problo
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: state(state_l1:state_h1, &
                    state_l2:state_h2, &
                    state_l3:state_h3,NVAR)

  
  real(rt) :: xcen, ycen, zcen, r2
  integer :: i,j,k

  do k = lo(3), hi(3)
     zcen = problo(3) + delta(3)*(dble(k) + 0.5e0_rt)

     do j = lo(2), hi(2)
        ycen = problo(2) + delta(2)*(dble(j) + 0.5e0_rt)

        do i = lo(1), hi(1)
           xcen = problo(1) + delta(1)*(dble(i) + 0.5e0_rt)

           state(i,j,k,URHO) = rho
           state(i,j,k,UMX) = rho*u_x
           state(i,j,k,UMY) = rho*u_y
           state(i,j,k,UMZ) = rho*u_z
           state(i,j,k,UEDEN) = rhoe + 0.5e0_rt*rho*(u_x**2+u_y**2+u_z**2) + 0.5e0_rt * (B_x**2 + B_y**2 + B_z**2)
           state(i,j,k,UEINT) = rhoe
           state(i,j,k,UTEMP) = T
           
           r2 = ((xcen-0.5d0)**2 + (ycen-0.5d0)**2 + (zcen-0.5d0)**2) / 0.01d0
           state(i,j,k,UFS:UFS-1+nspec) = ZERO
           state(i,j,k,UFS)  = exp(-r2)
           state(i,j,k,UFS+1)= 1.d0 - exp(-r2)
           
        enddo
     enddo
  enddo

 end subroutine ca_initdata


! :::
! ::: --------------------------------------------------------------------
! ::: 
subroutine ca_initmag(level, time, lo, hi, &
                      nbx, mag_x, bx_lo, bx_hi, &
                      nby, mag_y, by_lo, by_hi, &
                      nbz, mag_z, bz_lo, bz_hi, &
                      delta, xlo, xhi)

  use probdata_module
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => amrex_real
  
  implicit none
  
  integer :: level, nbx, nby, nbz
  integer :: lo(3), hi(3)
  integer :: bx_lo(3), bx_hi(3)
  integer :: by_lo(3), by_hi(3)
  integer :: bz_lo(3), bz_hi(3)
  real(rt) :: xlo(3), xhi(3), time, delta(3)

  real(rt) :: mag_x(bx_lo(1):bx_hi(1), bx_lo(2):bx_hi(2), bx_lo(3):bx_hi(3), nbx)
  real(rt) :: mag_y(by_lo(1):by_hi(1), by_lo(2):by_hi(2), by_lo(3):by_hi(3), nby)
  real(rt) :: mag_z(bz_lo(1):bz_hi(1), bz_lo(2):bz_hi(2), bz_lo(3):bz_hi(3), nbz)

  real(rt) :: xcen, ycen, zcen
  integer  :: i, j, k
  
  print *, "Initializing magnetic field!!"

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
                 mag_x(i,j,k,1) = B_x
                 mag_y(i,j,k,1) = B_y
                 mag_z(i,j,k,1) = B_z
           enddo
        enddo
     enddo

end subroutine ca_initmag
