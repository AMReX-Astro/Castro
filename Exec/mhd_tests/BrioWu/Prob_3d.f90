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

  namelist /fortin/ p_l, u_l, rho_l, p_r, u_r, rho_r, frac, &
       B_x_l, B_y_l, B_z_l, B_x_r, B_y_r, B_z_r, idir

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
  p_l = 1.0               ! left pressure (erg/cc)
  u_l = 0.0               ! left velocity (cm/s)
  rho_l = 1.0             ! left density (g/cc)
  B_x_l = 0.75
  B_y_l = 0.
  B_z_l = 1.

  p_r = 0.1               ! right pressure (erg/cc)
  u_r = 0.0               ! right velocity (cm/s)
  rho_r = 0.125           ! right density (g/cc)
  B_x_r = 0.75
  B_y_r = 0.
  B_z_r = -1.

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

  ! compute the internal energy (erg/cc) for the left and right state
  eos_state%rho = rho_l
  eos_state%p = p_l
  eos_state%T = 100000.e0_rt  ! initial guess
  eos_state%xn(:) = xn_zone(:)

  call eos(eos_input_rp, eos_state)

  rhoe_l = rho_l*eos_state%e
  T_l = eos_state%T

  eos_state%rho = rho_r
  eos_state%p = p_r
  eos_state%T = 100000.e0_rt  ! initial guess
  eos_state%xn(:) = xn_zone(:)

  call eos(eos_input_rp, eos_state)

  rhoe_r = rho_r*eos_state%e
  T_r = eos_state%T

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

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: state(state_l1:state_h1, &
                    state_l2:state_h2, &
                    state_l3:state_h3,NVAR)

  
  real(rt) :: xcen, ycen, zcen
  integer :: i,j,k

  do k = lo(3), hi(3)
     zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5e0_rt)

     do j = lo(2), hi(2)
        ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5e0_rt)

        do i = lo(1), hi(1)
           xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt)

           if (idir == 1) then
              if (xcen <= center(1)) then
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = rho_l*u_l
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*u_l*u_l + 0.5e0_rt * (B_x_l**2 + B_y_l**2 + B_z_l**2)
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = rho_r*u_r
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_r + 0.5e0_rt*rho_r*u_r*u_r + 0.5e0_rt *(B_x_r**2 + B_y_r**2 + B_z_r**2)
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              endif
           endif


           if (idir == 2) then
              if (ycen <= center(2)) then
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = rho_l*u_l
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*u_l*u_l + 0.5e0_rt * (B_x_l**2 + B_y_l**2 + B_z_l**2)
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = rho_r*u_r
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_r + 0.5e0_rt*rho_r*u_r*u_r + 0.5e0_rt *(B_x_r**2 + B_y_r**2 + B_z_r**2)
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              endif
           endif


           if (idir == 3) then
              if (zcen <= center(3)) then
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = rho_l*u_l
                 state(i,j,k,UEDEN) = rhoe_l + 0.5e0_rt*rho_l*u_l*u_l + 0.5e0_rt * (B_x_l**2 + B_y_l**2 + B_z_l**2)
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = rho_r*u_r
                 state(i,j,k,UEDEN) = rhoe_r + 0.5e0_rt*rho_r*u_r*u_r + 0.5e0_rt *(B_x_r**2 + B_y_r**2 + B_z_r**2)
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              endif
           endif
!              state(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
              state(i,j,k,UFS  ) = state(i,j,k,URHO)

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

  if (idir .eq. 1) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
              if (xcen <= center(1)+1) then
                 mag_x(i,j,k,1) = B_x_l
              else
                 mag_x(i,j,k,1) = B_x_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
              if (xcen <= center(1)) then
                 mag_y(i,j,k,1) = B_y_l
              else
                 mag_y(i,j,k,1) = B_y_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
              if (xcen <= center(1)) then
                 mag_z(i,j,k,1) = B_z_l
              else
                 mag_z(i,j,k,1) = B_z_r
              endif
           enddo
        enddo
     enddo

  endif




  if (idir .eq. 2) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
           do i = lo(1), hi(1)+1
              if (ycen <= center(2)) then
                 mag_x(i,j,k,1) = B_x_l
              else
                 mag_x(i,j,k,1) = B_x_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)+1
           ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
           do i = lo(1), hi(1)
              if (ycen <= center(2)+1) then
                 mag_y(i,j,k,1) = B_y_l
              else
                 mag_y(i,j,k,1) = B_y_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        do j = lo(2), hi(2)
           ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
           do i = lo(1), hi(1)
              if (ycen <= center(2)) then
                 mag_z(i,j,k,1) = B_z_l
              else
                 mag_z(i,j,k,1) = B_z_r
              endif
           enddo
        enddo
     enddo

  endif



  if (idir .eq. 3) then
     do k = lo(3), hi(3)
        zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              if (zcen <= center(3)) then
                 mag_x(i,j,k,1) = B_x_l
              else
                 mag_x(i,j,k,1) = B_x_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)
        zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              if (zcen <= center(3)) then
                 mag_y(i,j,k,1) = B_y_l
              else
                 mag_y(i,j,k,1) = B_y_r
              endif
           enddo
        enddo
     enddo

     do k = lo(3), hi(3)+1
        zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              if (zcen <= center(3)+1) then
                 mag_z(i,j,k,1) = B_z_l
              else
                 mag_z(i,j,k,1) = B_z_r
              endif
           enddo
        enddo
     enddo

  endif

end subroutine ca_initmag
