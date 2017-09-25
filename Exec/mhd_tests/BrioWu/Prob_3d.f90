subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_constants_module
  use probdata_module
  use prob_params_module, only : center
  use bl_error_module
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_rp
  use eos_module, only: eos
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i

  type(eos_t) :: eos_state

  namelist /fortin/ p_ambient, dens_ambient, exp_energy, &
       r_init, nsub, temp_ambient

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  p_ambient = 1.e-5_rt        ! ambient pressure (in erg/cc)
  dens_ambient = 1.e0_rt      ! ambient density (in g/cc)
  exp_energy = 1.e0_rt        ! absolute energy of the explosion (in erg)
  r_init = 0.05e0_rt          ! initial radius of the explosion (in cm)
  nsub = 4
  temp_ambient = -1.e2_rt     ! Set original temp. to negative, which is overwritten in the probin file

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

  ! override the pressure iwth the temperature
  if (temp_ambient > ZERO) then

     eos_state % rho = dens_ambient
     eos_state % xn(:) = xn_zone(:)
     eos_state % T = temp_ambient

     call eos(eos_input_rt, eos_state)

     p_ambient = eos_state % p

  endif

  ! Calculate ambient state data

  eos_state % rho = dens_ambient
  eos_state % p   = p_ambient
  eos_state % T   = 1.d5 ! Initial guess for iterations
  eos_state % xn  = xn_zone

  call eos(eos_input_rp, eos_state)

  e_ambient = eos_state % e

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
  use bl_constants_module, only: M_PI, FOUR3RD, ZERO, ONE
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp, eos_input_re

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: state(state_l1:state_h1, &
                    state_l2:state_h2, &
                    state_l3:state_h3,NVAR)

  real(rt) :: xmin,ymin,zmin
  real(rt) :: xx, yy, zz
  real(rt) :: dist
  real(rt) :: eint, p_zone
  real(rt) :: vctr, p_exp

  integer :: i,j,k, ii, jj, kk
  integer :: npert, nambient
  real(rt) :: e_zone
  type(eos_t) :: eos_state

  ! set explosion pressure -- we will convert the point-explosion energy into
  ! a corresponding pressure distributed throughout the perturbed volume
  vctr  = FOUR3RD*M_PI*r_init**3

  e_zone = exp_energy/vctr/dens_ambient

  eos_state % e = e_zone
  eos_state % rho = dens_ambient
  eos_state % xn(:) = xn_zone(:)
  eos_state % T = 100.0  ! initial guess

  call eos(eos_input_re, eos_state)

  p_exp = eos_state % p

  do k = lo(3), hi(3)
     zmin = xlo(3) + delta(3)*dble(k-lo(3))

     do j = lo(2), hi(2)
        ymin = xlo(2) + delta(2)*dble(j-lo(2))

        do i = lo(1), hi(1)
           xmin = xlo(1) + delta(1)*dble(i-lo(1))

           npert = 0
           nambient = 0

           do kk = 0, nsub-1
              zz = zmin + (delta(3)/dble(nsub))*(kk + 0.5e0_rt)

              do jj = 0, nsub-1
                 yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5e0_rt)

                 do ii = 0, nsub-1
                    xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5e0_rt)

                    dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2

                    if(dist <= r_init**2) then
                       npert = npert + 1
                    else
                       nambient = nambient + 1
                    endif

                 enddo
              enddo
           enddo

           p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient)/  &
                dble(nsub*nsub*nsub)

           eos_state % p = p_zone
           eos_state % rho = dens_ambient
           eos_state % xn(:) = xn_zone(:)

           call eos(eos_input_rp, eos_state)

           eint = dens_ambient * eos_state % e

           state(i,j,k,URHO) = dens_ambient
           state(i,j,k,UMX) = 0.e0_rt
           state(i,j,k,UMY) = 0.e0_rt
           state(i,j,k,UMZ) = 0.e0_rt

           state(i,j,k,UEDEN) = eint + &
                0.5e0_rt*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMZ)**2/state(i,j,k,URHO))

           state(i,j,k,UEINT) = eint

           state(i,j,k,UFS) = state(i,j,k,URHO)

        enddo
     enddo
  enddo

end subroutine ca_initdata
