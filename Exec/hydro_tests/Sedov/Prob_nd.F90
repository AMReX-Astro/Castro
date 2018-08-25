subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use prob_params_module, only : center, coord_type
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

  namelist /fortin/ p_ambient, dens_ambient, exp_energy, &
       r_init, nsub, temp_ambient

  real(rt) :: vctr

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

  p_ambient = 1.e-5_rt        ! ambient pressure (in erg/cc)
  dens_ambient = 1.e0_rt      ! ambient density (in g/cc)
  exp_energy = 1.e0_rt        ! absolute energy of the explosion (in erg)
  r_init = 0.05e0_rt          ! initial radius of the explosion (in cm)
  nsub = 4
  temp_ambient = -1.e2_rt     ! Set original temp. to negative, which is overwritten in the probin file

  ! set local variable defaults
  if (coord_type == 1 .or. coord_type == 2) then
     center = ZERO
  else
     center = HALF * (problo + probhi)
  end if

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

  xn_zone(:) = ZERO
  xn_zone(1) = ONE

  ! override the pressure with the temperature
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
  eos_state % T   = 1.d9 ! Initial guess for iterations
  eos_state % xn  = xn_zone

  call eos(eos_input_rp, eos_state)

  e_ambient = eos_state % e

  ! set explosion pressure -- we will convert the point-explosion energy into
  ! a corresponding pressure distributed throughout the perturbed volume

  if (coord_type == 0) then

#if AMREX_SPACEDIM == 1

#ifndef AMREX_USE_CUDA
     call amrex_error("Sedov problem unsupported in 1D Cartesian geometry.")
#endif

#elif AMREX_SPACEDIM == 2

     ! Cylindrical problem in Cartesian coordinates

     vctr = M_PI*r_init**2

#else

     ! Spherical problem in Cartesian coordinates

     vctr = FOUR3RD*M_PI*r_init**3

#endif

  else if (coord_type == 1) then

#if AMREX_SPACEDIM == 1

     vctr = M_PI*r_init**2

#elif AMREX_SPACEDIM == 2

     vctr = FOUR3RD*M_PI*r_init**3

#else

#ifndef AMREX_USE_CUDA
     call amrex_error("Sedov problem unsupported in 3D axisymmetric geometry.")
#endif

#endif

  else if (coord_type == 2) then

     ! Must have AMREX_SPACEDIM == 1 for this coord_type.

     vctr = FOUR3RD*M_PI*r_init**3

  end if

  e_exp = exp_energy/vctr/dens_ambient

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
  use amrex_constants_module, only: ZERO, ONE, HALF, FOUR3RD, M_PI
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS
  use prob_params_module, only : center, coord_type, dim
  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rp, eos_input_re
  use amrex_error_module, only: amrex_error

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt) :: xlo(3), xhi(3), time, delta(3)
  real(rt) :: state(state_lo(1):state_hi(1), &
                    state_lo(2):state_hi(2), &
                    state_lo(3):state_hi(3),NVAR)
  real(rt) :: xmin,ymin,zmin
  real(rt) :: xx, yy, zz
  real(rt) :: dist
  real(rt) :: eint, e_zone

  integer :: i,j,k, ii, jj, kk
  type(eos_t) :: eos_state
  integer :: nsubx, nsuby, nsubz
  real(rt) :: ds(3), xl, xr
  real(rt) :: vol_pert, vol_ambient
  logical :: pert

  if (dim == 1) then
     nsubx = nsub
     nsuby = 1
     nsubz = 1
  else if (dim == 2) then
     nsubx = nsub
     nsuby = nsub
     nsubz = 1
  else
     nsubx = nsub
     nsuby = nsub
     nsubz = nsub
  end if

  ds(1) = delta(1) / nsubx
  ds(2) = delta(2) / nsuby
  ds(3) = delta(3) / nsubz

  do k = lo(3), hi(3)
     zmin = xlo(3) + delta(3)*dble(k-lo(3))

     do j = lo(2), hi(2)
        ymin = xlo(2) + delta(2)*dble(j-lo(2))

        do i = lo(1), hi(1)
           xmin = xlo(1) + delta(1)*dble(i-lo(1))

           vol_pert = ZERO
           vol_ambient = ZERO

           do kk = 0, nsubz-1
              zz = zmin + ds(3) * (kk + HALF)

              do jj = 0, nsuby-1
                 yy = ymin + ds(2) * (jj + HALF)

                 do ii = 0, nsubx-1
                    xx = xmin + ds(1) * (ii + HALF)

                    dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2

                    if (dist <= r_init**2) then
                       pert =  .true.
                    else
                       pert = .false.
                    end if

                    if (coord_type == 1) then

                       ! The volume of a cell is a annular cylindrical region.
                       ! The main thing that matters is the distance from the
                       ! symmetry axis.
                       !   V = pi*dy*(x_r**2 - x_l**2) = pi*dy*dx*HALF*xx
                       ! (where x_r is the coordinate of the x right edge,
                       !        x_l is the coordinate of the x left edge,
                       !    and xx  is the coordinate of the x center of the cell)
                       !
                       ! since dx and dy are constant, they cancel out

                       if (pert) then
                          vol_pert = vol_pert + xx
                       else
                          vol_ambient = vol_ambient + xx
                       end if

                    else if (coord_type == 2) then

                       xl = xx - HALF * ds(1)
                       xr = xx + HALF * ds(1)

                       ! the volume of a subzone is (4/3) pi (xr^3 - xl^3).
                       ! we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)
                       ! The (4/3) pi dr factor is common, so we can neglect it.

                       if (pert) then
                          vol_pert = vol_pert + (xr*xr + xl*xr + xl*xl)
                       else
                          vol_ambient = vol_ambient + (xr*xr + xl*xr + xl*xl)
                       endif

                    else

                       ! Cartesian -- equal volume of dx * dy * dz is factored out.

                       if (pert) then
                          vol_pert = vol_pert + ONE
                       else
                          vol_ambient = vol_ambient + ONE
                       end if

                    end if

                 enddo
              enddo
           enddo

           e_zone = (vol_pert * e_exp + vol_ambient * e_ambient) / (vol_pert + vol_ambient)

           eint = dens_ambient * e_zone

           state(i,j,k,URHO) = dens_ambient
           state(i,j,k,UMX) = 0.e0_rt
           state(i,j,k,UMY) = 0.e0_rt
           state(i,j,k,UMZ) = 0.e0_rt

           state(i,j,k,UEDEN) = eint + &
                0.5e0_rt*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMZ)**2/state(i,j,k,URHO))

           state(i,j,k,UEINT) = eint

           ! The temperature initialization will be done later
           ! in the call to clean_state. We want to avoid
           ! EOS calls in ca_initdata.

           state(i,j,k,UFS) = state(i,j,k,URHO)

        enddo
     enddo
  enddo

end subroutine ca_initdata
