subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use prob_params_module, only : center
  use bl_error_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(2), probhi(2)

  integer :: untin,i

  namelist /fortin/ probtype, p_ambient, dens_ambient, exp_energy, &
       r_init, nsub

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Set namelist defaults

  p_ambient = 1.e-5_rt        ! ambient pressure (in erg/cc)
  dens_ambient = 1.e0_rt      ! ambient density (in g/cc)
  exp_energy = 1.e0_rt        ! absolute energy of the explosion (in erg)
  r_init = 0.05e0_rt          ! initial radius of the explosion (in cm)
  nsub = 4

  !     Set explosion center
  center(1) = (problo(1)+probhi(1))/2.e0_rt
  center(2) = (problo(2)+probhi(2))/2.e0_rt

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
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
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module, only : gamma_const
  use bl_constants_module, only: M_PI, FOUR3RD
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => c_real
  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  real(rt)         :: xlo(2), xhi(2), time, delta(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt)         :: xmin,ymin
  real(rt)         :: xx, yy
  real(rt)         :: dist
  real(rt)         :: eint, p_zone
  real(rt)         :: vctr, p_exp

  integer :: i,j, ii, jj
  real(rt)         :: vol_pert, vol_ambient

  ! Cylindrical problem in Cartesian coordinates
  if (probtype .eq. 21) then

     ! set explosion pressure -- we will convert the point-explosion
     ! energy into a corresponding pressure distributed throughout the
     ! perturbed volume
     vctr = M_PI*r_init**2
     p_exp = (gamma_const - 1.e0_rt)*exp_energy/vctr

     do j = lo(2), hi(2)
        ymin = xlo(2) + delta(2)*dble(j-lo(2))

        do i = lo(1), hi(1)
           xmin = xlo(1) + delta(1)*dble(i-lo(1))

           vol_pert    = 0.e0_rt
           vol_ambient = 0.e0_rt

           do jj = 0, nsub-1
              yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5e0_rt)

              do ii = 0, nsub-1
                 xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5e0_rt)

                 dist = (center(1)-xx)**2 + (center(2)-yy)**2

                 if(dist <= r_init**2) then
                    vol_pert    = vol_pert    + 1.e0_rt
                 else
                    vol_ambient = vol_ambient + 1.e0_rt
                 endif

              enddo
           enddo

           p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/ (vol_pert + vol_ambient)

           eint = p_zone/(gamma_const - 1.e0_rt)

           state(i,j,URHO) = dens_ambient
           state(i,j,UMX:UMZ) = 0.e0_rt

           state(i,j,UEDEN) = eint +  &
                0.5e0_rt*(sum(state(i,j,UMX:UMZ)**2)/state(i,j,URHO))

           state(i,j,UEINT) = eint

           state(i,j,UFS) = state(i,j,URHO)

        enddo
     enddo


  ! Cylindrical problem in cylindrical coordinates
  else if (probtype .eq. 22) then

     !  set explosion pressure -- we will convert the point-explosion
     !  energy into a corresponding pressure distributed throughout
     !  the perturbed volume
     vctr = M_PI*r_init**2
     p_exp = (gamma_const - 1.e0_rt)*exp_energy/vctr

     j = lo(2)

     do i = lo(1), hi(1)
        xmin = xlo(1) + delta(1)*dble(i-lo(1))

        vol_pert    = 0.e0_rt
        vol_ambient = 0.e0_rt

        do ii = 0, nsub-1
           xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5e0_rt)

           dist = xx

           if (dist <= r_init) then
              vol_pert    = vol_pert    + dist
           else
              vol_ambient = vol_ambient + dist
           endif

        enddo

        p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/ (vol_pert + vol_ambient)

        eint = p_zone/(gamma_const - 1.e0_rt)

        state(i,j,URHO) = dens_ambient
        state(i,j,UMX:UMZ) = 0.e0_rt

        state(i,j,UEDEN) = eint + &
             0.5e0_rt*(sum(state(i,j,UMX:UMZ)**2)/state(i,j,URHO))

        state(i,j,UEINT) = eint

        state(i,j,UFS) = state(i,j,URHO)

     enddo

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           state(i,j,URHO ) = state(i,lo(2),URHO)
           state(i,j,UMX:UMZ) = 0.e0_rt
           state(i,j,UEDEN) = state(i,lo(2),UEDEN)
           state(i,j,UEINT) = state(i,lo(2),UEINT)
           state(i,j,UFS  ) = state(i,lo(2),UFS)
        end do
     enddo


  ! Spherical problem in cylindrical (axisymmetric) coordinates
  else if (probtype .eq. 23) then

     ! set explosion pressure -- we will convert the point-explosion
     ! energy into a corresponding pressure distributed throughout the
     ! perturbed volume
     vctr = FOUR3RD*M_PI*r_init**3
     p_exp = (gamma_const - 1.e0_rt)*exp_energy/vctr

     do j = lo(2), hi(2)
        ymin = xlo(2) + delta(2)*dble(j-lo(2))

        do i = lo(1), hi(1)
           xmin = xlo(1) + delta(1)*dble(i-lo(1))

           vol_pert    = 0.e0_rt
           vol_ambient = 0.e0_rt

           do jj = 0, nsub-1
              yy = ymin + (delta(2)/dble(nsub))*(dble(jj) + 0.5e0_rt)

              do ii = 0, nsub-1
                 xx = xmin + (delta(1)/dble(nsub))*(dble(ii) + 0.5e0_rt)

                 dist = sqrt(xx**2 + yy**2)

                 ! The volume of a cell is a annular cylindrical region.
                 ! The main thing that matters is the distance from the
                 ! symmetry axis.
                 !   V = pi*dy*(x_r**2 - x_l**2) = pi*dy*dx*HALF*xx
                 ! (where x_r is the coordinate of the x right edge,
                 !        x_l is the coordinate of the x left edge,
                 !    and xx  is the coordinate of the x center of the cell)
                 !
                 ! since dx and dy are constant, they cancel out
                 if (dist <= r_init) then
                    vol_pert    = vol_pert    + xx
                 else
                    vol_ambient = vol_ambient + xx
                 endif

              enddo
           enddo

           p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/ (vol_pert + vol_ambient)

           eint = p_zone/(gamma_const - 1.e0_rt)

           state(i,j,URHO) = dens_ambient
           state(i,j,UMX:UMZ) = 0.e0_rt

           state(i,j,UEDEN) = eint + &
                0.5e0_rt*(sum(state(i,j,UMX:UMZ)**2)/state(i,j,URHO))

           state(i,j,UEINT) = eint

           state(i,j,UFS) = state(i,j,URHO)

        enddo
     enddo

  else
     call bl_abort('Dont know this probtype')
  end if

end subroutine ca_initdata
