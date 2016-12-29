subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use network, only : network_init
  use eos_module
  use network, only : nspec
  use bl_error_module

  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(1), probhi(1)

  type(eos_t) :: eos_state

  integer :: untin,i

  namelist /fortin/ rho_0, T_0

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  call network_init()

  if (namlen .gt. maxlen) then
     call bl_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  rho_0 = 1.0             ! not used -- no hydro
  T_0 = 5.797e5           ! 50 eV

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  eos_state % rho = rho_0
  eos_state % T   = T_0
  eos_state % xn  = 0.d0
  eos_state % xn(1) = 1.d0

  call eos(eos_input_rt, eos_state)

  rhoe_0 = rho_0 * eos_state % e

  !     domain extrema
  xmin = problo(1)
  xmax = probhi(1)

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
                       state,state_l1,state_h1,delta,xlo,xhi)

  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux

  implicit none

  integer level, nscal
  integer lo(1), hi(1)
  integer state_l1,state_h1
  double precision state(state_l1:state_h1,NVAR)
  double precision time, delta(1)
  double precision xlo(1), xhi(1)

  double precision xcen
  integer i

  do i = lo(1), hi(1)
     xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

     state(i,URHO ) = rho_0
     state(i,UMX  ) = 0.d0
     state(i,UEDEN) = rhoe_0 
     state(i,UEINT) = rhoe_0

     ! set the composition to be all in the first species
     state(i,UFS:UFS-1+nspec) = 0.d0
     state(i,UFS) = state(i,URHO)
     if (naux > 0) then
        state(i,UFX) = state(i,URHO)
     end if

     state(i,UTEMP) = T_0

  enddo

end subroutine ca_initdata

! ::: 
! ::: -----------------------------------------------------------
! :::
subroutine ca_initrad(level,time,lo,hi,nrad, &
                      rad_state,rad_state_l1,rad_state_h1, &
                      delta,xlo,xhi)

  use probdata_module
  use rad_params_module, only : dnugroup, nugroup
  use fundamental_constants_module, only : hplanck, c_light, k_B
  use bl_constants_module, only : ZERO, ONE, HALF, M_PI

  implicit none

  integer level, nrad
  integer lo(1), hi(1)
  integer rad_state_l1,rad_state_h1

  real(kind=8) xlo(1), xhi(1), time, delta(1)
  real(kind=8) rad_state(rad_state_l1:rad_state_h1,0:nrad-1)

  ! local variables
  double precision xcen, nu, xx
  integer i, igroup

  !        if (nrad .ne. 64) then
  !           print *, 'ERROR: initial radiation field assumes 64 groups'
  !           stop
  !        endif

  do i = lo(1), hi(1)  
     xcen = xmin + delta(1)*(float(i) + HALF)

     do igroup=0,nrad-1
        nu = nugroup(igroup)

        xx = hplanck*nu/(k_B*T_0)
        if (xx > 708.d0) then
           ! exp(708+eps) will give 1.e308 -- the limits of IEEE floating point, 
           ! and generate an overflow exception
           rad_state(i,igroup) = ZERO
        else
           rad_state(i,igroup) =  &
                (8.d0*M_PI*hplanck*nu**3/c_light**3)/(exp(xx) - ONE) &
                * dnugroup(igroup)
        endif
     end do

     ! Planck function * dnu as generated using analytic.f90
     ! with t_obs = 1.e-20 (then the analytic solution is just the
     ! Planck function).

     ! use cat analytic.out | awk '{print "rad_state(i,", $1, ") = ", $3}'
     ! to get the output from analytic.f90 into this form
     ! rad_state(i, 1 ) =  21.10636364
     ! rad_state(i, 2 ) =  71.33050571
     ! rad_state(i, 3 ) =  132.9955053
     ! rad_state(i, 4 ) =  247.8779431
     ! rad_state(i, 5 ) =  461.7858265
     ! rad_state(i, 6 ) =  859.8027717
     ! rad_state(i, 7 ) =  1599.762187
     ! rad_state(i, 8 ) =  2973.990893
     ! rad_state(i, 9 ) =  5522.858284
     ! rad_state(i, 10 ) =  10242.82341
     ! rad_state(i, 11 ) =  18965.84768
     ! rad_state(i, 12 ) =  35047.25997
     ! rad_state(i, 13 ) =  64603.54109
     ! rad_state(i, 14 ) =  118718.4320
     ! rad_state(i, 15 ) =  217326.4792
     ! rad_state(i, 16 ) =  395939.3330
     ! rad_state(i, 17 ) =  717045.7734
     ! rad_state(i, 18 ) =  1288867.597
     ! rad_state(i, 19 ) =  2294942.906
     ! rad_state(i, 20 ) =  4037920.159
     ! rad_state(i, 21 ) =  6997937.720
     ! rad_state(i, 22 ) =  11895771.81
     ! rad_state(i, 23 ) =  19726334.77
     ! rad_state(i, 24 ) =  31680948.94
     ! rad_state(i, 25 ) =  48809796.05
     ! rad_state(i, 26 ) =  71235973.01
     ! rad_state(i, 27 ) =  96866897.68
     ! rad_state(i, 28 ) =  120107308.8
     ! rad_state(i, 29 ) =  132103959.4
     ! rad_state(i, 30 ) =  124530830.3
     ! rad_state(i, 31 ) =  96479770.59
     ! rad_state(i, 32 ) =  58416883.53
     ! rad_state(i, 33 ) =  26028725.98
     ! rad_state(i, 34 ) =  7937733.816
     ! rad_state(i, 35 ) =  1516781.179
     ! rad_state(i, 36 ) =  162962.0252
     ! rad_state(i, 37 ) =  8615.969045
     ! rad_state(i, 38 ) =  190.3824212
     ! rad_state(i, 39 ) =  1.682720322
     ! rad_state(i, 40 ) =  0.2871907407E-02
     ! rad_state(i, 41 ) =  0.1126867091E-05
     ! rad_state(i, 42 ) =  0.5935611467E-10
     ! rad_state(i, 43 ) =  0.2636914282E-15
     ! rad_state(i, 44 ) =  0.5574322657E-22
     ! rad_state(i, 45 ) =  0.2771077089E-30
     ! rad_state(i, 46 ) =  0.1359913596E-40
     ! rad_state(i, 47 ) =  0.2262448350E-53
     ! rad_state(i, 48 ) =  0.3421433381E-69
     ! rad_state(i, 49 ) =  0.9299258736E-89
     ! rad_state(i, 50 ) =  0.000000000
     ! rad_state(i, 51 ) =  0.000000000
     ! rad_state(i, 52 ) =  0.000000000
     ! rad_state(i, 53 ) =  0.000000000
     ! rad_state(i, 54 ) =  0.000000000
     ! rad_state(i, 55 ) =  0.000000000
     ! rad_state(i, 56 ) =  0.000000000
     ! rad_state(i, 57 ) =  0.000000000
     ! rad_state(i, 58 ) =  0.000000000
     ! rad_state(i, 59 ) =  0.000000000
     ! rad_state(i, 60 ) =  0.000000000
     ! rad_state(i, 61 ) =  0.000000000
     ! rad_state(i, 62 ) =  0.000000000
     ! rad_state(i, 63 ) =  0.000000000
     ! rad_state(i, 64 ) =  0.000000000

  enddo

end subroutine ca_initrad

