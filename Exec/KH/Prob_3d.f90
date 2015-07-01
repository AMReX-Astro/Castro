subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_error_module
  use probdata_module
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ rho1, rho2, u1, u2, pres, L, vfac, vmode, idir

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  rho1 = 1.0d0
  rho2 = 2.0d0
  u1 = 0.5d0
  u2 = -0.5d0
  pres = 2.5d0

  L = 0.025d0

  vfac = 0.01d0
  vmode = 4.0d0

  idir     = 1

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

end subroutine PROBINIT


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
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, &
                                 UEINT, UFS, UTEMP, small_temp

  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1, &
                         state_l2:state_h2, &
                         state_l3:state_h3,NVAR)

  double precision xcen,ycen,zcen
  double precision rhom, um
  double precision e,rho,u
  double precision pi
  integer i,j,k

  type (eos_t) :: eos_state

  rhom = 0.5d0 * (rho1 - rho2) ! McNally+ Eq. 2
  um = 0.5d0 * (u1 - u2)       ! McNally+ Eq. 4

  pi = 2.0d0*asin(1.0d0)

  state = 0.0d0

  ! single species 
  state(:,:,:,UFS) = 1.0d0

  do k = lo(3), hi(3)
     zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
     do j = lo(2), hi(2)
        ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)
        do i = lo(1), hi(1)
           xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)
           
           if (idir == 1) then
              ! flow in x-direction; use z-direction as shear; perturb vy,vz
              
              ! this assumes 0 <= z <= 1
              ! McNally Eqns. 1 and 3
              if (zcen < 0.25d0) then
                 rho = smooth(rho1,-rhom,-0.25d0,L,zcen)
                 u   = smooth(u1  ,-um  ,-0.25d0,L,zcen)
              else if (zcen < 0.5d0) then
                 rho = smooth(rho2,rhom,0.25d0,L,-zcen)
                 u   = smooth(u2  ,um  ,0.25d0,L,-zcen)
              else if (zcen < 0.75d0) then
                 rho = smooth(rho2,rhom,-0.75d0,L,zcen)
                 u   = smooth(u2  ,um  ,-0.75d0,L,zcen)
              else
                 rho = smooth(rho1,-rhom,0.75d0,L,-zcen)
                 u   = smooth(u1  ,-um  ,0.75d0,L,-zcen)
              end if
              
              ! momentum field
              state(i,j,k,UMX) = u * rho
              state(i,j,k,UMY) = rho * vfac*cos(vmode*pi*xcen) ! out of phase
              state(i,j,k,UMZ) = rho * vfac*sin(vmode*pi*xcen) ! McNally+ Eq 5
              
           else if (idir == 2) then
              ! flow in y-direction; use z-direction as shear; perturb vx,vz

              ! this assumes 0 <= z <= 1
              ! McNally Eqns. 1 and 3
              if (zcen < 0.25d0) then
                 rho = smooth(rho1,-rhom,-0.25d0,L,zcen)
                 u   = smooth(u1  ,-um  ,-0.25d0,L,zcen)
              else if (zcen < 0.5d0) then
                 rho = smooth(rho2,rhom,0.25d0,L,-zcen)
                 u   = smooth(u2  ,um  ,0.25d0,L,-zcen)
              else if (zcen < 0.75d0) then
                 rho = smooth(rho2,rhom,-0.75d0,L,zcen)
                 u   = smooth(u2  ,um  ,-0.75d0,L,zcen)
              else
                 rho = smooth(rho1,-rhom,0.75d0,L,-zcen)
                 u   = smooth(u1  ,-um  ,0.75d0,L,-zcen)
              end if
              
              ! momentum field
              state(i,j,k,UMX) = rho * vfac*cos(vmode*pi*xcen) ! out of phase
              state(i,j,k,UMY) = u * rho
              state(i,j,k,UMZ) = rho * vfac*sin(vmode*pi*xcen) ! McNally+ Eq 5
              
           else if (idir == 3) then
              ! flow in z-direction; use y-direction as shear; perturb vx,vy
              
              ! this assumes 0 <= y <= 1
              ! McNally Eqns. 1 and 3
              if (ycen < 0.25d0) then
                 rho = smooth(rho1,-rhom,-0.25d0,L,ycen)
                 u   = smooth(u1  ,-um  ,-0.25d0,L,ycen)
              else if (ycen < 0.5d0) then
                 rho = smooth(rho2,rhom,0.25d0,L,-ycen)
                 u   = smooth(u2  ,um  ,0.25d0,L,-ycen)
              else if (ycen < 0.75d0) then
                 rho = smooth(rho2,rhom,-0.75d0,L,ycen)
                 u   = smooth(u2  ,um  ,-0.75d0,L,ycen)
              else
                 rho = smooth(rho1,-rhom,0.75d0,L,-ycen)
                 u   = smooth(u1  ,-um  ,0.75d0,L,-ycen)
              end if
              
              ! momentum field
              state(i,j,k,UMX) = rho * vfac*cos(vmode*pi*xcen) ! out of phase
              state(i,j,k,UMY) = rho * vfac*sin(vmode*pi*xcen) ! McNally+ Eq 5
              state(i,j,k,UMZ) = u * rho
              
           else

              call bl_abort("invalid idir")
           end if

           ! The rest of the variables
           state(i,j,k,URHO) = rho
           
           eos_state % rho = rho
           eos_state % p   = pres
           eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1)
           eos_state % T   = small_temp ! Initial guess for the EOS
           
           ! Get the temperature and internal energy assuming fixed pressure
           call eos(eos_input_rp, eos_state)
           state(i,j,k,UEINT) = eos_state % e * rho
           state(i,j,k,UTEMP) = eos_state % T
           
           !  Total energy
           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                0.5d0 * ( state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + &
                          state(i,j,k,UMZ)**2 ) / state(i,j,k,URHO)

           ! Convert mass fractions to conserved quantity
           state(i,j,k,UFS:) = state(i,j,k,UFS:) * rho
           
        enddo
     enddo
  enddo
  
contains
  function smooth(a,b,offset,L,pos) result(r)
    ! this is the general form of the smooth function in McNally+ Eq 1 or 3
    double precision :: a, b, offset, L, pos
    double precision :: r

    r = a + b*exp((pos+offset)/L)
  end function smooth

end subroutine ca_initdata

! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only : NVAR

  implicit none
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

  integer n

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
          adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                      domlo,domhi,delta,xlo,time,bc)

  implicit none
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill
