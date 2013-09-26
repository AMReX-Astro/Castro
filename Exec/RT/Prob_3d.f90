subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use eos_module, only : gamma_const
  use bl_error_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ denerr,dengrad,max_denerr_lev,max_dengrad_lev, &
       presserr,pressgrad,max_presserr_lev,max_pressgrad_lev,frac, &
       rho_1, rho_2, p0_base

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error('probin file name too long')
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults here
  frac = 0.5d0
  rho_1 = 1.0d0
  rho_2 = 2.0d0
  p0_base = 5.0d0

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! set local variable defaults
  center(1) = frac*(problo(1)+probhi(1))
  center(2) = frac*(problo(2)+probhi(2))
  center(3) = frac*(problo(3)+probhi(3))
  
  L_x = probhi(1) - problo(1)

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
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
       UEDEN, UEINT, UFS, UTEMP
  use bl_constants_module, only: ZERO, HALF, M_PI
  use eos_module, only : gamma_const
  
  implicit none
        
  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  
  integer :: i,j,k
  double precision :: x,y,z,r2d,pres,presmid,pertheight
  
  presmid  = p0_base - rho_1*center(3)
        
  state(:,:,:,UMX)   = ZERO
  state(:,:,:,UMY)   = ZERO
  state(:,:,:,UMZ)   = ZERO
  state(:,:,:,UTEMP) = ZERO

  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           
           if (z .lt. center(3)) then
              pres = p0_base - rho_1*z
              state(i,j,k,UEDEN) = pres / (gamma_const - 1.0d0)
              state(i,j,k,UEINT) = pres / (gamma_const - 1.0d0)
           else
              pres = presmid - rho_2*(z-center(3))
              state(i,j,k,UEDEN) = pres / (gamma_const - 1.0d0)
              state(i,j,k,UEINT) = pres / (gamma_const - 1.0d0)
           end if
           
        enddo
     enddo
  enddo

  do k = lo(3), hi(3)
     z = (k+HALF)*delta(3)
        
     do j = lo(2), hi(2)
        y = (j+HALF)*delta(2)
        
        do i = lo(1), hi(1)
           x = (i+HALF)*delta(1)
     
           r2d = min(sqrt((x-center(1))**2+(y-center(2))**2), 0.5d0*L_x)
           pertheight = 0.5d0 - 0.01d0*cos(2.0d0*M_PI*r2d/L_x)
           state(i,j,k,URHO) = rho_1 + ((rho_2-rho_1)/2.0d0)* &
                (1+tanh((z-pertheight)/0.005d0))
           state(i,j,k,UFS) = state(i,j,k,URHO)
           
        enddo
     enddo
  enddo

end subroutine ca_initdata


! ::: -----------------------------------------------------------
subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                      adv_h3,domlo,domhi,delta,xlo,time,bc)

  use bl_error_module
  use meth_params_module, only : NVAR
  implicit none

  include 'bc_types.fi'
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

  do n = 1,NVAR

     !        XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        call bl_error('SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) ')
     end if

     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        call bl_error('SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) ')
     end if

     !        YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        call bl_error('SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) ')
     end if

     !        YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        call bl_error('SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) ')
     end if

     !        ZLO
     if ( bc(3,1,n).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
        call bl_error('SHOULD NEVER GET HERE bc(3,1,n) .eq. EXT_DIR) ')
     end if

     !        ZHI
     if ( bc(3,2,n).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
        call bl_error('SHOULD NEVER GET HERE bc(3,2,n) .eq. EXT_DIR) ')
     end if

  end do

end subroutine ca_hypfill


subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                      adv_h3,domlo,domhi,delta,xlo,time,bc)

  use bl_error_module
  implicit none

  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer bc(3,2,*)
  integer domlo(3), domhi(3)
  double precision delta(3), xlo(3), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

  !     Note: this function should not be needed, technically, but is provided
  !     to filpatch because there are many times in the algorithm when just
  !     the density is needed.  We try to rig up the filling so that the same
  !     function is called here and in hypfill where all the states are filled.
  
  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
             domlo,domhi,delta,xlo,bc)

  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error('SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) ')
  end if

  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error('SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) ')
  end if

  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error('SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) ')
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error('SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) ')
  end if

  !     ZLO
  if ( bc(3,1,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     call bl_error('SHOULD NEVER GET HERE bc(3,1,1) .eq. EXT_DIR) ')
  endif

  !     ZHI
  if ( bc(3,2,1).eq.EXT_DIR .and. adv_l3.lt.domlo(3)) then
     call bl_error('SHOULD NEVER GET HERE bc(3,2,1) .eq. EXT_DIR) ')
  end if

end subroutine ca_denfill


! ::: -----------------------------------------------------------
subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer :: bc(3,2,*)
  integer :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
  
  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
       domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill


! ::: -----------------------------------------------------------
subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer :: bc(3,2,*)
  integer :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
  
  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
       domlo,domhi,delta,xlo,bc)

end subroutine ca_gravyfill


! ::: -----------------------------------------------------------
subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  include 'bc_types.fi'

  integer :: grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3
  integer :: bc(3,2,*)
  integer :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)
  
  call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3, &
       domlo,domhi,delta,xlo,bc)

end subroutine ca_gravzfill

