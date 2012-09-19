subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use eos_module
  use network, only : nspec
  implicit none 

  integer :: init,namlen,untin,i,k
  integer :: name(namlen)

  double precision :: center_x, center_y, center_z
  double precision :: problo(1), probhi(1)

  double precision :: eint

  namelist /fortin/ &
       rho_0, r_0, p_0, rho_ambient, smooth_delta, &
       denerr, dengrad, max_denerr_lev, max_dengrad_lev, &
       velerr, velgrad, max_velerr_lev, max_velgrad_lev, &
       presserr, pressgrad, max_presserr_lev, max_pressgrad_lev, &
       center_x, center_y, center_z

!
!     Build "probin" filename -- the name of file containing fortin namelist.
!     
  integer, parameter :: maxlen = 127
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
! set namelist defaults

  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1
  
  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1
  
  velgrad = 1.d20
  max_velgrad_lev = -1
  
  rho_0 = 1.d9
  r_0 = 6.5d8
  p_0 = 1.d10
  rho_ambient = 1.d0
  smooth_delta = 1.d-5

!     Read namelists in probin file
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! in 1-d spherical, the lower domain boundary should be the origin
  center(1) = 0.0d0

  xmin = problo(1)
  if (xmin /= 0.d0) then
     print *, 'ERROR: xmin should be 0!'
     stop
  endif

  xmax = probhi(1)

  ! set the composition to be uniform
  allocate(X_0(nspec))
  
  X_0(:) = 0.0
  X_0(1) = 1.0

  ! get the ambient temperature and sphere temperature, T_0
  call eos_e_given_RPX(eint, T_0,       rho_0,       p_0, X_0)
  call eos_e_given_RPX(eint, T_ambient, rho_ambient, p_0, X_0)

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
                       state,state_l1,state_h1,delta,xlo,xhi)

  use probdata_module
  use eos_module
  use network, only : nspec
  use interpolate_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS
  
  implicit none
  
  integer          :: level, nscal
  integer          :: lo(1), hi(1)
  integer          :: state_l1,state_h1
  double precision :: state(state_l1:state_h1,NVAR)
  double precision :: time, delta(1)
  double precision :: xlo(1), xhi(1)
  
  double precision :: xl,xx,dist,pres,eint,temp,avg_rho, rho_n
  double precision :: dx_sub
  integer          :: i,ii,n

  integer, parameter :: nsub = 5

  dx_sub = delta(1)/dble(nsub)
  
  do i = lo(1), hi(1)
     xl = dble(i) * delta(1)
      
     avg_rho = 0.d0

     do ii = 0, nsub-1

        xx = xl + (dble(ii) + 0.5d0) * dx_sub
        
        dist = (xx-center(1))
        
        ! use a tanh profile to smooth the transition between rho_0                                    
        ! and rho_ambient                                                                              
        rho_n = rho_0 - (rho_0 - rho_ambient)* &
             0.5d0 * (1.d0 + tanh((dist - r_0)/smooth_delta))
                
        avg_rho = avg_rho + rho_n
        
     enddo

     state(i,URHO) = avg_rho/dble(nsub)
     
     call eos_e_given_RPX(eint, temp, state(i,URHO), p_0, X_0)
     
     state(i,UTEMP) = temp
     state(i,UMX) = 0.d0
     state(i,UEDEN) = state(i,URHO) * eint
     state(i,UEINT) = state(i,URHO) * eint
     state(i,UFS:UFS+nspec-1) = state(i,URHO) * X_0(1:nspec)

  enddo
    
end subroutine ca_initdata


! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                      domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UEINT, UFS
  implicit none

  include 'bc_types.fi'

  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1,NVAR)

  integer n, i
  double precision :: vel

  ! call the generic ghostcell filling routine
  do n = 1,NVAR
     call filcc(adv(adv_l1,n), adv_l1,adv_h1, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill


! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_h1, &
                            domlo,domhi,delta,xlo,time,bc)
  implicit none
    
  integer          :: adv_l1,adv_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  logical          :: rho_only
  
  call filcc(adv,adv_l1,adv_h1, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill


! ::: -----------------------------------------------------------

subroutine ca_gravxfill(grav,grav_l1,grav_h1, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  
  integer          :: grav_l1,grav_h1
  integer          :: bc(1,2,*)
  integer          :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: grav(grav_l1:grav_h1)

  call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill

