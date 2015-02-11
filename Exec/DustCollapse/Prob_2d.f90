subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use eos_module
  use meth_params_module, only: small_temp
  use prob_params_module, only : center
  use network, only : nspec
  implicit none 

  integer :: init,namlen,untin,i,k
  integer :: name(namlen)

  double precision :: center_x, center_y, center_z
  double precision :: problo(2), probhi(2)

  type (eos_t) :: eos_state

  namelist /fortin/ &
       rho_0, r_0, r_old, p_0, rho_ambient, smooth_delta, &
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

  rho_0 = 1.d9
  r_0 = 6.5d8
  r_old = r_0
  p_0 = 1.d10
  rho_ambient = 1.d0
  smooth_delta = 1.d-5

!     Read namelists in probin file
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! we are going to enforce that the lower left corner of the domain is 0,0
  ! (i.e., we only model a quadrant)

  center(1) = center_x
  center(2) = center_y

  xmin = problo(1)
  if (xmin /= 0.d0) then
     print *, 'ERROR: xmin should be 0!'
     stop
  endif

  xmax = probhi(1)

  ymin = problo(2)
  if (ymin /= 0.d0) then
     print *, 'ERROR: ymin should be 0!'
     stop
  endif

  ymax = probhi(2)

  ! set the composition to be uniform
  allocate(X_0(nspec))
  
  X_0(:) = 0.0
  X_0(1) = 1.0

  ! get the ambient temperature and sphere temperature, T_0

  eos_state % rho = rho_0
  eos_state % p   = p_0
  eos_state % xn  = x_0
  eos_state % T   = small_temp ! Initial guess for the EOS

  call eos(eos_input_rp, eos_state)

  T_0 = eos_state % T
  
  eos_state % rho = rho_ambient
  
  call eos(eos_input_rp, eos_state)

  T_ambient = eos_state % T

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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module
  use network, only : nspec
  use interpolate_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, small_temp
  use prob_params_module, only : center

  implicit none

  integer          :: level, nscal
  integer          :: lo(2), hi(2)
  integer          :: state_l1,state_l2,state_h1,state_h2
  double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)
  double precision :: time, delta(2)
  double precision :: xlo(2), xhi(2)
  
  double precision :: xl,yl,xx,yy,dist,pres,eint,temp,avg_rho, rho_n
  double precision :: dx_sub,dy_sub
  integer          :: i,j,ii,jj,n

  type (eos_t) :: eos_state

  integer, parameter :: nsub = 5

  dx_sub = delta(1)/dble(nsub)
  dy_sub = delta(2)/dble(nsub)

  do j = lo(2), hi(2)
     yl = ymin + dble(j) * delta(2)

     do i = lo(1), hi(1)
        xl = xmin + dble(i) * delta(1)

        avg_rho = 0.d0

        do jj = 0, nsub-1
           yy = yl + (dble(jj) + 0.5d0) * dy_sub

           do ii = 0, nsub-1
              xx = xl + (dble(ii) + 0.5d0) * dx_sub

              dist = sqrt((xx-center(1))**2 + (yy-center(2))**2)

              ! use a tanh profile to smooth the transition between rho_0 
              ! and rho_ambient
              rho_n = rho_0 - 0.5d0*(rho_0 - rho_ambient)* &
                   (1.d0 + tanh((dist - r_0)/smooth_delta))

              avg_rho = avg_rho + rho_n

           enddo
        enddo
        
        state(i,j,URHO) = avg_rho/(nsub*nsub)
 
        eos_state % rho = state(i,j,URHO)
        eos_state % p   = p_0
        eos_state % T   = small_temp ! Initial guess for the EOS
        eos_state % xn  = X_0

        call eos(eos_input_rp, eos_state)

        temp = eos_state % T
        eint = eos_state % e

        state(i,j,UTEMP) = temp
        state(i,j,UMX) = 0.d0
        state(i,j,UMY) = 0.d0
        state(i,j,UEDEN) = state(i,j,URHO) * eint
        state(i,j,UEINT) = state(i,j,URHO) * eint
        state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO) * X_0(1:nspec)

     enddo
  enddo

end subroutine ca_initdata


! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : NVAR
  implicit none

  include 'bc_types.fi'

  integer          :: adv_l1,adv_l2,adv_h1,adv_h2
  integer          :: bc(2,2,*)
  integer          :: domlo(2), domhi(2)
  double precision :: delta(2), xlo(2), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

  integer n, i
  double precision :: vel

  ! call the generic ghostcell filling routine
  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n), &
          adv_l1,adv_l2,adv_h1,adv_h2, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)
  implicit none
    
  integer          :: adv_l1,adv_l2,adv_h1,adv_h2
  integer          :: bc(2,2,*)
  integer          :: domlo(2), domhi(2)
  double precision :: delta(2), xlo(2), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2)
  logical          :: rho_only
  
  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
             domlo,domhi,delta,xlo,bc)

end subroutine ca_denfill

! ::: -----------------------------------------------------------

subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                        domlo,domhi,delta,xlo,time,bc)

  use probdata_module
  implicit none
  
  integer          :: grav_l1,grav_l2,grav_h1,grav_h2
  integer          :: bc(2,2,*)
  integer          :: domlo(2), domhi(2)
  double precision :: delta(2), xlo(2), time
  double precision :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

  call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

end subroutine ca_gravxfill

! ::: -----------------------------------------------------------

      subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                              domlo,domhi,delta,xlo,time,bc)

      use probdata_module
      implicit none
      include 'bc_types.fi'

      integer :: grav_l1,grav_l2,grav_h1,grav_h2
      integer :: bc(2,2,*)
      integer :: domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

      call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

      end subroutine ca_gravyfill

