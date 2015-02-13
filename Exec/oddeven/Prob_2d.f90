subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use bl_error_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ p_ambient, dens_ambient, dens_pert_factor, vel_pert, &
                    denerr, dengrad, max_denerr_lev, max_dengrad_lev, &
                    presserr, pressgrad, max_presserr_lev, max_pressgrad_lev


  ! Build "probin" filename -- the name of file containing fortin namelist.

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! Set namelist defaults

  p_ambient = 1.d0         ! ambient pressure (in erg/cc)
  dens_ambient = 1.d0      ! ambient density (in g/cc)
  dens_pert_factor = 1.01  ! density enhancement for central zone
  vel_pert = 20.d0

  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1

  ! set center, domain extrema
  center(1) = (problo(1)+probhi(1))/2.d0
  center(2) = (problo(2)+probhi(2))/2.d0

  xmin = problo(1)
  xmax = probhi(1)

  ymin = problo(2)
  ymax = probhi(2)

  
  ! Read namelists
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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module
  use network, only: nspec
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UTEMP
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  double precision xlo(2), xhi(2), time, delta(2)
  double precision state(state_l1:state_h1,state_l2:state_h2,NVAR)

  double precision xcen, ycen
  double precision dens, eint, xvel, X(nspec), temp
  
  integer i,j, icen, jcen

  type (eos_t) :: eos_state

  ! compute the integer location of the center of the domain
  icen = center(1)/delta(1)
  jcen = center(2)/delta(2)

    
  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0)

     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0)

        if (i == icen .and. j == jcen) then
           dens = dens_ambient*dens_pert_factor
        else
           dens = dens_ambient
        endif

        state(i,j,URHO) = dens

        ! velocity perturbation
        if (xcen < center(1)) then 
           xvel = vel_pert
        else if (xcen > center(1)) then
           xvel = -vel_pert
        else
           xvel = 0.d0
        endif

        state(i,j,UMX) = dens*xvel
        state(i,j,UMY) = 0.d0


        ! set the composition
        X(:) = 0.d0
        X(1) = 1.d0
        
        
        ! compute the internal energy and temperature
        eos_state%T = 1.d0 ! initial guess
        eos_state%rho = dens
        eos_state%p = p_ambient
        eos_state%xn(:) = X

        call eos(eos_input_rp, eos_state)

        temp = eos_state%T
        eint = eos_state%e

        state(i,j,UEDEN) = dens*eint +  &
             0.5d0*(state(i,j,UMX)**2/state(i,j,URHO) + &
                    state(i,j,UMY)**2/state(i,j,URHO))

        state(i,j,UEINT) = dens*eint
        state(i,j,UTEMP) = temp

        state(i,j,UFS:UFS-1+nspec) = dens*X(:)

     enddo
  enddo

end subroutine ca_initdata


! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)
 
  use meth_params_module, only : NVAR
  use bl_error_module

  implicit none

  include 'bc_types.fi'
  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

  integer n

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,n), &
                adv_l1,adv_l2,adv_h1,adv_h2, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  do n = 1,NVAR

     !        XLO
     if ( bc(1,1,n).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
        call bl_error("SHOULD NEVER GET HERE bc(1,1,n) .eq. EXT_DIR) ")
     end if

     !        XHI
     if ( bc(1,2,n).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
        call bl_error("SHOULD NEVER GET HERE bc(1,2,n) .eq. EXT_DIR) ")
     end if

     !        YLO
     if ( bc(2,1,n).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
        call bl_error("SHOULD NEVER GET HERE bc(2,1,n) .eq. EXT_DIR) ")
     end if

     !        YHI
     if ( bc(2,2,n).eq.EXT_DIR .and. adv_l2.gt.domhi(2)) then
        call bl_error("SHOULD NEVER GET HERE bc(2,2,n) .eq. EXT_DIR) ")
     end if

  end do

end subroutine ca_hypfill

! ::: -----------------------------------------------------------

subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                      domlo,domhi,delta,xlo,time,bc)

  use bl_error_module

  implicit none
  include 'bc_types.fi'

  integer adv_l1,adv_l2,adv_h1,adv_h2
  integer bc(2,2,*)
  integer domlo(2), domhi(2)
  double precision delta(2), xlo(2), time
  double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

  ! Note: this function should not be needed, technically, but is
  ! provided to filpatch because there are many times in the algorithm
  ! when just the density is needed.  We try to rig up the filling so
  ! that the same function is called here and in hypfill where all the
  ! states are filled.

  call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     call bl_error("SHOULD NEVER GET HERE bc(1,1,1) .eq. EXT_DIR) ")
  end if

  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.lt.domhi(1)) then
     call bl_error("SHOULD NEVER GET HERE bc(1,2,1) .eq. EXT_DIR) ")
  end if

  !     YLO
  if ( bc(2,1,1).eq.EXT_DIR .and. adv_l2.lt.domlo(2)) then
     call bl_error("SHOULD NEVER GET HERE bc(2,1,1) .eq. EXT_DIR) ")
  end if

  !     YHI
  if ( bc(2,2,1).eq.EXT_DIR .and. adv_h2.lt.domhi(2)) then
     call bl_error("SHOULD NEVER GET HERE bc(2,2,1) .eq. EXT_DIR) ")
  end if

end subroutine ca_denfill

