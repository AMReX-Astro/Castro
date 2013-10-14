subroutine PROBINIT (init,name,namlen,problo,probhi)

  use eos_module
  use eos_type_module
  use bl_error_module 
  use network
  use probdata_module

  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(3), probhi(3)
  double precision xn(nspec)

  integer untin,i

  type (eos_t) :: eos_state

  namelist /fortin/ p_l, u_l, rho_l, p_r, u_r, rho_r, T_l, T_r, frac, idir, &
       use_Tinit, &
       denerr,  dengrad,  max_denerr_lev,  max_dengrad_lev, &
       velgrad,  max_velgrad_lev, &
       presserr,pressgrad,max_presserr_lev,max_pressgrad_lev, xcloud, cldradius

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  p_l = 1.0               ! left pressure (erg/cc)
  u_l = 0.0               ! left velocity (cm/s)
  rho_l = 1.0             ! left density (g/cc)
  T_l = 1.0

  p_r = 0.1               ! right pressure (erg/cc)
  u_r = 0.0               ! right velocity (cm/s)
  rho_r = 0.125           ! right density (g/cc)
  T_r = 1.0

  idir = 1                ! direction across which to jump
  frac = 0.5              ! fraction of the domain for the interface

  use_Tinit = .false.     ! optionally use T_l/r instead of p_l/r for initialization

  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1

  velgrad = 1.d20
  max_velgrad_lev = -1

  xcloud = -1.0
  cldradius = -1.0

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = frac*(problo(1)+probhi(1))
  center(2) = frac*(problo(2)+probhi(2))
  center(3) = frac*(problo(3)+probhi(3))
  
  ! compute the internal energy (erg/cc) for the left and right state
  xn(:) = 0.0d0
  xn(1) = 1.0d0

  if (use_Tinit) then

     eos_state%rho = rho_l
     eos_state%T = T_l
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state, .false.)
 
     rhoe_l = rho_l*eos_state%e
     p_l = eos_state%p

     eos_state%rho = rho_r
     eos_state%T = T_r
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state, .false.)
 
     rhoe_r = rho_r*eos_state%e
     p_r = eos_state%p

  else

     eos_state%rho = rho_l
     eos_state%p = p_l
     eos_state%T = 100000.d0  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state, .false.)
 
     rhoe_l = rho_l*eos_state%e
     T_l = eos_state%T

     eos_state%rho = rho_r
     eos_state%p = p_r
     eos_state%T = 100000.d0  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state, .false.)
 
     rhoe_r = rho_r*eos_state%e
     T_r = eos_state%T

  endif

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

  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision xlo(3), xhi(3), time, delta(3)
  double precision state(state_l1:state_h1,state_l2:state_h2, &
                         state_l3:state_h3,NVAR)

  double precision xcen,ycen,zcen
  integer i,j,k

  double precision rsq, dist, ycloud, zcloud, ecent, denfact, rhocld

  ycloud = 0.0
  zcloud = 0.0
  !  ecent  = 1.96
  ecent  = 1.00
  denfact = 3.0246
  rhocld  = rho_r * denfact

  do k = lo(3), hi(3)
     zcen = xlo(3) + delta(3)*(float(k-lo(3)) + 0.5d0)
     
     do j = lo(2), hi(2)
        ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5d0)

        do i = lo(1), hi(1)
           xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5d0)

           if (idir == 1) then
              if (xcen <= center(1)) then      !! left of shock
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = rho_l*u_l
                 state(i,j,k,UMY) = 0.d0
                 state(i,j,k,UMZ) = 0.d0
                 state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l*u_l
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else

               rsq = (((xcen-xcloud)**2) + ((ycen-ycloud)**2) + ((zcen-zcloud)**2))/ecent
               !  rsq = ((xcen-xcloud)**2) + ((ycen-ycloud)**2) + ((zcen-zcloud)**2)/ecent
               dist = sqrt(rsq)
               if (dist .le. cldradius) then   !! inside bubble
                 state(i,j,k,URHO) = rhocld
                 state(i,j,k,UMX) = 0.d0
                 state(i,j,k,UMY) = 0.d0
                 state(i,j,k,UMZ) = 0.d0
                 state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
               else                            !! right of shock
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = rho_r*u_r
                 state(i,j,k,UMY) = 0.d0
                 state(i,j,k,UMZ) = 0.d0
                 state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r*u_r
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
               endif
              endif
              
           else
              call bl_abort('invalid idir')
           endif
 
           state(i,j,k,UFS:UFS-1+nspec) = 0.0d0
           state(i,j,k,UFS  ) = state(i,j,k,URHO)

           
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
