
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use network, only : network_init
  use eos_module, only : gamma_const
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)
  
  integer untin,i
  
  namelist /fortin/ rho0, T0, v0, rho1, T1, v1, &
       denerr,   dengrad,  max_denerr_lev,  max_dengrad_lev, &
       velerr,   velgrad,  max_velerr_lev,  max_pressgrad_lev, &
       presserr, pressgrad,max_presserr_lev,max_pressgrad_lev, &
       temperr,  tempgrad, max_temperr_lev, max_tempgrad_lev, &
       raderr,   radgrad,  max_raderr_lev,  max_radgrad_lev
  
  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)
  
  call network_init()
  
  if (namlen .gt. maxlen) then
     write(6,*) 'probin file name too long'
     stop
  end if
  
  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
  
  ! set namelist defaults
  
  rho0 = 5.4588672836d-13
  T0 = 100.d0
  v0 = 235435.371882d0
  rho1 = 1.24793794736d-12
  T1 = 207.756999533d0
  v1 = 102986.727159d0
  
  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1
  
  velerr = 1.d20
  velgrad = 1.d20
  max_velerr_lev = -1
  max_velgrad_lev = -1
  
  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1
  
  temperr = 1.d20
  tempgrad = 1.d20
  max_temperr_lev = -1
  max_tempgrad_lev = -1
  
  raderr = 1.d20
  radgrad = 1.d20
  max_raderr_lev = -1
  max_radgrad_lev = -1
  
  center(1) = (problo(1)+probhi(1))/2.d0
  xmin = problo(1)
  xmax = probhi(1)

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
     state,state_l1,state_h1, &
     delta,xlo,xhi)
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use eos_module
  
  implicit none
  
  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  double precision :: state(state_l1:state_h1,NVAR)
  double precision :: time, delta(1)
  double precision :: xlo(1), xhi(1)
  
  integer :: i
  double precision :: X(nspec)
  double precision :: p, eint, xcell

  do i = lo(1), hi(1)
  
     xcell = xlo(1) + delta(1) * (float(i-lo(1)) + 0.5d0)

     if (xcell < 0.d0) then
        state(i,URHO) = rho0
        
        ! set the composition to be all in the first species
        state(i,UFS:UFS-1+nspec) = 0.d0
        state(i,UFS  ) = state(i,URHO)

        state(i,UTEMP) = T0
        state(i,UMX) = rho0*v0
     else
        state(i,URHO) = rho1
        
        ! set the composition to be all in the first species
        state(i,UFS:UFS-1+nspec) = 0.d0
        state(i,UFS  ) = state(i,URHO)

        state(i,UTEMP) = T1
        state(i,UMX) = rho1*v1
     end if
 
     if (naux > 0) then
        state(i,UFX) = state(i,URHO)        
     end if

     ! set the internal energy via the EOS
     X(:) = state(i,UFS:UFS-1+nspec)/state(i,URHO)
     call eos_given_RTX(eint, p, state(i,URHO), state(i,UTEMP), X)
     
     state(i,UEINT) = state(i,URHO)*eint
     state(i,UEDEN) = state(i,URHO)*eint + &
          0.5*(state(i,UMX)**2)/state(i,URHO)                
  enddo
  
end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! :::
subroutine ca_initrad(level,time,lo,hi,nrad, &
     rad_state,rad_state_l1, &
     rad_state_h1, &
     delta,xlo,xhi)

  use probdata_module
  use fundamental_constants_module, only: a_rad
  use rad_params_module, only : xnu
  use blackbody_module, only : BGroup
  
  implicit none
  integer :: level, nrad
  integer :: lo(1), hi(1)
  integer :: rad_state_l1,rad_state_h1
  double precision :: xlo(1), xhi(1), time, delta(1)
  double precision ::  rad_state(rad_state_l1:rad_state_h1, 0:nrad-1)

  ! local variables
  integer :: i, igroup
  double precision xcell, t

  do i = lo(1), hi(1)

     xcell = xlo(1) + delta(1) * (float(i-lo(1)) + 0.5d0)
   
     if (xcell < 0.d0) then
        T = T0
     else
        T = T1
     end if
     
     if (nrad == 1) then
        rad_state(i,:) = a_rad*T**4
     else
        do igroup=0,nrad-1
           rad_state(i,igroup) = BGroup(T, xnu(igroup), xnu(igroup+1))
        end do
     end if

  enddo

end subroutine ca_initrad


! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_h1, &
     domlo,domhi,delta,xlo,time,bc)
  
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UFS, UFX, UTEMP
  use network, only : nspec, naux
  use eos_module
  use probdata_module
  
  implicit none
  
  include 'bc_types.fi'
  integer adv_l1,adv_h1
  integer bc(1,2,*)
  integer domlo(1), domhi(1)
  double precision delta(1), xlo(1), time
  double precision adv(adv_l1:adv_h1,NVAR)
  
  integer i, n
  
  double precision :: X(nspec), eint, p
  double precision, save :: eint0, etot0, eint1, etot1
  logical, save :: first_call = .true.
  
  if (first_call) then
     first_call = .false.
     
     x(:) = 1.d0
     call eos_given_RTX(eint, p, rho0, T0, X)
     eint0 = rho0 * eint
     etot0 = eint0 + 0.5*rho0*v0**2
     
     call eos_given_RTX(eint, p, rho1, T1, X)
     eint1 = rho1 * eint
     etot1 = eint1 + 0.5*rho1*v1**2         
  end if
  
  do n = 1,NVAR
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,bc(1,1,n))
  enddo
  
  !     The strategy here is to set Dirichlet condition for inflow and
  !     outflow boundaries, and let the Riemann solver sort out the
  !     proper upwinding.  However, this decision makes this routine
  !     look somewhat non-orthodox, in that we need to set external
  !     values in either case....how do we know it's Outflow?  We have
  !     to assume that the setup routines converted Outflow to FOEXTRAP.
  
  !     XLO
  if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        adv(i,URHO) = rho0
        adv(i,UMX) = rho0*v0
        adv(i,UFS) = adv(i,URHO)
        if (naux>0) then
           adv(i,UFX) = adv(i,URHO)           
        end if
        adv(i,UEINT) = eint0
        adv(i,UEDEN) = etot0
        adv(i,UTEMP) = T0
     end do
  end if
  
  !     XHI
  if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        adv(i,URHO) = rho1
        adv(i,UMX) = rho1*v1
        adv(i,UFS) = adv(i,URHO)
        if (naux>0) then
           adv(i,UFX) = adv(i,URHO)           
        end if
        adv(i,UEINT) = eint1
        adv(i,UEDEN) = etot1
        adv(i,UTEMP) = T1
     end do
  end if
  
end subroutine ca_hypfill

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)

      use probdata_module

      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)
      integer i

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

!     XLO
      if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         do i = adv_l1, domlo(1)-1
            adv(i) = rho0
         end do
      end if            

!     XHI
      if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
         do i = domhi(1)+1, adv_h1
            adv(i) = rho1
         end do
      end if            

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_radfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)

      use probdata_module

      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)
      integer i

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

!     XLO
      if ( bc(1,1,1).eq.EXT_DIR .and. adv_l1.lt.domlo(1)) then
         do i = adv_l1, domlo(1)-1
            adv(i) = adv(domlo(1))
         end do
      end if            

!     XHI
      if ( bc(1,2,1).eq.EXT_DIR .and. adv_h1.gt.domhi(1)) then
         do i = domhi(1)+1, adv_h1
            adv(i) = adv(domhi(1))
         end do
      end if            

      end

