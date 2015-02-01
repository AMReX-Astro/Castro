subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)
  
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
         
  ! set namelist defaults

  p_ambient = 1.d-5        ! ambient pressure (in erg/cc)
  dens_ambient = 1.d0      ! ambient density (in g/cc)
  exp_energy = 1.d0        ! absolute energy of the explosion (in erg)
  r_init = 0.05d0          ! initial radius of the explosion (in cm)
  nsub = 4

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  center(1) = 0.d0

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
  use eos_module, only: gamma_const
  use bl_constants_module, only: M_PI, FOUR3RD
  use meth_params_module , only: NVAR, URHO, UMX, UEDEN, UEINT, UFS

  implicit none

  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  double precision :: xlo(1), xhi(1), time, delta(1)
  double precision :: state(state_l1:state_h1,NVAR)
  
  double precision :: xmin
  double precision :: xx, xl, xr
  double precision :: dx_sub,dist
  double precision :: eint, p_zone
  double precision :: vctr, p_exp

  integer :: i,ii
  double precision :: vol_pert, vol_ambient

  dx_sub = delta(1)/dble(nsub)

  ! Cylindrical coordinates
  if (probtype .eq. 11) then

     ! set explosion pressure -- we will convert the point-explosion energy into
     ! a corresponding pressure distributed throughout the perturbed volume
     vctr = M_PI*r_init**2
     p_exp = (gamma_const - 1.d0)*exp_energy/vctr

     do i = lo(1), hi(1)
        xmin = xlo(1) + delta(1)*dble(i-lo(1))
        vol_pert    = 0.d0
        vol_ambient = 0.d0
        do ii = 0, nsub-1
           xx = xmin + (dble(ii) + 0.5d0) * dx_sub
           dist = xx
           if(dist <= r_init) then
              vol_pert = vol_pert + dist
           else
              vol_ambient = vol_ambient + dist
           endif
        enddo

        p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/(vol_pert+vol_ambient)
  
        eint = p_zone/(gamma_const - 1.d0)
   
        state(i,URHO) = dens_ambient
        state(i,UMX) = 0.d0
  
        state(i,UEDEN) = eint + 0.5d0 * state(i,UMX)**2 / state(i,URHO)
        state(i,UEINT) = eint

        state(i,UFS) = state(i,URHO)
  
     end do

     ! Spherical coordinates
  else if (probtype .eq. 12) then

     ! set explosion pressure -- we will convert the point-explosion energy into
     ! a corresponding pressure distributed throughout the perturbed volume
     vctr = FOUR3RD*M_PI*r_init**3
     p_exp = (gamma_const - 1.d0)*exp_energy/vctr
     
     do i = lo(1), hi(1)
        xmin = xlo(1) + delta(1)*dble(i-lo(1))
        vol_pert    = 0.d0
        vol_ambient = 0.d0
        do ii = 0, nsub-1
           xl = xmin + (dble(ii)        ) * dx_sub
           xr = xmin + (dble(ii) + 1.0d0) * dx_sub
           xx = 0.5d0*(xl + xr)
           
           ! the volume of a subzone is (4/3) pi (xr^3 - xl^3).
           ! we can factor this as: (4/3) pi dr (xr^2 + xl*xr + xl^2)
           ! The (4/3) pi dr factor is common, so we can neglect it. 
           if(xx <= r_init) then
              vol_pert = vol_pert + (xr*xr + xl*xr + xl*xl)
           else
              vol_ambient = vol_ambient + (xr*xr + xl*xr + xl*xl)
           endif
        enddo
        
        p_zone = (vol_pert*p_exp + vol_ambient*p_ambient)/(vol_pert+vol_ambient)
  
        eint = p_zone/(gamma_const - 1.d0)
   
        state(i,URHO) = dens_ambient
        state(i,UMX) = 0.d0
        
        state(i,UEDEN) = eint + 0.5d0 * state(i,UMX)**2 / state(i,URHO)
        state(i,UEINT) = eint

        state(i,UFS) = state(i,URHO)
  
     end do
  else

     call bl_error('dont know this probtype in initdata')

  end if

end subroutine ca_initdata


! ::: 
! ::: -----------------------------------------------------------
! ::: 
subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                      domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only : NVAR

  implicit none

  include 'bc_types.fi'

  integer :: adv_l1,adv_h1
  integer :: bc(1,2,*)
  integer :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1,NVAR)
  
  double precision :: state(NVAR)
  double precision :: staten(NVAR)
  
  integer :: i, n
  logical rho_only

  do n = 1,NVAR
     call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  ! The strategy here is to set Dirichlet condition for inflow and
  ! outflow boundaries, and let the Riemann solver sort out the proper
  ! upwinding.  However, this decision makes this routine look
  ! somewhat non-orthodox, in that we need to set external values in
  ! either case....how do we know it's Outflow?  We have to assume
  ! that the setup routines converted Outflow to FOEXTRAP.

  ! Set flag for bc function
  rho_only = .FALSE.

  !     XLO
  if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        do n=1,NVAR
           state(n) = adv(domlo(1),n)
        enddo
        call bcnormal(state,staten,rho_only)
        do n=1,NVAR
           adv(i,n) = staten(n)
        enddo
     end do
  end if

  !     XHI
  if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        do n=1,NVAR
           state(n) = adv(domhi(1),n)
        enddo
        call bcnormal(state,staten,rho_only)
        do n=1,NVAR
           adv(i,n) = staten(n)
        enddo
     end do
  end if
  
end subroutine ca_hypfill


! ::: 
! ::: -----------------------------------------------------------
! ::: 
subroutine ca_denfill(adv,adv_l1,adv_h1, &
                      domlo,domhi,delta,xlo,time,bc)
  implicit none

  include 'bc_types.fi'

  integer :: adv_l1,adv_h1
  integer :: bc(1,2,*)
  integer :: domlo(1), domhi(1)
  double precision :: delta(1), xlo(1), time
  double precision :: adv(adv_l1:adv_h1)
  logical rho_only
  integer :: i

  ! Note: this function should not be needed, technically, but is
  ! provided to filpatch because there are many times in the algorithm
  ! when just the density is needed.  We try to rig up the filling so
  ! that the same function is called here and in hypfill where all the
  ! states are filled.

  call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

  rho_only = .TRUE.

  !     XLO
  if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        call bcnormal(adv(domlo(1)),adv(i),rho_only)
     end do
  end if
  
  !     XHI
  if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        call bcnormal(adv(domhi(1)),adv(i),rho_only)
     end do
  end if
  
end subroutine ca_denfill


! ::: 
! ::: -----------------------------------------------------------
! ::: 
subroutine bcnormal(u_int,u_ext,rho_only)

  use probdata_module
  use eos_module, only: gamma_const
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT

  implicit none

  double precision, intent(in   ) ::  u_int(*)
  double precision, intent(  out) ::  u_ext(*)
  logical         , intent(in   ) :: rho_only

  ! Local variables
  integer :: i

  ! For the Sedov problem, we will always set the state to the ambient conditions
  
  if (rho_only .EQV. .TRUE. ) then

     u_ext(1) = dens_ambient

  else

     ! First set everything from internal data (this is probably a bad
     ! thing to do...)  That is, we should have explicit boundary data
     ! for advected fields and species

     do i=1,NVAR
        u_ext(i) = u_int(i)
     enddo

     u_ext(URHO ) = dens_ambient
     u_ext(UMX  ) = 0.d0
     u_ext(UEDEN) = p_ambient/(gamma_const-1.d0)
     u_ext(UEINT) = u_ext(UEDEN)
     
  endif

end subroutine bcnormal

