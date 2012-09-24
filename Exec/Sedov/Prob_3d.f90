subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use bl_error_module

  implicit none
  integer :: init, namlen
  integer :: name(namlen)
  double precision :: problo(3), probhi(3)

  integer :: untin,i

  namelist /fortin/ probtype, p_ambient, dens_ambient, exp_energy, &
       r_init, nsub,  &
       denerr,dengrad,max_denerr_lev,max_dengrad_lev, &
       presserr,pressgrad,max_presserr_lev,max_pressgrad_lev
  
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

  denerr = 1.d20
  dengrad = 1.d20
  max_denerr_lev = -1
  max_dengrad_lev = -1

  presserr = 1.d20
  pressgrad = 1.d20
  max_presserr_lev = -1
  max_pressgrad_lev = -1

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! set local variable defaults
  center(1) = (problo(1)+probhi(1))/2.d0
  center(2) = (problo(2)+probhi(2))/2.d0
  center(3) = (problo(3)+probhi(3))/2.d0
  
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
  use eos_module, only : gamma_const
  use bl_constants_module, only: M_PI, FOUR3RD
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS

  implicit none
  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_l1:state_h1, &
                            state_l2:state_h2, &
                            state_l3:state_h3,NVAR)

  double precision :: xmin,ymin,zmin
  double precision :: xx, yy, zz
  double precision :: dist
  double precision :: eint, p_zone
  double precision :: vctr, p_exp

  integer :: i,j,k, ii, jj, kk
  integer :: npert, nambient
  
  if (probtype .eq. 32) then

     ! set explosion pressure -- we will convert the point-explosion
     ! energy into a corresponding pressure distributed throughout the
     ! perturbed volume
     vctr  = M_PI*r_init**2
     p_exp = (gamma_const - 1.d0)*exp_energy/vctr
     
     do k = lo(3), hi(3)
        zmin = xlo(3) + delta(3)*dble(k-lo(3)) 
        
        do j = lo(2), hi(2)
           ymin = xlo(2) + delta(2)*dble(j-lo(2))
           
           do i = lo(1), hi(1)
              xmin = xlo(1) + delta(1)*dble(i-lo(1))
              
              npert = 0
              nambient = 0
              
              do jj = 0, nsub-1
                 yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)
                 
                 do ii = 0, nsub-1
                    xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)
                    
                    dist = (center(1)-xx)**2 + (center(2)-yy)**2
                    
                    if(dist <= r_init**2) then
                       npert = npert + 1
                    else
                       nambient = nambient + 1
                    endif
                    
                 enddo
              enddo
              
              p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient) / &
                       (dble(npert) + dble(nambient))

              eint = p_zone/(gamma_const - 1.d0)

              state(i,j,k,URHO) = dens_ambient
              state(i,j,k,UMX) = 0.d0
              state(i,j,k,UMY) = 0.d0
              state(i,j,k,UMZ) = 0.d0
              
              state(i,j,k,UEDEN) = eint +  &
                   0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMZ)**2/state(i,j,k,URHO))

              state(i,j,k,UEINT) = eint

              state(i,j,k,UFS) = state(i,j,k,URHO)

           enddo
        enddo
     enddo
     
  else if (probtype .eq. 33) then

     ! set explosion pressure -- we will convert the point-explosion energy into
     ! a corresponding pressure distributed throughout the perturbed volume
     vctr  = FOUR3RD*M_PI*r_init**3
     p_exp = (gamma_const - 1.d0)*exp_energy/vctr

     do k = lo(3), hi(3)
        zmin = xlo(3) + delta(3)*dble(k-lo(3)) 

        do j = lo(2), hi(2)
           ymin = xlo(2) + delta(2)*dble(j-lo(2))
           
           do i = lo(1), hi(1)
              xmin = xlo(1) + delta(1)*dble(i-lo(1))

              npert = 0
              nambient = 0

              do kk = 0, nsub-1
                 zz = zmin + (delta(3)/dble(nsub))*(kk + 0.5d0)
                 
                 do jj = 0, nsub-1
                    yy = ymin + (delta(2)/dble(nsub))*(jj + 0.5d0)
                    
                    do ii = 0, nsub-1
                       xx = xmin + (delta(1)/dble(nsub))*(ii + 0.5d0)
                       
                       dist = (center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2
                       
                       if(dist <= r_init**2) then
                          npert = npert + 1
                       else
                          nambient = nambient + 1
                       endif
                       
                    enddo
                 enddo
              enddo
              
              p_zone = (dble(npert)*p_exp + dble(nambient)*p_ambient)/  &
                   dble(nsub*nsub*nsub)

              eint = p_zone/(gamma_const - 1.d0)

              state(i,j,k,URHO) = dens_ambient
              state(i,j,k,UMX) = 0.d0
              state(i,j,k,UMY) = 0.d0
              state(i,j,k,UMZ) = 0.d0
              
              state(i,j,k,UEDEN) = eint + &
                   0.5d0*(state(i,j,k,UMX)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMY)**2/state(i,j,k,URHO) + &
                          state(i,j,k,UMZ)**2/state(i,j,k,URHO))

              state(i,j,k,UEINT) = eint

              state(i,j,k,UFS) = state(i,j,k,URHO)

           enddo
        enddo
     enddo
     
  else 

     call bl_error('Dont know this probtype in initdata')

  end if
  
end subroutine ca_initdata

      
! ::: 
! ::: -----------------------------------------------------------
! ::: 
subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                      adv_h3,domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only : NVAR

  implicit none
  include 'bc_types.fi'

  integer :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer :: bc(3,2,*)
  integer :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)
  
  double precision :: state(NVAR)
  double precision :: staten(NVAR)
  
  integer :: i, j, k, n, lo(3), hi(3)
  double precision :: x, y, z
  logical rho_only
  
  lo(1) = adv_l1
  lo(2) = adv_l2
  lo(3) = adv_l3
  hi(1) = adv_h1
  hi(2) = adv_h2
  hi(3) = adv_h3

  do n = 1,NVAR
     call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
                adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                domlo,domhi,delta,xlo,bc(1,1,n))
  enddo

  ! The strategy here is to set Dirichlet condition for inflow and
  ! outflow boundaries, and let the Riemann solver sort out the proper
  ! upwinding.  However, this decision makes this routine look
  ! somewhat non-orthodox, in that we need to set external values in
  ! either case....how do we know it's Outflow?  We have to assume
  ! that the setup routines converted Outflow to FOEXTRAP.

  !     Set flag for bc function
  rho_only = .FALSE.

  !     XLO
  if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        do j = adv_l2, adv_h2
           do k = adv_l3, adv_h3
              do n=1,NVAR
                 state(n) = adv(domlo(1),j,k,n)
              enddo
              call bcnormal(state,staten,1,+1,rho_only)
              do n=1,NVAR
                 adv(i,j,k,n) = staten(n)
              enddo
           end do
        end do
     end do
  end if
  
  !     XHI
  if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        do j = adv_l2, adv_h2
           do k = adv_l3, adv_h3
              do n=1,NVAR
                 state(n) = adv(domhi(1),j,k,n)
              enddo
              call bcnormal(state,staten,1,-1,rho_only)
              do n=1,NVAR
                 adv(i,j,k,n) = staten(n)
              enddo
           end do
        end do
     end do
  end if
  
  !     YLO
  if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and.adv_l2.lt.domlo(2)) then
     do i = adv_l1, adv_h1
        do j = adv_l2, domlo(2)-1
           do k = adv_l3, adv_h3
              do n=1,NVAR
                 state(n) = adv(i,domlo(2),k,n)
              enddo
              call bcnormal(state,staten,2,+1,rho_only)
              do n=1,NVAR
                 adv(i,j,k,n) = staten(n)
              enddo
           end do
        end do
     end do
  end if
  
  !     YHI
  if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and.adv_h2.gt.domhi(2)) then
     do i = adv_l1, adv_h1
        do j = domhi(2)+1, adv_h2
           do k = adv_l3, adv_h3
              do n=1,NVAR
                 state(n) = adv(i,domhi(2),k,n)
              enddo
              call bcnormal(state,staten,2,-1,rho_only)
              do n=1,NVAR
                 adv(i,j,k,n) = staten(n)
              enddo
           end do
        end do
     end do
  end if

  !     ZLO
  if ( (bc(3,1,1).eq.EXT_DIR.or.bc(3,1,1).eq.FOEXTRAP).and.adv_l3.lt.domlo(3)) then
     do i = adv_l1, adv_h1
        do j = adv_l2, adv_h2
           do k = adv_l3, domlo(3)-1
              do n=1,NVAR
                 state(n) = adv(i,j,domlo(3),n)
              enddo
              call bcnormal(state,staten,3,+1,rho_only)
              do n=1,NVAR
                 adv(i,j,k,n) = staten(n)
              enddo
           end do
        end do
     end do
  end if

  !     ZHI
  if ( (bc(3,2,1).eq.EXT_DIR.or.bc(3,2,1).eq.FOEXTRAP).and.adv_h3.gt.domhi(3)) then
     do i = adv_l1, adv_h1
        do j = adv_l2, adv_h2
           do k = domhi(3)+1, adv_h3
              do n=1,NVAR
                 state(n) = adv(i,j,domhi(3),n)
              enddo
              call bcnormal(state,staten,3,-1,rho_only)
              do n=1,NVAR
                 adv(i,j,k,n) = staten(n)
              enddo
           end do
        end do
     end do
  end if
  
end subroutine ca_hypfill

! ::: 
! ::: -----------------------------------------------------------
! ::: 
subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                      adv_h3,domlo,domhi,delta,xlo,time,bc)

  implicit none
  include 'bc_types.fi'
  integer :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
  integer :: bc(3,2,*)
  integer :: domlo(3), domhi(3)
  double precision :: delta(3), xlo(3), time
  double precision :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
  logical rho_only
  integer :: i,j,k
  
  ! Note: this function should not be needed, technically, but is
  ! provided to filpatch because there are many times in the algorithm
  ! when just the density is needed.  We try to rig up the filling so
  ! that the same function is called here and in hypfill where all the
  ! states are filled.

  call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
             domlo,domhi,delta,xlo,bc)

  rho_only = .TRUE.

  !     XLO
  if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
     do i = adv_l1, domlo(1)-1
        do j = adv_l2, adv_h2
           do k = adv_l3, adv_h3
              call bcnormal(adv(domlo(1),j,k),adv(i,j,k),1,+1,rho_only)
           end do
        end do
     end do
  end if
  
  !     XHI
  if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
     do i = domhi(1)+1, adv_h1
        do j = adv_l2, adv_h2
           do k = adv_l3, adv_h3
              call bcnormal(adv(domhi(1),j,k),adv(i,j,k),1,-1,rho_only)
           end do
        end do
     end do
  end if
  
  !     YLO
  if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and.adv_l2.lt.domlo(2)) then
     do i = adv_l1, adv_h1
        do j = adv_l2, domlo(2)-1
           do k = adv_l3, adv_h3
              call bcnormal(adv(i,domlo(2),k),adv(i,j,k),2,+1,rho_only)
           end do
        end do
     end do
  end if
  
  !     YHI
  if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and.adv_h2.gt.domhi(2)) then
     do i = adv_l1, adv_h1
        do j = domhi(2)+1, adv_h2
           do k = adv_l3, adv_h3
              call bcnormal(adv(i,domhi(2),k),adv(i,j,k),2,-1,rho_only)
           end do
        end do
     end do
  end if
  
  !     ZLO
  if ( (bc(3,1,1).eq.EXT_DIR.or.bc(3,1,1).eq.FOEXTRAP).and.adv_l3.lt.domlo(3)) then
     do i = adv_l1, adv_h1
        do j = adv_l2, adv_h2
           do k = adv_l3, domlo(3)-1
              call bcnormal(adv(i,j,domlo(3)),adv(i,j,k),3,+1,rho_only)
           end do
        end do
     end do
  end if
  
  !     ZHI
  if ( (bc(3,2,1).eq.EXT_DIR.or.bc(3,2,1).eq.FOEXTRAP).and.adv_h3.gt.domhi(3)) then
     do i = adv_l1, adv_h1
        do j = adv_l2, adv_h2
           do k = domhi(3)+1, adv_h3
              call bcnormal(adv(i,j,domhi(3)),adv(i,j,k),3,-1,rho_only)
           end do
        end do
     end do
  end if
end subroutine ca_denfill


! ::: 
! ::: -----------------------------------------------------------
! ::: 
subroutine bcnormal(u_int,u_ext,dir,sgn,rho_only)

  use probdata_module
  use eos_module, only : gamma_const
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT
  implicit none
  
  double precision :: u_int(*),u_ext(*)
  logical rho_only
  integer :: dir,sgn
  double precision :: rho, rhou(3), eden, T, Y
  integer :: n,t1,t2,i
  
  ! for the Sedov problem, we will always set the state to the ambient conditions

  if (rho_only .EQV. .TRUE. ) then
     
     u_ext(1) = dens_ambient
     
  else
     
     ! First set everything from internal data (this is probably a bad
     ! thing to do...)  That is, we should have explicit boundary data
     ! for advected fields and species

     do i=1,NVAR
        u_ext(i) = u_int(i)
     enddo
     
     u_ext(URHO)   = dens_ambient
     u_ext(UMX)    = 0.d0
     u_ext(UMY)    = 0.d0
     u_ext(UMZ)    = 0.d0
     u_ext(UEDEN)  = p_ambient/(gamma_const-1.d0)
     u_ext(UEINT)  = u_ext(UEDEN)
     
  endif

end subroutine bcnormal
      
