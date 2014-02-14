subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use eos_module
  use network, only : nspec
  use meth_params_module, only : small_temp

  implicit none 

  integer :: init,namlen,untin,i,k
  integer :: name(namlen)

  double precision :: center_x, center_y, center_z
  double precision :: problo(3), probhi(3)

  type (eos_t) :: eos_state

  namelist /fortin/ &
       rho_0, r_0, r_old, p_0, rho_ambient, smooth_delta, &
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

  is_3d_fullstar = .false.

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
  r_old = r_0
  p_0 = 1.d10
  rho_ambient = 1.d0
  smooth_delta = 1.d-5

!     Read namelists in probin file
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! in 3-d, we center the sphere at (center_x, center_y, center_z)
  center(1) = center_x
  center(2) = center_y
  center(3) = center_z

  xmin = problo(1)
  xmax = probhi(1)

  ymin = problo(2)
  ymax = probhi(2)

  zmin = problo(3)
  zmax = probhi(3)

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
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)

  use probdata_module
  use eos_module
  use network, only : nspec
  use interpolate_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEDEN, UEINT, UFS, small_temp

  implicit none

  integer          :: level, nscal
  integer          :: lo(3), hi(3)
  integer          :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  double precision :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
  double precision :: time, delta(3)
  double precision :: xlo(3), xhi(3)
  
  double precision :: xl,yl,zl,xx,yy,zz,dist,pres,eint,temp,avg_rho,rho_n,volinv
  double precision :: dx_sub,dy_sub,dz_sub
  integer          :: i,j,k,ii,jj,kk,n

  type (eos_t) :: eos_state

  integer, parameter :: nsub = 5

  volinv = 1.d0/dble(nsub*nsub*nsub)

  dx_sub = delta(1)/dble(nsub)
  dy_sub = delta(2)/dble(nsub)
  dz_sub = delta(3)/dble(nsub)

  do k = lo(3), hi(3)
    zl = zmin + dble(k) * delta(3)

    do j = lo(2), hi(2)
    yl = ymin + dble(j) * delta(2)

     do i = lo(1), hi(1)
        xl = xmin + dble(i) * delta(1)

        avg_rho = 0.d0

        do kk = 0, nsub-1
           zz = zl + (dble(kk) + 0.5d0) * dz_sub

        do jj = 0, nsub-1
           yy = yl + (dble(jj) + 0.5d0) * dy_sub

           do ii = 0, nsub-1
              xx = xl + (dble(ii) + 0.5d0) * dx_sub

              dist = sqrt((xx-center(1))**2 + (yy-center(2))**2 + (zz-center(3))**2)

              ! use a tanh profile to smooth the transition between rho_0 
              ! and rho_ambient
              rho_n = rho_0 - 0.5d0*(rho_0 - rho_ambient)* &
                   (1.d0 + tanh((dist - r_0)/smooth_delta))

              avg_rho = avg_rho + rho_n

           enddo
        enddo
        enddo
        
        state(i,j,k,URHO) = avg_rho * volinv

        eos_state % rho = state(i,j,k,URHO)
        eos_state % p   = p_0
        eos_state % T   = small_temp ! Initial guess for the EOS
        eos_state % xn  = X_0

        call eos(eos_input_rp, eos_state)

        temp = eos_state % T
        eint = eos_state % e

        state(i,j,k,UTEMP) = temp
        state(i,j,k,UMX) = 0.d0
        state(i,j,k,UMY) = 0.d0
        state(i,j,k,UMZ) = 0.d0
        state(i,j,k,UEDEN) = state(i,j,k,URHO) * eint
        state(i,j,k,UEINT) = state(i,j,k,URHO) * eint
        state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO) * X_0(1:nspec)

     enddo
    enddo
  enddo

end subroutine ca_initdata

! ::: -----------------------------------------------------------

      subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
                            domlo,domhi,delta,xlo,time,bc)

      use probdata_module, only : center
      use meth_params_module, only : NVAR,UMX,UMY,UMZ

      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

      integer          :: i, j, k, n
      integer          :: ic,jc,kc
      double precision :: x,y,z,r
      double precision :: xc,yc,zc,rc
      double precision :: mom,momc

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
              adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      ! Do this for all the variables, but we will overwrite the momenta below
      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
              adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      if ( (bc(1,1,1).eq. EXT_DIR .or. bc(1,2,1).eq. EXT_DIR) .or.  &
           (bc(2,1,1).eq. EXT_DIR .or. bc(2,2,1).eq. EXT_DIR) .or. &
           (bc(3,1,1).eq. EXT_DIR .or. bc(3,2,1).eq. EXT_DIR) ) then
         print *,'NOT SET UP FOR EXT_DIR BCs IN HYPFILL'
         stop
      end if

!     XLO
      if ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1)) then
         do k = adv_l3, adv_h3
         do j = adv_l2, adv_h2

            y = (dble(j) + 0.5d0) * delta(2) - center(2)
            z = (dble(k) + 0.5d0) * delta(3) - center(3)

            ic = domlo(1)
            xc = (dble(ic) + 0.5d0) * delta(1) - center(1)
            rc = sqrt(xc**2 + y**2 + z**2)

            momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

            do i = adv_l1, domlo(1)-1
               x = (dble(i) + 0.5d0) * delta(1) - center(1)
               r = sqrt(x**2 + y**2 + z**2)

               mom  = momc * (rc/r)**2

               ! Project along the normal
               adv(i,j,k,UMX) =  sign(1.d0,adv(ic,j,k,UMX)) * mom * (x/r)
               adv(i,j,k,UMY) =  sign(1.d0,adv(ic,j,k,UMY)) * mom * (y/r)
               adv(i,j,k,UMZ) =  sign(1.d0,adv(ic,j,k,UMZ)) * mom * (z/r)

            end do
         end do
         end do
      end if            

!     XHI
      if ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1)) then
         do k = adv_l3, adv_h3
         do j = adv_l2, adv_h2

            y = (dble(j) + 0.5d0) * delta(2) - center(2)
            z = (dble(k) + 0.5d0) * delta(3) - center(3)

            ic = domhi(1)
            xc = (dble(ic) + 0.5d0) * delta(1) - center(1)
            rc = sqrt(xc**2 + y**2 + z**2)

            momc = sqrt(adv(ic,j,k,UMX)**2 + adv(ic,j,k,UMY)**2 + adv(ic,j,k,UMZ)**2)

            do i = domhi(1)+1, adv_h1
               x = (dble(i) + 0.5d0) * delta(1) - center(1)
               r = sqrt(x**2 + y**2 + z**2)

               mom  = momc * (rc/r)**2

               ! Project along the normal
               adv(i,j,k,UMX) =  sign(1.d0,adv(ic,j,k,UMX)) * mom * (x/r)
               adv(i,j,k,UMY) =  sign(1.d0,adv(ic,j,k,UMY)) * mom * (y/r)
               adv(i,j,k,UMZ) =  sign(1.d0,adv(ic,j,k,UMZ)) * mom * (z/r)

            end do
	 end do
	 end do
      end if            

!     YLO
      if ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2)) then
         do k = adv_l3, adv_h3
         do i = adv_l1, adv_h1

            x = (dble(i) + 0.5d0) * delta(1) - center(1)
            z = (dble(k) + 0.5d0) * delta(3) - center(3)

            jc = domlo(2)
            yc = (dble(jc) + 0.5d0) * delta(2) - center(2)
            rc = sqrt(x**2 + yc**2 + z**2)

            momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

            do j = adv_l2, domlo(2)-1

               y = (dble(j) + 0.5d0) * delta(2) - center(2)
               r = sqrt(x**2 + y**2 + z**2)

               mom  = momc * (rc/r)**2

               ! Project along the normal
               adv(i,j,k,UMX) =  sign(1.d0,adv(i,jc,k,UMX)) * mom * (x/r)
               adv(i,j,k,UMY) =  sign(1.d0,adv(i,jc,k,UMY)) * mom * (y/r)
               adv(i,j,k,UMZ) =  sign(1.d0,adv(i,jc,k,UMZ)) * mom * (z/r)

            end do
         end do
         end do
      end if            

!     YHI
      if ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2)) then
         do k = adv_l3, adv_h3
         do i = adv_l1, adv_h1

            x = (dble(i) + 0.5d0) * delta(1) - center(1)
            z = (dble(k) + 0.5d0) * delta(3) - center(3)

            jc = domhi(2)
            yc = (dble(jc) + 0.5d0) * delta(2) - center(2)
            rc = sqrt(x**2 + yc**2 + z**2)

            momc = sqrt(adv(i,jc,k,UMX)**2 + adv(i,jc,k,UMY)**2 + adv(i,jc,k,UMZ)**2)

            do j = domhi(2)+1, adv_h2
               y = (dble(j) + 0.5d0) * delta(2) - center(2)
               r = sqrt(x**2 + y**2 + z**2)

               mom  = momc * (rc/r)**2

               ! Project along the normal
               adv(i,j,k,UMX) =  sign(1.d0,adv(i,jc,k,UMX)) * mom * (x/r)
               adv(i,j,k,UMY) =  sign(1.d0,adv(i,jc,k,UMY)) * mom * (y/r)
               adv(i,j,k,UMZ) =  sign(1.d0,adv(i,jc,k,UMZ)) * mom * (z/r)

            end do
	 end do
	 end do
      end if            

!     ZLO
      if ( bc(3,1,1).eq. FOEXTRAP .and. adv_l3.lt.domlo(3)) then
         do j = adv_l2, adv_h2
         do i = adv_l1, adv_h1

            x = (dble(i) + 0.5d0) * delta(1) - center(1)
            y = (dble(j) + 0.5d0) * delta(2) - center(2)

            kc = domlo(3)
            zc = (dble(kc) + 0.5d0) * delta(3) - center(3)
            rc = sqrt(x**2 + y**2 + zc**2)

            momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

            do k = adv_l3, domlo(3)-1
               z = (dble(k) + 0.5d0) * delta(3) - center(3)
               r = sqrt(x**2 + y**2 + z**2)

               mom  = momc * (rc/r)**2

               ! Project along the normal
               adv(i,j,k,UMX) =  sign(1.d0,adv(i,j,kc,UMX)) * mom * (x/r)
               adv(i,j,k,UMY) =  sign(1.d0,adv(i,j,kc,UMY)) * mom * (y/r)
               adv(i,j,k,UMZ) =  sign(1.d0,adv(i,j,kc,UMZ)) * mom * (z/r)

            end do
         end do
         end do
      end if            

!     ZHI
      if ( bc(3,2,1).eq. FOEXTRAP .and. adv_h3.gt.domhi(3)) then
         do j = adv_l2, adv_h2
         do i = adv_l1, adv_h1

            x = (dble(i) + 0.5d0) * delta(1) - center(1)
            y = (dble(j) + 0.5d0) * delta(2) - center(2)

            kc = domhi(3)
            zc = (dble(kc) + 0.5d0) * delta(3) - center(3)
            rc = sqrt(x**2 + y**2 + zc**2)

            momc = sqrt(adv(i,j,kc,UMX)**2 + adv(i,j,kc,UMY)**2 + adv(i,j,kc,UMZ)**2)

            do k = domhi(3)+1, adv_h3
               z = (dble(k) + 0.5d0) * delta(3) - center(3)
               r = sqrt(x**2 + y**2 + z**2)

               mom  = momc * (rc/r)**2

               ! Project along the normal
               adv(i,j,k,UMX) =  sign(1.d0,adv(i,j,kc,UMX)) * mom * (x/r)
               adv(i,j,k,UMY) =  sign(1.d0,adv(i,j,kc,UMY)) * mom * (y/r)
               adv(i,j,k,UMZ) =  sign(1.d0,adv(i,j,kc,UMZ)) * mom * (z/r)

            end do
	 end do
	 end do
      end if            

      end subroutine ca_hypfill

! ::: -----------------------------------------------------------

      subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                            adv_h3,domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      double precision delta(3), xlo(3), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
      logical rho_only
      integer i,j,k

      call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,delta,xlo,bc)

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
      double precision delta(3), xlo(3), time
      double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

      call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

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
      double precision delta(3), xlo(3), time
      double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

      call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

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
      double precision delta(3), xlo(3), time
      double precision grav(grav_l1:grav_h1,grav_l2:grav_h2,grav_l3:grav_h3)

      call filcc(grav,grav_l1,grav_l2,grav_l3,grav_h1,grav_h2,grav_h3,domlo,domhi,delta,xlo,bc)

      end subroutine ca_gravzfill
