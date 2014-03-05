
subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use network, only : network_init
  implicit none

  integer init, namlen
  integer name(namlen)
  double precision problo(1), probhi(1)
  
  integer untin,i
  character(1) dummy
  
  namelist /fortin/ & 
       denerr,   dengrad,  max_denerr_lev,  max_dengrad_lev, &
       velerr,   velgrad,  max_velerr_lev,  max_pressgrad_lev, &
       presserr, pressgrad,max_presserr_lev,max_pressgrad_lev, &
       temperr,  tempgrad, max_temperr_lev, max_tempgrad_lev, &
       raderr,   radgrad,  max_raderr_lev,  max_radgrad_lev, &
       rwind0, rwind1, rhowind1, Twind1, rbasefac, filter_rhomax, filter_timemax
  
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
  
  rbasefac = 0.99d0
  rwind0 = 0.7d14
  rwind1 = 1.d14
  rhowind1 = 1.d-14
  Twind1 = 1.1d3

  filter_rhomax = -1.d20
  filter_timemax = -1.d20

  center(1) = 0.d0
  xmin = problo(1)
  xmax = probhi(1)

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  open(unit=9,file="model.input")
  print*,'reading model inputs'
  read(9,*) dummy
  read(9,*) npts_model
  read(9,*) dummy
  if (npts_model > npts_max) then
     write(6,*) 'npts_max in probdata.f90 is too small'
     stop
  end if
  do i = 1, npts_model
     read(9,*)model_r(i), model_rho(i), model_v(i), &
          model_T(i), model_Ye(i), model_Abar(i)
  enddo
  print*,'done reading model inputs'
  
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
  use interpolate_module
  use fundamental_constants_module, only: k_B, n_A
  
  implicit none
  
  integer :: level, nscal
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  double precision :: state(state_l1:state_h1,NVAR)
  double precision :: time, delta(1)
  double precision :: xlo(1), xhi(1)
  
  integer :: i, ii
  double precision :: rhoInv, xcl, xx, vtot, vsub, rho, T, u, rhosub, Tsub, usub, dx_sub
  double precision :: rhowind0, rlast, rholast, Twind0, Tlast
  integer, parameter :: nsub = 16
  double precision :: rho_tmp, rbase, T_tmp
  double precision :: Ye, Abar, invmu
  double precision, parameter :: Tindex=0.5
  type(eos_t) :: eos_state

  if (naux .ne. 2) then
     write(6,*) 'naux in network is not equal to 2'
     stop
  end if
  if (nspec .ne. 1) then
     write(6,*) 'nspec in network is not equal to 1'
     stop
  end if

  dx_sub = delta(1) / dble(nsub)

  rhowind0 = rhowind1 * (rwind1/rwind0)**2
  rlast = model_r(npts_model)
  rholast = model_rho(npts_model)
  rbase = rlast * rbasefac

  Twind0 = Twind1 * (rwind1/rwind0)**Tindex
  Tlast = model_T(npts_model)
  
  do i = lo(1), hi(1)
  
     xcl = xlo(1) + delta(1) * float(i-lo(1))
     
     vtot = 0.d0
     rho = 0.d0
     T = 0.d0
     u = 0.d0
     Ye = 0.d0
     Abar = 0.d0
     do ii=0,nsub-1
        xx = xcl + (dble(ii)+0.5d0) * dx_sub
        vsub = xx**2
        vtot = vtot + vsub
        if (xx .ge. model_r(npts_model)) then
           if (xx .ge. rwind0 ) then 
              rho_tmp = rhowind1 * (rwind1/xx)**2
           else
              rho_tmp = rholast * (rhowind0/rholast)** &
                   ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
           end if
           rho = rho + rho_tmp * vsub

           if (xx .ge. rwind0 ) then 
              T_tmp = Twind1 * (rwind1/xx)**Tindex
           else
              T_tmp = Tlast * (Twind0/Tlast)** &
                   ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
           end if
           T = T + vsub * T_tmp

           Ye = Ye + model_Ye(npts_model) * vsub
           Abar = Abar + model_Abar(npts_model) * vsub

           u = u + 0.d0
        else if (xx .le. model_r(1)) then
           rho = rho + model_rho(1) * vsub
           T = T + model_T(1) * vsub
           u = u + 0.d0
           Ye = Ye + model_Ye(1) * vsub
           Abar = Abar + model_Abar(1) * vsub
        else
           rho = rho + interpolate(xx,npts_model,model_r,model_rho) * vsub
           T   = T   + interpolate(xx,npts_model,model_r,model_T  ) * vsub
           u   = u   + interpolate(xx,npts_model,model_r,model_v  ) * vsub
           Ye  = Ye  + interpolate(xx,npts_model,model_r,model_Ye ) * vsub
           Abar=Abar + interpolate(xx,npts_model,model_r,model_Abar) * vsub
        end if
     end do
     rho = rho / vtot
     T = T / vtot
     u = u / vtot
     Ye = Ye / vtot
     Abar = Abar / vtot

     invmu = (1.d0+Abar*Ye)/Abar

     state(i,URHO)  = rho
     state(i,UTEMP)  = T
     state(i,UMX)   = rho * u

     ! set the composition to be all in the first species
     state(i,UFS:UFS-1+nspec) = 0.d0
     state(i,UFS  ) = state(i,URHO)
     state(i,UFX) = Ye*rho
     state(i,UFX+1) = invmu*rho

     ! set the internal energy via the EOS
     rhoInv = 1.d0 / state(i,URHO)
     eos_state % rho = state(i,URHO)
     eos_state % T   = state(i,UTEMP)
     eos_state % xn  = state(i,UFS:UFS+nspec-1) * rhoInv
     eos_state % aux = state(i,UFX:UFX+naux-1) * rhoInv

     call eos(eos_input_rt, eos_state)
     
     state(i,UEINT) = state(i,URHO) * eos_state % e
     state(i,UEDEN) = state(i,UEINT) + & 
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
  use interpolate_module

  implicit none
  integer :: level, nrad
  integer :: lo(1), hi(1)
  integer :: rad_state_l1,rad_state_h1
  double precision :: xlo(1), xhi(1), time, delta(1)
  double precision ::  rad_state(rad_state_l1:rad_state_h1, nrad)

  ! local variables
  integer :: i,ii
  double precision xcl, T, Tsub, xx, vtot, vsub, dx_sub
  integer, parameter :: nsub = 16

  double precision :: rlast, rbase, T_tmp, Twind0, Tlast
  double precision, parameter :: Tindex=0.5

  dx_sub = delta(1) / dble(nsub)

  rlast = model_r(npts_model) 
  rbase = rlast * rbasefac

  Twind0 = Twind1 * (rwind1/rwind0)**Tindex
  Tlast = model_T(npts_model)
            
  do i = lo(1), hi(1)

     xcl = xlo(1) + delta(1) * float(i-lo(1))
     
     vtot = 0.d0
     T = 0.d0
     do ii=0,nsub-1
        xx = xcl + (dble(ii)+0.5d0) * dx_sub
        vsub = xx**2
        vtot = vtot + vsub
        if (xx .ge. model_r(npts_model)) then
           if (xx .ge. rwind0 ) then 
              T_tmp = Twind1 * (rwind1/xx)**Tindex
           else
              T_tmp = Tlast * (Twind0/Tlast)** &
                   ((log(xx-rbase)-log(rlast-rbase))/(log(rwind0-rbase)-log(rlast-rbase)))
           end if
           T = T + vsub * T_tmp

        else if (xx .le. model_r(1)) then
           T = T + model_T(1) * vsub
        else
           T   = T   + interpolate(xx,npts_model,model_r,model_T  ) * vsub
        end if
     end do
     T = T / vtot
     ! set radiation energy density to a T**4
     rad_state(i,:) = a_rad*T**4

  enddo
        
end subroutine ca_initrad


! ::: -----------------------------------------------------------

      subroutine ca_hypfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)

      use meth_params_module, only : NVAR

      implicit none

      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1,NVAR)

      integer n

      do n = 1,NVAR
         call filcc(adv(adv_l1,n),adv_l1,adv_h1, &
             domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_denfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

      end

! ::: 
! ::: -----------------------------------------------------------
! :::

      subroutine ca_radfill(adv,adv_l1,adv_h1, &
          domlo,domhi,delta,xlo,time,bc)
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_h1
      integer bc(1,2,*)
      integer domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision adv(adv_l1:adv_h1)

      call filcc(adv,adv_l1,adv_h1,domlo,domhi,delta,xlo,bc)

      end

!-----------------------------------------------------------------------

      subroutine ca_gravxfill(grav,grav_l1,grav_h1,&
                              domlo,domhi,delta,xlo,time,bc)
 
      use probdata_module
      implicit none
      include 'bc_types.fi'
 
      integer :: grav_l1,grav_h1
      integer :: bc(1,2,*)
      integer :: domlo(1), domhi(1)
      double precision delta(1), xlo(1), time
      double precision grav(grav_l1:grav_h1)
 
      double precision :: ri,rim1
      integer          :: i
 
      call filcc(grav,grav_l1,grav_h1,domlo,domhi,delta,xlo,bc)
 
!     Outflow boundary condition for gravity at hi r
      if ( bc(1,2,1).eq. FOEXTRAP .and. grav_h1.gt.domhi(1)) then
         do i = domhi(1)+1, grav_h1
            ri   = (dble(i  )+0.5d0) * delta(1)
            rim1 = (dble(i-1)+0.5d0) * delta(1)
            grav(i) = grav(i-1) * (ri/rim1)**2
         end do
      end if
 
      end subroutine ca_gravxfill

! ::: -----------------------------------------------------------
