subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_error_module
  use probdata_module
  use prob_params_module, only: center
  use eos_module
  use eos_data_module
  use eos_type_module

  use network, only : nspec
  implicit none

  integer :: init,namlen,untin,i,k
  integer :: name(namlen)

  double precision problo(1), probhi(1)

  type (eos_t) :: eos_state

  namelist /fortin/ &
       model_name,  &
       heating_time, heating_rad, heating_peak, heating_sigma, &
       prob_type, &
       sponge_weighting, sponge_start_density, sponge_width_factor


  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer   :: ipos
  integer, parameter :: maxlen=127
  character probin*(maxlen)
  character (len=256) :: header_line

  integer :: nvars_model_file

  nvars_model_file = 3 + nspec

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  heating_time = 0.5d0
  heating_rad = 0.0d0
  heating_peak = 1.d16
  heating_sigma = 1.d7
  prob_type = 1

  sponge_weighting = 1.d3
  sponge_start_density = 1.d4
  sponge_width_factor = 10.0d0


  !     Read namelists in probin file
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  !     Open file of initial profile
  open(unit=99,file=model_name)

  ! the model file is assumed to be of the follow form:
  ! # npts = 896
  ! # num of variables = 6
  ! # density
  ! # temperature
  ! # pressure
  ! # carbon-12
  ! # oxygen-16
  ! # magnesium-24
  ! 195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
  ! 585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0

  ! the first line has the number of points in the model
  read (99, '(a256)') header_line

  ipos = index(header_line, '=') + 1
  read (header_line(ipos:),*) npts_model

  print *, npts_model, '    points found in the initial model file'

  ! now read in the number of variables
  read (99, '(a256)') header_line
  ipos = index(header_line, '=') + 1
  read (header_line(ipos:),*) nvars_model_file

  print *, nvars_model_file, 'variables found in the initial model'

  ! now read in the names of the variables
  do i = 1, nvars_model_file
     read (99, '(a256)') header_line
  enddo

  allocate(hse_r(npts_model),hse_rho(npts_model), &
           hse_t(npts_model),hse_p(npts_model))
  allocate(hse_s(nspec,npts_model))

  do k = 1,npts_model
     read(99,*) hse_r(k), hse_rho(k), hse_t(k), hse_p(k), &
                hse_s(1,k), hse_s(2,k), hse_s(3,k)
  end do


  center(1) = 0.0d0

  xmin = problo(1)
  if (xmin /= 0.d0) then
     call bl_error("ERROR: xmin should be 0!")
  endif

  xmax = probhi(1)

  ! store the state at the very top of the model for the boundary
  ! conditions
  allocate (hse_X_top(nspec))

  hse_rho_top  = hse_rho(npts_model)
  hse_t_top    = hse_t(npts_model)
  hse_X_top(:) = hse_s(:,npts_model)

  ! set hse_eint_top and hse_p_top via the EOS
  eos_state%rho   = hse_rho_top
  eos_state%T     = hse_T_top
  eos_state%xn(:) = hse_X_top

  call eos(eos_input_rt, eos_state)

  hse_eint_top = eos_state%e
  hse_p_top = eos_state%p

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
  use eos_data_module
  use eos_type_module
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

  double precision :: x,dist
  integer          :: i,n

  type (eos_t) :: eos_state

  ! Interpolate rho, T and X
  do i = lo(1), hi(1)

     dist = (dble(i) + 0.5d0) * delta(1)

     state(i,URHO ) = interpolate(dist,npts_model,hse_r,hse_rho)
     state(i,UTEMP) = interpolate(dist,npts_model,hse_r,hse_t)

     do n= 1, nspec
        state(i,UFS+n-1) = interpolate(dist,npts_model,hse_r, hse_s(n,:))
     enddo

  enddo

  ! Compute energy from rho,T and X
  do i = lo(1), hi(1)
     eos_state%rho = state(i,URHO)
     eos_state%T = state(i,UTEMP)
     eos_state%xn = state(i,UFS:UFS+nspec-1)

     call eos(eos_input_rt, eos_state)

     ! we'll add the density weighting shortly
     state(i,UEINT) = eos_state%e

  enddo

  do i = lo(1), hi(1)
     state(i,UMX) = 0.d0
     state(i,UEINT) = state(i,URHO) * state(i,UEINT)
     state(i,UEDEN) = state(i,UEINT)
     state(i,UFS:UFS+nspec-1) = state(i,URHO) * state(i,UFS:UFS+nspec-1)
  enddo

end subroutine ca_initdata

! ::: -----------------------------------------------------------

subroutine ca_hypfill(adv,adv_l1,adv_h1, &
                      domlo,domhi,delta,xlo,time,bc)

  use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UTEMP, UEINT, UFS
  use probdata_module, only: hse_rho_top, hse_t_top, hse_X_top, &
       hse_eint_top, hse_p_top
  use network, only: nspec
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

  ! override the generic routine at the top physical boundary
  ! by resetting the velocity to zero there.
  if (adv_h1.gt.domhi(1)) then
     if (bc(1,2,UMX).eq.FOEXTRAP) then
        do i = domhi(1)+1,adv_h1
           !adv(i,UMX) = adv(domhi(1),UMX)
           vel = max(adv(i,UMX)/adv(i,URHO),0.d0)
           adv(i,URHO)  = hse_rho_top
           adv(i,UMX)   = adv(i,URHO)*vel
           adv(i,UTEMP) = hse_T_top
           adv(i,UEINT) = hse_rho_top*hse_eint_top
           adv(i,UEDEN) = hse_rho_top*hse_eint_top + &
                0.5*adv(i,UMX)**2/adv(i,URHO)
           adv(i,UFS:UFS+nspec-1) = hse_rho_top*hse_X_top(:)
        enddo
     end if
  end if

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
