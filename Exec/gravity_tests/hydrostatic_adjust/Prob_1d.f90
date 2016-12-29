subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_error_module
  use probdata_module
  use prob_params_module, only: center
  use eos_module
  use eos_type_module
  use model_parser_module

  use network, only : nspec
  implicit none

  integer :: init,namlen,untin,i,k
  integer :: name(namlen)

  double precision problo(1), probhi(1)

  type (eos_t) :: eos_state

  namelist /fortin/ &
       model_name,  &
       heating_time, heating_rad, heating_peak, heating_sigma, &
       prob_type


  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !
  integer   :: ipos
  integer, parameter :: maxlen=127
  character probin*(maxlen)
  character (len=256) :: header_line

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


  !     Read namelists in probin file
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  ! Read in the initial model

  call read_model_file(model_name)

  ! Save some of the data locally

  allocate(hse_r(npts_model),hse_rho(npts_model), &
           hse_t(npts_model),hse_p(npts_model))
  allocate(hse_s(nspec,npts_model))

  hse_r   = model_r(:)
  hse_rho = model_state(:,idens_model)
  hse_t   = model_state(:,itemp_model)
  hse_p   = model_state(:,ipres_model)
  do i = 1, nspec
     hse_s(i,:) = model_state(:,ispec_model+i-1)
  enddo

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

end subroutine amrex_probinit

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
  use eos_type_module
  use network, only : nspec
  use interpolate_module
  use model_parser_module, only: npts_model
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
     state(i,UMX:UMZ) = 0.d0
     state(i,UEINT) = state(i,URHO) * state(i,UEINT)
     state(i,UEDEN) = state(i,UEINT)
     state(i,UFS:UFS+nspec-1) = state(i,URHO) * state(i,UFS:UFS+nspec-1)
  enddo

end subroutine ca_initdata

