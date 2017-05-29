subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use bl_types
  use bl_constants_module
  use bl_error_module
  use model_parser_module
  use probdata_module
  use prob_params_module, only: center
  use eos_type_module
  use eos_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(2), probhi(2)

  type (eos_t) :: eos_state

  integer :: untin, i

  namelist /fortin/ model_name, &
       pert_factor, x_pert_loc, pert_width, &
       cutoff_density, refine_cutoff_height, &
       zero_vels

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  if (namlen > maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults here
  zero_vels = .false.
  x_pert_loc = ONE
  pert_width = 0.1_rt
  pert_factor = ONE
  refine_cutoff_height = HALF*(problo(2)+probhi(2))
  
  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(untin)

  ! read the initial model
  call read_model_file(model_name)


  ! set the ambient conditions -- these are used for the upper boundary
  rho_ambient = model_state(npts_model,idens_model)
  T_ambient = model_state(npts_model,itemp_model)
  xn_ambient(:) = model_state(npts_model,ispec_model:ispec_model-1+nspec)

  eos_state%rho = rho_ambient
  eos_state%T = T_ambient
  eos_state%xn(:) = xn_ambient(:)

  call eos(eos_input_rt, eos_state)

  e_ambient = eos_state%e

  ! set center variable in prob_params_module
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))

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
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use bl_constants_module, only: HALF, ZERO, ONE
  use probdata_module
  use model_parser_module
  use interpolate_module
  use prob_params_module, only: problo
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, &
                                 UFS, UTEMP
  use eos_module
  use eos_type_module
  use network, only: nspec, network_species_index

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1, state_l2, state_h1, state_h2
  real(rt)         :: xlo(2), xhi(2), time, delta(2)
  real(rt)         :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  integer :: i, j, n
  real(rt)         :: x, y

  real(rt)         :: dens, temp, pres

  type (eos_t) :: eos_state
  
  integer :: ifuel

  ifuel = network_species_index("fuel")

  ! initialize from the model file and add an isobaric perturbation
  do j=lo(2),hi(2)
     y = problo(2) + (dble(j)+HALF)*delta(2)

     do i=lo(1),hi(1)
        x = problo(1) + (dble(i)+HALF)*delta(1)

        temp = interpolate(y,npts_model,model_r, &
                           model_state(:,itemp_model))

        dens = interpolate(y,npts_model,model_r, &
                           model_state(:,idens_model))

        pres = interpolate(y,npts_model,model_r, &
                           model_state(:,ipres_model))

        do n = 1, nspec
           state(i,j,UFS-1+n) = &
                interpolate(y,npts_model,model_r, model_state(:,ispec_model-1+n))
        enddo

        if (dens > cutoff_density .and. state(i,j,UFS-1+ifuel) > 0.99e0_rt) then
           state(i,j,UTEMP) = temp * (ONE + (pert_factor * &
                (ONE + tanh((x_pert_loc-x)/pert_width)) ) )
        else
           state(i,j,UTEMP) = temp
        endif

        eos_state%T = state(i,j,UTEMP)
        eos_state%rho = dens
        eos_state%p = pres
        eos_state%xn(:) = state(i,j,UFS:UFS-1+nspec)

        call eos(eos_input_tp, eos_state)

        state(i,j,URHO) = eos_state%rho
        state(i,j,UEINT) = eos_state%e

        ! make state conservative
        state(i,j,UFS:UFS-1+nspec) = state(i,j,UFS:UFS-1+nspec)*state(i,j,URHO)
        state(i,j,UEINT) = state(i,j,UEINT)*state(i,j,URHO)

        ! assumes ke=0
        state(i,j,UEDEN) = state(i,j,UEINT)

        state(i,j,UMX:UMZ) = ZERO

     end do
  end do

end subroutine ca_initdata

