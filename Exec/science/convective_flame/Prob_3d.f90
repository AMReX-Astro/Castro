subroutine PROBINIT (init,name,namlen,problo,probhi)

  use bl_types
  use bl_constants_module
  use bl_error_module
  use model_parser_module
  use probdata_module
  use prob_params_module, only: center

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: init, namlen
  integer :: name(namlen)
  real(rt)         :: problo(3), probhi(3)

  integer :: untin,i

  namelist /fortin/ model_name, &
       pert_factor, x_pert_loc, pert_width, &
       cutoff_density, &
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

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(untin)

  ! read the initial model
  call read_model_file(model_name)

  ! set center variable in prob_params_module
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))
  center(3) = HALF*(problo(3)+probhi(3))

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
  use model_parser_module
  use interpolate_module
  use prob_params_module, only: problo
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, &
                                 UFS, UTEMP
  use eos_module
  use eos_type_module
  use network, only: nspec, network_species_index

  use bl_fort_module, only : rt => c_real
  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt)         xlo(3), xhi(3), time, delta(3)
  real(rt)         state(state_l1:state_h1, &
                         state_l2:state_h2, &
                         state_l3:state_h3,NVAR)

  integer :: i, j, k, n
  real(rt)         :: x, y, z
  real(rt)         :: dens, temp, pres

  type (eos_t) :: eos_state

  integer :: ifuel

  ifuel = network_species_index("fuel")

  do k = lo(3), hi(3)
     z = problo(3) + (dble(k)+HALF)*delta(3)

     do j = lo(2), hi(2)
        y = problo(2) + (dble(j)+HALF)*delta(2)

        do i = lo(1), hi(1)
           x = problo(1) + (dble(i)+HALF)*delta(1)

           temp = interpolate(z,npts_model,model_r, &
                              model_state(:,itemp_model))

           dens = interpolate(z,npts_model,model_r, &
                              model_state(:,idens_model))

           pres = interpolate(z,npts_model,model_r, &
                              model_state(:,ipres_model))

           do n = 1, nspec
              state(i,j,k,UFS-1+n) = &
                   interpolate(z,npts_model,model_r, model_state(:,ispec_model-1+n))
           enddo

           if (dens > cutoff_density .and. state(i,j,k,UFS-1+ifuel) > 0.99e0_rt) then
              state(i,j,k,UTEMP) = temp * (ONE + (pert_factor * &
                   (ONE + tanh((x_pert_loc-x)/pert_width)) ) )
           else
              state(i,j,k,UTEMP) = temp
           endif

           eos_state%T = state(i,j,k,UTEMP)
           eos_state%rho = dens
           eos_state%p = pres
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)
           
           call eos(eos_input_tp, eos_state)

           state(i,j,k,URHO) = eos_state%rho
           state(i,j,k,UEINT) = eos_state%e

           ! make state conservative
           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,UFS:UFS-1+nspec)*state(i,j,k,URHO)
           state(i,j,k,UEINT) = state(i,j,k,UEINT)*state(i,j,k,URHO)

           ! assumes ke=0
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)
           
           state(i,j,k,UMX:UMZ) = ZERO

        enddo
     enddo
  enddo
  
end subroutine ca_initdata

