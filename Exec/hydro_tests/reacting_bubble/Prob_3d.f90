subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use model_parser_module
  use amrex_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(3), probhi(3)

  integer untin,i

  namelist /fortin/ model_name, pert_temp_factor, pert_rad_factor

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  ! Build "probin" filename from C++ land --
  ! the name of file containing fortin namelist.

  if (namlen .gt. maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! Namelist defaults

  ! Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)


  ! Read initial model
  call read_model_file(model_name)

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
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       delta,xlo,xhi)

  use probdata_module
  use interpolate_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt, eos_input_tp
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module
  use amrex_constants_module, only : ONE, HALF, TWO

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(3), hi(3)
  integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real(rt)         xlo(3), xhi(3), time, delta(3)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2, &
                         state_l3:state_h3,NVAR)

  real(rt)         dist, x, y, z
  integer i, j, k, n

  real(rt)         t0,x1,y1,z1,r1,x2,y2,z2,r2,x3,y3,z3,r3,x4,y4,r4,temp

  real(rt)         temppres(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3)

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

        do i = lo(1), hi(1)

           state(i,j,k,URHO)  = interpolate(z,npts_model,model_r, &
                                            model_state(:,idens_model))
           state(i,j,k,UTEMP) = interpolate(z,npts_model,model_r, &
                                            model_state(:,itemp_model))
           do n = 1, nspec
              state(i,j,k,UFS-1+n) = interpolate(z,npts_model,model_r, &
                                                 model_state(:,ispec_model-1+n))
           enddo

        enddo
     enddo
  enddo

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = eos_state%e
           temppres(i,j,k) = eos_state%p
        enddo
     enddo
  enddo

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           state(i,j,k,UEDEN) = state(i,j,k,URHO) * state(i,j,k,UEINT)
           state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           enddo

        enddo
     enddo
  enddo

  ! Initial velocities = 0
  state(:,:,:,UMX:UMZ) = 0.e0_rt

  ! Now add the perturbation
  do k = lo(3), hi(3)
     z = xlo(3) + delta(3)*(dble(k-lo(3)) + HALF)

     do j = lo(2), hi(2)
        y = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)

        do i = lo(1), hi(1)
           x = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)

           t0 = state(i,j,k,UTEMP)

           x1 = 5.0e7_rt
           y1 = 5.0e7_rt
           z1 = 6.5e7_rt
           r1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 ) / (2.5e6_rt*pert_rad_factor)

           x2 = 1.2e8_rt
           y2 = 1.2e8_rt
           z2 = 8.5e7_rt
           r2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2 ) / (2.5e6_rt*pert_rad_factor)

           x3 = 2.0e8_rt
           y3 = 2.0e8_rt
           z3 = 7.5e7_rt
           r3 = sqrt( (x-x3)**2 + (y-y3)**2 + (z-z3)**2 ) / (2.5e6_rt*pert_rad_factor)

           state(i,j,k,UTEMP) = t0 * (ONE + pert_temp_factor* &
                (0.150e0_rt * (ONE + tanh(TWO-r1)) + &
                 0.300e0_rt * (ONE + tanh(TWO-r2)) + &
                 0.225e0_rt * (ONE + tanh(TWO-r3))))

           state(i,j,k,UEINT) = state(i,j,k,UEINT) / state(i,j,k,URHO)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) =  state(i,j,k,UFS+n-1) / state(i,j,k,URHO)
           enddo

           eos_state%T = state(i,j,k,UTEMP)
           eos_state%p = temppres(i,j,k)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_tp, eos_state)

           state(i,j,k,UEINT) = eos_state%e
           state(i,j,k,URHO) = eos_state%rho

           state(i,j,k,UEDEN) = state(i,j,k,UEINT)*state(i,j,k,URHO)
           state(i,j,k,UEINT) = state(i,j,k,UEINT)*state(i,j,k,URHO)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           enddo

        enddo
     enddo
  enddo

end subroutine ca_initdata
