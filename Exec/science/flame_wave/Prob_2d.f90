subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)

  use probdata_module
  use model_parser_module
  use bl_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)

  integer untin,i

  namelist /fortin/ model_name, interp_BC, zero_vels, &
                    dtemp, x_half_max, x_half_width, &
                    X_min, cutoff_density

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  ! Build "probin" filename from C++ land --
  ! the name of file containing fortin namelist.

  if (namlen > maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Namelist defaults
  X_min = 1.e-4_rt
  cutoff_density = 500.e0_rt


  dtemp = 3.81e8_rt
  x_half_max = 1.2e5_rt
  x_half_width = 3.6e4_rt

  interp_BC = .false.
  zero_vels = .false.

  open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
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
! ::: -----------------------------------------------------------
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, state_l1, state_l2, state_h1, state_h2, &
                       delta, xlo, xhi)

  use bl_constants_module
  use probdata_module
  use interpolate_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt, eos_input_tp
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only: problo
  use network, only: nspec
  use model_parser_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt)         dist,x,y
  integer i,j,n

  real(rt)         t0,x1,y1,r1,temp

  real(rt)         temppres(state_l1:state_h1,state_l2:state_h2)

  type (eos_t) :: eos_state


  do j = lo(2), hi(2)
     y = problo(2) + (dble(j)+HALF)*delta(2)

     do i = lo(1), hi(1)

        state(i,j,URHO)  = interpolate(y,npts_model,model_r, &
                                       model_state(:,idens_model))
        state(i,j,UTEMP) = interpolate(y,npts_model,model_r, &
                                       model_state(:,itemp_model))

        state(i,j,UFS:UFS-1+nspec) = ZERO

        do n = 1, nspec
           state(i,j,UFS-1+n) = interpolate(y,npts_model,model_r, &
                                            model_state(:,ispec_model-1+n))
        enddo

     enddo
  enddo

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        eos_state%rho = state(i,j,URHO)
        eos_state%T = state(i,j,UTEMP)
        eos_state%xn(:) = state(i,j,UFS:UFS-1+nspec)

        call eos(eos_input_rt, eos_state)

        state(i,j,UEINT) = eos_state%e
        temppres(i,j) = eos_state%p

     end do
  end do

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        state(i,j,UEDEN) = state(i,j,URHO) * state(i,j,UEINT)
        state(i,j,UEINT) = state(i,j,URHO) * state(i,j,UEINT)

        do n = 1, nspec
           state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
        end do

     enddo
  enddo

  ! Initial velocities = 0
  state(:,:,UMX:UMZ) = 0.e0_rt

  ! Now add the perturbation
  do j = lo(2), hi(2)
     y = problo(2) + (dble(j)+HALF)*delta(2)

     do i = lo(1), hi(1)
        x = problo(1) + (dble(i)+HALF)*delta(1)

        if (state(i,j,UFS) > 0.1 .and. state(i,j,URHO) > 1.0e5_rt) then
           state(i,j,UTEMP)=state(i,j,UTEMP) + dtemp / &
                (ONE + exp((x-x_half_max)/x_half_width))
        end if

        do n = 1,nspec
           state(i,j,UFS+n-1) = state(i,j,UFS+n-1) / state(i,j,URHO)
        end do

        eos_state%T = state(i,j,UTEMP)
        eos_state%p = temppres(i,j)
        eos_state%xn(:) = state(i,j,UFS:UFS-1+nspec)

        call eos(eos_input_tp, eos_state)

        state(i,j,UEINT) = eos_state%e
        state(i,j,URHO) = eos_state%rho

        state(i,j,UEDEN) = state(i,j,UEINT)*state(i,j,URHO)
        state(i,j,UEINT) = state(i,j,UEINT)*state(i,j,URHO)

        do n = 1,nspec
           state(i,j,UFS+n-1) = state(i,j,URHO) * state(i,j,UFS+n-1)
        end do

     end do
  end do

end subroutine ca_initdata
