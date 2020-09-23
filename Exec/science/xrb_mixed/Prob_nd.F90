subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use model_parser_module
  use castro_error_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  real(rt) :: offset

  integer :: i

  if (num_vortices > max_num_vortices) then
     call castro_error("num_vortices too large, please increase max_num_vortices and the size of xloc_vortices")
  end if

  ! Read initial model
  call read_model_file(model_name)

  do i = 1, npts_model
     print *, i, model_r(i), model_state(i,idens_model)
  enddo

  ! velocity perturbation stuff
  offset = (probhi(1) - problo(1)) / (num_vortices + 1)

  allocate(xloc_vortices(num_vortices))

  do i = 1, num_vortices
     xloc_vortices(i) = dble(i) * offset
  enddo

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
subroutine ca_initdata(level,time, lo, hi, nscal, &
                       state, state_lo, state_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module
  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module
  use amrex_constants_module, only : ZERO, HALF
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : problo

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3),NVAR)

  real(rt) :: xdist, ydist, zdist, x, y, z, r, upert(3), height
  integer i, j, k, n, vortex

  real(rt) :: temppres(state_lo(1):state_hi(1), &
                       state_lo(2):state_hi(2), &
                       state_lo(3):state_hi(3))

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 2
           height = y
#elif AMREX_SPACEDIM == 3
           height = z
#else
           call castro_error("invalid dimensionality")
#endif

           call interpolate_sub(state(i,j,k,URHO), height, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), height, itemp_model)
           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), height, ispec_model-1+n)
           end do

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rt, eos_state)

           temppres(i,j,k) = eos_state % p

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state % e

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        end do
     end do
  end do

  ! Initial velocities = 0
  state(:,:,:,UMX:UMZ) = 0.e0_rt

  ! Now add the velocity perturbation
  if (apply_vel_field) then

     do k = lo(3), hi(3)
        z = problo(3) + delta(3)*(dble(k) + HALF)
        zdist = z - velpert_height_loc

        do j = lo(2), hi(2)
           y = problo(2) + delta(2)*(dble(j) + HALF)
           ydist = y - velpert_height_loc

           do i = lo(1), hi(1)
              x = problo(1) + delta(1)*(dble(i) + HALF)

              upert = ZERO

              ! loop over each vortex
              do vortex = 1, num_vortices

                 xdist = x - xloc_vortices(vortex)

#if AMREX_SPACEDIM == 2
                 r = sqrt(xdist**2 + ydist**2)

                 upert(1) = upert(1) - ydist * &
                      velpert_amplitude * exp( -r**2/(TWO*velpert_scale)) &
                      * (-ONE)**vortex

                 upert(2) = upert(2) + xdist * &
                      velpert_amplitude * exp(-r**2/(TWO*velpert_scale)) &
                      * (-ONE)**vortex

                 upert(3) = ZERO
#else
                 r = sqrt(xdist**2 + zdist**2)

                 upert(1) = upert(1) - zdist * &
                      velpert_amplitude * exp( -r**2/(TWO*velpert_scale)) &
                      * (-ONE)**vortex

                 upert(3) = upert(3) + xdist * &
                      velpert_amplitude * exp(-r**2/(TWO*velpert_scale)) &
                      * (-ONE)**vortex

                 upert(2) = ZERO
#endif

              end do

              state(i,j,k,UMX) = state(i,j,k,URHO) * upert(1)
              state(i,j,k,UMY) = state(i,j,k,URHO) * upert(2)
              state(i,j,k,UMZ) = state(i,j,k,URHO) * upert(3)

           end do
        end do
     end do
  endif

end subroutine ca_initdata
