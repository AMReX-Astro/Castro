subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use model_parser_module
  use castro_error_module
  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  integer :: untin, i
  real(rt) :: offset

  namelist /fortin/ model_name, apply_vel_field, &
                    velpert_scale, velpert_amplitude, velpert_height_loc, num_vortices

  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  ! Build "probin" filename from C++ land --
  ! the name of file containing fortin namelist.


  if (namlen .gt. maxlen) call castro_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! Namelist defaults
  apply_vel_field = .false.
  velpert_scale = 1.0e2_rt
  velpert_amplitude = 1.0e2_rt
  velpert_height_loc = 6.5e3_rt
  num_vortices = 1

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin, fortin)
  close(unit=untin)

  ! Read initial model
  call read_model_file(model_name)

  if (parallel_IOProcessor()) then
     do i = 1, npts_model
        print *, i, model_r(i), model_state(i,idens_model)
     enddo
  endif

  ! velocity perturbation stuff
  offset = (probhi(1) - problo(1)) / (num_vortices)

  allocate(xloc_vortices(num_vortices))

  do i = 1, num_vortices
     xloc_vortices(i) = (dble(i-1) + 0.5e0_rt) * offset + problo(1)
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
subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module, only : ZERO, HALF, ONE
  use prob_params_module, only : problo
  use probdata_module
  use eos_module
  use eos_type_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module
  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xdist, ydist, x, y, r, upert(2)
  integer i, j, k, n, vortex

  real(rt) :: temppres(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3))

  type (eos_t) :: eos_state

#if AMREX_SPACEDIM == 3
  call castro_error("Error: 3-d initialization not implemented")
#endif

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)

           call interpolate_sub(state(i,j,k,URHO), y, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), y, itemp_model)
           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), y, ispec_model-1+n)
           end do

           eos_state%rho = state(i,j,k,URHO)
           eos_state%T = state(i,j,k,UTEMP)
           eos_state%xn(:) = state(i,j,k,UFS:UFS-1+nspec)

           call eos(eos_input_rt, eos_state)

           state(i,j,k,UEINT) = eos_state%e
           temppres(i,j,k) = eos_state%p

           state(i,j,k,UEDEN) = state(i,j,k,URHO) * state(i,j,k,UEINT)
           state(i,j,k,UEINT) = state(i,j,k,URHO) * state(i,j,k,UEINT)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

           state(i,j,k,UMX:UMZ) = ZERO

        end do
     end do
  end do


  ! Now add the velocity perturbation
  if (apply_vel_field) then

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           y = problo(2) + delta(2)*(dble(j) + HALF)
           ydist = y - velpert_height_loc

           do i = lo(1), hi(1)
              x = problo(1) + delta(1)*(dble(i) + HALF)

              upert = ZERO

              ! loop over each vortex
              do vortex = 1, num_vortices

                 xdist = x - xloc_vortices(vortex)

                 r = sqrt(xdist**2 + ydist**2)

                 upert(1) = upert(1) - (ydist/velpert_scale) * &
                      velpert_amplitude * exp( -r**2/(2.e0_rt*velpert_scale**2)) &
                      * (-ONE)**vortex

                 upert(2) = upert(2) + (xdist/velpert_scale) * &
                      velpert_amplitude * exp(-r**2/(2.e0_rt*velpert_scale**2)) &
                      * (-ONE)**vortex
              end do

              state(i,j,k,UMX) = state(i,j,k,URHO) * upert(1)
              state(i,j,k,UMY) = state(i,j,k,URHO) * upert(2)

           end do
        end do
     end do
  end if

end subroutine ca_initdata
