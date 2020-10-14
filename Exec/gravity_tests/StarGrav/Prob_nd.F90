subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module
  use probdata_module
  use model_parser_module
  use castro_error_module
  use prob_params_module, only : center

  use amrex_fort_module, only : rt => amrex_real
  implicit none
  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  ! read initial model
  call read_model_file(model_name)

#if AMREX_SPACEDIM == 1
  center(1) = ZERO

#elif AMREX_SPACEDIM == 2
  ! assume axisymmetric
  center(1) = ZERO
  center(2) = HALF*(problo(2)+probhi(2))

#else
  center(1) = HALF*(problo(1)+probhi(1))
  center(2) = HALF*(problo(2)+probhi(2))
  center(3) = HALF*(problo(3)+probhi(3))
#endif

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
                       state, state_lo, state_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module
  use probdata_module
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP,&
                                 UEDEN, UEINT, UFS
  use network, only : nspec
  use model_parser_module
  use prob_params_module, only : center, problo
  use eos_type_module
  use eos_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(rt), intent(in) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1), &
                                   state_lo(2):state_hi(2), &
                                   state_lo(3):state_hi(3), NVAR)

  real(rt) :: xcen, ycen, zcen, dist, pres
  integer :: i, j, k, n

  type(eos_t) :: eos_state

  do k = lo(3), hi(3)
     zcen = problo(3) + delta(3)*(dble(k) + 0.5e0_rt) - center(3)

     do j = lo(2), hi(2)
        ycen = problo(2) + delta(2)*(dble(j) + 0.5e0_rt) - center(2)

        do i = lo(1), hi(1)
           xcen = problo(1) + delta(1)*(dble(i) + 0.5e0_rt) - center(1)

           dist = sqrt(xcen**2 + ycen**2 + zcen**2)

           call interpolate_sub(state(i,j,k,URHO), dist, idens_model)
           call interpolate_sub(state(i,j,k,UTEMP), dist, itemp_model)

           do n = 1, nspec
              call interpolate_sub(state(i,j,k,UFS-1+n), dist, ispec_model-1+n)
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

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state%e
           state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state%e

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           end do

        enddo
     enddo
  enddo

  ! Initial velocities = 0
  state(:,:,:,UMX:UMZ) = ZERO

end subroutine ca_initdata
