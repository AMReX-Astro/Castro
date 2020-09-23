subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use probdata_module
  use castro_error_module
  use initial_model_module, only : generate_initial_model, model_t
  use network, only : nspec, network_species_index
  use amrex_constants_module, only : ONE, HALF, ZERO
  use extern_probin_module, only : small_x
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type(model_t) :: model_params

  integer :: ihe4
  integer :: nbuf

  ihe4 = network_species_index("helium-4")
  if (ihe4 < 0) then
     call castro_error("Error: helium-4 not present")
  end if

  model_params % T_base = temp_base
  model_params % dens_base = dens_base
  model_params % xn(:) = 100*small_x
  model_params % xn(ihe4) = ONE - (nspec - 1) * 100*small_x

  ! we add some buffer to the model so we can use it to fill ghost cells in the boundary conditions
  nbuf = 8
  call generate_initial_model(nx_model, problo(AMREX_SPACEDIM), probhi(AMREX_SPACEDIM), model_params, nbuf)

  center(:) = HALF * (problo(:) + probhi(:))

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

  use probdata_module
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt, eos_input_tp
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use network, only: nspec
  use model_parser_module
  use prob_params_module, only : problo
  use amrex_constants_module, only : ZERO, ONE, HALF, TWO, FOUR
  use amrex_fort_module, only : rt => amrex_real
  use prob_params_module, only : center

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: dist, x, y, z, height, r
  integer :: i, j, k, n

  real(rt) :: t0

  type (eos_t) :: eos_state

  do k = lo(3), hi(3)
     z = problo(3) + delta(3)*(dble(k) + HALF)

     do j = lo(2), hi(2)
        y = problo(2) + delta(2)*(dble(j) + HALF)

        do i = lo(1), hi(1)
           x = problo(1) + delta(1)*(dble(i) + HALF)

#if AMREX_SPACEDIM == 1
           height = x
#elif AMREX_SPACEDIM == 2
           height = y
#else
           height = z
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

           state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state % e
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           do n = 1,nspec
              state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * state(i,j,k,UFS+n-1)
           enddo

           ! Initial velocities = 0
           state(i,j,k,UMX:UMZ) = ZERO

        end do
     end do
  end do

end subroutine ca_initdata
