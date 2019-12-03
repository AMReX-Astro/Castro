subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use castro_error_module
  use amrex_fort_module, only : rt => amrex_real
  use model_parser_module
  use probdata_module
  use burner_module
  use burn_type_module

  implicit none

  integer          :: init, namlen
  integer          :: name(namlen)
  real(rt)         :: problo(3), probhi(3)

  type(burn_t) :: burn_state_in, burn_state_out

  integer :: k, n

  call probdata_init(name, namlen)

  ! read the initial model
  call read_model_file(model_name)

  ! do the burning
  do k = 1, npts_model

     burn_state_in % rho = model_state(k, idens_model)
     burn_state_in % T = model_state(k, itemp_model)
     burn_state_in % e = ZERO

     do n = 1, nspec
        burn_state_in % xn(n) = model_state(k, ispec_model-1+n)
     enddo

     burn_state_in % i = k
     burn_state_in % j = -1
     burn_state_in % k = -1

     burn_state_in % success = .true.

     call burner(burn_state_in, burn_state_out, burn_dt, ZERO)

     if (.not. burn_state_out % success) then
        print *, "burning failed for zone ", k
     else
        print *, "zone, enuc = ", k, burn_state_out % e / burn_dt
     end if

  end do

end subroutine amrex_probinit


subroutine ca_initdata(level,time,lo,hi,nvar, &
     state,state_lo,state_hi, &
     dx,xlo,xhi) 

  use castro_error_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only: ZERO
  use network, only: nspec
  use eos_module, only : eos
  use eos_type_module, only : eos_t, eos_input_rt
  use meth_params_module, only: URHO, UTEMP, UFS, UEDEN, UEINT
  use model_parser_module, only: idens_model, itemp_model, ispec_model, npts_model, model_state, model_initialized
  implicit none

  integer :: level, nvar
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  real(rt)         :: xlo(3), xhi(3), dx(3), time
  real(rt)         :: state(state_lo(1):state_hi(1), &
       state_lo(2):state_hi(2), &
       state_lo(3):state_hi(3), NVAR)
  integer :: nzones_state, i, j, k, n, ii, m
  type (eos_t) :: eos_state

  ! Check to make sure model is initialized
  if (.not. model_initialized) then
#if !(defined(CUDA) || defined(ACC))
     call castro_error("Model has not been initialized")
#endif
  endif

  ! Zero the state
  state(:,:,:,:) = ZERO

  ! Fill state with model data
  m = 0
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           m = m + 1
           n = mod(m, npts_model) + 1
           
           state(i, j, k, URHO) = model_state(n, idens_model)
           state(i, j, k, UTEMP) = model_state(n, itemp_model)
           state(i, j, k, UFS:UFS-1+nspec) = model_state(n, ispec_model:ispec_model-1+nspec)

           ! Set the energy via an EOS call.
           eos_state % rho = state(i, j, k, URHO)
           eos_state % T   = state(i, j, k, UTEMP)
           eos_state % xn  = state(i, j, k, UFS:UFS-1+nspec)

           call eos(eos_input_rt, eos_state)

           state(i, j, k, UEINT) = state(i, j, k, URHO) * eos_state % e
           state(i, j, k, UEDEN) = state(i, j, k, URHO) * eos_state % e
           do ii = UFS, UFS-1+nspec
              state(i, j, k, ii) = state(i, j, k, URHO) * state(i, j, k, ii)
           enddo
        enddo
     enddo
  enddo

end subroutine ca_initdata

