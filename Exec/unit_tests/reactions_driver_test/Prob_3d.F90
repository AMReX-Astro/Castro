subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use model_parser_module, only: read_model_file
  use probdata_module
  implicit none

  integer          :: init, namlen
  integer          :: name(namlen)
  real(rt)         :: problo(3), probhi(3)

  integer :: untin,i

  namelist /fortin/ model_name

  ! Build "probin" filename -- the name of file containing fortin namelist.
  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  if (namlen > maxlen) call amrex_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(untin)

  ! read the initial model
  call read_model_file(model_name)

end subroutine amrex_probinit


module initdata_module

contains

  AMREX_DEVICE subroutine ca_initdata(lo, hi, &
                                      state, state_lo, state_hi, &
                                      dx, problo, domlo) bind(c, name='ca_initdata')

    use amrex_error_module

    use amrex_fort_module, only : rt => amrex_real
    use amrex_constants_module, only: ZERO
    use network, only: nspec
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rt
    use meth_params_module, only: URHO, UTEMP, UFS, UEDEN, UEINT, NVAR
    use model_parser_module, only: idens_model, itemp_model, ispec_model, npts_model, model_state, model_initialized
    implicit none

    integer :: lo(3), hi(3)
    integer :: state_lo(3), state_hi(3)
    real(rt)         :: problo(3), domlo(3), dx(3)
    real(rt)         :: state(state_lo(1):state_hi(1), &
                              state_lo(2):state_hi(2), &
                              state_lo(3):state_hi(3), NVAR)
    integer :: nzones_state, i, j, k, n, ii
    type (eos_t) :: eos_state

    ! Check to make sure model is initialized
    if (.not. model_initialized) then
#if !(defined(CUDA) || defined(ACC))
       call amrex_error("Model has not been initialized")
#endif
    endif

    ! Zero the state
    state(:,:,:,:) = ZERO

    ! Check to make sure the number of zones in state and model are the same
    nzones_state = (hi(1)-lo(1)+1) * (hi(2)-lo(2)+1) * (hi(3)-lo(3)+1)
    if (.not. nzones_state .eq. npts_model) then
#if !(defined(CUDA) || defined(ACC))
       write(*,*) "nzones_state: ", nzones_state
       write(*,*) "npts_model: ", npts_model
       call amrex_error("Number of zones in state not equal to number of zones in model! Make sure you are running with one box.")
#endif
    endif

    ! Fill state with model data
    n = 0
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             n = n + 1
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

end module initdata_module
