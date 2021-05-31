subroutine ca_extern_init(name,namlen) bind(C, name="ca_extern_init")
    ! initialize the external runtime parameters in
    ! extern_probin_module
    !
    ! Binds to C function `ca_extern_init`

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  call runtime_init(name,namlen)

end subroutine ca_extern_init

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine ca_microphysics_init() bind(C, name="ca_microphysics_init")

  use microphysics_module
  use meth_params_module, only: small_dens, small_temp
  implicit none

  call microphysics_init(small_dens=small_dens, small_temp=small_temp)

end subroutine ca_microphysics_init



#ifdef REACTIONS
subroutine ca_set_abort_on_failure(abort_on_failure_in) bind(C, name="ca_set_abort_on_failure")

  use extern_probin_module, only : abort_on_failure

  implicit none

  integer, intent(inout) :: abort_on_failure_in

  if (abort_on_failure_in >= 1) then
     abort_on_failure = .true.
  else
     abort_on_failure = .false.
  endif

end subroutine ca_set_abort_on_failure
#endif

! :::
! ::: ----------------------------------------------------------------
! :::

subroutine ca_set_method_params(dm) &
                                bind(C, name="ca_set_method_params")

  use meth_params_module
  use network, only : nspec, naux
  use amrex_constants_module, only : ZERO, ONE
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: dm

  integer :: ioproc


  !---------------------------------------------------------------------
  ! other initializations
  !---------------------------------------------------------------------

  ! This is a routine which links to the C++ ParallelDescriptor class

  call bl_pd_is_ioproc(ioproc)

  !---------------------------------------------------------------------
  ! safety checks
  !---------------------------------------------------------------------

  if (small_dens <= 0.e0_rt) then
     if (ioproc == 1) then
        call bl_warning("Warning:: small_dens has not been set, defaulting to 1.e-200_rt.")
     endif
     small_dens = 1.e-200_rt
  endif

  if (small_temp <= 0.e0_rt) then
     if (ioproc == 1) then
        call bl_warning("Warning:: small_temp has not been set, defaulting to 1.e-200_rt.")
     endif
     small_temp = 1.e-200_rt
  endif

  if (small_pres <= 0.e0_rt) then
     small_pres = 1.e-200_rt
  endif

  if (small_ener <= 0.e0_rt) then
     small_ener = 1.e-200_rt
  endif

end subroutine ca_set_method_params


! :::
! ::: ----------------------------------------------------------------
! :::



subroutine ca_set_problem_params(dm, &
                                 coord_type_in, &
                                 problo_in, probhi_in) &
                                 bind(C, name="ca_set_problem_params")
     ! Passing data from C++ into f90
     !
     ! Binds to C function `ca_set_problem_params`

  use amrex_constants_module, only: ZERO
  use castro_error_module
  use prob_params_module
  use meth_params_module, only: UMX, UMY, UMZ
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: dm
  integer,  intent(in) :: coord_type_in
  real(rt), intent(in) :: problo_in(dm), probhi_in(dm)

  allocate(dim)

  dim = dm

  allocate(coord_type)
  allocate(problo(3))
  allocate(probhi(3))

  coord_type = coord_type_in

  problo = ZERO
  probhi = ZERO

  problo(1:dm) = problo_in(1:dm)
  probhi(1:dm) = probhi_in(1:dm)

end subroutine ca_set_problem_params
