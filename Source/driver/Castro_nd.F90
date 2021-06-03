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

subroutine ca_microphysics_init(small_dens_in, small_temp_in) bind(C, name="ca_microphysics_init")

  use microphysics_module
  implicit none

  real (kind=rt), intent(in), value :: small_dens_in, small_temp_in

  call microphysics_init(small_dens=small_dens_in, small_temp=small_temp_in)

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

