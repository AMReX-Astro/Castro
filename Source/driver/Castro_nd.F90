#ifdef MICROPHYSICS_FORT
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
#endif


#ifdef REACTIONS
subroutine ca_set_abort_on_failure(abort_on_failure_in) bind(C, name="ca_set_abort_on_failure")

#ifdef MICROPHYSICS_FORT
  use extern_probin_module, only : abort_on_failure
#endif

  implicit none

  integer, intent(inout) :: abort_on_failure_in

#ifdef MICROPHYSICS_FORT
  if (abort_on_failure_in >= 1) then
     abort_on_failure = .true.
  else
     abort_on_failure = .false.
  endif
#endif

end subroutine ca_set_abort_on_failure
#endif

