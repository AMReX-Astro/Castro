subroutine ca_microphysics_init(small_dens_in, small_temp_in) bind(C, name="ca_microphysics_init")

  use amrex_fort_module, only: rt => amrex_real
#ifdef MICROPHYSICS_FORT
  use microphysics_module
#endif

  implicit none

  real (kind=rt), intent(in), value :: small_dens_in, small_temp_in

#ifdef MICROPHYSICS_FORT
  call microphysics_init(small_dens=small_dens_in, small_temp=small_temp_in)
#endif

end subroutine ca_microphysics_init



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

