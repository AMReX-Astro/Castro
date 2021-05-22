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

  allocate(dg(3))

  dg(:) = 1

  if (dim .lt. 2) then
     dg(2) = 0
  endif

  if (dim .lt. 3) then
     dg(3) = 0
  endif

end subroutine ca_set_problem_params


! :::
! ::: ----------------------------------------------------------------
! :::



subroutine ca_get_tagging_params(name, namlen) &
     bind(C, name="ca_get_tagging_params")
     ! Initialize the tagging parameters
     !
     ! Binds to C function `ca_get_tagging_params`

  use tagging_module
  use castro_error_module
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /tagging/ &
       denerr, dengrad, dengrad_rel, &
       max_denerr_lev, max_dengrad_lev, max_dengrad_rel_lev, &
       velerr, velgrad, velgrad_rel, &
       max_velerr_lev, max_velgrad_lev, max_velgrad_rel_lev, &
       presserr, pressgrad, pressgrad_rel, &
       max_presserr_lev, max_pressgrad_lev, max_pressgrad_rel_lev, &
       temperr, tempgrad, tempgrad_rel, &
       max_temperr_lev, max_tempgrad_lev, max_tempgrad_rel_lev, &
       raderr, radgrad, radgrad_rel, &
       max_raderr_lev, max_radgrad_lev, max_radgrad_rel_lev, &
       enucerr, max_enucerr_lev, &
       dxnuc_min, dxnuc_max, max_dxnuc_lev

  allocate(denerr)
  allocate(dengrad)
  allocate(dengrad_rel)
  allocate(max_denerr_lev)
  allocate(max_dengrad_lev)
  allocate(max_dengrad_rel_lev)

  ! Set namelist defaults
  denerr = 1.e20_rt
  dengrad = 1.e20_rt
  dengrad_rel = 1.e20_rt
  max_denerr_lev = -1
  max_dengrad_lev = -1
  max_dengrad_rel_lev = -1

  allocate(presserr)
  allocate(pressgrad)
  allocate(pressgrad_rel)
  allocate(max_presserr_lev)
  allocate(max_pressgrad_lev)
  allocate(max_pressgrad_rel_lev)

  presserr = 1.e20_rt
  pressgrad = 1.e20_rt
  pressgrad_rel = 1.e20_rt
  max_presserr_lev = -1
  max_pressgrad_lev = -1
  max_pressgrad_rel_lev = -1

  allocate(velerr)
  allocate(velgrad)
  allocate(velgrad_rel)
  allocate(max_velerr_lev)
  allocate(max_velgrad_lev)
  allocate(max_velgrad_rel_lev)

  velerr  = 1.e20_rt
  velgrad = 1.e20_rt
  velgrad_rel = 1.e20_rt
  max_velerr_lev = -1
  max_velgrad_lev = -1
  max_velgrad_rel_lev = -1

  allocate(temperr)
  allocate(tempgrad)
  allocate(tempgrad_rel)
  allocate(max_temperr_lev)
  allocate(max_tempgrad_lev)
  allocate(max_tempgrad_rel_lev)

  temperr  = 1.e20_rt
  tempgrad = 1.e20_rt
  tempgrad_rel = 1.e20_rt
  max_temperr_lev = -1
  max_tempgrad_lev = -1
  max_tempgrad_rel_lev = -1

  allocate(raderr)
  allocate(radgrad)
  allocate(radgrad_rel)
  allocate(max_raderr_lev)
  allocate(max_radgrad_lev)
  allocate(max_radgrad_rel_lev)

  raderr  = 1.e20_rt
  radgrad = 1.e20_rt
  radgrad_rel = 1.e20_rt
  max_raderr_lev = -1
  max_radgrad_lev = -1
  max_radgrad_rel_lev = -1

  allocate(enucerr)
  allocate(max_enucerr_lev)

  enucerr = 1.e200_rt
  max_enucerr_lev = -1

  allocate(dxnuc_min)
  allocate(dxnuc_max)
  allocate(max_dxnuc_lev)

  dxnuc_min = 1.e200_rt
  dxnuc_max = 1.e200_rt
  max_dxnuc_lev = -1

  ! create the filename
  if (namlen > maxlen) then
     call castro_error('probin file name too long')
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=tagging, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     call castro_error('ERROR: problem in the tagging namelist')
  endif

  close (unit=un)

end subroutine ca_get_tagging_params

subroutine ca_get_ambient_params(name, namlen) bind(C, name="ca_get_ambient_params")
    ! Initialize the ambient parameters

  use ambient_module, only: ambient_state
  use amrex_error_module, only: amrex_error
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR, URHO, UMX, UMZ, UTEMP, UEINT, UEDEN, UFS, &
                                small_dens, small_temp, small_ener
  use amrex_constants_module, only: ZERO
  use actual_network, only: nspec

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  real(rt) :: ambient_density, ambient_temp, ambient_energy

  namelist /ambient/ ambient_density, ambient_temp, ambient_energy

  ! Set namelist defaults

  ambient_density = small_dens
  ambient_temp    = small_temp
  ambient_energy  = small_ener

  ! create the filename
  if (namlen > maxlen) then
     call amrex_error('probin file name too long')
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=ambient, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     call amrex_error('ERROR: problem in the ambient namelist')
  endif

  close (unit=un)

  allocate(ambient_state(NVAR))

  ! Set some initial data in the state for safety, though the
  ! intent is that any problems using this may override these.

  ambient_state(:) = ZERO

  ambient_state(URHO)    = ambient_density
  ambient_state(UMX:UMZ) = ZERO
  ambient_state(UTEMP)   = ambient_temp
  ambient_state(UEINT)   = ambient_density * ambient_energy
  ambient_state(UEDEN)   = ambient_density * ambient_energy
  ambient_state(UFS:UFS+nspec-1) = ambient_density * (1.0e0_rt / nspec)

end subroutine ca_get_ambient_params



subroutine get_ambient_data(ambient_state_out) bind(C, name='get_ambient_data')

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR
  use ambient_module, only: ambient_state

  implicit none

  real(rt), intent(out) :: ambient_state_out(NVAR)

  ambient_state_out(:) = ambient_state(:)

end subroutine get_ambient_data
