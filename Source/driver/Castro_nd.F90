

subroutine ca_network_init() bind(C, name="ca_network_init")
    !
    ! Binds to C function `ca_network_init`
    !

  use network, only: network_init
#ifdef REACTIONS
  use actual_rhs_module, only: actual_rhs_init
#endif

  call network_init()

#ifdef REACTIONS
  call actual_rhs_init()
#endif

end subroutine ca_network_init


subroutine ca_network_finalize() bind(C, name="ca_network_finalize")
    !
    ! Binds to C function `ca_network_finalize`
    !

  use network, only: network_finalize

  call network_finalize()

end subroutine ca_network_finalize


subroutine ca_eos_init() bind(C, name="ca_eos_init")

  use meth_params_module, only: small_dens, small_temp
  use eos_module, only: eos_init

  implicit none

  call eos_init(small_dens=small_dens, small_temp=small_temp)

  !$acc update device(small_dens, small_temp)

end subroutine ca_eos_init


subroutine ca_eos_finalize() bind(C, name="ca_eos_finalize")
    !
    ! Binds to C function `ca_eos_finalize`
    !

  use eos_module, only: eos_finalize

  call eos_finalize()

end subroutine ca_eos_finalize


! :::
! ::: ----------------------------------------------------------------
! :::



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



subroutine ca_amrinfo_init() bind(C, name="ca_amrinfo_init")
    !
    ! Binds to C function `ca_amrinfo_init`
    !

  use amrinfo_module, only: amr_level, amr_iteration, amr_ncycle, amr_time, amr_dt
  use amrex_fort_module, only: rt => amrex_real

  allocate(amr_level)
  amr_level = 0
  allocate(amr_iteration)
  amr_iteration = 0
  allocate(amr_ncycle)
  amr_ncycle = 0
  allocate(amr_time)
  amr_time = 0.0e0_rt
  allocate(amr_dt)
  amr_dt = 0.0e0_rt

end subroutine ca_amrinfo_init



subroutine ca_amrinfo_finalize() bind(C, name="ca_amrinfo_finalize")
    !
    ! Binds to C function `ca_amrinfo_finalize`
    !

  use amrinfo_module, only: amr_level, amr_iteration, amr_ncycle, amr_time, amr_dt

  deallocate(amr_level)
  deallocate(amr_iteration)
  deallocate(amr_ncycle)
  deallocate(amr_time)
  deallocate(amr_dt)

end subroutine ca_amrinfo_finalize




subroutine ca_set_amr_info(level_in, iteration_in, ncycle_in, time_in, dt_in) &
     bind(C, name="ca_set_amr_info")
     !
     ! Binds to C function `ca_set_amr_info`
     !
     ! level_in integer
     ! iteration_in integer
     ! ncycle_in integer
     ! time_in real(rt)
     ! dt_in real(rt)
     !

  use amrinfo_module, only: amr_level, amr_iteration, amr_ncycle, amr_time, amr_dt
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: level_in, iteration_in, ncycle_in
  real(rt), intent(in) :: time_in, dt_in

  if (level_in .ge. 0) then
     amr_level = level_in
  endif

  if (iteration_in .ge. 0) then
     amr_iteration = iteration_in
  endif

  if (ncycle_in .ge. 0) then
     amr_ncycle = ncycle_in
  endif

  if (time_in .ge. 0.0e0_rt) then
     amr_time = time_in
  endif

  if (dt_in .ge. 0.0e0_rt) then
     amr_dt = dt_in
  endif


end subroutine ca_set_amr_info


! :::
! ::: ----------------------------------------------------------------
! :::



subroutine allocate_outflow_data(np,nc) &
     bind(C, name="allocate_outflow_data")
     !
     ! Binds to C function `allocate_outflow_data`

  use meth_params_module, only: outflow_data_old, outflow_data_new, outflow_data_allocated
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: np,nc

  if (.not. outflow_data_allocated) then
     allocate(outflow_data_old(nc,np))
     allocate(outflow_data_new(nc,np))
  end if

  outflow_data_allocated = .true.

end subroutine allocate_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::


subroutine set_old_outflow_data(radial,time,np,nc) &
     bind(C, name="set_old_outflow_data")
     ! Passing data from C++ to f90
     !
     ! Binds to C function `set_old_outflow_data`

  use meth_params_module, only: outflow_data_old, outflow_data_old_time
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: radial(nc,np)
  real(rt), intent(in) :: time
  integer,  intent(in) :: np,nc

  ! Do this so the routine has the right size
  deallocate(outflow_data_old)
  allocate(outflow_data_old(nc,np))

  outflow_data_old(1:nc,1:np) = radial(1:nc,1:np)

  outflow_data_old_time = time

end subroutine set_old_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::


subroutine set_new_outflow_data(radial,time,np,nc) &
     bind(C, name="set_new_outflow_data")
     ! Passing data from C++ to f90
     !
     ! Binds to C function `set_new_outflow_data`

  use meth_params_module, only: outflow_data_new, outflow_data_new_time
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: radial(nc,np)
  real(rt), intent(in) :: time
  integer,  intent(in) :: np,nc

  ! Do this so the routine has the right size
  deallocate(outflow_data_new)
  allocate(outflow_data_new(nc,np))

  outflow_data_new(1:nc,1:np) = radial(1:nc,1:np)

  outflow_data_new_time = time

end subroutine set_new_outflow_data

! :::
! ::: ----------------------------------------------------------------
! :::


subroutine swap_outflow_data() bind(C, name="swap_outflow_data")
    !
    ! Binds to C function `swap_outflow_data`
    !

  use meth_params_module, only: outflow_data_new, outflow_data_new_time, &
       outflow_data_old, outflow_data_old_time
  use castro_error_module
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer :: np, nc

  nc = size(outflow_data_new,dim=1)
  np = size(outflow_data_new,dim=2)

  if (size(outflow_data_old,dim=2) .ne. size(outflow_data_new,dim=2)) then
     ! Do this so the routine has the right size
     deallocate(outflow_data_old)
     allocate(outflow_data_old(nc,np))
  end if

#ifndef AMREX_USE_CUDA
  if (size(outflow_data_old,dim=2) .ne. size(outflow_data_new,dim=2)) then
     print *,'size of old and new dont match in swap_outflow_data '
     call castro_error("Error:: Castro_nd.f90 :: swap_outflow_data")
  end if
#endif

  outflow_data_old(1:nc,1:np) = outflow_data_new(1:nc,1:np)

  if (outflow_data_new_time .ge. 0.e0_rt) &
       outflow_data_old_time = outflow_data_new_time
  outflow_data_new_time = -1.e0_rt

end subroutine swap_outflow_data

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

  integer :: iadv, ispec, iaux, ipassive

  integer :: ioproc


  ! easy indexing for the passively advected quantities.  This lets us
  ! loop over all groups (advected, species, aux) in a single loop.
  ! Note: these sizes are the maximum size we expect for passives.
  allocate(qpass_map(npassive))
  allocate(upass_map(npassive))

  ipassive = 1

  if (nadv > 0) then
     do iadv = 1, nadv
        upass_map(ipassive) = UFA + iadv - 1
        qpass_map(ipassive) = QFA + iadv - 1
        ipassive = ipassive + 1
     end do
  end if

  if (nspec > 0) then
     do ispec = 1, nspec
        upass_map(ipassive) = UFS + ispec - 1
        qpass_map(ipassive) = QFS + ispec - 1
        ipassive = ipassive + 1
     end do
  endif

  if (naux > 0) then
     do iaux = 1, naux
        upass_map(ipassive) = UFX + iaux - 1
        qpass_map(ipassive) = QFX + iaux - 1
        ipassive = ipassive + 1
     end do
  endif


  !---------------------------------------------------------------------
  ! other initializations
  !---------------------------------------------------------------------

  ! This is a routine which links to the C++ ParallelDescriptor class

  call bl_pd_is_ioproc(ioproc)

#ifdef ROTATION
  rot_vec = ZERO
  rot_vec(rot_axis) = ONE
#endif


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

  !$acc update device(small_dens, small_temp, small_pres, small_ener)

end subroutine ca_set_method_params


! :::
! ::: ----------------------------------------------------------------
! :::



subroutine ca_set_problem_params(dm,physbc_lo_in,physbc_hi_in,&
     Interior_in, Inflow_in, Outflow_in, &
     Symmetry_in, SlipWall_in, NoSlipWall_in, &
     coord_type_in, &
     problo_in, probhi_in, center_in) &
     bind(C, name="ca_set_problem_params")
     ! Passing data from C++ into f90
     !
     ! Binds to C function `ca_set_problem_params`

  use amrex_constants_module, only: ZERO
  use castro_error_module
  use prob_params_module
  use meth_params_module, only: UMX, UMY, UMZ
#ifdef ROTATION
  use meth_params_module, only: rot_axis
#endif
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: dm
  integer,  intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
  integer,  intent(in) :: Interior_in, Inflow_in, Outflow_in, Symmetry_in, SlipWall_in, NoSlipWall_in
  integer,  intent(in) :: coord_type_in
  real(rt), intent(in) :: problo_in(dm), probhi_in(dm), center_in(dm)

  allocate(dim)

  dim = dm

  allocate(physbc_lo(3))
  allocate(physbc_hi(3))

  physbc_lo(:) = 0
  physbc_hi(:) = 0

  physbc_lo(1:dm) = physbc_lo_in(1:dm)
  physbc_hi(1:dm) = physbc_hi_in(1:dm)

  allocate(Interior)
  allocate(Inflow)
  allocate(Outflow)
  allocate(Symmetry)
  allocate(SlipWall)
  allocate(NoSlipWall)

  allocate(coord_type)
  allocate(center(3))
  allocate(problo(3))
  allocate(probhi(3))

  Interior   = Interior_in
  Inflow     = Inflow_in
  Outflow    = Outflow_in
  Symmetry   = Symmetry_in
  SlipWall   = SlipWall_in
  NoSlipWall = NoSlipWall_in

  coord_type = coord_type_in

  problo = ZERO
  probhi = ZERO
  center = ZERO

  problo(1:dm) = problo_in(1:dm)
  probhi(1:dm) = probhi_in(1:dm)
  center(1:dm) = center_in(1:dm)

  allocate(dg(3))

  dg(:) = 1

  if (dim .lt. 2) then
     dg(2) = 0
  endif

  if (dim .lt. 3) then
     dg(3) = 0
  endif

#ifdef ROTATION
  if (coord_type == 1) then
     rot_axis = 2
  endif
#endif

  allocate(mom_flux_has_p(3))

  ! sanity check on our allocations
#ifndef AMREX_USE_CUDA
  if (UMZ > MAX_MOM_INDEX) then
     call castro_error("ERROR: not enough space in comp in mom_flux_has_p")
  endif
#endif

  ! keep track of which components of the momentum flux have pressure
  if (dim == 1 .or. (dim == 2 .and. coord_type == 1)) then
     mom_flux_has_p(1)%comp(UMX) = .false.
  else
     mom_flux_has_p(1)%comp(UMX) = .true.
  endif
  mom_flux_has_p(1)%comp(UMY) = .false.
  mom_flux_has_p(1)%comp(UMZ) = .false.

  mom_flux_has_p(2)%comp(UMX) = .false.
  mom_flux_has_p(2)%comp(UMY) = .true.
  mom_flux_has_p(2)%comp(UMZ) = .false.

  mom_flux_has_p(3)%comp(UMX) = .false.
  mom_flux_has_p(3)%comp(UMY) = .false.
  mom_flux_has_p(3)%comp(UMZ) = .true.

  !$acc update device(physbc_lo, physbc_hi)
  !$acc update device(Interior, Inflow, Outflow, Symmetry, Slipwall, NoSlipWall)
  !$acc update device(dim)
  !$acc update device(dg)
  !$acc update device(coord_type)
  !$acc update device(center, problo, probhi)
  !$acc update device(domlo_level, domhi_level, dx_level)
  !$acc update device(ref_ratio, n_error_buf, blocking_factor)
  !$acc update device(mom_flux_has_p)

end subroutine ca_set_problem_params

! :::
! ::: ----------------------------------------------------------------
! :::



subroutine ca_set_grid_info(max_level_in, dx_level_in, domlo_in, domhi_in, &
     ref_ratio_in, n_error_buf_in, blocking_factor_in) &
     bind(C, name="ca_set_grid_info")
     !
     ! Sometimes this routine can get called multiple
     ! times upon initialization; in this case, just to
     ! be safe, we'll deallocate and start again.
     !
     ! Binds to C function `ca_set_grid_info`

  use prob_params_module, only: max_level, dx_level, domlo_level, domhi_level, n_error_buf, ref_ratio, blocking_factor
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: max_level_in
  real(rt), intent(in) :: dx_level_in(3*(max_level_in+1))
  integer,  intent(in) :: domlo_in(3*(max_level_in+1)), domhi_in(3*(max_level_in+1))
  integer,  intent(in) :: ref_ratio_in(3*(max_level_in+1))
  integer,  intent(in) :: n_error_buf_in(0:max_level_in)
  integer,  intent(in) :: blocking_factor_in(0:max_level_in)

  integer :: lev, dir

  if (allocated(dx_level)) then
     deallocate(dx_level)
  endif
  if (allocated(domlo_level)) then
     deallocate(domlo_level)
  endif
  if (allocated(domhi_level)) then
     deallocate(domhi_level)
  endif

  if (allocated(ref_ratio)) then
     deallocate(ref_ratio)
  endif
  if (allocated(n_error_buf)) then
     deallocate(n_error_buf)
  endif
  if (allocated(blocking_factor)) then
     deallocate(blocking_factor)
  endif

  max_level = max_level_in

  allocate(dx_level(1:3, 0:max_level))
  allocate(domlo_level(1:3, 0:max_level))
  allocate(domhi_level(1:3, 0:max_level))
  allocate(ref_ratio(1:3, 0:max_level))
  allocate(n_error_buf(0:max_level))
  allocate(blocking_factor(0:max_level))

  do lev = 0, max_level
     do dir = 1, 3
        dx_level(dir,lev) = dx_level_in(3*lev + dir)
        domlo_level(dir,lev) = domlo_in(3*lev + dir)
        domhi_level(dir,lev) = domhi_in(3*lev + dir)
        ref_ratio(dir,lev) = ref_ratio_in(3*lev + dir)
     enddo
     n_error_buf(lev) = n_error_buf_in(lev)
     blocking_factor(lev) = blocking_factor_in(lev)
  enddo

end subroutine ca_set_grid_info


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
#ifndef AMREX_USE_CUDA
  if (namlen > maxlen) then
     call castro_error('probin file name too long')
  endif
#endif

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
#ifndef AMREX_USE_CUDA
     call castro_error('ERROR: problem in the tagging namelist')
#endif
  endif

  close (unit=un)

end subroutine ca_get_tagging_params

#ifdef SPONGE
! :::
! ::: ----------------------------------------------------------------
! :::



subroutine ca_read_sponge_params(name, namlen) bind(C, name="ca_read_sponge_params")
    ! Initialize the sponge parameters
    !
    ! Binds to C function `ca_get_sponge_params`

  use sponge_module
  use castro_error_module
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  real(rt) :: sponge_target_x_velocity
  real(rt) :: sponge_target_y_velocity
  real(rt) :: sponge_target_z_velocity

  namelist /sponge/ &
       sponge_lower_factor, sponge_upper_factor, &
       sponge_lower_radius, sponge_upper_radius, &
       sponge_lower_density, sponge_upper_density, &
       sponge_lower_pressure, sponge_upper_pressure, &
       sponge_target_x_velocity, &
       sponge_target_y_velocity, &
       sponge_target_z_velocity, &
       sponge_timescale

  ! Set namelist defaults

  sponge_lower_factor = 0.e0_rt
  sponge_upper_factor = 1.e0_rt

  sponge_lower_radius = -1.e0_rt
  sponge_upper_radius = -1.e0_rt

  sponge_lower_density = -1.e0_rt
  sponge_upper_density = -1.e0_rt

  sponge_lower_pressure = -1.e0_rt
  sponge_upper_pressure = -1.e0_rt

  sponge_target_x_velocity = 0.e0_rt
  sponge_target_y_velocity = 0.e0_rt
  sponge_target_z_velocity = 0.e0_rt

  sponge_timescale    = -1.e0_rt

  ! create the filename
#ifndef AMREX_USE_CUDA
  if (namlen > maxlen) then
     call castro_error('probin file name too long')
  endif
#endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! read in the namelist
  open (newunit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=sponge, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
#ifndef AMREX_USE_CUDA
     call castro_error('ERROR: problem in the sponge namelist')
#endif
  endif

  close (unit=un)

  sponge_target_velocity = [sponge_target_x_velocity, &
       sponge_target_y_velocity, &
       sponge_target_z_velocity]

  ! Sanity check

#ifndef AMREX_USE_CUDA
  if (sponge_lower_factor < 0.e0_rt .or. sponge_lower_factor > 1.e0_rt) then
     call castro_error('ERROR: sponge_lower_factor cannot be outside of [0, 1].')
  endif

  if (sponge_upper_factor < 0.e0_rt .or. sponge_upper_factor > 1.e0_rt) then
     call castro_error('ERROR: sponge_upper_factor cannot be outside of [0, 1].')
  endif
#endif

end subroutine ca_read_sponge_params


subroutine ca_get_sponge_params(sponge_lower_factor_in, sponge_upper_factor_in, &
                                sponge_lower_radius_in, sponge_upper_radius_in, &
                                sponge_lower_density_in, sponge_upper_density_in, &
                                sponge_lower_pressure_in, sponge_upper_pressure_in, &
                                sponge_target_velocity_in, &
                                sponge_timescale_in) bind(C, name="ca_get_sponge_params")

  use sponge_module

  implicit none

  real(rt), intent(inout) :: sponge_lower_factor_in
  real(rt), intent(inout) :: sponge_upper_factor_in
  real(rt), intent(inout) :: sponge_lower_radius_in
  real(rt), intent(inout) :: sponge_upper_radius_in
  real(rt), intent(inout) :: sponge_lower_density_in
  real(rt), intent(inout) :: sponge_upper_density_in
  real(rt), intent(inout) :: sponge_lower_pressure_in
  real(rt), intent(inout) :: sponge_upper_pressure_in
  real(rt), intent(inout) :: sponge_target_velocity_in(3)
  real(rt), intent(inout) :: sponge_timescale_in

  sponge_lower_factor_in = sponge_lower_factor
  sponge_upper_factor_in = sponge_upper_factor

  sponge_lower_radius_in = sponge_lower_radius
  sponge_upper_radius_in = sponge_upper_radius

  sponge_lower_density_in = sponge_lower_density
  sponge_upper_density_in = sponge_upper_density

  sponge_lower_pressure_in = sponge_lower_pressure
  sponge_upper_pressure_in = sponge_upper_pressure

  sponge_target_velocity_in(:) = sponge_target_velocity(:)

  sponge_timescale_in = sponge_timescale

end subroutine ca_get_sponge_params


subroutine ca_allocate_sponge_params() bind(C, name="ca_allocate_sponge_params")
    ! allocate sponge parameters
    !
    ! Binds to C function `ca_allocate_sponge_params`
    !

  use sponge_module
  allocate(sponge_lower_factor, sponge_upper_factor)
  allocate(sponge_lower_radius, sponge_upper_radius)
  allocate(sponge_lower_density, sponge_upper_density)
  allocate(sponge_lower_pressure, sponge_upper_pressure)
  allocate(sponge_target_velocity(3))
  allocate(sponge_timescale)

end subroutine ca_allocate_sponge_params



subroutine ca_deallocate_sponge_params() bind(C, name="ca_deallocate_sponge_params")
    ! deallocate sponge parameters
    !
    ! Binds to C function `ca_deallocate_sponge_params`
    !

  use sponge_module

  deallocate(sponge_lower_factor, sponge_upper_factor)
  deallocate(sponge_lower_radius, sponge_upper_radius)
  deallocate(sponge_lower_density, sponge_upper_density)
  deallocate(sponge_lower_pressure, sponge_upper_pressure)
  deallocate(sponge_target_velocity)
  deallocate(sponge_timescale)

end subroutine ca_deallocate_sponge_params
#endif

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



#ifdef GRAVITY

subroutine set_pointmass(pointmass_in) bind(C, name='set_pointmass')
    !
    ! pointmass_in real(rt)
    !

  use meth_params_module, only: point_mass
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: pointmass_in

  point_mass = pointmass_in

end subroutine set_pointmass
#endif
