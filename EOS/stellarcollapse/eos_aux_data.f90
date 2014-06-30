module eos_aux_data_module

  use bl_types

  implicit none

  ! for reading HDF5 table
  integer, save :: nrho,ntemp,nye
  real(kind=dp_t), save, allocatable :: eos_table(:,:,:,:)
  real(kind=dp_t), save, allocatable :: eos_logrho(:),eos_logtemp(:),eos_ye(:)
  real(kind=dp_t), save :: energy_shift = 0.0d0
  ! these are the indices in the eos_table
  ! we probably don't need all of these, but oh well
  integer, parameter :: eos_nvars = 19
  integer, parameter :: ilogpress = 1
  integer, parameter :: ilogenergy = 2
  integer, parameter :: ientropy = 3
  integer, parameter :: imunu = 4
  integer, parameter :: ics2 = 5
  integer, parameter :: idedt = 6
  integer, parameter :: idpdrhoe = 7
  integer, parameter :: idpderho = 8
  integer, parameter :: imuhat = 9
  integer, parameter :: imu_e = 10
  integer, parameter :: imu_p = 11
  integer, parameter :: imu_n = 12
  integer, parameter :: ixa = 13
  integer, parameter :: ixh = 14
  integer, parameter :: ixn = 15
  integer, parameter :: ixp = 16
  integer, parameter :: iabar = 17
  integer, parameter :: izbar = 18
  integer, parameter :: igamma = 19
  real(kind=dp_t), save :: eos_minrho, eos_maxrho
  real(kind=dp_t), save :: eos_mintemp, eos_maxtemp
  real(kind=dp_t), save :: eos_minye, eos_maxye
  ! this is a convenience for grabbing the ye auxilliary variable
  ! it should be initialized in eos_init
  integer, save :: iye_eos = -1

contains

  subroutine read_stellarcollapse_file(eos_input_file,use_energy_shift)

    use bl_error_module
    use hdf5
    use parallel

    implicit none

   character(len=256), intent(in) :: eos_input_file
   logical, intent(in) :: use_energy_shift

   integer(HID_T) :: file_id,dset_id,dspace_id
   integer(HSIZE_T) :: dims1(1),dims3(3)
   integer :: error,total_error


    if (trim(eos_input_file) == "") then
       call bl_error("EOS: eos_file not specified in probin!")
    endif

    ! initialize HDF5 library
    call h5open_f(error)
    if (error .ne. 0) &
         call bl_error("EOS: couldn't initialize HDF5 library")

    ! open the HDF5 eos file, and get the file id
    call h5fopen_f(trim(eos_input_file),H5F_ACC_RDONLY_F,file_id,error)
    if (error .ne. 0) &
         call bl_error("EOS: couldn't open eos_file for reading")
    
    ! get the number of density points
    ! open the dataset and get it's id
    dims1(1) = 1
    call h5dopen_f(file_id,"pointsrho",dset_id,error)
    ! read the data
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,nrho,dims1,error)
    ! close the dataset
    call h5dclose_f(dset_id,error)
    if (error .ne. 0) call bl_error("EOS: couldn't read nrho")

    ! get the number of temperature points
    call h5dopen_f(file_id,'pointstemp',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,ntemp,dims1,error)
    call h5dclose_f(dset_id,error)
    if (error .ne. 0) call bl_error("EOS: couldn't read ntemp")

    ! get the number of ye points
    call h5dopen_f(file_id,'pointsye',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,nye,dims1,error)
    call h5dclose_f(dset_id,error)
    if (error .ne. 0) call bl_error("EOS: couldn't read nye")

    allocate(eos_table(nrho,ntemp,nye,eos_nvars))

    ! thermo variables
    total_error = 0
    dims3(1) = nrho
    dims3(2) = ntemp
    dims3(3) = nye
    ! log of pressure (cgs)
    call h5dopen_f(file_id,'logpress',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ilogpress),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! log of energy (cgs)
    call h5dopen_f(file_id,'logenergy',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ilogenergy),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! entropy (k_B / baryon)
    call h5dopen_f(file_id,'entropy',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ientropy),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! square of (non-relativistic) sound speed (cgs)
    call h5dopen_f(file_id,'cs2',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ics2),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! derivative of energy wrt temperature (cgs)
    call h5dopen_f(file_id,'dedt',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,idedt),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! derivative of pressure wrt density at constant energy (cgs)
    call h5dopen_f(file_id,'dpdrhoe',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,idpdrhoe),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! derivative of pressure wrt energy at constant density (cgs)
    call h5dopen_f(file_id,'dpderho',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,idpderho),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error
    
    ! gamma_1
    call h5dopen_f(file_id,'gamma',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,igamma),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error


    ! chemical potentials (including rest mass)
    ! electron (MeV / baryon)
    call h5dopen_f(file_id,'mu_e',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,imu_e),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! proton (MeV / baryon)
    call h5dopen_f(file_id,'mu_p',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,imu_p),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! neutron (MeV / baryon)
    call h5dopen_f(file_id,'mu_n',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,imu_n),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! muhat = mu_n - mu_p  (MeV / baryon)
    call h5dopen_f(file_id,'muhat',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,imuhat),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! munu = mu_e - muhat  (MeV / baryon)
    call h5dopen_f(file_id,'munu',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,imunu),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error


    ! composition
    ! alpha particle mass fraction
    call h5dopen_f(file_id,'Xa',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ixa),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! 'heavy nucleus' mass fraction
    call h5dopen_f(file_id,'Xh',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ixh),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! neutron mass fraction
    call h5dopen_f(file_id,'Xn',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ixn),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! proton mass fraction
    call h5dopen_f(file_id,'Xp',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,ixp),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! 'heavy nucleus' nuclear properties
    ! average heavy nucleas A
    call h5dopen_f(file_id,'Abar',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,iabar),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! average heavy nucleus Z
    call h5dopen_f(file_id,'Zbar',dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE, &
                   eos_table(:,:,:,izbar),dims3,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    
    ! read in rho,t,ye grid
    ! log density (cgs)
    dims1(1) = nrho
    allocate(eos_logrho(nrho))
    call h5dopen_f(file_id,"logrho",dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,eos_logrho,dims1,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! log temperature (MeV)
    dims1(1) = ntemp
    allocate(eos_logtemp(ntemp))
    call h5dopen_f(file_id,"logtemp",dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,eos_logtemp,dims1,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error

    ! electron fraction
    dims1(1) = nye
    allocate(eos_ye(nye))
    call h5dopen_f(file_id,"ye",dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,eos_ye,dims1,error)
    call h5dclose_f(dset_id,error)
    total_error = total_error + error


    ! get the energy shift that makes energy consistent amongst different EOS's 
    ! that make up the table; this is in erg/g
    if (use_energy_shift) then
       call h5dopen_f(file_id,"energy_shift",dset_id,error)
       call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,energy_shift,dims1,error)
       call h5dclose_f(dset_id,error)
       total_error = total_error + error
    endif

    if (parallel_IOProcessor()) print *, 'energy_shift', energy_shift

    if (total_error .ne. 0) call bl_error("EOS: Error reading EOS table")
    
    ! close the file
    call h5fclose_f(file_id,error)
    ! cleanup library
    call h5close_f(error)

    eos_minrho = eos_logrho(1)
    eos_maxrho = eos_logrho(nrho)
    
    eos_mintemp = eos_logtemp(1)
    eos_maxtemp = eos_logtemp(ntemp)

    eos_minye = eos_ye(1)
    eos_maxye = eos_ye(nye)


  end subroutine read_stellarcollapse_file


! this converts from the units used in Castro to the units of the table
  subroutine convert_to_table_format(state)

    use eos_type_module
    use fundamental_constants_module, only: k_B, ev2erg, MeV2eV, n_A
    use bl_constants_module, only: ONE

    implicit none

    type(eos_t), intent(inout) :: state

    ! the stellarcollapse.org tables use some log10 variables, as well as 
    ! units of MeV for temperature and chemical potential, and k_B / baryon 
    ! for entropy
    state%rho = dlog10(state%rho)
    state%p = dlog10(state%p)
    state%e = dlog10(state%e - energy_shift)
    state%T = dlog10(state%T * k_B / ev2erg / MeV2eV)
    ! assuming baryon mass to be ~ 1 amu = 1/N_A
    state%s = state%s * k_B / n_A

  end subroutine convert_to_table_format

! this converts from the units used in the table to the units of Castro
! and builds some quantities, such as enthalpy
  subroutine convert_from_table_format(state)

    use fundamental_constants_module
    use eos_type_module
    use bl_constants_module, only: TEN
    use parallel
    
    implicit none

    type(eos_t), intent(inout) :: state

    state%rho = TEN**state%rho
    state%p   = TEN**state%p
    state%e   = TEN**state%e + energy_shift
    state%T = (TEN**state%T) * MeV2eV * ev2erg / k_B
    state%s = state%s * n_A / k_B

    ! construct enthalpy
    state%h = state%e + state%p/state%rho

  end subroutine convert_from_table_format

  subroutine table_lookup(state)
    ! this routine will populate the (known) state parameters by interpolating
    ! the EOS table, assuming density, temperature and ye are the independent
    ! variables
    
    use interpolate_module
    use eos_type_module

    implicit none

    type(eos_t), intent(inout) :: state

    double precision :: rho,temp,ye
    double precision :: derivs(3),cs2
    logical :: err
    character(len=128) :: errstring

    rho = state%rho
    temp = state%T
    ye = state%aux(iye_eos)

    write(errstring,'(3(e12.5,x))') rho,temp,ye

    ! by definition
    state%mu_e = 1./ye

    ! this should probably be vectorized...
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,ilogpress), &
                        state%p,derivs,err)
    ! this should only fail once because the same inputs are used throughout
    if (err) call bl_error('table_lookup: tri-interpolate failure:',trim(errstring))
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,ilogenergy), &
                        state%e,derivs,err)
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,ientropy), &
                        state%s,derivs,err)
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,ics2), &
                        cs2,derivs,err)
    state%cs = sqrt(cs2)
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,igamma), &
                        state%gam1,derivs,err)
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,idedt), &
                        state%dedT,derivs,err)
    call tri_interpolate(rho,temp,ye,nrho,ntemp,nye, &
                        eos_logrho,eos_logtemp,eos_ye, &
                        eos_table(:,:,:,idpdrhoe), &
                        state%dpdr_e,derivs,err)
    
  end subroutine table_lookup

end module eos_aux_data_module
