module model_parser_module

  ! read in an initial model and return arrays with the model data.
  ! take care to match the species available in the model file to
  ! those defined by the network

  ! the model file is assumed to be of the follow form:
  ! # npts = 896
  ! # num of variables = 6
  ! # density
  ! # temperature
  ! # pressure
  ! # carbon-12
  ! # oxygen-16
  ! # magnesium-24
  ! 195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
  ! 585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0
  !
  ! we read in the number of variables and their order and use this to
  ! map them into the model_state array.  We ignore anything other than
  ! density, temperature, pressure and composition.
  !
  ! composition is assumed to be in terms of mass fractions

  use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
  use network
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! integer keys for indexing the model_state array
  integer, parameter :: idens_model = 1
  integer, parameter :: itemp_model = 2
  integer, parameter :: ipres_model = 3

#ifdef INTERPOLATE_MODEL_VELOCITY
  integer, parameter :: ivelr_model = 4
  integer, parameter :: ispec_model = 5
#else
  integer, parameter :: ispec_model = 4
#endif

  integer, parameter :: nvars_model = ispec_model + nspec - 1

  ! number of points in the model file
  integer,   allocatable, save :: npts_model

  ! arrays for storing the model data
  real (rt), allocatable, save :: model_state(:,:)
  real (rt), allocatable, save :: model_r(:)

  ! model_initialized will be .true. once the model is read in and the
  ! model data arrays are initialized and filled
  logical, save :: model_initialized = .false.

  integer, parameter :: MAX_VARNAME_LENGTH=80

  public :: read_model_file, close_model_file, interpolate_sub, locate_sub

#ifdef AMREX_USE_CUDA
  attributes(managed) :: model_state, model_r, npts_model
#endif

contains

  subroutine read_model_file(model_file)

    use amrex_constants_module
    use castro_error_module

    character(len=*), intent(in   ) :: model_file

    ! local variables
    integer :: nvars_model_file
    integer :: ierr

    integer :: i, j, comp

    real(rt), allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found_model, found_dens, found_temp, found_pres, found_velr
    logical :: found_spec(nspec)
    integer :: ipos
    character (len=256) :: header_line

    allocate(npts_model)

    ! open the model file
    open(99,file=trim(model_file),status='old',iostat=ierr)

    if (ierr .ne. 0) then
       print *,'Couldnt open model_file: ',model_file
       call castro_error('Aborting now -- please supply model_file')
    end if

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model

    ! now read in the number of variables
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

    allocate (vars_stored(nvars_model_file))
    allocate (varnames_stored(nvars_model_file))

    ! now read in the names of the variables
    do i = 1, nvars_model_file
       read (99, '(a256)') header_line
       ipos = index(header_line, '#') + 1
       varnames_stored(i) = trim(adjustl(header_line(ipos:)))
    enddo

    ! allocate storage for the model data
    allocate (model_state(npts_model, nvars_model))
    allocate (model_r(npts_model))

887 format(78('-'))
889 format(a60)

    if ( amrex_pd_ioprocessor() ) then
       write (*,889) ' '
       write (*,887)
       write (*,*)   'reading initial model'
       write (*,*)   npts_model, 'points found in the initial model file'
       write (*,*)   nvars_model_file, ' variables found in the initial model file'
    endif


    ! start reading in the data
    do i = 1, npts_model
       read(99,*) model_r(i), (vars_stored(j), j = 1, nvars_model_file)

       model_state(i,:) = ZERO

       ! make sure that each of the variables that MAESTRO cares about
       ! are found
       found_dens = .false.
       found_temp = .false.
       found_pres = .false.
#ifdef INTERPOLATE_MODEL_VELOCITY
       found_velr = .false.
#endif
       found_spec(:) = .false.

       do j = 1,nvars_model_file

          ! keep track of whether the current variable from the model
          ! file is one that MAESTRO cares about
          found_model = .false.


          if (varnames_stored(j) == "density") then
             model_state(i,idens_model) = vars_stored(j)
             found_model = .true.
             found_dens  = .true.

          else if (varnames_stored(j) == "temperature") then
             model_state(i,itemp_model) = vars_stored(j)
             found_model = .true.
             found_temp  = .true.

          else if (varnames_stored(j) == "pressure") then
             model_state(i,ipres_model) = vars_stored(j)
             found_model = .true.
             found_pres  = .true.

#ifdef INTERPOLATE_MODEL_VELOCITY
          else if (varnames_stored(j) == "radial velocity") then
             model_state(i,ivelr_model) = vars_stored(j)
             found_model = .true.
             found_velr  = .true.
#endif

          else
             do comp = 1, nspec
                if (varnames_stored(j) == spec_names(comp)) then
                   model_state(i,ispec_model-1+comp) = vars_stored(j)
                   found_model = .true.
                   found_spec(comp) = .true.
                   exit
                endif
             enddo
          endif

          ! is the current variable from the model file one that we
          ! care about?
          if (.NOT. found_model .and. i == 1) then
             if ( amrex_pd_ioprocessor() ) then
                print *, 'WARNING: variable not found: ', &
                     trim(varnames_stored(j))
             end if
          endif

       enddo   ! end loop over nvars_model_file

       ! were all the variables we care about provided?
       if (i == 1) then
          if (.not. found_dens) then
             if ( amrex_pd_ioprocessor() ) then
                print *, 'WARNING: density not provided in inputs file'
             end if
          endif

          if (.not. found_temp) then
             if ( amrex_pd_ioprocessor() ) then
                print *, 'WARNING: temperature not provided in inputs file'
             end if
          endif

          if (.not. found_pres) then
             if ( amrex_pd_ioprocessor() ) then
                print *, 'WARNING: pressure not provided in inputs file'
             end if
          endif

#ifdef INTERPOLATE_MODEL_VELOCITY
          if (.not. found_velr) then
             if ( amrex_pd_ioprocessor() ) then
                print *, 'WARNING: radial velocity not provided in inputs file'
             end if
          endif
#endif

          do comp = 1, nspec
             if (.not. found_spec(comp)) then
                if ( amrex_pd_ioprocessor() ) then
                   print *, 'WARNING: ', trim(spec_names(comp)), &
                        ' not provided in inputs file'
                end if
             endif
          enddo
       endif

    end do   ! end loop over npts_model

    close(99)

    model_initialized = .true.

    deallocate(vars_stored,varnames_stored)

  end subroutine read_model_file


  function get_model_npts(model_file)

    integer :: get_model_npts

    ! look in the model file and return the number of points
    character(len=256), intent(in   ) :: model_file

    character (len=256) :: header_line
    integer :: ipos

    open(99,file=model_file)

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) get_model_npts

    close(99)

    return

  end function get_model_npts

  subroutine close_model_file

    if (model_initialized) then
       deallocate(model_r)
       deallocate(model_state)
       deallocate(npts_model)
       npts_model = -1
       model_initialized = .false.
    endif
  end subroutine close_model_file

  ! Initialize the model by passing in the number of points, r and state
  ! directly (rather than reading from a file)
  subroutine model_parser_init(npts, r, state)

      integer, intent(in) :: npts
      real(rt), intent(in) :: r(npts)
      real(rt), intent(in) :: state(npts, nvars_model)

      ! don't reinitialize
      if (.not. model_initialized) then

          allocate(npts_model)

          npts_model = npts

          ! allocate storage for the model data
          allocate (model_state(npts_model, nvars_model))
          allocate (model_r(npts_model))

          model_r(:) = r(:)
          model_state(:,:) = state(:,:)

          model_initialized = .true.

      endif

  end subroutine model_parser_init

  subroutine interpolate_sub(interp, r, var_index, iloc)
    ! find the value of model_state component var_index at point r
    ! using linear interpolation.  Eventually, we can do something
    ! fancier here.

    use amrex_fort_module, only : rt => amrex_real
    real(rt), intent(  out) :: interp
    real(rt), intent(in   ) :: r
    integer, intent(in) :: var_index
    integer, intent(in), optional   :: iloc

    ! Local variables
    integer                         :: id
    real(rt)                        :: slope,minvar,maxvar

    !$gpu

    !     find the location in the coordinate array where we want to interpolate
    if (present(iloc)) then
       id = iloc
    else
       call locate_sub(r, npts_model, model_r, id)
    end if

    if (id .eq. 1) then

       slope = (model_state(id+1, var_index) - model_state(id, var_index)) / &
            (model_r(id+1) - model_r(id))
       interp = slope*(r - model_r(id)) + model_state(id, var_index)

       ! safety check to make sure interp lies within the bounding points
       minvar = min(model_state(id+1, var_index), model_state(id, var_index))
       maxvar = max(model_state(id+1, var_index), model_state(id, var_index))
       interp = max(interp, minvar)
       interp = min(interp, maxvar)

    else if (id .eq. npts_model) then

       slope = (model_state(id, var_index) - model_state(id-1, var_index)) / &
            (model_r(id) - model_r(id-1))
       interp = slope*(r - model_r(id)) + model_state(id, var_index)

       ! safety check to make sure interp lies within the bounding points
       minvar = min(model_state(id, var_index), model_state(id-1, var_index))
       maxvar = max(model_state(id, var_index), model_state(id-1, var_index))
       interp = max(interp, minvar)
       interp = min(interp, maxvar)

    else

       if (r .ge. model_r(id)) then

          slope = (model_state(id+1, var_index) - model_state(id, var_index)) / &
               (model_r(id+1) - model_r(id))
          interp = slope*(r - model_r(id)) + model_state(id, var_index)

       else

          slope = (model_state(id, var_index) - model_state(id-1, var_index)) / &
               (model_r(id) - model_r(id-1))
          interp = slope*(r - model_r(id)) + model_state(id, var_index)

       end if

    endif

  end subroutine interpolate_sub




  subroutine locate_sub(x, n, xs, loc)

    use amrex_fort_module, only : rt => amrex_real
    integer,  intent(in   ) :: n
    real(rt), intent(in   ) :: x, xs(n)
    integer,  intent(  out) :: loc

    integer :: ilo, ihi, imid

    !$gpu

    if (x .le. xs(1)) then
       loc = 1
    else if (x .gt. xs(n-1)) then
       loc = n
    else

       ilo = 1
       ihi = n-1

       do while (ilo+1 .ne. ihi)
          imid = (ilo+ihi)/2
          if (x .le. xs(imid)) then
             ihi = imid
          else
             ilo = imid
          end if
       end do

       loc = ihi

    end if

  end subroutine locate_sub



end module model_parser_module
