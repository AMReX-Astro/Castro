module conservative_map_module

  ! read in an initial model and return arrays with the model data.
  ! take care to match the species available in the model file to
  ! those defined by the network.

  ! This version differs from the model_parser in that we expect
  ! the inputs file to have the conserved variables stored, with the
  ! same names that Castro would output in the plotfile.  We use the
  ! same indexing as with the Castro conserved variable state array.

  use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
  use network
  use amrex_fort_module, only : rt => amrex_real

  implicit none


  ! number of points in the model file
  integer,   allocatable, save :: npts_model

  ! arrays for storing the model data
  real (rt), allocatable, save :: model_state(:,:)
  real (rt), allocatable, save :: model_r(:)

  ! model_initialized will be .true. once the model is read in and the
  ! model data arrays are initialized and filled
  logical, save :: model_initialized = .false.

  integer, parameter :: MAX_VARNAME_LENGTH=80

  public :: read_conserved_model_file, close_conserved_model_file, interpolate_conservative

#ifdef AMREX_USE_CUDA
  attributes(managed) :: model_state, model_r, npts_model
#endif

contains

  subroutine read_conserved_model_file(model_file)

    use meth_params_module, only : NVAR, URHO, UMX, UTEMP, UEDEN, UEINT, UFS
    use amrex_constants_module
    use amrex_error_module

    character(len=*), intent(in   ) :: model_file

    ! local variables
    integer :: nvars_model_file
    integer :: ierr

    integer :: i, j, comp
    integer :: un

    real(rt), allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found_model, found_dens, found_xmom, found_rho_E, found_rho_eint, found_temp
    logical :: found_spec(nspec)
    integer :: ipos
    character (len=256) :: header_line

    allocate(npts_model)

    ! open the model file
    open(newunit=un, file=trim(model_file), status='old', iostat=ierr)

    if (ierr .ne. 0) then
       print *,'Couldnt open model_file: ',model_file
       call amrex_error('Aborting now -- please supply model_file')
    end if

    ! the first line has the number of points in the model
    read (un, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model

    ! now read in the number of variables
    read (un, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

    allocate (vars_stored(nvars_model_file))
    allocate (varnames_stored(nvars_model_file))

    ! now read in the names of the variables
    do i = 1, nvars_model_file
       read (un, '(a256)') header_line
       ipos = index(header_line, '#') + 1
       varnames_stored(i) = trim(adjustl(header_line(ipos:)))
    enddo

    ! allocate storage for the model data -- we use NVAR here, since
    ! it maps 1-to-1 to the conservative state
    allocate (model_state(npts_model, NVAR))
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
       read(un,*) model_r(i), (vars_stored(j), j = 1, nvars_model_file)

       model_state(i,:) = ZERO

       ! make sure that each of the variables that MAESTRO cares about
       ! are found
       found_dens = .false.
       found_xmom = .false.
       found_rho_E = .false.
       found_rho_eint = .false.
       found_temp = .false.
       found_spec(:) = .false.

       do j = 1, nvars_model_file

          ! keep track of whether the current variable from the model
          ! file is one that MAESTRO cares about
          found_model = .false.

          if (varnames_stored(j) == "density") then
             model_state(i,URHO) = vars_stored(j)
             found_model = .true.
             found_dens  = .true.

          else if (varnames_stored(j) == "Temp") then
             model_state(i,UTEMP) = vars_stored(j)
             found_model = .true.
             found_temp  = .true.

          else if (varnames_stored(j) == "xmom") then
             model_state(i,UMX) = vars_stored(j)
             found_model = .true.
             found_xmom  = .true.

          else if (varnames_stored(j) == "rho_E") then
             model_state(i,UEDEN) = vars_stored(j)
             found_model = .true.
             found_rho_E  = .true.

          else if (varnames_stored(j) == "rho_e") then
             model_state(i,UEINT) = vars_stored(j)
             found_model = .true.
             found_rho_eint  = .true.

          else
             do comp = 1, nspec
                if (varnames_stored(j) == "rho_" // trim(short_spec_names(comp))) then
                   model_state(i,UFS-1+comp) = vars_stored(j)
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
             call amrex_error("ERROR: density not provided in inputs file")
          end if

          if (.not. found_xmom) then
             call amrex_error("ERROR: x-momentum not provided in inputs file")
          end if

          if (.not. found_rho_E) then
             call amrex_error("ERROR: rho_E not provided in inputs file")
          end if

          if (.not. found_rho_eint) then
             call amrex_error("ERROR: rho_e not provided in inputs file")
          end if

          if (.not. found_temp) then
             call amrex_error("ERROR: Temp not provided in inputs file")
          end if

          do comp = 1, nspec
             if (.not. found_spec(comp)) then
                call amrex_error("ERROR: " // trim(spec_names(comp)), &
                     ' not provided in inputs file')
             end if
          end do
       end if

    end do   ! end loop over npts_model

    close(unit=un)

    model_initialized = .true.

    deallocate(vars_stored,varnames_stored)

  end subroutine read_conserved_model_file


  function get_model_npts(model_file)

    integer :: get_model_npts

    ! look in the model file and return the number of points
    character(len=256), intent(in   ) :: model_file

    character (len=256) :: header_line
    integer :: ipos, un

    open(newunit=un, file=model_file)

    ! the first line has the number of points in the model
    read (un, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) get_model_npts

    close(unit=un)

    return

  end function get_model_npts


  subroutine interpolate_conservative(interp, xl, xr, var_index)
    ! This interpolation routine is for conservative data at a
    ! different (assumed coarser) resolution than the model data.  We
    ! come in with the left and right edges of our zone, xl and xr,
    ! and we will average all of the model zones that fall inbetween
    ! as the interpolated value.
    !
    ! This assumes that the model data is uniformly spaced and that we
    ! are properly nested between the grid and the model.

    use amrex_constants_module, only : ZERO, HALF
    use amrex_error_module, only : amrex_error
    use amrex_fort_module, only : rt => amrex_real

    real(rt), intent(out) :: interp
    real(rt), intent(in) :: xl, xr
    integer, intent(in) :: var_index

    ! Local variables
    integer :: n
    integer :: ileft, iright, npts
    real(rt) :: x_model, x_model_l, x_model_r, dx, xscale
    real(rt), parameter :: tol = 1.e-12_rt
    !$gpu

    ! we assume that the model is uniformly spaced
    dx = model_r(2) - model_r(1)

    ! find the range of zones in the model that fit into our grid zone xl:xr
    ileft = -1
    iright = -1

    do n = 1, npts_model
       x_model = model_r(n)
       x_model_l = x_model - HALF*dx
       x_model_r = x_model + HALF*dx

       if (xl == ZERO) then
          xscale = 1
       else
          xscale = abs(xl)
       endif

       if (abs(x_model_l - xl) < tol*xscale) then
          if (ileft > 0) then
             call amrex_error("Error: ileft already set")
          else
             ileft = n
          end if
       end if

       if (abs(x_model_r - xr) < tol*abs(xr)) then
          if (iright > 0) then
             call amrex_error("Error: iright already set")
          else
             iright = n
          end if
       end if

       if (ileft >= 0 .and. iright >= 0) then
          exit
       end if
    end do

    if (ileft == -1 .or. iright == -1) then
       call amrex_error("Error: ileft or iright not set")
    end if

    if (iright < ileft) then
       call amrex_error("Error: iright < ileft")
    end if

    npts = iright - ileft + 1

    interp = sum(model_state(ileft:iright, var_index))/npts

  end subroutine interpolate_conservative


  subroutine close_conserved_model_file

    if (model_initialized) then
       deallocate(model_r)
       deallocate(model_state)
       deallocate(npts_model)
       npts_model = -1
       model_initialized = .false.
    endif
  end subroutine close_conserved_model_file

end module conservative_map_module
