module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module, only: ZERO, HALF, TWO
  use eos_type_module
  use eos_data_module
  use eos_aux_data_module

  implicit none

  integer,         save, private :: max_newton = 100

  double precision, save, private :: ttol = 1.0d-8
  double precision, save, private :: dtol = 1.0d-8

  character(len=15) :: errfmt = '(3(e12.5,x))'

  public eos_init, eos, get_munu

contains

  ! EOS initialization routine 
  ! this reads in the HDF5 file containing the tabulated data
  subroutine eos_init(small_temp, small_dens)

    use parallel
    use extern_probin_module, only: eos_file, use_energy_shift
    use network, only: network_species_index

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
 
    if (parallel_IOProcessor()) print *, 'Reading HDF5 file', eos_file
    call read_stellarcollapse_file(eos_file,use_energy_shift)

    ! make sure we have iye_eos
    iye_eos = network_species_index('ye')
    if (iye_eos .lt. 0) call bl_error("EOS_INIT: couldn't find iye_eos in network")
    
    ! note that the small* values are NOT in the units used in the table
    ! but we write them here as other places in the code may refer to them
    smallt = eos_mintemp
    if (present(small_temp)) then
      if (small_temp > ZERO) then
       smallt = small_temp
      end if
    endif

    smalld = eos_minrho
    if (present(small_dens)) then
       if (small_dens > ZERO) then
         smalld = small_dens
       endif
    endif

    initialized = .true.
 
  end subroutine eos_init


  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag, pt_index)

    ! A generic wrapper for the stellarcollapse EOS
    ! 
    ! The stellarcollapse tables are indexed by log(density), 
    ! log(temperature), and electron fraction.  As such, the usual
    ! 'composition' variable passed to the EOS is the electron fraction.
    !
    ! Make sure you use a network that uses ye as a species!

    implicit none

    ! Input arguments
    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state 
    logical, optional, intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables and arrays
    double precision :: e_want, p_want, s_want, h_want
    double precision, parameter :: tol = 1.0d-8

    logical :: eos_diag

    integer :: n, ierr

    logical :: input_is_constant

    input_is_constant = .true.

    if (.not. initialized) call bl_error('EOS: not initialized')

    eos_diag = .false.

    if (present(do_eos_diag)) eos_diag = do_eos_diag

    ! Check to make sure the electron fraction was set properly
    if (state%aux(iye_eos) .lt. init_test) call eos_error(ierr_general, input, pt_index)
    ! Make sure it fits in the table.
    if ((state%aux(iye_eos) .lt. eos_minye) .or. (state%aux(iye_eos) .gt. eos_maxye)) &
         call eos_error(ierr_out_of_bounds, input, pt_index)

    ! Convert to the units used by the table.
    call convert_to_table_format(state)

    ierr = 0

    select case (input)

!---------------------------------------------------------------------------
! dens, temp, and ye are inputs; 
! this is direct table interpolation, so fall through after bounds checking
!---------------------------------------------------------------------------

    case (eos_input_rt)

       ! Make sure these were set properly
       if (state % rho .lt. init_test .or. state % T .lt. init_test) call eos_error(ierr_init, input, pt_index)
       ! Make sure they fit in the table.
       if ((state%rho .lt. eos_minrho) .or. (state%rho .gt. eos_maxrho) .or. &
            (state%T .lt. eos_mintemp) .or. (state%T .gt. eos_maxtemp)) &
            call eos_error(ierr_out_of_bounds, input, pt_index)


!---------------------------------------------------------------------------
! dens, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_rh)
       ! NOT CURRENTLY IMPLEMENTED
       call eos_error(ierr_not_implemented, input, pt_index)

!---------------------------------------------------------------------------
! temp, pres, and ye are inputs; iterate to find density
!---------------------------------------------------------------------------

    case (eos_input_tp)

       ! check inputs
       if (state % T .lt. init_test .or. state % p .lt. init_test) call eos_error(ierr_init, input, pt_index)
       ! make sure the input temp is within table
       if ((state%T .lt. eos_mintemp) .or. (state%T .gt. eos_maxtemp)) &
            call eos_error(ierr_out_of_bounds, input, pt_index)
       ! make sure the initial density guess is within table
       state%rho = max(eos_minrho, min(eos_maxrho, state%rho))

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given pressure
       p_want = state % p
       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

       if (eos_diag) print *, 'WANT p ', p_want       
       
       call newton_iter(state, eos_diag, ierr, ipres, idens, p_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)




!---------------------------------------------------------------------------
! dens, pres, and ye are inputs; iterate to find the temperature
!---------------------------------------------------------------------------

    case (eos_input_rp)

       ! check inputs
       if (state % rho .lt. init_test .or. state % p .lt. init_test) call eos_error(ierr_init, input, pt_index)
       ! make sure the input density is within the table
       if ((state%rho .lt. eos_minrho) .or. (state%rho .gt. eos_maxrho)) &
            call eos_error(ierr_out_of_bounds, input, pt_index)
       ! make sure the initial temperature guess is within the table
       state%T = max(eos_mintemp, min(eos_maxtemp, state%T))

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given pressure
       p_want = state % p
       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

       if (eos_diag) print *, 'WANT p ', p_want       

       call newton_iter(state, eos_diag, ierr, ipres, itemp, p_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)



!---------------------------------------------------------------------------
! dens, energy, and ye are inputs; iterate to find temperature
!---------------------------------------------------------------------------

    case (eos_input_re)

       ! check inputs
       if (state % rho .lt. init_test .or. state % e .lt. init_test) call eos_error(ierr_init, input, pt_index)
       ! make sure the input density is within the table
       if ((state%rho .lt. eos_minrho) .or. (state%rho .gt. eos_maxrho)) &
            call eos_error(ierr_out_of_bounds, input, pt_index)
       ! make sure the initial guess for temperature is within the table
       state%T = max(eos_mintemp, min(state%T, eos_maxtemp))

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given energy
       e_want = state % e
       if (e_want < ZERO) call eos_error(ierr_neg_e, input, pt_index)

       if (eos_diag) print *, 'WANT e ', e_want

       ! iterate to get the temperature
       call newton_iter(state, eos_diag, ierr, iener, itemp, e_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)



!---------------------------------------------------------------------------
! pres, entropy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_ps)
       ! NOT CURRENTLY IMPLEMENTED
       call eos_error(ierr_not_implemented, input, pt_index)


!---------------------------------------------------------------------------
! pres, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_ph)
       ! NOT CURRENTLY IMPLEMENTED
       call eos_error(ierr_not_implemented, input, pt_index)


!---------------------------------------------------------------------------
! temp, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_th)
       ! NOT CURRENTLY IMPLEMENTED
       call eos_error(ierr_not_implemented, input, pt_index)


!---------------------------------------------------------------------------
! The EOS input doesn't match any of the available options.
!---------------------------------------------------------------------------

    case default 

       call eos_error(ierr_input, input, pt_index)

    end select

    ! do a final lookup - by now we should have a consistent density and temperature
    call table_lookup(state)


    ! convert back to hydro units from table units
    ! also builds some quantities, like enthalpy
    call convert_from_table_format(state)


  end subroutine eos



  subroutine newton_iter(state, eos_diag, ierr, var, dvar, f_want)

    use interpolate_module

     implicit none

     type (eos_t),       intent(inout) :: state
     integer,            intent(in   ) :: var, dvar
     double precision,   intent(in   ) :: f_want
     logical,            intent(in   ) :: eos_diag
     integer,            intent(inout) :: ierr

     integer          :: iter, ivar
     double precision :: smallx, error, xnew, xtol
     double precision :: f, x, dfdx, df(3)

     logical :: converged, err
     character(len=128) :: errstring

     if (.not. (dvar .eq. itemp .or. dvar .eq. idens) ) then
       ierr = ierr_iter_var
       return
     endif

     converged = .false.
     err = .false.

     ! find out which table variable we are interpolating for
     select case(var)
     case (ipres)
        ivar = ilogpress
     case (iener)
        ivar = ilogenergy
     case (ientr)
        ivar = ientropy
     case default
        call bl_error("newton_iter: don't know how to handle var",var)
     end select

     do iter = 1, max_newton

        ! If we're converged, exit the loop
        if (converged) return

        ! interpolate the table for var; df is dfdrho,dfdT,dfdye
        call tri_interpolate(state%rho,state%T,state%aux(iye_eos), &
                             nrho,ntemp,nye, &
                             eos_logrho,eos_logtemp,eos_ye, &
                             eos_table(:,:,:,ivar), &
                             f, df, err)
        if (err) then
           write(errstring,trim(errfmt)) state%rho,state%T,state%aux(iye_eos)
           call bl_error('newton iter: failure to interpolate',trim(errstring))
        endif

        ! Figure out what variable we're working with
        if (dvar .eq. itemp) then

          x = state % T
          ! note that we are in table units
          smallx = eos_mintemp
          xtol = ttol

          dfdx = df(2)

        else ! dvar == density

          x = state % rho
          ! note that we are in table units
          smallx = eos_minrho
          xtol = dtol

          dfdx = df(1)

        endif

        ! Now do the calculation for the next guess for T/rho
        if (eos_diag) then
          print *, 'VAR  = ', var , iter, ' f    = ', f
          print *, 'DVAR = ', dvar, iter, ' dfdx = ', dfdx
        endif

        xnew = x - (f - f_want) / dfdx

        if (eos_diag) then
          print *, 'XNEW FIRST ', x, ' - ', f - f_want, ' / ', dfdx,' = ', xnew
        endif

        ! Don't let the temperature/density change by more than a factor of two
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NOTE that temperature is in log(MeV) and can be negative
        ! account for this by sign check; hopefully we aren't crossing zero...
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (x .lt. ZERO) then
           xnew = max(TWO * x, min(xnew, HALF * x))
        else
           xnew = max(HALF * x, min(xnew, TWO * x))
        endif

        ! Don't let us freeze/evacuate
        xnew = max(smallx, xnew)

        if (eos_diag) then
          print *, 'XNEW AFTER ', iter, xnew
        endif

        ! Store the new temperature/density
        if (dvar .eq. itemp) then
          state % T    = xnew
        else
          state % rho  = xnew
        endif

        ! Compute the error from the last iteration
        error = abs( (xnew - x) / x )

        if (error .lt. xtol) converged = .true.
               
     enddo

     ! Call error if too many iterations were needed
     if (.not. converged) ierr = ierr_iter_conv

  end subroutine newton_iter



  subroutine eos_error(err, input, pt_index)

    implicit none

    integer,           intent(in) :: err
    integer,           intent(in) :: input
    integer, optional, intent(in) :: pt_index(:)

    integer :: dim_ptindex

    character (len=128) :: err_string, zone_string, eos_input_str

    write(eos_input_str, '(A13, I1)') ' EOS input = ', input

    select case (err)

      case (ierr_general)

        err_string = 'EOS: error in the EOS.'

      case (ierr_input)

        err_string = 'EOS: invalid input.'

      case (ierr_iter_conv)

        err_string = 'EOS: Newton-Raphson iterations failed to converge.'

      case (ierr_neg_e)

        err_string = 'EOS: energy < 0 in the EOS.'

      case (ierr_neg_p)

        err_string = 'EOS: pressure < 0 in the EOS.'

      case (ierr_neg_h)

        err_string = 'EOS: enthalpy < 0 in the EOS.'

      case (ierr_neg_s)

        err_string = 'EOS: entropy < 0 in the EOS.'

      case (ierr_init)
 
        err_string = 'EOS: the input variables were not initialized.'

      case (ierr_init_xn)

        err_string = 'EOS: the species abundances were not initialized.'

      case (ierr_iter_var)

        err_string = 'EOS: the variable you are iterating over was not recognized.'

     case (ierr_out_of_bounds)
        
        err_string = 'EOS: an input variable value  is out of bounds of the table.'

     case (ierr_not_implemented)
        
        err_string = 'EOS: input method is currently not implemented for this EOS.'

      case default

        err_string = 'EOS: invalid input to error handler.'

    end select

    err_string = trim(err_string) // trim(eos_input_str)

    ! this format statement is for writing into zone_string -- make sure that
    ! the len of z_err can accomodate this format specifier
1001 format(1x,"zone index info: i = ", i5)
1002 format(1x,"zone index info: i = ", i5, '  j = ', i5)
1003 format(1x,"zone index info: i = ", i5, '  j = ', i5, '  k = ', i5)

    if (present(pt_index)) then
 
       dim_ptindex = size(pt_index,dim=1)

       if (dim_ptindex .eq. 1) then 
          write (zone_string,1001) pt_index(1)
       else if (dim_ptindex .eq. 2) then 
          write (zone_string,1002) pt_index(1), pt_index(2)
       else if (dim_ptindex .eq. 3) then 
          write (zone_string,1003) pt_index(1), pt_index(2), pt_index(3)
       end if

    else

      zone_string = ''

    endif

    call bl_error(err_string, zone_string)

  end subroutine eos_error


  function get_munu(rho,T,ye) result(munu)
    use interpolate_module 

    real(kind=dp_t), intent(in   ) :: rho, T, ye
    real(kind=dp_t)                :: munu

    type(eos_t) :: state
    real(kind=dp_t) :: derivs(3)
    logical :: err
    character(len=128) :: errstring

    ! convert our values to table format before interpolating
    state%rho = rho
    state%T = T
    state%aux(iye_eos) = ye
    call convert_to_table_format(state)

    ! look it up
    call tri_interpolate(state%rho,state%T,state%aux(iye_eos),&
                         nrho,ntemp,nye, &
                         eos_logrho,eos_logtemp,eos_ye, &
                         eos_table(:,:,:,imunu),&
                         munu,derivs,err)

    ! check return
    if (err) then
       write(errstring,trim(errfmt)) state%rho,state%T,state%aux(iye_eos)
       call bl_error('get_munu: tri-interpolate failure:',trim(errstring))
    endif
    
  end function get_munu


end module eos_module
