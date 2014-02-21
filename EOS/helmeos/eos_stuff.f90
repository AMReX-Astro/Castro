module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  logical,         save, private :: do_coulomb
  integer,         save, private :: max_newton = 100
  logical,         save, private :: input_is_constant

  double precision, save, private :: ttol = 1.0d-8
  double precision, save, private :: dtol = 1.0d-8

  public eos_init, eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and CASTRO
  ! For this general EOS, this calls helmeos_init() which reads in the 
  ! table with the electron component's properties.
  subroutine eos_init(small_temp, small_dens)

    use parallel
    use extern_probin_module, only: use_eos_coulomb, eos_input_is_constant

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
 
    do_coulomb = use_eos_coulomb
    input_is_constant = eos_input_is_constant 

    smallt = 1.d4

    if (present(small_temp)) then
      if (small_temp > ZERO) then
       smallt = small_temp
      end if
    endif

    smalld = 1.d-5
 
    if (present(small_dens)) then
       if (small_dens > ZERO) then
         smalld = small_dens
       endif
    endif

    if (parallel_IOProcessor()) print *, 'Initializing helmeos... Coulomb corrections = ', do_coulomb

    ! Call the helmeos initialization routine and read in the table 
    ! containing the electron contribution.

    call helmeos_init()

    initialized = .true.
 
  end subroutine eos_init



  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag, pt_index)

    ! A generic wrapper for the Helmholtz electron/positron degenerate EOS.  

    implicit none

    ! Input arguments

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state 
    logical, optional, intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables and arrays
    
    double precision :: ymass(nspec), ysum, yzsum
    double precision :: e_want, p_want, s_want, h_want

    logical eosfail, eos_diag

    integer :: n, ierr



    if (.not. initialized) call bl_error('EOS: not initialized')

    eos_diag = .false.

    if (present(do_eos_diag)) eos_diag = do_eos_diag

    ! Check to make sure the composition was set properly.

    do n = 1, nspec
      if (state % xn(n) .lt. init_test) call eos_error(ierr_init_xn, input, pt_index)
    enddo

    ! Get abar, zbar, etc.

    call composition(state, .false.)

    eosfail = .false.

    ierr = 0

    select case (input)

!---------------------------------------------------------------------------
! dens, temp, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_rt)

       if (state % rho .lt. init_test .or. state % T .lt. init_test) call eos_error(ierr_init, input, pt_index)

       ! Call the EOS.

       call helmeos(do_coulomb, eosfail, state)

       if (eosfail) call eos_error(ierr_general, input, pt_index)



!---------------------------------------------------------------------------
! dens, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_rh)

       if (state % rho .lt. init_test .or. state % h .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given enthalpy.

       h_want = state % h

       if (h_want < ZERO) call eos_error(ierr_neg_h, input, pt_index)

       if (eos_diag) print *, 'WANT h ', h_want

       call newton_iter(state, eos_diag, ierr, ienth, itemp, h_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) state % h = h_want



!---------------------------------------------------------------------------
! temp, pres, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_tp)

       if (state % T .lt. init_test .or. state % p .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given pressure
       p_want = state % p

       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

       if (eos_diag) print *, 'WANT p ', p_want       

       call newton_iter(state, eos_diag, ierr, ipres, idens, p_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) state % p = p_want



!---------------------------------------------------------------------------
! dens, pres, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_rp)

       if (state % rho .lt. init_test .or. state % p .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given pressure
       p_want = state % p

       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

       if (eos_diag) print *, 'WANT p ', p_want       

       call newton_iter(state, eos_diag, ierr, ipres, itemp, p_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) state % p = p_want



!---------------------------------------------------------------------------
! dens, energy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_re)

       if (state % rho .lt. init_test .or. state % e .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given energy
       e_want = state % e

       if (e_want < ZERO) call eos_error(ierr_neg_e, input, pt_index)

       if (eos_diag) print *, 'WANT e ', e_want

       call newton_iter(state, eos_diag, ierr, iener, itemp, e_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) state % e = e_want
              


!---------------------------------------------------------------------------
! pres, entropy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_ps)

       if (state % p .lt. init_test .or. state % s .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given entropy and pressure
       s_want = state % s
       p_want = state % p

       if (s_want < ZERO) call eos_error(ierr_neg_s, input, pt_index)
       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

       if (eos_diag) then
          print *, 'WANT s ', s_want
          print *, 'WANT p ', p_want
       endif

       call newton_iter2(state, eos_diag, ierr, ipres, p_want, ientr, s_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) then
          state % s = s_want
          state % p = p_want
       endif



!---------------------------------------------------------------------------
! pres, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_ph)

       if (state % p .lt. init_test .or. state % h .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given enthalpy and pressure
       s_want = state % s
       h_want = state % h

       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)
       if (h_want < ZERO) call eos_error(ierr_neg_h, input, pt_index)

       if (eos_diag) then
          print *, 'WANT p ', p_want
          print *, 'WANT h ', h_want
       endif

       call newton_iter2(state, eos_diag, ierr, ipres, p_want, ienth, h_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) then
          state % p = p_want
          state % h = h_want
       endif



!---------------------------------------------------------------------------
! temp, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

    case (eos_input_th)

       if (state % t .lt. init_test .or. state % h .lt. init_test) call eos_error(ierr_init, input, pt_index)

       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given enthalpy
       h_want = state % h

       if (eos_diag) print *, 'WANT h ', h_want

       if (h_want < ZERO) call eos_error(ierr_neg_h, input, pt_index)

       call newton_iter(state, eos_diag, ierr, ienth, idens, h_want)

       if (ierr > 0) call eos_error(ierr, input, pt_index)

       if (input_is_constant) state % h = h_want



!---------------------------------------------------------------------------
! The EOS input doesn't match any of the available options.
!---------------------------------------------------------------------------

    case default 

       call eos_error(ierr_input, input, pt_index)

    end select



    ! Take care of final housekeeping.

    ! Count the positron contribution in the electron quantities.
    state % xne  = state % xne  + state % xnp
    state % pele = state % pele + state % ppos

    ! Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
    ! This replaces the relativistic version that comes out of helmeos.

    state % cs = sqrt(state % gam1 * state % p / state % rho)

    ! Get dpdX, dedX, dhdX.

    call composition_derivatives(state, .false.)

  end subroutine eos



  subroutine newton_iter(state, eos_diag, ierr, var, dvar, f_want)

     implicit none

     type (eos_t),       intent(inout) :: state
     integer,            intent(in   ) :: var, dvar
     double precision,   intent(in   ) :: f_want
     logical,            intent(in   ) :: eos_diag
     integer,            intent(inout) :: ierr

     integer          :: iter
     double precision :: smallx, error, xnew, xtol
     double precision :: f, x, dfdx

     logical :: converged, eosfail

     if (.not. (dvar .eq. itemp .or. dvar .eq. idens) ) then
       ierr = ierr_iter_var
       return
     endif

     converged = .false.
     xnew = ZERO

     do iter = 1, max_newton

        ! For each iteration, start by filling the state with the EOS

        call helmeos(do_coulomb, eosfail, state)

        if (eosfail) then
          ierr = ierr_general
          return
        endif

        ! First, figure out what variable we're working with

        if (dvar .eq. itemp) then

          x = state % T

          smallx = smallt
          xtol = ttol

          select case (var)

            case (ipres)
              f    = state % p
              dfdx = state % dpdT
            case (iener)
              f    = state % e
              dfdx = state % dedT
            case (ientr)
              f    = state % s
              dfdx = state % dsdT
            case (ienth)
              f    = state % h
              dfdx = state % dhdT
            case default
              ierr = ierr_iter_var
              return

          end select

        else ! dvar == density

          x = state % rho

          smallx = smalld
          xtol = dtol

          select case (var)

            case (ipres)
              f    = state % p
              dfdx = state % dpdr
            case (iener)
              f    = state % e
              dfdx = state % dedr
            case (ientr)
              f    = state % s
              dfdx = state % dsdr
            case (ienth)
              f    = state % h
              dfdx = state % dhdr
            case default
              ierr = ierr_iter_var
              return
 
          end select

        endif

        ! Now do the calculation for the next guess for T/rho
 
        if (eos_diag) then
          print *, 'VAR  = ', var , iter, ' f    = ', f
          print *, 'DVAR = ', dvar, iter, ' dfdx = ', dfdx
        endif

        xnew = x - (f - f_want) / dfdx

        if (eos_diag) then
          print *, 'XNEW FIRST ', x, ' - ', f - f_want, ' / ', dfdx
        endif

        ! Don't let the temperature/density change by more than a factor of two
        xnew = max(HALF * x, min(xnew, TWO * x))

        ! Don't let us freeze/evacuate
        xnew = max(smallx, xnew)

        if (eos_diag) then
          print *, 'XNEW AFTER ', iter, xnew
        endif

        ! Compute the error

        error = abs( (xnew - x) / x )

        if (error .lt. xtol) then
          converged = .true.
          return
        endif

        ! Store the new temperature/density if we're still iterating

        if (dvar .eq. itemp) then
          state % T    = xnew
        else
          state % rho  = xnew
        endif
               
     enddo

     ! Call error if too many iterations are needed

     if (.not. converged) ierr = ierr_iter_conv

  end subroutine newton_iter



  subroutine newton_iter2(state, eos_diag, ierr, var1, f_want, var2, g_want)

     implicit none

     type (eos_t),       intent(inout) :: state
     integer,            intent(in   ) :: var1, var2
     double precision,   intent(in   ) :: f_want, g_want
     logical,            intent(in   ) :: eos_diag
     integer,            intent(inout) :: ierr

     integer          :: iter
     double precision :: error1, error2, fi, gi, rnew, tnew, delr
     double precision :: f, dfdt, dfdr
     double precision :: g, dgdt, dgdr
     double precision :: temp, dens

     logical :: converged, eosfail

     converged = .false.     

     ! First pass

     rnew = ZERO
     tnew = ZERO

     do iter = 1, max_newton

        ! Start each iteration by filling the state with the EOS

        call helmeos(do_coulomb, eosfail, state)

        if (eosfail) then
          ierr = ierr_general
          return
        endif

        ! First, figure out which variables we're using
 
        temp = state % T
        dens = state % rho

        select case (var1)

           case (ipres)
             f    = state % p
             dfdt = state % dpdT
             dfdr = state % dpdr
           case (iener)
             f    = state % e
             dfdt = state % dedT
             dfdr = state % dedr
           case (ientr)
             f    = state % s
             dfdt = state % dsdT
             dfdr = state % dsdr
           case (ienth)
             f    = state % h
             dfdT = state % dhdT
             dfdr = state % dhdr
           case default
             ierr = ierr_iter_var
             return

         end select

         select case (var2)

           case (ipres)
             g    = state % p
             dgdt = state % dpdT
             dgdr = state % dpdr
           case (iener)
             g    = state % e
             dgdt = state % dedT
             dgdr = state % dedr
           case (ientr)
             g    = state % s
             dgdt = state % dsdT
             dgdr = state % dsdr
           case (ienth)
             g    = state % h
             dgdt = state % dhdT
             dgdr = state % dhdr
           case default
             ierr = ierr_iter_var
             return

         end select

         if (eos_diag) then
           print *, 'VAR1 ', var1, iter, f
           print *, 'VAR2 ', var2, iter, g
         end if

        ! Two functions, f and g, to iterate over
        fi = f_want - f
        gi = g_want - g

        !
        ! 0 = f + dfdr * delr + dfdt * delt
        ! 0 = g + dgdr * delr + dgdt * delt
        !

        delr = (fi*dgdt - gi*dfdt) / (dgdr*dfdt - dgdt*dfdr)

        rnew = dens + delr

        tnew = temp - (fi + dfdr*delr) / dfdt

        if (eos_diag) then
           print *, 'RNEW FIRST ', dens, ' + ', &
                fi*dgdt - gi*dfdt, ' / ', dgdr*dfdt - dgdt*dfdr
           print *, 'TNEW FIRST ', temp, ' - ', &
                fi + dfdr*delr, ' / ', dfdt
        endif

        ! Don't let the temperature or density change by more
        ! than a factor of two
        tnew = max(HALF * temp, min(tnew, TWO * temp))
        rnew = max(HALF * dens, min(rnew, TWO * dens))

        ! Don't let us freeze or evacuate
        tnew = max(smallt, tnew)
        rnew = max(smalld, rnew)

        if (eos_diag) then
           print *, 'RNEW AFTER ', iter, rnew
           print *, 'TNEW AFTER ', iter, tnew
        endif

        ! Compute the errors
        error1 = abs( (rnew - dens) / dens )
        error2 = abs( (tnew - temp) / temp )

        if (error1 .LT. dtol .and. error2 .LT. ttol) then
          converged = .true.
          exit
        endif
     
        ! Store the new temperature and density if we're still iterating
        state % rho = rnew
        state % T   = tnew
                
     enddo

     ! Call error if too many iterations are needed

     if (.not. converged) ierr = ierr_iter_conv

  end subroutine newton_iter2



  subroutine eos_error(err, input, pt_index)

    implicit none

    integer,           intent(in) :: err
    integer,           intent(in) :: input
    integer, optional, intent(in) :: pt_index(:)

    integer :: dim_ptindex

    character (len=64) :: err_string, zone_string, eos_input_str

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

      case default

        err_string = 'EOS: invalid input to error handler.'

    end select

    err_string = err_string // eos_input_str

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


end module eos_module
