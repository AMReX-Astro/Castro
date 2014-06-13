module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module
  use helmeos_module

  implicit none

  logical,         save, private :: do_coulomb
  logical,         save, private :: input_is_constant

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

    call helmeos_init(initialized)

  end subroutine eos_init



  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state_in, do_eos_diag, pt_index, state_len)

    ! A generic wrapper for the Helmholtz electron/positron degenerate EOS.  

    implicit none

    ! Input arguments

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state_in
    logical, optional, intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)
    integer, optional, intent(in   ) :: state_len

    integer :: j, N

    type (eos_t), allocatable :: state(:)

    ! Local variables and arrays
    
    double precision :: ymass(nspec), ysum, yzsum
    double precision :: e_want, p_want, s_want, h_want

    logical eosfail, eos_diag

    integer :: ns, ierr


    if(present(state_len)) then
      N = state_len
    else
      N = 1
    endif

    allocate(state(N))

    state(1) = state_in

    if (.not. initialized) call bl_error('EOS: not initialized')

    eos_diag = .false.

    if (present(do_eos_diag)) eos_diag = do_eos_diag

    do j = 1, N
       ! Check to make sure the composition was set properly.

       do ns = 1, nspec
         if (state(j) % xn(ns) .lt. init_test) call eos_error(ierr_init_xn, input, pt_index)
       enddo

       ! Get abar, zbar, etc.

       call composition(state(j), .false.)
    enddo

    eosfail = .false.

    ierr = 0

!---------------------------------------------------------------------------
! dens, temp, and xmass are inputs
!---------------------------------------------------------------------------


       ! Call the EOS.

       call helmeos(do_coulomb, eosfail, state, N, input)

!       if (eosfail) call eos_error(ierr_general, input, pt_index)



!---------------------------------------------------------------------------
! dens, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_rh) then

!       if (state % rho .lt. init_test .or. state % h .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given enthalpy.

!       h_want = state % h

!       if (h_want < ZERO) call eos_error(ierr_neg_h, input, pt_index)

!       if (eos_diag) print *, 'WANT h ', h_want

!       call newton_iter(state, eos_diag, ierr, ienth, itemp, h_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) state % h = h_want



!---------------------------------------------------------------------------
! temp, pres, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_tp) then

!       if (state % T .lt. init_test .or. state % p .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given pressure
!       p_want = state % p

!       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

!       if (eos_diag) print *, 'WANT p ', p_want       

!       call newton_iter(state, eos_diag, ierr, ipres, idens, p_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) state % p = p_want



!---------------------------------------------------------------------------
! dens, pres, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_rp) then

!       if (state % rho .lt. init_test .or. state % p .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given pressure
!       p_want = state % p

!       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

!       if (eos_diag) print *, 'WANT p ', p_want       

!       call newton_iter(state, eos_diag, ierr, ipres, itemp, p_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) state % p = p_want



!---------------------------------------------------------------------------
! dens, energy, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_re) then

!       if (state % rho .lt. init_test .or. state % e .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given energy
!       e_want = state % e

!       if (e_want < ZERO) call eos_error(ierr_neg_e, input, pt_index)

!       if (eos_diag) print *, 'WANT e ', e_want

!       call newton_iter(state, eos_diag, ierr, iener, itemp, e_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) state % e = e_want
              


!---------------------------------------------------------------------------
! pres, entropy, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_ps) then

!       if (state % p .lt. init_test .or. state % s .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given entropy and pressure
!       s_want = state % s
!       p_want = state % p

!       if (s_want < ZERO) call eos_error(ierr_neg_s, input, pt_index)
!       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)

!       if (eos_diag) then
!          print *, 'WANT s ', s_want
!          print *, 'WANT p ', p_want
!       endif

!       call newton_iter2(state, eos_diag, ierr, ipres, p_want, ientr, s_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) then
!          state % s = s_want
!          state % p = p_want
!       endif



!---------------------------------------------------------------------------
! pres, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_ph) then

!       if (state % p .lt. init_test .or. state % h .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given enthalpy and pressure
!       p_want = state % p
!       h_want = state % h

!       if (p_want < ZERO) call eos_error(ierr_neg_p, input, pt_index)
!       if (h_want < ZERO) call eos_error(ierr_neg_h, input, pt_index)

!       if (eos_diag) then
!          print *, 'WANT p ', p_want
!          print *, 'WANT h ', h_want
!       endif

!       call newton_iter2(state, eos_diag, ierr, ipres, p_want, ienth, h_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) then
!          state % p = p_want
!          state % h = h_want
!       endif



!---------------------------------------------------------------------------
! temp, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

!    elseif (input .eq. eos_input_th) then

!       if (state % t .lt. init_test .or. state % h .lt. init_test) call eos_error(ierr_init, input, pt_index)

!       if (eos_diag) print *, 'T/D INIT ', state % T, state % rho

       ! We want to converge to the given enthalpy
!       h_want = state % h

!       if (eos_diag) print *, 'WANT h ', h_want

!       if (h_want < ZERO) call eos_error(ierr_neg_h, input, pt_index)

!       call newton_iter(state, eos_diag, ierr, ienth, idens, h_want)

!       if (ierr > 0) call eos_error(ierr, input, pt_index)

!       if (input_is_constant) state % h = h_want



!---------------------------------------------------------------------------
! The EOS input doesn't match any of the available options.
!---------------------------------------------------------------------------

!    else

!       call eos_error(ierr_input, input, pt_index)

!    endif



    ! Take care of final housekeeping.

    ! Count the positron contribution in the electron quantities.
    state(:) % xne  = state(:) % xne  + state(:) % xnp
    state(:) % pele = state(:) % pele + state(:) % ppos

    ! Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
    ! This replaces the relativistic version that comes out of helmeos.

    state(:) % cs = sqrt(state(:) % gam1 * state(:) % p / state(:) % rho)

    ! Get dpdX, dedX, dhdX.

    do j = 1, N
       call composition_derivatives(state(j), .false.)
    enddo

    state_in = state(1)

  end subroutine eos



  subroutine eos_error(err, input, pt_index)

    implicit none

    integer,           intent(in) :: err
    integer,           intent(in) :: input
    integer, optional, intent(in) :: pt_index(:)

    integer :: dim_ptindex

    character (len=64) :: err_string, zone_string, eos_input_str

    write(eos_input_str, '(A13, I1)') ' EOS input = ', input

    if (err .eq. ierr_general) then

      err_string = 'EOS: error in the EOS.'

    elseif (err .eq. ierr_input) then

      err_string = 'EOS: invalid input.'

    elseif (err .eq. ierr_iter_conv) then

      err_string = 'EOS: Newton-Raphson iterations failed to converge.'

    elseif (err .eq. ierr_neg_e) then

      err_string = 'EOS: energy < 0 in the EOS.'

    elseif (err .eq. ierr_neg_p) then

      err_string = 'EOS: pressure < 0 in the EOS.'

    elseif (err .eq. ierr_neg_h) then

      err_string = 'EOS: enthalpy < 0 in the EOS.'

    elseif (err .eq. ierr_neg_s) then

      err_string = 'EOS: entropy < 0 in the EOS.'

    elseif (err .eq. ierr_init) then

      err_string = 'EOS: the input variables were not initialized.'

    elseif (err .eq. ierr_init_xn) then

      err_string = 'EOS: the species abundances were not initialized.'

    elseif (err .eq. ierr_iter_var) then

      err_string = 'EOS: the variable you are iterating over was not recognized.'

    else

      err_string = 'EOS: invalid input to error handler.'

    endif

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
