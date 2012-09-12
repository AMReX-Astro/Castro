module eos_module

  use bl_space, only: MAX_SPACEDIM
  use bl_types
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module

  implicit none

  private 

  real(kind=dp_t), public :: xn_eos(nspec)
  real(kind=dp_t), public :: temp_eos
  real(kind=dp_t), public :: den_eos
  real(kind=dp_t), public :: abar_eos
  real(kind=dp_t), public :: zbar_eos
  real(kind=dp_t), public :: e_eos
  real(kind=dp_t), public :: p_eos
  real(kind=dp_t), public :: h_eos
  real(kind=dp_t), public :: cv_eos
  real(kind=dp_t), public :: cp_eos
  real(kind=dp_t), public :: xne_eos
  real(kind=dp_t), public :: eta_eos
  real(kind=dp_t), public :: pele_eos
  real(kind=dp_t), public :: dpdt_eos
  real(kind=dp_t), public :: dpdr_eos
  real(kind=dp_t), public :: dedr_eos
  real(kind=dp_t), public :: dedt_eos
  real(kind=dp_t), public :: gam1_eos
  real(kind=dp_t), public ::   cs_eos
  real(kind=dp_t), public ::    s_eos
  real(kind=dp_t), public :: dsdt_eos
  real(kind=dp_t), public :: dsdr_eos
  real(kind=dp_t), public :: dpdX_eos(nspec)
  real(kind=dp_t), public :: dhdX_eos(nspec)
  real(kind=dp_t), public :: conduct_eos

  integer, public         :: pt_index_eos(MAX_SPACEDIM)

  common /eos_common/ xn_eos,temp_eos,den_eos,abar_eos,zbar_eos,e_eos,p_eos,h_eos
  common /eos_common/ cv_eos,cp_eos,xne_eos,eta_eos,pele_eos,dpdt_eos,dpdr_eos,dedr_eos
  common /eos_common/ dedt_eos,gam1_eos,cs_eos,s_eos,dsdt_eos,dsdr_eos,dpdX_eos,dhdX_eos
  common /eos_common/ conduct_eos,pt_index_eos
  SAVE /eos_common/
!$omp threadprivate(/eos_common/)

  integer, parameter, public :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter, public :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter, public :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter, public :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter, public :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter, public :: eos_input_ps = 6  ! p, s are inputs
 

  logical, save, private :: do_coulomb
  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld

  logical, save, private :: initialized = .false.

  private nspec, aion, zion

  public eos_init, eos_get_small_temp, eos_get_small_dens, eos_given_ReX, &
         eos_e_given_RPX, eos_S_given_ReX, eos_given_RTX, eos_dpdr_given_RTX, &
         eos_given_TPX, eos_given_PSX, eos

  interface eos
     module procedure eos_old
     module procedure eos_new
  end interface eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  ! For this general EOS, this calls helmeos_init() which reads in the 
  ! table with the electron component's properties.
  subroutine eos_init(small_temp, small_dens, gamma_in)

    use parallel
    use extern_probin_module, only: use_eos_coulomb

    implicit none
 
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens

    ! gamma_in is a dummy variable -- it is needed in a generic interface
    ! for an EOS, but only used in a gamma-law EOS, not this general EOS
    real(kind=dp_t), intent(in), optional :: gamma_in
 
    do_coulomb = use_eos_coulomb
 
    if (present(small_temp)) then
      if (small_temp > 0.d0) then
       smallt = small_temp
      else
       smallt = 1.d4
      end if
    else
       smallt = 1.d4
    endif
 
    if (present(small_dens)) then
       if (small_dens > 0.d0) then
         smalld = small_dens
       else
         smalld = 1.d-5
       end if
    else
       smalld = 1.d-5
    endif

    if (parallel_IOProcessor()) print *, 'Initializing helmeos... Coulomb corrections = ', do_coulomb
    ! call the helmeos initialization routine and read in the table 
    ! containing the electron contribution.
    call helmeos_init()
    initialized = .true.
 
  end subroutine eos_init


  !---------------------------------------------------------------------------
  ! Castro interfaces 
  !---------------------------------------------------------------------------
  subroutine eos_get_small_temp(small_temp_out)
 
    real(kind=dp_t), intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
 
    real(kind=dp_t), intent(out) :: small_dens_out
 
    small_dens_out = smalld
 
  end subroutine eos_get_small_dens

  subroutine eos_given_ReX(G, P, C, T, dpdr_e, dpde, R, e, X, pt_index)

    ! note: here, dpdr_e is partial p / partial rho at constant e 
    !       and   dpde   is partial p / partial e   at constant rho


     ! In/out variables
     real(kind=dp_t), intent(  out) :: G, P, C, dpdr_e, dpde
     real(kind=dp_t), intent(inout) :: T
     real(kind=dp_t), intent(in   ) :: R, e, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos = T
      den_eos = R
        e_eos = e
      xn_eos(1:nspec) = X(1:nspec)

     call eos(eos_input_re, den_eos, temp_eos, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag, pt_index)

    G  = gam1_eos
    P  =    p_eos
    C  =   cs_eos
    T  = temp_eos
    dpdr_e = dpdr_eos - dpdt_eos*dedr_eos/dedt_eos
    dpde = dpdt_eos / dedt_eos

  end subroutine eos_given_ReX

  subroutine eos_e_given_RPX(e, T, R, P, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e
     real(kind=dp_t), intent(in   ) :: R, p, X(:)
     real(kind=dp_t), intent(inout) :: T
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos = T
      den_eos = R
        p_eos = P
      xn_eos(1:nspec) = X(1:nspec)

     call eos(eos_input_rp, den_eos, temp_eos, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag, pt_index)

    e  =    e_eos
    T  = temp_eos

  end subroutine eos_e_given_RPX

  subroutine eos_S_given_ReX(S, R, e, T, X, pt_index)

     implicit none

     ! In/out variables
     real(kind=dp_t), intent(  out) :: S
     real(kind=dp_t), intent(in   ) :: R, e, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos = T
      den_eos = R
        e_eos = e
      xn_eos(1:nspec) = X(1:nspec)

     call eos(eos_input_re, den_eos, temp_eos, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag, pt_index)

    S  = s_eos

  end subroutine eos_S_given_ReX

  subroutine eos_given_RTX(e, P, R, T, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e, P
     real(kind=dp_t), intent(in   ) :: R, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

      den_eos = R
     temp_eos = T
      xn_eos(1:nspec) = X(1:nspec)

     call eos(eos_input_rt, den_eos, temp_eos, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag, pt_index)

    P  =    p_eos
    e  =    e_eos

  end subroutine eos_given_RTX

  subroutine eos_dpdr_given_RTX(e, P, R, T, X, dpdr, pt_index)

    ! note: here, dpdr is partial p / partial rho at constant T
    ! this is different than the dpdr_e that Castro uses for source
    ! terms in the primitive variable formulation.

    ! In/out variables
    real(kind=dp_t), intent(  out) :: e, P, dpdr
    real(kind=dp_t), intent(in   ) :: R, T, X(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    logical :: do_diag

    do_diag = .false.

    den_eos = R
    temp_eos = T
    xn_eos(1:nspec) = X(1:nspec)

    call eos(eos_input_rt, den_eos, temp_eos, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag, pt_index)
    
    P  =    p_eos
    e  =    e_eos
    dpdr =  dpdr_eos

  end subroutine eos_dpdr_given_RTX

  subroutine eos_given_TPX(e, P, R, T, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e
     real(kind=dp_t), intent(inout) :: R
     real(kind=dp_t), intent(in   ) :: P, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     ! An initial guess of density needs to be given
     den_eos = R 
     p_eos = P
     temp_eos = T
     xn_eos(1:nspec) = X(1:nspec)

     call eos(eos_input_tp, den_eos, temp_eos, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag, pt_index)

    R  =    den_eos
    e  =    e_eos

  end subroutine eos_given_TPX

  subroutine eos_given_PSX(S, P, X, e, R, T, pt_index)
    
    ! In/out variables
    real(kind=dp_t), intent(  out) :: e
    real(kind=dp_t), intent(inout) :: R, T
    real(kind=dp_T), intent(in   ) :: S, P, X(:)
    integer, optional, intent(in  ) :: pt_index(:)

    ! Local variables
    logical :: do_diag

    do_diag = .false.

    s_eos    = S
    p_eos    = P
    den_eos  = R
    temp_eos = T
    xn_eos(1:nspec) = X(1:nspec)

    call eos(eos_input_ps, den_eos, temp_eos, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag, pt_index)

    R = den_eos
    T = temp_eos
    e = e_eos

  end subroutine eos_given_PSX

  !---------------------------------------------------------------------------
  ! new interface
  !---------------------------------------------------------------------------
  subroutine eos_new(input, eos_state, do_eos_diag, pt_index)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: eos_state
    logical,           intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

    call eos_old(input, eos_state%rho, eos_state%T, &
                 eos_state%xn, &
                 eos_state%p, eos_state%h, eos_state%e, &
                 eos_state%cv, eos_state%cp, eos_state%xne, &
                 eos_state%eta, eos_state%pele, &
                 eos_state%dpdT, eos_state%dpdr, &
                 eos_state%dedT, eos_state%dedr, &
                 eos_state%dpdX, eos_state%dhdX, &
                 eos_state%gam1, eos_state%cs, eos_state%s, &
                 eos_state%dsdT, eos_state%dsdr, &
                 do_eos_diag, pt_index)

  end subroutine eos_new


  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos_old(input, dens, temp, &
                     xmass, &
                     pres, enthalpy, eint, &
                     c_v, c_p, ne, eta, pele, &
                     dPdT, dPdR, dEdT, dEdR, &
                     dPdX, dhdX, &
                     gam1, cs, entropy, &
                     dsdT, dsdR, &
                     do_eos_diag, &
                     pt_index)

    use bl_error_module

! a generic wrapper for the Helmholtz based electron/positron degenerate
! EOS.  
!
! dens     -- mass density (g/cc)
! temp     -- temperature (K)
! xmass    -- the mass fractions of the individual isotopes
! pres     -- the pressure (dyn/cm**2)
! enthalpy -- the enthalpy (erg/g)
! eint     -- the internal energy (erg/g)
! c_v      -- specific heat at constant volume
! c_p      -- specific heat at constant pressure
! ne       -- number density of electrons + positrons
! eta      -- degeneracy parameter
! pele     -- electron pressure + positron pressure
! dPdT     -- d pressure/ d temperature
! dPdR     -- d pressure/ d density
! dEdT     -- d energy/ d temperature
! dEdR     -- d energy/ d density
! dPdX     -- d pressure / d xmass
! dhdX     -- d enthalpy / d xmass  -- AT CONSTANT PRESSURE!!!
! gam1     -- first adiabatic index (d log P/ d log rho) |_s
! cs       -- sound speed -- note that this is the non-relativistic one
!             (we compute it in this wrapper as sqrt(gam1 p /rho) instead
!             of taking the relativistic version from helmeos.
! entropy  -- entropy (erg/g/K)
! dsdT     -- d entropy / d temperature
! dsdR     -- d entropy / d density
!
! input = 1 means dens, temp    , and xmass are inputs, return enthalpy, eint
!       = 2 means dens, enthalpy, and xmass are inputs, return temp    , eint
!                (note, temp should be filled with an initial guess)
!       = 3 means temp, pres    , and xmass are inputs, return dens    , etc
!       = 4 means dens, pres    , and xmass are inputs, return temp    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc
!       = 6 means pres, entropy , and xmass are inputs, return dens    , etc
!
!
! derivatives wrt X_k:
!
!   The EOS does not return the thermodynamic derivatives with respect
!   to the mass fractions, but rather, only due to the average atomic
!   mass (abar) and average proton number (zbar):
!
!     abar = ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
!
!     zbar = ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )
!
!   using the chain rule:
!
!   dp/dX_k = dp/d(abar) d(abar)/dX_k  +  dp/d(zbar) d(zbar)/dX_k
!
!   and using the above definitions of abar and zbar and sum_k {X_k} = 1
!
!   d(abar)/dX_k = abar * (1 - abar/A_k)
!   d(zbar)/dX_k = (Z_k - zbar) / ( A_k * sum_i {X_i/A_i} )
!
    implicit none

!    include 'vector_eos.dek'


!     ::::: Arguments
    logical             :: do_eos_diag
    integer, intent(in) :: input

    ! some of these quantites can be inputs or outputs
    real(kind=dp_t), intent(inout) :: dens, temp
    real(kind=dp_t), intent(in)    :: xmass(nspec)
    real(kind=dp_t), intent(inout) :: pres, enthalpy, &
                                      eint, entropy

    ! these quantities are always outputs
    real(kind=dp_t), intent(out) :: c_v, c_p
    real(kind=dp_t), intent(out) :: ne, eta, pele
    real(kind=dp_t), intent(out) :: dPdT, dPdR, &
                                    dedT, dedR
    real(kind=dp_t), intent(out) :: gam1
    real(kind=dp_t), intent(out) :: cs
    real(kind=dp_t), intent(out) :: dPdX(nspec), &
                                    dhdX(nspec)
    real(kind=dp_t), intent(out) :: dsdT, dsdR

    integer, optional, intent(in   ) :: pt_index(:)

    
!     ::::: Local variables and arrays

    integer :: i, n, iter, niter, max_newton
    parameter (max_newton = 100)
    
    real(kind=dp_t) :: error, error2
    real(kind=dp_t) :: ymass(nspec)
    real(kind=dp_t) :: abar, zbar
    real(kind=dp_t) :: energy_want
    real(kind=dp_t) :: enthalpy_want
    real(kind=dp_t) :: pres_want
    real(kind=dp_t) :: entropy_want
    real(kind=dp_t) :: dhdt
    real(kind=dp_t) :: tnew
    real(kind=dp_t) :: dnew
    real(kind=dp_t) :: enth1
    real(kind=dp_t) :: ener1
    real(kind=dp_t) :: dedX(nspec)

    real(kind=dp_t) :: dpdd, pres1, entr1
    real(kind=dp_t) :: f, g, dfdd, dfdt, dgdd, dgdt, deld

    real(kind=dp_t) :: ttol
    parameter (ttol = 1.0d-8)
    real(kind=dp_t) :: dtol
    parameter (dtol = 1.0d-8)
    real(kind=dp_t) :: stol
    parameter (stol = 1.0d-8)

    logical eosfail
    integer dim_ptindex

    ! err_string is used to convert the pt_index information into a string
    character (len=64) :: err_string  

!     ::::: Input/Output arrays for call to helmeos
    real(kind=dp_t) :: temp_row, den_row, abar_row, &
                     zbar_row, etot_row, ptot_row, &
                     cv_row, cp_row, &
                     xne_row, xnp_row, etaele_row, &
                     pele_row, ppos_row, dpd_row, &
                     dpt_row, dpa_row, dpz_row, &
                     ded_row, det_row, dea_row, &
                     dez_row, &
                     stot_row, dsd_row, dst_row
    real(kind=dp_t) :: gam1_row, cs_row

    if (present(pt_index)) dim_ptindex = size(pt_index,dim=1)
      
    if (.not. initialized) call bl_error('EOS: not initialized')

    ! this format statement is for writing into err_string -- make sure that
    ! the len of err_string can accomodate this format specifier
1001 format(1x,"zone index info: i = ", i5)
1002 format(1x,"zone index info: i = ", i5, '  j = ', i5)
1003 format(1x,"zone index info: i = ", i5, '  j = ', i5, '  k = ', i5)

    tnew  = 0.0d0
    dnew   = 0.0d0

    do i=1,nspec
       ymass(i) = xmass(i)/aion(i)
       dnew    = dnew + ymass(i)
       tnew    = tnew + zion(i) * ymass(i)
    enddo

    abar = 1.0d0/dnew
    zbar = tnew * abar

    if (input .EQ. eos_input_rt) then

!---------------------------------------------------------------------------
! input = 1: dens, temp, and xmass are inputs
!---------------------------------------------------------------------------

! we are taking density, temperature, and composition as given
       temp_row = temp
       den_row  = dens
       abar_row = abar
       zbar_row = zbar
         
! call the eos
       call helmeos(do_coulomb,eosfail, &
                   temp_row, den_row, abar_row, zbar_row, &
                   etot_row, ptot_row, &
                   cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                   pele_row, ppos_row, &
                   dpd_row, dpt_row, dpa_row, dpz_row, &
                   ded_row, det_row, dea_row, dez_row, & 
                   gam1_row, cs_row, stot_row, &
                   dsd_row, dst_row)

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

       ! fill the outputs
       pres = ptot_row
       eint = etot_row
         
       enthalpy = eint + pres/dens

       c_v = cv_row
       c_p = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       ne   = xne_row + xnp_row
       eta  = etaele_row
       pele = pele_row + ppos_row

       dPdR = dpd_row
       dPdT = dpt_row
       dEdR = ded_row
       dEdT = det_row
       gam1 = gam1_row
!       cs =   cs_row
       cs =   sqrt(gam1*pres/dens)
       entropy = stot_row
       dsdT = dst_row
       dsdR = dsd_row

       do n = 1, nspec
          dpdX(n) = dpa_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          dhdX(n) = dEdX(n) + &
               (pres/dens**2 - dEdR)*dpdX(n)/dPdr

       enddo

    else if (input .EQ. eos_input_rh) then

!---------------------------------------------------------------------------
! input = 2: dens, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = temp
       den_row  = dens
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',temp,dens

       ! we want to converge to the given enthalpy
       enthalpy_want = enthalpy

       if (do_eos_diag) print*,'WANT H ',enthalpy

       call helmeos(do_coulomb,eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

       do iter = 1, max_newton

          niter = iter

          ! recompute the enthalpy and it's temperature derivative
          enth1 = etot_row + ptot_row/dens
          if (do_eos_diag) print*,'ENTH1 ',iter,enth1

          dhdt = det_row + dpt_row/dens
          if (do_eos_diag) print*,'DHDT ',iter,dhdt

          tnew = temp_row - &
               (enth1 - enthalpy_want)/dhdt

          if (do_eos_diag) then
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  enth1 - enthalpy_want, ' / ', dhdt
          endif

          ! don't let the temperature change by more than a factor of two
          tnew = max(.5d0*temp_row, &
                        min(tnew, 2.d0*temp_row))

          ! don't let us freeze
          tnew = max(smallt, tnew)

          if (do_eos_diag) print*,'TNEW AFTER ',iter,tnew

          ! compute the error
          error = 0.0d0
          error = max(error,abs(tnew - temp_row)/temp_row)

          ! store the new temperature
          temp_row = tnew

          if (error .LT. ttol) goto 70
        
          call helmeos(do_coulomb,eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       enddo

       ! Land here if too many iterations are needed

       continue

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations')
       endif

70     continue

       ! store the end result
       temp = tnew
       pres = ptot_row
       eint = etot_row
       
       c_v = cv_row
       c_p = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       ne   = xne_row + xnp_row
       eta  = etaele_row
       pele = pele_row + ppos_row
       
       dPdR = dpd_row
       dPdT = dpt_row
       dEdR = ded_row
       dEdT = det_row   ! c_v
       gam1 = gam1_row
!      cs =   cs_row
       cs =   sqrt(gam1*pres/dens)
       entropy = stot_row
       dsdT = dst_row
       dsdR = dsd_row

       do n = 1, nspec
          dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          dhdX(n) = dEdX(n) + &
               (pres/dens**2 - dEdR)*dpdX(n)/dPdr

       enddo

    else if (input .EQ. eos_input_tp ) then

!---------------------------------------------------------------------------
! input = 3: temp, pres, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = temp
       den_row  = dens
       abar_row = abar
       zbar_row = zbar

       ! we want to converge to the given pressure
       pres_want = pres

       if (pres_want < ZERO) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('ERROR: pressure < 0 in the EOS', err_string)
          else
             call bl_error('ERROR: pressure < 0 in the EOS')
          endif
       endif
         
       call helmeos(do_coulomb,eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)
       
       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

       do iter = 1, max_newton

          niter = iter

          ! recompute the density and it's temperature derivative
          pres1 = ptot_row
          dpdd  = dpd_row
          
          dnew = den_row - &
               (pres1 - pres_want)/dpdd

          ! don't let the density change by more than an order of magnitude
          dnew = max(.5d0*den_row, &
                        min(dnew, 2.d0*den_row))

          ! compute the error
          error = 0.0d0
          error = max(error,abs(dnew - den_row)/den_row)

          ! store the new density
          den_row = dnew

          ! check if we are evacuating, if so, set the density to smalld, and adjust
          ! the error so we iterate on this one
          if (den_row .LT. smalld) then
             den_row = smalld
             error = 1.1d0*dtol
          endif

          if (error .LT. dtol) goto 170
        
          call helmeos(do_coulomb,eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)
        
          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif

       end do

       ! Land here if too many iterations are needed

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:3: too many iterations', err_string)         
       else
          call bl_error('EOS: Newton-Raphson failed:3: too many iterations')         
       endif

170    continue

       ! store the end result
       dens = dnew
       temp = temp_row
       eint = etot_row
       enthalpy = eint + ptot_row/dens

       c_v = cv_row
       c_p = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       ne   = xne_row + xnp_row
       eta  = etaele_row
       pele = pele_row + ppos_row
       
       dPdR = dpd_row
       dPdT = dpt_row
       dEdR = ded_row
       dEdT = det_row   ! c_v
       gam1 = gam1_row
!      cs =   cs_row
       cs =   sqrt(gam1*pres/dens)
       entropy = stot_row
       dsdT = dst_row
       dsdR = dsd_row

       do n = 1, nspec
          dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          dhdX(n) = dEdX(n) + &
               (pres/dens**2 - dEdR)*dpdX(n)/dPdr

       enddo

    else if (input .EQ. eos_input_rp ) then

!---------------------------------------------------------------------------
! input = 4: dens, pres, and xmass are inputs
!---------------------------------------------------------------------------

       ! Load the initial guess
       temp_row = temp
       den_row  = dens
       abar_row = abar
       zbar_row = zbar
       if (do_eos_diag) print*,'T/D INIT ',temp,dens

       ! We want to converge to the given pressure
       pres_want = pres
       if (do_eos_diag) print*,'P WANT ',pres

       if (pres_want < ZERO) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('ERROR: pressure < 0 in the EOS', err_string)
          else
             call bl_error('ERROR: pressure < 0 in the EOS')
          endif
       endif

       ! First pass
       call helmeos(do_coulomb,eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif
       
       do iter = 1, max_newton

          niter = iter

          ! recompute the temperature but dont allow it to change too much
          tnew = temp_row - &
               (ptot_row - pres_want)/dpt_row

          if (do_eos_diag) print*,'PRES ',ptot_row,pres_want
          if (do_eos_diag) print*,'PRES DIFF ',ptot_row-pres_want
          
          if (do_eos_diag) print*,'DPDT FAC ', 1.0/dpt_row

          if (do_eos_diag) print*,'TNEW BEFORE MAX ',iter,tnew

          ! don't let the temperature change by more than a factor of 2
          tnew = max(.5d0*temp_row, &
                        min(tnew, 2.d0*temp_row))

          ! don't let us freeze
          tnew = max(smallt, tnew)

          if (do_eos_diag) print*,'TNEW AFTER MAX ',iter,tnew
          if (do_eos_diag) print*,' '

          ! compute the error and store the new temperature
          error = 0.0d0
          error = max(error,abs(tnew - temp_row)/temp_row)
          if (do_eos_diag) print *,'ERROR  ',iter,error
          if (do_eos_diag) print*,' '
          temp_row = tnew

          if (error .LT. ttol) goto 870
        
          call helmeos(do_coulomb,eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       end do

       ! Land here if too many iterations are needed

       print *, 'helmeos input==4 failed to converge, iter = ',niter
       print *, 'error, temp_row, den_row = ', error, temp_row, den_row

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:4: too many iterations', err_string)          
       else
          call bl_error('EOS: Newton-Raphson failed:4: too many iterations')          
       endif

870    continue

       ! store the end result
       ! jbb
       ! temp = tnew
       temp = temp_row
       dens = den_row
       eint = etot_row
       c_v = cv_row
       c_p = cp_row

       enthalpy = eint + ptot_row/dens

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       ne   = xne_row + xnp_row
       eta  = etaele_row
       pele = pele_row + ppos_row
       
       dPdR = dpd_row
       dPdT = dpt_row
       dEdR = ded_row
       dEdT = det_row   ! c_v
       gam1 = gam1_row
!      cs =   cs_row
       cs =   sqrt(gam1*pres/dens)
       entropy = stot_row
       dsdT = dst_row
       dsdR = dsd_row

       do n = 1, nspec
          dpdX(n) = dpa_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          dhdX(n) = dEdX(n) + &
               (pres/dens**2 - dEdR)*dpdX(n)/dPdr

       enddo

    else if (input .EQ. eos_input_re) then

!---------------------------------------------------------------------------
! input = 5: dens, energy, and xmass are inputs
!---------------------------------------------------------------------------

!      do_eos_diag = .true.
       ! load the initial guess
       temp_row = temp
       den_row  = dens
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',temp,dens

       ! we want to converge to the given energy
       energy_want = eint

       if (energy_want < ZERO) then
          print *,'BAD HERE ',pt_index(1)
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('ERROR: energy < 0 in the EOS', err_string)
          else
             call bl_error('ERROR: energy < 0 in the EOS')
          endif
       endif

       if (do_eos_diag) print*,'WANT e ',energy_want

       call helmeos(do_coulomb,eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

       do iter = 1, max_newton

          niter = iter

          ! recompute the energy and its temperature derivative
          ener1 = etot_row

          if (do_eos_diag) then
             print*,'ENER1 ',iter,ener1
             print*,'DEDT ',iter,det_row
          end if

          tnew = temp_row - (ener1-energy_want)/det_row

          if (do_eos_diag) then
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  ener1 - energy_want, ' / ', det_row
          endif

          ! don't let the temperature change by more than a factor of two
          tnew = max(.5d0*temp_row, &
                        min(tnew, 2.d0*temp_row))

          ! don't let us freeze
          tnew = max(smallt, tnew)

          if (do_eos_diag) print*,'TNEW AFTER ',iter,tnew

          ! compute the error
          error = 0.0d0
          error = max(error,abs(tnew - temp_row)/temp_row)

          ! store the new temperature
          temp_row = tnew

          if (error .LT. ttol) goto 270
        
          call helmeos(do_coulomb,eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       enddo

       ! Land here if too many iterations are needed

       continue

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations')
       endif

270     continue

       ! store the end result
       temp = tnew
       pres = ptot_row
       enthalpy = eint + ptot_row/dens
       
       c_v = cv_row
       c_p = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       ne   = xne_row + xnp_row
       eta  = etaele_row
       pele = pele_row + ppos_row
       
       dPdR = dpd_row
       dPdT = dpt_row
       dEdR = ded_row
       dEdT = det_row   ! c_v
       gam1 = gam1_row
!      cs =   cs_row
       cs =   sqrt(gam1*pres/dens)
       entropy = stot_row
       dsdT = dst_row
       dsdR = dsd_row

       do n = 1, nspec
          dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          dhdX(n) = dEdX(n) + &
               (pres/dens**2 - dEdR)*dpdX(n)/dPdr

       enddo

    else if (input .EQ. eos_input_ps) then
!---------------------------------------------------------------------------
! input = 6: pres, entropy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = temp
       den_row  = dens
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',temp,dens

       ! we want to converge to the given entropy and pressure
       entropy_want = entropy
       pres_want    = pres

       if (entropy_want < ZERO) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('ERROR: entropy < 0 in the EOS', err_string)
          else
             call bl_error('ERROR: entropy < 0 in the EOS')
          endif
       endif

       if (pres_want < ZERO) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('ERROR: pressure < 0 in the EOS', err_string)
          else
             call bl_error('ERROR: pressure < 0 in the EOS')
          endif
       endif

       if (do_eos_diag) then
          print*,'WANT s ',entropy_want
          print*,'WANT pres', pres_want
       endif

       call helmeos(do_coulomb,eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

       do iter = 1, max_newton

          niter = iter

          ! correct density and temperature
          pres1 = ptot_row
          entr1 = stot_row

          if (do_eos_diag) then
             print*,'PRES1 ',iter,pres1
             print*,'ENTR1 ',iter,entr1
          end if

          ! two functions, f and g, to iterate over
          f = pres_want - pres1
          dfdd = -dpd_row
          dfdt = -dpt_row
          
          g = entropy_want - entr1
          dgdd = -dsd_row
          dgdt = -dst_row
          !
          ! 0 = f + dfdd * deld + dfdt * delt
          ! 0 = g + dgdd * deld + dgdt * delt
          !
          deld = (f*dgdt - g*dfdt) / (dgdd*dfdt - dgdt*dfdd)

          dnew = den_row + deld

          tnew = temp_row - (f + dfdd*deld) / dfdt

          if (do_eos_diag) then
             print *, 'DNEW FIRST ', den_row, ' + ', &
                  f*dgdt - g*dfdt, ' / ', dgdd*dfdt - dgdt*dfdd
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  f + dfdd*deld, ' / ', dfdt
          endif

          ! don't let the temperature or density change by more
          ! than a factor of two
          tnew = max(HALF*temp_row, &
                        min(tnew, TWO*temp_row))
          dnew = max(HALF*den_row, &
                        min(dnew, TWO*den_row))

          ! don't let us freeze or evacuate
          tnew = max(smallt, tnew)
          dnew = max(smalld, dnew)

          if (do_eos_diag) then
             print*,'DNEW AFTER ',iter,dnew
             print*,'TNEW AFTER ',iter,tnew
          endif

          ! compute the errors
          error = ZERO
          error2 = ZERO
          error  = max(error ,abs(dnew - den_row)/den_row)
          error2 = max(error2,abs(tnew - temp_row)/temp_row)

          ! store the new temperature and density
          den_row = dnew
          temp_row = tnew

          if (error .LT. dtol .and. error2 .LT. ttol) goto 370
     
          call helmeos(do_coulomb,eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       enddo

       ! Land here if too many iterations are needed

       continue

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations')
       endif

370     continue

       ! store the end result
       dens = dnew
       temp = tnew
       eint = etot_row
       enthalpy = eint + ptot_row/dens
          
       c_v = cv_row
       c_p = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       ne   = xne_row + xnp_row
       eta  = etaele_row
       pele = pele_row + ppos_row
       
       dPdR = dpd_row
       dPdT = dpt_row
       dEdR = ded_row
       dEdT = det_row   ! c_v
       gam1 = gam1_row
!      cs =   cs_row
       cs =   sqrt(gam1*pres/dens)
       dsdT = dst_row
       dsdR = dsd_row

       do n = 1, nspec
          dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          dhdX(n) = dEdX(n) + &
               (pres/dens**2 - dEdR)*dpdX(n)/dPdr

       enddo

    else 

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: invalid input', err_string)
       else
          call bl_error('EOS: invalid input')
       endif

    endif


    return
  end subroutine eos_old

end module eos_module
