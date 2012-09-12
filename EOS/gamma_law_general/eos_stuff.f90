! This is a constant gamma equation of state, using and ideal gas.
!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species.
!
!   NOTE: in the helmholtz EOS, we use Abar, which is the weighted
!   ion mass.  Abar does not include any free electron contribution.
!
! There are 2 ways to compute mu, depending on whether the gas is
! completely ionized or not.
!
! For a neutral gas, the mean molecular weight is:
!
!   1/mu = sum_k { X_k / A_k }
!
! For a completely ionized gas, the mean molecular weight is:
!
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).  

module eos_module

  use bl_types
  use bl_space
  use bl_constants_module, only: M_PI, ONE
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

  integer, parameter, public :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter, public :: eos_input_rh = 2   ! density, enthalpy are inputs
  integer, parameter, public :: eos_input_tp = 3   ! temperature, pressure are inputs
  integer, parameter, public :: eos_input_rp = 4   ! density, pressure are inputs
  integer, parameter, public :: eos_input_re = 5   ! density, internal energy are inputs
  integer, parameter, public :: eos_input_ps = 6   ! pressure, entropy are inputs
 

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  real(kind=dp_t), save, public  :: gamma_const

  logical, save, private :: initialized = .false.

  private nspec, aion, zion

  public eos_init, eos_get_small_temp, eos_get_small_dens, eos_given_ReX, &
       eos_e_given_RPX, eos_S_given_ReX, eos_given_RTX, eos_dpdr_given_RTX, &
       eos_given_TPX, eos_given_PSX, eos_get_cv, eos

  interface eos
     module procedure eos_old
     module procedure eos_new
  end interface eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  subroutine eos_init(small_temp, small_dens, gamma_in)

    use extern_probin_module

    implicit none
 
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens
    real(kind=dp_t), intent(in), optional :: gamma_in
 
    ! constant ratio of specific heats
    if (present(gamma_in)) then
       gamma_const = gamma_in
    else
       gamma_const = 5.d0/3.d0
    end if
 
    ! small temperature and density parameters
    if (present(small_temp)) then
       smallt = small_temp
    else
       smallt = 0.d0
    endif
 
    if (present(small_dens)) then
       smalld = small_dens
    else
       smalld = 0.d0
    endif

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
    !       and   dpde is partial p / partial e   at constant rho 


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
              do_diag)

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
              do_diag)

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
              do_diag)

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
              do_diag)

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
             do_diag)

    P  =    p_eos
    e  =    e_eos
    dpdr =  dpdr_eos

  end subroutine eos_dpdr_given_RTX

  subroutine eos_given_TPX(e, P, R, T, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e, R
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
              do_diag)

    R  =    den_eos
    e  =    e_eos

  end subroutine eos_given_TPX

  subroutine eos_given_PSX(P, S, X, R, T, e, pt_index)

    ! In/out variables
    real(kind=dp_t), intent(  out) :: e, R, T
    real(kind=dp_t), intent(in   ) :: P, S, X(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    logical :: do_diag

    do_diag = .false.

    p_eos = P
    s_eos = S
    xn_eos(1:nspec) = X(1:nspec)

    call eos(eos_input_ps, den_eos, temp_eos, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    R = den_eos
    T = temp_eos
    e = e_eos

  end subroutine eos_given_PSX


  subroutine eos_get_cv(cv, R, T, X, pt_index)

! input/output variables
    real(kind=dp_t), intent(out) :: cv
    real(kind=dp_t), intent(in)  :: R, T, X(:)
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
         do_diag)
    
    cv = cv_eos

  end subroutine eos_get_cv


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
    use fundamental_constants_module, only: k_B, n_A, hbar
    use extern_probin_module, only: eos_assume_neutral

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
! dPdX     -- d pressure / d xmass(k)
! dhdX     -- d enthalpy / d xmass(k)  -- AT CONSTANT PRESSURE!!!
! gam1     -- first adiabatic index (d log P/ d log rho) |_s
! cs       -- sound speed -- note that this is the non-relativistic one
!             (we compute it in this wrapper as sqrt(gam1 p /rho) instead
!             of taking the relativistic version from helmeos.
! entropy  -- entropy (erg/g/K)  NOTE: presently the entropy expression is 
!             valid only for an ideal MONATOMIC gas (gamma = 5/3).
!
! input = 1 means dens, temp    , and xmass are inputs, return enthalpy, eint
!       = 2 means dens, enthalpy, and xmass are inputs, return temp    , eint
!                (note, temp should be filled with an initial guess)
!       = 3 means temp, pres    , and xmass are inputs, return dens    , etc
!       = 4 means dens, pres    , and xmass are inputs, return temp    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc
!       = 6 means pres, entr    , and xmass are inputs, return temp    , etc
!
!
! derivatives wrt X_k:
!
!   For an ideal gas, the thermodynamic quantities only depend on composition
!   through the mean molecular weight, mu.
!
!   Using the chain rule:
!
!   dp/dX_k = dp/d(mu) d(mu)/dX_k
!

    implicit none

    logical do_eos_diag
    integer, intent(in) :: input

    real(kind=dp_t) :: dens, temp
    real(kind=dp_t) :: xmass(nspec)
    real(kind=dp_t) :: pres, enthalpy, eint
    real(kind=dp_t) :: c_v, c_p
    real(kind=dp_t) :: ne, eta, pele
    real(kind=dp_t) :: dPdT, dPdR, dedT, dedR
    real(kind=dp_t) :: gam1, entropy, cs
    real(kind=dp_t) :: dPdX(nspec), dedX(nspec), dhdX(nspec)
    real(kind=dp_t) :: dsdT, dsdR

    integer, optional, intent(in   ) :: pt_index(:)


    ! local variables
    real(kind=dp_t) :: ymass(nspec)    
    real(kind=dp_t) :: mu
    real(kind=dp_t) :: dmudX, sum_y

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')
      

    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------
    if (eos_assume_neutral) then
       ! assume completely neutral atoms

       sum_y  = 0.d0
          
       do n = 1, nspec
          ymass(n) = xmass(n)/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
          
       mu = 1.d0/sum_y

    else
       ! assume completely ionized species

       sum_y  = 0.d0
          
       do n = 1, nspec
          ymass(n) = xmass(n)*(1.d0 + zion(n))/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
          
       mu = 1.d0/sum_y

    endif

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------
    if (input .EQ. eos_input_rh) then

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
       temp = (enthalpy*mu*m_nucleon/k_B)*(gamma_const - 1.0_dp_t)/gamma_const


    else if (input .EQ. eos_input_tp ) then

       ! temp, pres, and xmass are inputs
          
       ! Solve for the density:
       ! p = rho k T / (mu m_nucleon)
       dens = pres*mu*m_nucleon/(k_B*temp)


    else if (input .EQ. eos_input_rp ) then

       ! dens, pres, and xmass are inputs

       ! Solve for the temperature:
       ! p = rho k T / (mu m_nucleon)
       temp = pres*mu*m_nucleon/(k_B*dens)


    else if (input .EQ. eos_input_re) then

       ! dens, energy, and xmass are inputs

       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]
       temp = eint*mu*m_nucleon*(gamma_const-1.0_dp_t)/k_B


    else if (input .EQ. eos_input_ps) then
       
       ! pressure and entropy are inputs

       ! Solve for the temperature
       ! Invert Sackur-Tetrode eqn (below) using 
       ! rho = p mu m_nucleon / (k T)
       temp = pres**(2.0_dp_t/5.0_dp_t) * &
            ( 2.0_dp_t*M_PI*hbar*hbar/(mu*m_nucleon) )**(3.0_dp_t/5.0_dp_t) * &
            dexp(2.0_dp_t*mu*m_nucleon*entropy/(5.0_dp_t*k_B) - 1.0_dp_t) / &
            k_B

       ! Solve for the density
       ! rho = p mu m_nucleon / (k T)
       dens = pres*mu*m_nucleon/(k_B*temp)

    endif

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    ! compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation
    pres = dens*k_B*temp/(mu*m_nucleon)
    eint = pres/(gamma_const - 1.0_dp_t)/dens

    ! enthalpy is h = e + p/rho
    enthalpy = eint + pres/dens

    ! entropy (per gram) of an ideal monoatomic gas (the Sactur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.
    entropy = (k_B/(mu*m_nucleon))*(2.5_dp_t + &
         log( ( (mu*m_nucleon)**2.5/dens )*(k_B*temp)**1.5_dp_t / (2.0_dp_t*M_PI*hbar*hbar)**1.5_dp_t ) )

    ! compute the thermodynamic derivatives and specific heats 
    dPdT = pres/temp
    dPdR = pres/dens
    dedT = eint/temp
    dedR = 0.d0
    dsdT = 0.d0
    dsdR = 0.d0

    c_v = dedT
    c_p = gamma_const*c_v

    gam1 = gamma_const

    do n = 1, nspec

       ! the species only come into p and e (and therefore h)
       ! through mu, so first compute dmu/dX
       !
       ! NOTE: an extra, constant term appears in dmudx, which
       ! results from writing mu = sum { X_k} / sum {X_k / A_k}
       ! (for the neutral, analogous for the ionized).  The
       ! numerator is simply 1, but we can differentiate
       ! wrt it, giving the constant mu(k) term in dmudx.  Since
       ! dPdX only appears in a sum over species creation rate 
       ! (omegadot) and sum{omegadot} = 0, this term has no effect.
       ! If is added simply for completeness.

       if (eos_assume_neutral) then
          dmudX =  (mu/aion(n))*(aion(n) - mu)
       else
          dmudX =  (mu/aion(n))*(aion(n) - mu*(1.0_dp_t + zion(n)))
       endif

       dPdX(n) = -(pres/mu)*dmudX
       dedX(n) = -(eint/mu)*dmudX
          
       ! dhdX is at constant pressure -- see paper III for details
       dhdX(n) = dedX(n) + &
            (pres/dens**2 - dedR)*dPdX(n)/dPdr
    enddo

    ! electron-specific stuff (really for the degenerate EOS)
    ne = 0.d0
    eta = 0.d0
    pele = 0.d0

    ! sound speed
    cs = sqrt(gamma_const*pres/dens)

    return
  end subroutine eos_old

end module eos_module
