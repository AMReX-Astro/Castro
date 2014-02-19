module eos_type_module
  
  use bl_types
  use network

  implicit none

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  ! rho      -- mass density (g/cm**3)
  ! T        -- temperature (K)
  ! xn       -- the mass fractions of the individual isotopes
  ! p        -- the pressure (dyn/cm**2)
  ! h        -- the enthalpy (erg/g)
  ! e        -- the internal energy (erg/g)
  ! s        -- the entropy (erg/g/K)
  ! c_v      -- specific heat at constant volume
  ! c_p      -- specific heat at constant pressure
  ! ne       -- number density of electrons + positrons
  ! np       -- number density of positrons only
  ! eta      -- degeneracy parameter
  ! pele     -- electron pressure + positron pressure
  ! ppos     -- position pressure only
  ! mu       -- mean molecular weight
  ! mu_e     -- mean number of nucleons per electron
  ! dPdT     -- d pressure/ d temperature
  ! dPdr     -- d pressure/ d density
  ! dedT     -- d energy/ d temperature
  ! dedr     -- d energy/ d density
  ! dsdT     -- d entropy/ d temperature
  ! dsdr     -- d entropy/ d density
  ! dhdT     -- d enthalpy/ d temperature
  ! dhdr     -- d enthalpy/ d density
  ! dPdX     -- d pressure / d xmass
  ! dhdX     -- d enthalpy / d xmass at constant pressure
  ! gam1     -- first adiabatic index (d log P/ d log rho) |_s
  ! cs       -- sound speed
  ! abar     -- average atomic number ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
  ! zbar     -- average proton number ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )

  ! Initialize the main quantities to an unphysical number
  ! so that we know if the user forgot to initialize them
  ! when calling the EOS in a particular mode.

  real (kind=dp_t), parameter :: init_num = -1.0d200

  type eos_t

     double precision :: rho         = init_num
     double precision :: T           = init_num
     double precision :: p           = init_num
     double precision :: e           = init_num
     double precision :: h           = init_num
     double precision :: s           = init_num
     double precision :: dpdT        
     double precision :: dpdr        
     double precision :: dedT        
     double precision :: dedr        
     double precision :: dhdT        
     double precision :: dhdr        
     double precision :: dsdT        
     double precision :: dsdr        
     double precision :: dpde        
     double precision :: dpdr_e      

     double precision :: xn(nspec)   = init_num
     double precision :: cv          
     double precision :: cp          
     double precision :: xne         
     double precision :: xnp         
     double precision :: eta         
     double precision :: pele        
     double precision :: ppos        
     double precision :: mu
     double precision :: mu_e
     double precision :: dedX(nspec) 
     double precision :: dpdX(nspec) 
     double precision :: dhdX(nspec) 
     double precision :: gam1        
     double precision :: cs          

     double precision :: abar        
     double precision :: zbar        
     double precision :: dpa          
     double precision :: dpz         
     double precision :: dea         
     double precision :: dez         

  end type eos_t

contains



  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state, assume_neutral)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(in   ) :: assume_neutral

    double precision :: ymass, ysum, yzsum, ysumi
    integer :: n

    ysum = ZERO
    yzsum = ZERO

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight, and
    ! mu_e, the mean number of nucleons per electron.

    do n = 1, nspec
       ymass = state % xn(n) / aion(n)
       ysumi = ysumi + (ONE + zion(n)) * ymass
       ysum  = ysum  + ymass
       yzsum = yzsum + zion(n) * ymass
    enddo

    state % abar = ONE / ysum
    state % zbar = yzsum * state % abar
    state % mu_e = ONE / yzsum

    if (assume_neutral) then

      state % mu = state % abar

    else

      state % mu = ONE / ysumi

    endif       


  end subroutine composition      



  subroutine composition_derivatives(state, assume_neutral)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(in   ) :: assume_neutral

    double precision :: dmudX
    integer :: n

    ! The species only come into p and e (and therefore h)
    ! through mu, so first compute dmu/dX.
    !
    ! NOTE: an extra, constant term appears in dmudX, which
    ! results from writing mu = sum {X_k} / sum {X_k / A_k}
    ! (for the neutral, analogous for the ionized).  The
    ! numerator is simply 1, but we can differentiate
    ! wrt it, giving the constant mu(k) term in dmudX.  Since
    ! dPdX only appears in a sum over species creation rate 
    ! (omegadot) and sum{omegadot} = 0, this term has no effect.
    ! It is added simply for completeness.

    do n = 1, nspec       

       state % dpdX(n) = state % dpa * (state % abar/aion(n)) &
                       * (aion(n) - state % abar)             &
                       + state % dpz * (state % abar/aion(n)) &
                       * (zion(n) - state % zbar)

       state % dEdX(n) = state % dea * (state % abar/aion(n)) &
                       * (aion(n) - state % abar)             &
                       + state % dez * (state % abar/aion(n)) &
                       * (zion(n) - state % zbar)

       ! dhdX is at constant pressure -- see paper III for details.
       state % dhdX(n) = state % dedX(n) + (state % p / state % rho**2 - state % dedr) &
                                         *  state % dPdX(n) / state % dPdr

    enddo

  end subroutine composition_derivatives

end module eos_type_module
