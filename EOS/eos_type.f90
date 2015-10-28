module eos_type_module
  
  use bl_types
  use network
  use eos_data_module
  use mempool_module
  use bl_constants_module
  
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
  ! y_e      -- electron fraction == 1 / mu_e
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
  ! dpdA     -- d pressure/ d abar
  ! dpdZ     -- d pressure/ d zbar
  ! dedA     -- d energy/ d abar
  ! dedZ     -- d energy/ d zbar

  ! Initialize the main quantities to an unphysical number
  ! so that we know if the user forgot to initialize them
  ! when calling the EOS in a particular mode.

  double precision, parameter :: init_num  = -1.0d200
  double precision, parameter :: init_test = -1.0d199

  type eos_type

     integer :: N, width(1), spec_width(2), aux_width(2)
     
  end type eos_type
  
  type, extends(eos_type) :: eos_t

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
    double precision :: aux(naux)   = init_num
    double precision :: cv
    double precision :: cp
    double precision :: xne
    double precision :: xnp
    double precision :: eta
    double precision :: pele
    double precision :: ppos
    double precision :: mu
    double precision :: mu_e
    double precision :: y_e
    double precision :: dedX(nspec)
    double precision :: dpdX(nspec)
    double precision :: dhdX(nspec)
    double precision :: gam1
    double precision :: cs

    double precision :: abar
    double precision :: zbar
    double precision :: dpdA

    double precision :: dpdZ
    double precision :: dedA
    double precision :: dedZ

    logical :: reset
    
  end type eos_t



  type, extends(eos_type) :: eos_t_vector

    double precision, pointer :: rho(:)
    double precision, pointer :: T(:)
    double precision, pointer :: p(:)
    double precision, pointer :: e(:)
    double precision, pointer :: h(:)
    double precision, pointer :: s(:)
    double precision, pointer :: dpdT(:)
    double precision, pointer :: dpdr(:)
    double precision, pointer :: dedT(:)
    double precision, pointer :: dedr(:)
    double precision, pointer :: dhdT(:)
    double precision, pointer :: dhdr(:)
    double precision, pointer :: dsdT(:)
    double precision, pointer :: dsdr(:)
    double precision, pointer :: dpde(:)
    double precision, pointer :: dpdr_e(:)

    double precision, pointer :: xn(:,:)
    double precision, pointer :: aux(:,:)
    double precision, pointer :: cv(:)
    double precision, pointer :: cp(:)
    double precision, pointer :: xne(:)
    double precision, pointer :: xnp(:)
    double precision, pointer :: eta(:)
    double precision, pointer :: pele(:)
    double precision, pointer :: ppos(:)
    double precision, pointer :: mu(:)
    double precision, pointer :: mu_e(:)
    double precision, pointer :: y_e(:)
    double precision, pointer :: dedX(:,:)
    double precision, pointer :: dpdX(:,:)
    double precision, pointer :: dhdX(:,:)
    double precision, pointer :: gam1(:)
    double precision, pointer :: cs(:)

    double precision, pointer :: abar(:)
    double precision, pointer :: zbar(:)
    double precision, pointer :: dpdA(:)
    double precision, pointer :: dpdZ(:)
    double precision, pointer :: dedA(:)
    double precision, pointer :: dedZ(:)

    logical :: reset
    
  end type eos_t_vector

  
 
  type, extends(eos_type) :: eos_t_1D

    integer :: lo(1), hi(1)
    
    double precision, pointer :: rho(:)
    double precision, pointer :: T(:)
    double precision, pointer :: p(:)
    double precision, pointer :: e(:)
    double precision, pointer :: h(:)
    double precision, pointer :: s(:)
    double precision, pointer :: dpdT(:)
    double precision, pointer :: dpdr(:)
    double precision, pointer :: dedT(:)
    double precision, pointer :: dedr(:)
    double precision, pointer :: dhdT(:)        
    double precision, pointer :: dhdr(:)        
    double precision, pointer :: dsdT(:)        
    double precision, pointer :: dsdr(:)        
    double precision, pointer :: dpde(:)        
    double precision, pointer :: dpdr_e(:)      

    double precision, pointer :: xn(:,:)
    double precision, pointer :: aux(:,:)
    double precision, pointer :: cv(:)          
    double precision, pointer :: cp(:)          
    double precision, pointer :: xne(:)         
    double precision, pointer :: xnp(:)         
    double precision, pointer :: eta(:)         
    double precision, pointer :: pele(:)        
    double precision, pointer :: ppos(:)        
    double precision, pointer :: mu(:)
    double precision, pointer :: mu_e(:)
    double precision, pointer :: y_e(:)
    double precision, pointer :: dedX(:,:) 
    double precision, pointer :: dpdX(:,:) 
    double precision, pointer :: dhdX(:,:) 
    double precision, pointer :: gam1(:)        
    double precision, pointer :: cs(:)          

    double precision, pointer :: abar(:)        
    double precision, pointer :: zbar(:)        
    double precision, pointer :: dpdA(:)          
    double precision, pointer :: dpdZ(:)        
    double precision, pointer :: dedA(:)         
    double precision, pointer :: dedZ(:)         

    logical :: reset
    
  end type eos_t_1D

 
  
  type, extends(eos_type) :: eos_t_2D

    integer :: lo(2), hi(2)
    
    double precision, pointer :: rho(:,:)
    double precision, pointer :: T(:,:)
    double precision, pointer :: p(:,:)
    double precision, pointer :: e(:,:)
    double precision, pointer :: h(:,:)
    double precision, pointer :: s(:,:)
    double precision, pointer :: dpdT(:,:)
    double precision, pointer :: dpdr(:,:)
    double precision, pointer :: dedT(:,:)
    double precision, pointer :: dedr(:,:)
    double precision, pointer :: dhdT(:,:)        
    double precision, pointer :: dhdr(:,:)        
    double precision, pointer :: dsdT(:,:)        
    double precision, pointer :: dsdr(:,:)        
    double precision, pointer :: dpde(:,:)        
    double precision, pointer :: dpdr_e(:,:)      

    double precision, pointer :: xn(:,:,:)
    double precision, pointer :: aux(:,:,:)
    double precision, pointer :: cv(:,:)          
    double precision, pointer :: cp(:,:)          
    double precision, pointer :: xne(:,:)         
    double precision, pointer :: xnp(:,:)         
    double precision, pointer :: eta(:,:)         
    double precision, pointer :: pele(:,:)        
    double precision, pointer :: ppos(:,:)        
    double precision, pointer :: mu(:,:)
    double precision, pointer :: mu_e(:,:)
    double precision, pointer :: y_e(:,:)
    double precision, pointer :: dedX(:,:,:) 
    double precision, pointer :: dpdX(:,:,:) 
    double precision, pointer :: dhdX(:,:,:) 
    double precision, pointer :: gam1(:,:)        
    double precision, pointer :: cs(:,:)          

    double precision, pointer :: abar(:,:)        
    double precision, pointer :: zbar(:,:)        
    double precision, pointer :: dpdA(:,:)          
    double precision, pointer :: dpdZ(:,:)        
    double precision, pointer :: dedA(:,:)         
    double precision, pointer :: dedZ(:,:)         

    logical :: reset
    
  end type eos_t_2D


  
  type, extends(eos_type) :: eos_t_3D

    integer :: lo(3), hi(3)
     
    double precision, pointer :: rho(:,:,:)
    double precision, pointer :: T(:,:,:)
    double precision, pointer :: p(:,:,:)
    double precision, pointer :: e(:,:,:)
    double precision, pointer :: h(:,:,:)
    double precision, pointer :: s(:,:,:)
    double precision, pointer :: dpdT(:,:,:)
    double precision, pointer :: dpdr(:,:,:)
    double precision, pointer :: dedT(:,:,:)
    double precision, pointer :: dedr(:,:,:)
    double precision, pointer :: dhdT(:,:,:)        
    double precision, pointer :: dhdr(:,:,:)        
    double precision, pointer :: dsdT(:,:,:)        
    double precision, pointer :: dsdr(:,:,:)        
    double precision, pointer :: dpde(:,:,:)        
    double precision, pointer :: dpdr_e(:,:,:)      

    double precision, pointer :: xn(:,:,:,:)
    double precision, pointer :: aux(:,:,:,:)
    double precision, pointer :: cv(:,:,:)          
    double precision, pointer :: cp(:,:,:)          
    double precision, pointer :: xne(:,:,:)         
    double precision, pointer :: xnp(:,:,:)         
    double precision, pointer :: eta(:,:,:)         
    double precision, pointer :: pele(:,:,:)        
    double precision, pointer :: ppos(:,:,:)        
    double precision, pointer :: mu(:,:,:)
    double precision, pointer :: mu_e(:,:,:)
    double precision, pointer :: y_e(:,:,:)
    double precision, pointer :: dedX(:,:,:,:) 
    double precision, pointer :: dpdX(:,:,:,:) 
    double precision, pointer :: dhdX(:,:,:,:) 
    double precision, pointer :: gam1(:,:,:)        
    double precision, pointer :: cs(:,:,:)          

    double precision, pointer :: abar(:,:,:)        
    double precision, pointer :: zbar(:,:,:)        
    double precision, pointer :: dpdA(:,:,:)          
    double precision, pointer :: dpdZ(:,:,:)        
    double precision, pointer :: dedA(:,:,:)         
    double precision, pointer :: dedZ(:,:,:)         

    logical :: reset
    
  end type eos_t_3D

contains

  subroutine eos_allocate(state, lo, hi)

    implicit none

    class (eos_type) :: state
    integer          :: lo(:), hi(:)

    integer :: ndim

    if (size(hi) .ne. size(lo)) then
       call bl_error("Error: cannot have different lo and hi dimensions in EOS type constructor.")
    endif

    ndim = size(lo)
    
    select type (state)

    type is (eos_t_1D)

       if (ndim .ne. 1) then
          call bl_error("Error: dimensionality of lo and hi must match EOS dimensionality.")
       endif

       state % lo = lo
       state % hi = hi
       state % reset = .false.
       
       call bl_allocate(state % rho, lo(1), hi(1))
       call bl_allocate(state % T, lo(1), hi(1))
       call bl_allocate(state % p, lo(1), hi(1))
       call bl_allocate(state % e, lo(1), hi(1))
       call bl_allocate(state % h, lo(1), hi(1))
       call bl_allocate(state % s, lo(1), hi(1))
       call bl_allocate(state % dpdT, lo(1), hi(1))
       call bl_allocate(state % dpdr, lo(1), hi(1))
       call bl_allocate(state % dedT, lo(1), hi(1))
       call bl_allocate(state % dedr, lo(1), hi(1))
       call bl_allocate(state % dhdT, lo(1), hi(1))
       call bl_allocate(state % dhdr, lo(1), hi(1))
       call bl_allocate(state % dsdT, lo(1), hi(1))
       call bl_allocate(state % dsdr, lo(1), hi(1))
       call bl_allocate(state % dpde, lo(1), hi(1))
       call bl_allocate(state % dpdr_e, lo(1), hi(1))
       call bl_allocate(state % xn, lo(1), hi(1), 1, nspec)
       call bl_allocate(state % aux, lo(1), hi(1), 1, naux)
       call bl_allocate(state % cv, lo(1), hi(1))
       call bl_allocate(state % cp, lo(1), hi(1))
       call bl_allocate(state % xne, lo(1), hi(1))
       call bl_allocate(state % xnp, lo(1), hi(1))
       call bl_allocate(state % eta, lo(1), hi(1))
       call bl_allocate(state % pele, lo(1), hi(1))
       call bl_allocate(state % ppos, lo(1), hi(1))
       call bl_allocate(state % mu, lo(1), hi(1))
       call bl_allocate(state % mu_e, lo(1), hi(1))
       call bl_allocate(state % y_e, lo(1), hi(1))
       call bl_allocate(state % dedX, lo(1), hi(1), 1, nspec)
       call bl_allocate(state % dpdX, lo(1), hi(1), 1, nspec)
       call bl_allocate(state % dhdX, lo(1), hi(1), 1, nspec)
       call bl_allocate(state % gam1, lo(1), hi(1))
       call bl_allocate(state % cs, lo(1), hi(1))
       call bl_allocate(state % abar, lo(1), hi(1))
       call bl_allocate(state % zbar, lo(1), hi(1))
       call bl_allocate(state % dpdA, lo(1), hi(1))
       call bl_allocate(state % dpdZ, lo(1), hi(1))
       call bl_allocate(state % dedA, lo(1), hi(1))
       call bl_allocate(state % dedZ, lo(1), hi(1))

       state % rho(:)   = init_num
       state % T(:)     = init_num
       state % p(:)     = init_num
       state % e(:)     = init_num
       state % h(:)     = init_num
       state % s(:)     = init_num
       state % xn(:,:)  = init_num
       state % aux(:,:) = init_num

    type is (eos_t_2D)

       if (ndim .ne. 2) then
          call bl_error("Error: dimensionality of lo and hi must match EOS dimensionality.")
       endif       

       state % lo = lo
       state % hi = hi
       state % reset = .false.
       
       call bl_allocate(state % rho, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % T, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % p, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % e, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % h, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % s, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dpdT, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dpdr, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dedT, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dedr, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dhdT, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dhdr, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dsdT, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dsdr, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dpde, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dpdr_e, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % xn, lo(1), hi(1), lo(2), hi(2), 1, nspec)
       call bl_allocate(state % aux, lo(1), hi(1), lo(2), hi(2), 1, naux)
       call bl_allocate(state % cv, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % cp, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % xne, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % xnp, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % eta, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % pele, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % ppos, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % mu, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % mu_e, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % y_e, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dedX, lo(1), hi(1), lo(2), hi(2), 1, nspec)
       call bl_allocate(state % dpdX, lo(1), hi(1), lo(2), hi(2), 1, nspec)
       call bl_allocate(state % dhdX, lo(1), hi(1), lo(2), hi(2), 1, nspec)
       call bl_allocate(state % gam1, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % cs, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % abar, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % zbar, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dpdA, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dpdZ, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dedA, lo(1), hi(1), lo(2), hi(2))
       call bl_allocate(state % dedZ, lo(1), hi(1), lo(2), hi(2))

       state % rho(:,:)   = init_num
       state % T(:,:)     = init_num
       state % p(:,:)     = init_num
       state % e(:,:)     = init_num
       state % h(:,:)     = init_num
       state % s(:,:)     = init_num
       state % xn(:,:,:)  = init_num
       state % aux(:,:,:) = init_num

    type is (eos_t_3D)

       if (ndim .ne. 3) then
          call bl_error("Error: dimensionality of lo and hi must match EOS dimensionality.")
       endif       

       state % lo = lo
       state % hi = hi
       state % reset = .false.
       
       call bl_allocate(state % rho, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % T, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % p, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % e, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % h, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % s, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dpdT, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dpdr, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dedT, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dedr, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dhdT, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dhdr, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dsdT, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dsdr, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dpde, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dpdr_e, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % xn, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3), 1, nspec)
       call bl_allocate(state % aux, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3), 1, naux)
       call bl_allocate(state % cv, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % cp, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % xne, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % xnp, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % eta, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % pele, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % ppos, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % mu, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % mu_e, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % y_e, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dedX, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3), 1, nspec)
       call bl_allocate(state % dpdX, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3), 1, nspec)
       call bl_allocate(state % dhdX, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3), 1, nspec)
       call bl_allocate(state % gam1, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % cs, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % abar, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % zbar, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dpdA, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dpdZ, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dedA, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))
       call bl_allocate(state % dedZ, lo(1), hi(1), lo(2), hi(2), lo(3), hi(3))

       state % rho(:,:,:)   = init_num
       state % T(:,:,:)     = init_num
       state % p(:,:,:)     = init_num
       state % e(:,:,:)     = init_num
       state % h(:,:,:)     = init_num
       state % s(:,:,:)     = init_num
       state % xn(:,:,:,:)  = init_num
       state % aux(:,:,:,:) = init_num

    end select

    state % N = product(hi - lo + 1)
    state % width(1) = state % N
    state % spec_width(1) = state % N
    state % spec_width(2) = nspec
    state % aux_width(1) = state % N
    state % aux_width(2) = naux    
    
  end subroutine eos_allocate


  
  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i
    
    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

    do i = 1, state % N

       state % mu_e(i) = ONE / (sum(state % xn(i,:) * zion(:) / aion(:)))       
       state % y_e(i) = ONE / state % mu_e(i)
       
       state % abar(i) = ONE / (sum(state % xn(i,:) / aion(:)))
       state % zbar(i) = state % abar(i) / state % mu_e(i)

    enddo

  end subroutine composition

  

  ! Normalize the mass fractions: they must be individually positive and less than one,
  ! and they must all sum to unity.

  subroutine normalize_abundances(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i    

    do i = 1, state % N

       state % xn(i,:) = max(smallx, min(ONE, state % xn(i,:)))

       state % xn(i,:) = state % xn(i,:) / sum(state % xn(i,:))

    enddo
  
  end subroutine normalize_abundances



  ! Given an index i, return a scalar eos_t type corresponding to that index.
  
  subroutine get_eos_t(state_vector, state, i)

    implicit none

    type (eos_t)        :: state
    type (eos_t_vector) :: state_vector
    integer             :: i

    if (i < 1) then
       call bl_error("Error: cannot access index < 1 in the EOS.")
    else if (i > state_vector % N) then
       call bl_error("Error: cannot access index > length of EOS vector.")
    endif

    state % rho = state_vector % rho(i)
    state % T   = state_vector % T(i)
    state % p   = state_vector % p(i)
    state % h   = state_vector % h(i)
    state % s   = state_vector % s(i)
    state % dpdT = state_vector % dpdT(i)
    state % dpdr = state_vector % dpdr(i)
    state % rho = state_vector % rho(i)
    state % T   = state_vector % T(i)
    state % p   = state_vector % p(i)
    state % e   = state_vector % e(i)
    state % h   = state_vector % h(i)
    state % s   = state_vector % s(i)
    state % dpdT = state_vector % dpdT(i)
    state % dpdr = state_vector % dpdr(i)
    state % dedT = state_vector % dedT(i)
    state % dedr = state_vector % dedr(i)
    state % dhdT = state_vector % dhdT(i)
    state % dhdr = state_vector % dhdr(i)
    state % dsdT = state_vector % dsdT(i)
    state % dsdr = state_vector % dsdr(i)
    state % dpde = state_vector % dpde(i)
    state % dpdr_e = state_vector % dpdr_e(i)
    state % xn(:) = state_vector % xn(i,:)
    if (naux > 0) then
       state % aux(:) = state_vector % aux(i,:)
    endif
    state % cv = state_vector % cv(i)
    state % cp = state_vector % cp(i)
    state % xne = state_vector % xne(i)
    state % xnp = state_vector % xnp(i)
    state % eta = state_vector % eta(i)
    state % pele = state_vector % pele(i)
    state % ppos = state_vector % ppos(i)
    state % mu = state_vector % mu(i)
    state % mu_e = state_vector % mu(i)
    state % y_e = state_vector % y_e(i)
    state % dedX(:) = state_vector % dedX(i,:)
    state % dpdX(:) = state_vector % dpdX(i,:)
    state % dhdX(:) = state_vector % dhdX(i,:)
    state % gam1 = state_vector % gam1(i)
    state % cs = state_vector % cs(i)
    state % abar = state_vector % abar(i)
    state % zbar = state_vector % zbar(i)
    state % dpdA = state_vector % dpdA(i)
    state % dpdZ = state_vector % dpdZ(i)
    state % dedA = state_vector % dedA(i)
    state % dedZ = state_vector % dedZ(i)
    state % reset = state_vector % reset
    
  end subroutine get_eos_t




  ! Given an index i, merge a scalar eos_t type corresponding to that index.
  
  subroutine put_eos_t(state_vector, state, i)

    implicit none

    type (eos_t)        :: state
    type (eos_t_vector) :: state_vector
    integer             :: i

    if (i < 1) then
       call bl_error("Error: cannot access index < 1 in the EOS.")
    else if (i > state_vector % N) then
       call bl_error("Error: cannot access index > length of EOS vector.")
    endif

    state_vector % rho(i) = state % rho
    state_vector % T(i)   = state % T
    state_vector % p(i)   = state % p
    state_vector % e(i)   = state % e
    state_vector % h(i)   = state % h
    state_vector % s(i)   = state % s
    state_vector % dpdT(i) = state % dpdT
    state_vector % dpdr(i) = state % dpdr
    state_vector % dedT(i) = state % dedT
    state_vector % dedr(i) = state % dedr
    state_vector % dhdT(i) = state % dhdT
    state_vector % dhdr(i) = state % dhdr
    state_vector % dsdT(i) = state % dsdT
    state_vector % dsdr(i) = state % dsdr
    state_vector % dpde(i) = state % dpde
    state_vector % dpdr_e(i) = state % dpdr_e
    state_vector % xn(i,:) = state % xn(:)
    if (naux > 0) then
       state_vector % aux(i,:) = state % aux
    endif
    state_vector % cv(i) = state % cv
    state_vector % cp(i) = state % cp
    state_vector % xne(i) = state % xne
    state_vector % xnp(i) = state % xnp
    state_vector % eta(i) = state % eta
    state_vector % pele(i) = state % pele
    state_vector % ppos(i) = state % ppos
    state_vector % mu(i) = state % mu
    state_vector % mu_e(i) = state % mu
    state_vector % y_e(i) = state % y_e
    state_vector % dedX(i,:) = state % dedX(:)
    state_vector % dpdX(i,:) = state % dpdX(:)
    state_vector % dhdX(i,:) = state % dhdX(:)
    state_vector % gam1(i) = state % gam1
    state_vector % cs(i) = state % cs
    state_vector % abar(i) = state % abar
    state_vector % zbar(i) = state % zbar
    state_vector % dpdA(i) = state % dpdA
    state_vector % dpdZ(i) = state % dpdZ
    state_vector % dedA(i) = state % dedA
    state_vector % dedZ(i) = state % dedZ
    state_vector % reset = state % reset
    
  end subroutine put_eos_t
  

  recursive subroutine eos_vector_in(state, state_in)

    use iso_c_binding
    
    implicit none

    class (eos_type),    intent(in   ) :: state_in
    type (eos_t_vector), intent(inout) :: state

    integer :: lo(3), hi(3)

    select type (state_in)

    type is (eos_t_vector)

       state = state_in
       
    type is (eos_t)

       lo(1) = 1
       hi(1) = 1
       
       state % N = 1
       state % width(1) = state % N
       state % spec_width(1) = state % N
       state % spec_width(2) = nspec
       state % aux_width(1) = state % N
       state % aux_width(2) = naux       
       
       call c_f_pointer(c_loc(state_in % rho), state % rho, state % width)
       call c_f_pointer(c_loc(state_in % T), state % T, state % width)
       call c_f_pointer(c_loc(state_in % p), state % p, state % width)
       call c_f_pointer(c_loc(state_in % e), state % e, state % width)
       call c_f_pointer(c_loc(state_in % h), state % h, state % width)
       call c_f_pointer(c_loc(state_in % s), state % s, state % width)
       call c_f_pointer(c_loc(state_in % dpdT), state % dpdT, state % width)
       call c_f_pointer(c_loc(state_in % dpdr), state % dpdr, state % width)
       call c_f_pointer(c_loc(state_in % dedT), state % dedT, state % width)
       call c_f_pointer(c_loc(state_in % dedr), state % dedr, state % width)
       call c_f_pointer(c_loc(state_in % dhdT), state % dhdT, state % width)
       call c_f_pointer(c_loc(state_in % dhdr), state % dhdr, state % width)
       call c_f_pointer(c_loc(state_in % dsdT), state % dsdT, state % width)
       call c_f_pointer(c_loc(state_in % dsdr), state % dsdr, state % width)
       call c_f_pointer(c_loc(state_in % dpde), state % dpde, state % width)
       call c_f_pointer(c_loc(state_in % dpdr_e), state % dpdr_e, state % width)
       call c_f_pointer(c_loc(state_in % xn(1)), state % xn, state % spec_width)
       if (naux > 0) then
          call c_f_pointer(c_loc(state_in % aux(1)), state % aux, state % aux_width)
       endif
       call c_f_pointer(c_loc(state_in % cv), state % cv, state % width)
       call c_f_pointer(c_loc(state_in % cp), state % cp, state % width)
       call c_f_pointer(c_loc(state_in % xne), state % xne, state % width)
       call c_f_pointer(c_loc(state_in % xnp), state % xnp, state % width)
       call c_f_pointer(c_loc(state_in % eta), state % eta, state % width)
       call c_f_pointer(c_loc(state_in % pele), state % pele, state % width)
       call c_f_pointer(c_loc(state_in % ppos), state % ppos, state % width)
       call c_f_pointer(c_loc(state_in % mu), state % mu, state % width)
       call c_f_pointer(c_loc(state_in % mu_e), state % mu_e, state % width)
       call c_f_pointer(c_loc(state_in % y_e), state % y_e, state % width)
       call c_f_pointer(c_loc(state_in % dedX(1)), state % dedX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dpdX(1)), state % dpdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dhdX(1)), state % dhdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % gam1), state % gam1, state % width)
       call c_f_pointer(c_loc(state_in % cs), state % cs, state % width)
       call c_f_pointer(c_loc(state_in % abar), state % abar, state % width)
       call c_f_pointer(c_loc(state_in % zbar), state % zbar, state % width)
       call c_f_pointer(c_loc(state_in % dpdA), state % dpdA, state % width)
       call c_f_pointer(c_loc(state_in % dpdZ), state % dpdZ, state % width)
       call c_f_pointer(c_loc(state_in % dedA), state % dedA, state % width)
       call c_f_pointer(c_loc(state_in % dedZ), state % dedZ, state % width)

       state % reset = state_in % reset       
       
    type is (eos_t_1D)

       lo(1) = state_in % lo(1)
       hi(1) = state_in % hi(1)

       state % N = state_in % N
       state % width(1) = state % N
       state % spec_width(1) = state % N
       state % spec_width(2) = nspec
       state % aux_width(1) = state % N
       state % aux_width(2) = naux
    
       call c_f_pointer(c_loc(state_in % rho(lo(1))), state % rho, state % width)
       call c_f_pointer(c_loc(state_in % T(lo(1))), state % T, state % width)
       call c_f_pointer(c_loc(state_in % p(lo(1))), state % p, state % width)
       call c_f_pointer(c_loc(state_in % e(lo(1))), state % e, state % width)
       call c_f_pointer(c_loc(state_in % h(lo(1))), state % h, state % width)
       call c_f_pointer(c_loc(state_in % s(lo(1))), state % s, state % width)
       call c_f_pointer(c_loc(state_in % dpdT(lo(1))), state % dpdT, state % width)
       call c_f_pointer(c_loc(state_in % dpdr(lo(1))), state % dpdr, state % width)
       call c_f_pointer(c_loc(state_in % dedT(lo(1))), state % dedT, state % width)
       call c_f_pointer(c_loc(state_in % dedr(lo(1))), state % dedr, state % width)
       call c_f_pointer(c_loc(state_in % dhdT(lo(1))), state % dhdT, state % width)
       call c_f_pointer(c_loc(state_in % dhdr(lo(1))), state % dhdr, state % width)
       call c_f_pointer(c_loc(state_in % dsdT(lo(1))), state % dsdT, state % width)
       call c_f_pointer(c_loc(state_in % dsdr(lo(1))), state % dsdr, state % width)
       call c_f_pointer(c_loc(state_in % dpde(lo(1))), state % dpde, state % width)
       call c_f_pointer(c_loc(state_in % dpdr_e(lo(1))), state % dpdr_e, state % width)
       call c_f_pointer(c_loc(state_in % xn(lo(1),1)), state % xn, state % spec_width)
       if (naux > 0) then
          call c_f_pointer(c_loc(state_in % aux(lo(1),1)), state % aux, state % aux_width)
       endif
       call c_f_pointer(c_loc(state_in % cv(lo(1))), state % cv, state % width)
       call c_f_pointer(c_loc(state_in % cp(lo(1))), state % cp, state % width)
       call c_f_pointer(c_loc(state_in % xne(lo(1))), state % xne, state % width)
       call c_f_pointer(c_loc(state_in % xnp(lo(1))), state % xnp, state % width)
       call c_f_pointer(c_loc(state_in % eta(lo(1))), state % eta, state % width)
       call c_f_pointer(c_loc(state_in % pele(lo(1))), state % pele, state % width)
       call c_f_pointer(c_loc(state_in % ppos(lo(1))), state % ppos, state % width)
       call c_f_pointer(c_loc(state_in % mu(lo(1))), state % mu, state % width)
       call c_f_pointer(c_loc(state_in % mu_e(lo(1))), state % mu_e, state % width)
       call c_f_pointer(c_loc(state_in % y_e(lo(1))), state % y_e, state % width)
       call c_f_pointer(c_loc(state_in % dedX(lo(1),1)), state % dedX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dpdX(lo(1),1)), state % dpdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dhdX(lo(1),1)), state % dhdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % gam1(lo(1))), state % gam1, state % width)
       call c_f_pointer(c_loc(state_in % cs(lo(1))), state % cs, state % width)
       call c_f_pointer(c_loc(state_in % abar(lo(1))), state % abar, state % width)
       call c_f_pointer(c_loc(state_in % zbar(lo(1))), state % zbar, state % width)
       call c_f_pointer(c_loc(state_in % dpdA(lo(1))), state % dpdA, state % width)
       call c_f_pointer(c_loc(state_in % dpdZ(lo(1))), state % dpdZ, state % width)
       call c_f_pointer(c_loc(state_in % dedA(lo(1))), state % dedA, state % width)
       call c_f_pointer(c_loc(state_in % dedZ(lo(1))), state % dedZ, state % width)

       state % reset = state_in % reset
       
    type is (eos_t_2D)
    
       lo(1:2) = state_in % lo
       hi(1:2) = state_in % hi

       state % N = state_in % N
       state % width(1) = state % N
       state % spec_width(1) = state % N
       state % spec_width(2) = nspec
       state % aux_width(1) = state % N
       state % aux_width(2) = naux
       
       call c_f_pointer(c_loc(state_in % rho(lo(1),lo(2))), state % rho, state % width)
       call c_f_pointer(c_loc(state_in % T(lo(1),lo(2))), state % T, state % width)
       call c_f_pointer(c_loc(state_in % p(lo(1),lo(2))), state % p, state % width)
       call c_f_pointer(c_loc(state_in % e(lo(1),lo(2))), state % e, state % width)
       call c_f_pointer(c_loc(state_in % h(lo(1),lo(2))), state % h, state % width)
       call c_f_pointer(c_loc(state_in % s(lo(1),lo(2))), state % s, state % width)
       call c_f_pointer(c_loc(state_in % dpdT(lo(1),lo(2))), state % dpdT, state % width)
       call c_f_pointer(c_loc(state_in % dpdr(lo(1),lo(2))), state % dpdr, state % width)
       call c_f_pointer(c_loc(state_in % dedT(lo(1),lo(2))), state % dedT, state % width)
       call c_f_pointer(c_loc(state_in % dedr(lo(1),lo(2))), state % dedr, state % width)
       call c_f_pointer(c_loc(state_in % dhdT(lo(1),lo(2))), state % dhdT, state % width)
       call c_f_pointer(c_loc(state_in % dhdr(lo(1),lo(2))), state % dhdr, state % width)
       call c_f_pointer(c_loc(state_in % dsdT(lo(1),lo(2))), state % dsdT, state % width)
       call c_f_pointer(c_loc(state_in % dsdr(lo(1),lo(2))), state % dsdr, state % width)
       call c_f_pointer(c_loc(state_in % dpde(lo(1),lo(2))), state % dpde, state % width)
       call c_f_pointer(c_loc(state_in % dpdr_e(lo(1),lo(2))), state % dpdr_e, state % width)
       call c_f_pointer(c_loc(state_in % xn(lo(1),lo(2),1)), state % xn, state % spec_width)
       if (naux > 0) then
          call c_f_pointer(c_loc(state_in % aux(lo(1),lo(2),1)), state % aux, state % aux_width)
       endif
       call c_f_pointer(c_loc(state_in % cv(lo(1),lo(2))), state % cv, state % width)
       call c_f_pointer(c_loc(state_in % cp(lo(1),lo(2))), state % cp, state % width)
       call c_f_pointer(c_loc(state_in % xne(lo(1),lo(2))), state % xne, state % width)
       call c_f_pointer(c_loc(state_in % xnp(lo(1),lo(2))), state % xnp, state % width)
       call c_f_pointer(c_loc(state_in % eta(lo(1),lo(2))), state % eta, state % width)
       call c_f_pointer(c_loc(state_in % pele(lo(1),lo(2))), state % pele, state % width)
       call c_f_pointer(c_loc(state_in % ppos(lo(1),lo(2))), state % ppos, state % width)
       call c_f_pointer(c_loc(state_in % mu(lo(1),lo(2))), state % mu, state % width)
       call c_f_pointer(c_loc(state_in % mu_e(lo(1),lo(2))), state % mu_e, state % width)
       call c_f_pointer(c_loc(state_in % y_e(lo(1),lo(2))), state % y_e, state % width)
       call c_f_pointer(c_loc(state_in % dedX(lo(1),lo(2),1)), state % dedX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dpdX(lo(1),lo(2),1)), state % dpdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dhdX(lo(1),lo(2),1)), state % dhdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % gam1(lo(1),lo(2))), state % gam1, state % width)
       call c_f_pointer(c_loc(state_in % cs(lo(1),lo(2))), state % cs, state % width)
       call c_f_pointer(c_loc(state_in % abar(lo(1),lo(2))), state % abar, state % width)
       call c_f_pointer(c_loc(state_in % zbar(lo(1),lo(2))), state % zbar, state % width)
       call c_f_pointer(c_loc(state_in % dpdA(lo(1),lo(2))), state % dpdA, state % width)
       call c_f_pointer(c_loc(state_in % dpdZ(lo(1),lo(2))), state % dpdZ, state % width)
       call c_f_pointer(c_loc(state_in % dedA(lo(1),lo(2))), state % dedA, state % width)
       call c_f_pointer(c_loc(state_in % dedZ(lo(1),lo(2))), state % dedZ, state % width)

       state % reset = state_in % reset
       
    type is (eos_t_3D)
       
       lo = state_in % lo
       hi = state_in % hi

       state % N = state_in % N
       state % width(1) = state % N
       state % spec_width(1) = state % N
       state % spec_width(2) = nspec
       state % aux_width(1) = state % N
       state % aux_width(2) = naux
       
       call c_f_pointer(c_loc(state_in % rho(lo(1),lo(2),lo(3))), state % rho, state % width)
       call c_f_pointer(c_loc(state_in % T(lo(1),lo(2),lo(3))), state % T, state % width)
       call c_f_pointer(c_loc(state_in % p(lo(1),lo(2),lo(3))), state % p, state % width)
       call c_f_pointer(c_loc(state_in % e(lo(1),lo(2),lo(3))), state % e, state % width)
       call c_f_pointer(c_loc(state_in % h(lo(1),lo(2),lo(3))), state % h, state % width)
       call c_f_pointer(c_loc(state_in % s(lo(1),lo(2),lo(3))), state % s, state % width)
       call c_f_pointer(c_loc(state_in % dpdT(lo(1),lo(2),lo(3))), state % dpdT, state % width)
       call c_f_pointer(c_loc(state_in % dpdr(lo(1),lo(2),lo(3))), state % dpdr, state % width)
       call c_f_pointer(c_loc(state_in % dedT(lo(1),lo(2),lo(3))), state % dedT, state % width)
       call c_f_pointer(c_loc(state_in % dedr(lo(1),lo(2),lo(3))), state % dedr, state % width)
       call c_f_pointer(c_loc(state_in % dhdT(lo(1),lo(2),lo(3))), state % dhdT, state % width)
       call c_f_pointer(c_loc(state_in % dhdr(lo(1),lo(2),lo(3))), state % dhdr, state % width)
       call c_f_pointer(c_loc(state_in % dsdT(lo(1),lo(2),lo(3))), state % dsdT, state % width)
       call c_f_pointer(c_loc(state_in % dsdr(lo(1),lo(2),lo(3))), state % dsdr, state % width)
       call c_f_pointer(c_loc(state_in % dpde(lo(1),lo(2),lo(3))), state % dpde, state % width)
       call c_f_pointer(c_loc(state_in % dpdr_e(lo(1),lo(2),lo(3))), state % dpdr_e, state % width)
       call c_f_pointer(c_loc(state_in % xn(lo(1),lo(2),lo(3),1)), state % xn, state % spec_width)
       if (naux > 0) then
          call c_f_pointer(c_loc(state_in % aux(lo(1),lo(2),lo(3),1)), state % aux, state % aux_width)
       endif
       call c_f_pointer(c_loc(state_in % cv(lo(1),lo(2),lo(3))), state % cv, state % width)
       call c_f_pointer(c_loc(state_in % cp(lo(1),lo(2),lo(3))), state % cp, state % width)
       call c_f_pointer(c_loc(state_in % xne(lo(1),lo(2),lo(3))), state % xne, state % width)
       call c_f_pointer(c_loc(state_in % xnp(lo(1),lo(2),lo(3))), state % xnp, state % width)
       call c_f_pointer(c_loc(state_in % eta(lo(1),lo(2),lo(3))), state % eta, state % width)
       call c_f_pointer(c_loc(state_in % pele(lo(1),lo(2),lo(3))), state % pele, state % width)
       call c_f_pointer(c_loc(state_in % ppos(lo(1),lo(2),lo(3))), state % ppos, state % width)
       call c_f_pointer(c_loc(state_in % mu(lo(1),lo(2),lo(3))), state % mu, state % width)
       call c_f_pointer(c_loc(state_in % mu_e(lo(1),lo(2),lo(3))), state % mu_e, state % width)
       call c_f_pointer(c_loc(state_in % y_e(lo(1),lo(2),lo(3))), state % y_e, state % width)
       call c_f_pointer(c_loc(state_in % dedX(lo(1),lo(2),lo(3),1)), state % dedX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dpdX(lo(1),lo(2),lo(3),1)), state % dpdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % dhdX(lo(1),lo(2),lo(3),1)), state % dhdX, state % spec_width)
       call c_f_pointer(c_loc(state_in % gam1(lo(1),lo(2),lo(3))), state % gam1, state % width)
       call c_f_pointer(c_loc(state_in % cs(lo(1),lo(2),lo(3))), state % cs, state % width)
       call c_f_pointer(c_loc(state_in % abar(lo(1),lo(2),lo(3))), state % abar, state % width)
       call c_f_pointer(c_loc(state_in % zbar(lo(1),lo(2),lo(3))), state % zbar, state % width)
       call c_f_pointer(c_loc(state_in % dpdA(lo(1),lo(2),lo(3))), state % dpdA, state % width)
       call c_f_pointer(c_loc(state_in % dpdZ(lo(1),lo(2),lo(3))), state % dpdZ, state % width)
       call c_f_pointer(c_loc(state_in % dedA(lo(1),lo(2),lo(3))), state % dedA, state % width)
       call c_f_pointer(c_loc(state_in % dedZ(lo(1),lo(2),lo(3))), state % dedZ, state % width)

       state % reset = state_in % reset
       
    end select
       
  end subroutine eos_vector_in

  

  subroutine eos_vector_out(state, state_out)

    use iso_c_binding
    
    implicit none

    class (eos_type),    intent(inout) :: state_out
    type (eos_t_vector), intent(in)    :: state

    ! At present there are no EOS types that need to have any special handling.
    
    select type (state_out)

    end select

    
  end subroutine eos_vector_out



  subroutine eos_copy(state_in, state_out)

    class (eos_type) :: state_in, state_out

    select type (state_in)

    type is (eos_t_1D)
       select type (state_out)
       type is (eos_t_1D)
          state_out % rho = state_in % rho
          state_out % T   = state_in % T
          state_out % p   = state_in % p
          state_out % e   = state_in % e
          state_out % h   = state_in % h
          state_out % s   = state_in % s
          state_out % dpdT = state_in % dpdT
          state_out % dpdr = state_in % dpdr
          state_out % dedT = state_in % dedT
          state_out % dedr = state_in % dedr
          state_out % dhdT = state_in % dhdT
          state_out % dhdr = state_in % dhdr
          state_out % dsdT = state_in % dsdT
          state_out % dsdr = state_in % dsdr
          state_out % dpde = state_in % dpde
          state_out % dpdr_e = state_in % dpdr_e
          state_out % xn = state_in % xn
          state_out % aux = state_in % aux
          state_out % cv = state_in % cv
          state_out % cp = state_in % cp
          state_out % xne = state_in % xne
          state_out % xnp = state_in % xnp
          state_out % eta = state_in % eta
          state_out % pele = state_in % pele
          state_out % ppos = state_in % ppos
          state_out % mu = state_in % mu
          state_out % mu_e = state_in % mu
          state_out % y_e = state_in % y_e
          state_out % dedX = state_in % dedX
          state_out % dpdX = state_in % dpdX
          state_out % dhdX = state_in % dhdX
          state_out % gam1 = state_in % gam1
          state_out % cs = state_in % cs
          state_out % abar = state_in % abar
          state_out % zbar = state_in % zbar
          state_out % dpdA = state_in % dpdA
          state_out % dpdZ = state_in % dpdZ
          state_out % dedA = state_in % dedA
          state_out % dedZ = state_in % dedZ

          state_out % reset = state_in % reset
       class default
          call bl_error("Error: Cannot copy a eos_t_1D to a different EOS type.")
       end select

    type is (eos_t_2D)
       select type (state_out)
       type is (eos_t_2D)
          state_out % rho = state_in % rho
          state_out % T   = state_in % T
          state_out % p   = state_in % p
          state_out % e   = state_in % e
          state_out % h   = state_in % h
          state_out % s   = state_in % s
          state_out % dpdT = state_in % dpdT
          state_out % dpdr = state_in % dpdr
          state_out % dedT = state_in % dedT
          state_out % dedr = state_in % dedr
          state_out % dhdT = state_in % dhdT
          state_out % dhdr = state_in % dhdr
          state_out % dsdT = state_in % dsdT
          state_out % dsdr = state_in % dsdr
          state_out % dpde = state_in % dpde
          state_out % dpdr_e = state_in % dpdr_e
          state_out % xn = state_in % xn
          state_out % aux = state_in % aux
          state_out % cv = state_in % cv
          state_out % cp = state_in % cp
          state_out % xne = state_in % xne
          state_out % xnp = state_in % xnp
          state_out % eta = state_in % eta
          state_out % pele = state_in % pele
          state_out % ppos = state_in % ppos
          state_out % mu = state_in % mu
          state_out % mu_e = state_in % mu
          state_out % y_e = state_in % y_e
          state_out % dedX = state_in % dedX
          state_out % dpdX = state_in % dpdX
          state_out % dhdX = state_in % dhdX
          state_out % gam1 = state_in % gam1
          state_out % cs = state_in % cs
          state_out % abar = state_in % abar
          state_out % zbar = state_in % zbar
          state_out % dpdA = state_in % dpdA
          state_out % dpdZ = state_in % dpdZ
          state_out % dedA = state_in % dedA
          state_out % dedZ = state_in % dedZ

          state_out % reset = state_in % reset
       class default
          call bl_error("Error: Cannot copy a eos_t_2D to a different EOS type.")
       end select

    type is (eos_t_3D)
       select type (state_out)
       type is (eos_t_3D)
          state_out % rho = state_in % rho
          state_out % T   = state_in % T
          state_out % p   = state_in % p
          state_out % e   = state_in % e
          state_out % h   = state_in % h
          state_out % s   = state_in % s
          state_out % dpdT = state_in % dpdT
          state_out % dpdr = state_in % dpdr
          state_out % dedT = state_in % dedT
          state_out % dedr = state_in % dedr
          state_out % dhdT = state_in % dhdT
          state_out % dhdr = state_in % dhdr
          state_out % dsdT = state_in % dsdT
          state_out % dsdr = state_in % dsdr
          state_out % dpde = state_in % dpde
          state_out % dpdr_e = state_in % dpdr_e
          state_out % xn = state_in % xn
          state_out % aux = state_in % aux
          state_out % cv = state_in % cv
          state_out % cp = state_in % cp
          state_out % xne = state_in % xne
          state_out % xnp = state_in % xnp
          state_out % eta = state_in % eta
          state_out % pele = state_in % pele
          state_out % ppos = state_in % ppos
          state_out % mu = state_in % mu
          state_out % mu_e = state_in % mu
          state_out % y_e = state_in % y_e
          state_out % dedX = state_in % dedX
          state_out % dpdX = state_in % dpdX
          state_out % dhdX = state_in % dhdX
          state_out % gam1 = state_in % gam1
          state_out % cs = state_in % cs
          state_out % abar = state_in % abar
          state_out % zbar = state_in % zbar
          state_out % dpdA = state_in % dpdA
          state_out % dpdZ = state_in % dpdZ
          state_out % dedA = state_in % dedA
          state_out % dedZ = state_in % dedZ

          state_out % reset = state_in % reset
       class default
          call bl_error("Error: Cannot copy a eos_t_3D to a different EOS type.")
       end select

    class default
       call bl_error("Error: eos_copy does not recognize the input EOS state type.")

    end select

  end subroutine eos_copy



  subroutine eos_type_error(err, input, pt_index)
    
    implicit none

    integer,           intent(in) :: err
    integer,           intent(in) :: input
    integer, optional, intent(in) :: pt_index(3)

    integer :: dim_ptindex

    character (len=64) :: err_string, zone_string, eos_input_str

    write(eos_input_str, '(A13, I1)') ' EOS input = ', input

    if (err .eq. ierr_input) then

      err_string = 'EOS: invalid input.'

    elseif (err .eq. ierr_init) then

      err_string = 'EOS: input variables were not initialized.'

    elseif (err .eq. ierr_init_xn) then

      err_string = 'EOS: species abundances were not initialized.'

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

       dim_ptindex = 3

       if (pt_index(3) .eq. -99) then
          dim_ptindex = 2
          if (pt_index(2) .eq. -99) then
             dim_ptindex = 1
             if (pt_index(1) .eq. -99) then
                dim_ptindex = 0
             endif
          endif
       endif

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

    err_string = err_string // zone_string
    
    call bl_error(err_string)

  end subroutine eos_type_error



  subroutine eos_deallocate(state)

    implicit none
    
    class (eos_type) :: state

    select type (state)

    type is (eos_t)

       ! Nothing to do here since no arrays are allocated.

    type is (eos_t_vector)

       ! Nothing to do here since we don't allocate new space
       ! for this type, we only point to existing data.

    type is (eos_t_1D)

       call bl_deallocate(state % rho)
       call bl_deallocate(state % T)
       call bl_deallocate(state % p)
       call bl_deallocate(state % e)
       call bl_deallocate(state % h)
       call bl_deallocate(state % s)
       call bl_deallocate(state % dpdT)
       call bl_deallocate(state % dpdr)
       call bl_deallocate(state % dedT)
       call bl_deallocate(state % dedr)
       call bl_deallocate(state % dhdT)
       call bl_deallocate(state % dhdr)
       call bl_deallocate(state % dsdT)
       call bl_deallocate(state % dsdr)
       call bl_deallocate(state % dpde)
       call bl_deallocate(state % dpdr_e)
       call bl_deallocate(state % xn)
       call bl_deallocate(state % aux)
       call bl_deallocate(state % cv)
       call bl_deallocate(state % cp)
       call bl_deallocate(state % xne)
       call bl_deallocate(state % xnp)
       call bl_deallocate(state % eta)
       call bl_deallocate(state % pele)
       call bl_deallocate(state % ppos)
       call bl_deallocate(state % mu)
       call bl_deallocate(state % mu_e)
       call bl_deallocate(state % y_e)
       call bl_deallocate(state % dedX)
       call bl_deallocate(state % dpdX)
       call bl_deallocate(state % dhdX)
       call bl_deallocate(state % gam1)
       call bl_deallocate(state % cs)
       call bl_deallocate(state % abar)
       call bl_deallocate(state % zbar)
       call bl_deallocate(state % dpdA)
       call bl_deallocate(state % dpdZ)
       call bl_deallocate(state % dedA)
       call bl_deallocate(state % dedZ)

    type is (eos_t_2D)

       call bl_deallocate(state % rho)
       call bl_deallocate(state % T)
       call bl_deallocate(state % p)
       call bl_deallocate(state % e)
       call bl_deallocate(state % h)
       call bl_deallocate(state % s)
       call bl_deallocate(state % dpdT)
       call bl_deallocate(state % dpdr)
       call bl_deallocate(state % dedT)
       call bl_deallocate(state % dedr)
       call bl_deallocate(state % dhdT)
       call bl_deallocate(state % dhdr)
       call bl_deallocate(state % dsdT)
       call bl_deallocate(state % dsdr)
       call bl_deallocate(state % dpde)
       call bl_deallocate(state % dpdr_e)
       call bl_deallocate(state % xn)
       call bl_deallocate(state % aux)
       call bl_deallocate(state % cv)
       call bl_deallocate(state % cp)
       call bl_deallocate(state % xne)
       call bl_deallocate(state % xnp)
       call bl_deallocate(state % eta)
       call bl_deallocate(state % pele)
       call bl_deallocate(state % ppos)
       call bl_deallocate(state % mu)
       call bl_deallocate(state % mu_e)
       call bl_deallocate(state % y_e)
       call bl_deallocate(state % dedX)
       call bl_deallocate(state % dpdX)
       call bl_deallocate(state % dhdX)
       call bl_deallocate(state % gam1)
       call bl_deallocate(state % cs)
       call bl_deallocate(state % abar)
       call bl_deallocate(state % zbar)
       call bl_deallocate(state % dpdA)
       call bl_deallocate(state % dpdZ)
       call bl_deallocate(state % dedA)
       call bl_deallocate(state % dedZ)

    type is (eos_t_3D)

       call bl_deallocate(state % rho)
       call bl_deallocate(state % T)
       call bl_deallocate(state % p)
       call bl_deallocate(state % e)
       call bl_deallocate(state % h)
       call bl_deallocate(state % s)
       call bl_deallocate(state % dpdT)
       call bl_deallocate(state % dpdr)
       call bl_deallocate(state % dedT)
       call bl_deallocate(state % dedr)
       call bl_deallocate(state % dhdT)
       call bl_deallocate(state % dhdr)
       call bl_deallocate(state % dsdT)
       call bl_deallocate(state % dsdr)
       call bl_deallocate(state % dpde)
       call bl_deallocate(state % dpdr_e)
       call bl_deallocate(state % xn)
       call bl_deallocate(state % aux)
       call bl_deallocate(state % cv)
       call bl_deallocate(state % cp)
       call bl_deallocate(state % xne)
       call bl_deallocate(state % xnp)
       call bl_deallocate(state % eta)
       call bl_deallocate(state % pele)
       call bl_deallocate(state % ppos)
       call bl_deallocate(state % mu)
       call bl_deallocate(state % mu_e)
       call bl_deallocate(state % y_e)
       call bl_deallocate(state % dedX)
       call bl_deallocate(state % dpdX)
       call bl_deallocate(state % dhdX)
       call bl_deallocate(state % gam1)
       call bl_deallocate(state % cs)
       call bl_deallocate(state % abar)
       call bl_deallocate(state % zbar)
       call bl_deallocate(state % dpdA)
       call bl_deallocate(state % dpdZ)
       call bl_deallocate(state % dedA)
       call bl_deallocate(state % dedZ)

    end select
    
  end subroutine eos_deallocate


  
end module eos_type_module
