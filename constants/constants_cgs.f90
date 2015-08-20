! Fundamental constants taken from NIST's 2010 CODATA recommended values

module fundamental_constants_module

  use bl_types
  use bl_constants_module, only: M_PI
  
  implicit none

  ! newton's gravitational constant
  real(kind=dp_t), parameter :: Gconst = 6.67428e-8_dp_t      ! cm^3/g/s^2
! new value; if uncommented initial models will need to be re-HSE'ed
!  real(kind=dp_t), parameter :: Gconst = 6.67384e-8_dp_t      ! cm^3/g/s^2

  ! boltzmann's constant
  real(kind=dp_t), parameter :: k_B    = 1.3806488e-16_dp_t   ! erg/K

  ! planck's constant over 2pi
  real(kind=dp_t), parameter :: hbar   = 1.054571726e-27_dp_t ! erg s

  ! planck's constant 
  real(kind=dp_t), parameter :: hplanck = 6.62606957e-27_dp_t ! erg s

  ! avogradro's Number
  real(kind=dp_t), parameter :: n_A    = 6.02214129e23_dp_t   ! mol^-1

  ! convert eV to erg
  real(kind=dp_t), parameter :: ev2erg = 1.602176487e-12_dp_t

  ! convert MeV to eV
  real(kind=dp_t), parameter :: MeV2eV = 1.0e6_dp_t

  ! mass of proton
  real(kind=dp_t), parameter :: m_p     = 1.672621777e-24_dp_t ! g

  ! mass of neutron
  real(kind=dp_t), parameter :: m_n      = 1.674927351e-24_dp_t ! g

  ! mass of electron
  real(kind=dp_t), parameter :: m_e     = 9.10938291e-28_dp_t  ! g

  ! speed of light in vacuum
  real(kind=dp_t), parameter :: c_light = 2.99792458e10_dp_t   ! cm/s

  ! electron charge
  ! NIST: q_e = 1.602176565e-19 C
  !
  ! C is the SI unit Coulomb; in cgs we have the definition:
  !     1 C = 0.1 * |c_light| * 1 statC
  ! where statC is the cgs unit statCoulomb; 1 statC = 1 erg^1/2 cm^1/2
  ! and |c_light| is the speed of light in cgs (but without units)
  real(kind=dp_t), parameter :: q_e     = 4.80320451e-10_dp_t  ! erg^1/2 cm^1/2

  ! stefan-boltzmann constant
  real(kind=dp_t), parameter :: sigma_SB = 5.670373e-5_dp_t    ! erg/s/cm^2/K^4

  ! radiation constant
  real(kind=dp_t), parameter :: a_rad = 4.0_dp_t*sigma_SB/c_light

  ! Number of centimeters in a parsec and an AU.
  ! Note that since the length of an AU is defined exactly
  ! by IAU convention, the length of a parsec is also
  ! defined exactly as (6.48e5 / pi) AU.
  real(kind=dp_t), parameter :: AU = 1.49597871e13_dp_t            ! cm
  real(kind=dp_t), parameter :: parsec = 3.085677587679311e18_dp_t ! cm

  ! Hubble constant (in s^{-1}, converted from 100 (km/s)/Mpc by dividing by 3.08568025e19km/Mpc)
  real(kind=dp_t), parameter :: Hubble_const = 32.407764868e-19_dp_t

  ! solar mass (from http://asa.usno.navy.mil/SecK/Constants.html)
  real(kind=dp_t), parameter :: M_solar = 1.9884e33_dp_t

end module fundamental_constants_module
