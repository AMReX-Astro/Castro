subroutine amrex_probinit (init, name, namlen, problo, probhi) bind(c)

  use eos_module, only: eos
  use eos_type_module, only: eos_input_rt, eos_input_rp, eos_t
  use castro_error_module, only: castro_error
  use network, only: nspec
  use probdata_module, only: p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, T_l, T_r, &
                             frac, idir, use_Tinit, xcloud, cldradius, split
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: init, namlen
  integer,  intent(in   ) :: name(namlen)
  real(rt), intent(in   ) :: problo(3), probhi(3)

  real(rt) :: xn(nspec)
  integer  :: untin, i

  type (eos_t) :: eos_state

  namelist /fortin/ p_l, u_l, rho_l, p_r, u_r, rho_r, T_l, T_r, frac, idir, &
                    use_Tinit, xcloud, cldradius

  integer, parameter :: maxlen = 256
  character :: probin*(maxlen)

  if (namlen .gt. maxlen) then
     call castro_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults

  p_l = 1.0               ! left pressure (erg/cc)
  u_l = 0.0               ! left velocity (cm/s)
  rho_l = 1.0             ! left density (g/cc)
  T_l = 1.0

  p_r = 0.1               ! right pressure (erg/cc)
  u_r = 0.0               ! right velocity (cm/s)
  rho_r = 0.125           ! right density (g/cc)
  T_r = 1.0

  idir = 1                ! direction across which to jump
  frac = 0.5              ! fraction of the domain for the interface

  use_Tinit = .false.     ! optionally use T_l/r instead of p_l/r for initialization

  xcloud = -1.0
  cldradius = -1.0

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  if (idir /= 1) then
     call castro_error('invalid idir')
  end if

  split(:) = frac * (problo(:) + probhi(:))

  ! compute the internal energy (erg/cc) for the left and right state
  xn(:) = 0.0e0_rt
  xn(1) = 1.0e0_rt

  if (use_Tinit) then

     eos_state%rho = rho_l
     eos_state%T = T_l
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)
 
     rhoe_l = rho_l*eos_state%e
     p_l = eos_state%p

     eos_state%rho = rho_r
     eos_state%T = T_r
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)
 
     rhoe_r = rho_r*eos_state%e
     p_r = eos_state%p

  else

     eos_state%rho = rho_l
     eos_state%p = p_l
     eos_state%T = 100000.e0_rt  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)
 
     rhoe_l = rho_l*eos_state%e
     T_l = eos_state%T

     eos_state%rho = rho_r
     eos_state%p = p_r
     eos_state%T = 100000.e0_rt  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)
 
     rhoe_r = rho_r*eos_state%e
     T_r = eos_state%T

  endif

end subroutine amrex_probinit

subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec
  use probdata_module, only: p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, T_l, T_r, &
                             frac, idir, use_Tinit, xcloud, cldradius, split
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xcen, ycen, zcen
  integer  :: i, j, k

  real(rt) :: rsq, dist, ycloud, zcloud, ecent, denfact, rhocld

  ycloud = 0.0
  zcloud = 0.0
  ecent  = 1.00
  denfact = 3.0246
  rhocld  = rho_r * denfact

  do k = lo(3), hi(3)
     zcen = xlo(3) + dx(3) * (dble(k-lo(3)) + 0.5e0_rt)

     do j = lo(2), hi(2)
        ycen = xlo(2) + dx(2) * (dble(j-lo(2)) + 0.5e0_rt)

        do i = lo(1), hi(1)
           xcen = xlo(1) + dx(1) * (dble(i-lo(1)) + 0.5e0_rt)

           if (xcen <= split(1)) then      !! left of shock

              state(i,j,k,URHO) = rho_l
              state(i,j,k,UMX) = rho_l * u_l
              state(i,j,k,UMY) = 0.e0_rt
              state(i,j,k,UMZ) = 0.e0_rt
              state(i,j,k,UEDEN) = rhoe_l + 0.5 * rho_l * u_l * u_l
              state(i,j,k,UEINT) = rhoe_l
              state(i,j,k,UTEMP) = T_l

           else

              rsq = (((xcen-xcloud)**2) + ((ycen-ycloud)**2) + ((zcen-zcloud)**2)) / ecent

              dist = sqrt(rsq)
              if (dist .le. cldradius) then   !! inside bubble
                 state(i,j,k,URHO) = rhocld
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_r
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              else                            !! right of shock
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = rho_r * u_r
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_r + 0.5 * rho_r * u_r * u_r
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              end if

           end if
 
           state(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
           state(i,j,k,UFS) = state(i,j,k,URHO)
           
        end do
     end do
  end do

end subroutine ca_initdata
