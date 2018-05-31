subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use eos_module
  use eos_type_module
  use bl_error_module
  use network
  use probdata_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(2), probhi(2)
  real(rt)         xn(nspec)

  integer untin,i

  type (eos_t) :: eos_state

  namelist /fortin/ p_l, u_l, v_l, rho_l, rhoe_l, p_r, u_r, v_r, rho_r, rhoe_r, T_l, T_r, frac, idir, &
       use_Tinit

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer maxlen
  parameter (maxlen=256)
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call bl_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do
         
  ! set namelist defaults

  p_l = 116.5             ! left pressure (erg/cc)
  u_l = 7.1447096          ! left u (cm/s)
  v_l = -4.125          ! left v (cm/s)
  rho_l = 8.0             ! left density (g/cc)
  T_l = 1.0

  p_r = 1.0               ! right pressure (erg/cc)
  u_r = 0.0               ! right u (cm/s)
  v_r = 0.0               ! right v (cm/s)
  rho_r = 1.4             ! right density (g/cc)
  T_r = 1.0

  use_Tinit = .false.     ! optionally use T_l/r instead of p_l/r for initialization

  !     Read namelists
  untin = 9
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  !     set local variable defaults -- the 'center' variables are the location of the

  !     compute the internal energy (erg/cc) for the left and right state
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
     eos_state%T = 10.e0_rt  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)
 
     rhoe_l = rho_l*eos_state%e
     T_l = eos_state%T

     eos_state%rho = rho_r
     eos_state%p = p_r
     eos_state%T = 10.e0_rt  ! initial guess
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)
 
     rhoe_r = rho_r*eos_state%e
     T_r = eos_state%T

  endif

end subroutine amrex_probinit


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nstate    => number of state components.  You should know
! :::		   this already!
! ::: state     <=  Scalar array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS

  use amrex_constants_module, only : M_PI, sixth
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer level, nscal
  integer lo(2), hi(2)
  integer state_l1,state_l2,state_h1,state_h2
  real(rt)         xlo(2), xhi(2), time, delta(2)
  real(rt)         state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt), parameter :: gp(2) = (/ -1.d0/sqrt(3.d0), 1.d0/sqrt(3.d0) /)
  real(rt), parameter :: pi_over_3 = M_pi / 3.d0
  real(rt), parameter :: ff = 0.25d0

  real(rt) :: x,y,xcen,ycen,shockfront
  integer  :: i,j,ii,jj

  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(float(j-lo(2))+ 0.5e0_rt)
     do i = lo(1), hi(1)
        xcen = xlo(1) + delta(1)*(float(i-lo(1))+ 0.5e0_rt)

        state(i,j,URHO   ) = 0.d0 
        state(i,j,UMX:UMY) = 0.d0 
        state(i,j,UEDEN  ) = 0.d0 
        state(i,j,UEINT  ) = 0.d0 

        do jj = 1, 2
          y = ycen + 0.5d0*delta(2)*gp(jj)

          do ii = 1, 2
            x = xcen + 0.5d0*delta(1)*gp(ii)

           shockfront = tan(pi_over_3)*(x - sixth) ! initial shock front

           if (y .ge. shockfront ) then
              state(i,j,URHO) = state(i,j,URHO) + ff*rho_l
              state(i,j,UMX ) = state(i,j,UMX ) + ff*rho_l*u_l
              state(i,j,UMY ) = state(i,j,UMY ) + ff*rho_l*v_l

              state(i,j,UEDEN) = state(i,j,UEDEN) + ff*(rhoe_l + 0.5*rho_l*(u_l*u_l + v_l*v_l))
              state(i,j,UEINT) = state(i,j,UEINT) + ff*rhoe_l 
              state(i,j,UTEMP) = state(i,j,UTEMP) + ff*T_l
           else
              state(i,j,URHO) = state(i,j,URHO) + ff*rho_r
              state(i,j,UMX ) = state(i,j,UMX ) + ff*rho_r*u_r
              state(i,j,UMY ) = state(i,j,UMY ) + ff*rho_r*v_r

              state(i,j,UEDEN) = state(i,j,UEDEN) + ff*(rhoe_r + 0.5*rho_r*(u_r*u_r + v_r*v_r))
              state(i,j,UEINT) = state(i,j,UEINT) + ff*rhoe_r
              state(i,j,UTEMP) = state(i,j,UTEMP) + ff*T_r
           endif

          end do
        end do

        state(i,j,UFS:UFS-1+nspec) = 0.0e0_rt
        state(i,j,UFS  ) = state(i,j,URHO)
        
     enddo
  enddo
  
end subroutine ca_initdata

