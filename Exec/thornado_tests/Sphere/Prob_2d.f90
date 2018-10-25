subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use eos_module
  use eos_type_module
  use amrex_constants_module, only: half
  use amrex_error_module 
  use network
  use probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: init, namlen
  integer, intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(2), probhi(2)

  integer untin,i

  type (eos_t) :: eos_state

  namelist /fortin/ centx, centy

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 256
  character probin*(maxlen)

  if (namlen .gt. maxlen) then
     call amrex_error("probin file name too long")
  end if

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  !     Read namelists
  open(newunit=untin, file=probin(1:namlen), form='formatted', status='old')
  read(untin,fortin)
  close(unit=untin)

  xn(:) = 0.0e0_rt
  xn(1) = 1.0e0_rt

  centx = half * (problo(1)+probhi(1))
  centy = half * (problo(2)+probhi(2))

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

  use amrex_constants_module, only: zero, half, one
  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX
  
  use UnitsModule
  use EquationOfStateModule_TABLE, only: ComputeThermodynamicStates_Primitive_TABLE

  implicit none

  integer, intent(in) :: level, nscal
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: state_l1,state_l2,state_h1,state_h2
  real(rt), intent(in) :: xlo(2), xhi(2), time, delta(2)
  real(rt), intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  real(rt), allocatable :: rho_in(:), T_in(:), Ye_in(:), Epervol_out(:), Epermass_out(:), Ne_out(:)

  real(rt) :: x,y,radius
  real(rt) :: rho_min,rho_max,r_rho,H_rho,T_min,T_max,r_T,H_T,Ye_min,Ye_max,r_Ye,H_Ye
  real(rt) :: tanh_r, tanh_t, tanh_y
  integer  :: i,j

  allocate(rho_in(lo(1):hi(1)))
  allocate(  T_in(lo(1):hi(1)))
  allocate( Ye_in(lo(1):hi(1)))
  allocate(Epervol_out(lo(1):hi(1)))
  allocate(Epermass_out(lo(1):hi(1)))
  allocate(Ne_out(lo(1):hi(1)))

  if (UFX .lt. 0.d0) &
     call amrex_abort("Must have UFX defined to run this problem!")

  ! ************************ Min and max values of rho, T, Ye ************************
  rho_min = 1.0e8 
  rho_max = 4.0e14

  T_min = 5.0e9 
  T_max = 2.6e11

  Ye_min = 0.3
  Ye_max = 0.46

  ! ************************ Radii and widths for rho, T, Ye ************************
  r_rho = 2.0e6
  H_rho = 1.0e6

  r_T = 2.5e6
  H_T = 2.0e6

  r_Ye = 4.5e6
  H_Ye = 1.0e6

  ! ************************ ************************ ************************

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        x = xlo(1) + delta(1)*(dble(i-lo(1))+half) - centx
        y = xlo(2) + delta(2)*(dble(j-lo(2))+half) - centy

        radius = sqrt(x*x+y*y)
 
        ! These go from near-zero at radius = 0 to near-one at large radius
        tanh_r = half * (one + tanh((radius - r_rho)/H_rho))
        tanh_t = half * (one + tanh((radius - r_T  )/H_T))
        tanh_y = half * (one + tanh((radius - r_Ye )/H_Ye))

        ! The profile has max values at radius = 0 and min values at large radius
        state(i,j,URHO ) = ( rho_max - (rho_max - rho_min) * tanh_r )
        state(i,j,UTEMP) = (   T_max - (  T_max -   T_min) * tanh_t )

        state(i,j,UMX:UMY) = zero

        state(i,j,UFS:UFS-1+nspec) = 0.0e0_rt
        state(i,j,UFS            ) = state(i,j,URHO)

        rho_in(i) = state(i,j,URHO) * (Gram/Centimeter**3)
          T_in(i) = state(i,j,UTEMP) * Kelvin

         Ye_in(i) = (  Ye_min - ( Ye_min -  Ye_max) * tanh_y )
           
     enddo

     call ComputeThermodynamicStates_Primitive_TABLE(rho_in, T_in, Ye_in, Epervol_out, Epermass_out, Ne_out )

     do i = lo(1), hi(1)

        state(i,j,UEINT) =  Epervol_out(i) / (Erg/Centimeter**3)    ! UEINT = (rho e) 
        state(i,j,UEDEN) =  state(i,j,UEINT)   ! (rho E) = (rho e) since momentum = 0
        state(i,j,UFX  ) =  Ne_out(i) * Centimeter**3
        state(i,j,UFS  ) =  Ne_out(i) * Centimeter**3 * AtomicMassUnit / Gram

     enddo

  enddo

  deallocate(rho_in,T_in,Ye_in)
  deallocate(Epervol_out,Epermass_out,Ne_out)

end subroutine ca_initdata

! hardwired assuming 4 moments
subroutine get_rad_ncomp(rad_ncomp) bind(C,name="ca_get_rad_ncomp")

  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule, only : nE, nDOF, nNodesX, nNodesE

  integer :: rad_ncomp
  integer :: n_moments = 4

  rad_ncomp =  nSpecies * n_moments * nE * nNodesX(1) *  nNodesX(2) * nNodesE

end subroutine get_rad_ncomp

! hardwired assuming 4 moments
! streaming sine wave, J = H_x = 1 + sin(2*pi*x)
subroutine ca_init_thornado_data(level,time,lo,hi, &
                                 nrad_comp,rad_state, &
                                 rad_state_l1,rad_state_l2, &
                                 rad_state_h1,rad_state_h2, &
                                 state,state_l1,state_l2,state_h1,state_h2, &
                                 delta,xlo,xhi) bind(C,name="ca_init_thornado_data")

  use probdata_module
  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule, only : nE, nDOF, nNodesX, nNodesE
  use amrex_fort_module, only : rt => amrex_real
  use amrex_error_module
  use amrex_constants_module, only : M_PI
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX
  use MeshModule, only: MeshE, MeshX, NodeCoordinate
  use UnitsModule
  use EquationOfStateModule_TABLE, only: ComputeThermodynamicStates_Auxiliary_TABLE, ComputeElectronChemicalPotential_TABLE, &
                                         ComputeProtonChemicalPotential_TABLE, ComputeNeutronChemicalPotential_TABLE 

  implicit none

  integer , intent(in) :: level, nrad_comp
  integer , intent(in) :: lo(2), hi(2)
  integer , intent(in) :: rad_state_l1,rad_state_h1
  integer , intent(in) :: rad_state_l2,rad_state_h2
  integer, intent(in) :: state_l1,state_l2,state_h1,state_h2
  real(rt), intent(in) :: xlo(2), xhi(2), time, delta(2)
  real(rt), intent(inout) ::  rad_state(rad_state_l1:rad_state_h1,rad_state_l2:rad_state_h2,&
                                        0:nrad_comp-1)
  real(rt), intent(inout) :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  ! Local parameter
  integer, parameter :: n_moments = 4

  ! local variables
  integer :: i,j,ixnode,iynode,ienode
  integer :: ii,ii_0,is,im,ie,id
  integer :: nx,ny
  real(rt) :: xcen, ycen, xnode, ynode
  real(rt) :: rho_in(1), T_in(1), Ye_in(1), Evol(1), Ne_loc(1), Em_in(1), M_e(1), M_p(1), M_n(1), M_nu(1), E(1)
  ! zero it out, just in case
  rad_state = 0.0e0_rt

  print *,'nrad_comp ',nrad_comp
  print *,'nSpecies  ',nSpecies
  print *,'n_moments ',n_moments
  print *,'nE        ',nE
  print *,'nNodesX   ',nNodesX(:)
  print *,'nNodesE   ',nNodesE
  print *,'MULT ', nSpecies * n_moments * nE * nNodesX(1) *  nNodesX(2) * nNodesE

  print *,'nDOF      ',nDOF

  ny = nNodesE*nNodesX(1)
  nx = nNodesE

  if (nDOF .ne. ny*nNodesX(2)) then
     print *,'nDOF is ', nDOF
     print *,'nNodesX(1)*nNodesX(2)*nNodesE is ', nNodesX(1)*nNodesX(2)*nNodesE
     call amrex_abort("nDOF ne nNodesX(1)*nNodesX(2)*nNodesE")
  end if
     
  do j = lo(2), hi(2)
     ycen = xlo(2) + delta(2)*(float(j-lo(2)) + 0.5e0_rt) - centy

     do i = lo(1), hi(1)

        xcen = xlo(1) + delta(1)*(float(i-lo(1)) + 0.5e0_rt) - centx
        
        ! Get Castro fluid variables unit convert to thornado units
        rho_in(1) = state(i,j,URHO) * Gram / Centimeter**3
        T_in(1) = state(i,j,UTEMP) * Kelvin
        Evol(1) = state(i,j,UEINT) * (Erg/Centimeter**3)
        Ne_loc(1) = state(i,j,UFX) / Centimeter**3  
        
        ! Calculate chemical potentials via thornado subroutines
        call ComputeThermodynamicStates_Auxiliary_TABLE( rho_in, Evol, Ne_loc, T_in, Em_in, Ye_in) 
        call ComputeElectronChemicalPotential_TABLE(rho_in,T_in,Ye_in,M_e)        
        call ComputeProtonChemicalPotential_TABLE(rho_in,T_in,Ye_in,M_p)        
        call ComputeNeutronChemicalPotential_TABLE(rho_in,T_in,Ye_in,M_n)        

        M_nu = M_e + M_p - M_n

        do is = 1, nSpecies
        do im = 1, n_moments
        do ie = 1, nE

           ii_0 = (is-1)*(n_moments*nE*nDOF) + (im-1)*(nE*nDOF) + (ie-1)*nDOF

           do iynode = 1, nNodesX(2)
           do ixnode = 1, nNodesX(1)
           do ienode = 1, nNodesE
 
              ! Calculate the indices
              id = (ienode-1) + nx*(ixnode-1) + ny*(iynode-1)
              ii = ii_0 + id

              ! Calculate actual positions of the nodes used for the gaussian quadrature
              xnode = xcen + ( float(ixnode)-1.5e0_rt )*delta(1)/sqrt(3.0e0_rt)
              ynode = ycen + ( float(iynode)-1.5e0_rt )*delta(2)/sqrt(3.0e0_rt)

              ! Get energy at given node coordinate via thornado subroutine
              E = NodeCoordinate( MeshE, ie, ienode)

              ! J moment, im = 1
              if (im .eq. 1) then 
                      rad_state(i,j,ii) = max(1.0e0_rt / (exp( (E(1)-M_nu(1)) / T_in(1))  + 1.0e0_rt), 1.0d-99)
              endif 
              ! H_x moment, im = 2
              if (im .eq. 2) rad_state(i,j,ii) = 0.0e0_rt  

              ! H_y moment, im = 3
              if (im .eq. 3) rad_state(i,j,ii) = 0.0e0_rt
   
              ! H_z moment, im = 4
              if (im .eq. 4) rad_state(i,j,ii) = 0.0e0_rt
   
           end do
           end do
           end do

        end do
        end do
        end do

     enddo
  enddo

end subroutine ca_init_thornado_data


