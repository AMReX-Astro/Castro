
! :: ----------------------------------------------------------
! :: summass
! ::             MASS = sum{ vol(i)*rho(i) }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho aray
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  r          => radius at cell center
! ::  irlo,hi    => index limits of r array
! ::  tmp        => temp column array
! :: ----------------------------------------------------------
! ::
      subroutine ca_summass(rho,r_l1,r_h1,lo,hi, &
                            dx,mass,r,irlo,irhi)
 
      use prob_params_module, only: coord_type 

      implicit none
      integer irlo, irhi
      integer r_l1,r_h1
      integer lo(1), hi(1)
      double precision mass, dx(1)
      double precision rho(r_l1:r_h1)
      double precision r(irlo:irhi)
 
      integer i
      double precision vol,ro,ri
      double precision, parameter :: spherical_factor =4.d0/3.d0*3.14159265358979323846d0

      mass = 0.d0
      if (coord_type .eq. 2) then
        do i = lo(1), hi(1)
           ro = r(i)+0.5d0*dx(1)
           ri = r(i)-0.5d0*dx(1)
           vol = spherical_factor * (ro**3 - ri**3)
           mass = mass + vol*rho(i)
        enddo
      else
        do i = lo(1), hi(1)
           mass = mass + rho(i)
        enddo
        mass = mass * dx(1)
      end if

      end subroutine ca_summass

! :: ----------------------------------------------------------
! :: sumsquared
! ::             MASS = sum{ vol(i)*rho(i)^2 }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho aray
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  r          => radius at cell center
! ::  irlo,hi    => index limits of r array
! ::  tmp        => temp column array
! :: ----------------------------------------------------------
! ::
      subroutine ca_sumsquared(rho,r_l1,r_h1,lo,hi, &
                               dx,mass,r,irlo,irhi)
 
      use prob_params_module, only: coord_type 

      implicit none
      integer irlo, irhi
      integer r_l1,r_h1
      integer lo(1), hi(1)
      double precision mass, dx(1)
      double precision rho(r_l1:r_h1)
      double precision r(irlo:irhi)
 
      integer i
      double precision vol,ro,ri
      double precision, parameter :: spherical_factor =4.d0/3.d0*3.14159265358979323846d0

      mass = 0.d0
      if (coord_type .eq. 2) then
        do i = lo(1), hi(1)
           ro = r(i)+0.5d0*dx(1)
           ri = r(i)-0.5d0*dx(1)
           vol = spherical_factor * (ro**3 - ri**3)
           mass = mass + vol*rho(i)*rho(i)
        enddo
      else
        do i = lo(1), hi(1)
           mass = mass + rho(i)*rho(i)
        enddo
        mass = mass * dx(1)
      end if

      end subroutine ca_sumsquared

! :: ----------------------------------------------------------
! :: sumlocmass
! ::             MASS = sum{ vol(i)*rho(i)*loc(i) }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho aray
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  r          => radius at cell center
! ::  irlo,hi    => index limits of r array
! ::  idir       => direction of location weighting
! ::  tmp        => temp column array
! :: ----------------------------------------------------------
! ::
      subroutine ca_sumlocmass(rho,r_l1,r_h1,lo,hi, &
                               dx,mass,r,irlo,irhi,idir) 
 
      use prob_params_module, only: coord_type 

      implicit none
      integer irlo, irhi, idir
      integer r_l1,r_h1
      integer lo(1), hi(1)
      double precision mass, dx(1)
      double precision rho(r_l1:r_h1)
      double precision r(irlo:irhi)
 
      integer i
      double precision vol,ro,ri
      double precision, parameter :: spherical_factor =4.d0/3.d0*3.14159265358979323846d0

      mass = 0.d0
      if (coord_type .eq. 2) then
        do i = lo(1), hi(1)
           ro = r(i)+0.5d0*dx(1)
           ri = r(i)-0.5d0*dx(1)
           vol = spherical_factor * (ro**3 - ri**3)
           mass = mass + vol*rho(i)*r(i)
        enddo
      else
        do i = lo(1), hi(1)
           mass = mass + rho(i)*r(i)
        enddo
        mass = mass * dx(1)
      end if

      end subroutine ca_sumlocmass

! ::
! :: ----------------------------------------------------------
! :: sumproduct
! ::
! :: INPUTS / OUTPUTS:
! ::  f1,f2      => fields whose product is summed
! ::  flo,fhi    => index limits of field arrays
! ::  lo,hi      => index limits of grid interior
! ::  problo     => lower bound of domain in physical units
! ::  dx         => cell size
! ::  product    <= sum of product of two quantities
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumproduct(f1, f1_l1,f1_h1,&
                                f2, f2_l1,f2_h1,&
                                lo,hi,dx,product)
       implicit none

       integer          :: f1_l1,f1_h1
       integer          :: f2_l1,f2_h1
       integer          :: lo(1), hi(1)
       double precision :: product
       double precision :: dx(1), vol
       double precision :: f1(f1_l1:f1_h1)
       double precision :: f2(f2_l1:f2_h1)

       integer          :: i

       product = 0.d0
       vol = dx(1)
 
       !$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:product)
       do i = lo(1), hi(1)
          product = product + f1(i)*f2(i)
       enddo
       !$OMP END PARALLEL DO

       product = product * vol

       end subroutine ca_sumproduct
