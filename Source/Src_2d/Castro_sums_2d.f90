
! :: ----------------------------------------------------------
! :: summass
! ::             MASS = sum{ vol(i,j)*rho(i,j) }
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

       subroutine ca_summass(rho,r_l1,r_l2,r_h1,r_h2,lo,hi, &
                             dx,mass,r,irlo,irhi)

       use prob_params_module, only : coord_type
  
       implicit none
       integer irlo, irhi
       integer r_l1,r_l2,r_h1,r_h2
       integer lo(2), hi(2)
       double precision mass, dx(2)
       double precision rho(r_l1:r_h1,r_l2:r_h2)
       double precision r(irlo:irhi)
       double precision  tmp(lo(2):hi(2))

       double precision, parameter ::  Pi = 3.1415926535897932384d0

       integer i, j
       double precision vol

       if (coord_type .eq. 0) then
          vol = dx(1) * dx(2)
       else
          vol = 2.d0 * Pi * dx(1) * dx(2)
       endif

       do j = lo(2),hi(2)
          tmp(j) = 0.d0
       enddo

       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             tmp(j) = tmp(j) + vol * r(i) * rho(i,j)
          enddo
       enddo

       mass = 0.d0
       do j = lo(2), hi(2)
          mass = mass + tmp(j)
       enddo

       end subroutine ca_summass

! :: ----------------------------------------------------------
! :: sumsquared
! ::             MASS = sum{ vol(i,j)*rho(i,j)^2 }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho aray
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  sum of quantity squared
! ::  r          => radius at cell center
! ::  irlo,hi    => index limits of r array
! ::  tmp        => temp column array
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumsquared(rho,r_l1,r_l2,r_h1,r_h2,lo,hi, &
                                dx,mass,r,irlo,irhi)

       use prob_params_module, only : coord_type
  
       implicit none
       integer irlo, irhi
       integer r_l1,r_l2,r_h1,r_h2
       integer lo(2), hi(2)
       double precision mass, dx(2)
       double precision rho(r_l1:r_h1,r_l2:r_h2)
       double precision r(irlo:irhi)
       double precision  tmp(lo(2):hi(2))

       double precision, parameter ::  Pi = 3.1415926535897932384d0

       integer i, j
       double precision vol

       if (coord_type .eq. 0) then
          vol = dx(1) * dx(2)
       else
          vol = 2.d0 * Pi * dx(1) * dx(2)
       endif

       do j = lo(2),hi(2)
          tmp(j) = 0.d0
       enddo

       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             tmp(j) = tmp(j) + vol * r(i) * rho(i,j) * rho(i,j)
          enddo
       enddo

       mass = 0.d0
       do j = lo(2), hi(2)
          mass = mass + tmp(j)
       enddo

       end subroutine ca_sumsquared

! :: ----------------------------------------------------------
! :: sumlocmass
! ::             MASS = sum{ vol(i,j)*rho(i,j) }
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
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumlocmass(rho,r_l1,r_l2,r_h1,r_h2,lo,hi, &
                                problo,dx,mass,r,irlo,irhi,idir)

       use prob_params_module, only : coord_type
       use probdata_module, only : center

       implicit none

       integer          :: irlo, irhi, idir
       integer          :: r_l1,r_l2,r_h1,r_h2
       integer          :: lo(2), hi(2)
       double precision :: mass, problo(2), dx(2)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2)
       double precision :: r(irlo:irhi)

       double precision, parameter ::  Pi = 3.1415926535897932384d0

       integer i, j
       double precision x,y,vol_weighting

       if (coord_type .eq. 0) then
          vol_weighting = dx(1) * dx(2)
       else
          vol_weighting = 2.d0 * Pi * dx(1) * dx(2)
       endif

       mass = 0.d0
       if (idir .eq. 0) then
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+0.5d0) * dx(1) - center(1)
             do j = lo(2), hi(2)
                mass = mass + (vol_weighting * r(i)) * rho(i,j) * x
             enddo
          enddo
       else if (idir .eq. 1) then
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+0.5d0) * dx(2) - center(2)
             do i = lo(1), hi(1)
                mass = mass + (vol_weighting * r(i)) * rho(i,j) * y
             enddo
          enddo
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

       subroutine ca_sumproduct(f1, f1_l1,f1_l2,f1_h1,f1_h2,&
                                f2, f2_l1,f2_l2,f2_h1,f2_h2,&
                                lo,hi,dx,product)
       implicit none

       integer          :: f1_l1,f1_l2,f1_h1,f1_h2
       integer          :: f2_l1,f2_l2,f2_h1,f2_h2
       integer          :: lo(2), hi(2)
       double precision :: product
       double precision :: dx(2), vol
       double precision :: f1(f1_l1:f1_h1,f1_l2:f1_h2)
       double precision :: f2(f2_l1:f2_h1,f2_l2:f2_h2)

       integer          :: i,j

       product = 0.d0
       vol = dx(1) * dx(2)
 
       !$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:product)
       do i = lo(1), hi(1)
           do j = lo(2), hi(2)
             product = product + f1(i,j)*f2(i,j)
           enddo
       enddo
       !$OMP END PARALLEL DO

       product = product * vol

       end subroutine ca_sumproduct
