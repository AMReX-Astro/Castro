
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
       use bl_constants_module
  
       implicit none
       integer irlo, irhi
       integer r_l1,r_l2,r_h1,r_h2
       integer lo(2), hi(2)
       double precision mass, dx(2)
       double precision rho(r_l1:r_h1,r_l2:r_h2)
       double precision r(irlo:irhi)
       double precision  tmp(lo(2):hi(2))

       integer i, j
       double precision vol

       if (coord_type .eq. 0) then
          vol = dx(1) * dx(2)
       else
          vol = TWO * M_PI * dx(1) * dx(2)
       endif

       do j = lo(2),hi(2)
          tmp(j) = ZERO
       enddo

       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             tmp(j) = tmp(j) + vol * r(i) * rho(i,j)
          enddo
       enddo

       mass = ZERO
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
       use bl_constants_module
  
       implicit none
       integer irlo, irhi
       integer r_l1,r_l2,r_h1,r_h2
       integer lo(2), hi(2)
       double precision mass, dx(2)
       double precision rho(r_l1:r_h1,r_l2:r_h2)
       double precision r(irlo:irhi)
       double precision  tmp(lo(2):hi(2))

       integer i, j
       double precision vol

       if (coord_type .eq. 0) then
          vol = dx(1) * dx(2)
       else
          vol = TWO * M_PI * dx(1) * dx(2)
       endif

       do j = lo(2),hi(2)
          tmp(j) = ZERO
       enddo

       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             tmp(j) = tmp(j) + vol * r(i) * rho(i,j) * rho(i,j)
          enddo
       enddo

       mass = ZERO
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
       use bl_constants_module

       implicit none

       integer          :: irlo, irhi, idir
       integer          :: r_l1,r_l2,r_h1,r_h2
       integer          :: lo(2), hi(2)
       double precision :: mass, problo(2), dx(2)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2)
       double precision :: r(irlo:irhi)

       integer i, j
       double precision x,y,vol_weighting

       if (coord_type .eq. 0) then
          vol_weighting = dx(1) * dx(2)
       else
          vol_weighting = TWO * M_PI * dx(1) * dx(2)
       endif

       mass = ZERO
       if (idir .eq. 0) then
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             do j = lo(2), hi(2)
                mass = mass + (vol_weighting * r(i)) * rho(i,j) * x
             enddo
          enddo
       else if (idir .eq. 1) then
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             do i = lo(1), hi(1)
                mass = mass + (vol_weighting * r(i)) * rho(i,j) * y
             enddo
          enddo
       end if

       end subroutine ca_sumlocmass

! :: ----------------------------------------------------------
! :: sumlocmass2d
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
! ::  idir1,2    => directions of location weighting
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumlocmass2d(rho,r_l1,r_l2,r_h1,r_h2,lo,hi, &
                                  problo,dx,mass,r,irlo,irhi,idir1,idir2)

       use prob_params_module, only : coord_type
       use probdata_module, only : center
       use bl_constants_module

       implicit none

       integer          :: irlo, irhi, idir1, idir2
       integer          :: r_l1,r_l2,r_h1,r_h2
       integer          :: lo(2), hi(2)
       double precision :: mass, problo(2), dx(2)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2)
       double precision :: r(irlo:irhi)

       integer i, j
       double precision x,y,vol_weighting

       if (coord_type .eq. 0) then
          vol_weighting = dx(1) * dx(2)
       else
          vol_weighting = TWO * M_PI * dx(1) * dx(2)
       endif

       mass = ZERO
       if (idir1 .eq. 0) then
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             if (idir2 .eq. 0) then
               do j = lo(2), hi(2)
                  mass = mass + (vol_weighting * r(i)) * rho(i,j) * x * x
               enddo
             else
               do j = lo(2), hi(2)
                  y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
                  mass = mass + (vol_weighting * r(i)) * rho(i,j) * x * y
               enddo
             endif
          enddo
       else 
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             if (idir2 .eq. 0 ) then
               do i = lo(1), hi(1)
                  x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
                  mass = mass + (vol_weighting * r(i)) * rho(i,j) * y * x
               enddo
             else
               do i = lo(1), hi(1)
                  mass = mass + (vol_weighting * r(i)) * rho(i,j) * y * y
               enddo
             endif
          enddo
       end if

       end subroutine ca_sumlocmass2d

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

       use prob_params_module, only : coord_type
       use probdata_module, only : center
       use bl_constants_module

       implicit none

       integer          :: f1_l1,f1_l2,f1_h1,f1_h2
       integer          :: f2_l1,f2_l2,f2_h1,f2_h2
       integer          :: lo(2), hi(2)
       double precision :: product
       double precision :: dx(2), vol
       double precision :: f1(f1_l1:f1_h1,f1_l2:f1_h2)
       double precision :: f2(f2_l1:f2_h1,f2_l2:f2_h2)

       integer          :: i,j

       product = ZERO

       if (coord_type .eq. 0) then
          vol = dx(1) * dx(2)
       else
          vol = TWO * M_PI * dx(1) * dx(2)
       endif
 
       do i = lo(1), hi(1)
           do j = lo(2), hi(2)
             product = product + f1(i,j)*f2(i,j)
           enddo
       enddo

       product = product * vol

       end subroutine ca_sumproduct
