 
! ::
! :: ----------------------------------------------------------
! :: summass
! ::             MASS = sum{ vol(i,j,k)*rho(i,j,k) }
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! :: ----------------------------------------------------------
! ::

      subroutine ca_summass(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,mass)

        use bl_constants_module

        implicit none

        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: lo(3), hi(3)
        double precision :: mass, dx(3)
        double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        
        integer          :: i, j, k
        double precision :: vol

        vol  = dx(1)*dx(2)*dx(3)
        mass = ZERO

        do k = lo(3), hi(3)
           do i = lo(1), hi(1)
              do j = lo(2), hi(2)
                 mass = mass + rho(i,j,k)
              enddo
           enddo
        enddo

        mass = mass * vol

      end subroutine ca_summass

! ::
! :: ----------------------------------------------------------
! :: sumsquared
! ::             MASS = sum{ vol(i,j,k)*rho(i,j,k)**2) }
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  sum of squared quantity
! :: ----------------------------------------------------------
! ::

      subroutine ca_sumsquared(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi,dx,mass)

        use bl_constants_module

        implicit none

        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: lo(3), hi(3)
        double precision :: mass, dx(3)
        double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        
        integer          :: i, j, k
        double precision :: vol

        vol  = dx(1)*dx(2)*dx(3)
        mass = ZERO

        do k = lo(3), hi(3)
           do i = lo(1), hi(1)
              do j = lo(2), hi(2)
                 mass = mass + rho(i,j,k)*rho(i,j,k)
              enddo
           enddo
        enddo

        mass = mass * vol

      end subroutine ca_sumsquared

! ::
! :: ----------------------------------------------------------
! :: sumlocmass
! ::             MASS = sum{ vol(i,j,k)*rho(i,j,k)*x_idir }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  r          => radius at cell center
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumlocmass(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi, &
                                problo,dx,mass,idir)

       use prob_params_module, only : center
       use bl_constants_module

       implicit none
       integer          :: idir
       integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
       integer          :: lo(3), hi(3)
       double precision :: mass, problo(3), dx(3)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

       integer          :: i, j, k
       double precision :: x,y,z,vol

       vol  = dx(1)*dx(2)*dx(3)
       mass = ZERO

       if (idir .eq. 0) then
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   mass = mass + rho(i,j,k) * x
                enddo
             enddo
          enddo
       else if (idir .eq. 1) then
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * y
                enddo
             enddo
          enddo
       else
          do k = lo(3), hi(3)
             z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * z
                enddo
             enddo
          enddo
       end if

       mass = mass * vol

       end subroutine ca_sumlocmass




! ::
! :: ----------------------------------------------------------
! :: sumlocmass2d
! ::             MASS = sum{ vol(i,j,k)*rho(i,j,k)*x_idir1*x_idir2 }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  problo     => physical location of the origin
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  idir1      => first direction to weight by
! ::  idir2      => second direction to weight by
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumlocmass2d(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi, &
                                  problo,dx,mass,idir1,idir2)

       use prob_params_module, only : center
       use bl_constants_module

       implicit none
       integer          :: idir1, idir2
       integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
       integer          :: lo(3), hi(3)
       double precision :: mass, problo(3), dx(3)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

       integer          :: i, j, k
       double precision :: x,y,z,vol

       vol  = dx(1)*dx(2)*dx(3)
       mass = ZERO

       if (idir1 .eq. 0) then
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             if (idir2 .eq. 0) then
                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      mass = mass + rho(i,j,k) * x * x
                   enddo
                enddo
             elseif (idir2 .eq. 1) then
                do j = lo(2), hi(2)
                   y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
                   do k = lo(3), hi(3)
                      mass = mass + rho(i,j,k) * x * y
                   enddo
                enddo
             else
                do k = lo(3), hi(3)
                   z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
                   do j = lo(2), hi(2)
                      mass = mass + rho(i,j,k) * x * z
                   enddo
                enddo
             endif
          enddo
       else if (idir1 .eq. 1) then
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             if (idir2 .eq. 0) then
                do i = lo(1), hi(1)
                   x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
                   do k = lo(3), hi(3)
                      mass = mass + rho(i,j,k) * y * x
                   enddo
                enddo
             elseif (idir2 .eq. 1) then
                do i = lo(1), hi(1)
                   do k = lo(3), hi(3)
                      mass = mass + rho(i,j,k) * y * y
                   enddo
                enddo
             else
                do k = lo(3), hi(3)
                   z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
                   do i = lo(1), hi(1)
                      mass = mass + rho(i,j,k) * y * z
                   enddo
                enddo
             endif
          enddo
       else
          do k = lo(3), hi(3)
             z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
             if (idir2 .eq. 0) then
                do i = lo(1), hi(1)
                   x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
                   do j = lo(2), hi(2)
                      mass = mass + rho(i,j,k) * z * x
                   enddo
                enddo
             elseif (idir2 .eq. 1) then
                do j = lo(2), hi(2)
                   y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
                   do i = lo(1), hi(1)
                      mass = mass + rho(i,j,k) * z * y
                   enddo
                enddo
             else
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      mass = mass + rho(i,j,k) * z * z
                   enddo
                enddo
             endif
          enddo
       end if

       mass = mass * vol

       end subroutine ca_sumlocmass2d



! ::
! :: ----------------------------------------------------------
! :: sumlocsquaredmass
! ::             MASS = sum{ vol(i,j,k)*rho(i,j,k)*x_idir }
! ::
! :: INPUTS / OUTPUTS:
! ::  rho        => density field
! ::  rlo,rhi    => index limits of rho array
! ::  lo,hi      => index limits of grid interior
! ::  dx         => cell size
! ::  mass      <=  total mass
! ::  r          => radius at cell center
! ::  irlo,hi    => index limits of r array
! ::  rz_flag    => == 1 if R_Z coords
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumlocsquaredmass(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi, &
                                       problo,dx,mass,idir)

       use prob_params_module, only : center
       use bl_constants_module

       implicit none
       integer          :: idir
       integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
       integer          :: lo(3), hi(3)
       double precision :: mass, problo(3), dx(3)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

       integer          :: i, j, k
       double precision :: x,y,z,vol

       vol  = dx(1)*dx(2)*dx(3)
       mass = ZERO

       if (idir .eq. 0) then
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+HALF) * dx(1) - center(1)
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   mass = mass + rho(i,j,k) * (x**2)
                enddo
             enddo
          enddo
       else if (idir .eq. 1) then
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+HALF) * dx(2) - center(2)
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * (y**2)
                enddo
             enddo
          enddo
       else
          do k = lo(3), hi(3)
             z = problo(3) + (dble(k)+HALF) * dx(3) - center(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * (z**2)
                enddo
             enddo
          enddo
       end if

       mass = mass * vol

       end subroutine ca_sumlocsquaredmass

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

       subroutine ca_sumproduct(f1, f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3, &
                                f2, f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3, &
                                lo,hi,dx,product)

       use bl_constants_module

       implicit none

       integer          :: f1_l1,f1_l2,f1_l3,f1_h1,f1_h2,f1_h3
       integer          :: f2_l1,f2_l2,f2_l3,f2_h1,f2_h2,f2_h3
       integer          :: lo(3), hi(3)
       double precision :: product
       double precision :: dx(3), vol
       double precision :: f1(f1_l1:f1_h1,f1_l2:f1_h2,f1_l3:f1_h3)
       double precision :: f2(f2_l1:f2_h1,f2_l2:f2_h2,f2_l3:f2_h3)

       integer          :: i,j,k

       product = ZERO
       vol = dx(1) * dx(2) * dx(3)
 
       do i = lo(1), hi(1)
         do k = lo(3), hi(3)
           do j = lo(2), hi(2)
             product = product + f1(i,j,k)*f2(i,j,k)
           enddo
         enddo
       enddo

       product = product * vol

       end subroutine ca_sumproduct
