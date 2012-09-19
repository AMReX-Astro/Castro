
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
        implicit none
        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: lo(3), hi(3)
        double precision :: mass, dx(3)
        double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        
        integer          :: i, j, k
        double precision :: vol

        vol  = dx(1)*dx(2)*dx(3)
        mass = 0.d0

        !$OMP PARALLEL DO PRIVATE(i,j,k) reduction(+:mass)
        do k = lo(3), hi(3)
           do i = lo(1), hi(1)
              do j = lo(2), hi(2)
                 mass = mass + rho(i,j,k)
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

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
        implicit none
        integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
        integer          :: lo(3), hi(3)
        double precision :: mass, dx(3)
        double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)
        
        integer          :: i, j, k
        double precision :: vol

        vol  = dx(1)*dx(2)*dx(3)
        mass = 0.d0

        !$OMP PARALLEL DO PRIVATE(i,j,k) reduction(+:mass)
        do k = lo(3), hi(3)
           do i = lo(1), hi(1)
              do j = lo(2), hi(2)
                 mass = mass + rho(i,j,k)*rho(i,j,k)
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

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
! ::  irlo,hi    => index limits of r array
! ::  rz_flag    => == 1 if R_Z coords
! :: ----------------------------------------------------------
! ::

       subroutine ca_sumlocmass(rho,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,lo,hi, &
                                problo,dx,mass,idir)

       use probdata_module, only : center

       implicit none
       integer          :: idir
       integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
       integer          :: lo(3), hi(3)
       double precision :: mass, problo(3), dx(3)
       double precision :: rho(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

       integer          :: i, j, k
       double precision :: x,y,z,vol

       vol  = dx(1)*dx(2)*dx(3)
       mass = 0.d0

       if (idir .eq. 0) then
          !$OMP PARALLEL DO PRIVATE(i,j,k,x) REDUCTION(+:mass)
          do i = lo(1), hi(1)
             x = problo(1) + (dble(i)+0.5d0) * dx(1) - center(1)
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   mass = mass + rho(i,j,k) * x
                enddo
             enddo
          enddo
          !$OMP END PARALLEL DO
       else if (idir .eq. 1) then
          !$OMP PARALLEL DO PRIVATE(i,j,k,y) REDUCTION(+:mass)
          do j = lo(2), hi(2)
             y = problo(2) + (dble(j)+0.5d0) * dx(2) - center(2)
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * y
                enddo
             enddo
          enddo
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(i,j,k,z) REDUCTION(+:mass)
          do k = lo(3), hi(3)
             z = problo(3) + (dble(k)+0.5d0) * dx(3) - center(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   mass = mass + rho(i,j,k) * z
                enddo
             enddo
          enddo
          !$OMP END PARALLEL DO
       end if

       mass = mass * vol

       end subroutine ca_sumlocmass
