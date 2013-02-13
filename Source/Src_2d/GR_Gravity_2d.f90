! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_avgpres (lo,hi,dx,dr,&
                                     var,r_l1,r_l2,r_h1,r_h2, &
                                     radial_pres,radial_vol,problo, &
                                     n1d,drdxfac,level)
      use probdata_module
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS
      use eos_module
      use network, only : nspec

      implicit none

      integer          :: lo(2),hi(2)
      double precision :: dx(2),dr
      double precision :: problo(2)

      integer          :: n1d,drdxfac,level
      double precision :: radial_pres(0:n1d-1)
      double precision :: radial_vol (0:n1d-1)

      integer          :: r_l1,r_l2,r_h1,r_h2
      double precision :: var(r_l1:r_h1,r_l2:r_h2,NVAR)

      integer          :: i,j,n,index
      integer          :: ii,jj
      integer          :: pt_index(2)
      double precision :: xc,yc,r
      double precision :: fac,xx,yy,dx_frac,dy_frac,vol_frac
      double precision :: lo_i,lo_j,rlo,rhi
      double precision :: rho, e, G, P, C, T, dpdr, dpde, X(nspec)

      fac  = dble(drdxfac)
      dx_frac = dx(1) / fac
      dy_frac = dx(2) / fac

      do j = lo(2), hi(2)
         yc = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)

         do i = lo(1), hi(1)
            xc = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)

            r = sqrt(xc**2  + yc**2)
            index = int(r/dr)

            if (index .gt. n1d-1) then

               if (level .eq. 0) then
                  print *,'   '
                  print *,'>>> Error: Gravity_2d::ca_compute_avgpres ',i,j
                  print *,'>>> ... index too big: ', index,' > ',n1d-1
                  print *,'>>> ... at (i,j)     : ',i,j
                  print *,'    '
                  call bl_error("Error:: Gravity_2d.f90 :: ca_compute_avgpres")
               end if

            else

               rho =  var(i,j,URHO)
               e   =  var(i,j,UEINT) / rho
               T   =  var(i,j,UTEMP)
               do n = 1, nspec
                  X(n)= var(i,j,UFS+n-1)/rho
               enddo

               ! Compute pressure from the EOS
               pt_index(1) = i
               pt_index(2) = j
               call eos_given_ReX(G, P, C, T, dpdr, dpde, rho, e, X, pt_index=pt_index)

               ! Note that we assume we are in r-z coordinates in 2d or we wouldn't be 
               !      doing monopole gravity
               lo_i = problo(1) + dble(i)*dx(1) - center(1)
               lo_j = problo(2) + dble(j)*dx(2) - center(2)
               do ii = 0,drdxfac-1
                  xx  = lo_i + (dble(ii  )+0.5d0)*dx_frac
                  rlo = lo_i +  dble(ii  )       *dx_frac
                  rhi = lo_i +  dble(ii+1)       *dx_frac
                  vol_frac = (rhi**2 - rlo**2)  * dy_frac
                  do jj = 0,drdxfac-1
                     yy = lo_j + (dble(jj)+0.5d0)*dy_frac
                      r = sqrt(xx**2  + yy**2)
                     index = int(r/dr)

                     if (index .le. n1d-1) then
                        radial_pres(index) = radial_pres(index) + vol_frac * P
                        radial_vol (index) = radial_vol (index) + vol_frac
                     end if

                  end do
               end do

            end if
         enddo
      enddo

      end subroutine ca_compute_avgpres

