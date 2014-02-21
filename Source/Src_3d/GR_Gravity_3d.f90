! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_avgpres (lo,hi,dx,dr,&
                                    var,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3,&
                                    radial_pres,problo,&
                                    n1d,drdxfac,level)
      use probdata_module
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS
      use eos_module
      use network, only : nspec

      implicit none

      integer          :: lo(3),hi(3)
      double precision :: dx(3),dr
      double precision :: problo(3)

      integer          :: n1d,drdxfac,level
      double precision :: radial_pres(0:n1d-1)

      integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
      double precision :: var(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3,NVAR)

      integer          :: i,j,k,n,index
      integer          :: ii,jj,kk
      integer          :: pt_index(3)
      double precision :: xc,yc,zc,r
      double precision :: fac,xx,yy,zz,dx_frac,dy_frac,dz_frac
      double precision :: lo_i,lo_j,lo_k

      type (eos_t) :: eos_state

      fac     = dble(drdxfac)
      dx_frac = dx(1) / fac
      dy_frac = dx(2) / fac
      dz_frac = dx(3) / fac
      !
      ! Don't OMP this.
      !
      do k = lo(3), hi(3)
         zc = problo(3) + (dble(k)+0.50d0) * dx(3) - center(3)

         do j = lo(2), hi(2)
            yc = problo(2) + (dble(j)+0.50d0) * dx(2) - center(2)

            do i = lo(1), hi(1)
               xc    = problo(1) + (dble(i)+0.50d0) * dx(1) - center(1)
               r     = sqrt(xc**2 + yc**2 + zc**2)
               index = int(r/dr)

               if (index .gt. n1d-1) then

                  if (level .eq. 0) then
                     print *,'   '  
                     print *,'>>> Error: Gravity_3d::ca_compute_avgpres ',i,j,k
                     print *,'>>> ... index too big: ', index,' > ',n1d-1
                     print *,'>>> ... at (i,j,k)   : ',i,j,k
                     call bl_error("Error:: Gravity_3d.f90 :: ca_compute_avgpres")
                  end if

               else

                  eos_state % rho = var(i,j,k,URHO)
                  eos_state % e   = var(i,j,k,UEINT) / rho
                  eos_state % T   = var(i,j,k,UTEMP)
                  eos_state % xn  = var(i,j,k,UFS:UFS+nspec-1) / rho
   
                  ! Compute pressure from the EOS
                  pt_index(1) = i
                  pt_index(2) = j
                  pt_index(3) = k
                  call eos(eos_input_re, eos_state, pt_index = pt_index)

                  lo_i =  problo(1) + dble(i)*dx(1) - center(1)
                  lo_j =  problo(2) + dble(j)*dx(2) - center(2)
                  lo_k =  problo(3) + dble(k)*dx(3) - center(3)

                  do kk = 0,drdxfac-1
                     zz = lo_k + (dble(kk)+0.5d0)*dz_frac
                     do jj = 0,drdxfac-1
                        yy = lo_j + (dble(jj)+0.5d0)*dy_frac
                        do ii = 0,drdxfac-1

                           xx    = lo_i + (dble(ii)+0.5d0)*dx_frac
                           r     = sqrt(xx**2  + yy**2 + zz**2)
                           index = int(r/dr)

                           if (index .le. n1d-1) then
                              radial_pres(index) = radial_pres(index) + eos_state % P
                           end if

                        end do
                     end do
                  end do

               end if
            enddo
         enddo
      enddo

      end subroutine ca_compute_avgpres
