! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_compute_avgpres (lo,hi,dx,dr,&
                                     var,r_l1,r_l2,r_h1,r_h2, &
                                     radial_pres,problo, &
                                     n1d,drdxfac,level) bind(C,name='ca_compute_avgpres')
      use prob_params_module, only : center
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, UFX
      use eos_module
      use network, only : nspec, naux
      use amrex_constants_module
      use castro_error_module

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer , intent(in   ) :: lo(2),hi(2)
      real(rt), intent(in   ) :: dx(2),dr
      real(rt), intent(in   ) :: problo(2)

      integer , intent(in   ) :: n1d,drdxfac,level
      real(rt), intent(inout) :: radial_pres(0:n1d-1)

      integer , intent(in   ) :: r_l1,r_l2,r_h1,r_h2
      real(rt), intent(in   ) :: var(r_l1:r_h1,r_l2:r_h2,NVAR)

      integer          :: i,j,n,index
      integer          :: ii,jj
      real(rt)         :: xc,yc,r
      real(rt)         :: fac,xx,yy,dx_frac,dy_frac,vol_frac
      real(rt)         :: lo_i,lo_j,rlo,rhi

      type (eos_t) :: eos_state

      fac  = dble(drdxfac)
      dx_frac = dx(1) / fac
      dy_frac = dx(2) / fac

      do j = lo(2), hi(2)
         yc = problo(2) + (dble(j)+HALF) * dx(2) - center(2)

         do i = lo(1), hi(1)
            xc = problo(1) + (dble(i)+HALF) * dx(1) - center(1)

            r = sqrt(xc**2  + yc**2)
            index = int(r/dr)

            if (index .gt. n1d-1) then

               if (level .eq. 0) then
                  print *,'   '
                  print *,'>>> Error: Gravity_2d::ca_compute_avgpres ',i,j
                  print *,'>>> ... index too big: ', index,' > ',n1d-1
                  print *,'>>> ... at (i,j)     : ',i,j
                  print *,'    '
                  call castro_error("Error:: Gravity_2d.f90 :: ca_compute_avgpres")
               end if

            else

               eos_state % rho = var(i,j,URHO)
               eos_state % e   = var(i,j,UEINT) / eos_state % rho
               eos_state % T   = var(i,j,UTEMP)
               eos_state % xn  = var(i,j,UFS:UFS+nspec-1) / eos_state % rho 
               eos_state % aux = var(i,j,UFX:UFX+naux-1) / eos_state % rho

               ! Compute pressure from the EOS
               call eos(eos_input_re, eos_state)

               ! Note that we assume we are in r-z coordinates in 2d or we wouldn't be 
               !      doing monopole gravity
               lo_i = problo(1) + dble(i)*dx(1) - center(1)
               lo_j = problo(2) + dble(j)*dx(2) - center(2)
               do ii = 0,drdxfac-1
                  xx  = lo_i + (dble(ii  )+HALF)*dx_frac
                  rlo = lo_i +  dble(ii  )      *dx_frac
                  rhi = lo_i +  dble(ii+1)      *dx_frac
                  vol_frac = (rhi**2 - rlo**2)  *dy_frac
                  do jj = 0,drdxfac-1
                     yy = lo_j + (dble(jj)+HALF)*dy_frac
                     r = sqrt(xx**2  + yy**2)
                     index = int(r/dr)

                     if (index .le. n1d-1) then
                        radial_pres(index) = radial_pres(index) + vol_frac * eos_state % P
                     end if

                  end do
               end do

            end if
         enddo
      enddo

      end subroutine ca_compute_avgpres

