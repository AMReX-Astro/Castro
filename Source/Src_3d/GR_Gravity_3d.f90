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
      double precision :: rho, e, G, P, C, T, dpdr, dpde, X(nspec)

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

                  rho =  var(i,j,k,URHO)
                  e   =  var(i,j,k,UEINT) / rho
                  T   =  var(i,j,k,UTEMP)
                  do n = 1, nspec
                     X(n)= var(i,j,k,UFS+n-1)/rho
                  enddo
   
                  ! Compute pressure from the EOS
                  pt_index(1) = i
                  pt_index(2) = j
                  pt_index(3) = k
                  call eos_given_ReX(G, P, C, T, dpdr, dpde, rho, e, X, pt_index=pt_index)

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
                              radial_pres(index) = radial_pres(index) + P
                           end if

                        end do
                     end do
                  end do

               end if
            enddo
         enddo
      enddo

      end subroutine ca_compute_avgpres

! ::
! :: ----------------------------------------------------------
! ::

      subroutine ca_integrate_gr_grav (rho,pres,grav,dr,numpts_1d)

      use bl_constants_module,          only : M_PI
      use fundamental_constants_module, only : Gconst

      implicit none
      integer          :: numpts_1d
      double precision ::  rho(0:numpts_1d-1)
      double precision :: pres(0:numpts_1d-1)
      double precision :: grav(0:numpts_1d-1)
      double precision :: dr

      integer          :: i
      double precision :: mass_encl,rc,rlo,halfdr
      double precision :: ga, gb, gc, P,R

      double precision, parameter ::  fourpi       = 4.d0 * M_PI
      double precision, parameter ::  fourthirdspi = fourpi / 3.d0
      double precision, parameter ::  sqvc         = 29979245800.d0**2

      halfdr    = 0.5d0 * dr
      mass_encl = 0.d0

      do i = 0,numpts_1d-1

         rc  = (dble(i)+0.5d0) * dr
         rlo = (dble(i)      ) * dr
         if (i.eq.0) then
            mass_encl = fourthirdspi * rc**3 * rho(i)
         else
            mass_encl = mass_encl + fourthirdspi * halfdr * (rlo**2 + rlo*(rlo-halfdr) + (rlo-halfdr)**2) * rho(i-1) + &
                                    fourthirdspi * halfdr * ( rc**2 +  rc* rlo         +  rlo**2        ) * rho(i  )
         end if

         grav(i) = -Gconst * mass_encl / rc**2

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This adds the post-Newtonian correction
!!       Added by Ken Chen, 6/9 2010
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!       Tolman-Oppenheimer-Volkoff(TOV) case

         if (rho(i) .gt. 0.d0) then

            P  =  pres(i)
            R  =   rho(i)
            ga = (1.d0 + P/(R*sqvc))
            gb = (1.d0 + fourpi * rc**3 * P / (mass_encl*sqvc))
            gc = 1.d0 / (1.d0 - 2.d0 * Gconst * mass_encl / (rc*sqvc))

            grav(i) = grav(i)*ga*gb*gc

         end if

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This ends the post-Newtonian correction
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

      end subroutine ca_integrate_gr_grav

