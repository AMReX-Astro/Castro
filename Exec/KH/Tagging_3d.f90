! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the Laplacian.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_laplac_error(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                 set,clear, &
                                 var,varl1,varl2,varl3,varh1,varh2,varh3, &
                                 lo,hi,nd,domlo,domhi, &
                                 delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision :: var(varl1:varh1,varl2:varh2,varl3:varh3)
      double precision :: delta(3), xlo(3), problo(3), time
      integer          :: i,j,k

      double precision ::  delu(3,varl1:varh1,varl2:varh2,varl3:varh3)
      double precision :: delua(3,varl1:varh1,varl2:varh2,varl3:varh3)
      double precision :: delu2(9), delu3(9), delu4(9)
      double precision :: num, denom, error

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8
      double precision, parameter:: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(1,i,j,k) =     var(i+1,j,k) -      var(i-1,j,k) 
         delua(1,i,j,k) = abs(var(i+1,j,k)) + abs(var(i-1,j,k))
      end do
      end do
      end do

      ! d/dy
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(2,i,j,k) =     var(i,j+1,k)  -     var(i,j-1,k) 
         delua(2,i,j,k) = abs(var(i,j+1,k)) + abs(var(i,j-1,k))
      end do
      end do
      end do

      ! d/dz
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(3,i,j,k) =     var(i,j,k+1)  -     var(i,j,k-1)
         delua(3,i,j,k) = abs(var(i,j,k+1)) + abs(var(i,j,k-1))
      end do
      end do
      end do

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         

         ! d/dxdx
          delu2(1) =     delu(1,i+1,j,k)  -     delu(1,i-1,j,k)
          delu3(1) = abs(delu(1,i+1,j,k)) + abs(delu(1,i-1,j,k))
          delu4(1) =    delua(1,i+1,j,k)  +    delua(1,i-1,j,k)

          ! d/dydx
          delu2(2) =     delu(1,i,j+1,k)  -     delu(1,i,j-1,k)
          delu3(2) = abs(delu(1,i,j+1,k)) + abs(delu(1,i,j-1,k))
          delu4(2) =    delua(1,i,j+1,k)  +    delua(1,i,j-1,k)

          ! d/dxdy
          delu2(3) =     delu(2,i+1,j,k)  -     delu(2,i-1,j,k)
          delu3(3) = abs(delu(2,i+1,j,k)) + abs(delu(2,i-1,j,k))
          delu4(3) =    delua(2,i+1,j,k)  +    delua(2,i-1,j,k)

          ! d/dydy
          delu2(4) =     delu(2,i,j+1,k)  -     delu(2,i,j-1,k)
          delu3(4) = abs(delu(2,i,j+1,k)) + abs(delu(2,i,j-1,k))
          delu4(4) =    delua(2,i,j+1,k)  +    delua(2,i,j-1,k)

          ! d/dzdx
          delu2(5) =     delu(1,i,j,k+1)  -     delu(1,i,j,k-1)
          delu3(5) = abs(delu(1,i,j,k+1)) + abs(delu(1,i,j,k-1))
          delu4(5) =    delua(1,i,j,k+1)  +    delua(1,i,j,k-1)

          ! d/dzdy
          delu2(6) =     delu(2,i,j,k+1)  -     delu(2,i,j,k-1)
          delu3(6) = abs(delu(2,i,j,k+1)) + abs(delu(2,i,j,k-1))
          delu4(6) =    delua(2,i,j,k+1)  +    delua(2,i,j,k-1)

          ! d/dxdz
          delu2(7) =     delu(3,i+1,j,k)  -     delu(3,i-1,j,k)
          delu3(7) = abs(delu(3,i+1,j,k)) + abs(delu(3,i-1,j,k))
          delu4(7) =    delua(3,i+1,j,k)  +    delua(3,i-1,j,k)

          ! d/dydz
          delu2(8) =     delu(3,i,j+1,k)  -     delu(3,i,j-1,k)
          delu3(8) = abs(delu(3,i,j+1,k)) + abs(delu(3,i,j-1,k))
          delu4(8) =    delua(3,i,j+1,k)  +    delua(3,i,j-1,k)

          ! d/dzdz
          delu2(9) =     delu(3,i,j,k+1)  -     delu(3,i,j,k-1)
          delu3(9) = abs(delu(3,i,j,k+1)) + abs(delu(3,i,j,k-1))
          delu4(9) =    delua(3,i,j,k+1)  +    delua(3,i,j,k-1)

         ! compute the error
         num   =  delu2(1)**2 + delu2(2)**2 + delu2(3)**2 + delu2(4)**2 &
                 +delu2(5)**2 + delu2(6)**2 + delu2(7)**2 + delu2(8)**2 &
                 +delu2(9)**2

         denom = (delu3(1) + (epsil*delu4(1)+1.d-99))**2 + &
                 (delu3(2) + (epsil*delu4(2)+1.d-99))**2 + &
                 (delu3(3) + (epsil*delu4(3)+1.d-99))**2 + &
                 (delu3(4) + (epsil*delu4(4)+1.d-99))**2 + &
                 (delu3(5) + (epsil*delu4(5)+1.d-99))**2 + &
                 (delu3(6) + (epsil*delu4(6)+1.d-99))**2 + &
                 (delu3(7) + (epsil*delu4(7)+1.d-99))**2 + &
                 (delu3(8) + (epsil*delu4(8)+1.d-99))**2 + &
                 (delu3(9) + (epsil*delu4(9)+1.d-99))**2

         error = sqrt(num/denom)

         if (error .gt. ctore) tag(i,j,k)=set

      end do
      end do
      end do

      end subroutine ca_laplac_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_denerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             den,denl1,denl2,denl3,denh1,denh2,denh3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision den(denl1:denh1,denl2:denh2,denl3:denh3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
                  ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
                  az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
                  ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
                  ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
                  az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. dengrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the velocity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vel       => velocity array
! ::: nv        => number of components in vel array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_velerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             vel,vell1,vell2,vell3,velh1,velh2,velh3, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, nv, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer vell1,vell2,vell3,velh1,velh2,velh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision vel(vell1:velh1,vell2:velh2,vell3:velh3,nv)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

!     Tag on regions of high velocity gradient
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(vel(i+1,j,k,1) - vel(i,j,k,1))
                  ay = ABS(vel(i,j+1,k,1) - vel(i,j,k,1))
                  az = ABS(vel(i,j,k+1,1) - vel(i,j,k,1))
                  ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1,j,k,1)))
                  ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1,k,1)))
                  az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. velgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo

      end

