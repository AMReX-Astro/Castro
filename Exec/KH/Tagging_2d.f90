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
      subroutine ca_laplac_error(tag,tagl1,tagl2,tagh1,tagh2, &
                                 set,clear, &
                                 var,varl1,varl2,varh1,varh2, &
                                 lo,hi,nd,domlo,domhi, &
                                 delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagh1,tagh2
      integer          :: varl1,varl2,varh1,varh2
      integer          :: lo(2), hi(2), domlo(2), domhi(2)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2)
      double precision :: var(varl1:varh1,varl2:varh2)
      double precision :: delta(2), xlo(2), problo(2), time
      integer          :: i,j

      double precision ::  delu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)
      double precision :: delua(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)
      double precision :: delu2(4), delu3(4), delu4(4)
      double precision :: num, denom, error

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8
      double precision, parameter :: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,1) =     var(i+1,j) -      var(i-1,j)
         delua(i,j,1) = abs(var(i+1,j)) + abs(var(i-1,j))
      end do
      end do

      ! d/dy
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,2) =     var(i,j+1) -      var(i,j-1)
         delua(i,j,2) = abs(var(i,j+1)) + abs(var(i,j-1))
      end do
      end do

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         ! d/dxdx
         delu2(1) =     delu(i+1,j,1)  -     delu(i-1,j,1)
         delu3(1) = abs(delu(i+1,j,1)) + abs(delu(i-1,j,1))
         delu4(1) =    delua(i+1,j,1)  +    delua(i-1,j,1)
                                                         
         ! d/dydx                                        
         delu2(2) =     delu(i,j+1,1)  -     delu(i,j-1,1)
         delu3(2) = abs(delu(i,j+1,1)) + abs(delu(i,j-1,1))
         delu4(2) =    delua(i,j+1,1)  +    delua(i,j-1,1)
                                                         
         ! d/dxdy                                        
         delu2(3) =     delu(i+1,j,2)  -     delu(i-1,j,2)
         delu3(3) = abs(delu(i+1,j,2)) + abs(delu(i-1,j,2))
         delu4(3) =    delua(i+1,j,2)  +    delua(i-1,j,2)
                                                         
         ! d/dydy                                        
         delu2(4) =     delu(i,j+1,2)  -     delu(i,j-1,2)
         delu3(4) = abs(delu(i,j+1,2)) + abs(delu(i,j-1,2))
         delu4(4) =    delua(i,j+1,2)  +    delua(i,j-1,2)

         ! compute the error
         num   =  delu2(1)**2 + delu2(2)**2 + delu2(3)**2 + delu2(4)**2
         denom = (delu3(1) + (epsil*delu4(1)+1.d-99))**2 + &
                 (delu3(2) + (epsil*delu4(2)+1.d-99))**2 + &
                 (delu3(3) + (epsil*delu4(3)+1.d-99))**2 + &
                 (delu3(4) + (epsil*delu4(4)+1.d-99))**2

         error = sqrt(num/denom)

         if (error .gt. ctore) tag(i,j)=set

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
      subroutine ca_denerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             den,denl1,denl2,denh1,denh2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer denl1,denl2,denh1,denh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision den(denl1:denh1,denl2:denh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high density gradient
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ax = ABS(den(i+1,j,1) - den(i,j,1))
            ay = ABS(den(i,j+1,1) - den(i,j,1))
            ax = MAX(ax,ABS(den(i,j,1) - den(i-1,j,1)))
            ay = MAX(ay,ABS(den(i,j,1) - den(i,j-1,1)))
            if ( MAX(ax,ay) .ge. dengrad) then
               tag(i,j) = set
            endif
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
! ::: np        => number of components in velocity array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine ca_velerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             vel,vell1,vell2,velh1,velh2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer vell1,vell2,velh1,velh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision vel(vell1:velh1,vell2:velh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high velocity gradient
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            ax = ABS(vel(i+1,j,1) - vel(i,j,1))
            ay = ABS(vel(i,j+1,1) - vel(i,j,1))
            ax = MAX(ax,ABS(vel(i,j,1) - vel(i-1,j,1)))
            ay = MAX(ay,ABS(vel(i,j,1) - vel(i,j-1,1)))
            if ( MAX(ax,ay) .ge. velgrad) then
               tag(i,j) = set
            endif
         enddo
      enddo
      
      end

