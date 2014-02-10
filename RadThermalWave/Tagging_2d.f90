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

      double precision ::  delu(2,varl1:varh1,varl2:varh2)
      double precision :: delua(2,varl1:varh1,varl2:varh2)
      double precision :: delu2(4), delu3(4), delu4(4)
      double precision :: num, denom, error

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8
      double precision, parameter :: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(1,i,j) =     var(i+1,j) -      var(i-1,j)
         delua(1,i,j) = abs(var(i+1,j)) + abs(var(i-1,j))
      end do
      end do

      ! d/dy
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(2,i,j) =     var(i,j+1) -      var(i,j-1)
         delua(2,i,j) = abs(var(i,j+1)) + abs(var(i,j-1))
      end do
      end do

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         ! d/dxdx
         delu2(1) =     delu(1,i+1,j)  -     delu(1,i-1,j)
         delu3(1) = abs(delu(1,i+1,j)) + abs(delu(1,i-1,j))
         delu4(1) =    delua(1,i+1,j)  +    delua(1,i-1,j)

         ! d/dydx
         delu2(2) =     delu(1,i,j+1)  -     delu(1,i,j-1)
         delu3(2) = abs(delu(1,i,j+1)) + abs(delu(1,i,j-1))
         delu4(2) =    delua(1,i,j+1)  +    delua(1,i,j-1)

         ! d/dxdy
         delu2(3) =     delu(2,i+1,j)  -     delu(2,i-1,j)
         delu3(3) = abs(delu(2,i+1,j)) + abs(delu(2,i-1,j))
         delu4(3) =    delua(2,i+1,j)  +    delua(2,i-1,j)

         ! d/dydy
         delu2(4) =     delu(2,i,j+1)  -     delu(2,i,j-1)
         delu3(4) = abs(delu(2,i,j+1)) + abs(delu(2,i,j-1))
         delu4(4) =    delua(2,i,j+1)  +    delua(2,i,j-1)

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

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (den(i,j,1) .ge. denerr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
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
      endif
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the temperature
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => temperature array
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine ca_temperror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             temp,templ1,templ2,temph1,temph2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer templ1,templ2,temph1,temph2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision temp(templ1:temph1,templ2:temph2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay, gradT
      integer i, j

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (temp(i,j,1) .ge. temperr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(temp(i+1,j,1) - temp(i,j,1))
               ay = ABS(temp(i,j+1,1) - temp(i,j,1))
               ax = MAX(ax,ABS(temp(i,j,1) - temp(i-1,j,1)))
               ay = MAX(ay,ABS(temp(i,j,1) - temp(i,j-1,1)))
               gradT = sqrt(ax**2 + ay**2)
               if (gradT / temp(i,j,1) .ge. tempgrad) then
              
!               if ( MAX(ax*delta(1),ay*delta(2))/temp(i,j,1) &
!                    .ge. tempgrad * 2.d0**(level-1) ) then
!                    .ge. tempgrad ) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the pressure
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: press     => pressure array
! ::: np        => number of components in press array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_presserror(tag,tagl1,tagl2,tagh1,tagh2, &
                               set,clear, &
                               press,pressl1,pressl2, &
                                     pressh1,pressh2, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagh1,tagh2
      integer pressl1,pressl2,pressh1,pressh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision press(pressl1:pressh1,pressl2:pressh2,np)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (press(i,j,1) .ge. presserr) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(press(i+1,j,1) - press(i,j,1))
                  ay = ABS(press(i,j+1,1) - press(i,j,1))
                  ax = MAX(ax,ABS(press(i,j,1) - press(i-1,j,1)))
                  ay = MAX(ay,ABS(press(i,j,1) - press(i,j-1,1)))
                  if ( MAX(ax,ay) .ge. pressgrad) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

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
      subroutine ca_velerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             vel,vell1,vell2,velh1,velh2, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nv, level
      integer tagl1,tagl2,tagh1,tagh2
      integer vell1,vell2,velh1,velh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision vel(vell1:velh1,vell2:velh2,nv)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
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
      endif

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the radiation
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: rad       => radiation array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_raderror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             rad,radl1,radl2, &
                                 radh1,radh2, &
                             lo,hi,nr,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, nr, level
      integer tagl1,tagl2,tagh1,tagh2
      integer radl1,radl2,radh1,radh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision rad(radl1:radh1,radl2:radh2,nr)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high radiation
      if (level .lt. max_raderr_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (rad(i,j,1) .ge. raderr) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

!     Tag on regions of high radiation gradient
      if (level .lt. max_radgrad_lev) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(rad(i+1,j,1) - rad(i,j,1))
                  ay = ABS(rad(i,j+1,1) - rad(i,j,1))
                  ax = MAX(ax,ABS(rad(i,j,1) - rad(i-1,j,1)))
                  ay = MAX(ay,ABS(rad(i,j,1) - rad(i,j-1,1)))
                  if ( MAX(ax,ay) .ge. radgrad) then
                     tag(i,j) = set
                  endif
               enddo
            enddo
      endif

      end
