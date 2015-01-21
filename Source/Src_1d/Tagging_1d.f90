
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the Laplacian.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_laplac_error(tag,tagl1,tagh1, &
                                 set,clear, &
                                 var,varl1,varh1, &
                                 lo,hi,nd,domlo,domhi, &
                                 delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagh1
      integer          :: varl1,varh1
      integer          :: lo(1), hi(1), domlo(1), domhi(1)
      integer          :: tag(tagl1:tagh1)
      double precision :: var(varl1:varh1)
      double precision :: delta(1), xlo(1), problo(1), time
      integer          :: i

      double precision :: error, num, denom
      double precision ::  delu(lo(1)-1:hi(1)+1)
      double precision :: delua(lo(1)-1:hi(1)+1)
      double precision :: delu2, delu3, delu4

      double precision, parameter :: ctore=0.8
      double precision, parameter :: eps=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do i=lo(1)-1,hi(1)+1
         delu(i)  =     var(i+1)  -     var(i-1) 
         delua(i) = abs(var(i+1)) + abs(var(i-1))
      end do

      ! d/dxdx
      do i = lo(1),hi(1)

         delu2 =     delu(i+1)  -     delu(i-1)
         delu3 = abs(delu(i+1)) + abs(delu(i-1))
         delu4 =    delua(i+1)  +    delua(i-1)

         ! compute the error
         num   = abs(delu2)
         denom = abs(delu3 + (eps*delu4+1.d-99))

         error = num/denom

         if (error .gt. ctore) tag(i)=set

      end do
      
      end subroutine ca_laplac_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_denerror(tag,tagl1,tagh1, &
                             set,clear, &
                             den,denl1,denh1, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagh1
      integer denl1,denh1
      integer lo(1), hi(1), domlo(1), domhi(1)
      integer tag(tagl1:tagh1)
      double precision den(denl1:denh1,nd)
      double precision delta(1), xlo(1), problo(1), time

      double precision ax
      integer i

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do i = lo(1), hi(1)
            if (den(i,1) .ge. denerr) then
               tag(i) = set
            endif
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         do i = lo(1), hi(1)
            ax = ABS(den(i+1,1) - den(i,1))
            ax = MAX(ax,ABS(den(i,1) - den(i-1,1)))
            if ( ax .ge. dengrad) then
               tag(i) = set
            endif
         enddo
      endif
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the pressure
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => temperature array
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine ca_temperror(tag,tagl1,tagh1,&
                              set,clear, &
                              temp,templ1,temph1,&
                              lo,hi,nd,domlo,domhi, &
                              delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagh1
      integer templ1,temph1
      integer lo(1), hi(1), domlo(1), domhi(1)
      integer tag(tagl1:tagh1)
      double precision temp(templ1:temph1,nd)
      double precision delta(1), xlo(1), problo(1), time

      double precision ax
      integer i

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         do i = lo(1), hi(1)
            if (temp(i,1) .ge. temperr) &
               tag(i) = set
         enddo
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         do i = lo(1), hi(1)
            ax = ABS(temp(i+1,1) - temp(i,1))
            ax = MAX(ax,ABS(temp(i,1) - temp(i-1,1)))
            if ( ax .ge. tempgrad) &
               tag(i) = set
         enddo
      endif
      
      end
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the pressure
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: press     => pressure array
! ::: np        => number of components in press array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_presserror(tag,tagl1,tagh1,set,clear, &
                               press,pressl1,pressh1, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)
      use probdata_module
      implicit none
      integer set, clear, np, level
      integer tagl1,tagh1
      integer pressl1,pressh1
      integer lo(1), hi(1), domlo(1), domhi(1)
      integer tag(tagl1:tagh1)
      double precision press(pressl1:pressh1,np)
      double precision delta(1), xlo(1), problo(1), time

      double precision ax
      integer i

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
            do i = lo(1), hi(1)
               if (press(i,1) .ge. presserr) then
                  tag(i) = set
               endif
            enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
            do i = lo(1), hi(1)
               ax = ABS(press(i+1,1) - press(i,1))
               ax = MAX(ax,ABS(press(i,1) - press(i-1,1)))
               if ( ax .ge. pressgrad) then
                  tag(i) = set
               endif
            enddo
      endif

      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the velocity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vel       => velocity array
! ::: nv        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_velerror(tag,tagl1,tagh1, &
                             set,clear, &
                             vel,vell1,velh1, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none
      integer set, clear, nv, level
      integer tagl1,tagh1
      integer vell1,velh1
      integer lo(1), hi(1), domlo(1), domhi(1)
      integer tag(tagl1:tagh1)
      double precision vel(vell1:velh1,nv)
      double precision delta(1), xlo(1), problo(1), time

      double precision ax
      integer i

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         do i = lo(1), hi(1)
            ax = ABS(vel(i+1,1) - vel(i,1))
            ax = MAX(ax,ABS(vel(i,1) - vel(i-1,1)))
            if ( ax .ge. velgrad) then
               tag(i) = set
            endif
         enddo
      endif

      end



! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the radiation
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: rad       => radiation array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_raderror(tag,tagl1,tagh1, &
                             set,clear, &
                             rad,radl1,radh1, &
                             lo,hi,nr,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none
      integer set, clear, nr, level
      integer tagl1,tagh1
      integer radl1,radh1
      integer lo(1), hi(1), domlo(1), domhi(1)
      integer tag(tagl1:tagh1)
      double precision rad(radl1:radh1,nr)
      double precision delta(1), xlo(1), problo(1), time

      double precision ax
      integer i

!     Tag on regions of high radiation
      if (level .lt. max_raderr_lev) then
         do i = lo(1), hi(1)
            if (rad(i,1) .ge. raderr) then
               tag(i) = set
            endif
         enddo
      endif

!     Tag on regions of high radiation gradient
      if (level .lt. max_radgrad_lev) then
         do i = lo(1), hi(1)
            ax = ABS(rad(i+1,1) - rad(i,1))
            ax = MAX(ax,ABS(rad(i,1) - rad(i-1,1)))
            if ( ax .ge. radgrad) then
               tag(i) = set
            endif
         enddo
      endif

      end


