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
      subroutine ca_laplac_error(tag,taglo,taghi, &
                                 set,clear, &
                                 var,varlo,varhi, &
                                 lo,hi,nd,domlo,domhi, &
                                 delta,xlo,problo,time,level)

      use tagging_params_module
      use prob_params_module, only: dg, dim
      
      implicit none

      integer          :: set, clear, nd, level
      integer          :: taglo(3), taghi(3)
      integer          :: varlo(3), varhi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: var(varlo(1):varhi(1),varlo(2):varhi(2),varlo(3):varhi(3))
      double precision :: delta(3), xlo(3), problo(3), time

      integer          :: i, j, k
      double precision ::  delu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
      double precision :: delua(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
      double precision :: delu2(9), delu3(9), delu4(9)
      double precision :: num, denom, error

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8
      double precision, parameter :: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      delu = 0.0
      delua = 0.0
      
      ! d/dx
      do k=lo(3)-1*dg(3),hi(3)+1*dg(3)
         do j=lo(2)-1*dg(2),hi(2)+1*dg(2)
            do i=lo(1)-1*dg(1),hi(1)+1*dg(1)
               delu(i,j,k,1)  =     var(i+1*dg(1),j,k)  -     var(i-1*dg(1),j,k) 
               delua(i,j,k,1) = abs(var(i+1*dg(1),j,k)) + abs(var(i-1*dg(1),j,k))
            end do
         end do
      end do

      ! d/dy
      if (dim .ge. 2) then
         do k=lo(3)-1*dg(3),hi(3)+1*dg(3)
            do j=lo(2)-1*dg(2),hi(2)+1*dg(2)
               do i=lo(1)-1*dg(1),hi(1)+1*dg(1)
                  delu(i,j,k,2)  =     var(i,j+1*dg(2),k)  -     var(i,j-1*dg(2),k) 
                  delua(i,j,k,2) = abs(var(i,j+1*dg(2),k)) + abs(var(i,j-1*dg(2),k))
               end do
            end do
         end do
      endif
         
      ! d/dz
      if (dim .eq. 3) then
         do k=lo(3)-1*dg(3),hi(3)+1*dg(3)
            do j=lo(2)-1*dg(2),hi(2)+1*dg(2)
               do i=lo(1)-1*dg(1),hi(1)+1*dg(1)
                  delu(i,j,k,3)  =     var(i,j,k+1*dg(3))  -     var(i,j,k-1*dg(3))
                  delua(i,j,k,3) = abs(var(i,j,k+1*dg(3))) + abs(var(i,j,k-1*dg(3)))
               end do
            end do
         end do
      endif

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
         
               ! d/dxdx
               delu2(1) =     delu(i+1,j,k,1)  -     delu(i-1,j,k,1)
               delu3(1) = abs(delu(i+1,j,k,1)) + abs(delu(i-1,j,k,1))
               delu4(1) =    delua(i+1,j,k,1)  +    delua(i-1,j,k,1)

               ! d/dydx
               delu2(2) =     delu(i,j+1,k,1)  -     delu(i,j-1,k,1)
               delu3(2) = abs(delu(i,j+1,k,1)) + abs(delu(i,j-1,k,1))
               delu4(2) =    delua(i,j+1,k,1)  +    delua(i,j-1,k,1)

               ! d/dxdy
               delu2(3) =     delu(i+1,j,k,2)  -     delu(i-1,j,k,2)
               delu3(3) = abs(delu(i+1,j,k,2)) + abs(delu(i-1,j,k,2))
               delu4(3) =    delua(i+1,j,k,2)  +    delua(i-1,j,k,2)
                                       
               ! d/dydy                     
               delu2(4) =     delu(i,j+1,k,2)  -     delu(i,j-1,k,2)
               delu3(4) = abs(delu(i,j+1,k,2)) + abs(delu(i,j-1,k,2))
               delu4(4) =    delua(i,j+1,k,2)  +    delua(i,j-1,k,2)
                                                              
               ! d/dzdx                                            
               delu2(5) =     delu(i,j,k+1,1)  -     delu(i,j,k-1,1)
               delu3(5) = abs(delu(i,j,k+1,1)) + abs(delu(i,j,k-1,1))
               delu4(5) =    delua(i,j,k+1,1)  +    delua(i,j,k-1,1)
                                                              
               ! d/dzdy                                            
               delu2(6) =     delu(i,j,k+1,2)  -     delu(i,j,k-1,2)
               delu3(6) = abs(delu(i,j,k+1,2)) + abs(delu(i,j,k-1,2))
               delu4(6) =    delua(i,j,k+1,2)  +    delua(i,j,k-1,2)
                                                              
               ! d/dxdz                                            
               delu2(7) =     delu(i+1,j,k,3)  -     delu(i-1,j,k,3)
               delu3(7) = abs(delu(i+1,j,k,3)) + abs(delu(i-1,j,k,3))
               delu4(7) =    delua(i+1,j,k,3)  +    delua(i-1,j,k,3)
                                                              
               ! d/dydz                                            
               delu2(8) =     delu(i,j+1,k,3)  -     delu(i,j-1,k,3)
               delu3(8) = abs(delu(i,j+1,k,3)) + abs(delu(i,j-1,k,3))
               delu4(8) =    delua(i,j+1,k,3)  +    delua(i,j-1,k,3)
               
               ! d/dzdz                                            
               delu2(9) =     delu(i,j,k+1,3)  -     delu(i,j,k-1,3)
               delu3(9) = abs(delu(i,j,k+1,3)) + abs(delu(i,j,k-1,3))
               delu4(9) =    delua(i,j,k+1,3)  +    delua(i,j,k-1,3)

               ! compute the error
               num   = sum(delu2**2)

               denom = sum((delu3 + (epsil*delu4+1.d-99))**2)

               error = sqrt(num/denom)

               if (error .gt. ctore) tag(i,j,k) = set

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
      subroutine ca_denerror(tag,taglo,taghi, &
                             set,clear, &
                             den,denlo,denhi, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, nd, level
      integer          :: taglo(3), taghi(3)
      integer          :: denlo(3), denhi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: den(denlo(1):denhi(1),denlo(2):denhi(2),denlo(3):denhi(3),nd)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (den(i,j,k,1) .ge. denerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(den(i+1*dg(1),j,k,1) - den(i,j,k,1))
                  ay = ABS(den(i,j+1*dg(2),k,1) - den(i,j,k,1))
                  az = ABS(den(i,j,k+1*dg(3),1) - den(i,j,k,1))
                  ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. dengrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_denerror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the temperature
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

      subroutine ca_temperror(tag,taglo,taghi, &
                              set,clear, &
                              temp,templo,temphi, &
                              lo,hi,np,domlo,domhi, &
                              delta,xlo,problo,time,level)

      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, np, level
      integer          :: taglo(3), taghi(3)
      integer          :: templo(3), temphi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: temp(templo(1):temphi(1),templo(2):temphi(2),templo(3):temphi(3),np)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (temp(i,j,k,1) .ge. temperr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(temp(i+1*dg(1),j,k,1) - temp(i,j,k,1))
                  ay = ABS(temp(i,j+1*dg(2),k,1) - temp(i,j,k,1))
                  az = ABS(temp(i,j,k+1*dg(3),1) - temp(i,j,k,1))
                  ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. tempgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_temperror

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
      subroutine ca_presserror(tag,taglo,taghi, &
                               set,clear, &
                               press,presslo,presshi, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)

      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, np, level
      integer          :: taglo(3), taghi(3)
      integer          :: presslo(3), presshi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: press(presslo(1):presshi(1),presslo(2):presshi(2),presslo(3):presshi(3),np)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (press(i,j,k,1) .ge. presserr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(press(i+1*dg(1),j,k,1) - press(i,j,k,1))
                  ay = ABS(press(i,j+1*dg(2),k,1) - press(i,j,k,1))
                  az = ABS(press(i,j,k+1*dg(3),1) - press(i,j,k,1))
                  ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. pressgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_presserror

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
! ::: nv        => number of components in vel array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_velerror(tag,taglo,taghi, &
                             set,clear, &
                             vel,vello,velhi, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, nv, level
      integer          :: taglo(3), taghi(3)
      integer          :: vello(3), velhi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: vel(vello(1):velhi(1),vello(2):velhi(2),vello(3):velhi(3),nv)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high velocity
      if (level .lt. max_velerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (vel(i,j,k,1) .ge. velerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif      
      
!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(vel(i+1*dg(1),j,k,1) - vel(i,j,k,1))
                  ay = ABS(vel(i,j+1*dg(2),k,1) - vel(i,j,k,1))
                  az = ABS(vel(i,j,k+1*dg(3),1) - vel(i,j,k,1))
                  ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. velgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_velerror


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
      subroutine ca_raderror(tag,taglo,taghi, &
                             set,clear, &
                             rad,radlo,radhi, &
                             lo,hi,nr,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, nr, level
      integer          :: taglo(3), taghi(3)
      integer          :: radlo(3), radhi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: rad(radlo(1):radhi(1),radlo(2):radhi(2),radlo(3):radhi(3),nr)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high radiation
      if (level .lt. max_raderr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (rad(i,j,k,1) .ge. raderr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high radiation gradient
      if (level .lt. max_radgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(rad(i+1*dg(1),j,k,1) - rad(i,j,k,1))
                  ay = ABS(rad(i,j+1*dg(2),k,1) - rad(i,j,k,1))
                  az = ABS(rad(i,j,k+1*dg(3),1) - rad(i,j,k,1))
                  ax = MAX(ax,ABS(rad(i,j,k,1) - rad(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(rad(i,j,k,1) - rad(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(rad(i,j,k,1) - rad(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. radgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_raderror

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the entropy
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: ent       => entropy array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_enterror(tag,taglo,taghi, &
                             set,clear, &
                             ent,entlo,enthi, &
                             lo,hi,ne,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, ne, level
      integer          :: taglo(3), taghi(3)
      integer          :: entlo(3), enthi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: ent(entlo(1):enthi(1),entlo(2):enthi(2),entlo(3):enthi(3),ne)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high entropy
      if (level .lt. max_enterr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (ent(i,j,k,1) .ge. enterr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high entropy gradient
      if (level .lt. max_entgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(ent(i+1*dg(1),j,k,1) - ent(i,j,k,1))
                  ay = ABS(ent(i,j+1*dg(2),k,1) - ent(i,j,k,1))
                  az = ABS(ent(i,j,k+1*dg(3),1) - ent(i,j,k,1))
                  ax = MAX(ax,ABS(ent(i,j,k,1) - ent(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(ent(i,j,k,1) - ent(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(ent(i,j,k,1) - ent(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. entgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_enterror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the Ye
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: ye        => ye array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_yeerror(tag,taglo,taghi, &
                            set,clear, &
                            ye,yelo,yehi, &
                            lo,hi,nr,domlo,domhi, &
                            delta,xlo,problo,time,level)

      use probdata_module
      use tagging_params_module
      use prob_params_module, only: dg
      
      implicit none

      integer          :: set, clear, nr, level
      integer          :: taglo(3), taghi(3)
      integer          :: yelo(3), yehi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
      double precision :: ye ( yelo(1): yehi(1), yelo(2): yehi(2), yelo(3): yehi(3),nr)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision :: ax, ay, az
      integer          :: i, j, k

!     Tag on regions of high Ye
      if (level .lt. max_yeerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (ye(i,j,k,1) .ge. yeerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high ye gradient
      if (level .lt. max_yegrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(ye(i+1*dg(1),j,k,1) - ye(i,j,k,1))
                  ay = ABS(ye(i,j+1*dg(2),k,1) - ye(i,j,k,1))
                  az = ABS(ye(i,j,k+1*dg(3),1) - ye(i,j,k,1))
                  ax = MAX(ax,ABS(ye(i,j,k,1) - ye(i-1*dg(1),j,k,1)))
                  ay = MAX(ay,ABS(ye(i,j,k,1) - ye(i,j-1*dg(2),k,1)))
                  az = MAX(az,ABS(ye(i,j,k,1) - ye(i,j,k-1*dg(3),1)))
                  if ( MAX(ax,ay,az) .ge. yegrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine ca_yeerror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the mass
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of work region
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: grav      => grav array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of work region
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_graverror(tag,taglo,taghi, &
                              set,clear, &
                              grav,gravlo,gravhi, &
                              lo,hi,nr,domlo,domhi, &
                              delta,xlo,problo,time,level)

      use probdata_module
      use tagging_params_module
      use prob_params_module, only: dim, center
      use fundamental_constants_module, only : Gconst
      
      implicit none

      integer          :: set, clear, nr, level
      integer          :: taglo(3), taghi(3)
      integer          :: gravlo(3), gravhi(3)
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag ( taglo(1): taghi(1), taglo(2): taghi(2), taglo(3): taghi(3))
      double precision :: grav(gravlo(1):gravhi(1),gravlo(2):gravhi(2),gravlo(3):gravhi(3),nr)
      double precision :: delta(3), xlo(3), problo(3), time

      double precision, parameter :: Msun = 1.98892d33 
      double precision mass, x, y, z, r2
      integer          :: i, j, k

      if (level .lt. max_masserr_lev) then
         if (dim .eq. 1) then
            j = lo(2)
            k = lo(3)
            do i = lo(1), hi(1)
               x = xlo(1) + (i-lo(1)+0.5d0)*delta(1)
               mass = -grav(i,j,k,1) * x**2 / (Gconst * Msun)
               if (mass .lt. masserr .and. mass .gt. 0.d0 .or. time < 1.d-10) then
                  tag(i,j,k) = set
               endif
            enddo
         else if (dim .eq. 2) then
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  x = xlo(1) + (i-lo(1)+0.5d0)*delta(1) - center(1)
                  y = xlo(2) + (j-lo(2)+0.5d0)*delta(2) - center(2)
                  r2 = x**2 + y**2
                  mass = grav(i,j,k,1) * r2 / (Gconst * Msun) ! grav is mag. grav
                  if (mass .lt. masserr .and. mass .gt. 0.d0) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         else
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     x = xlo(1) + (i-lo(1)+0.5d0)*delta(1) - center(1)
                     y = xlo(2) + (j-lo(2)+0.5d0)*delta(2) - center(2)
                     z = xlo(3) + (k-lo(3)+0.5d0)*delta(3) - center(3)
                     r2 = x**2 + y**2 + z**2
                     mass = grav(i,j,k,1) * r2 / (Gconst * Msun)
                     if (mass .lt. masserr .and. mass .gt. 0.d0) then
                        tag(i,j,k) = set
                     endif
                  enddo
               enddo
            enddo
         endif
      end if

    end subroutine ca_graverror
