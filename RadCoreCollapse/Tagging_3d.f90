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

      double precision ::  delu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
      double precision :: delua(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
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
          delu(i,j,k,1) =     var(i+1,j,k)  -     var(i-1,j,k) 
         delua(i,j,k,1) = abs(var(i+1,j,k)) + abs(var(i-1,j,k))
      end do
      end do
      end do

      ! d/dy
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,k,2) =     var(i,j+1,k)  -     var(i,j-1,k) 
         delua(i,j,k,2) = abs(var(i,j+1,k)) + abs(var(i,j-1,k))
      end do
      end do
      end do

      ! d/dz
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,k,3) =     var(i,j,k+1)  -     var(i,j,k-1)
         delua(i,j,k,3) = abs(var(i,j,k+1)) + abs(var(i,j,k-1))
      end do
      end do
      end do

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

      subroutine ca_temperror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              temp,templ1,templ2,templ3,temph1, &
                              temph2,temph3, &
                              lo,hi,np,domlo,domhi, &
                              delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer templ1,templ2,templ3,temph1,temph2,temph3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision temp(templ1:temph1,templ2:temph2, &
                            templ3:temph3,np)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

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
                  ax = ABS(temp(i+1,j,k,1) - temp(i,j,k,1))
                  ay = ABS(temp(i,j+1,k,1) - temp(i,j,k,1))
                  az = ABS(temp(i,j,k+1,1) - temp(i,j,k,1))
                  ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1,j,k,1)))
                  ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1,k,1)))
                  az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. tempgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
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
      subroutine ca_presserror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             press,pressl1,pressl2,pressl3,pressh1, &
                             pressh2,pressh3, &
                             lo,hi,np,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer pressl1,pressl2,pressl3,pressh1,pressh2,pressh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision press(pressl1:pressh1,pressl2:pressh2, &
                             pressl3:pressh3,np)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

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
                  ax = ABS(press(i+1,j,k,1) - press(i,j,k,1))
                  ay = ABS(press(i,j+1,k,1) - press(i,j,k,1))
                  az = ABS(press(i,j,k+1,1) - press(i,j,k,1))
                  ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1,j,k,1)))
                  ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1,k,1)))
                  az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. pressgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
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
      if (level .lt. max_velgrad_lev) then
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
      subroutine ca_raderror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             rad,radl1,radl2,radl3,radh1,radh2,radh3, &
                             lo,hi,nr,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, nr, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer radl1,radl2,radl3,radh1,radh2,radh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision rad(radl1:radh1,radl2:radh2,radl3:radh3,nr)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

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
                  ax = ABS(rad(i+1,j,k,1) - rad(i,j,k,1))
                  ay = ABS(rad(i,j+1,k,1) - rad(i,j,k,1))
                  az = ABS(rad(i,j,k+1,1) - rad(i,j,k,1))
                  ax = MAX(ax,ABS(rad(i,j,k,1) - rad(i-1,j,k,1)))
                  ay = MAX(ay,ABS(rad(i,j,k,1) - rad(i,j-1,k,1)))
                  az = MAX(az,ABS(rad(i,j,k,1) - rad(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. radgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the entropy
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: ent       => entropy array
! ::: nd        => number of components in ent array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_enterror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             ent,entl1,entl2,entl3,enth1,enth2,enth3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer entl1,entl2,entl3,enth1,enth2,enth3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision ent(entl1:enth1,entl2:enth2,entl3:enth3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

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
            if (k.eq.domhi(3) .or. k.eq.domlo(3)) then
               cycle
            end if
            do j = lo(2), hi(2)
               if (j.eq.domhi(2) .or. j.eq.domlo(2)) then
                  cycle
               end if
               do i = lo(1), hi(1)
                  if (i.eq.domhi(1) .or. i.eq.domlo(1)) then
                     cycle
                  end if
                  ax = ABS(ent(i+1,j,k,1) - ent(i,j,k,1))
                  ay = ABS(ent(i,j+1,k,1) - ent(i,j,k,1))
                  az = ABS(ent(i,j,k+1,1) - ent(i,j,k,1))
                  ax = MAX(ax,ABS(ent(i,j,k,1) - ent(i-1,j,k,1)))
                  ay = MAX(ay,ABS(ent(i,j,k,1) - ent(i,j-1,k,1)))
                  az = MAX(az,ABS(ent(i,j,k,1) - ent(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. entgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on Ye
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: Ye        => Ye array
! ::: nd        => number of components in ent array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_yeerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             Ye,Yel1,Yel2,Yel3,Yeh1,Yeh2,Yeh3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer Yel1,Yel2,Yel3,Yeh1,Yeh2,Yeh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision Ye(Yel1:Yeh1,Yel2:Yeh2,Yel3:Yeh3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

!     Tag on regions of low Ye
      if (level .lt. max_Yeerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (Ye(i,j,k,1) .lt. Yeerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high Ye gradient
      if (level .lt. max_Yegrad_lev) then
         do k = lo(3), hi(3)
            if (k.eq.domhi(3) .or. k.eq.domlo(3)) then
               cycle
            end if
            do j = lo(2), hi(2)
               if (j.eq.domhi(2) .or. j.eq.domlo(2)) then
                  cycle
               end if
               do i = lo(1), hi(1)
                  if (i.eq.domhi(1) .or. i.eq.domlo(1)) then
                     cycle
                  end if
                  ax = ABS(Ye(i+1,j,k,1) - Ye(i,j,k,1))
                  ay = ABS(Ye(i,j+1,k,1) - Ye(i,j,k,1))
                  az = ABS(Ye(i,j,k+1,1) - Ye(i,j,k,1))
                  ax = MAX(ax,ABS(Ye(i,j,k,1) - Ye(i-1,j,k,1)))
                  ay = MAX(ay,ABS(Ye(i,j,k,1) - Ye(i,j-1,k,1)))
                  az = MAX(az,ABS(Ye(i,j,k,1) - Ye(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. Yegrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on mass
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: grav      => grav array
! ::: nd        => number of components in ent array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_graverror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             grav,gravl1,gravl2,gravl3,gravh1,gravh2, gravh3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      use fundamental_constants_module, only : Gconst
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagh1,tagl2,tagh2,tagl3,tagh3
      integer gravl1,gravh1,gravl2,gravh2,gravl3,gravh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision grav(gravl1:gravh1,gravl2:gravh2,gravl3:gravh3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision, parameter :: Msun = 1.98892d33 

      double precision mass, x, y, z, r2
      integer i, j, k

      if (level .lt. max_masserr_lev) then
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
      
      end

