      subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,&
                            domlo,domhi,delta,xlo,time,bc)
 
      use eos_module
      use probdata_module
      use meth_params_module, only: NVAR,URHO,UMX,UMY,UFS,UFX,UEINT,UEDEN,UTEMP,&
                                    outflow_data_old, outflow_data_new, &
                                    outflow_data_old_time, outflow_data_new_time
 
      implicit none
      include 'bc_types.fi'
      integer adv_l1,adv_l2,adv_h1,adv_h2
      integer bc(2,2,*)
      integer domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

      integer          :: i, j, index, n
      integer          :: numpts
      double precision :: x,y,r
      double precision :: alpha,omalpha,dt,eps
      double precision :: dr,cen,xi,ghi,gmd,glo,den,maxvar,minvar
      double precision, allocatable :: radial(:,:)
      double precision staten(NVAR)
      type(eos_t) :: eos_state
 
      eps = 1.d-10

      if ( ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1) ) .or. &
           ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1) ) .or. &
           ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2) ) .or. &
           ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2) ) ) then
 
         dt = outflow_data_new_time - outflow_data_old_time

         if (time == outflow_data_old_time) then
            numpts = size(outflow_data_old,dim=2)
            allocate(radial(NVAR,0:numpts-1))
            do i = 1, numpts
               do n = 1,NVAR
                  radial(n,i-1) = outflow_data_old(n,i)
               end do
            end do
         else if (time == outflow_data_new_time) then
            numpts = size(outflow_data_new,dim=2)
            allocate(radial(NVAR,0:numpts-1))
            do i = 1, numpts
               do n = 1,NVAR
                  radial(n,i-1) = outflow_data_new(n,i)
               end do
            end do
         else if (time > outflow_data_new_time*(1.d0+eps) .or. time < outflow_data_old_time*(1.d0-eps)) then
               print *,'BAD TIME IN HYPFILL ',time
               print *,'OLD_TIME IS ',outflow_data_old_time
               print *,'NEW_TIME IS ',outflow_data_new_time
               stop
         else
            numpts = size(outflow_data_new,dim=2)
            if (numpts .ne. size(outflow_data_old,dim=2)) then
               print *,'SIZE OF OUTFLOW_DATA_OLD .NE. SIZE OF OUTFLOW_DATA_NEW IN HYPFILL'
               print *,'SIZE OF OUTFLOW_DATA_OLD ',size(outflow_data_old,dim=2)
               print *,'SIZE OF OUTFLOW_DATA_NEW ',size(outflow_data_new,dim=2)
               stop
            end if
            allocate(radial(NVAR,0:numpts-1))
            alpha = (time - outflow_data_old_time) / dt
            omalpha = 1.d0 - alpha 
            do i = 1, numpts
               do n = 1,NVAR
                  radial(n,i-1) = omalpha*outflow_data_old(n,i) + alpha*outflow_data_new(n,i)
               end do
            end do
         end if

         ! For the purposes of the outflow boundary fill, dr = dx
         dr = delta(1)

      end if

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,n), &
              adv_l1,adv_l2,adv_h1,adv_h2,&
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

      if ( (bc(1,1,1).eq. EXT_DIR .or. bc(1,2,1).eq. EXT_DIR) .or.  &
           (bc(2,1,1).eq. EXT_DIR .or. bc(2,2,1).eq. EXT_DIR) ) then
         print *,'NOT SET UP FOR EXT_DIR BCs IN HYPFILL'
         stop
      end if

      staten(URHO ) = rho_bndry
      staten(UMX  ) = rho_bndry * v_bndry
      staten(UMY  ) = rho_bndry * v_bndry

      ! Convert Kelvin to MeV
!      staten(UTEMP) = T_bndry / 1.160445d10
      staten(UTEMP) = T_bndry

      ! Set default species concentration to 1.
      staten(UFS  ) = 1.d0

      ! Use Ye from initial model
      staten(UFX  ) = Ye_bndry

      eos_state % rho = staten(URHO)
      eos_state % T   = staten(UTEMP)
      eos_state % xn  = staten(UFS)
      eos_state % aux = staten(UFX)

      call eos(eos_input_rt, eos_state)

      ! Convert e to rho*e
      staten(UEINT) = staten(URHO) * eos_state % e

      ! Convert e to rho*E = rho*(e + 1/2 u^2)
      staten(UEDEN) = staten(UEINT) + 0.5d0*staten(URHO)*(v_bndry**2)

      ! Convert Y to (rho*Y)
      staten(UFS)   = staten(URHO) * staten(UFS)

      ! Convert Ye to (rho*Ye)
      staten(UFX)   = staten(URHO) * staten(UFX)

!     XLO
      if ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1)) then
         do j = adv_l2, adv_h2
            do i = adv_l1, domlo(1)-1
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'HYPFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  do n = 1,NVAR
                     ghi = radial(n,index+1)
                     gmd = radial(n,index  )
                     glo = radial(n,index-1)
                     den = &
                      ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                      ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                      (-ghi + 26.d0*gmd - glo)/24.d0
      
                     minvar = min(gmd, min(glo,ghi))
                     maxvar = max(gmd, max(glo,ghi))
                     den = max(den,minvar)
                     den = min(den,maxvar)
                     adv(i,j,n) = den
                  end do
               else 
!                  index = numpts-1
!                  adv(i,j,1:NVAR) = radial(1:NVAR,index)
                   adv(i,j,1:NVAR) = staten(1:NVAR)
               end if

               ! The radial momentum was stored in both the UMX and UMY components --
               !   now we project it back along the normal
               adv(i,j,UMX) = adv(i,j,UMX) * (x/r)
               adv(i,j,UMY) = adv(i,j,UMY) * (y/r)

               if (adv(i,j,1) <= 0.0) then
                   print *,"XLO:BAD DENSITY AT ",i,j 
                   print *,"ADV ",adv(i,j,1)
                   print *,"INDEX ",index, radial(1,index)
               end if

            end do
         end do
      end if            

!     XHI
      if ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1)) then
         do j = adv_l2, adv_h2
            do i = domhi(1)+1, adv_h1
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'HYPFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  do n = 1,NVAR
                     ghi = radial(n,index+1)
                     gmd = radial(n,index  )
                     glo = radial(n,index-1)
                     den = &
                      ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                      ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                      (-ghi + 26.d0*gmd - glo)/24.d0
      
                     minvar = min(gmd, min(glo,ghi))
                     maxvar = max(gmd, max(glo,ghi))
                     den = max(den,minvar)
                     den = min(den,maxvar)
                     adv(i,j,n) = den
                  end do
               else 
!                  index = numpts-1
!                  adv(i,j,1:NVAR) = radial(1:NVAR,index)
                   adv(i,j,1:NVAR) = staten(1:NVAR)
               end if

               ! The radial momentum was stored in both the UMX and UMY components --
               !   now we project it back along the normal
               adv(i,j,UMX) = adv(i,j,UMX) * (x/r)
               adv(i,j,UMY) = adv(i,j,UMY) * (y/r)

               if (adv(i,j,1) <= 0.0) then
                   print *,"XHI:BAD DENSITY AT ",i,j 
                   print *,"ADV ",adv(i,j,1)
                   print *,"INDEX ",index, radial(1,index)
               end if

            end do
         end do
      end if            

!     YLO
      if ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, domlo(2)-1
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'HYPFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  do n = 1,NVAR
                     ghi = radial(n,index+1)
                     gmd = radial(n,index  )
                     glo = radial(n,index-1)
                     den = &
                      ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                      ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                      (-ghi + 26.d0*gmd - glo)/24.d0
      
                     minvar = min(gmd, min(glo,ghi))
                     maxvar = max(gmd, max(glo,ghi))
                     den = max(den,minvar)
                     den = min(den,maxvar)
                     adv(i,j,n) = den
                  end do
               else 
!                  index = numpts-1
!                  adv(i,j,1:NVAR) = radial(1:NVAR,index)
                   adv(i,j,1:NVAR) = staten(1:NVAR)
               end if

               ! The radial momentum was stored in both the UMX and UMY components --
               !   now we project it back along the normal
               adv(i,j,UMX) = adv(i,j,UMX) * (x/r)
               adv(i,j,UMY) = adv(i,j,UMY) * (y/r)

               if (adv(i,j,1) <= 0.0) then
                   print *,"YLO:BAD DENSITY AT ",i,j 
                   print *,"ADV ",adv(i,j,1)
                   print *,"INDEX ",index, radial(1,index)
               end if
            end do
         end do
      end if            

!     YHI
      if ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2)) then
         do i = adv_l1, adv_h1
            do j = domhi(2)+1, adv_h2
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'HYPFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  do n = 1,NVAR
                     ghi = radial(n,index+1)
                     gmd = radial(n,index  )
                     glo = radial(n,index-1)
                     den = &
                      ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                      ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                      (-ghi + 26.d0*gmd - glo)/24.d0
      
                     minvar = min(gmd, min(glo,ghi))
                     maxvar = max(gmd, max(glo,ghi))
                     den = max(den,minvar)
                     den = min(den,maxvar)
                     adv(i,j,n) = den
                  end do
               else 
!                  index = numpts-1
!                  adv(i,j,1:NVAR) = radial(1:NVAR,index)
                   adv(i,j,1:NVAR) = staten(1:NVAR)
               end if

               ! The radial momentum was stored in both the UMX and UMY components --
               !   now we project it back along the normal
               adv(i,j,UMX) = adv(i,j,UMX) * (x/r)
               adv(i,j,UMY) = adv(i,j,UMY) * (y/r)

               if (adv(i,j,1) <= 0.0) then
                   print *,"YHI:BAD DENSITY AT ",i,j 
                   print *,"ADV ",adv(i,j,1)
                   print *,"INDEX ",index, radial(1,index)
               end if
            end do
         end do
      end if            


      if ( ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1) ) .or. &
           ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1) ) .or. &
           ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2) ) .or. &
           ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2) ) ) then
         deallocate(radial)
      end if

      end subroutine ca_hypfill

! ::: -----------------------------------------------------------

      subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,&
                            domlo,domhi,delta,xlo,time,bc)

      use probdata_module
      use meth_params_module, only: outflow_data_old, outflow_data_new, &
                                    outflow_data_old_time, outflow_data_new_time
      implicit none
      include 'bc_types.fi'

      integer :: adv_l1,adv_l2,adv_h1,adv_h2
      integer :: bc(2,2,*)
      integer :: domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)

      integer          :: i,j,index
      integer          :: numpts
      double precision :: x,y,r
      double precision :: alpha,omalpha,dt,eps
      double precision :: dr,cen,xi,ghi,gmd,glo,den,maxvar,minvar
      double precision, allocatable :: radial(:)

      eps = 1.d-10

      if ( ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1) ) .or. &
           ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1) ) .or. &
           ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2) ) .or. &
           ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2) ) ) then

         dt = outflow_data_new_time - outflow_data_old_time

         if (time == outflow_data_old_time) then
            numpts = size(outflow_data_old,dim=2)
            allocate(radial(0:numpts-1))
            do i = 1, numpts
               radial(i-1) = outflow_data_old(1,i)
            end do
         else if (time == outflow_data_new_time) then
            numpts = size(outflow_data_new,dim=2)
            allocate(radial(0:numpts-1))
            do i = 1, numpts
               radial(i-1) = outflow_data_new(1,i)
            end do
         else if (time > outflow_data_new_time*(1.d0+eps) .or. time < outflow_data_old_time*(1.d0-eps)) then
               print *,'BAD TIME IN DENFILL ',time
               print *,'OLD_TIME IS ',outflow_data_old_time
               print *,'NEW_TIME IS ',outflow_data_new_time
               stop
         else
            numpts = size(outflow_data_new,dim=2)
            if (numpts .ne. size(outflow_data_old,dim=2)) then
               print *,'SIZE OF OUTFLOW_DATA_OLD .NE. SIZE OF OUTFLOW_DATA_NEW IN HYPFILL'
               print *,'SIZE OF OUTFLOW_DATA_OLD ',size(outflow_data_old,dim=2)
               print *,'SIZE OF OUTFLOW_DATA_NEW ',size(outflow_data_new,dim=2)
               stop
            end if
            allocate(radial(0:numpts-1))
            alpha   = (time - outflow_data_old_time) / (outflow_data_new_time - outflow_data_old_time)
            omalpha = 1.d0 - alpha
            do i = 1, numpts
               radial(i-1) = omalpha * outflow_data_old(1,i) + alpha * outflow_data_new(1,i)
            end do
         end if

         ! For the purposes of the outflow boundary fill, dr = dx
         dr = delta(1)

      end if

      call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,bc)

      if ( (bc(1,1,1).eq. EXT_DIR .or. bc(1,2,1).eq. EXT_DIR) .or.  &
           (bc(2,1,1).eq. EXT_DIR .or. bc(2,2,1).eq. EXT_DIR) ) then
         print *,'NOT SET UP FOR EXT_DIR BCs IN DENFILL '
         stop
      end if

!     XLO
      if ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1)) then
         do j = adv_l2, adv_h2
            do i = adv_l1, domlo(1)-1
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'DENFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  ghi = radial(index+1)
                  gmd = radial(index  )
                  glo = radial(index-1)
                  den = &
                   ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                   ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                   (-ghi + 26.d0*gmd - glo)/24.d0
     
                  minvar = min(gmd, min(glo,ghi))
                  maxvar = max(gmd, max(glo,ghi))
                  den = max(den,minvar)
                  den = min(den,maxvar)
                  adv(i,j) = den
               else 
!                  index = numpts-1
!                  adv(i,j) = radial(index)
                   adv(i,j) = rho_bndry
               end if
            end do
         end do
      end if            

!     XHI
      if ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1)) then
         do j = adv_l2, adv_h2
            do i = domhi(1)+1, adv_h1
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'DENFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  ghi = radial(index+1)
                  gmd = radial(index  )
                  glo = radial(index-1)
                  den = &
                   ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                   ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                   (-ghi + 26.d0*gmd - glo)/24.d0
     
                  minvar = min(gmd, min(glo,ghi))
                  maxvar = max(gmd, max(glo,ghi))
                  den = max(den,minvar)
                  den = min(den,maxvar)
                  adv(i,j) = den
               else 
!                  index = numpts-1
!                  adv(i,j) = radial(index)
                   adv(i,j) = rho_bndry
               end if
            end do
         end do
      end if            

!     YLO
      if ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, domlo(2)-1
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'DENFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  ghi = radial(index+1)
                  gmd = radial(index  )
                  glo = radial(index-1)
                  den = &
                   ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                   ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                   (-ghi + 26.d0*gmd - glo)/24.d0
     
                  minvar = min(gmd, min(glo,ghi))
                  maxvar = max(gmd, max(glo,ghi))
                  den = max(den,minvar)
                  den = min(den,maxvar)
                  adv(i,j) = den
               else 
!                  index = numpts-1
!                  adv(i,j) = radial(index)
                   adv(i,j) = rho_bndry
               end if
            end do
         end do
      end if            

!     YHI
      if ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2)) then
         do i = adv_l1, adv_h1
            do j = domhi(2)+1, adv_h2
               x = (dble(i) + 0.5d0) * delta(1) + corner(1) - center(1)
               y = (dble(j) + 0.5d0) * delta(2) + corner(2) - center(2)
               r = sqrt(x**2 + y**2)
               index = int(r/dr)
               if (index < 1) then
                   print *,'DENFILL: INDEX TOO SMALL ',index,' < 1'
                   print *,'AT (i,j) ',i,j
                   print *,'R / DR ',r,dr
                   stop
               else if (index < numpts-1) then
                  cen = (dble(index)+0.5d0)*dr
                   xi = r - cen
                  ! Quadratic interpolation
                  ghi = radial(index+1)
                  gmd = radial(index  )
                  glo = radial(index-1)
                  den = &
                   ( ghi -  2.d0*gmd + glo)*xi**2/(2.d0*dr**2) + &
                   ( ghi       - glo      )*xi   /(2.d0*dr   ) + &
                   (-ghi + 26.d0*gmd - glo)/24.d0
     
                  minvar = min(gmd, min(glo,ghi))
                  maxvar = max(gmd, max(glo,ghi))
                  den = max(den,minvar)
                  den = min(den,maxvar)
                  adv(i,j) = den
               else 
!                  index = numpts-1
!                  adv(i,j) = radial(index)
                   adv(i,j) = rho_bndry
               end if
            end do
         end do
      end if            

      if ( ( bc(1,1,1).eq. FOEXTRAP .and. adv_l1.lt.domlo(1) ) .or. &
           ( bc(1,2,1).eq. FOEXTRAP .and. adv_h1.gt.domhi(1) ) .or. &
           ( bc(2,1,1).eq. FOEXTRAP .and. adv_l2.lt.domlo(2) ) .or. &
           ( bc(2,2,1).eq. FOEXTRAP .and. adv_h2.gt.domhi(2) ) ) then

         deallocate(radial)
      end if

      end subroutine ca_denfill

! ::: -----------------------------------------------------------

      subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
           domlo,domhi,delta,xlo,time,bc)

      use probdata_module
      implicit none
      include 'bc_types.fi'

      integer :: grav_l1,grav_l2,grav_h1,grav_h2
      integer :: bc(2,2,*)
      integer :: domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

      call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

      end subroutine ca_gravxfill

! ::: -----------------------------------------------------------

      subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
           domlo,domhi,delta,xlo,time,bc)

      use probdata_module
      implicit none
      include 'bc_types.fi'

      integer :: grav_l1,grav_l2,grav_h1,grav_h2
      integer :: bc(2,2,*)
      integer :: domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision grav(grav_l1:grav_h1,grav_l2:grav_h2)

      call filcc(grav,grav_l1,grav_l2,grav_h1,grav_h2,domlo,domhi,delta,xlo,bc)

      end subroutine ca_gravyfill

! ::: -----------------------------------------------------------

      subroutine ca_radfill(rad,rad_l1,rad_l2,rad_h1,rad_h2, &
           domlo,domhi,delta,xlo,time,bc)

      use probdata_module
      implicit none
      include 'bc_types.fi'

      integer :: rad_l1,rad_l2,rad_h1,rad_h2
      integer :: bc(2,2,*)
      integer :: domlo(2), domhi(2)
      double precision delta(2), xlo(2), time
      double precision rad(rad_l1:rad_h1,rad_l2:rad_h2)

      call filcc(rad,rad_l1,rad_l2,rad_h1,rad_h2,domlo,domhi,delta,xlo,bc)

    end subroutine ca_radfill
